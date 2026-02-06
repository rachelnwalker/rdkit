/// This is a simple test file for reference on how MonomerMol and SCSR could be compatible.
#include "RDGeneral/test.h"
#include <catch2/catch_all.hpp>

#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/MonomerMol/MonomerLibrary.h>
#include <GraphMol/MonomerMol/MonomerMol.h>
#include <GraphMol/MonomerMol/Conversions.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <fstream>
#include <sstream>
#include <regex>
#include <iostream>

using namespace RDKit;

namespace {

// Structure to hold template info extracted from SCSR
struct TemplateInfo {
    std::string symbol;           // e.g., "A", "G", "C"
    std::string molblock;         // V3000 molblock for the template
};

// Convert attachment type to atom map number: Al=1, Br=2, Cx=3
int attachTypeToMapNum(const std::string& attachType) {
    if (attachType == "Al") return 1;
    if (attachType == "Br") return 2;
    if (attachType == "Cx") return 3;
    throw std::runtime_error("Unknown attachment type: " + attachType);
}

// Extract templates from SCSR file content (just molblocks, not SAP parsing)
std::vector<TemplateInfo> extractTemplates(const std::string& fileContent) {
    std::vector<TemplateInfo> templates;

    // Find TEMPLATE section
    size_t templateStart = fileContent.find("M  V30 BEGIN TEMPLATE");
    size_t templateEnd = fileContent.find("M  V30 END TEMPLATE");
    if (templateStart == std::string::npos || templateEnd == std::string::npos) {
        return templates;
    }

    std::string templateSection = fileContent.substr(templateStart, templateEnd - templateStart);

    // Pattern to match template headers like "TEMPLATE 1 AA/Ala/A/"
    std::regex templateHeaderRe(R"(M\s+V30\s+TEMPLATE\s+\d+\s+\w+/\w+/(\w+)/)");

    // Find each template
    std::istringstream iss(templateSection);
    std::string line;
    TemplateInfo currentTemplate;
    bool inCtab = false;
    std::ostringstream molblockStream;

    while (std::getline(iss, line)) {
        std::smatch match;

        if (std::regex_search(line, match, templateHeaderRe)) {
            // New template header
            if (!currentTemplate.symbol.empty() && !currentTemplate.molblock.empty()) {
                templates.push_back(currentTemplate);
            }
            currentTemplate = TemplateInfo();
            currentTemplate.symbol = match[1].str();
            inCtab = false;
            molblockStream.str("");
            molblockStream.clear();
        }
        else if (line.find("M  V30 BEGIN CTAB") != std::string::npos) {
            inCtab = true;
            // Start fresh molblock with proper header
            molblockStream.str("");
            molblockStream.clear();
            molblockStream << "\n  RDKit          3D\n\n";
            molblockStream << "  0  0  0  0  0  0  0  0  0  0999 V3000\n";
            molblockStream << "M  V30 BEGIN CTAB\n";
        }
        else if (line.find("M  V30 END CTAB") != std::string::npos && inCtab) {
            molblockStream << "M  V30 END CTAB\n";
            molblockStream << "M  END\n";
            currentTemplate.molblock = molblockStream.str();
            inCtab = false;
        }
        else if (inCtab) {
            molblockStream << line << "\n";
        }
    }

    // Don't forget the last template
    if (!currentTemplate.symbol.empty() && !currentTemplate.molblock.empty()) {
        templates.push_back(currentTemplate);
    }

    return templates;
}

// Set atom map numbers on a molecule based on SubstanceGroup attach points
// Finds the SUP type SGroup with CLASS=AA and uses its attach points
void setAtomMapsFromSubstanceGroups(RWMol& mol) {
    for (const auto& sg : getSubstanceGroups(mol)) {
        // Look for SUP (superatom) type with CLASS=AA (the monomer itself)
        if (sg.getProp<std::string>("TYPE") != "SUP") {
            continue;
        }

        std::string sgClass;
        if (!sg.getPropIfPresent("CLASS", sgClass) || sgClass != "AA") {
            continue;
        }

        // Found the main monomer SGroup - get its attach points
        for (const auto& ap : sg.getAttachPoints()) {
            // ap.aIdx is the core atom index (0-based)
            // ap.id is the attachment type ("Al", "Br", "Cx")
            int mapNum = attachTypeToMapNum(ap.id);
            mol.getAtomWithIdx(ap.lvIdx)->setAtomMapNum(mapNum);
        }
    }
}

// Read entire file into string
std::string readFile(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + path);
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

} // anonymous namespace

TEST_CASE("SCSR to MonomerMol Conversion", "[MonomerMol][SCSR]") {
    std::string fName = std::string(getenv("RDBASE")) +
        "/Code/GraphMol/FileParsers/test_data/macromols/cyclic.mol";

    SECTION("Build MonomerMol from SCSR templates") {
        // Step 1: Read the SCSR file and extract template molblocks
        std::string fileContent = readFile(fName);
        auto templates = extractTemplates(fileContent);

        REQUIRE(templates.size() == 5);  // A, G, C, L, P

        // Step 2: Create empty custom MonomerLibrary
        auto customLib = std::make_shared<MonomerLibrary>();

        // Step 3: Parse each template molblock, set atom maps from SGroups, add to library
        for (const auto& tmpl : templates) {
            INFO("Processing template: " << tmpl.symbol);

            // Parse the V3000 molblock (preserves SubstanceGroups)
            bool sanitize = false;
            bool removeHs = false;
            bool strictParsing = false;
            std::unique_ptr<RWMol> mol(MolBlockToMol(tmpl.molblock, sanitize, removeHs, strictParsing));
            REQUIRE(mol != nullptr);

            // Set atom map numbers based on SubstanceGroup attach points
            // The parsed mol contains SGROUPs with SAP info
            // Al -> :1, Br -> :2, Cx -> :3
            setAtomMapsFromSubstanceGroups(*mol);

            std::cout << "Template " << tmpl.symbol << " has " << mol->getNumAtoms() << " atoms and "
                      << mol->getNumBonds() << " bonds." << std::endl;
            std::cout << "SMILES with atom maps: " << MolToSmiles(*mol) << std::endl;

            // Add to library using addMonomer with the parsed mol
            customLib->addMonomer(
                std::shared_ptr<ROMol>(mol.release()),
                tmpl.symbol,
                "AA",
                tmpl.symbol  // Use symbol as PDB code for simplicity
            );

            CHECK(customLib->hasMonomer(tmpl.symbol, "AA"));
        }

        // Step 4: Build MonomerMol with the sequence from SCSR
        // From cyclic.mol, sequence by SEQID: G, G, C, L, P, C, A, A
        MonomerMol monomerMol(customLib);

        // note that order this is added does not matter
        monomerMol.addMonomer("G", 1, "AA", "AA1");  // SEQID=1
        monomerMol.addMonomer("G");   // SEQID=2
        monomerMol.addMonomer("C");   // SEQID=3
        monomerMol.addMonomer("L");   // SEQID=4
        monomerMol.addMonomer("P");   // SEQID=5
        monomerMol.addMonomer("C");   // SEQID=6
        monomerMol.addMonomer("A");   // SEQID=7
        monomerMol.addMonomer("A");   // SEQID=8

        CHECK(monomerMol.getNumAtoms() == 8);

        // Step 5: Add connections based on ATTCHORD
        // Backbone (R2-R1)
        monomerMol.addConnection(0, 1, BACKBONE_LINKAGE);  // G(1) -> G(2)
        monomerMol.addConnection(1, 2, BACKBONE_LINKAGE);  // G(2) -> C(3)
        monomerMol.addConnection(2, 3, BACKBONE_LINKAGE);  // C(3) -> L(4)
        monomerMol.addConnection(3, 4, BACKBONE_LINKAGE);  // L(4) -> P(5)
        monomerMol.addConnection(4, 5, BACKBONE_LINKAGE);  // P(5) -> C(6)
        monomerMol.addConnection(5, 6, BACKBONE_LINKAGE);  // C(6) -> A(7)
        monomerMol.addConnection(6, 7, BACKBONE_LINKAGE);  // A(7) -> A(8)

        // Cyclic closure
        monomerMol.addConnection(7, 0, BACKBONE_LINKAGE);  // A(8) -> G(1)

        // Disulfide bridge (R3-R3)
        monomerMol.addConnection(2, 5, CROSS_LINKAGE);     // C(3) <-> C(6)

        CHECK(monomerMol.getNumBonds() == 9);

        std::cout << "MonomerMol as SMILES: " << MolToSmiles(monomerMol) << std::endl;

        // Step 6: Convert to atomistic
        auto atomisticFromMonomer = toAtomistic(monomerMol);
        REQUIRE(atomisticFromMonomer != nullptr);
        CHECK(atomisticFromMonomer->getNumAtoms() > 0);

        // Step 7: Compare with MolFromSCSRFile
        v2::FileParsers::MolFileParserParams pp;
        pp.strictParsing = false;
        v2::FileParsers::MolFromSCSRParams scsrParams;

        auto atomisticFromSCSR = v2::FileParsers::MolFromSCSRFile(fName, pp, scsrParams);
        REQUIRE(atomisticFromSCSR != nullptr);

        // Compare atom and bond counts
        CHECK(atomisticFromMonomer->getNumAtoms() == atomisticFromSCSR->getNumAtoms());
        CHECK(atomisticFromMonomer->getNumBonds() == atomisticFromSCSR->getNumBonds());

        // Compare canonical SMILES
        std::string smilesFromMonomer = MolToSmiles(*atomisticFromMonomer);
        std::string smilesFromSCSR = MolToSmiles(*atomisticFromSCSR);

        INFO("SMILES from MonomerMol: " << smilesFromMonomer);
        INFO("SMILES from SCSR:       " << smilesFromSCSR);

        CHECK(smilesFromMonomer == smilesFromSCSR);
    }
}
