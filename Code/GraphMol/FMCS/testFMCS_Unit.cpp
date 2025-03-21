//
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifdef _WIN32
#include <RDGeneral/test.h>
#else
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <cstdio>
#include <cstring>
#include <ctime>
#include <string>
#include <iostream>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include "../RDKitBase.h"
#include "../FileParsers/FileParsers.h"  //MOL single molecule !
#include "../FileParsers/MolSupplier.h"  //SDF
#include "../SmilesParse/SmilesParse.h"
#include "../SmilesParse/SmilesWrite.h"
#include "../SmilesParse/SmartsWrite.h"
#include "FMCS.h"
#include "DebugTrace.h"  //#ifdef VERBOSE_STATISTICS_ON

#include "../Substruct/SubstructMatch.h"

#ifdef _WIN32
#include <windows.h>
#endif

using namespace RDKit;

unsigned long long T0;
unsigned long long t0;

void printTime() {
  unsigned long long t1 = nanoClock();
  double sec = double(t1 - t0) / 1000000.;
  printf("Time elapsed %.3lf seconds\n", sec);
  t0 = nanoClock();
}

std::string getSmilesOnly(
    const char *smiles,
    std::string *id = nullptr) {  // remove label, because RDKit parse FAILED
  const char *sp = strchr(smiles, ' ');
  unsigned int n = (sp ? sp - smiles + 1 : strlen(smiles));
  if (id) {
    *id = std::string(smiles + n);
  }
  return std::string(smiles, n);
}

// UNIT Test Set:
//=========================================================================

void test1Basics() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "FMCS test1Basics()" << std::endl;

  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      "CC1CCC(N)CC1", "CC1CC(C)CC(C)C1",  // OK test.sdf
  };

  for (auto &i : smi) {
    std::string id;
    mols.emplace_back(SmilesToMol(getSmilesOnly(i, &id)));
  }
  t0 = nanoClock();
  MCSResult res = findMCS(mols);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 7 && res.NumBonds == 7);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test32() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS test32" << std::endl;
  std::cout << "\n test32()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      // #  Using CHEMBL1515359 CHEMBL1590658 CHEMBL1447567 CHEMBL1384017
      // CHEMBL1456416 CHEMBL1308819 CHEMBL1703007 CHEMBL1707819 CHEMBL1500793
      // CHEMBL1334715
      // 32 . 1 31 33 0.82
      // S(-N1-C-C-O-C-C-1)(-c1:c:c:c(-N(-C-C)-C-C):c(-N-C(-C=C-c2:c:c:c:c:c:2)=O):c:1)(=O)=O
      "O=C(Nc1cc(S(N2CCOCC2)(=O)=O)ccc1N1CCOCC1)C=Cc1ccc(Cl)cc1  CHEMBL1515359",
      "c1ccc(C=CC(Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCOCC2)=O)cc1  CHEMBL1590658",
      "Cc1ccc(C=CC(=O)Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCOCC2)cc1  CHEMBL1447567",
      "c1ccc(C=CC(Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCCC2)=O)cc1  CHEMBL1384017",
      "O=C(C=Cc1ccc(F)cc1)Nc1cc(S(N2CCOCC2)(=O)=O)ccc1N1CCCC1  CHEMBL1456416",
      "c1cc(F)cc(C=CC(=O)Nc2c(N3CCCC3)ccc(S(N3CCOCC3)(=O)=O)c2)c1  "
      "CHEMBL1308819",
      "CCN1CCN(c2ccc(S(N3CCOCC3)(=O)=O)cc2NC(=O)C=Cc2ccc(C)cc2)CC1  "
      "CHEMBL1703007",
      "c1cc(C=CC(=O)Nc2cc(S(N3CCOCC3)(=O)=O)ccc2N2CCOCC2)c([N+]([O-])=O)cc1  "
      "CHEMBL1707819",
      "N#CC(=Cc1ccccc1)C(=O)Nc1cc(S(N2CCOCC2)(=O)=O)ccc1N1CCCC1  CHEMBL1500793",
      "C(=Cc1ccc2c(c1)OCO2)C(Nc1cc(S(=O)(=O)N2CCOCC2)ccc1N1CCOCC1)=O  "
      "CHEMBL1334715",
      // 31 33 0.35 sec MCS: CCN(CC)c1ccc(cc1NC(=O)C=Cc1ccccc1)S(=O)(=O)N1CCOCC1
  };
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i)));
  }
  t0 = nanoClock();
  MCSResult res = findMCS(mols);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 31 && res.NumBonds == 33);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test190() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS test190" << std::endl;
  std::cout << "\n test190()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      // # 190
      "COc1cc2nc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)ccc3)oc2cc1  CHEMBL1479679",
      "COc1cc2nc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)c(C)cc3)oc2cc1  CHEMBL1333382",
      "Cc1cc2oc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)ccc3)nc2cc1  CHEMBL1437584",
      "COc1c(NC(=O)CSc2ccc(Cl)cc2)cc(-c2nc3ccccc3o2)cc1  CHEMBL1601350",
      "Cc1cc2nc(-c3cccc(NC(=O)CSc4ccc(Cl)cc4)c3)oc2cc1C  CHEMBL1398008",
      "Cc1cc2oc(-c3cc(NC(=O)CSc4ccc(Cl)cc4)c(C)cc3)nc2cc1  CHEMBL1612903",
      "COc1cc2nc(-c3cc(NC(=O)Cc4ccc(Cl)cc4)c(C)cc3)oc2cc1  CHEMBL1316483",
      "Cc1c(NC(=O)CSc2ccc(Cl)cc2)cccc1-c1nc2cc(Cl)ccc2o1  CHEMBL1568754",
      "COc1ccc2oc(-c3ccc(C)c(NC(=O)COc4cc(C)cc(C)c4)c3)nc2c1  CHEMBL1436972",
      "Cc1ccc(SCC(=O)Nc2cc(-c3nc4cc(C)ccc4o3)c(O)cc2)cc1  CHEMBL1611932",
      // # 19 21 1.37 sec MCS: CC(=O)Nc1cccc(c1)-c1nc2ccccc2o1
      //  19 21 2.36 sec MCS: CC(=O)Nc1cccc(c1)-c1nc2ccccc2o1 19 atoms, 21 bonds
  };
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i)));
  }
  t0 = nanoClock();
  MCSResult res = findMCS(mols);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 19 && res.NumBonds == 21);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test3() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS test3" << std::endl;
  std::cout << "\n test3()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      // TEST 3
      "CN(C)c1ccc(CC(=O)NCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL153934",
      "CN(C)c1ccc(CC(=O)NCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL152361",
      "CN(C)c1ccc(CC(=O)NCCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL157336",
      "CN(C)c1ccc(CC(=O)NCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL157429",
      "CN(C)c1ccc(CC(=O)NCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL357551",
      "CN(C)c1ccc(CC(=O)NCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL421974",
      "CN(C)c1ccc(CC(NCCCCCC(NO)=O)=O)cc1 CHEMBL484488",
      "CC(C)Cc1ccc(C(C)C(=O)NC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL564780",
      "c1cc([N+]([O-])=O)ccc1CC(=O)NC1CCCCCC1 CHEMBL1553142",
      "CC1(C)NC(C)(C)CC(NC(=O)Cc2ccccc2)C1 CHEMBL1703640",
      // # 3 . 1 14 14 0.08 sec MCS: CCCCNC(=O)Cc1ccccc1
  };
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i)));
  }
  t0 = nanoClock();
  MCSResult res = findMCS(mols);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 14 && res.NumBonds == 14);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testRing1() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testRing1" << std::endl;
  std::cout << "\ntestRing1()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      "COCc1c(ncc2[nH]c3cccc(Oc4ccc(Cl)cc4)c3c12)C(=O)OC(C)C",
      //      "COCc1cnc(C(=O)OC(C)C)c2[nH]c3cc(Oc4ccc(Cl)cc4)ccc3c12", //
      //      original molecule
      "COCc1cnc(C(=O)OC(C)C)c2[nH]ccc(Oc4ccc(Cl)cc4)cccc12",  // ring 3 removed
  };
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i)));  // with RING INFO
  }

  {
    MCSParameters p;
    p.BondCompareParameters.RingMatchesRingOnly = true;
    p.BondCompareParameters.CompleteRingsOnly = true;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
              << " atoms, " << res.NumBonds << " bonds\n";
    printTime();
    TEST_ASSERT(res.NumAtoms == 16 && res.NumBonds == 17);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRings = true;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS MatchFusedRings: " << res.SmartsString << " "
              << res.NumAtoms << " atoms, " << res.NumBonds << " bonds\n";
    printTime();
    TEST_ASSERT(res.NumAtoms == 16 && res.NumBonds == 17);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRingsStrict = true;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS MatchFusedRingsStrict: " << res.SmartsString << " "
              << res.NumAtoms << " atoms, " << res.NumBonds << " bonds\n";
    printTime();
    TEST_ASSERT(res.NumAtoms == 12 && res.NumBonds == 12);
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test504() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS test504" << std::endl;
  std::cout << "\ntest504()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      // TEST 504
      "C(CCNC(C1CC1[c:1]1[c:2]c(Cl)c(Cl)c[c:3]1)=O)CCN1CCC(NC(Nc2ccc(Cl)cc2)=O)"
      "C1 CHEMBL545864",  // - QUERY

      "FC(F)(F)c1cc(NC(N2CCCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)ccc1Cl "
      "CHEMBL528228",
      "FC(F)(F)c1cc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)C2)=O)ccc1Cl "
      "CHEMBL525875",
      "Fc1ccc(NC(N2CCCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)cc1C(F)(F)F "
      "CHEMBL527277",
      "FC(F)(F)c1cc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)ccc1Cl "
      "CHEMBL537333",
      "Fc1ccc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)C2)=O)cc1C(F)(F)F "
      "CHEMBL588077",
      "FC(F)(F)c1ccc(NC(NC2CCN(CCCCCNC(C3CC3c3cc(Cl)c(Cl)cc3)=O)C2)=O)cc1 "
      "CHEMBL525307",
      "Fc1ccc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)cc1C(F)(F)F "
      "CHEMBL581847",
      "FC(F)(F)c1ccc(NC(NC2CCN(CCCCCNC(C3CC3c3cc(Cl)c(Cl)cc3)=O)CC2)=O)cc1 "
      "CHEMBL579547",
      "N#Cc1cccc(NC(NC2CCN(CCCCCNC(C3CC3c3ccc(Cl)c(Cl)c3)=O)CC2)=O)c1 "
      "CHEMBL529994",
  };
  RWMol *qm = SmilesToMol(getSmilesOnly(smi[0]));
  unsigned int nq = qm->getNumAtoms();
  for (size_t ai = 0; ai < nq; ai++) {
    Atom *atom = qm->getAtomWithIdx(ai);
    atom->setProp("molAtomMapNumber", (int)ai);
  }
  std::cout << "Query +MAP " << MolToSmiles(*qm) << "\n";
  mols.emplace_back(qm);  // with RING INFO
  for (size_t i = 1; i < sizeof(smi) / sizeof(smi[0]); i++) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(smi[i])));  // with RING INFO
  }
  t0 = nanoClock();
  MCSResult res = findMCS(mols);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 34 && res.NumBonds == 36);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test18() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS test18" << std::endl;
  std::cout << "\ntest18()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      // TEST 18
      "Cc1nc(CN(C(C)c2ncccc2)CCCCN)ccc1 CHEMBL1682991",  //-- QUERY
      "Cc1ccc(CN(C(C)c2ccccn2)CCCCN)nc1 CHEMBL1682990",
      "Cc1cccnc1CN(C(C)c1ccccn1)CCCCN CHEMBL1682998",
      "CC(N(CCCCN)Cc1c(N)cccn1)c1ccccn1 CHEMBL1682987",
      "Cc1cc(C)c(CN(C(C)c2ccccn2)CCCCN)nc1 CHEMBL1682992",
      "Cc1cc(C(C)N(CCCCN)Cc2c(C)cccn2)ncc1 CHEMBL1682993",
      "Cc1nc(C(C)N(CCCCN)Cc2nc3c([nH]2)cccc3)ccc1 CHEMBL1682878",
      "CC(c1ncccc1)N(CCCCN)Cc1nc2c([nH]1)cccc2 CHEMBL1682867",
      "CC(N(CCCCN)Cc1c(C(C)(C)C)cccn1)c1ccccn1 CHEMBL1682989",
      "CC(N(CCCCN)Cc1c(C(F)(F)F)cccn1)c1ccccn1 CHEMBL1682988",
      // # 18 .  20 20 0.04 sec. Python MCS: CC(c1ccccn1)N(CCCCN)Ccnccc
  };
  RWMol *qm = SmilesToMol(getSmilesOnly(smi[0]));
  unsigned int nq = qm->getNumAtoms();
  for (size_t ai = 0; ai < nq; ai++) {
    Atom *atom = qm->getAtomWithIdx(ai);
    atom->setProp("molAtomMapNumber", (int)ai);
  }
  std::cout << "Query +MAP " << MolToSmiles(*qm) << "\n";
  mols.emplace_back(qm);  // with RING INFO
  for (size_t i = 1; i < sizeof(smi) / sizeof(smi[0]); i++) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(smi[i])));  // with RING INFO
  }
  t0 = nanoClock();
  MCSResult res = findMCS(mols);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 21 && res.NumBonds == 21);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testThreshold() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testThreshold" << std::endl;
  std::cout << "\ntestThreshold()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      "CCC", "CCCO", "CCCN", "CC",
      //        "CCC", "CC", //th=0.5
  };
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i)));
  }
  findMCS(mols);
  MCSParameters p;
  p.Threshold = 0.7;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  p.Threshold = 1.0;  // restore default value
  TEST_ASSERT(res.NumAtoms == 3 && res.NumBonds == 2);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testTarget_no_10188_30149() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testTarget_no_10188_30149" << std::endl;
  std::cout << "\ntestTarget_no_10188_30149()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      // Target_no_10188_30149.txt // VERY SLOWER than Python
      "CN(C)CCNC(=O)c1ccc(-c2n[nH]c3cc(Nc4ccccc4Cl)ccc32)cc1 CHEMBL399167",
      "O=C(O)c1cccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL197613",
      "c1ccc(Nc2ccc3c(c2)[nH]nc3-c2ccccc2)cc1 CHEMBL383177",  /// == QUERY
      "NC(=O)c1cccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL199136",
      "Clc1ccccc1Nc1ccc2c(c1)n[nH]c2-c1ccccc1 CHEMBL440566",
      "O=C(NCCCN1CCOCC1)c1cccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL198687",
      "O=C(O)c1ccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)cc1 CHEMBL197698",
      "O=C(NC1CCNCC1)c1cccc(-c2n[nH]c3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL194806",
      "COc1ccccc1Nc1ccc2c(c1)[nH]nc2-c1ccccc1 CHEMBL254443",
      "CN(C)CCNC(=O)c1cccc(-c2[nH]nc3cc(Nc4ccccc4Cl)ccc32)c1 CHEMBL198821",
  };
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i)));
  }
  MCSParameters p;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 15 && res.NumBonds == 14);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testTarget_no_10188_49064() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testTarget_no_10188_49064" << std::endl;
  std::cout << "\ntestTarget_no_10188_49064()\n";
  std::vector<ROMOL_SPTR> mols;
  // clang-format off
  const char* smi[] = {
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3cccc(O)c3)nc21",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(F)cc3)nc21",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(OC4OC(CO)C(O)C(O)C4O)cc3)nc21",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(NC(=O)CCl)cc3)nc21",
      "Cn1c2nc(Nc3ccc(NC(=O)CCN)cc3)ncc2cc(-c2c(Cl)cccc2Cl)c1=O",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3cccc(COCC(O)CO)c3)nc21",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(O)cc3)nc21",
      "CC(=O)Nc1ccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)cc1",
      "Cn1c2nc(Nc3ccc(N)cc3)ncc2cc(-c2c(Cl)cccc2Cl)c1=O",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(NC(=O)CCNC(=O)OC(C)(C)C)cc3)nc21",
      "Cc1ccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)cc1",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3cccc(CO)c3)nc21",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(NCC(O)CO)cc3)nc21",
      "CCc1cccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)c1",
      "Cn1c2nc(Nc3cccc(N)c3)ncc2cc(-c2c(Cl)cccc2Cl)c1=O",
      "CC(=O)Nc1cccc(Nc2ncc3cc(-c4c(Cl)cccc4Cl)c(=O)n(C)c3n2)c1",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(CCO)cc3)nc21",
      "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(I)cc3)nc21",
      "CN1CCN(C(=O)c2ccc(Nc3ncc4cc(-c5c(Cl)cccc5Cl)c(=O)n(C)c4n3)cc2)CC1",
  };
  // clang-format on
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i)));
  }
  MCSParameters p;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 15 && res.NumBonds == 14);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

#define MCSTESTREPEATS 0  // To run MCS repeatedly to measure performance
MCSResult checkMCS(const std::vector<ROMOL_SPTR> mols, const MCSParameters p,
                   unsigned int expectedAtoms, unsigned int expectedBonds) {
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  // std::shared_ptr<RWMol>
  std::unique_ptr<ROMol> mcsMol(SmartsToMol(res.SmartsString));
  std::unique_ptr<ROMol> cleanMCSMol(SmilesToMol(MolToSmiles(*mcsMol)));
  std::cout << "MCS: " << res.SmartsString << " " << MolToSmiles(*cleanMCSMol)
            << " " << res.NumAtoms << " atoms, " << res.NumBonds << " bonds"
            << std::endl;
#ifdef MCSTESTREPEATS
  for (int i = 0; i < MCSTESTREPEATS; i++) {
    res = findMCS(mols, &p);
  }
#endif
  printTime();
  if (res.NumAtoms != expectedAtoms || res.NumBonds != expectedBonds) {
    std::cerr << "testMaxDistance failed, expected " << expectedAtoms
              << " atoms, " << expectedBonds << " bonds"
              << " but got " << res.NumAtoms << " atoms and " << res.NumBonds
              << " bonds" << std::endl;
    TEST_ASSERT(res.NumAtoms == expectedAtoms && res.NumBonds == expectedBonds);
  }
  return res;
}

/* TODO: best practice on where to put a test data file into the repo? */
void testJnk1LigandsDistance() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testJnk1LigandsDistance" << std::endl;
  std::cout << "\ntestJnk1LigandsDistance()\n";
  std::vector<ROMOL_SPTR> mols;
  std::string rdbase = getenv("RDBASE");
  const char *jnk1sdf = "/Code/GraphMol/FMCS/testData/Jnk1_ligands.sdf";
  std::string fn(rdbase + jnk1sdf);

  std::unique_ptr<RDKit::MolSupplier> suppl;
  try {
    suppl.reset(new RDKit::SDMolSupplier(fn));
  } catch (...) {
    std::cerr << "ERROR: RDKit could not load input file" << std::endl;
    TEST_ASSERT(false);
  }
  ROMol *m1 = nullptr;
  ROMol *m2 = nullptr;
  while (!suppl->atEnd()) {
    ROMol *m = suppl->next();
    if (m) {
      if (m->getProp<std::string>("_Name") == "17124-1") {
        m1 = m;
      } else if (m->getProp<std::string>("_Name") == "18629-1") {
        m2 = m;
      } else {
        ROMOL_SPTR cleanupMol(m);  // don't leak memory
      }
    }
  }
  mols.emplace_back(m1);
  mols.emplace_back(m2);

  MCSParameters p;
  p.AtomTyper = MCSAtomCompareAnyHeavyAtom;
  p.BondTyper = MCSBondCompareOrderExact;
  p.AtomCompareParameters.MaxDistance = 3.0;
  MCSResult res = checkMCS(mols, p, 22, 23);

  SubstructMatchParameters smp;
  smp.useChirality = true;
  smp.uniquify = false;
  std::vector<MatchVectType> mvt1 = SubstructMatch(*m1, *(res.QueryMol), smp);
  std::vector<MatchVectType> mvt2 = SubstructMatch(*m2, *(res.QueryMol), smp);
  if (mvt1.size() != 2 || mvt2.size() != 2) {
    std::cerr << "jnk match atoms expected 2, 2: " << mvt1.size() << ","
              << mvt2.size() << std::endl;
    TEST_ASSERT(mvt1.size() == 2);
    TEST_ASSERT(mvt2.size() == 2);
  }

  std::list<int> forbidden1 = {18, 19, 25, 26};
  std::list<int> forbidden2 = {19};
  for (auto &matchVect : mvt1) {
    for (auto &matchPair : matchVect) {
      auto isPresent =
          std::find(forbidden1.begin(), forbidden1.end(), matchPair.second);
      if (isPresent != forbidden1.end()) {
        std::cerr << "mol1 index forbidden: " << matchPair.second << std::endl;
        TEST_ASSERT(isPresent == forbidden1.end());
      }
    }
  }
  for (auto &matchVect : mvt2) {
    for (auto &matchPair : matchVect) {
      auto isPresent =
          std::find(forbidden2.begin(), forbidden2.end(), matchPair.second);
      if (isPresent != forbidden2.end()) {
        std::cerr << "mol2 index forbidden: " << matchPair.second << std::endl;
        TEST_ASSERT(isPresent == forbidden2.end());
      }
    }
  }
  p.AtomCompareParameters.MaxDistance = -1.0;
  // Should match the flipped N if we don't filter on max distance
  checkMCS(mols, p, 23, 24);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testSegFault() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testSegFault" << std::endl;
  std::cout << "\ntestSegFault()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      "CN(CCCN(C)CCc1ccccc1)CCOC(c1ccccc1)c1ccccc1",
      "CN(CCCc1ccccc1)CCCN(C)CCOC(c1ccccc1)c1ccccc1",
      "Fc1ccc(C(OCCNCCCNCCc2ccccc2)c2ccc(F)cc2)cc1",
      "O=C(Cc1ccccc1)NCCCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "O=C(Cc1ccccc1)NCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "O=C(Cc1ccc(Br)cc1)NC=CNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "O=C(Cc1ccc(F)cc1)NCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "O=C(Cc1ccccc1)NCCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "CN(CCOC(c1ccc(F)cc1)c1ccc(F)cc1)CCN(C)CCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "COC(=O)C1C2CCC(CC1C(=O)Oc1ccccc1)N2C",
      "O=C1CN(CCc2ccccc2)CCN1CCOC(c1ccc(F)cc1)c1ccc(F)cc1",
      "CN(CCOC(c1ccccc1)c1ccccc1)CCN(C)CCc1ccc(F)cc1",
  };
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i)));
  }
  {
    MCSParameters p;
    p.BondCompareParameters.RingMatchesRingOnly = true;
    p.BondCompareParameters.CompleteRingsOnly = true;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
              << " atoms, " << res.NumBonds << " bonds\n";
    printTime();
    TEST_ASSERT(res.NumAtoms == 6 && res.NumBonds == 6);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRings = true;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS MatchFusedRings: " << res.SmartsString << " "
              << res.NumAtoms << " atoms, " << res.NumBonds << " bonds\n";
    printTime();
    TEST_ASSERT(res.NumAtoms == 6 && res.NumBonds == 6);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRingsStrict = true;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS MatchFusedRingsStrict: " << res.SmartsString << " "
              << res.NumAtoms << " atoms, " << res.NumBonds << " bonds\n";
    printTime();
    TEST_ASSERT(res.NumAtoms == 6 && res.NumBonds == 6);
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testAtomCompareIsotopes() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testAtomCompareIsotopes" << std::endl;
  std::cout << "\ntestAtomCompareIsotopes()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      "CC[13NH2]",
      "CC[13CH3]",
  };
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i)));
  }
  MCSParameters p;
  p.AtomTyper = MCSAtomCompareIsotopes;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 3 && res.NumBonds == 2);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testAtomCompareAnyAtom() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testAtomCompareAnyAtom" << std::endl;
  std::cout << "\ntestAtomCompareAnyAtom()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      "c1ccccc1C", "c1ccccc1O", "c1ccccc1Cl",
      "c1ccccc1F",  // opt
      "c1ccccc1N",  // opt
  };
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i)));
  }
  MCSParameters p;
  p.AtomTyper = MCSAtomCompareAny;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 7 && res.NumBonds == 7);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testAtomCompareAnyAtomBond() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testAtomCompareAnyAtomBond"
                       << std::endl;
  std::cout << "\ntestAtomCompareAnyAtom()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      "C1CCCCC1=C", "c1ccccc1O", "c1ccccc1Cl",
      "c1ccccc1F",  // opt
      "c1ccccc1N",  // opt
  };
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i)));
  }
  t0 = nanoClock();
  MCSParameters p;
  p.AtomTyper = MCSAtomCompareAny;
  p.BondTyper = MCSBondCompareAny;
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 7 && res.NumBonds == 7);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testAtomCompareAnyHeavyAtom() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testAtomCompareAnyAtom" << std::endl;
  std::cout << "\ntestAtomCompareAnyAtom()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      "[H]c1ccccc1C", "[H]c1ccccc1O",  // H matches H, O matches C
  };
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i), 0, false));
  }
  MCSParameters p;
  p.AtomTyper = MCSAtomCompareAnyHeavyAtom;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 8 && res.NumBonds == 8);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testAtomCompareAnyHeavyAtom1() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testAtomCompareAnyAtom" << std::endl;
  std::cout << "\ntestAtomCompareAnyAtom()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      "[H]c1ccccc1C", "Oc1ccccc1O",  // O matches C, H does not match O
  };
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i), 0, false));
  }
  MCSParameters p;
  p.AtomTyper = MCSAtomCompareAnyHeavyAtom;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 7 && res.NumBonds == 7);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testSimple() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testSimple" << std::endl;
  std::cout << "\ntestSimple()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      // LONG TIME TEST for performance analisis (about 30 sec)
      "CC(C)CC(NC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(C("
      "CO)NC(C(NC(c1ccncc1)=O)NC(=O)C(Cc1ccc(Cl)cc1)NC=O)=O)=O)C(NC("
      "CCCCNC(C)C)C(N1C(C(=O)NC(C)C(N)=O)CCC1)=O)=O CHEMBL439258 modified "
      "QUERY",  // CHEMBL439258
      "CC(C)CC(NC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(=O)C(Cc1ccccc1)NC(C(CO)NC(C("
      "NC(c1ccncc1)=O)NC(=O)C(Cc1ccc(Cl)cc1)NC(C(NC(C)=O)Cc1cc2ccccc2cc1)="
      "O)=O)=O)C(NC(CCCCNC(C)C)C(N1C(C(=O)NC(C)C(N)=O)CCC1)=O)=O "
      "CHEMBL439258",  // CHEMBL439258
      "CC(C)CC(NC(=O)CNC(=O)C(Cc1ccc(NC(C)=O)cc1)NC(C(CO)NC(C(NC(c1ccncc1)"
      "=O)NC(=O)C(Cc1ccc(Cl)cc1)NC(C(NC(C)=O)Cc1cc2ccccc2cc1)=O)=O)=O)C("
      "NC(CCCCNC(C)C)C(N1C(C(=O)NC(C)C(N)=O)CCC1)=O)=O CHEMBL439258 "
      "modified",  // CHEMBL439258

      "CCCCC(NC(C(CCC(O)=O)NC(C(CC(C)C)NC(C(C(C)C)NC(=O)C(CCC(O)=O)NC(C("
      "CCCN=C(N)N)NC(C(NC(=O)C(NC(C(NC(C1CCCNC(=O)CCC(N)C(=O)NC(CC(C)C)C(="
      "O)NC(C(C)O)C(=O)N1)=O)Cc1c[nH]cn1)=O)CC(C)C)CC(C)C)=O)=O)=O)=O)=O)"
      "C(NC(C)C(NC(CCCN=C(N)N)C(NC(C)C(NC(CCC(O)=O)C(NC(CCC(N)=O)C(NC(CC("
      "C)C)C(NC(C)C(NC(CCC(N)=O)C(NC(CCC(N)=O)C(NC(C)C(NC(Cc1c[nH]cn1)C("
      "NC(CO)C(NC(CC(N)=O)C(NC(CCCN=C(N)N)C(NC(CCCCN)C(NC(CC(C)C)C(NC("
      "CCCC)C(NC(C(NC(C(C)CC)C(NC(C(N)=O)C(C)CC)=O)=O)CCC(O)=O)=O)=O)=O)="
      "O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O CHEMBL438567",
      "CCC(C)C(NC(CNC(=O)C(C)NC(=O)C(C)NC(C(Cc1nc[nH]c1)NC(C(CC(N)=O)NC("
      "CNC(C(CO)NC(=O)C(C)NC(=O)C(CCC(N)=O)NC(C(NC(=O)C(NC(C(CCCN=C(N)N)"
      "NC(C(CCC(N)=O)NC(C(NC(C(CCCN=C(N)N)NC(CNC(C(CCC(N)=O)NC(C(CC(C)C)"
      "NC(C(C)N)=O)=O)=O)=O)=O)CC(C)C)=O)=O)=O)CC(C)C)CC(C)C)=O)=O)=O)=O)="
      "O)=O)C(NC(CC(C)C)C(NC(C(O)C)C(NC(CCSC)C(O)=O)=O)=O)=O CHEMBL429374",
      "CC(C)CC1NC(=O)C(CCCCN)NC(=O)C(Cc2ccc(O)cc2)NC(=O)CNC(=O)C2NC(=O)C("
      "NC(C(C(C)C)NC(CNC(C3NC(=O)CC3)=O)=O)=O)CSSCC(C(O)=O)NC(=O)C3N("
      "CCC3O)C(=O)C(Cc3ccccc3)NC(=O)C(CSSC2)NC1=O CHEMBL1076370",
  };
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i)));
  }
  {
    MCSParameters p;
    p.BondCompareParameters.RingMatchesRingOnly = true;
    p.BondCompareParameters.CompleteRingsOnly = true;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms
              << " atoms, " << res.NumBonds << " bonds\n";
    printTime();
    TEST_ASSERT(res.NumAtoms == 15 && res.NumBonds == 14);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRings = true;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS MatchFusedRings: " << res.SmartsString << " "
              << res.NumAtoms << " atoms, " << res.NumBonds << " bonds\n";
    printTime();
    TEST_ASSERT(res.NumAtoms == 15 && res.NumBonds == 14);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRingsStrict = true;
    t0 = nanoClock();
    MCSResult res = findMCS(mols, &p);
    std::cout << "MCS MatchFusedRingsStrict: " << res.SmartsString << " "
              << res.NumAtoms << " atoms, " << res.NumBonds << " bonds\n";
    printTime();
    TEST_ASSERT(res.NumAtoms == 15 && res.NumBonds == 14);
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testSimpleFast() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testSimpleFast" << std::endl;
  std::cout << "\ntestSimpleFast()\n";
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      // SHORT TEST for 26 bonds.
      // Python MCS = 26 bonds : COCc1cncc(c1):n:c1cccc(Oc2ccc(Cl)cc2)c1
      // MCS 26: COCc1c-ncc(c1)nc1cccc(c1)Oc1ccc(Cl)cc1 24 atoms, 26 bonds
      ///            "COCC1=C(N=CC2=C1C1=C(OC3=CC=C(Cl)C=C3)C=CC=C1N2)C(=O)OC(C)C",
      ///            "COCC1=CN=C(C(=O)OC(C)C)C2=C1C1=CC=C(OC3=CC=C(Cl)C=C3)C=C1N2",
      // The SAME, but pre AROMATIZATED (else PRECONDITION Exception with
      // Implicit Hs / 16 bonds only)
      "COCc1c(ncc2[nH]c3cccc(Oc4ccc(Cl)cc4)c3c12)C(=O)OC(C)C",
      "COCc1cnc(C(=O)OC(C)C)c2[nH]c3cc(Oc4ccc(Cl)cc4)ccc3c12",
  };
  for (auto &i : smi) {
    mols.emplace_back(SmilesToMol(getSmilesOnly(i)));
  }
  MCSParameters p;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 24 && res.NumBonds == 26);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void compareChirality(const char *target, const char *query,
                      bool useChirality) {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  std::cout << "test target = " << target << " query = " << query
            << " useChirality = " << useChirality << std::endl;
  ROMOL_SPTR target_ptr(SmilesToMol(target));
  ROMOL_SPTR query_ptr(SmartsToMol(query));

  std::vector<std::pair<int, int>> vect;
  bool sub_res = SubstructMatch(*target_ptr.get(), *query_ptr.get(), vect, true,
                                useChirality);

  std::cout << "res sub = " << sub_res << " size=" << vect.size() << std::endl;

  std::vector<ROMOL_SPTR> mols;
  mols.push_back(target_ptr);
  mols.push_back(query_ptr);

  MCSParameters p;
  p.AtomCompareParameters.MatchChiralTag = useChirality;
  p.BondCompareParameters.MatchStereo = useChirality;

  MCSResult mcs_res = findMCS(mols, &p);
  std::cout << "MCS: " << mcs_res.SmartsString << " " << mcs_res.NumAtoms
            << " atoms, " << mcs_res.NumBonds << " bonds\n";

  bool mcs_resb = query_ptr->getNumAtoms() == mcs_res.NumAtoms &&
                  query_ptr->getNumBonds() == mcs_res.NumBonds;
  TEST_ASSERT(sub_res == mcs_resb);
  if (sub_res != mcs_resb) {  // || vect.size() != mcs_res.NumAtoms) {
    BOOST_LOG(rdInfoLog) << "mcs_resb=" << mcs_resb
                         << "\t*** TEST FAILED ***\n";  // exit(1);
  }
}

void testSubMcsChirality(const char *target, const char *query) {
  compareChirality(target, query, false);
  compareChirality(target, query, true);
}

void testChirality() {
  BOOST_LOG(rdInfoLog) << "\n-------------------------------------"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testChirality" << std::endl;

  testSubMcsChirality("O[C@H](F)CCl", "C(F)C");  // MCS = CCF
  testSubMcsChirality("CC[C@H](F)Cl", "CCC");    // MCS = CCC
  testSubMcsChirality("CCC(F)Cl", "CC");         // MCS = CC

  // GREG's tests:
  testSubMcsChirality("C[C@H](F)CCl", "C[C@H](F)C");
  testSubMcsChirality("C[C@H](F)CCl", "C[C@@H](F)C");
  testSubMcsChirality(
      "O[C@H](F)CCl",
      "CC(F)C");  // DISMATCH. but actually it has non chiral MCS = CCF
  testSubMcsChirality("O[C@H](F)CCl", "OC(F)C");
  testSubMcsChirality("O[C@H](F)CCl", "O[C@H](F)C");
  testSubMcsChirality("O[C@H](F)CCl", "O[C@@H](F)C");

  testSubMcsChirality("OC(F)CCl", "OC(F)C");
  testSubMcsChirality("OC(F)CCl", "O[C@H](F)C");
  testSubMcsChirality("OC(F)CCl", "O[C@@H](F)C");

  testSubMcsChirality(
      "O[C@H](F)CCl",
      "O[C@H]C");  // Is OCC = MCS ? Degree of [C@H] is different
  testSubMcsChirality(
      "O[C@H](F)CCl",
      "O[C@@H]C");  // Is OCC = MCS ? Degree of [C@H] is different

  // ADD-IN TESTS:
  std::cout
      << "\n0. <<<<<<<<< actual MCS: [#6]-[#6]-[#6] 3 atoms, 2 bonds >>>>>>>\n";
  testSubMcsChirality(
      "CC[C@H](F)Cl",
      "CC[C@H]");  //  actual MCS is CCC always. We lost last [C@H]

  /* //TEMP (asymmetric implementation of match algorithm, that is incorrect for
  FMCS task):
     std::cout << "\nTEMP query/targ==MCS[asymmetric implementation of match
  algorithm]\n";
     std::cout << "1. <<<<<<<<< MCS: [#8]-[#6]-[#9] 3 atoms, 2 bonds >>>>>>>\n";
     testSubMcsChirality ("O[C@H](F)", "OC(F)"    ); // PASSED
  / /   testSubMcsChirality ("OC(F)"    , "O[C@H](F)"); // FAILED !!!
     std::cout << "\n2.1<<<<<<<<< MCS: [#8]-[#6]-[#6] 3 atoms, 2 bonds
  >>>>>>>\n";
     testSubMcsChirality ("OC[C@H](F)C", "OC(F)C"    ); // PASSED but incorrect
     std::cout << "\n2.2<<<<<<<<< MCS: [#8]-[#6](-[#9])-[#6] 4 atoms, 3 bonds
  >>>>>>>\n";
     testSubMcsChirality ("OC(F)C"     , "O[C@H](F)C"); // FAILED !!!
  */
  std::cout << "\tdone" << std::endl;
}

void testJSONParameters() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing FMCS testJSONParameters" << std::endl;
  MCSParameters pj;

  parseMCSParametersJSON(nullptr, &pj);
  parseMCSParametersJSON("{}", nullptr);

  pj = MCSParameters();  // parsing of empty string keeps default values
  parseMCSParametersJSON("", &pj);
  TEST_ASSERT(pj.MaximizeBonds && pj.Threshold == 1.0 && pj.Timeout == 0 &&
              !pj.AtomCompareParameters.MatchValences &&
              !pj.AtomCompareParameters.MatchChiralTag &&
              !pj.BondCompareParameters.MatchStereo &&
              !pj.BondCompareParameters.RingMatchesRingOnly &&
              !pj.BondCompareParameters.CompleteRingsOnly &&
              !pj.BondCompareParameters.MatchFusedRings &&
              !pj.BondCompareParameters.MatchFusedRingsStrict);

  {
    pj = MCSParameters();
    const char json[] =
        "{\"MaximizeBonds\": false, \"Threshold\": 0.7, \"Timeout\": 3"
        ", \"MatchValences\": true, \"MatchChiralTag\": true"
        ", \"MatchStereo\": true, \"RingMatchesRingOnly\": true"
        ", \"CompleteRingsOnly\": true"
        ", \"MatchFusedRings\": true"
        ", \"MatchFusedRingsStrict\": true"
        ", \"InitialSeed\": \"CNC\""
        ", \"AtomCompare\": \"Isotopes\""
        ", \"BondCompare\": \"OrderExact\""
        "}";
    parseMCSParametersJSON(json, &pj);
    TEST_ASSERT(!pj.MaximizeBonds && pj.Threshold == 0.7 && pj.Timeout == 3 &&
                pj.AtomCompareParameters.MatchValences &&
                pj.AtomCompareParameters.MatchChiralTag &&
                pj.AtomCompareParameters.RingMatchesRingOnly &&
                !pj.AtomCompareParameters.CompleteRingsOnly &&
                pj.BondCompareParameters.MatchStereo &&
                pj.BondCompareParameters.RingMatchesRingOnly &&
                pj.BondCompareParameters.CompleteRingsOnly &&
                pj.BondCompareParameters.MatchFusedRings &&
                pj.BondCompareParameters.MatchFusedRingsStrict &&
                pj.AtomTyper == MCSAtomCompareIsotopes &&
                pj.BondTyper == MCSBondCompareOrderExact &&
                0 == strcmp(pj.InitialSeed.c_str(), "CNC"));
  }
  {
    // Atom* and Bond* versions override simple
    // RingMatchesRingOnly and CompleteRingsOnly
    pj = MCSParameters();
    const char json[] =
        "{\"MaximizeBonds\": false, \"Threshold\": 0.7, \"Timeout\": 3"
        ", \"MatchValences\": true, \"MatchChiralTag\": true"
        ", \"MatchStereo\": true, \"RingMatchesRingOnly\": false"
        ", \"AtomRingMatchesRingOnly\": true"
        ", \"BondRingMatchesRingOnly\": true"
        ", \"CompleteRingsOnly\": false"
        ", \"AtomCompleteRingsOnly\": true"
        ", \"BondCompleteRingsOnly\": true"
        ", \"MatchFusedRings\": true"
        ", \"MatchFusedRingsStrict\": true"
        ", \"InitialSeed\": \"CNC\""
        "}";
    parseMCSParametersJSON(json, &pj);
    TEST_ASSERT(!pj.MaximizeBonds && pj.Threshold == 0.7 && pj.Timeout == 3 &&
                pj.AtomCompareParameters.MatchValences &&
                pj.AtomCompareParameters.MatchChiralTag &&
                pj.AtomCompareParameters.RingMatchesRingOnly &&
                pj.AtomCompareParameters.CompleteRingsOnly &&
                pj.BondCompareParameters.MatchStereo &&
                pj.BondCompareParameters.RingMatchesRingOnly &&
                pj.BondCompareParameters.CompleteRingsOnly &&
                pj.BondCompareParameters.MatchFusedRings &&
                pj.BondCompareParameters.MatchFusedRingsStrict &&
                0 == strcmp(pj.InitialSeed.c_str(), "CNC"));
  }

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithubIssue481() {
  BOOST_LOG(rdInfoLog) << "\n-------------------------------------"
                       << std::endl;
  BOOST_LOG(rdInfoLog)
      << "Testing github #481 : order dependence in FMCS with chirality"
      << std::endl;

  {
    std::string s1 = "O[C@H](C)Cl";
    std::string s2 = "O[C@@H](C)Cl";
    ROMOL_SPTR ptr1(SmilesToMol(s1));
    ROMOL_SPTR ptr2(SmilesToMol(s2));
    std::vector<ROMOL_SPTR> mols;
    mols.push_back(ptr1);
    mols.push_back(ptr2);

    BOOST_LOG(rdInfoLog) << "**** mols:" << s1 << "   " << s2 << "\n";
    // #688:
    std::vector<std::pair<int, int>> vect;
    ROMOL_SPTR mcs_mol(SmilesToMol("OCCl"));
    bool sub_res = SubstructMatch(*ptr2, *mcs_mol, vect, true, false);
    std::cout << "query=OCCl "
              << "SubstructMatch(useChirality=false) res =" << sub_res
              << " size=" << vect.size() << std::endl;
    sub_res = SubstructMatch(*ptr2, *mcs_mol, vect, true, true);
    std::cout << "query=OCCl "
              << "SubstructMatch(useChirality=true ) res =" << sub_res
              << " size=" << vect.size() << std::endl;

    {
      MCSParameters p;
      p.AtomCompareParameters.MatchChiralTag = false;
      p.BondCompareParameters.MatchStereo = false;
      MCSResult mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      TEST_ASSERT(mcs_res.NumAtoms == 4);
      TEST_ASSERT(mcs_res.NumBonds == 3);

      p.AtomCompareParameters.MatchChiralTag = true;
      p.BondCompareParameters.MatchStereo = true;
      mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      TEST_ASSERT(mcs_res.NumAtoms == 3);  // #688 was 2);
      TEST_ASSERT(mcs_res.NumBonds == 2);  // #688 was 1);
    }

    BOOST_LOG(rdInfoLog) << "------ REVERSE mols -------- \n";

    mols.clear();
    mols.push_back(ptr2);
    mols.push_back(ptr1);
    {
      MCSParameters p;
      p.AtomCompareParameters.MatchChiralTag = false;
      p.BondCompareParameters.MatchStereo = false;
      MCSResult mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      TEST_ASSERT(mcs_res.NumAtoms == 4);
      TEST_ASSERT(mcs_res.NumBonds == 3);

      p.AtomCompareParameters.MatchChiralTag = true;
      p.BondCompareParameters.MatchStereo = true;
      mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      TEST_ASSERT(mcs_res.NumAtoms == 3);  // #688 was 2);
      TEST_ASSERT(mcs_res.NumBonds == 2);  // #688 was 1);
    }
  }

  {
    std::string s1 = "OC(C)Cl";
    std::string s2 = "O[C@H](C)Cl";
    ROMOL_SPTR ptr1(SmilesToMol(s1));
    ROMOL_SPTR ptr2(SmilesToMol(s2));
    std::vector<ROMOL_SPTR> mols;
    mols.push_back(ptr1);
    mols.push_back(ptr2);

    BOOST_LOG(rdInfoLog) << "**** mols:" << s1 << "   " << s2 << "\n";

    {
      MCSParameters p;
      p.AtomCompareParameters.MatchChiralTag = false;
      p.BondCompareParameters.MatchStereo = false;
      MCSResult mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      TEST_ASSERT(mcs_res.NumAtoms == 4);
      TEST_ASSERT(mcs_res.NumBonds == 3);

      p.AtomCompareParameters.MatchChiralTag = true;
      p.BondCompareParameters.MatchStereo = true;
      mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      std::vector<std::pair<int, int>> vect;
      bool sub_res =
          SubstructMatch(*mols[1].get(), *mols[0].get(), vect, true, true);
      if (!sub_res) {  // actually == true & 4, 3 !!!
        TEST_ASSERT(mcs_res.NumAtoms == 0);
        TEST_ASSERT(mcs_res.NumBonds == 0);
      }
    }

    BOOST_LOG(rdInfoLog) << "------ REVERSE mols -------- \n";

    mols.clear();
    mols.push_back(ptr2);
    mols.push_back(ptr1);
    {
      MCSParameters p;
      p.AtomCompareParameters.MatchChiralTag = false;
      p.BondCompareParameters.MatchStereo = false;
      MCSResult mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      TEST_ASSERT(mcs_res.NumAtoms == 4);
      TEST_ASSERT(mcs_res.NumBonds == 3);

      p.AtomCompareParameters.MatchChiralTag = true;
      p.BondCompareParameters.MatchStereo = true;
      mcs_res = findMCS(mols, &p);
      BOOST_LOG(rdInfoLog) << "MCS: " << mcs_res.SmartsString << " "
                           << mcs_res.NumAtoms << " atoms, " << mcs_res.NumBonds
                           << " bonds\n";
      std::vector<std::pair<int, int>> vect;
      bool sub_res =
          SubstructMatch(*mols[1].get(), *mols[0].get(), vect, true, true);
      if (!sub_res) {
        TEST_ASSERT(mcs_res.NumAtoms == 1);
        TEST_ASSERT(mcs_res.NumBonds == 0);
        TEST_ASSERT(mcs_res.SmartsString == "[#17]");
      }
    }
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testInitialSeed() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "FMCS testInitialSeed()" << std::endl;

  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      "CC1CCC(N)CC1", "CC1CC(C)CC(C)C1",  // OK test.sdf
  };

  for (auto &i : smi) {
    std::string id;
    mols.emplace_back(SmilesToMol(getSmilesOnly(i, &id)));
  }
  MCSParameters p;
  p.InitialSeed = "CC";
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms == 7 && res.NumBonds == 7);
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testInitialSeed2() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "FMCS testInitialSeed2()" << std::endl;

  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      "Cc1c(F)c(N2CCNC(C)C2)cc2c1c(=O)c(C(=O)O)cn2C1CC1",
      "COc1c(N2CCNC(C)C2)c(F)cc2c(=O)c(C(=O)O)cn(C3CC3)c12",
  };
  const char *initial_smarts = "CCNCCNcccccccnC1CC1";
  BOOST_LOG(rdInfoLog) << "initial_smarts: " << initial_smarts << std::endl;

  for (auto &i : smi) {
    std::string id;
    mols.emplace_back(SmilesToMol(getSmilesOnly(i, &id)));
    std::unique_ptr<ROMol> seed(SmartsToMol(initial_smarts));
    MatchVectType match;
    bool matched = SubstructMatch(*mols.back(), *seed, match);
    BOOST_LOG(rdInfoLog) << (matched ? "RDKit MATCHED " : "RDKit DISmatched ")
                         << i << std::endl;
  }
  MCSParameters p;
  p.InitialSeed = initial_smarts;
  t0 = nanoClock();
  MCSResult res = findMCS(mols, &p);
  std::cout << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms != 0);

  // Make Initial Seed from MCS
  p.Verbose = true;
  p.InitialSeed =
      "[#6]1-[#6]-[#7]-[#6](-[#6]-[#7]-1-[#6]1:[#6](:[#6]:[#6]2:[#6](:[#6](:[#"
      "6]:[#7](-[#6]3-[#6]-[#6]-3):[#6]:2:[#6]:1)-[#6](=[#8])-[#8])=[#8])-[#9])"
      "-[#6]";  // 25 atoms, 28 bonds
  BOOST_LOG(rdInfoLog)
      << "\n\nFound MCS as the only initial seed (25 atoms, 28 bonds): \n"
      << p.InitialSeed << std::endl;
  t0 = nanoClock();
  res = findMCS(mols, &p);
  BOOST_LOG(rdInfoLog) << "MCS: " << res.SmartsString << " " << res.NumAtoms
                       << " atoms, " << res.NumBonds << " bonds\n";
  printTime();
  TEST_ASSERT(res.NumAtoms != 0);

  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub631() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing github issue 631: FindMCS "
                          "matchChiralTag=True does not match self"
                       << std::endl;
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {
      // all examples derived from the bug report
      "CN(C)[C@@H]1CCCNC1",  // == MCS
      "C1C=CCN1[C@@H]1CCCNC1",
      "Cc1cc2c(cc1C)C(=O)N([C@@H]1CCC(=O)NC1=O)C2=O",
  };

  for (int pass = 0; pass < 2; ++pass, mols.clear()) {
    for (auto &i : smi) {
      RWMol *m = SmilesToMol(getSmilesOnly(i));
      TEST_ASSERT(m);

      if (0 == pass) {
        mols.clear();  // use a pair of the same molecules only. On the second
      }
      // pass use all.

      mols.emplace_back(m);
      mols.push_back(ROMOL_SPTR(new ROMol(*m)));
      {
        MCSParameters p;
        p.AtomCompareParameters.MatchChiralTag = false;
        //          p.Verbose = true;
        MCSResult res = findMCS(mols, &p);
        BOOST_LOG(rdInfoLog)
            << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n"
            << std::endl;
        ;

        TEST_ASSERT(res.NumAtoms == mols[0]->getNumAtoms());
        TEST_ASSERT(res.NumBonds == mols[0]->getNumBonds());
      }

      {
        MCSParameters p;
        p.AtomCompareParameters.MatchChiralTag = true;
        // p.BondCompareParameters.MatchStereo = useChirality;
        //        p.Verbose = true;
        MCSResult res = findMCS(mols, &p);
        BOOST_LOG(rdInfoLog)
            << "MCS: " << res.SmartsString << " " << res.NumAtoms << " atoms, "
            << res.NumBonds << " bonds\n"
            << std::endl;
        ;

        TEST_ASSERT(res.NumAtoms == mols[0]->getNumAtoms());
        TEST_ASSERT(res.NumBonds == mols[0]->getNumBonds());
      }
      BOOST_LOG(rdInfoLog) << "============================================"
                           << std::endl;
    }
  }
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testFormalChargeMatch() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing including formal charge in matches "
                       << std::endl;
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {"CCNC", "CCN(C)C", "CC[N+](C)C"};

  for (auto &i : smi) {
    RWMol *m = SmilesToMol(getSmilesOnly(i));
    TEST_ASSERT(m);

    mols.emplace_back(m);
  }
  {
    // by default charge is not used.
    MCSParameters p;
    MCSResult res = findMCS(mols, &p);
    BOOST_LOG(rdInfoLog) << "MCS: " << res.SmartsString << " " << res.NumAtoms
                         << " atoms, " << res.NumBonds << " bonds\n"
                         << std::endl;
    ;

    TEST_ASSERT(res.NumAtoms == 4);
    TEST_ASSERT(res.NumBonds == 3);
  }
  {
    // check charge
    MCSParameters p;
    p.AtomCompareParameters.MatchFormalCharge = true;
    MCSResult res = findMCS(mols, &p);
    BOOST_LOG(rdInfoLog) << "MCS: " << res.SmartsString << " " << res.NumAtoms
                         << " atoms, " << res.NumBonds << " bonds\n"
                         << std::endl;
    ;

    TEST_ASSERT(res.NumAtoms == 2);
    TEST_ASSERT(res.NumBonds == 1);
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub2034() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github #2034: add test for ring--non-ring "
                          "matches at the atom level"
                       << std::endl;
  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {"C1CC1N2CC2", "C1CC1N"};

  for (auto &i : smi) {
    auto m = SmilesToMol(getSmilesOnly(i));
    TEST_ASSERT(m);

    mols.emplace_back(m);
  }
  {
    // by default we're not doing ringMatchesRingOnly.
    MCSResult res = findMCS(mols);
    BOOST_LOG(rdInfoLog) << "MCS: " << res.SmartsString << " " << res.NumAtoms
                         << " atoms, " << res.NumBonds << " bonds\n"
                         << std::endl;

    TEST_ASSERT(res.NumAtoms == 4);
    TEST_ASSERT(res.NumBonds == 4);
  }

  {
    MCSParameters p;
    p.AtomCompareParameters.RingMatchesRingOnly = true;
    //          p.Verbose = true;
    MCSResult res = findMCS(mols, &p);
    BOOST_LOG(rdInfoLog) << "MCS: " << res.SmartsString << " " << res.NumAtoms
                         << " atoms, " << res.NumBonds << " bonds\n"
                         << std::endl;

    TEST_ASSERT(res.NumAtoms == 3);
    TEST_ASSERT(res.NumBonds == 3);
  }
  {
    // set it:
    bool maximizeBonds = true, verbose = false, matchValences = false,
         ringMatchesRingOnly = true;
    double threshold = 1.0;
    unsigned int timeout = 3000;
    MCSResult res = findMCS(mols, maximizeBonds, threshold, timeout, verbose,
                            matchValences, ringMatchesRingOnly);
    BOOST_LOG(rdInfoLog) << "MCS: " << res.SmartsString << " " << res.NumAtoms
                         << " atoms, " << res.NumBonds << " bonds\n"
                         << std::endl;

    TEST_ASSERT(res.NumAtoms == 3);
    TEST_ASSERT(res.NumBonds == 3);
  }
  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub945() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github #945: MCS returning partial rings "
                          "with completeRingsOnly=True"
                       << std::endl;
  {
    std::vector<ROMOL_SPTR> mols;
    const char *smi[] = {"c1cc2ccccc2s1", "c1cc2ccccc2o1"};

    for (auto &i : smi) {
      auto m = SmilesToMol(getSmilesOnly(i));
      TEST_ASSERT(m);

      mols.emplace_back(m);
    }
    {
      MCSParameters p;
      // p.Verbose = true;
      p.BondCompareParameters.CompleteRingsOnly = true;
      MCSResult res = findMCS(mols, &p);
      // std::cerr << "MCS CompleteRingsOnly: "
      //           << res.SmartsString << " " << res.NumAtoms
      //           << " atoms, " << res.NumBonds << " bonds\n"
      //           << std::endl;

      TEST_ASSERT(res.NumAtoms == 6);
      TEST_ASSERT(res.NumBonds == 6);
    }
    {
      MCSParameters p;
      p.BondCompareParameters.CompleteRingsOnly = true;
      p.BondCompareParameters.MatchFusedRings = true;
      MCSResult res = findMCS(mols, &p);
      // std::cerr << "MCS MatchFusedRings: "
      //           << res.SmartsString << " " << res.NumAtoms
      //           << " atoms, " << res.NumBonds << " bonds\n"
      //           << std::endl;
      TEST_ASSERT(res.NumAtoms == 6);
      TEST_ASSERT(res.NumBonds == 6);
    }
    {
      MCSParameters p;
      p.BondCompareParameters.CompleteRingsOnly = true;
      p.BondCompareParameters.MatchFusedRingsStrict = true;
      MCSResult res = findMCS(mols, &p);
      // std::cerr << "MCS MatchFusedRingsStrict: "
      //           << res.SmartsString << " " << res.NumAtoms
      //           << " atoms, " << res.NumBonds << " bonds\n"
      //           << std::endl;
      TEST_ASSERT(res.NumAtoms == 6);
      TEST_ASSERT(res.NumBonds == 6);
    }
  }
  {
    std::vector<ROMOL_SPTR> mols;
    const char *smi[] = {"c1cc2ccc(C)cc2s1", "c1cc2c(cccc2s1)C"};

    for (auto &i : smi) {
      auto m = SmilesToMol(getSmilesOnly(i));
      TEST_ASSERT(m);

      mols.emplace_back(m);
    }
    {
      MCSParameters p;
      // p.Verbose = true;
      p.BondCompareParameters.CompleteRingsOnly = true;
      MCSResult res = findMCS(mols, &p);
      // std::cerr << "MCS CompleteRingsOnly: "
      //           << res.SmartsString << " " << res.NumAtoms
      //           << " atoms, " << res.NumBonds << " bonds\n"
      //           << std::endl;
      TEST_ASSERT(res.NumAtoms == 9);
      TEST_ASSERT(res.NumBonds == 10);
    }
    {
      MCSParameters p;
      p.BondCompareParameters.CompleteRingsOnly = true;
      p.BondCompareParameters.MatchFusedRings = true;
      MCSResult res = findMCS(mols, &p);
      // std::cerr << "MCS MatchFusedRings: "
      //           << res.SmartsString << " " << res.NumAtoms
      //           << " atoms, " << res.NumBonds << " bonds\n"
      //           << std::endl;
      TEST_ASSERT(res.NumAtoms == 9);
      TEST_ASSERT(res.NumBonds == 10);
    }
    {
      MCSParameters p;
      p.BondCompareParameters.CompleteRingsOnly = true;
      p.BondCompareParameters.MatchFusedRingsStrict = true;
      MCSResult res = findMCS(mols, &p);
      // std::cerr << "MCS MatchFusedRingsStrict: "
      //           << res.SmartsString << " " << res.NumAtoms
      //           << " atoms, " << res.NumBonds << " bonds\n"
      //           << std::endl;
      TEST_ASSERT(res.NumAtoms == 9);
      TEST_ASSERT(res.NumBonds == 10);
    }

    // if we don't require rings to be closed, then we get the solution that
    // has more atoms:
    {
      MCSParameters p;
      p.BondCompareParameters.CompleteRingsOnly = false;
      MCSResult res = findMCS(mols, &p);
      // std::cerr << "MCS: " << res.SmartsString << " " << res.NumAtoms
      //           << " atoms, " << res.NumBonds << " bonds\n"
      //           << std::endl;
      TEST_ASSERT(res.NumAtoms == 10);
      TEST_ASSERT(res.NumBonds == 10);
    }
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub2420() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github #2420: MCS code doesn't return "
                          "envelope MCS if CompleteRingsOnly is true"
                       << std::endl;

  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {"C1CC2C(CC1)CCCC2", "C1CCCCCCCCC1"};

  for (auto &i : smi) {
    auto m = SmilesToMol(getSmilesOnly(i));
    TEST_ASSERT(m);

    mols.emplace_back(m);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS CompleteRingsOnly: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 10);
    TEST_ASSERT(res.NumBonds == 10);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRings = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS MatchFusedRings: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 10);
    TEST_ASSERT(res.NumBonds == 10);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRingsStrict = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS MatchFusedRingsStrict: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 0);
    TEST_ASSERT(res.NumBonds == 0);
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub2663() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github #2663: The C++ MCS code "
                          "generates ambiguous SMARTS strings"
                       << std::endl;

  {
    std::vector<ROMOL_SPTR> mols;
    const char *smi[] = {"C1C(C)CC2CC12", "CC1CCCC2CCCC12"};

    for (auto &i : smi) {
      auto m = SmilesToMol(getSmilesOnly(i));
      TEST_ASSERT(m);

      mols.emplace_back(m);
    }
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    MCSResult res = findMCS(mols, &p);
    TEST_ASSERT(res.NumAtoms == 7);
    TEST_ASSERT(res.NumBonds == 7);
    TEST_ASSERT(res.SmartsString ==
                "[#6]1-&@[#6](-&!@[#6])-&@[#6]-&@[#6]-&@[#6]-&@[#6]-&@1");
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub2662() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github #2662: C++ MCS code returns a null "
                          "MCS between methylcyclopentane and methylcyclohexane"
                       << std::endl;

  {
    std::vector<ROMOL_SPTR> mols;
    const char *smi[] = {"CC1CCCC1", "CC1CCCCC1"};

    for (auto &i : smi) {
      auto m = SmilesToMol(getSmilesOnly(i));
      TEST_ASSERT(m);

      mols.emplace_back(m);
    }
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    MCSResult res = findMCS(mols, &p);
    TEST_ASSERT(res.NumAtoms == 2);
    TEST_ASSERT(res.NumBonds == 1);
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testNaphthalenes() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing naphthalenes with different substitution "
                          "patterns"
                       << std::endl;

  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {"Cc1cccc2ccccc12", "Cc1ccc2ccccc2c1"};

  for (auto &i : smi) {
    auto m = SmilesToMol(getSmilesOnly(i));
    TEST_ASSERT(m);

    mols.emplace_back(m);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS CompleteRingsOnly: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 11);
    TEST_ASSERT(res.NumBonds == 11);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRings = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS MatchFusedRings: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 10);
    TEST_ASSERT(res.NumBonds == 11);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRingsStrict = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS MatchFusedRingsStrict: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 10);
    TEST_ASSERT(res.NumBonds == 11);
  }
  // add cyclodecapentaene: then we'll get the envelope
  // unless we specify MatchFusedRingsStrict
  char cyclodecapentaene[] = "Cc1ccccccccc1";
  {
    auto m = SmilesToMol(cyclodecapentaene);
    TEST_ASSERT(m);
    mols.emplace_back(m);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS CompleteRingsOnly: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 11);
    TEST_ASSERT(res.NumBonds == 11);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRings = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS MatchFusedRings: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 11);
    TEST_ASSERT(res.NumBonds == 11);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRingsStrict = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS MatchFusedRingsStrict: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 2);
    TEST_ASSERT(res.NumBonds == 1);
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testBicycles() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing bicycles with different substitution "
                          "patterns"
                       << std::endl;

  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {"CC1CCCC2CCCC12", "CC1CCC2CC12", "C1C(C)CC2CC12"};

  for (auto &i : smi) {
    auto m = SmilesToMol(getSmilesOnly(i));
    TEST_ASSERT(m);

    mols.emplace_back(m);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS CompleteRingsOnly: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 7);
    TEST_ASSERT(res.NumBonds == 7);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRings = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS MatchFusedRings: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 5);
    TEST_ASSERT(res.NumBonds == 5);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRingsStrict = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS MatchFusedRingsStrict: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 5);
    TEST_ASSERT(res.NumBonds == 5);
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testBicyclesTricycles() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing bicycles and tricycles with "
                          "different substitution patterns"
                       << std::endl;

  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {"CC1CCCC2CCCC12", "C1C(C)CCC2CCCC12",
                       "CC1CCCC(C3)2CCCC123", "C1C(C)CCC(C3)2CCCC123"};

  for (auto &i : smi) {
    auto m = SmilesToMol(getSmilesOnly(i));
    TEST_ASSERT(m);

    mols.emplace_back(m);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS CompleteRingsOnly: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 10);
    TEST_ASSERT(res.NumBonds == 10);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRings = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS MatchFusedRings: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 9);
    TEST_ASSERT(res.NumBonds == 10);
  }
  {
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.MatchFusedRingsStrict = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS MatchFusedRingsStrict: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 9);
    TEST_ASSERT(res.NumBonds == 10);
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void test_p38() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing p38 ligands" << std::endl;

  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {"C1COCCC1Nc1ncc2cc(Cc3c(F)cccc3)c(=O)n(C)c2n1",
                       "O(c1cc2c(n(C)c1=O)nc(NC1CCOCC1)nc2)c1ccccc1"};

  for (auto &i : smi) {
    auto m = SmilesToMol(getSmilesOnly(i));
    TEST_ASSERT(m);

    mols.emplace_back(m);
  }
  {
    MCSResult res = findMCS(mols, true);
    // std::cerr << "MCS CompleteRingsOnly: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 19);
    TEST_ASSERT(res.NumBonds == 21);
  }
  {
    MCSResult res =
        findMCS(mols, true, 1.0, 3600, false, false, true, true, false,
                AtomCompareElements, BondCompareOrder, PermissiveRingFusion);
    // std::cerr << "MCS MatchFusedRings: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 19);
    TEST_ASSERT(res.NumBonds == 21);
  }
  {
    MCSResult res =
        findMCS(mols, true, 1.0, 3600, false, false, true, true, false,
                AtomCompareElements, BondCompareOrder, StrictRingFusion);
    // std::cerr << "MCS MatchFusedRings: "
    //           << res.SmartsString << " " << res.NumAtoms
    //           << " atoms, " << res.NumBonds << " bonds\n"
    //           << std::endl;
    TEST_ASSERT(res.NumAtoms == 19);
    TEST_ASSERT(res.NumBonds == 21);
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGithub2714() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github #2714" << std::endl;

  {
    std::vector<ROMOL_SPTR> mols;
    const char *smi[] = {"CC1CCC1", "CCC1CC1"};

    for (auto &i : smi) {
      auto m = SmilesToMol(getSmilesOnly(i));
      TEST_ASSERT(m);

      mols.emplace_back(m);
    }
    MCSParameters p;
    p.BondCompareParameters.CompleteRingsOnly = true;
    p.BondCompareParameters.RingMatchesRingOnly = true;
    p.AtomCompareParameters.RingMatchesRingOnly = true;
    MCSResult res = findMCS(mols, &p);
    // std::cerr << "MCS: " << res.SmartsString << " " << res.NumAtoms
    //          << " atoms, " << res.NumBonds << " bonds\n"
    //          << std::endl;
    TEST_ASSERT(res.NumAtoms == 2);
    TEST_ASSERT(res.NumBonds == 1);
    TEST_ASSERT(res.SmartsString == "[#6&!R]-&!@[#6&R]");
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGitHub2731_comment546175466() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing Github #2731 comment 546175466" << std::endl;

  {
    std::vector<ROMOL_SPTR> mols;
    const char *smi[] = {"C1=CC2=CC3=CC=CC=C3N=C2C=C1",
                         "C1=CC=C2N=C3C=NC=CC3=CC2=C1"};

    for (auto &i : smi) {
      auto m = SmilesToMol(getSmilesOnly(i));
      TEST_ASSERT(m);

      mols.emplace_back(m);
    }
    {
      MCSParameters p;
      p.BondCompareParameters.MatchFusedRingsStrict = true;
      MCSResult res = findMCS(mols, &p);
      // std::cerr << "MCS: " << res.SmartsString << " " << res.NumAtoms
      //          << " atoms, " << res.NumBonds << " bonds\n"
      //          << std::endl;
      TEST_ASSERT(res.NumAtoms == 13);
      TEST_ASSERT(res.NumBonds == 14);
      TEST_ASSERT(
          res.SmartsString ==
          "[#6;r6,!R1]1:&@[#6;r6,!R1]:&@[#6]2:&@[#6;r6,!R1]:&@[#6](:&@[#6]:&@[#6]):&@[#6](:&@[#7;r6,!R1]:&@[#6]:&@2:&@[#6;r6,!R1]:&@[#6;r6,!R1]:&@1):&@[#6]");
    }
    {
      MCSParameters p;
      p.BondCompareParameters.CompleteRingsOnly = true;
      p.BondCompareParameters.MatchFusedRingsStrict = true;
      MCSResult res = findMCS(mols, &p);
      // std::cerr << "MCS: " << res.SmartsString << " " << res.NumAtoms
      //          << " atoms, " << res.NumBonds << " bonds\n"
      //          << std::endl;
      TEST_ASSERT(res.NumAtoms == 10);
      TEST_ASSERT(res.NumBonds == 11);
      TEST_ASSERT(
          res.SmartsString ==
          "[#6;r6,!R1]1:&@[#6;r6,!R1]:&@[#6]2:&@[#6;r6,!R1]:&@[#6]:&@[#6]:&@[#7;r6,!R1]:&@[#6]:&@2:&@[#6;r6,!R1]:&@[#6;r6,!R1]:&@1");
    }
  }
  {
    std::vector<ROMOL_SPTR> mols;
    const char *smi[] = {"C12CCC1C1CCCCC12", "C12CC1CC1CCC3CC123"};

    for (auto &i : smi) {
      auto m = SmilesToMol(getSmilesOnly(i));
      TEST_ASSERT(m);

      mols.emplace_back(m);
    }
    {
      MCSParameters p;
      MCSResult res = findMCS(mols, &p);
      // std::cerr << "MCS: " << res.SmartsString << " " << res.NumAtoms
      //          << " atoms, " << res.NumBonds << " bonds\n"
      //          << std::endl;
      TEST_ASSERT(res.NumAtoms == 10);
      TEST_ASSERT(res.NumBonds == 11);
      TEST_ASSERT(res.SmartsString ==
                  "[#6]1-[#6]-[#6]-[#6]-[#6]2-[#6]-1-[#6]-[#6]-[#6]-[#6]-2");
    }
    {
      MCSParameters p;
      p.BondCompareParameters.MatchFusedRings = true;
      MCSResult res = findMCS(mols, &p);
      // std::cerr << "MCS: " << res.SmartsString << " " << res.NumAtoms
      //          << " atoms, " << res.NumBonds << " bonds\n"
      //          << std::endl;
      TEST_ASSERT(res.NumAtoms == 10);
      TEST_ASSERT(res.NumBonds == 10);
      TEST_ASSERT(res.SmartsString ==
                  "[#6](-[#6]-[#6])(-[#6])-[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1");
    }
    {
      MCSParameters p;
      p.BondCompareParameters.MatchFusedRingsStrict = true;
      MCSResult res = findMCS(mols, &p);
      // std::cerr << "MCS: " << res.SmartsString << " " << res.NumAtoms
      //          << " atoms, " << res.NumBonds << " bonds\n"
      //          << std::endl;
      TEST_ASSERT(res.NumAtoms == 10);
      TEST_ASSERT(res.NumBonds == 9);
      TEST_ASSERT(
          res.SmartsString ==
          "[#6](-&@[#6]-&@[#6])(-&@[#6])-&@[#6](-&@[#6]-&@[#6]-&@[#6])-&@[#6]-&@[#6]");
    }
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testQueryMolVsSmarts() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Testing QueryMol vs SmartsString" << std::endl;

  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {"C1CCC2CCCCC12", "C12CC1C1C3CCCC3CCC12"};

  for (auto &i : smi) {
    auto m = SmilesToMol(getSmilesOnly(i));
    TEST_ASSERT(m);

    mols.emplace_back(m);
  }
  MCSParameters p;
  p.BondCompareParameters.MatchFusedRingsStrict = true;
  MCSResult res = findMCS(mols, &p);
  // std::cerr << "MCS: " << res.SmartsString << " " << res.NumAtoms
  //          << " atoms, " << res.NumBonds << " bonds\n"
  //          << std::endl;
  TEST_ASSERT(res.NumAtoms == 9);
  TEST_ASSERT(res.NumBonds == 10);
  ROMOL_SPTR smartsMol(SmartsToMol(res.SmartsString));
  std::vector<MatchVectType> matchVectFromQueryMol;
  std::vector<MatchVectType> matchVectFromSmartsMol;
  TEST_ASSERT(SubstructMatch(*mols[0], *res.QueryMol, matchVectFromQueryMol) ==
              1);
  TEST_ASSERT(SubstructMatch(*mols[0], *smartsMol, matchVectFromSmartsMol) ==
              1);
  TEST_ASSERT(SubstructMatch(*mols[1], *res.QueryMol, matchVectFromQueryMol) ==
              1);
  TEST_ASSERT(SubstructMatch(*mols[1], *smartsMol, matchVectFromSmartsMol) ==
              1);

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testCompareNonExistent() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testAtomCompareNonExistent" << std::endl;

  std::vector<ROMOL_SPTR> mols;
  const char *smi[] = {"C", "CC"};

  for (auto &i : smi) {
    auto m = SmilesToMol(getSmilesOnly(i));
    TEST_ASSERT(m);

    mols.emplace_back(m);
  }
  {
    MCSParameters p;
    bool hasThrown = false;
    try {
      p.setMCSAtomTyperFromEnum(static_cast<AtomComparator>(99));
    } catch (const std::runtime_error &e) {
      BOOST_LOG(rdInfoLog) << e.what() << std::endl;
      hasThrown = true;
    }
    TEST_ASSERT(hasThrown);
  }
  {
    MCSParameters p;
    bool hasThrown = false;
    try {
      p.setMCSAtomTyperFromConstChar("hello");
    } catch (const std::runtime_error &e) {
      BOOST_LOG(rdInfoLog) << e.what() << std::endl;
      hasThrown = true;
    }
    TEST_ASSERT(!hasThrown);
  }
  {
    MCSParameters p;
    bool hasThrown = false;
    try {
      p.setMCSBondTyperFromEnum(static_cast<BondComparator>(99));
    } catch (const std::runtime_error &e) {
      BOOST_LOG(rdInfoLog) << e.what() << std::endl;
      hasThrown = true;
    }
    TEST_ASSERT(hasThrown);
  }
  {
    MCSParameters p;
    bool hasThrown = false;
    try {
      p.setMCSBondTyperFromConstChar("hello");
    } catch (const std::runtime_error &e) {
      BOOST_LOG(rdInfoLog) << e.what() << std::endl;
      hasThrown = true;
    }
    TEST_ASSERT(!hasThrown);
  }
}

void testGitHub3095() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testGitHub3095" << std::endl;

  {
    std::vector<ROMOL_SPTR> mols;
    const char *smi[] = {"C", "OC"};

    for (auto &i : smi) {
      auto m = SmilesToMol(getSmilesOnly(i));
      TEST_ASSERT(m);

      mols.emplace_back(m);
    }
    MCSParameters p;
    MCSResult res = findMCS(mols, &p);
    TEST_ASSERT(res.NumAtoms == 1);
    TEST_ASSERT(res.NumBonds == 0);
    TEST_ASSERT(res.SmartsString == "[#6]");
  }
  {
    std::vector<ROMOL_SPTR> mols;
    const char *smi[] = {"C1CC1", "OC"};

    for (auto &i : smi) {
      auto m = SmilesToMol(getSmilesOnly(i));
      TEST_ASSERT(m);

      mols.emplace_back(m);
    }
    {
      MCSParameters p;
      p.AtomCompareParameters.RingMatchesRingOnly = true;
      MCSResult res = findMCS(mols, &p);
      TEST_ASSERT(res.NumAtoms == 0);
      TEST_ASSERT(res.NumBonds == 0);
      TEST_ASSERT(res.SmartsString.empty());
    }
    {
      MCSParameters p;
      p.BondCompareParameters.RingMatchesRingOnly = true;
      MCSResult res = findMCS(mols, &p);
      TEST_ASSERT(res.NumAtoms == 1);
      TEST_ASSERT(res.NumBonds == 0);
      TEST_ASSERT(res.SmartsString == "[#6]");
    }
  }
  {
    std::vector<ROMOL_SPTR> mols;
    const char *smi[] = {"C1CC1", "C1CCC1"};

    for (auto &i : smi) {
      auto m = SmilesToMol(getSmilesOnly(i));
      TEST_ASSERT(m);

      mols.emplace_back(m);
    }
    {
      MCSParameters p;
      p.BondCompareParameters.CompleteRingsOnly = true;
      MCSResult res = findMCS(mols, &p);
      TEST_ASSERT(res.NumAtoms == 0);
      TEST_ASSERT(res.NumBonds == 0);
      TEST_ASSERT(res.SmartsString.empty());
    }
  }
  {
    std::vector<ROMOL_SPTR> mols;
    const char *smi[] = {"CC1CC1", "CN1CCC1"};

    for (auto &i : smi) {
      auto m = SmilesToMol(getSmilesOnly(i));
      TEST_ASSERT(m);

      mols.emplace_back(m);
    }
    {
      MCSParameters p;
      p.AtomCompareParameters.RingMatchesRingOnly = true;
      p.BondCompareParameters.RingMatchesRingOnly = true;
      p.BondCompareParameters.CompleteRingsOnly = true;
      MCSResult res = findMCS(mols, &p);
      TEST_ASSERT(res.NumAtoms == 1);
      TEST_ASSERT(res.NumBonds == 0);
      TEST_ASSERT(res.SmartsString == "[#6&!R]");
    }
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGitHub3458() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testGitHub3458" << std::endl;

  {
    std::vector<ROMOL_SPTR> mols;
    const char *smi[] = {
        "Brc1cccc(Nc2ncnc3cc4ccccc4cc23)c1",
        "CCOc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OCC",
        "CN(C)c1cc2c(Nc3cccc(Br)c3)ncnc2cn1",
        "CNc1cc2c(Nc3cccc(Br)c3)ncnc2cn1",
        "Brc1cccc(Nc2ncnc3cc4[nH]cnc4cc23)c1",
        "Cn1cnc2cc3ncnc(Nc4cccc(Br)c4)c3cc21",
        "Cn1cnc2cc3c(Nc4cccc(Br)c4)ncnc3cc21",
        "COc1cc2ncnc(Nc3cccc(Br)c3)c2cc1OC",
        "C#CCNC/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)c(C#N)cnc2cc1OCC",
        "C=CC(=O)Nc1ccc2ncnc(Nc3cc(Cl)c(Cl)cc3F)c2c1"};

    for (auto &i : smi) {
      auto m = SmilesToMol(getSmilesOnly(i));
      TEST_ASSERT(m);

      mols.emplace_back(m);
    }
    {
      MCSParameters p;
      p.AtomCompareParameters.RingMatchesRingOnly = true;
      p.BondCompareParameters.RingMatchesRingOnly = true;
      p.Threshold = 0.8;
      MCSResult res = findMCS(mols, &p);
      TEST_ASSERT(res.NumAtoms == 17);
      TEST_ASSERT(res.NumBonds == 18);
      TEST_ASSERT(
          res.SmartsString ==
          "[#6&R]:&@[#6&R]:&@[#6]1:&@[#6](-&!@[#7&!R]-&!@[#6]2:&@[#6]:&@[#6]:&@[#6]:&@[#6](:&@[#6]:&@2)-&!@[#35&!R]):&@[#7]:&@[#6]:&@[#7]:&@[#6]:&@1:&@[#6&R]");
    }
    {
      MCSParameters p;
      p.AtomCompareParameters.RingMatchesRingOnly = true;
      p.BondCompareParameters.RingMatchesRingOnly = true;
      p.Threshold = 1.0;
      MCSResult res = findMCS(mols, &p);
      TEST_ASSERT(res.NumAtoms == 15);
      TEST_ASSERT(res.NumBonds == 15);
      TEST_ASSERT(
          res.SmartsString ==
          "[#6&R]:&@[#6&R]:&@[#6&R](:&@[#6&R]-&!@[#7&!R]-&!@[#6]1:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@1):&@[#6&R](:&@[#7&R]:&@[#6&R]):&@[#6&R]");
    }
  }

  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGitHub3693() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testGitHub3693" << std::endl;

  {
    std::vector<ROMOL_SPTR> mols = {
        "Nc1ccc(O)cc1c1ccc2ccccc2c1"_smiles,
        "Oc1cnc(NC2CCC2)c(c1)c1ccc2ccccc2c1"_smiles};

    {
      MCSParameters p;
      MCSResult res = findMCS(mols, &p);
      TEST_ASSERT(res.NumAtoms == 17);
      TEST_ASSERT(res.NumBonds == 18);
      TEST_ASSERT(res.SmartsString ==
                  "[#7]-,:[#6]:[#6](:[#6]:[#6](:[#6])-[#8])-[#6]1:[#6]:[#6]:[#"
                  "6]2:[#6](:[#6]:1):[#6]:[#6]:[#6]:[#6]:2");
    }
    {
      MCSParameters p;
      p.BondCompareParameters.CompleteRingsOnly = true;
      MCSResult res = findMCS(mols, &p);
      TEST_ASSERT(res.NumAtoms == 11);
      TEST_ASSERT(res.NumBonds == 12);
      TEST_ASSERT(res.SmartsString ==
                  "[#6]-&!@[#6]1:&@[#6]:&@[#6]:&@[#6]2:&@[#6](:&@[#6]:&@1):&@[#"
                  "6]:&@[#6]:&@[#6]:&@[#6]:&@2");
    }
    {
      MCSParameters p;
      p.AtomCompareParameters.CompleteRingsOnly = true;
      MCSResult res = findMCS(mols, &p);
      TEST_ASSERT(res.NumAtoms == 10);
      TEST_ASSERT(res.NumBonds == 11);
      TEST_ASSERT(
          res.SmartsString ==
          "[#6]1:&@[#6]:&@[#6]:&@[#6]2:&@[#6](:&@[#6]:&@1):&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@2");
    }
    {
      MCSParameters p;
      p.AtomCompareParameters.CompleteRingsOnly = true;
      p.BondCompareParameters.CompleteRingsOnly = true;
      MCSResult res = findMCS(mols, &p);
      TEST_ASSERT(res.NumAtoms == 10);
      TEST_ASSERT(res.NumBonds == 11);
      TEST_ASSERT(
          res.SmartsString ==
          "[#6]1:&@[#6]:&@[#6]:&@[#6]2:&@[#6](:&@[#6]:&@1):&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@2");
    }
    {
      MCSParameters p;
      p.AtomCompareParameters.CompleteRingsOnly = true;
      // this will automatically be set to true
      p.BondCompareParameters.CompleteRingsOnly = false;
      MCSResult res = findMCS(mols, &p);
      TEST_ASSERT(res.NumAtoms == 10);
      TEST_ASSERT(res.NumBonds == 11);
      TEST_ASSERT(
          res.SmartsString ==
          "[#6]1:&@[#6]:&@[#6]:&@[#6]2:&@[#6](:&@[#6]:&@1):&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@2");
    }
  }
  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGitHub3886() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "testGitHub3886" << std::endl;

  std::vector<ROMOL_SPTR> mols = {"c1cccnc1"_smiles, "Fc1ccccc1"_smiles};

  MCSParameters p;
  p.InitialSeed = "c1ccc*c1";
  std::stringstream captureLog;
  rdWarningLog->SetTee(captureLog);
  MCSResult res = findMCS(mols, &p);
  rdWarningLog->ClearTee();
  TEST_ASSERT(captureLog.str().find("The provided InitialSeed is not an MCS") !=
              std::string::npos);
  TEST_ASSERT(res.NumAtoms == 5);
  TEST_ASSERT(res.NumBonds == 4);
  TEST_ASSERT(res.SmartsString == "[#6](:[#6]:[#6]:[#6]):[#6]");
  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testAtomCompareCompleteRingsOnly() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "When AtomCompareParameters.CompleteRingsOnly is true single atoms "
         "which are part of a ring in one of the molecules should not be "
         "included in MCS"
      << std::endl;
  std::vector<ROMOL_SPTR> mols = {"C1CCCC1C"_smiles, "C1CCCC1C1CCCCC1"_smiles};
  {
    MCSParameters p;
    p.AtomCompareParameters.CompleteRingsOnly = true;
    MCSResult res = findMCS(mols, &p);
    TEST_ASSERT(res.NumAtoms == 5);
    TEST_ASSERT(res.NumBonds == 5);
    TEST_ASSERT(res.SmartsString == "[#6]1-&@[#6]-&@[#6]-&@[#6]-&@[#6]-&@1");
  }
  {
    MCSParameters p;
    p.AtomCompareParameters.CompleteRingsOnly = true;
    // this will automatically be set to true
    p.AtomCompareParameters.RingMatchesRingOnly = false;
    MCSResult res = findMCS(mols, &p);
    TEST_ASSERT(res.NumAtoms == 5);
    TEST_ASSERT(res.NumBonds == 5);
    TEST_ASSERT(res.SmartsString == "[#6]1-&@[#6]-&@[#6]-&@[#6]-&@[#6]-&@1");
  }
  BOOST_LOG(rdInfoLog) << "============================================"
                       << std::endl;
  BOOST_LOG(rdInfoLog) << "\tdone" << std::endl;
}

void testGitHub4498() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "GitHub Issue #4498: FindMCS may leave mols with "
                          "fake (empty) ring info"
                       << std::endl;
  {
    // Test with rings
    std::vector<ROMOL_SPTR> mols = {"NC1=CC(N)=C(N)C=C1"_smiles,
                                    "NC1=CC(N)=C(N)C=C1"_smiles};

    TEST_ASSERT(mols[0]);
    TEST_ASSERT(mols[1]);

    mols[1]->getRingInfo()->reset();
    TEST_ASSERT(mols[0]->getRingInfo()->isInitialized() == true);
    TEST_ASSERT(mols[1]->getRingInfo()->isInitialized() == false);

    MCSResult res = findMCS(mols);
    TEST_ASSERT(res.NumAtoms == 9);

    TEST_ASSERT(mols[0]->getRingInfo()->isInitialized() == true);
    TEST_ASSERT(mols[0]->getRingInfo()->numRings() == 1);

    TEST_ASSERT(mols[1]->getRingInfo()->isInitialized() == false);
  }
  {
    // Test without rings
    std::vector<ROMOL_SPTR> mols = {"NC=CC(N)=C(N)C=C"_smiles,
                                    "NC=CC(N)=C(N)C=C"_smiles};

    TEST_ASSERT(mols[0]);
    TEST_ASSERT(mols[1]);

    mols[1]->getRingInfo()->reset();
    TEST_ASSERT(mols[0]->getRingInfo()->isInitialized() == true);
    TEST_ASSERT(mols[1]->getRingInfo()->isInitialized() == false);

    MCSResult res = findMCS(mols);
    TEST_ASSERT(res.NumAtoms == 9);

    TEST_ASSERT(mols[0]->getRingInfo()->isInitialized() == true);
    TEST_ASSERT(mols[0]->getRingInfo()->numRings() == 0);

    TEST_ASSERT(mols[1]->getRingInfo()->isInitialized() == false);
  }
}

void testBondStereo() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "FindMCS should check bond stereo" << std::endl;
  {
    std::vector<ROMOL_SPTR> mols = {"CC\\C=C/CC"_smiles, "CC\\C=C\\CC"_smiles};

    MCSParameters p;
    p.BondCompareParameters.MatchStereo = false;
    MCSResult mcs_resf = findMCS(mols, &p);
    std::cout << "MCS MatchStereo false : " << mcs_resf.SmartsString << " "
              << mcs_resf.NumAtoms << " atoms, " << mcs_resf.NumBonds
              << " bonds\n";
    TEST_ASSERT(mcs_resf.NumAtoms == 6);
    TEST_ASSERT(mcs_resf.NumBonds == 5);

    p.BondCompareParameters.MatchStereo = true;
    MCSResult mcs_rest = findMCS(mols, &p);
    std::cout << "MCS MatchStereo true  : " << mcs_rest.SmartsString << " "
              << mcs_rest.NumAtoms << " atoms, " << mcs_rest.NumBonds
              << " bonds\n";
    TEST_ASSERT(mcs_rest.NumAtoms == 3);
    TEST_ASSERT(mcs_rest.NumBonds == 2);
    TEST_ASSERT(mcs_resf.SmartsString != mcs_rest.SmartsString);
  }
  {
    std::vector<ROMOL_SPTR> mols = {"CC\\C=C/CC"_smiles, "CCC=CCC"_smiles};
    MCSParameters p;
    p.BondCompareParameters.MatchStereo = false;
    MCSResult mcs_resf = findMCS(mols, &p);
    std::cout << "MCS MatchStereo false : " << mcs_resf.SmartsString << " "
              << mcs_resf.NumAtoms << " atoms, " << mcs_resf.NumBonds
              << " bonds\n";
    TEST_ASSERT(mcs_resf.NumAtoms == 6);
    TEST_ASSERT(mcs_resf.NumBonds == 5);

    p.BondCompareParameters.MatchStereo = true;
    MCSResult mcs_rest = findMCS(mols, &p);
    std::cout << "MCS MatchStereo true  : " << mcs_rest.SmartsString << " "
              << mcs_rest.NumAtoms << " atoms, " << mcs_rest.NumBonds
              << " bonds\n";
    TEST_ASSERT(mcs_rest.NumAtoms == 3);
    TEST_ASSERT(mcs_rest.NumBonds == 2);
    TEST_ASSERT(mcs_resf.SmartsString != mcs_rest.SmartsString);
  }
  {
    std::vector<ROMOL_SPTR> mols = {"CCC=CCC"_smiles, "CCC=CCC"_smiles};
    MCSParameters p;
    p.BondCompareParameters.MatchStereo = false;
    MCSResult mcs_resf = findMCS(mols, &p);
    std::cout << "MCS MatchStereo false : " << mcs_resf.SmartsString << " "
              << mcs_resf.NumAtoms << " atoms, " << mcs_resf.NumBonds
              << " bonds\n";
    TEST_ASSERT(mcs_resf.NumAtoms == 6);
    TEST_ASSERT(mcs_resf.NumBonds == 5);

    p.BondCompareParameters.MatchStereo = true;
    MCSResult mcs_rest = findMCS(mols, &p);
    std::cout << "MCS MatchStereo true  : " << mcs_rest.SmartsString << " "
              << mcs_rest.NumAtoms << " atoms, " << mcs_rest.NumBonds
              << " bonds\n";
    TEST_ASSERT(mcs_rest.NumAtoms == 6);
    TEST_ASSERT(mcs_rest.NumBonds == 5);
    TEST_ASSERT(mcs_resf.SmartsString == mcs_rest.SmartsString);
  }
}

void testSlowCompleteRingsOnly() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "FindMCS should be more efficient when" << std::endl;
  BOOST_LOG(rdInfoLog) << "CompleteRingsOnly is set" << std::endl;
  std::vector<ROMOL_SPTR> mols = {
      "CC(CC(C1CC2CCC(CC3CC(CC3C)C(C)(C)C)CC2C1)C1CCC(C)(C)CC1)C1CCCC1C"_smiles,
      "CC(CC(C1CC2CCC(CC3CC(CC3C)C(C)(C)C)CC2C1)C1CCC(C)(C)C1)C1CCCC1C"_smiles};

  MCSParameters p;
  p.BondCompareParameters.CompleteRingsOnly = true;
  p.AtomTyper = MCSAtomCompareAny;
  p.BondTyper = MCSBondCompareAny;
#ifdef NDEBUG
  p.Timeout = 10;
#else  // allow more time when building in debug mode
  p.Timeout = 60;
#endif
  auto mcs = findMCS(mols, &p);

  TEST_ASSERT(!mcs.Canceled);
  TEST_ASSERT(mcs.NumAtoms == 31);
  TEST_ASSERT(mcs.NumBonds == 34);
}

void testCustomShouldAcceptMCS() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "Use custom ShouldAcceptMCS function" << std::endl;
  std::vector<ROMOL_SPTR> mols = {"OC(=O)CCc1c(C)ccc2ccccc12"_smiles,
                                  "c1(C)cc(cccc2)c2c(*)c1"_smiles};
  struct customMCSFunction {
    static bool atomCompare(const MCSAtomCompareParameters &p,
                            const ROMol &mol1, unsigned int atom1,
                            const ROMol &mol2, unsigned int atom2, void *) {
      const auto a1 = mol1.getAtomWithIdx(atom1);
      const auto a2 = mol2.getAtomWithIdx(atom2);
      bool a1IsDummy = (a1->getAtomicNum() == 0);
      bool a2IsDummy = (a2->getAtomicNum() == 0);
      if (a1IsDummy ^ a2IsDummy) {
        const Atom *atoms[] = {a1, a2};
        unsigned int dummyAtomIdx = a1IsDummy ? 0 : 1;
        unsigned int otherAtomIdx = 1 - dummyAtomIdx;
        return atoms[dummyAtomIdx]->getDegree() == 1 &&
               atoms[otherAtomIdx]->getDegree() > 1;
      }
      if (a1->getAtomicNum() != a2->getAtomicNum()) {
        return false;
      }
      return checkAtomRingMatch(p, mol1, atom1, mol2, atom2);
    }

    static bool shouldAcceptMCS(const ROMol &query, const ROMol &target,
                                const MatchVectType &atomIdxMatch,
                                const MatchVectType &bondIdxMatch,
                                const MCSParameters *) {
      for (const auto &atomIdxPair : atomIdxMatch) {
        if (query.getAtomWithIdx(atomIdxPair.first)->getAtomicNum() == 0) {
          return true;
        }
        if (target.getAtomWithIdx(atomIdxPair.second)->getAtomicNum() == 0) {
          return true;
        }
      }
      boost::dynamic_bitset<> queryAtomIndices(query.getNumAtoms());
      std::for_each(atomIdxMatch.begin(), atomIdxMatch.end(),
                    [&queryAtomIndices](const auto &pair) {
                      queryAtomIndices.set(pair.first);
                    });
      boost::dynamic_bitset<> targetAtomIndices(target.getNumAtoms());
      std::for_each(atomIdxMatch.begin(), atomIdxMatch.end(),
                    [&targetAtomIndices](const auto &pair) {
                      targetAtomIndices.set(pair.second);
                    });
      for (const auto &bondIdxPair : bondIdxMatch) {
        const auto queryBond = query.getBondWithIdx(bondIdxPair.first);
        const auto targetBond = target.getBondWithIdx(bondIdxPair.second);
        TEST_ASSERT(queryBond);
        TEST_ASSERT(targetBond);
        TEST_ASSERT(queryAtomIndices.test(queryBond->getBeginAtomIdx()));
        TEST_ASSERT(queryAtomIndices.test(queryBond->getEndAtomIdx()));
        TEST_ASSERT(targetAtomIndices.test(targetBond->getBeginAtomIdx()));
        TEST_ASSERT(targetAtomIndices.test(targetBond->getEndAtomIdx()));
      }
      return false;
    }
  };

  MCSParameters p;
  p.AtomTyper = customMCSFunction::atomCompare;
  p.BondTyper = MCSBondCompareAny;
  p.ShouldAcceptMCS = customMCSFunction::shouldAcceptMCS;
  p.AtomCompareParameters.RingMatchesRingOnly = false;
  p.AtomCompareParameters.CompleteRingsOnly = false;
  p.BondCompareParameters.RingMatchesRingOnly = true;
  p.BondCompareParameters.CompleteRingsOnly = false;
  auto mcs = findMCS(mols, &p);
  TEST_ASSERT(mcs.NumAtoms == 11);
  TEST_ASSERT(mcs.NumBonds == 12);
  TEST_ASSERT(
      mcs.SmartsString ==
      "[#6]1:&@[#6]:&@[#6]2:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@2:&@[#6](:&@[#6]:&@1)-&!@[#0,#6]");
}

void testGitHub5510() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test GitHub issue #5510" << std::endl;
  {
    std::vector<ROMOL_SPTR> mols = {"C"_smiles, "O"_smiles, "N"_smiles};
    auto mcs = findMCS(mols);
    TEST_ASSERT(mcs.NumAtoms == 0);
  }
  {
    std::vector<ROMOL_SPTR> mols = {"C"_smiles, "O"_smiles, "N"_smiles,
                                    "CCC"_smiles};
    auto mcs = findMCS(mols);
    TEST_ASSERT(mcs.NumAtoms == 0);
  }
}

void testDegenerateMCS() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test retrieving all degenerate MCSs" << std::endl;
  std::string para("CC1:C:C:C(N):C:C:1");
  std::string ortho("CC1:C:C:C:C:C:1N");
  {
    std::vector<ROMOL_SPTR> mols1 = {
        "Nc1ccc(cc1)C-Cc1c(N)cccc1"_smiles,
        "Nc1ccc(cc1)C=Cc1c(N)cccc1"_smiles,
    };
    std::vector<ROMOL_SPTR> mols2 = {
        "Nc1ccccc1C-Cc1ccc(N)cc1"_smiles,
        "Nc1ccccc1C=Cc1ccc(N)cc1"_smiles,
    };
    MCSParameters p;
    auto mcs1 = findMCS(mols1);
    TEST_ASSERT(mcs1.NumAtoms == 8);
    TEST_ASSERT(mcs1.NumBonds == 8);
    TEST_ASSERT(!mcs1.SmartsString.empty());
    TEST_ASSERT(mcs1.QueryMol);
    auto mcs1Smiles = MolToSmiles(*ROMOL_SPTR(SmartsToMol(mcs1.SmartsString)));
    TEST_ASSERT(mcs1Smiles == para);
    TEST_ASSERT(mcs1.DegenerateSmartsQueryMolDict.empty());
    auto mcs2 = findMCS(mols2);
    TEST_ASSERT(mcs2.NumAtoms == 8);
    TEST_ASSERT(mcs2.NumBonds == 8);
    TEST_ASSERT(!mcs2.SmartsString.empty());
    TEST_ASSERT(mcs2.QueryMol);
    auto mcs2Smiles = MolToSmiles(*ROMOL_SPTR(SmartsToMol(mcs2.SmartsString)));
    TEST_ASSERT(mcs2Smiles == ortho);
    TEST_ASSERT(mcs2.DegenerateSmartsQueryMolDict.empty());
    TEST_ASSERT(mcs1Smiles != mcs2Smiles);
    TEST_ASSERT(MolToSmiles(*mols1[0]) == MolToSmiles(*mols2[0]));
    TEST_ASSERT(MolToSmiles(*mols1[1]) == MolToSmiles(*mols2[1]));
    p.StoreAll = true;
    mcs1 = findMCS(mols1, &p);
    TEST_ASSERT(mcs1.NumAtoms == 8);
    TEST_ASSERT(mcs1.NumBonds == 8);
    TEST_ASSERT(mcs1.SmartsString.empty());
    TEST_ASSERT(!mcs1.QueryMol);
    TEST_ASSERT(mcs1.DegenerateSmartsQueryMolDict.size() == 2);
    mcs2 = findMCS(mols2, &p);
    TEST_ASSERT(mcs2.NumAtoms == 8);
    TEST_ASSERT(mcs2.NumBonds == 8);
    TEST_ASSERT(mcs2.SmartsString.empty());
    TEST_ASSERT(!mcs2.QueryMol);
    TEST_ASSERT(mcs2.DegenerateSmartsQueryMolDict.size() == 2);
    auto smartsToCanonicalSmiles = [](const auto &pair) {
      return MolToSmiles(*ROMOL_SPTR(SmartsToMol(pair.first)));
    };
    std::set<std::string> degSmiles1;
    std::transform(mcs1.DegenerateSmartsQueryMolDict.begin(),
                   mcs1.DegenerateSmartsQueryMolDict.end(),
                   std::inserter(degSmiles1, degSmiles1.end()),
                   smartsToCanonicalSmiles);
    std::set<std::string> degSmiles2;
    std::transform(mcs2.DegenerateSmartsQueryMolDict.begin(),
                   mcs2.DegenerateSmartsQueryMolDict.end(),
                   std::inserter(degSmiles2, degSmiles2.end()),
                   smartsToCanonicalSmiles);
    TEST_ASSERT(degSmiles1.size() == 2);
    TEST_ASSERT(degSmiles1 == degSmiles2);
  }
  {
    std::vector<ROMOL_SPTR> mols = {
        "Nc1ccc(cc1)C-C"_smiles,
        "Nc1ccc(cc1)C=C"_smiles,
    };
    MCSParameters p;
    auto mcs1 = findMCS(mols);
    TEST_ASSERT(mcs1.NumAtoms == 8);
    TEST_ASSERT(mcs1.NumBonds == 8);
    TEST_ASSERT(mcs1.QueryMol);
    TEST_ASSERT(!mcs1.SmartsString.empty());
    TEST_ASSERT(MolToSmiles(*ROMOL_SPTR(SmartsToMol(mcs1.SmartsString))) ==
                para);
    TEST_ASSERT(mcs1.DegenerateSmartsQueryMolDict.empty());
    p.StoreAll = true;
    auto mcs2 = findMCS(mols, &p);
    TEST_ASSERT(mcs2.NumAtoms == 8);
    TEST_ASSERT(mcs2.NumBonds == 8);
    TEST_ASSERT(!mcs2.QueryMol);
    TEST_ASSERT(mcs2.SmartsString.empty());
    TEST_ASSERT(mcs2.DegenerateSmartsQueryMolDict.size() == 1);
    TEST_ASSERT(MolToSmiles(*ROMOL_SPTR(SmartsToMol(
                    mcs2.DegenerateSmartsQueryMolDict.begin()->first))) ==
                para);
  }
}

void testGitHub5411() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test GitHub issue #5411" << std::endl;
  ROMOL_SPTR mola(
      SmilesToMol("[NH3+]CCCCn1c(SCCc2c[nH]c3ccccc23)nnc1-c1cccc2ccccc12"));
  ROMOL_SPTR molb(
      SmilesToMol("[NH3+]CCCCn1c(SCCc2c[nH]c3ccccc23)nnc1-c1ccc2ccccc2n1"));
  std::string expectedSmarts(
      "[#7]-&!@[#6]-&!@[#6]-&!@[#6]-&!@[#6]-&!@[#7;r5,!R1]1:&@[#6;r5,!R1](-&!@[#16]-&!@[#6]-&!@[#6]-&!@[#6;r5,!R1]2:&@[#6;r5,!R1]:&@[#7;r5,!R1]:&@[#6]3:&@[#6]:&@2:&@[#6;r6,!R1]:&@[#6;r6,!R1]:&@[#6;r6,!R1]:&@[#6;r6,!R1]:&@3):&@[#7;r5,!R1]:&@[#7;r5,!R1]:&@[#6;r5,!R1]:&@1-&!@[#6]:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@[#6]");
  MCSParameters p;
  p.BondCompareParameters.MatchFusedRings = true;
  p.BondCompareParameters.MatchFusedRingsStrict = true;
  std::vector<ROMOL_SPTR> mola_molb{mola, molb};
  std::vector<ROMOL_SPTR> molb_mola{molb, mola};
  auto mcs1 = findMCS(mola_molb, &p);
  TEST_ASSERT(mcs1.NumAtoms == 31);
  TEST_ASSERT(mcs1.NumBonds == 33);
  TEST_ASSERT(mcs1.QueryMol);
  TEST_ASSERT(mcs1.SmartsString == expectedSmarts);
  auto mcs2 = findMCS(molb_mola, &p);
  TEST_ASSERT(mcs2.NumAtoms == mcs1.NumAtoms);
  TEST_ASSERT(mcs2.NumBonds == mcs1.NumBonds);
  TEST_ASSERT(mcs2.QueryMol);
  TEST_ASSERT(mcs2.SmartsString == mcs1.SmartsString);
}

void testGitHub5440() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test GitHub issue #5440" << std::endl;
  ROMOL_SPTR mola(
      SmilesToMol("c1cccc(c12)sc(C([O-])=O)c2CCCOc3cccc(c34)cccc4"));
  ROMOL_SPTR molb(
      SmilesToMol("c1cccc(c12)oc(C([O-])=O)c2CCCOc3cccc(c34)cccc4"));
  TEST_ASSERT(mola->getNumAtoms() == 26 && mola->getNumBonds() == 29);
  TEST_ASSERT(molb->getNumAtoms() == mola->getNumAtoms() &&
              molb->getNumBonds() == mola->getNumBonds());
  std::string seed(
      "[#6&!R](=&!@[#8&!R])(-&!@[#8&!R])-&!@[#6&R]1:&@[#6&R](-&!@[#6&!R]-&!@[#6&!R]-&!@[#6&!R]-&!@[#8&!R]-&!@[#6&R]2:&@[#6&R]:&@[#6&R]:&@[#6&R]:&@[#6&R]:&@[#6&R]:&@2):&@[#6&R]2:&@[#6&R](:&@[#7,#8,#16;R]:&@1):&@[#6&R]:&@[#6&R]:&@[#6&R]:&@[#6&R]:&@2");
  MCSParameters p;
  p.MaximizeBonds = false;
  p.AtomTyper = MCSAtomCompareAny;
  p.AtomCompareParameters.RingMatchesRingOnly = true;
  p.AtomCompareParameters.CompleteRingsOnly = true;
  p.BondTyper = MCSBondCompareAny;
  p.BondCompareParameters.CompleteRingsOnly = true;
  p.BondCompareParameters.MatchFusedRings = true;
  p.BondCompareParameters.MatchFusedRingsStrict = true;
  ROMOL_SPTR seedMol(SmartsToMol(seed));
  TEST_ASSERT(!SubstructMatch(*mola, *seedMol).empty());
  TEST_ASSERT(!SubstructMatch(*molb, *seedMol).empty());
  p.InitialSeed = seed;
  std::vector<ROMOL_SPTR> mola_molb{mola, molb};
  auto mcs1 = findMCS(mola_molb, &p);
  TEST_ASSERT(mcs1.NumAtoms == 26);
  TEST_ASSERT(mcs1.NumBonds == 29);
  TEST_ASSERT(!SubstructMatch(*mola, *mcs1.QueryMol).empty());
  TEST_ASSERT(SubstructMatch(*mola, *mcs1.QueryMol) ==
              SubstructMatch(*molb, *mcs1.QueryMol));
  p.InitialSeed = "";
  auto mcs2 = findMCS(mola_molb, &p);
  TEST_ASSERT(mcs2.NumAtoms == mcs1.NumAtoms);
  TEST_ASSERT(mcs2.NumBonds == mcs1.NumBonds);
  TEST_ASSERT(!SubstructMatch(*mola, *mcs2.QueryMol).empty());
  TEST_ASSERT(SubstructMatch(*mola, *mcs2.QueryMol) ==
              SubstructMatch(*molb, *mcs2.QueryMol));
  p.InitialSeed = "AAAA";
  auto mcs3 = findMCS(mola_molb, &p);
  TEST_ASSERT(mcs3.NumAtoms == mcs1.NumAtoms);
  TEST_ASSERT(mcs3.NumBonds == mcs1.NumBonds);
  TEST_ASSERT(!SubstructMatch(*mola, *mcs3.QueryMol).empty());
  TEST_ASSERT(SubstructMatch(*mola, *mcs3.QueryMol) ==
              SubstructMatch(*molb, *mcs3.QueryMol));
  p.InitialSeed = "a1aaaa1";
  auto mcs4 = findMCS(mola_molb, &p);
  TEST_ASSERT(mcs4.NumAtoms == mcs1.NumAtoms);
  TEST_ASSERT(mcs4.NumBonds == mcs1.NumBonds);
  TEST_ASSERT(!SubstructMatch(*mola, *mcs4.QueryMol).empty());
  TEST_ASSERT(SubstructMatch(*mola, *mcs4.QueryMol) ==
              SubstructMatch(*molb, *mcs4.QueryMol));
}

void testGitHub5457() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test GitHub issue #5457" << std::endl;
  ROMOL_SPTR mola(SmilesToMol(
      "c1cccc(c12)ccc(c2)-c3n(CCC[NH3+])c(nn3)SCCc4c[nH]c(c45)cccc5"));
  ROMOL_SPTR molb(SmilesToMol(
      "c1cccc(c12)sc(c2)-c3n(CCCC[NH3+])c(nn3)SCCc4c[nH]c(c45)cccc5"));
  std::string seed(
      "[#6&!R](-&!@[#6,#7;!R])-&!@[#6&!R]-&!@[#6&!R]-&!@[#7&R]1:&@[#6&R](-&!@[#16&!R]-&!@[#6&!R]-&!@[#6&!R]-&!@[#6&R]2:&@[#6&R]:&@[#7&R]:&@[#6&R]3:&@[#6&R]:&@2:&@[#6&R]:&@[#6&R]:&@[#6&R]:&@[#6&R]:&@3):&@[#7&R]:&@[#7&R]:&@[#6&R]:&@1");
  MCSParameters p;
  p.AtomTyper = MCSAtomCompareAny;
  p.AtomCompareParameters.RingMatchesRingOnly = true;
  p.AtomCompareParameters.CompleteRingsOnly = true;
  p.BondTyper = MCSBondCompareAny;
  p.BondCompareParameters.CompleteRingsOnly = true;
  ROMOL_SPTR seedMol(SmartsToMol(seed));
  TEST_ASSERT(!SubstructMatch(*mola, *seedMol).empty());
  TEST_ASSERT(!SubstructMatch(*molb, *seedMol).empty());
  p.InitialSeed = seed;
  std::vector<ROMOL_SPTR> mola_molb{mola, molb};
  auto mcs = findMCS(mola_molb, &p);
  TEST_ASSERT(mcs.NumAtoms == 21);
  TEST_ASSERT(mcs.NumBonds == 23);
  TEST_ASSERT(
      mcs.SmartsString ==
      "[#6&!R](-&!@[#7,#6;!R])-&!@[#6&!R]-&!@[#6&!R]-&!@[#7]1:&@[#6](-&!@[#16&!R]-&!@[#6&!R]-&!@[#6&!R]-&!@[#6]2:&@[#6]:&@[#7]:&@[#6]3:&@[#6]:&@2:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@3):&@[#7]:&@[#7]:&@[#6]:&@1");
  TEST_ASSERT(!SubstructMatch(*mola, *mcs.QueryMol).empty());
  TEST_ASSERT(!SubstructMatch(*molb, *mcs.QueryMol).empty());
  ;
}

void testGitHub6082() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test GitHub issue #6082" << std::endl;
  std::vector<ROMOL_SPTR> mols{"CCc1ccc2nc(CO)cn2n1"_smiles,
                               "OCc1cc2nc(C3CC3)ccn2n1"_smiles};
  auto mcs = findMCS(mols);
  TEST_ASSERT(mcs.NumAtoms == 11);
  TEST_ASSERT(mcs.NumBonds == 10);
  TEST_ASSERT(mcs.SmartsString ==
              "[#6]-[#6]-[#6]:,-[#6]:[#6]:[#6]:[#7](:[#6]:[#6]-,:[#6]):[#7]");
}

void testGitHub3965() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test GitHub issue #6082" << std::endl;
  // The molecules originally reported in Issue 3965 are larger
  // While the test completes successfully, it takes ~120s,
  // so we use a smaller version which triggered the same neverending
  // behavior with the old code, and now completes in a few seconds
  std::vector<ROMOL_SPTR> mols{
      "Nc1cc2cc(c1)C(=O)N[C@H](c1ccccc1)c1cccc(COCCNC2=O)c1"_smiles,
      "Nc1nc2cc(c1)C(=O)N[C@H](c1ccccc1)c1cccc(COCCNC2=O)c1"_smiles};
  MCSParameters p;
#ifdef NDEBUG
  p.Timeout = 30;
#else  // allow more time when building in debug mode
  p.Timeout = 120;
#endif
  p.BondCompareParameters.CompleteRingsOnly = true;
  auto mcs = findMCS(mols, &p);

  TEST_ASSERT(!mcs.Canceled);
  TEST_ASSERT(mcs.NumAtoms == 26);
  TEST_ASSERT(mcs.NumBonds == 28);
  TEST_ASSERT(
      mcs.SmartsString ==
      "[#6]1:&@[#6]:&@[#6]-&@[#6](=&!@[#8])-&@[#7]-&@[#6](-&!@[#6]2:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@2)-&@[#6]2:&@[#6]:&@[#6]:&@[#6]:&@[#6](-&@[#6]-&@[#8]-&@[#6]-&@[#6]-&@[#7]-&@[#6]-&@1=&!@[#8]):&@[#6]:&@2");
}

void testGitHub6578() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test GitHub issue #6578" << std::endl;
  {
    std::vector<ROMOL_SPTR> mols{
        "Clc1ccc(cc1)c2nn3ncccc3c2c4ccnc(Nc5ccc6OCCOc6c5)n4"_smiles,
        "Cc1ccc2c(c3ccnc(Nc4ccccc4)n3)c(nn2n1)c5ccccc5"_smiles,
        "Cc1ccc2c(c3ccnc(Nc4cccc(c4)C(F)(F)F)n3)c(nn2n1)c5ccccc5"_smiles,
        "c1[nH]ccc1"_smiles, ""_smiles};
    MCSParameters p;
    p.Threshold = 0.5;
    auto mcs = findMCS(mols, &p);
    TEST_ASSERT(mcs.NumAtoms == 29);
    TEST_ASSERT(mcs.NumBonds == 33);
    TEST_ASSERT(
        mcs.SmartsString ==
        "[#6]-[#6]1:[#6]:[#6]:[#6]2:[#7](:[#7]:1):[#7]:[#6](:[#6]:2-[#6]1:[#6]:[#6]:[#7]:[#6](:[#7]:1)-[#7]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1");
  }
  {
    std::vector<ROMOL_SPTR> mols{
        "Clc1ccc(cc1)c2nn3ncccc3c2c4ccnc(Nc5ccc6OCCOc6c5)n4"_smiles, ""_smiles,
        "Cc1ccc2c(c3ccnc(Nc4ccccc4)n3)c(nn2n1)c5ccccc5"_smiles, ""_smiles,
        ""_smiles};
    MCSParameters p;
    p.Threshold = 0.0;
    auto mcs = findMCS(mols, &p);
    TEST_ASSERT(mcs.NumAtoms == 29);
    TEST_ASSERT(mcs.NumBonds == 33);
    TEST_ASSERT(
        mcs.SmartsString ==
        "[#6]-[#6]1:[#6]:[#6]:[#6]2:[#7](:[#7]:1):[#7]:[#6](:[#6]:2-[#6]1:[#6]:[#6]:[#7]:[#6](:[#7]:1)-[#7]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1");
  }
  {
    std::vector<ROMOL_SPTR> mols{""_smiles, ""_smiles, ""_smiles, ""_smiles};
    MCSParameters p;
    p.Threshold = 0.0;
    auto mcs = findMCS(mols, &p);
    TEST_ASSERT(mcs.NumAtoms == 0);
    TEST_ASSERT(mcs.NumBonds == 0);
    TEST_ASSERT(mcs.SmartsString == "");
  }
}

void testGitHub6773() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test GitHub issue #6773" << std::endl;
  {
    // these are 2- and 1-hydroxydecalin
    std::vector<ROMOL_SPTR> mols{"OC1CCC2CCCCC2C1"_smiles,
                                 "OC1CCCC2CCCCC12"_smiles};
    {
      MCSParameters p;
      p.BondCompareParameters.RingMatchesRingOnly = true;
      p.AtomCompareParameters.RingMatchesRingOnly = true;
      p.BondCompareParameters.MatchFusedRingsStrict = true;
      auto mcs = findMCS(mols, &p);
      TEST_ASSERT(mcs.NumAtoms == 11);
      TEST_ASSERT(mcs.NumBonds == 11);
      TEST_ASSERT(
          mcs.SmartsString ==
          "[#8&!R]-&!@[#6;r6,!R1]1-&@[#6;r6,!R1]-&@[#6;r6,!R1]-&@[#6]-&@[#6](-&@[#6;r6,!R1]-&@1)-&@[#6&R]-&@[#6&R]-&@[#6&R]-&@[#6&R]");
    }
    {
      MCSParameters p;
      p.BondCompareParameters.RingMatchesRingOnly = true;
      p.AtomCompareParameters.RingMatchesRingOnly = true;
      p.BondCompareParameters.CompleteRingsOnly = true;
      p.AtomCompareParameters.CompleteRingsOnly = true;
      p.BondCompareParameters.MatchFusedRingsStrict = true;
      auto mcs = findMCS(mols, &p);
      TEST_ASSERT(mcs.NumAtoms == 10);
      TEST_ASSERT(mcs.NumBonds == 11);
      TEST_ASSERT(
          mcs.SmartsString ==
          "[#6;r6,!R1]1-&@[#6;r6,!R1]-&@[#6;r6,!R1]-&@[#6]2-&@[#6](-&@[#6;r6,!R1]-&@1)-&@[#6;r6,!R1]-&@[#6;r6,!R1]-&@[#6;r6,!R1]-&@[#6;r6,!R1]-&@2");
    }
    // this is 1-hydroxyperhydroazulene
    mols.push_back("OC1CCC2CCCCCC12"_smiles);
    {
      MCSParameters p;
      p.BondCompareParameters.RingMatchesRingOnly = true;
      p.AtomCompareParameters.RingMatchesRingOnly = true;
      p.BondCompareParameters.MatchFusedRingsStrict = true;
      auto mcs = findMCS(mols, &p);
      TEST_ASSERT(mcs.NumAtoms == 10);
      TEST_ASSERT(mcs.NumBonds == 9);
      TEST_ASSERT(
          mcs.SmartsString ==
          "[#8&!R]-&!@[#6&R](-&@[#6&R]-&@[#6&R])-&@[#6&R]-&@[#6&R]-&@[#6&R]-&@[#6&R]-&@[#6&R]-&@[#6&R]");
    }
    {
      MCSParameters p;
      p.BondCompareParameters.RingMatchesRingOnly = true;
      p.AtomCompareParameters.RingMatchesRingOnly = true;
      p.BondCompareParameters.CompleteRingsOnly = true;
      p.AtomCompareParameters.CompleteRingsOnly = true;
      p.BondCompareParameters.MatchFusedRingsStrict = true;
      auto mcs = findMCS(mols, &p);
      TEST_ASSERT(mcs.NumAtoms == 2);
      TEST_ASSERT(mcs.NumBonds == 1);
      TEST_ASSERT(mcs.SmartsString == "[#8&!R]-&!@[#6&R]");
    }
    mols = std::vector<ROMOL_SPTR>{"C1CC2CCC3CC2C(C1)CCCCCC3"_smiles,
                                   "C1CC2CCC3CCCC4CCC(C1)C2C34"_smiles};
    {
      MCSParameters p;
      p.BondCompareParameters.MatchFusedRingsStrict = true;
      auto mcs = findMCS(mols, &p);
      TEST_ASSERT(mcs.NumAtoms == 15);
      TEST_ASSERT(mcs.NumBonds == 16);
      TEST_ASSERT(
          mcs.SmartsString ==
          "[#6;r6,!R1]1-&@[#6;r6,!R1]-&@[#6]2-&@[#6;r6,!R1]-&@[#6;r6,!R1]-&@[#6](-&@[#6]-&@[#6]-&@2-&@[#6](-&@[#6;r6,!R1]-&@1)-&@[#6]-&@[#6])-&@[#6]-&@[#6]-&@[#6]");
    }
    {
      MCSParameters p;
      p.BondCompareParameters.MatchFusedRingsStrict = true;
      p.BondCompareParameters.CompleteRingsOnly = true;
      auto mcs = findMCS(mols, &p);
      TEST_ASSERT(mcs.NumAtoms == 10);
      TEST_ASSERT(mcs.NumBonds == 11);
      TEST_ASSERT(
          mcs.SmartsString ==
          "[#6;r6,!R1]1-&@[#6;r6,!R1]-&@[#6]2-&@[#6;r6,!R1]-&@[#6;r6,!R1]-&@[#6]-&@[#6]-&@[#6]-&@2-&@[#6]-&@[#6;r6,!R1]-&@1");
      p.AtomCompareParameters.CompleteRingsOnly = true;
      mcs = findMCS(mols, &p);
      TEST_ASSERT(mcs.NumAtoms == 10);
      TEST_ASSERT(mcs.NumBonds == 11);
      TEST_ASSERT(
          mcs.SmartsString ==
          "[#6;r6,!R1]1-&@[#6;r6,!R1]-&@[#6]2-&@[#6;r6,!R1]-&@[#6;r6,!R1]-&@[#6]-&@[#6]-&@[#6]-&@2-&@[#6]-&@[#6;r6,!R1]-&@1");
    }
  }
}

void testBondCompareCompleteRingsOnly() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog)
      << "test CompleteRingsOnly should work also on unsubstituted fused systems"
      << std::endl;
  std::vector<ROMOL_SPTR> mols{"C1CCCC2CCCCC12"_smiles, "C1CCC2CCCC12"_smiles};
  MCSParameters p;
  p.BondCompareParameters.CompleteRingsOnly = true;
  auto mcs = findMCS(mols, &p);
  TEST_ASSERT(mcs.NumAtoms == 0);
  TEST_ASSERT(mcs.NumBonds == 0);
}

void testAtomRingQueries() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test atom ring queries to enforce strict ring fusion"
                       << std::endl;
  std::vector<ROMOL_SPTR> mols{"NC(=O)CC1C2CC3CCCCC3C12"_smiles,
                               "NC(=O)CC1CCCC2CCCCC12"_smiles};
  MCSParameters p;
  p.BondCompareParameters.CompleteRingsOnly = true;
  auto mcs = findMCS(mols, &p);
  TEST_ASSERT(mcs.NumAtoms == 14);
  TEST_ASSERT(mcs.NumBonds == 14);
  TEST_ASSERT(
      mcs.SmartsString ==
      "[#7]-&!@[#6](=&!@[#8])-&!@[#6]-&!@[#6]1-&@[#6]-&@[#6]-&@[#6]-&@[#6]-&@[#6]-&@[#6]-&@[#6]-&@[#6]-&@[#6]-&@1");
  p.BondCompareParameters.MatchFusedRings = true;
  mcs = findMCS(mols, &p);
  TEST_ASSERT(mcs.NumAtoms == 10);
  TEST_ASSERT(mcs.NumBonds == 11);
  TEST_ASSERT(
      mcs.SmartsString ==
      "[#6]1-&@[#6]-&@[#6]-&@[#6]-&@[#6]2-&@[#6]-&@1-&@[#6]-&@[#6]-&@[#6]-&@[#6]-&@2");
  p.BondCompareParameters.MatchFusedRingsStrict = true;
  mcs = findMCS(mols, &p);
  TEST_ASSERT(mcs.NumAtoms == 6);
  TEST_ASSERT(mcs.NumBonds == 6);
  TEST_ASSERT(
      mcs.SmartsString ==
      "[#6;r6,!R1]1-&@[#6;r6,!R1]-&@[#6;r6,!R1]-&@[#6;r6,!R1]-&@[#6]-&@[#6]-&@1");
  ROMOL_SPTR smartsMol(SmartsToMol(mcs.SmartsString));
  std::vector<MatchVectType> matches;
  TEST_ASSERT(SubstructMatch(*mols[0], *mcs.QueryMol, matches) == 1);
  TEST_ASSERT(SubstructMatch(*mols[1], *mcs.QueryMol, matches) == 2);
  TEST_ASSERT(SubstructMatch(*mols[0], *smartsMol, matches) == 1);
  TEST_ASSERT(SubstructMatch(*mols[1], *smartsMol, matches) == 2);
}

//====================================================================================================
//====================================================================================================

int main(int argc, const char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  BOOST_LOG(rdInfoLog) << "FMCS Unit Test \n";

// use maximum CPU resources to increase time measuring accuracy and stability
// in multi process environment
#ifdef WIN32
  //    SetPriorityClass (GetCurrentProcess(), REALTIME_PRIORITY_CLASS );
  SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_HIGHEST);
#else
  setpriority(PRIO_PROCESS, getpid(), -20);
#endif

  T0 = nanoClock();
  t0 = nanoClock();

  testJSONParameters();

  test1Basics();

  test32();
  test190();
  test3();

  testSimpleFast();
  testSimple();
  testSegFault();
  testThreshold();
  testRing1();

  testAtomCompareIsotopes();
  testAtomCompareAnyAtom();
  testAtomCompareAnyAtomBond();
  testAtomCompareAnyHeavyAtom();
  testAtomCompareAnyHeavyAtom1();

  testJnk1LigandsDistance();

  test18();
  test504();

  testInitialSeed();
  testInitialSeed2();

  // chirality check:
  testGithubIssue481();
  testChirality();
  testGithub631();
  //---

  testFormalChargeMatch();
  testGithub2034();
  testGithub945();
  testGithub2420();
  testGithub2663();
  testGithub2662();
  testNaphthalenes();
  testBicycles();
  testBicyclesTricycles();
  test_p38();
  testGithub2714();
  testGitHub2731_comment546175466();
  testQueryMolVsSmarts();
  testCompareNonExistent();
  testGitHub3095();
  testGitHub3458();
  testGitHub3693();
  testGitHub3886();
  testAtomCompareCompleteRingsOnly();
  testGitHub4498();
  testBondStereo();
  testSlowCompleteRingsOnly();
  testCustomShouldAcceptMCS();
  testGitHub5510();
  testDegenerateMCS();
  testGitHub5411();
  testGitHub5440();
  testGitHub5457();
  testGitHub6082();
  testGitHub3965();
  testGitHub6578();
  testGitHub6773();
  testBondCompareCompleteRingsOnly();
  testAtomRingQueries();

  unsigned long long t1 = nanoClock();
  double sec = double(t1 - T0) / 1000000.;
  printf("TOTAL Time elapsed %.2lf seconds\n", sec);

  BOOST_LOG(rdInfoLog)
      << "*******************************************************\n";
  return 0;
}
