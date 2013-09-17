#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <map>
using namespace std;

#include "TTree.h"
#include "TFile.h"
#include "TDirectoryFile.h"

#include "SetStyle.h"
#include "PlotHelper2.h"

#include "TemplateClass.h"
#include "LeptSfProvider.h"
#include "XSecProvider.h"
#include "FakeRateCalculator.h"
#include "pu.h"
#include "scales2.h"

class Index
{
public:
   unsigned int run;
   unsigned int lumi;
   unsigned int event;
public:
   Index() {}
   bool operator <(const Index &other) const
   {
      if(run < other.run)   return true;
      if(run > other.run)   return false;
      if(event < other.event)   return true;
      if(event > other.event)   return false;
      if(lumi < other.lumi)   return true;
      if(lumi > other.lumi)   return false;
      return false;
   }
};

int main(int argc, char *argv[])
{
   SetStyle();

   if(argc != 3)
   {
      cerr << "Usage: " << argv[0] << " ID Energy" << endl;
      cerr << "ID == 0 means running ZX" << endl;
      cerr << "ID < 0 means running on data" << endl;
      cerr << "Energy is the collision energy, can be either 7 or 8" << endl;
      return -1;
   }

   bool Do8TeV = true;
   if(argv[2][0] == '7')
      Do8TeV = false;

   init(!Do8TeV);   // false = 8 TeV, true = 7 TeV

   double Lumi = 1;   // in fb^(-1)
   int ID = atoi(argv[1]);

   cout << "Starting job with ID = " << ID << endl;

   bool DoZX = false;
   bool DoData = false;

   char IDChar[999] = "";
   sprintf(IDChar, "%d", ID);
   string IDString = IDChar;
   string FileName = "hzzTree_id" + IDString + ".root";

   string CMString = "8";
   if(Do8TeV == false)
      CMString = "7";
   
   string BaseDir = "8TeV";
   if(Do8TeV == false)
      BaseDir = "7TeV";
   FileName = BaseDir + "/" + FileName;

   string DirectoryName = "zz4lTree";
   string TreeName = "probe_tree";

   string DataFileName = "hcp8.root";
   if(Do8TeV == false)
      DataFileName = "hcp7.root";

   if(ID < 0)
   {
      DoZX = false;
      DoData = true;
      DirectoryName = "zz4lTree";
      ID = -1;
      FileName = DataFileName;
      IDChar[0] = '\0';
      IDString = "DATA";
   }
   else if(ID == 0)
   {
      DoZX = true;
      DoData = false;
      DirectoryName = "zxTree";
      ID = 0;
      FileName = DataFileName;
      IDChar[0] = '\0';
      IDString = "ZX";
   }

   TFile F(FileName.c_str());
   TDirectoryFile *D = (TDirectoryFile *)F.Get(DirectoryName.c_str());
   TTree *Tree = (TTree *)D->Get(TreeName.c_str());

   FakeRateCalculator FR(DataFileName, true, 40, 120, 0, 0, true);   // whatever the parameters mean...

   double CrossSection = 1;
   if(Do8TeV == true && ID > 0)
      CrossSection = xsecProvider.get8TeVBkgXsec(ID);
   if(Do8TeV == false && ID > 0)
      CrossSection = xsecProvider.get7TeVBkgXsec(ID);

   if(ID < 0)
      CrossSection = 1;
   else if(ID > 1000 && ID <= 1999 && Do8TeV == false)   // ggH
      CrossSection = xsecProvider.get7TeVggHXsec(ID - 1000);
   else if(ID > 1000 && ID <= 1999 && Do8TeV == true)   // ggH
      CrossSection = xsecProvider.get8TeVggHXsec(ID - 1000);
   else if(ID > 2000 && ID <= 2999 && Do8TeV == false)   // VBF
      CrossSection = xsecProvider.get7TeVVBFXsec(ID - 2000);
   else if(ID > 2000 && ID <= 2999 && Do8TeV == true)   // VBF
      CrossSection = xsecProvider.get8TeVVBFXsec(ID - 2000);
   else if(ID > 3000 && ID <= 3999 && Do8TeV == false)   // VH
      CrossSection = xsecProvider.get7TeVWHiXsec(ID - 3000);
   else if(ID > 3000 && ID <= 3999 && Do8TeV == true)   // VH
      CrossSection = xsecProvider.get8TeVWHiXsec(ID - 3000);
   else if(ID > 8000 && ID <= 8999 && Do8TeV == true)   // Pseudoscalar
     CrossSection = xsecProvider.get8TeVggHXsec(ID - 8000) + xsecProvider.get8TeVVBFXsec(ID - 8000) + xsecProvider.get8TeVWHiXsec(ID - 8000) ;
   else if(ID > 9000 && ID <= 9999 && Do8TeV == true)   // SpinTwoMinimal
      CrossSection = xsecProvider.get8TeVggHXsec(ID - 9000) + xsecProvider.get8TeVVBFXsec(ID - 9000) + xsecProvider.get8TeVWHiXsec(ID - 9000) ;

   double EventCount = 1;
   if(Do8TeV == true && ID > 0)
      EventCount = evt_8TeV(ID);
   if(Do8TeV == false && ID > 0 && ID != 1125)
      EventCount = evt_7TeV(ID);
   if(Do8TeV == false && ID > 0 && ID == 1125)
      EventCount = 295356;   // roughly...

   double BaseWeight = CrossSection / EventCount;
   cout << BaseWeight << " " << CrossSection << " " << EventCount << endl;

   TemplateClass M;
   M.Init(Tree);

   TFile OutputFile(("Plots" + CMString + "TeV_" + IDString + ".root").c_str(), "recreate");
   TH1D HMass1("HMass1", "4l Mass;ZZ Mass (GeV);Event", 25, 80, 180);
   TH1D HMass2("HMass2", "4l Mass;ZZ Mass (GeV);Event", 25, 80, 280);
   TH1D HMass3("HMass3", "4l Mass;ZZ Mass (GeV);Event", 27, 100, 181);
   TH1D HZ1Mass("HZ1Mass", "Z1 Mass;Z1 Mass (GeV);Event", 25, 30, 100);
   TH1D HZ1MassRestricted1("HZ1MassRestricted1", "Z1 Mass (M4l = 110-140);Z1 Mass (GeV);Event", 25, 30, 100);
   TH1D HZ1MassRestricted2("HZ1MassRestricted2", "Z1 Mass (M4l = 120-130);Z1 Mass (GeV);Event", 25, 30, 100);
   TH1D HZ1MassRestricted3("HZ1MassRestricted3", "Z1 Mass (M4l = 110-125);Z1 Mass (GeV);Event", 25, 30, 100);
   TH1D HZ1MassRestricted4("HZ1MassRestricted4", "Z1 Mass (M4l = 125-140);Z1 Mass (GeV);Event", 25, 30, 100);
   TH1D HZ2Mass("HZ2Mass", "Z2 Mass;Z2 Mass (GeV);Event", 25, 10, 100);
   TH1D HZ2MassRestricted1("HZ2MassRestricted1", "Z2 Mass (M4l = 110-140);Z2 Mass (GeV);Event", 25, 10, 100);
   TH1D HZ2MassRestricted2("HZ2MassRestricted2", "Z2 Mass (M4l = 120-130);Z2 Mass (GeV);Event", 25, 10, 100);
   TH1D HPT1("HPT1", "PT of the 4l system;PT (GeV);Event", 50, 0, 160);
   TH1D HPT1Restricted1("HPT1Restricted1", "PT of the 4l system (M4l = 110-140);PT (GeV);Event", 50, 0, 160);
   TH1D HPT1Restricted2("HPT1Restricted2", "PT of the 4l system (M4l = 120-130);PT (GeV);Event", 50, 0, 160);
   TH1D HPT2("HPT2", "PT of the 4l system;PT (GeV);Event", 8, 0, 160);
   TH1D HPT2Restricted1("HPT2Restricted1", "PT of the 4l system (M4l = 110-140);PT (GeV);Event", 8, 0, 160);
   TH1D HPT2Restricted2("HPT2Restricted2", "PT of the 4l system (M4l = 120-130);PT (GeV);Event", 8, 0, 160);
   TH1D HRapidity("HRapidity", "Rapidity of the 4l system;Rapidity;Event", 50, -3, 3);
   TH1D HRapidityRestricted1("HRapidityRestricted1", "Rapidity of the 4l system (M4l = 110-140);Rapidity;Event", 50, -3, 3);
   TH1D HRapidityRestricted2("HRapidityRestricted2", "Rapidity of the 4l system (M4l = 120-130);Rapidity;Event", 50, -3, 3);
   TH1D HMELA("HMELA", "MELA;MELA;Event", 50, 0, 1);
   TH1D HMELARestricted1("HMELARestricted1", "MELA (M4l = 110-140);MELA;Event", 50, 0, 1);
   TH1D HMELARestricted2("HMELARestricted2", "MELA (M4l = 120-130);MELA;Event", 50, 0, 1);
   TH1D HPSMELA("HPSMELA", "PSMELA;PSMELA;Event", 50, 0, 1);
   TH1D HPSMELARestricted1("HPSMELARestricted1", "PSMELA (M4l = 110-140);PSMELA;Event", 50, 0, 1);
   TH1D HPSMELARestricted2("HPSMELARestricted2", "PSMELA (M4l = 120-130);PSMELA;Event", 50, 0, 1);

   TH1D HPSMELA_AfterMELACut_Mass110To140("HPSMELA_AfterMELACut_Mass110To140", "PSMELA (MELA > 0.5 && M4l = 110-140);PSMELA;Event", 20, 0, 1);
   TH1D HPSMELA_AfterMELACut_Mass120To130("HPSMELA_AfterMELACut_Mass120To130", "PSMELA (MELA > 0.5 && M4l = 120-130);PSMELA;Event", 20, 0, 1);

   const int nbinsX=21;
   float binsX[nbinsX+1]={0.000, 0.030, 0.060, 0.100, 0.200, 0.300, 0.400, 0.500, 0.550, 0.600, 
      0.633, 0.666, 0.700, 0.733, 0.766, 0.800, 0.833, 0.866, 0.900, 0.933,
      0.966, 1.000};
   const int nbinsY=25;
   float binsY[nbinsY+1]={0.000, 0.100, 0.150, 0.200, 0.233, 0.266, 0.300, 0.333, 0.366, 0.400, 
      0.433, 0.466, 0.500, 0.533, 0.566, 0.600, 0.633, 0.666, 0.700, 0.733, 
      0.766, 0.800, 0.850, 0.900, 0.950, 1.000};

   TH1D HSuperMELA_Mass105To140("HSuperMELA_Mass105To140", "SuperMELA, mass = 105-140;SuperMELA;Event", 50, 0, 1);
   TH2D HSuperMELAVsPSMELA_Mass105To140("HSuperMELAVsPSMELA_Mass105To140",
      "SuperMELA vs. psMELA, mass = 105-140;SuperMELA;psMELA", nbinsX, binsX, nbinsY, binsY);
   TH2D HSuperMELAVsPSMELA_Mass105To140_ee("HSuperMELAVsPSMELA_Mass105To140_ee",
      "SuperMELA vs. psMELA (ee), mass = 105-140;SuperMELA;psMELA", nbinsX, binsX, nbinsY, binsY);
   TH2D HSuperMELAVsPSMELA_Mass105To140_em("HSuperMELAVsPSMELA_Mass105To140_em",
      "SuperMELA vs. psMELA (emu), mass = 105-140;SuperMELA;psMELA", nbinsX, binsX, nbinsY, binsY);
   TH2D HSuperMELAVsPSMELA_Mass105To140_mm("HSuperMELAVsPSMELA_Mass105To140_mm",
      "SuperMELA vs. psMELA (mumu), mass = 105-140;SuperMELA;psMELA", nbinsX, binsX, nbinsY, binsY);

   map<Index, bool> AppearedEvents;

   int EntryCount = Tree->GetEntries();
   for(int iEntry = 0; iEntry < EntryCount; iEntry++)
   {
      Tree->GetEntry(iEntry);

      if(DoData == true || DoZX == true)
      {
         Index newindex;
         newindex.run = M.run;
         newindex.lumi = M.lumi;
         newindex.event = M.event;
         if(AppearedEvents.find(newindex) != AppearedEvents.end())
            continue;
         AppearedEvents.insert(pair<Index, bool>(newindex, true));
      }

      double Weight = 1;

      if(DoZX == false && DoData == false)
      {
         if(Do8TeV == false)
            Weight = BaseWeight * getPUWeight2011((int)M.numTrueInteractions);
         else
            Weight = BaseWeight * getPUWeight2012((int)M.numTrueInteractions, 1);
         Weight = Weight * getSF(M.l1pt, M.l1eta, M.l1pdgId);
         Weight = Weight * getSF(M.l2pt, M.l2eta, M.l2pdgId);
         Weight = Weight * getSF(M.l3pt, M.l3eta, M.l3pdgId);
         Weight = Weight * getSF(M.l4pt, M.l4eta, M.l4pdgId);
      }
      else if(DoData == false)
      {
         double F1 = FR.getFakeRate(M.l3pt, M.l3eta, M.l3pdgId);
         double F2 = FR.getFakeRate(M.l4pt, M.l4eta, M.l4pdgId);
         double P1 = getPR(M.l3pt, M.l3eta, M.l3pdgId);
         double P2 = getPR(M.l4pt, M.l4eta, M.l4pdgId);

         double Eps1 = F1 / (1 - F1);
         double Eps2 = F2 / (1 - F2);
         double Eta1 = (1 - P1) / P1;
         double Eta2 = (1 - P2) / P2;
         double Denominator = (1 - (Eps1 * Eta1)) * (1 - (Eps2 * Eta2));
         
         if(M.l3pdgId == M.l4pdgId)
         {
            double ScaleFactor = 1;
            if(M.channel == 0)
               ScaleFactor = 1.28;
            else if(M.channel == 1)
               ScaleFactor = 0.93;
            else
               ScaleFactor = 0.94;
            Weight = F1 * F2 * ScaleFactor;
         }
         /*
         else if(M.l3pdgId == -M.l4pdgId && (M.l3idNew == 0 || M.l3bdtIso / M.l3pt > 0.4)
            && (M.l4idNew == 0 || M.l4bdtIso / M.l4pt > 0.4))
            Weight = -Eps1 * Eps2 / Denominator;
         else if(M.l3pdgId == -M.l4pdgId && (M.l3idNew == 0 || M.l3bdtIso / M.l3pt > 0.4)
            && (M.l4idNew == 1 || M.l4bdtIso / M.l4pt <= 0.4))
            Weight = Eps1 / Denominator;
         else if(M.l3pdgId == -M.l4pdgId && (M.l3idNew == 1 || M.l3bdtIso / M.l3pt <= 0.4)
            && (M.l4idNew == 0 || M.l4bdtIso / M.l4pt > 0.4))
            Weight = Eps2 / Denominator;
         */
         else
            Weight = 0;   // Hmm?
      }
      else
         Weight = 1;

      HMass1.Fill(M.mass, Weight);
      HMass2.Fill(M.mass, Weight);
      HMass3.Fill(M.mass, Weight);
      HZ1Mass.Fill(M.z1mass, Weight);
      HZ2Mass.Fill(M.z2mass, Weight);
      HPT1.Fill(M.pt, Weight);
      HPT2.Fill(M.pt, Weight);
      HRapidity.Fill(M.rap, Weight);
      HMELA.Fill(M.melaLD, Weight);
      HPSMELA.Fill(M.melaPSLD, Weight);
      if(M.mass > 105 && M.mass < 140)
      {
         HSuperMELA_Mass105To140.Fill(M.SuperMELA, Weight);
         HSuperMELAVsPSMELA_Mass105To140.Fill(M.SuperMELA, M.melaPSLD, Weight);
         if(M.channel == 1)
            HSuperMELAVsPSMELA_Mass105To140_ee.Fill(M.SuperMELA, M.melaPSLD, Weight);
         if(M.channel == 0)
            HSuperMELAVsPSMELA_Mass105To140_mm.Fill(M.SuperMELA, M.melaPSLD, Weight);
         if(M.channel == 2 || M.channel == 3)
            HSuperMELAVsPSMELA_Mass105To140_em.Fill(M.SuperMELA, M.melaPSLD, Weight);
      }
      if(M.mass > 110 && M.mass < 140)
      {
         HZ1MassRestricted1.Fill(M.z1mass, Weight);
         HZ2MassRestricted1.Fill(M.z2mass, Weight);
         HPT1Restricted1.Fill(M.pt, Weight);
         HPT2Restricted1.Fill(M.pt, Weight);
         HRapidityRestricted1.Fill(M.rap, Weight);
         HMELARestricted1.Fill(M.melaLD, Weight);
         HPSMELARestricted1.Fill(M.melaPSLD, Weight);
         if(M.melaLD > 0.5) {
           HPSMELA_AfterMELACut_Mass110To140.Fill(M.melaPSLD, Weight);
         }
      }
      if(M.mass > 120 && M.mass < 130)
      {
         HZ1MassRestricted2.Fill(M.z1mass, Weight);
         HZ2MassRestricted2.Fill(M.z2mass, Weight);
         HPT1Restricted2.Fill(M.pt, Weight);
         HPT2Restricted2.Fill(M.pt, Weight);
         HRapidityRestricted2.Fill(M.rap, Weight);
         HMELARestricted2.Fill(M.melaLD, Weight);
         HPSMELARestricted2.Fill(M.melaPSLD, Weight);
         if(M.melaLD > 0.5) {
           HPSMELA_AfterMELACut_Mass120To130.Fill(M.melaPSLD, Weight);
         }
      }
      if(M.mass > 110 && M.mass < 125)
         HZ1MassRestricted3.Fill(M.z1mass, Weight);
      if(M.mass > 125 && M.mass < 140)
         HZ1MassRestricted4.Fill(M.z1mass, Weight);
      
   }

   HSuperMELA_Mass105To140.SetStats(0);
   HSuperMELAVsPSMELA_Mass105To140.SetStats(0);

   PsFileHelper PsFile(("Plots" + CMString + "TeV_" + IDString + ".ps").c_str());
   PsFile.AddTextPage(IDString);

   PsFile.AddPlot(HMass1, "", false);
   PsFile.AddPlot(HMass1, "", true);
   PsFile.AddPlot(HMass2, "", false);
   PsFile.AddPlot(HMass2, "", true);
   PsFile.AddPlot(HMass3, "", false);
   PsFile.AddPlot(HMass3, "", true);
   PsFile.AddPlot(HZ1Mass, "", false);
   PsFile.AddPlot(HZ1Mass, "", true);
   PsFile.AddPlot(HZ1MassRestricted1, "", false);
   PsFile.AddPlot(HZ1MassRestricted1, "", true);
   PsFile.AddPlot(HZ1MassRestricted2, "", false);
   PsFile.AddPlot(HZ1MassRestricted2, "", true);
   PsFile.AddPlot(HZ1MassRestricted3, "", false);
   PsFile.AddPlot(HZ1MassRestricted3, "", true);
   PsFile.AddPlot(HZ1MassRestricted4, "", false);
   PsFile.AddPlot(HZ1MassRestricted4, "", true);
   PsFile.AddPlot(HZ2Mass, "", false);
   PsFile.AddPlot(HZ2Mass, "", true);
   PsFile.AddPlot(HZ2MassRestricted1, "", false);
   PsFile.AddPlot(HZ2MassRestricted1, "", true);
   PsFile.AddPlot(HZ2MassRestricted2, "", false);
   PsFile.AddPlot(HZ2MassRestricted2, "", true);
   PsFile.AddPlot(HPT1, "", false);
   PsFile.AddPlot(HPT1, "", true);
   PsFile.AddPlot(HPT1Restricted1, "", false);
   PsFile.AddPlot(HPT1Restricted1, "", true);
   PsFile.AddPlot(HPT1Restricted2, "", false);
   PsFile.AddPlot(HPT1Restricted2, "", true);
   PsFile.AddPlot(HPT2, "", false);
   PsFile.AddPlot(HPT2, "", true);
   PsFile.AddPlot(HPT2Restricted1, "", false);
   PsFile.AddPlot(HPT2Restricted1, "", true);
   PsFile.AddPlot(HPT2Restricted2, "", false);
   PsFile.AddPlot(HPT2Restricted2, "", true);
   PsFile.AddPlot(HRapidity, "", false);
   PsFile.AddPlot(HRapidity, "", true);
   PsFile.AddPlot(HRapidityRestricted1, "", false);
   PsFile.AddPlot(HRapidityRestricted1, "", true);
   PsFile.AddPlot(HRapidityRestricted2, "", false);
   PsFile.AddPlot(HRapidityRestricted2, "", true);
   PsFile.AddPlot(HMELA, "", false);
   PsFile.AddPlot(HMELA, "", true);
   PsFile.AddPlot(HMELARestricted1, "", false);
   PsFile.AddPlot(HMELARestricted1, "", true);
   PsFile.AddPlot(HMELARestricted2, "", false);
   PsFile.AddPlot(HMELARestricted2, "", true);
   PsFile.AddPlot(HPSMELA, "", false);
   PsFile.AddPlot(HPSMELA, "", true);
   PsFile.AddPlot(HPSMELARestricted1, "", false);
   PsFile.AddPlot(HPSMELARestricted1, "", true);
   PsFile.AddPlot(HPSMELARestricted2, "", false);
   PsFile.AddPlot(HPSMELARestricted2, "", true);
   PsFile.AddPlot(HPSMELA_AfterMELACut_Mass110To140, "", false);
   PsFile.AddPlot(HPSMELA_AfterMELACut_Mass110To140, "", true);
   PsFile.AddPlot(HPSMELA_AfterMELACut_Mass120To130, "", false);
   PsFile.AddPlot(HPSMELA_AfterMELACut_Mass120To130, "", true);
   PsFile.AddPlot(HSuperMELA_Mass105To140, "", false);
   PsFile.AddPlot(HSuperMELA_Mass105To140, "", true);
   
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140, "colz", false, false);
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140, "colz", false, true);
   TH2D *HSuperMELAVsPSMELA_Mass105To140_Smoothed = (TH2D *)HSuperMELAVsPSMELA_Mass105To140.Clone("HSuperMELAVsPSMELA_Mass105To140_Smoothed");
   HSuperMELAVsPSMELA_Mass105To140_Smoothed->Smooth(1, "k5b");
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140_Smoothed, "colz", false, false);
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140_Smoothed, "colz", false, true);
   
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140_ee, "colz", false, false);
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140_ee, "colz", false, true);
   TH2D *HSuperMELAVsPSMELA_Mass105To140_ee_Smoothed = (TH2D *)HSuperMELAVsPSMELA_Mass105To140.Clone("HSuperMELAVsPSMELA_Mass105To140_ee_Smoothed");
   HSuperMELAVsPSMELA_Mass105To140_ee_Smoothed->Smooth(1, "k5b");
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140_ee_Smoothed, "colz", false, false);
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140_ee_Smoothed, "colz", false, true);
   
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140_em, "colz", false, false);
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140_em, "colz", false, true);
   TH2D *HSuperMELAVsPSMELA_Mass105To140_em_Smoothed = (TH2D *)HSuperMELAVsPSMELA_Mass105To140_em.Clone("HSuperMELAVsPSMELA_Mass105To140_em_Smoothed");
   HSuperMELAVsPSMELA_Mass105To140_em_Smoothed->Smooth(1, "k5b");
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140_em_Smoothed, "colz", false, false);
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140_em_Smoothed, "colz", false, true);
   
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140_mm, "colz", false, false);
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140_mm, "colz", false, true);
   TH2D *HSuperMELAVsPSMELA_Mass105To140_mm_Smoothed = (TH2D *)HSuperMELAVsPSMELA_Mass105To140_mm.Clone("HSuperMELAVsPSMELA_Mass105To140_mm_Smoothed");
   HSuperMELAVsPSMELA_Mass105To140_mm_Smoothed->Smooth(1, "k5b");
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140_mm_Smoothed, "colz", false, false);
   PsFile.AddPlot(HSuperMELAVsPSMELA_Mass105To140_mm_Smoothed, "colz", false, true);
   
   PsFile.AddTimeStampPage();
   PsFile.Close();

   HMass1.Write();
   HMass2.Write();
   HMass3.Write();
   HZ1Mass.Write();
   HZ1MassRestricted1.Write();
   HZ1MassRestricted2.Write();
   HZ1MassRestricted3.Write();
   HZ1MassRestricted4.Write();
   HZ2Mass.Write();
   HZ2MassRestricted1.Write();
   HZ2MassRestricted2.Write();
   HPT1.Write();
   HPT1Restricted1.Write();
   HPT1Restricted2.Write();
   HPT2.Write();
   HPT2Restricted1.Write();
   HPT2Restricted2.Write();
   HRapidity.Write();
   HRapidityRestricted1.Write();
   HRapidityRestricted2.Write();
   HMELA.Write();
   HMELARestricted1.Write();
   HMELARestricted2.Write();
   HPSMELA.Write();
   HPSMELARestricted1.Write();
   HPSMELARestricted2.Write();
   HPSMELA_AfterMELACut_Mass110To140.Write();
   HPSMELA_AfterMELACut_Mass120To130.Write();
   HSuperMELA_Mass105To140.Write();
   HSuperMELAVsPSMELA_Mass105To140.Write();
   HSuperMELAVsPSMELA_Mass105To140_Smoothed->Write();
   HSuperMELAVsPSMELA_Mass105To140_ee.Write();
   HSuperMELAVsPSMELA_Mass105To140_ee_Smoothed->Write();
   HSuperMELAVsPSMELA_Mass105To140_em.Write();
   HSuperMELAVsPSMELA_Mass105To140_em_Smoothed->Write();
   HSuperMELAVsPSMELA_Mass105To140_mm.Write();
   HSuperMELAVsPSMELA_Mass105To140_mm_Smoothed->Write();

   OutputFile.Close();

   F.Close();

   return 0;
}





