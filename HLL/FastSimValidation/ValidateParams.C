//================================================================================================
//					ValidateParams
// 
// Reads kinematic parameters from fullsim and fastsim ntuples, plots fullsim and fastsim for comparison
//
// USAGE
//
// ValidateParams(string fullsimFilename, string fastsimFilename, string outputDirectory, Int_t norm_option = 1,  Int_t error_option = 2, const string Label = "ZZ")
// 
// INPUT
// fullsimFilename and fastsimFilename		.root files containing the ntuples with fullsim and fastsim
// outputDirectory				directory to which all histograms will be outputted.
// norm_option					0 = do not normalize, 1 = normalize (DEFAULT)
// error_option 				0 = no error bars, 1 = on fastsim, 2 = on fullsim (DEFAULT), 3 = on both
// Label					label to be attached to the output files' name
//
// OUTPUT
// write .gif and .root files to the specified directory.
//________________________________________________________________________________________________


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TStyle.h>		    // access to gStyle
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <iostream>                 // standard I/O
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"

#include "CITHZZ/CommonCode/HZZEventTree.h"
#include "CITHZZ/CommonCode/CommonDefs.hh"

#include "TLorentzVector.h"
#include "CITHZZ/UtilityFunctions/AngleConversion.h"

#endif


//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}

//*************************************************************************************************
//Plotting function
//*************************************************************************************************
void PlotTwo(TH1* fullsim, TH1* fastsim, const char* outputFile, Int_t error_option = 2) {
  // error_option sets whether to plot error bars or not
  // 0	no error bars
  // 1	error bars on fastsim
  // 2	error bars on fullsim 
  // 3	error bars EVERYWHERE

  TCanvas* cv = new TCanvas("cv","cv", 800,600);
  TLegend* legend = new TLegend(0.7,0.8,0.9,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(fastsim,"FastSim", "LP");
  legend->AddEntry(fullsim,"FullSim", "LP");
  fastsim->SetLineColor(kBlue);
  fullsim->SetLineColor(kRed);
  if (error_option == 0) {
    fastsim->Draw("hist");
    fullsim->Draw("histsame");
  }
  if (error_option == 1) {
    fastsim->Draw("e1");
    fullsim->Draw("histsame");
  }
  if (error_option == 2) {
    fastsim->Draw("hist");
    fullsim->Draw("e1same");
  }
  if (error_option == 3) {
    fastsim->Draw("e1");
    fullsim->Draw("e1same");
  }

  legend->Draw();
  cv->SaveAs(outputFile);

  // Cleaning up
  delete cv;
  delete legend;
}


//*************************************************************************************************
//Main Function
//*************************************************************************************************
void ValidateParams(string fullsimFilename, string fastsimFilename, string outputDirectory, Int_t norm_option = 1,  Int_t error_option = 2, const string Label = "ZZ") {

  // Setting graphic style
  gStyle->SetOptStat(0);

  string label = Label;
  if (Label != "") label = "_" + Label;


  //*****************************************************************************************
  // Load Input
  //*****************************************************************************************
  HZZEventTree fullsimHZZEventTree;
  fullsimHZZEventTree.LoadTree(fullsimFilename.c_str());
  fullsimHZZEventTree.InitTree();
  HZZEventTree fastsimHZZEventTree;
  fastsimHZZEventTree.LoadTree(fastsimFilename.c_str());
  fastsimHZZEventTree.InitTree();
  cout << "Events in the ntuple: " << fastsimHZZEventTree.tree_->GetEntries() << endl;




  //*************************************************************************************************
  //Histograms
  //*************************************************************************************************

  // Saving all the histograms to file
  TFile outRootFile((outputDirectory + "ValidateParams_AllHistograms" + label + ".root").c_str(), "recreate");

  TH1F *histZ1Mass_fullsim_ee = new TH1F( "histZ1Mass_fullsim_ee", ";Z1 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ1Mass_fastsim_ee = new TH1F( "histZ1Mass_fastsim_ee", ";Z1 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ1Mass_fullsim_mm = new TH1F( "histZ1Mass_fullsim_mm", ";Z1 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ1Mass_fastsim_mm = new TH1F( "histZ1Mass_fastsim_mm", ";Z1 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);

  TH1F *histZ2Mass_fullsim_ee = new TH1F( "histZ2Mass_fullsim_ee", ";Z2 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ2Mass_fastsim_ee = new TH1F( "histZ2Mass_fastsim_ee", ";Z2 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ2Mass_fullsim_mm = new TH1F( "histZ2Mass_fullsim_mm", ";Z2 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);
  TH1F *histZ2Mass_fastsim_mm = new TH1F( "histZ2Mass_fastsim_mm", ";Z2 Mass [GeV/c^{2}]; Number of Events", 100, 0,200);

  TH1F *histM4l_fullsim_4e = new TH1F( "histM4l_fullsim_4e", ";m4l [GeV/c^{2}]; Number of Events", 100, 0,300);
  TH1F *histM4l_fastsim_4e = new TH1F( "histM4l_fastsim_4e", ";m4l [GeV/c^{2}]; Number of Events", 100, 0,300);
  TH1F *histM4l_fullsim_4m = new TH1F( "histM4l_fullsim_4m", ";m4l [GeV/c^{2}]; Number of Events", 100, 0,300);
  TH1F *histM4l_fastsim_4m = new TH1F( "histM4l_fastsim_4m", ";m4l [GeV/c^{2}]; Number of Events", 100, 0,300);
  TH1F *histM4l_fullsim_2e2m = new TH1F( "histM4l_fullsim_2e2m", ";m4l [GeV/c^{2}]; Number of Events", 100, 0,300);
  TH1F *histM4l_fastsim_2e2m = new TH1F( "histM4l_fastsim_2e2m", ";m4l [GeV/c^{2}]; Number of Events", 100, 0,300);

  TH1F *histLep1Eta_fullsim_e = new TH1F("histLep1Eta_fullsim_e", ";Lepton 1 eta; Number of Events", 50, -4, 4);
  TH1F *histLep1Eta_fastsim_e = new TH1F("histLep1Eta_fastsim_e", ";Lepton 1 eta; Number of Events", 50, -4, 4);
  TH1F *histLep1Phi_fullsim_e = new TH1F("histLep1Phi_fullsim_e", ";Lepton 1 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep1Phi_fastsim_e = new TH1F("histLep1Phi_fastsim_e", ";Lepton 1 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep1Pt_fullsim_e = new TH1F("histLep1Pt_fullsim_e", ";Lepton 1 Pt; Number of Events", 150, 0, 300);
  TH1F *histLep1Pt_fastsim_e = new TH1F("histLep1Pt_fastsim_e", ";Lepton 1 Pt; Number of Events", 150, 0, 300);

  TH1F *histLep1Eta_fullsim_m = new TH1F("histLep1Eta_fullsim_m", ";Lepton 1 eta; Number of Events", 50, -4, 4);
  TH1F *histLep1Eta_fastsim_m = new TH1F("histLep1Eta_fastsim_m", ";Lepton 1 eta; Number of Events", 50, -4, 4);
  TH1F *histLep1Phi_fullsim_m = new TH1F("histLep1Phi_fullsim_m", ";Lepton 1 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep1Phi_fastsim_m = new TH1F("histLep1Phi_fastsim_m", ";Lepton 1 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep1Pt_fullsim_m = new TH1F("histLep1Pt_fullsim_m", ";Lepton 1 Pt; Number of Events", 150, 0, 300);
  TH1F *histLep1Pt_fastsim_m = new TH1F("histLep1Pt_fastsim_m", ";Lepton 1 Pt; Number of Events", 150, 0, 300);

  TH1F *histLep2Eta_fullsim_e = new TH1F("histLep2Eta_fullsim_e", ";Lepton 2 eta; Number of Events", 50, -4, 4);
  TH1F *histLep2Eta_fastsim_e = new TH1F("histLep2Eta_fastsim_e", ";Lepton 2 eta; Number of Events", 50, -4, 4);
  TH1F *histLep2Phi_fullsim_e = new TH1F("histLep2Phi_fullsim_e", ";Lepton 2 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep2Phi_fastsim_e = new TH1F("histLep2Phi_fastsim_e", ";Lepton 2 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep2Pt_fullsim_e = new TH1F("histLep2Pt_fullsim_e", ";Lepton 2 Pt; Number of Events", 150, 0, 300);
  TH1F *histLep2Pt_fastsim_e = new TH1F("histLep2Pt_fastsim_e", ";Lepton 2 Pt; Number of Events", 150, 0, 300);

  TH1F *histLep2Eta_fullsim_m = new TH1F("histLep2Eta_fullsim_m", ";Lepton 2 eta; Number of Events", 50, -4, 4);
  TH1F *histLep2Eta_fastsim_m = new TH1F("histLep2Eta_fastsim_m", ";Lepton 2 eta; Number of Events", 50, -4, 4);
  TH1F *histLep2Phi_fullsim_m = new TH1F("histLep2Phi_fullsim_m", ";Lepton 2 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep2Phi_fastsim_m = new TH1F("histLep2Phi_fastsim_m", ";Lepton 2 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep2Pt_fullsim_m = new TH1F("histLep2Pt_fullsim_m", ";Lepton 2 Pt; Number of Events", 150, 0, 300);
  TH1F *histLep2Pt_fastsim_m = new TH1F("histLep2Pt_fastsim_m", ";Lepton 2 Pt; Number of Events", 150, 0, 300);

  TH1F *histLep3Eta_fullsim_e = new TH1F("histLep3Eta_fullsim_e", ";Lepton 3 eta; Number of Events", 50, -4, 4);
  TH1F *histLep3Eta_fastsim_e = new TH1F("histLep3Eta_fastsim_e", ";Lepton 3 eta; Number of Events", 50, -4, 4);
  TH1F *histLep3Phi_fullsim_e = new TH1F("histLep3Phi_fullsim_e", ";Lepton 3 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep3Phi_fastsim_e = new TH1F("histLep3Phi_fastsim_e", ";Lepton 3 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep3Pt_fullsim_e = new TH1F("histLep3Pt_fullsim_e", ";Lepton 3 Pt; Number of Events", 150, 0, 300);
  TH1F *histLep3Pt_fastsim_e = new TH1F("histLep3Pt_fastsim_e", ";Lepton 3 Pt; Number of Events", 150, 0, 300);

  TH1F *histLep3Eta_fullsim_m = new TH1F("histLep3Eta_fullsim_m", ";Lepton 3 eta; Number of Events", 50, -4, 4);
  TH1F *histLep3Eta_fastsim_m = new TH1F("histLep3Eta_fastsim_m", ";Lepton 3 eta; Number of Events", 50, -4, 4);
  TH1F *histLep3Phi_fullsim_m = new TH1F("histLep3Phi_fullsim_m", ";Lepton 3 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep3Phi_fastsim_m = new TH1F("histLep3Phi_fastsim_m", ";Lepton 3 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep3Pt_fullsim_m = new TH1F("histLep3Pt_fullsim_m", ";Lepton 3 Pt; Number of Events", 150, 0, 300);
  TH1F *histLep3Pt_fastsim_m = new TH1F("histLep3Pt_fastsim_m", ";Lepton 3 Pt; Number of Events", 150, 0, 300);

  TH1F *histLep4Eta_fullsim_e = new TH1F("histLep4Eta_fullsim_e", ";Lepton 4 eta; Number of Events", 50, -4, 4);
  TH1F *histLep4Eta_fastsim_e = new TH1F("histLep4Eta_fastsim_e", ";Lepton 4 eta; Number of Events", 50, -4, 4);
  TH1F *histLep4Phi_fullsim_e = new TH1F("histLep4Phi_fullsim_e", ";Lepton 4 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep4Phi_fastsim_e = new TH1F("histLep4Phi_fastsim_e", ";Lepton 4 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep4Pt_fullsim_e = new TH1F("histLep4Pt_fullsim_e", ";Lepton 4 Pt; Number of Events", 150, 0, 300);
  TH1F *histLep4Pt_fastsim_e = new TH1F("histLep4Pt_fastsim_e", ";Lepton 4 Pt; Number of Events", 150, 0, 300);

  TH1F *histLep4Eta_fullsim_m = new TH1F("histLep4Eta_fullsim_m", ";Lepton 4 eta; Number of Events", 50, -4, 4);
  TH1F *histLep4Eta_fastsim_m = new TH1F("histLep4Eta_fastsim_m", ";Lepton 4 eta; Number of Events", 50, -4, 4);
  TH1F *histLep4Phi_fullsim_m = new TH1F("histLep4Phi_fullsim_m", ";Lepton 4 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep4Phi_fastsim_m = new TH1F("histLep4Phi_fastsim_m", ";Lepton 4 phi; Number of Events", 50, -3.2, 3.2);
  TH1F *histLep4Pt_fullsim_m = new TH1F("histLep4Pt_fullsim_m", ";Lepton 4 Pt; Number of Events", 150, 0, 300);
  TH1F *histLep4Pt_fastsim_m = new TH1F("histLep4Pt_fastsim_m", ";Lepton 4 Pt; Number of Events", 150, 0, 300);

  TH1F *histAllLepsEta_fullsim_e = new TH1F("histAllLepsEta_fullsim_e", ";Leptons eta; Number of Events", 50, -4, 4);
  TH1F *histAllLepsEta_fastsim_e = new TH1F("histAllLepsEta_fastsim_e", ";Leptons eta; Number of Events", 50, -4, 4);
  TH1F *histAllLepsEta_fullsim_m = new TH1F("histAllLepsEta_fullsim_m", ";Leptons eta; Number of Events", 50, -4, 4);
  TH1F *histAllLepsEta_fastsim_m = new TH1F("histAllLepsEta_fastsim_m", ";Leptons eta; Number of Events", 50, -4, 4);

  TH1F *histAllLepsPhi_fullsim_e = new TH1F("histAllLepsPhi_fullsim_e", ";Leptons phi; Number of Events", 50, -4, 4);
  TH1F *histAllLepsPhi_fastsim_e = new TH1F("histAllLepsPhi_fastsim_e", ";Leptons phi; Number of Events", 50, -4, 4);
  TH1F *histAllLepsPhi_fullsim_m = new TH1F("histAllLepsPhi_fullsim_m", ";Leptons phi; Number of Events", 50, -4, 4);
  TH1F *histAllLepsPhi_fastsim_m = new TH1F("histAllLepsPhi_fastsim_m", ";Leptons phi; Number of Events", 50, -4, 4);

  TH1F *histAllLepsPt_fullsim_e = new TH1F("histAllLepsPt_fullsim_e", ";Leptons Pt; Number of Events", 150, 0, 300);
  TH1F *histAllLepsPt_fastsim_e = new TH1F("histAllLepsPt_fastsim_e", ";Leptons Pt; Number of Events", 150, 0, 300);
  TH1F *histAllLepsPt_fullsim_m = new TH1F("histAllLepsPt_fullsim_m", ";Leptons Pt; Number of Events", 150, 0, 300);
  TH1F *histAllLepsPt_fastsim_m = new TH1F("histAllLepsPt_fastsim_m", ";Leptons Pt; Number of Events", 150, 0, 300);

  //*****************************************************************************************
  // Loop over fullsim events
  //*****************************************************************************************
  for(UInt_t ientry=0; ientry < fullsimHZZEventTree.tree_->GetEntries(); ientry++) {       	
    fullsimHZZEventTree.tree_->GetEntry(ientry);
    if (ientry % 1000000 == 0) cout << "Event " << ientry << endl;

    //how to get lepton four momenta
    TLorentzVector lepton1FourVector;
    lepton1FourVector.SetPtEtaPhiM( fullsimHZZEventTree.fLep1Pt, fullsimHZZEventTree.fLep1Eta, fullsimHZZEventTree.fLep1Phi, ELECTRONMASS);


    // Setting up normalization options
    Float_t weight = fullsimHZZEventTree.fWeight; 

    //fill events passing selection

    // Filling Z masses
    if (fullsimHZZEventTree.fPassFullSelection) {
      if (abs(fullsimHZZEventTree.fLep1Type) == 11) {
	histZ1Mass_fullsim_ee->Fill(fullsimHZZEventTree.fZ1Mass, weight);
	histZ2Mass_fullsim_ee->Fill(fullsimHZZEventTree.fZ2Mass, weight);
      } else {
	histZ1Mass_fullsim_mm->Fill(fullsimHZZEventTree.fZ1Mass, weight);
	histZ2Mass_fullsim_mm->Fill(fullsimHZZEventTree.fZ2Mass, weight);
      }

      // Filling m4l masses (in 3 categories: 4e, 4m and 2e2m)
      if (abs(fullsimHZZEventTree.fLep1Type) == 11 && abs(fullsimHZZEventTree.fLep3Type) == 11) {
	histM4l_fullsim_4e->Fill(fullsimHZZEventTree.fFourLeptonMass, weight);
      } else if (abs(fullsimHZZEventTree.fLep1Type) == 13 && abs(fullsimHZZEventTree.fLep3Type) == 13) {
	histM4l_fullsim_4m->Fill(fullsimHZZEventTree.fFourLeptonMass, weight);
      } else {
	histM4l_fullsim_2e2m->Fill(fullsimHZZEventTree.fFourLeptonMass, weight);
      }

      // Filling all the lepton parameters 
      if (abs(fullsimHZZEventTree.fLep1Type) == 11) {
	histLep1Eta_fullsim_e->Fill(fullsimHZZEventTree.fLep1Eta, weight);
	histLep1Phi_fullsim_e->Fill(fullsimHZZEventTree.fLep1Phi, weight);
	histLep1Pt_fullsim_e->Fill(fullsimHZZEventTree.fLep1Pt, weight);
	histAllLepsEta_fullsim_e->Fill(fullsimHZZEventTree.fLep1Eta, weight);
	histAllLepsPhi_fullsim_e->Fill(fullsimHZZEventTree.fLep1Phi, weight);
	histAllLepsPt_fullsim_e->Fill(fullsimHZZEventTree.fLep1Pt, weight);
      } else {
	histLep1Eta_fullsim_m->Fill(fullsimHZZEventTree.fLep1Eta, weight);
	histLep1Phi_fullsim_m->Fill(fullsimHZZEventTree.fLep1Phi, weight);
	histLep1Pt_fullsim_m->Fill(fullsimHZZEventTree.fLep1Pt, weight);
	histAllLepsEta_fullsim_m->Fill(fullsimHZZEventTree.fLep1Eta, weight);
	histAllLepsPhi_fullsim_m->Fill(fullsimHZZEventTree.fLep1Phi, weight);
	histAllLepsPt_fullsim_m->Fill(fullsimHZZEventTree.fLep1Pt, weight);
      }

      if (abs(fullsimHZZEventTree.fLep2Type) == 11) {
	histLep2Eta_fullsim_e->Fill(fullsimHZZEventTree.fLep2Eta, weight);
	histLep2Phi_fullsim_e->Fill(fullsimHZZEventTree.fLep2Phi, weight);
	histLep2Pt_fullsim_e->Fill(fullsimHZZEventTree.fLep2Pt, weight);
	histAllLepsEta_fullsim_e->Fill(fullsimHZZEventTree.fLep2Eta, weight);
	histAllLepsPhi_fullsim_e->Fill(fullsimHZZEventTree.fLep2Phi, weight);
	histAllLepsPt_fullsim_e->Fill(fullsimHZZEventTree.fLep2Pt, weight);
      } else {
	histLep2Eta_fullsim_m->Fill(fullsimHZZEventTree.fLep2Eta, weight);
	histLep2Phi_fullsim_m->Fill(fullsimHZZEventTree.fLep2Phi, weight);
	histLep2Pt_fullsim_m->Fill(fullsimHZZEventTree.fLep2Pt, weight);
	histAllLepsEta_fullsim_m->Fill(fullsimHZZEventTree.fLep2Eta, weight);
	histAllLepsPhi_fullsim_m->Fill(fullsimHZZEventTree.fLep2Phi, weight);
	histAllLepsPt_fullsim_m->Fill(fullsimHZZEventTree.fLep2Pt, weight);
      }

      if (abs(fullsimHZZEventTree.fLep3Type) == 11) {
	histLep3Eta_fullsim_e->Fill(fullsimHZZEventTree.fLep3Eta, weight);
	histLep3Phi_fullsim_e->Fill(fullsimHZZEventTree.fLep3Phi, weight);
	histLep3Pt_fullsim_e->Fill(fullsimHZZEventTree.fLep3Pt, weight);
	histAllLepsEta_fullsim_e->Fill(fullsimHZZEventTree.fLep3Eta, weight);
	histAllLepsPhi_fullsim_e->Fill(fullsimHZZEventTree.fLep3Phi, weight);
	histAllLepsPt_fullsim_e->Fill(fullsimHZZEventTree.fLep3Pt, weight);
      } else {
	histLep3Eta_fullsim_m->Fill(fullsimHZZEventTree.fLep3Eta, weight);
	histLep3Phi_fullsim_m->Fill(fullsimHZZEventTree.fLep3Phi, weight);
	histLep3Pt_fullsim_m->Fill(fullsimHZZEventTree.fLep3Pt, weight);
	histAllLepsEta_fullsim_m->Fill(fullsimHZZEventTree.fLep3Eta, weight);
	histAllLepsPhi_fullsim_m->Fill(fullsimHZZEventTree.fLep3Phi, weight);
	histAllLepsPt_fullsim_m->Fill(fullsimHZZEventTree.fLep3Pt, weight);
      }

      if (abs(fullsimHZZEventTree.fLep4Type) == 11) {
	histLep4Eta_fullsim_e->Fill(fullsimHZZEventTree.fLep4Eta, weight);
	histLep4Phi_fullsim_e->Fill(fullsimHZZEventTree.fLep4Phi, weight);
	histLep4Pt_fullsim_e->Fill(fullsimHZZEventTree.fLep4Pt, weight);
	histAllLepsEta_fullsim_e->Fill(fullsimHZZEventTree.fLep4Eta, weight);
	histAllLepsPhi_fullsim_e->Fill(fullsimHZZEventTree.fLep4Phi, weight);
	histAllLepsPt_fullsim_e->Fill(fullsimHZZEventTree.fLep4Pt, weight);
      } else {
	histLep4Eta_fullsim_m->Fill(fullsimHZZEventTree.fLep4Eta, weight);
	histLep4Phi_fullsim_m->Fill(fullsimHZZEventTree.fLep4Phi, weight);
	histLep4Pt_fullsim_m->Fill(fullsimHZZEventTree.fLep4Pt, weight);
	histAllLepsEta_fullsim_m->Fill(fullsimHZZEventTree.fLep4Eta, weight);
	histAllLepsPhi_fullsim_m->Fill(fullsimHZZEventTree.fLep4Phi, weight);
	histAllLepsPt_fullsim_m->Fill(fullsimHZZEventTree.fLep4Pt, weight);
      }
    }
  }

  //*****************************************************************************************
  // Loop over fastsim events
  //*****************************************************************************************
  for(Int_t ientry=0; ientry < fastsimHZZEventTree.tree_->GetEntries(); ientry++) {       	
    fastsimHZZEventTree.tree_->GetEntry(ientry);
    if (ientry % 1000000 == 0) cout << "Event " << ientry << endl;

    //how to get lepton four momenta
    TLorentzVector lepton1FourVector;
    lepton1FourVector.SetPtEtaPhiM( fastsimHZZEventTree.fLep1Pt, fastsimHZZEventTree.fLep1Eta, fastsimHZZEventTree.fLep1Phi, ELECTRONMASS);

    //fill events passing selection

    // Setting up normalization options
    Float_t weight = fastsimHZZEventTree.fWeight; 

    // Filling Z masses
    if (fastsimHZZEventTree.fPassFullSelection) {
      if (abs(fastsimHZZEventTree.fLep1Type) == 11) {
	histZ1Mass_fastsim_ee->Fill(fastsimHZZEventTree.fZ1Mass, weight);
	histZ2Mass_fastsim_ee->Fill(fastsimHZZEventTree.fZ2Mass, weight);
      } else {
	histZ1Mass_fastsim_mm->Fill(fastsimHZZEventTree.fZ1Mass, weight);
	histZ2Mass_fastsim_mm->Fill(fastsimHZZEventTree.fZ2Mass, weight);
      }

      // Filling m4l masses (in 3 categories: 4e, 4m and 2e2m)
      if (abs(fastsimHZZEventTree.fLep1Type) == 11 && abs(fastsimHZZEventTree.fLep3Type) == 11) {
	histM4l_fastsim_4e->Fill(fastsimHZZEventTree.fFourLeptonMass, weight);
      } else if (abs(fastsimHZZEventTree.fLep1Type) == 13 && abs(fastsimHZZEventTree.fLep3Type) == 13) {
	histM4l_fastsim_4m->Fill(fastsimHZZEventTree.fFourLeptonMass, weight);
      } else {
	histM4l_fastsim_2e2m->Fill(fastsimHZZEventTree.fFourLeptonMass, weight);
      }

      // Filling all the lepton parameters 
      if (abs(fastsimHZZEventTree.fLep1Type) == 11) {
	histLep1Eta_fastsim_e->Fill(fastsimHZZEventTree.fLep1Eta, weight);
	histLep1Phi_fastsim_e->Fill(fastsimHZZEventTree.fLep1Phi, weight);
	histLep1Pt_fastsim_e->Fill(fastsimHZZEventTree.fLep1Pt, weight);
	histAllLepsEta_fastsim_e->Fill(fastsimHZZEventTree.fLep1Eta, weight);
	histAllLepsPhi_fastsim_e->Fill(fastsimHZZEventTree.fLep1Phi, weight);
	histAllLepsPt_fastsim_e->Fill(fastsimHZZEventTree.fLep1Pt, weight);
      } else {
	histLep1Eta_fastsim_m->Fill(fastsimHZZEventTree.fLep1Eta, weight);
	histLep1Phi_fastsim_m->Fill(fastsimHZZEventTree.fLep1Phi, weight);
	histLep1Pt_fastsim_m->Fill(fastsimHZZEventTree.fLep1Pt, weight);
	histAllLepsEta_fastsim_m->Fill(fastsimHZZEventTree.fLep1Eta, weight);
	histAllLepsPhi_fastsim_m->Fill(fastsimHZZEventTree.fLep1Phi, weight);
	histAllLepsPt_fastsim_m->Fill(fastsimHZZEventTree.fLep1Pt, weight);
      }

      if (abs(fastsimHZZEventTree.fLep2Type) == 11) {
	histLep2Eta_fastsim_e->Fill(fastsimHZZEventTree.fLep2Eta, weight);
	histLep2Phi_fastsim_e->Fill(fastsimHZZEventTree.fLep2Phi, weight);
	histLep2Pt_fastsim_e->Fill(fastsimHZZEventTree.fLep2Pt, weight);
	histAllLepsEta_fastsim_e->Fill(fastsimHZZEventTree.fLep2Eta, weight);
	histAllLepsPhi_fastsim_e->Fill(fastsimHZZEventTree.fLep2Phi, weight);
	histAllLepsPt_fastsim_e->Fill(fastsimHZZEventTree.fLep2Pt, weight);
      } else {
	histLep2Eta_fastsim_m->Fill(fastsimHZZEventTree.fLep2Eta, weight);
	histLep2Phi_fastsim_m->Fill(fastsimHZZEventTree.fLep2Phi, weight);
	histLep2Pt_fastsim_m->Fill(fastsimHZZEventTree.fLep2Pt, weight);
	histAllLepsEta_fastsim_m->Fill(fastsimHZZEventTree.fLep2Eta, weight);
	histAllLepsPhi_fastsim_m->Fill(fastsimHZZEventTree.fLep2Phi, weight);
	histAllLepsPt_fastsim_m->Fill(fastsimHZZEventTree.fLep2Pt, weight);
      }

      if (abs(fastsimHZZEventTree.fLep3Type) == 11) {
	histLep3Eta_fastsim_e->Fill(fastsimHZZEventTree.fLep3Eta, weight);
	histLep3Phi_fastsim_e->Fill(fastsimHZZEventTree.fLep3Phi, weight);
	histLep3Pt_fastsim_e->Fill(fastsimHZZEventTree.fLep3Pt, weight);
	histAllLepsEta_fastsim_e->Fill(fastsimHZZEventTree.fLep3Eta, weight);
	histAllLepsPhi_fastsim_e->Fill(fastsimHZZEventTree.fLep3Phi, weight);
	histAllLepsPt_fastsim_e->Fill(fastsimHZZEventTree.fLep3Pt, weight);
      } else {
	histLep3Eta_fastsim_m->Fill(fastsimHZZEventTree.fLep3Eta, weight);
	histLep3Phi_fastsim_m->Fill(fastsimHZZEventTree.fLep3Phi, weight);
	histLep3Pt_fastsim_m->Fill(fastsimHZZEventTree.fLep3Pt, weight);
	histAllLepsEta_fastsim_m->Fill(fastsimHZZEventTree.fLep3Eta, weight);
	histAllLepsPhi_fastsim_m->Fill(fastsimHZZEventTree.fLep3Phi, weight);
	histAllLepsPt_fastsim_m->Fill(fastsimHZZEventTree.fLep3Pt, weight);
      }

      if (abs(fastsimHZZEventTree.fLep4Type) == 11) {
	histLep4Eta_fastsim_e->Fill(fastsimHZZEventTree.fLep4Eta, weight);
	histLep4Phi_fastsim_e->Fill(fastsimHZZEventTree.fLep4Phi, weight);
	histLep4Pt_fastsim_e->Fill(fastsimHZZEventTree.fLep4Pt, weight);
	histAllLepsEta_fastsim_e->Fill(fastsimHZZEventTree.fLep4Eta, weight);
	histAllLepsPhi_fastsim_e->Fill(fastsimHZZEventTree.fLep4Phi, weight);
	histAllLepsPt_fastsim_e->Fill(fastsimHZZEventTree.fLep4Pt, weight);
      } else {
	histLep4Eta_fastsim_m->Fill(fastsimHZZEventTree.fLep4Eta, weight);
	histLep4Phi_fastsim_m->Fill(fastsimHZZEventTree.fLep4Phi, weight);
	histLep4Pt_fastsim_m->Fill(fastsimHZZEventTree.fLep4Pt, weight);
	histAllLepsEta_fastsim_m->Fill(fastsimHZZEventTree.fLep4Eta, weight);
	histAllLepsPhi_fastsim_m->Fill(fastsimHZZEventTree.fLep4Phi, weight);
	histAllLepsPt_fastsim_m->Fill(fastsimHZZEventTree.fLep4Pt, weight);
      }
    }
  }

  //********************************************************************
  // Now, normalizing if norm_option one was selected
  //********************************************************************
  if (norm_option) {
    NormalizeHist(histZ1Mass_fullsim_ee);
    NormalizeHist(histZ1Mass_fastsim_ee);
    NormalizeHist(histZ1Mass_fullsim_mm);
    NormalizeHist(histZ1Mass_fastsim_mm);
    NormalizeHist(histZ2Mass_fullsim_ee);
    NormalizeHist(histZ2Mass_fastsim_ee);
    NormalizeHist(histZ2Mass_fullsim_mm);
    NormalizeHist(histZ2Mass_fastsim_mm);

    NormalizeHist(histM4l_fastsim_4m);
    NormalizeHist(histM4l_fastsim_4e);
    NormalizeHist(histM4l_fastsim_2e2m);
    NormalizeHist(histM4l_fullsim_4m);
    NormalizeHist(histM4l_fullsim_4e);
    NormalizeHist(histM4l_fullsim_2e2m);

    NormalizeHist(histLep1Eta_fastsim_e);
    NormalizeHist(histLep1Phi_fastsim_e);
    NormalizeHist(histLep1Pt_fastsim_e);
    NormalizeHist(histLep1Eta_fullsim_e);
    NormalizeHist(histLep1Phi_fullsim_e);
    NormalizeHist(histLep1Pt_fullsim_e);

    NormalizeHist(histLep2Eta_fastsim_e);
    NormalizeHist(histLep2Phi_fastsim_e);
    NormalizeHist(histLep2Pt_fastsim_e);
    NormalizeHist(histLep2Eta_fullsim_e);
    NormalizeHist(histLep2Phi_fullsim_e);
    NormalizeHist(histLep2Pt_fullsim_e);

    NormalizeHist(histLep3Eta_fastsim_e);
    NormalizeHist(histLep3Phi_fastsim_e);
    NormalizeHist(histLep3Pt_fastsim_e);
    NormalizeHist(histLep3Eta_fullsim_e);
    NormalizeHist(histLep3Phi_fullsim_e);
    NormalizeHist(histLep3Pt_fullsim_e);

    NormalizeHist(histLep4Eta_fastsim_e);
    NormalizeHist(histLep4Phi_fastsim_e);
    NormalizeHist(histLep4Pt_fastsim_e);
    NormalizeHist(histLep4Eta_fullsim_e);
    NormalizeHist(histLep4Phi_fullsim_e);
    NormalizeHist(histLep4Pt_fullsim_e);

    NormalizeHist(histAllLepsEta_fastsim_e);
    NormalizeHist(histAllLepsPhi_fastsim_e);
    NormalizeHist(histAllLepsPt_fastsim_e);
    NormalizeHist(histAllLepsEta_fullsim_e);
    NormalizeHist(histAllLepsPhi_fullsim_e);
    NormalizeHist(histAllLepsPt_fullsim_e);

    NormalizeHist(histLep1Eta_fastsim_m);
    NormalizeHist(histLep1Phi_fastsim_m);
    NormalizeHist(histLep1Pt_fastsim_m);
    NormalizeHist(histLep1Eta_fullsim_m);
    NormalizeHist(histLep1Phi_fullsim_m);
    NormalizeHist(histLep1Pt_fullsim_m);

    NormalizeHist(histLep2Eta_fastsim_m);
    NormalizeHist(histLep2Phi_fastsim_m);
    NormalizeHist(histLep2Pt_fastsim_m);
    NormalizeHist(histLep2Eta_fullsim_m);
    NormalizeHist(histLep2Phi_fullsim_m);
    NormalizeHist(histLep2Pt_fullsim_m);

    NormalizeHist(histLep3Eta_fastsim_m);
    NormalizeHist(histLep3Phi_fastsim_m);
    NormalizeHist(histLep3Pt_fastsim_m);
    NormalizeHist(histLep3Eta_fullsim_m);
    NormalizeHist(histLep3Phi_fullsim_m);
    NormalizeHist(histLep3Pt_fullsim_m);

    NormalizeHist(histLep4Eta_fastsim_m);
    NormalizeHist(histLep4Phi_fastsim_m);
    NormalizeHist(histLep4Pt_fastsim_m);
    NormalizeHist(histLep4Eta_fullsim_m);
    NormalizeHist(histLep4Phi_fullsim_m);
    NormalizeHist(histLep4Pt_fullsim_m);

    NormalizeHist(histAllLepsEta_fastsim_m);
    NormalizeHist(histAllLepsPhi_fastsim_m);
    NormalizeHist(histAllLepsPt_fastsim_m);
    NormalizeHist(histAllLepsEta_fullsim_m);
    NormalizeHist(histAllLepsPhi_fullsim_m);
    NormalizeHist(histAllLepsPt_fullsim_m);
  }

  //*****************************************************************************************
  // Plot
  //*****************************************************************************************

  // Setting label offset
  histZ1Mass_fullsim_ee->GetYaxis()->SetTitleOffset(1.5);
  histZ1Mass_fastsim_ee->GetYaxis()->SetTitleOffset(1.5);
  histZ1Mass_fullsim_mm->GetYaxis()->SetTitleOffset(1.5);
  histZ1Mass_fastsim_mm->GetYaxis()->SetTitleOffset(1.5);
  histZ2Mass_fullsim_ee->GetYaxis()->SetTitleOffset(1.5);
  histZ2Mass_fastsim_ee->GetYaxis()->SetTitleOffset(1.5);
  histZ2Mass_fullsim_mm->GetYaxis()->SetTitleOffset(1.5);
  histZ2Mass_fastsim_mm->GetYaxis()->SetTitleOffset(1.5);

  histM4l_fastsim_4m->GetYaxis()->SetTitleOffset(1.5);
  histM4l_fastsim_4e->GetYaxis()->SetTitleOffset(1.5);
  histM4l_fastsim_2e2m->GetYaxis()->SetTitleOffset(1.5);
  histM4l_fullsim_4m->GetYaxis()->SetTitleOffset(1.5);
  histM4l_fullsim_4e->GetYaxis()->SetTitleOffset(1.5);
  histM4l_fullsim_2e2m->GetYaxis()->SetTitleOffset(1.5);

  histLep1Eta_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep1Phi_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep1Pt_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep1Eta_fullsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep1Phi_fullsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep1Pt_fullsim_e->GetYaxis()->SetTitleOffset(1.5);

  histLep2Eta_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep2Phi_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep2Pt_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep2Eta_fullsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep2Phi_fullsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep2Pt_fullsim_e->GetYaxis()->SetTitleOffset(1.5);

  histLep3Eta_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep3Phi_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep3Pt_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep3Eta_fullsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep3Phi_fullsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep3Pt_fullsim_e->GetYaxis()->SetTitleOffset(1.5);

  histLep4Eta_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep4Phi_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep4Pt_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep4Eta_fullsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep4Phi_fullsim_e->GetYaxis()->SetTitleOffset(1.5);
  histLep4Pt_fullsim_e->GetYaxis()->SetTitleOffset(1.5);

  histAllLepsEta_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histAllLepsPhi_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histAllLepsPt_fastsim_e->GetYaxis()->SetTitleOffset(1.5);
  histAllLepsEta_fullsim_e->GetYaxis()->SetTitleOffset(1.5);
  histAllLepsPhi_fullsim_e->GetYaxis()->SetTitleOffset(1.5);
  histAllLepsPt_fullsim_e->GetYaxis()->SetTitleOffset(1.5);

  histLep1Eta_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep1Phi_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep1Pt_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep1Eta_fullsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep1Phi_fullsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep1Pt_fullsim_m->GetYaxis()->SetTitleOffset(1.5);

  histLep2Eta_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep2Phi_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep2Pt_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep2Eta_fullsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep2Phi_fullsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep2Pt_fullsim_m->GetYaxis()->SetTitleOffset(1.5);

  histLep3Eta_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep3Phi_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep3Pt_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep3Eta_fullsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep3Phi_fullsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep3Pt_fullsim_m->GetYaxis()->SetTitleOffset(1.5);

  histLep4Eta_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep4Phi_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep4Pt_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep4Eta_fullsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep4Phi_fullsim_m->GetYaxis()->SetTitleOffset(1.5);
  histLep4Pt_fullsim_m->GetYaxis()->SetTitleOffset(1.5);

  histAllLepsEta_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histAllLepsPhi_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histAllLepsPt_fastsim_m->GetYaxis()->SetTitleOffset(1.5);
  histAllLepsEta_fullsim_m->GetYaxis()->SetTitleOffset(1.5);
  histAllLepsPhi_fullsim_m->GetYaxis()->SetTitleOffset(1.5);
  histAllLepsPt_fullsim_m->GetYaxis()->SetTitleOffset(1.5);

  // Saving all histograms to file
  outRootFile.Write();

  // Plotting all histograms
  PlotTwo(histZ1Mass_fullsim_ee, histZ1Mass_fastsim_ee, (outputDirectory + "ValidateParams_Z1Mass_ee" + label + ".gif").c_str(), error_option);
  PlotTwo(histZ1Mass_fullsim_mm, histZ1Mass_fastsim_mm, (outputDirectory + "ValidateParams_Z1Mass_mm" + label + ".gif").c_str(), error_option);

  PlotTwo(histZ2Mass_fullsim_ee, histZ2Mass_fastsim_ee, (outputDirectory + "ValidateParams_Z2Mass_ee" + label + ".gif").c_str(), error_option);
  PlotTwo(histZ2Mass_fullsim_ee, histZ2Mass_fastsim_mm, (outputDirectory + "ValidateParams_Z2Mass_mm" + label + ".gif").c_str(), error_option);

  PlotTwo(histM4l_fullsim_4e, histM4l_fastsim_4e, (outputDirectory + "ValidateParams_M4l_4e" + label + ".gif").c_str(), error_option);
  PlotTwo(histM4l_fullsim_4m, histM4l_fastsim_4m, (outputDirectory + "ValidateParams_M4l_4m" + label + ".gif").c_str(), error_option);
  PlotTwo(histM4l_fullsim_2e2m, histM4l_fastsim_2e2m, (outputDirectory + "ValidateParams_M4l_2e2m" + label + ".gif").c_str(), error_option);

  PlotTwo(histLep1Eta_fullsim_e, histLep1Eta_fastsim_e, (outputDirectory + "ValidateParams_Lep1Eta_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep1Eta_fullsim_m, histLep1Eta_fastsim_m, (outputDirectory + "ValidateParams_Lep1Eta_m" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep1Phi_fullsim_e, histLep1Phi_fastsim_e, (outputDirectory + "ValidateParams_Lep1Phi_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep1Phi_fullsim_m, histLep1Phi_fastsim_m, (outputDirectory + "ValidateParams_Lep1Phi_m" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep1Pt_fullsim_e, histLep1Pt_fastsim_e, (outputDirectory + "ValidateParams_Lep1Pt_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep1Pt_fullsim_m, histLep1Pt_fastsim_m, (outputDirectory + "ValidateParams_Lep1Pt_m" + label + ".gif").c_str(), error_option);

  PlotTwo(histLep2Eta_fullsim_e, histLep2Eta_fastsim_e, (outputDirectory + "ValidateParams_Lep2Eta_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep2Eta_fullsim_m, histLep2Eta_fastsim_m, (outputDirectory + "ValidateParams_Lep2Eta_m" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep2Phi_fullsim_e, histLep2Phi_fastsim_e, (outputDirectory + "ValidateParams_Lep2Phi_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep2Phi_fullsim_m, histLep2Phi_fastsim_m, (outputDirectory + "ValidateParams_Lep2Phi_m" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep2Pt_fullsim_e, histLep2Pt_fastsim_e, (outputDirectory + "ValidateParams_Lep2Pt_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep2Pt_fullsim_m, histLep2Pt_fastsim_m, (outputDirectory + "ValidateParams_Lep2Pt_m" + label + ".gif").c_str(), error_option);

  PlotTwo(histLep3Eta_fullsim_e, histLep3Eta_fastsim_e, (outputDirectory + "ValidateParams_Lep3Eta_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep3Eta_fullsim_m, histLep3Eta_fastsim_m, (outputDirectory + "ValidateParams_Lep3Eta_m" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep3Phi_fullsim_e, histLep3Phi_fastsim_e, (outputDirectory + "ValidateParams_Lep3Phi_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep3Phi_fullsim_m, histLep3Phi_fastsim_m, (outputDirectory + "ValidateParams_Lep3Phi_m" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep3Pt_fullsim_e, histLep3Pt_fastsim_e, (outputDirectory + "ValidateParams_Lep3Pt_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep3Pt_fullsim_m, histLep3Pt_fastsim_m, (outputDirectory + "ValidateParams_Lep3Pt_m" + label + ".gif").c_str(), error_option);

  PlotTwo(histLep4Eta_fullsim_e, histLep4Eta_fastsim_e, (outputDirectory + "ValidateParams_Lep4Eta_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep4Eta_fullsim_m, histLep4Eta_fastsim_m, (outputDirectory + "ValidateParams_Lep4Eta_m" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep4Phi_fullsim_e, histLep4Phi_fastsim_e, (outputDirectory + "ValidateParams_Lep4Phi_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep4Phi_fullsim_m, histLep4Phi_fastsim_m, (outputDirectory + "ValidateParams_Lep4Phi_m" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep4Pt_fullsim_e, histLep4Pt_fastsim_e, (outputDirectory + "ValidateParams_Lep4Pt_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLep4Pt_fullsim_m, histLep4Pt_fastsim_m, (outputDirectory + "ValidateParams_Lep4Pt_m" + label + ".gif").c_str(), error_option);

  PlotTwo(histAllLepsEta_fullsim_e, histAllLepsEta_fastsim_e, (outputDirectory + "ValidateParams_AllLepsEta_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histAllLepsEta_fullsim_m, histAllLepsEta_fastsim_m, (outputDirectory + "ValidateParams_AllLepsEta_m" + label + ".gif").c_str(), error_option);
  PlotTwo(histAllLepsPhi_fullsim_e, histAllLepsPhi_fastsim_e, (outputDirectory + "ValidateParams_AllLepsPhi_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histAllLepsPhi_fullsim_m, histAllLepsPhi_fastsim_m, (outputDirectory + "ValidateParams_AllLepsPhi_m" + label + ".gif").c_str(), error_option);
  PlotTwo(histAllLepsPt_fullsim_e, histAllLepsPt_fastsim_e, (outputDirectory + "ValidateParams_AllLepsPt_e" + label + ".gif").c_str(), error_option);
  PlotTwo(histAllLepsPt_fullsim_m, histAllLepsPt_fastsim_m, (outputDirectory + "ValidateParams_AllLepsPt_m" + label + ".gif").c_str(), error_option);
}
