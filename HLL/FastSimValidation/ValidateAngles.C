//================================================================================================
//					ValidateAngles
// 
// Reads kinematic parameters from fullsim and fastsim ntuples, calculates kinematic angles and
// plots fullsim and fastsim for comparison
//
// USAGE
//
// ValidateAngles(string fullsimFilename, string fastsimFilename, string outputDirectory, Int_t norm_option = 1,  Int_t error_option = 2, const string Label = "ZZ")
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
void ValidateAngles(string fullsimFilename, string fastsimFilename, string outputDirectory, Int_t norm_option = 1,  Int_t error_option = 2, const string Label = "ZZ") {

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
  TFile outRootFile((outputDirectory + "ValidateAngles_AllHistograms" + label + ".root").c_str(), "recreate");

  TH1F* histPhi0_4m_fullsim = new TH1F("histPhi0_4m_fullsim", "Phi0, 4m;Phi0;Number of Events", 50, 0, 6.3);
  TH1F* histPhi0_4m_fastsim = new TH1F("histPhi0_4m_fastsim", "Phi0, 4m;Phi0;Number of Events", 50, 0, 6.3);
  TH1F* histPhi0_4e_fullsim = new TH1F("histPhi0_4e_fullsim", "Phi0, 4e;Phi0;Number of Events", 50, 0, 6.3);
  TH1F* histPhi0_4e_fastsim = new TH1F("histPhi0_4e_fastsim", "Phi0, 4e;Phi0;Number of Events", 50, 0, 6.3);
  TH1F* histPhi0_2e2m_fullsim = new TH1F("histPhi0_2e2m_fullsim", "Phi0, 2e2m;Phi0;Number of Events", 50, 0, 6.3);
  TH1F* histPhi0_2e2m_fastsim = new TH1F("histPhi0_2e2m_fastsim", "Phi0, 2e2m;Phi0;Number of Events", 50, 0, 6.3);

  TH1F* histTheta0_4m_fullsim = new TH1F("histTheta0_4m_fullsim", "Theta0, 4m;Theta0;Number of Events", 50, 0, 3.2);
  TH1F* histTheta0_4m_fastsim = new TH1F("histTheta0_4m_fastsim", "Theta0, 4m;Theta0;Number of Events", 50, 0, 3.2);
  TH1F* histTheta0_4e_fullsim = new TH1F("histTheta0_4e_fullsim", "Theta0, 4e;Theta0;Number of Events", 50, 0, 3.2);
  TH1F* histTheta0_4e_fastsim = new TH1F("histTheta0_4e_fastsim", "Theta0, 4e;Theta0;Number of Events", 50, 0, 3.2);
  TH1F* histTheta0_2e2m_fullsim = new TH1F("histTheta0_2e2m_fullsim", "Theta0, 2e2m;Theta0;Number of Events", 50, 0, 3.2);
  TH1F* histTheta0_2e2m_fastsim = new TH1F("histTheta0_2e2m_fastsim", "Theta0, 2e2m;Theta0;Number of Events", 50, 0, 3.2);

  TH1F* histPhi_4m_fullsim = new TH1F("histPhi_4m_fullsim", "Phi, 4m;Phi;Number of Events", 50, 0, 6.3);
  TH1F* histPhi_4m_fastsim = new TH1F("histPhi_4m_fastsim", "Phi, 4m;Phi;Number of Events", 50, 0, 6.3);
  TH1F* histPhi_4e_fullsim = new TH1F("histPhi_4e_fullsim", "Phi, 4e;Phi;Number of Events", 50, 0, 6.3);
  TH1F* histPhi_4e_fastsim = new TH1F("histPhi_4e_fastsim", "Phi, 4e;Phi;Number of Events", 50, 0, 6.3);
  TH1F* histPhi_2e2m_fullsim = new TH1F("histPhi_2e2m_fullsim", "Phi, 2e2m;Phi;Number of Events", 50, 0, 6.3);
  TH1F* histPhi_2e2m_fastsim = new TH1F("histPhi_2e2m_fastsim", "Phi, 2e2m;Phi;Number of Events", 50, 0, 6.3);

  TH1F* histTheta1_4m_fullsim = new TH1F("histTheta1_4m_fullsim", "Theta1, 4m;Theta1;Number of Events", 50, 0, 3.2);
  TH1F* histTheta1_4m_fastsim = new TH1F("histTheta1_4m_fastsim", "Theta1, 4m;Theta1;Number of Events", 50, 0, 3.2);
  TH1F* histTheta1_4e_fullsim = new TH1F("histTheta1_4e_fullsim", "Theta1, 4e;Theta1;Number of Events", 50, 0, 3.2);
  TH1F* histTheta1_4e_fastsim = new TH1F("histTheta1_4e_fastsim", "Theta1, 4e;Theta1;Number of Events", 50, 0, 3.2);
  TH1F* histTheta1_2e2m_fullsim = new TH1F("histTheta1_2e2m_fullsim", "Theta1, 2e2m;Theta1;Number of Events", 50, 0, 3.2);
  TH1F* histTheta1_2e2m_fastsim = new TH1F("histTheta1_2e2m_fastsim", "Theta1, 2e2m;Theta1;Number of Events", 50, 0, 3.2);

  TH1F* histTheta2_4m_fullsim = new TH1F("histTheta2_4m_fullsim", "Theta2, 4m;Theta2;Number of Events", 50, 0, 3.2);
  TH1F* histTheta2_4m_fastsim = new TH1F("histTheta2_4m_fastsim", "Theta2, 4m;Theta2;Number of Events", 50, 0, 3.2);
  TH1F* histTheta2_4e_fullsim = new TH1F("histTheta2_4e_fullsim", "Theta2, 4e;Theta2;Number of Events", 50, 0, 3.2);
  TH1F* histTheta2_4e_fastsim = new TH1F("histTheta2_4e_fastsim", "Theta2, 4e;Theta2;Number of Events", 50, 0, 3.2);
  TH1F* histTheta2_2e2m_fullsim = new TH1F("histTheta2_2e2m_fullsim", "Theta2, 2e2m;Theta2;Number of Events", 50, 0, 3.2);
  TH1F* histTheta2_2e2m_fastsim = new TH1F("histTheta2_2e2m_fastsim", "Theta2, 2e2m;Theta2;Number of Events", 50, 0, 3.2);

  TH1F* histHMass_4m_fullsim = new TH1F("histHMass_4m_fullsim", "HMass, 4m;HMass;Number of Events", 100, 0, 300);
  TH1F* histHMass_4m_fastsim = new TH1F("histHMass_4m_fastsim", "HMass, 4m;HMass;Number of Events", 100, 0, 300);
  TH1F* histHMass_4e_fullsim = new TH1F("histHMass_4e_fullsim", "HMass, 4e;HMass;Number of Events", 100, 0, 300);
  TH1F* histHMass_4e_fastsim = new TH1F("histHMass_4e_fastsim", "HMass, 4e;HMass;Number of Events", 100, 0, 300);
  TH1F* histHMass_2e2m_fullsim = new TH1F("histHMass_2e2m_fullsim", "HMass, 2e2m;HMass;Number of Events", 100, 0, 300);
  TH1F* histHMass_2e2m_fastsim = new TH1F("histHMass_2e2m_fastsim", "HMass, 2e2m;HMass;Number of Events", 100, 0, 300);

  TH1F* histZMass_4m_fullsim = new TH1F("histZMass_4m_fullsim", "ZMass, 4m;ZMass;Number of Events", 100, 40, 140);
  TH1F* histZMass_4m_fastsim = new TH1F("histZMass_4m_fastsim", "ZMass, 4m;ZMass;Number of Events", 100, 40, 140);
  TH1F* histZMass_4e_fullsim = new TH1F("histZMass_4e_fullsim", "ZMass, 4e;ZMass;Number of Events", 100, 40, 140);
  TH1F* histZMass_4e_fastsim = new TH1F("histZMass_4e_fastsim", "ZMass, 4e;ZMass;Number of Events", 100, 40, 140);
  TH1F* histZMass_2e2m_fullsim = new TH1F("histZMass_2e2m_fullsim", "ZMass, 2e2m;ZMass;Number of Events", 100, 40, 140);
  TH1F* histZMass_2e2m_fastsim = new TH1F("histZMass_2e2m_fastsim", "ZMass, 2e2m;ZMass;Number of Events", 100, 40, 140);

  TH1F* histZ2Mass_4m_fullsim = new TH1F("histZ2Mass_4m_fullsim", "Z2Mass, 4m;Z2Mass;Number of Events", 150, 0, 150);
  TH1F* histZ2Mass_4m_fastsim = new TH1F("histZ2Mass_4m_fastsim", "Z2Mass, 4m;Z2Mass;Number of Events", 150, 0, 150);
  TH1F* histZ2Mass_4e_fullsim = new TH1F("histZ2Mass_4e_fullsim", "Z2Mass, 4e;Z2Mass;Number of Events", 150, 0, 150);
  TH1F* histZ2Mass_4e_fastsim = new TH1F("histZ2Mass_4e_fastsim", "Z2Mass, 4e;Z2Mass;Number of Events", 150, 0, 150);
  TH1F* histZ2Mass_2e2m_fullsim = new TH1F("histZ2Mass_2e2m_fullsim", "Z2Mass, 2e2m;Z2Mass;Number of Events", 150, 0, 150);
  TH1F* histZ2Mass_2e2m_fastsim = new TH1F("histZ2Mass_2e2m_fastsim", "Z2Mass, 2e2m;Z2Mass;Number of Events", 150, 0, 150);

  TH1F* histPhiH_4m_fullsim = new TH1F("histPhiH_4m_fullsim", "PhiH, 4m;PhiH;Number of Events", 100, 0, 0.5);
  TH1F* histPhiH_4m_fastsim = new TH1F("histPhiH_4m_fastsim", "PhiH, 4m;PhiH;Number of Events", 100, 0, 0.5);
  TH1F* histPhiH_4e_fullsim = new TH1F("histPhiH_4e_fullsim", "PhiH, 4e;PhiH;Number of Events", 100, 0, 0.5);
  TH1F* histPhiH_4e_fastsim = new TH1F("histPhiH_4e_fastsim", "PhiH, 4e;PhiH;Number of Events", 100, 0, 0.5);
  TH1F* histPhiH_2e2m_fullsim = new TH1F("histPhiH_2e2m_fullsim", "PhiH, 2e2m;PhiH;Number of Events", 100, 0, 0.5);
  TH1F* histPhiH_2e2m_fastsim = new TH1F("histPhiH_2e2m_fastsim", "PhiH, 2e2m;PhiH;Number of Events", 100, 0, 0.5);

  //*****************************************************************************************
  // Loop over fullsim events
  //*****************************************************************************************
  for(UInt_t ientry=0; ientry < fullsimHZZEventTree.tree_->GetEntries(); ientry++) {       	
    fullsimHZZEventTree.tree_->GetEntry(ientry);
    if (ientry % 1000000 == 0) cout << "Event " << ientry << endl;

    //how to get lepton four momenta
    FourVector lepton1FourVector, lepton2FourVector, lepton3FourVector, lepton4FourVector;
    if (abs(fullsimHZZEventTree.fLep1Type) == 11) 
      lepton1FourVector.SetPtEtaPhiMass( fullsimHZZEventTree.fLep1Pt, fullsimHZZEventTree.fLep1Eta, fullsimHZZEventTree.fLep1Phi, ELECTRONMASS);
    else if (abs(fullsimHZZEventTree.fLep1Type) == 13)
      lepton1FourVector.SetPtEtaPhiMass( fullsimHZZEventTree.fLep1Pt, fullsimHZZEventTree.fLep1Eta, fullsimHZZEventTree.fLep1Phi, MUONMASS);

    if (abs(fullsimHZZEventTree.fLep2Type) == 11) 
      lepton2FourVector.SetPtEtaPhiMass( fullsimHZZEventTree.fLep2Pt, fullsimHZZEventTree.fLep2Eta, fullsimHZZEventTree.fLep2Phi, ELECTRONMASS);
    else if (abs(fullsimHZZEventTree.fLep2Type) == 13)
      lepton2FourVector.SetPtEtaPhiMass( fullsimHZZEventTree.fLep2Pt, fullsimHZZEventTree.fLep2Eta, fullsimHZZEventTree.fLep2Phi, MUONMASS);

    if (abs(fullsimHZZEventTree.fLep3Type) == 11) 
      lepton3FourVector.SetPtEtaPhiMass( fullsimHZZEventTree.fLep3Pt, fullsimHZZEventTree.fLep3Eta, fullsimHZZEventTree.fLep3Phi, ELECTRONMASS);
    else if (abs(fullsimHZZEventTree.fLep3Type) == 13)
      lepton3FourVector.SetPtEtaPhiMass( fullsimHZZEventTree.fLep3Pt, fullsimHZZEventTree.fLep3Eta, fullsimHZZEventTree.fLep3Phi, MUONMASS);

    if (abs(fullsimHZZEventTree.fLep4Type) == 11) 
      lepton4FourVector.SetPtEtaPhiMass( fullsimHZZEventTree.fLep4Pt, fullsimHZZEventTree.fLep4Eta, fullsimHZZEventTree.fLep4Phi, ELECTRONMASS);
    else if (abs(fullsimHZZEventTree.fLep4Type) == 13)
      lepton4FourVector.SetPtEtaPhiMass( fullsimHZZEventTree.fLep4Pt, fullsimHZZEventTree.fLep4Eta, fullsimHZZEventTree.fLep4Phi, MUONMASS);


    // Getting correct weight for each event
    Float_t weight = fullsimHZZEventTree.fWeight; 

    // Using vector-angle converter to find the angles
    LeptonVectors leptons;
    leptons.Lepton11 = lepton1FourVector;
    leptons.Lepton12 = lepton2FourVector;
    leptons.Lepton21 = lepton3FourVector;
    leptons.Lepton22 = lepton4FourVector;
    EventParameters thisEvent = ConvertVectorsToAngles(leptons);

    // Filling event parameters
    if (fullsimHZZEventTree.fPassFullSelection) {

      if (abs(fullsimHZZEventTree.fLep1Type) == 11 && abs(fullsimHZZEventTree.fLep3Type) == 11) {
	histPhi0_4e_fullsim->Fill(thisEvent.Phi0, weight);
	histTheta0_4e_fullsim->Fill(thisEvent.Theta0, weight);
	histPhi_4e_fullsim->Fill(thisEvent.Phi, weight);
	histTheta1_4e_fullsim->Fill(thisEvent.Theta1, weight);
	histTheta2_4e_fullsim->Fill(thisEvent.Theta2, weight);
	histHMass_4e_fullsim->Fill(thisEvent.HMass, weight);
	histZMass_4e_fullsim->Fill(thisEvent.ZMass, weight);
	histZ2Mass_4e_fullsim->Fill(thisEvent.Z2Mass, weight);
	histPhiH_4e_fullsim->Fill(thisEvent.PhiH, weight);

      } else if (abs(fullsimHZZEventTree.fLep1Type) == 13 && abs(fullsimHZZEventTree.fLep3Type) == 13) {
	histPhi0_4m_fullsim->Fill(thisEvent.Phi0, weight);
	histTheta0_4m_fullsim->Fill(thisEvent.Theta0, weight);
	histPhi_4m_fullsim->Fill(thisEvent.Phi, weight);
	histTheta1_4m_fullsim->Fill(thisEvent.Theta1, weight);
	histTheta2_4m_fullsim->Fill(thisEvent.Theta2, weight);
	histHMass_4m_fullsim->Fill(thisEvent.HMass, weight);
	histZMass_4m_fullsim->Fill(thisEvent.ZMass, weight);
	histZ2Mass_4m_fullsim->Fill(thisEvent.Z2Mass, weight);
	histPhiH_4m_fullsim->Fill(thisEvent.PhiH, weight);

      } else {
	histPhi0_2e2m_fullsim->Fill(thisEvent.Phi0, weight);
	histTheta0_2e2m_fullsim->Fill(thisEvent.Theta0, weight);
	histPhi_2e2m_fullsim->Fill(thisEvent.Phi, weight);
	histTheta1_2e2m_fullsim->Fill(thisEvent.Theta1, weight);
	histTheta2_2e2m_fullsim->Fill(thisEvent.Theta2, weight);
	histHMass_2e2m_fullsim->Fill(thisEvent.HMass, weight);
	histZMass_2e2m_fullsim->Fill(thisEvent.ZMass, weight);
	histZ2Mass_2e2m_fullsim->Fill(thisEvent.Z2Mass, weight);
	histPhiH_2e2m_fullsim->Fill(thisEvent.PhiH, weight);
      }
    }
  }

  //*****************************************************************************************
  // Loop over fastsim events
  //*****************************************************************************************
  for(UInt_t ientry=0; ientry < fastsimHZZEventTree.tree_->GetEntries(); ientry++) {       	
    fastsimHZZEventTree.tree_->GetEntry(ientry);
    if (ientry % 1000000 == 0) cout << "Event " << ientry << endl;

    //how to get lepton four momenta
    FourVector lepton1FourVector, lepton2FourVector, lepton3FourVector, lepton4FourVector;
    if (abs(fastsimHZZEventTree.fLep1Type) == 11) 
      lepton1FourVector.SetPtEtaPhiMass( fastsimHZZEventTree.fLep1Pt, fastsimHZZEventTree.fLep1Eta, fastsimHZZEventTree.fLep1Phi, ELECTRONMASS);
    else if (abs(fastsimHZZEventTree.fLep1Type) == 13)
      lepton1FourVector.SetPtEtaPhiMass( fastsimHZZEventTree.fLep1Pt, fastsimHZZEventTree.fLep1Eta, fastsimHZZEventTree.fLep1Phi, MUONMASS);

    if (abs(fastsimHZZEventTree.fLep2Type) == 11) 
      lepton2FourVector.SetPtEtaPhiMass( fastsimHZZEventTree.fLep2Pt, fastsimHZZEventTree.fLep2Eta, fastsimHZZEventTree.fLep2Phi, ELECTRONMASS);
    else if (abs(fastsimHZZEventTree.fLep2Type) == 13)
      lepton2FourVector.SetPtEtaPhiMass( fastsimHZZEventTree.fLep2Pt, fastsimHZZEventTree.fLep2Eta, fastsimHZZEventTree.fLep2Phi, MUONMASS);

    if (abs(fastsimHZZEventTree.fLep3Type) == 11) 
      lepton3FourVector.SetPtEtaPhiMass( fastsimHZZEventTree.fLep3Pt, fastsimHZZEventTree.fLep3Eta, fastsimHZZEventTree.fLep3Phi, ELECTRONMASS);
    else if (abs(fastsimHZZEventTree.fLep3Type) == 13)
      lepton3FourVector.SetPtEtaPhiMass( fastsimHZZEventTree.fLep3Pt, fastsimHZZEventTree.fLep3Eta, fastsimHZZEventTree.fLep3Phi, MUONMASS);

    if (abs(fastsimHZZEventTree.fLep4Type) == 11) 
      lepton4FourVector.SetPtEtaPhiMass( fastsimHZZEventTree.fLep4Pt, fastsimHZZEventTree.fLep4Eta, fastsimHZZEventTree.fLep4Phi, ELECTRONMASS);
    else if (abs(fastsimHZZEventTree.fLep4Type) == 13)
      lepton4FourVector.SetPtEtaPhiMass( fastsimHZZEventTree.fLep4Pt, fastsimHZZEventTree.fLep4Eta, fastsimHZZEventTree.fLep4Phi, MUONMASS);


    // Weight for each event
    Float_t weight = fastsimHZZEventTree.fWeight; 

    // Using vector-angle converter to find the angles
    LeptonVectors leptons;
    leptons.Lepton11 = lepton1FourVector;
    leptons.Lepton12 = lepton2FourVector;
    leptons.Lepton21 = lepton3FourVector;
    leptons.Lepton22 = lepton4FourVector;
    EventParameters thisEvent = ConvertVectorsToAngles(leptons);

    // Filling event parameters
    if (fastsimHZZEventTree.fPassFullSelection) {

      if (abs(fastsimHZZEventTree.fLep1Type) == 11 && abs(fastsimHZZEventTree.fLep3Type) == 11) {
	histPhi0_4e_fastsim->Fill(thisEvent.Phi0, weight);
	histTheta0_4e_fastsim->Fill(thisEvent.Theta0, weight);
	histPhi_4e_fastsim->Fill(thisEvent.Phi, weight);
	histTheta1_4e_fastsim->Fill(thisEvent.Theta1, weight);
	histTheta2_4e_fastsim->Fill(thisEvent.Theta2, weight);
	histHMass_4e_fastsim->Fill(thisEvent.HMass, weight);
	histZMass_4e_fastsim->Fill(thisEvent.ZMass, weight);
	histZ2Mass_4e_fastsim->Fill(thisEvent.Z2Mass, weight);
	histPhiH_4e_fastsim->Fill(thisEvent.PhiH, weight);

      } else if (abs(fastsimHZZEventTree.fLep1Type) == 13 && abs(fastsimHZZEventTree.fLep3Type) == 13) {
	histPhi0_4m_fastsim->Fill(thisEvent.Phi0, weight);
	histTheta0_4m_fastsim->Fill(thisEvent.Theta0, weight);
	histPhi_4m_fastsim->Fill(thisEvent.Phi, weight);
	histTheta1_4m_fastsim->Fill(thisEvent.Theta1, weight);
	histTheta2_4m_fastsim->Fill(thisEvent.Theta2, weight);
	histHMass_4m_fastsim->Fill(thisEvent.HMass, weight);
	histZMass_4m_fastsim->Fill(thisEvent.ZMass, weight);
	histZ2Mass_4m_fastsim->Fill(thisEvent.Z2Mass, weight);
	histPhiH_4m_fastsim->Fill(thisEvent.PhiH, weight);

      } else {
	histPhi0_2e2m_fastsim->Fill(thisEvent.Phi0, weight);
	histTheta0_2e2m_fastsim->Fill(thisEvent.Theta0, weight);
	histPhi_2e2m_fastsim->Fill(thisEvent.Phi, weight);
	histTheta1_2e2m_fastsim->Fill(thisEvent.Theta1, weight);
	histTheta2_2e2m_fastsim->Fill(thisEvent.Theta2, weight);
	histHMass_2e2m_fastsim->Fill(thisEvent.HMass, weight);
	histZMass_2e2m_fastsim->Fill(thisEvent.ZMass, weight);
	histZ2Mass_2e2m_fastsim->Fill(thisEvent.Z2Mass, weight);
	histPhiH_2e2m_fastsim->Fill(thisEvent.PhiH, weight);
      }
    }
  }

  //********************************************************************
  // Now, normalizing if norm_option one was selected
  //********************************************************************
  if (norm_option) {
    NormalizeHist(histPhi0_4e_fullsim);
    NormalizeHist(histPhi0_4m_fullsim);
    NormalizeHist(histPhi0_2e2m_fullsim);
    NormalizeHist(histTheta0_4e_fullsim);
    NormalizeHist(histTheta0_4m_fullsim);
    NormalizeHist(histTheta0_2e2m_fullsim);
    NormalizeHist(histPhi_4e_fullsim);
    NormalizeHist(histPhi_4m_fullsim);
    NormalizeHist(histPhi_2e2m_fullsim);
    NormalizeHist(histTheta1_4e_fullsim);
    NormalizeHist(histTheta1_4m_fullsim);
    NormalizeHist(histTheta1_2e2m_fullsim);
    NormalizeHist(histTheta2_4e_fullsim);
    NormalizeHist(histTheta2_4m_fullsim);
    NormalizeHist(histTheta2_2e2m_fullsim);
    NormalizeHist(histHMass_4e_fullsim);
    NormalizeHist(histHMass_4m_fullsim);
    NormalizeHist(histHMass_2e2m_fullsim);
    NormalizeHist(histZMass_4e_fullsim);
    NormalizeHist(histZMass_4m_fullsim);
    NormalizeHist(histZMass_2e2m_fullsim);
    NormalizeHist(histZ2Mass_4e_fullsim);
    NormalizeHist(histZ2Mass_4m_fullsim);
    NormalizeHist(histZ2Mass_2e2m_fullsim);
    NormalizeHist(histPhiH_4e_fullsim);
    NormalizeHist(histPhiH_4m_fullsim);
    NormalizeHist(histPhiH_2e2m_fullsim);

    NormalizeHist(histPhi0_4e_fastsim);
    NormalizeHist(histPhi0_4m_fastsim);
    NormalizeHist(histPhi0_2e2m_fastsim);
    NormalizeHist(histTheta0_4e_fastsim);
    NormalizeHist(histTheta0_4m_fastsim);
    NormalizeHist(histTheta0_2e2m_fastsim);
    NormalizeHist(histPhi_4e_fastsim);
    NormalizeHist(histPhi_4m_fastsim);
    NormalizeHist(histPhi_2e2m_fastsim);
    NormalizeHist(histTheta1_4e_fastsim);
    NormalizeHist(histTheta1_4m_fastsim);
    NormalizeHist(histTheta1_2e2m_fastsim);
    NormalizeHist(histTheta2_4e_fastsim);
    NormalizeHist(histTheta2_4m_fastsim);
    NormalizeHist(histTheta2_2e2m_fastsim);
    NormalizeHist(histHMass_4e_fastsim);
    NormalizeHist(histHMass_4m_fastsim);
    NormalizeHist(histHMass_2e2m_fastsim);
    NormalizeHist(histZMass_4e_fastsim);
    NormalizeHist(histZMass_4m_fastsim);
    NormalizeHist(histZMass_2e2m_fastsim);
    NormalizeHist(histZ2Mass_4e_fastsim);
    NormalizeHist(histZ2Mass_4m_fastsim);
    NormalizeHist(histZ2Mass_2e2m_fastsim);
    NormalizeHist(histPhiH_4e_fastsim);
    NormalizeHist(histPhiH_4m_fastsim);
    NormalizeHist(histPhiH_2e2m_fastsim);
  }


  //*****************************************************************************************
  // Plot
  //*****************************************************************************************

  // Setting label offset
  histPhi0_4e_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histPhi0_4m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histPhi0_2e2m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta0_4e_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta0_4m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta0_2e2m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histPhi_4e_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histPhi_4m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histPhi_2e2m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta1_4e_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta1_4m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta1_2e2m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta2_4e_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta2_4m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta2_2e2m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histHMass_4e_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histHMass_4m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histHMass_2e2m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histZMass_4e_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histZMass_4m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histZMass_2e2m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histZ2Mass_4e_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histZ2Mass_4m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histZ2Mass_2e2m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histPhiH_4e_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histPhiH_4m_fullsim->GetYaxis()->SetTitleOffset(1.4);
  histPhiH_2e2m_fullsim->GetYaxis()->SetTitleOffset(1.4);

  histPhi0_4e_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histPhi0_4m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histPhi0_2e2m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta0_4e_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta0_4m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta0_2e2m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histPhi_4e_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histPhi_4m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histPhi_2e2m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta1_4e_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta1_4m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta1_2e2m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta2_4e_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta2_4m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histTheta2_2e2m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histHMass_4e_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histHMass_4m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histHMass_2e2m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histZMass_4e_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histZMass_4m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histZMass_2e2m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histZ2Mass_4e_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histZ2Mass_4m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histZ2Mass_2e2m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histPhiH_4e_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histPhiH_4m_fastsim->GetYaxis()->SetTitleOffset(1.4);
  histPhiH_2e2m_fastsim->GetYaxis()->SetTitleOffset(1.4);

  // Saving all histograms to file
  outRootFile.Write();

  // Plotting all histograms
  PlotTwo(histPhi0_4e_fullsim, histPhi0_4e_fastsim, (outputDirectory + "ValidateAngles_Phi0_4e" + label + ".gif").c_str(), error_option);
  PlotTwo(histPhi0_4m_fullsim, histPhi0_4m_fastsim, (outputDirectory + "ValidateAngles_Phi0_4m" + label + ".gif").c_str(), error_option);
  PlotTwo(histPhi0_2e2m_fullsim, histPhi0_2e2m_fastsim, (outputDirectory + "ValidateAngles_Phi0_2e2m" + label + ".gif").c_str(), error_option);
  PlotTwo(histTheta0_4e_fullsim, histTheta0_4e_fastsim, (outputDirectory + "ValidateAngles_Theta0_4e" + label + ".gif").c_str(), error_option);
  PlotTwo(histTheta0_4m_fullsim, histTheta0_4m_fastsim, (outputDirectory + "ValidateAngles_Theta0_4m" + label + ".gif").c_str(), error_option);
  PlotTwo(histTheta0_2e2m_fullsim, histTheta0_2e2m_fastsim, (outputDirectory + "ValidateAngles_Theta0_2e2m" + label + ".gif").c_str(), error_option);
  PlotTwo(histPhi_4e_fullsim, histPhi_4e_fastsim, (outputDirectory + "ValidateAngles_Phi_4e" + label + ".gif").c_str(), error_option);
  PlotTwo(histPhi_4m_fullsim, histPhi_4m_fastsim, (outputDirectory + "ValidateAngles_Phi_4m" + label + ".gif").c_str(), error_option);
  PlotTwo(histPhi_2e2m_fullsim, histPhi_2e2m_fastsim, (outputDirectory + "ValidateAngles_Phi_2e2m" + label + ".gif").c_str(), error_option);
  PlotTwo(histTheta1_4e_fullsim, histTheta1_4e_fastsim, (outputDirectory + "ValidateAngles_Theta1_4e" + label + ".gif").c_str(), error_option);
  PlotTwo(histTheta1_4m_fullsim, histTheta1_4m_fastsim, (outputDirectory + "ValidateAngles_Theta1_4m" + label + ".gif").c_str(), error_option);
  PlotTwo(histTheta1_2e2m_fullsim, histTheta1_2e2m_fastsim, (outputDirectory + "ValidateAngles_Theta1_2e2m" + label + ".gif").c_str(), error_option);
  PlotTwo(histTheta2_4e_fullsim, histTheta2_4e_fastsim, (outputDirectory + "ValidateAngles_Theta2_4e" + label + ".gif").c_str(), error_option);
  PlotTwo(histTheta2_4m_fullsim, histTheta2_4m_fastsim, (outputDirectory + "ValidateAngles_Theta2_4m" + label + ".gif").c_str(), error_option);
  PlotTwo(histTheta2_2e2m_fullsim, histTheta2_2e2m_fastsim, (outputDirectory + "ValidateAngles_Theta2_2e2m" + label + ".gif").c_str(), error_option);
  PlotTwo(histHMass_4e_fullsim, histHMass_4e_fastsim, (outputDirectory + "ValidateAngles_HMass_4e" + label + ".gif").c_str(), error_option);
  PlotTwo(histHMass_4m_fullsim, histHMass_4m_fastsim, (outputDirectory + "ValidateAngles_HMass_4m" + label + ".gif").c_str(), error_option);
  PlotTwo(histHMass_2e2m_fullsim, histHMass_2e2m_fastsim, (outputDirectory + "ValidateAngles_HMass_2e2m" + label + ".gif").c_str(), error_option);
  PlotTwo(histZMass_4e_fullsim, histZMass_4e_fastsim, (outputDirectory + "ValidateAngles_ZMass_4e" + label + ".gif").c_str(), error_option);
  PlotTwo(histZMass_4m_fullsim, histZMass_4m_fastsim, (outputDirectory + "ValidateAngles_ZMass_4m" + label + ".gif").c_str(), error_option);
  PlotTwo(histZMass_2e2m_fullsim, histZMass_2e2m_fastsim, (outputDirectory + "ValidateAngles_ZMass_2e2m" + label + ".gif").c_str(), error_option);
  PlotTwo(histZ2Mass_4e_fullsim, histZ2Mass_4e_fastsim, (outputDirectory + "ValidateAngles_Z2Mass_4e" + label + ".gif").c_str(), error_option);
  PlotTwo(histZ2Mass_4m_fullsim, histZ2Mass_4m_fastsim, (outputDirectory + "ValidateAngles_Z2Mass_4m" + label + ".gif").c_str(), error_option);
  PlotTwo(histZ2Mass_2e2m_fullsim, histZ2Mass_2e2m_fastsim, (outputDirectory + "ValidateAngles_Z2Mass_2e2m" + label + ".gif").c_str(), error_option);
  PlotTwo(histPhiH_4e_fullsim, histPhiH_4e_fastsim, (outputDirectory + "ValidateAngles_PhiH_4e" + label + ".gif").c_str(), error_option);
  PlotTwo(histPhiH_4m_fullsim, histPhiH_4m_fastsim, (outputDirectory + "ValidateAngles_PhiH_4m" + label + ".gif").c_str(), error_option);
  PlotTwo(histPhiH_2e2m_fullsim, histPhiH_2e2m_fastsim, (outputDirectory + "ValidateAngles_PhiH_2e2m" + label + ".gif").c_str(), error_option);
}
