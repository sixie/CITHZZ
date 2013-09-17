//================================================================================================
//			plotRegressionVarianceComparison_2D.C
// 
// Plots, for comparison, the deltaE / E obtained from the variance training versus the one obtained from the energy training
// To be used with trainElectronEnergyRegression_option.C and applyElectronEnergyRegression_option.C
//
// USAGE
// plotRegressionVarianceComparison_2D(string applyingFilename, string targetFilename1, string targetFilename2, string targetFilename3, string targetFilename4, 
//				       string versionOption1, string versionOption2, string versionOption3, string versionOption4,
//				       string outDirectoryName, Int_t Option)
//
// applyingFilename = file to which regression will be applied (must agree with file given to applyElectronEnergyRegression.C)
// targetFilename#  = target files outputted by applyElectronEnergyRegression.C
// versionOption#   = version ("V00" etc.) used for each target file outputted by applyElectronEnergyRegression.C
// outDirectoryName = name of the directory to which plot files (.gif and .root) will be outputted
// Option           = indicates which eta and pt bin to be plotted
//
// 0:	        eta < 0.8	 7 < pT < 10
// 1:	0.8   < eta < 1.479	 7 < pT < 10
// 2:	1.479 < eta  		 7 < pT < 10
// 3:	        eta < 0.8	10 < pT < 20
// 4:	0.8   < eta < 1.479	10 < pT < 20
// 5:	1.479 < eta  		10 < pT < 20
// 6:	        eta < 0.8	20 < pT
// 7:	0.8   < eta < 1.479	20 < pT
// 8:	1.479 < eta  		20 < pT
// 10:	any eta			7 < pT < 10
// -1:	any eta			7 < pT
//
// Bins based on tracker variables:
// 20:	barrel		brem < 0.2
// 21:	barrel		brem > 0.2	      EoP < 0.8
// 22:	barrel		brem > 0.2	0.8 < EoP < 1.2
// 23:	barrel		brem > 0.2	1.2 < EoP
//
// 30:	endcap		brem < 0.2
// 31:	endcap		brem > 0.2	      EoP < 0.8
// 32:	endcap		brem > 0.2	0.8 < EoP < 1.2
// 33:	endcap		brem > 0.2	1.2 < EoP
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TString.h>
#include <TGraphAsymmErrors.h>     
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "TLegend.h"
#include "TEfficiency.h"
#include "TProfile.h"
#include "CITCommon/CommonData/interface/ElectronTree.h"
#endif


//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void plotRegressionVarianceComparison_2D(string applyingFilename, string targetFilename1, string targetFilename2, string versionOption1, string versionOption2, string outDirectoryName, Int_t Option) {  

  // Setting up graphic style
  gStyle->SetOptStat(0);

  string targetFilenames_array[] = {targetFilename1, targetFilename2};
  string versionOptions_array[] = {versionOption1, versionOption2};
  vector<string> targetFilenames(targetFilenames_array, targetFilenames_array+2);
  vector<string> versionOptions(versionOptions_array, versionOptions_array+2);
  int nFiles = 2;

  // Preparing Option string to include in name
  string OptionStr; 
  ostringstream OptionConv;
  OptionConv << Option;
  OptionStr = OptionConv.str();

  vector<TH2F*> EleRegressionVariance_hists;		// Vector to store all the histograms

  cout << "Applying weights to file " << applyingFilename.c_str() << endl;

  for (int i = 0; i < nFiles; i++) {
    string targetFilename = targetFilenames[i];
    string versionOption  = versionOptions[i];

    cout << "Using target file " << targetFilename.c_str() << " with option " << versionOption.c_str() << " in bin " << Option << endl; 

    // Preparing strings for name
    string binNameStr;
    if (Option == 0)  binNameStr = "Energy regression variance, #eta < 0.8, 7 < p_{T} < 10";
    if (Option == 1)  binNameStr = "Energy regression variance, 0.8 < #eta < 1.479, 7 < p_{T} < 10";
    if (Option == 2)  binNameStr = "Energy regression variance, #eta > 1.479, 7 < p_{T} < 10";
    if (Option == 3)  binNameStr = "Energy regression variance, #eta < 0.8, 10 < p_{T} < 20";
    if (Option == 4)  binNameStr = "Energy regression variance, 0.8 < #eta < 1.479, 10 < p_{T} < 20";
    if (Option == 5)  binNameStr = "Energy regression variance, #eta > 1.479, 10 < p_{T} < 20";
    if (Option == 6)  binNameStr = "Energy regression variance, #eta < 0.8, p_{T} > 20";
    if (Option == 7)  binNameStr = "Energy regression variance, 0.8 < #eta < 1.479, p_{T} > 20";
    if (Option == 8)  binNameStr = "Energy regression variance, #eta > 1.479, p_{T} > 20";
    if (Option == 10) binNameStr = "Energy regression variance, 7 < p_{T} < 10";
    if (Option == -1) binNameStr = "Energy regression variance, p_{T} > 7";

    if (Option == 20) binNameStr = "Energy regression variance, barrel, brem < 0.2";
    if (Option == 21) binNameStr = "Energy regression variance, barrel, brem > 0.2, EoP < 0.8";
    if (Option == 22) binNameStr = "Energy regression variance, barrel, brem > 0.2, 0.8 < EoP < 1.2";
    if (Option == 23) binNameStr = "Energy regression variance, barrel, brem > 0.2, EoP > 1.2";

    if (Option == 30) binNameStr = "Energy regression variance, endcap, brem < 0.2";
    if (Option == 31) binNameStr = "Energy regression variance, endcap, brem > 0.2, EoP < 0.8";
    if (Option == 32) binNameStr = "Energy regression variance, endcap, brem > 0.2, 0.8 < EoP < 1.2";
    if (Option == 33) binNameStr = "Energy regression variance, endcap, brem > 0.2, EoP > 1.2";

    //--------------------------------------------------------------------------------------------------------------
    // Histograms
    //==============================================================================================================  
    TH2F *EleRegressionVariance = new TH2F(("EleRegressionVariance_" + versionOption).c_str(), (binNameStr + "; #Delta E / E (variance regression) ; #Delta E/ E (energy regression)").c_str(), 100, 0., 0.1, 100, 0., 0.1);

    // Adding histogram to vector
    EleRegressionVariance_hists.push_back(EleRegressionVariance);
    
    //*****************************************************************************************
    //RealEleTree
    //*****************************************************************************************
    citana::ElectronTree RealEleTree;
    RealEleTree.LoadTree(applyingFilename.c_str());

    // If versionOption is V00 or V01, there aren't any different treatment for low- and high-pT bins
    if (versionOption == "V00" || versionOption == "V01" || versionOption == "V5") {
      //*****************************************************************************************
      //Load Regression
      //*****************************************************************************************
      Float_t       fEleRegressionTargetEB;
      Float_t       fEleRegressionTargetEE;
      Float_t       fEleRegressionVarianceTargetEB;
      Float_t       fEleRegressionVarianceTargetEE;
      RealEleTree.InitTree();

      // Loading the weights
      TFile targetFile(targetFilename.c_str());
      TTree * targetee_tree = (TTree*) targetFile.Get("targetee_tree");
      TTree * targeteb_tree = (TTree*) targetFile.Get("targeteb_tree");
      TTree * targeteevar_tree = (TTree*) targetFile.Get("targeteevar_tree");
      TTree * targetebvar_tree = (TTree*) targetFile.Get("targetebvar_tree");

      targeteb_tree->SetBranchAddress( "targeteb", &fEleRegressionTargetEB); 
      targetee_tree->SetBranchAddress( "targetee", &fEleRegressionTargetEE); 
      targetebvar_tree->SetBranchAddress( "targetebvar", &fEleRegressionVarianceTargetEB); 
      targeteevar_tree->SetBranchAddress( "targeteevar", &fEleRegressionVarianceTargetEE); 

      for(UInt_t ientry=0; ientry < RealEleTree.tree_->GetEntries(); ientry++) {       	
	RealEleTree.tree_->GetEntry(ientry);
	targeteb_tree->GetEntry(ientry);
	targetee_tree->GetEntry(ientry);
	targetebvar_tree->GetEntry(ientry);
	targeteevar_tree->GetEntry(ientry);

	if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

	//don't evaluate performance using training events
	if (RealEleTree.fEventNumber % 2 == 0) continue;
	if (RealEleTree.fElePt < 7) continue;

	//classify by eta and pt bins
	Int_t subdet = 0;
	if (fabs(RealEleTree.fEleSCEta) < 0.8) subdet = 0;
	else if (fabs(RealEleTree.fEleSCEta) < 1.479) subdet = 1;
	else subdet = 2;

	Int_t ptBin = 0;
	if (RealEleTree.fElePt > 10.0) ptBin = 1;
	if (RealEleTree.fElePt > 20.0) ptBin = 2;

	// Also separating by fBrem and EoP
	Int_t fbremBin = 0;
	if (RealEleTree.fEleNBrem > 0.2) fbremBin = 1;

	Int_t EoPBin = 0;
	if (RealEleTree.fEleEOverP > 0.8) EoPBin = 1;
	if (RealEleTree.fEleEOverP > 1.2) EoPBin = 2;

	Bool_t passCuts = kFALSE;
	if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
	if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
	if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
	if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
	if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
	if (Option == 5) passCuts = (subdet == 2 && ptBin == 1);    
	if (Option == 6) passCuts = (subdet == 0 && ptBin == 2);
	if (Option == 7) passCuts = (subdet == 1 && ptBin == 2);
	if (Option == 8) passCuts = (subdet == 2 && ptBin == 2);    
	if (Option == 10) passCuts = (ptBin == 0 );
	if (Option == -1) passCuts = kTRUE;
	
	// Now for "category" bins
	if (Option == 20) passCuts = (subdet < 2 && fbremBin == 0);			// Barrel, low brem
	if (Option == 21) passCuts = (subdet < 2 && fbremBin == 1 && EoPBin == 0);	// Barrel, high brem and low EoP	
	if (Option == 22) passCuts = (subdet < 2 && fbremBin == 1 && EoPBin == 1);	// Barrel, high brem and middle EoP	
	if (Option == 23) passCuts = (subdet < 2 && fbremBin == 1 && EoPBin == 2);	// Barrel, high brem and high EoP	
	if (Option == 24) passCuts = (subdet < 2 && fbremBin == 1 && EoPBin != 1);	// Barrel, high brem and low+high EoP	

	if (Option == 30) passCuts = (subdet == 2 && fbremBin == 0);			// Endcap, low brem
	if (Option == 31) passCuts = (subdet == 2 && fbremBin == 1 && EoPBin == 0);	// Endcap, high brem and low EoP	
	if (Option == 32) passCuts = (subdet == 2 && fbremBin == 1 && EoPBin == 1);	// Endcap, high brem and middle EoP	
	if (Option == 33) passCuts = (subdet == 2 && fbremBin == 1 && EoPBin == 2);	// Endcap, high brem and high EoP	
	if (Option == 34) passCuts = (subdet == 2 && fbremBin == 1 && EoPBin != 1);	// Endcap, high brem and low+high EoP	

	if (!passCuts) continue;    

	// Debugging
	//cout << "Debug barrel: " << fEleRegressionTargetEB << " : " << fEleRegressionVarianceTargetEB << " : " << RealEleTree.fGeneratedEnergyStatus1 << " : " << RealEleTree.fSCRawEnergy << endl;

	// Delta E over E from the energy regression training
	Double_t EoE_energy;
	if (fabs(RealEleTree.fEleSCEta) < 1.479) {
	  EoE_energy = fabs((fEleRegressionTargetEB * RealEleTree.fSCRawEnergy - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1);
	} else {
	  EoE_energy = fabs((fEleRegressionTargetEE * RealEleTree.fSCRawEnergy*(1+RealEleTree.fElePreShowerOverRaw) - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1);
	}

	// Delta E over E from the variance regression
	Double_t EoE_variance;
	if (fabs(RealEleTree.fEleSCEta) < 1.479) {
	  EoE_variance = (fEleRegressionVarianceTargetEB / fEleRegressionTargetEB);
	} else {
	  EoE_variance = (fEleRegressionVarianceTargetEE / fEleRegressionTargetEE);
	}

	EleRegressionVariance->Fill(EoE_variance, EoE_energy, RealEleTree.fWeight);
      } 
    }

    // Now if version is V10 or V11
    if (versionOption == "V10" || versionOption == "V11") {
      //*****************************************************************************************
      //Load Regression
      //*****************************************************************************************
      Float_t       fEleRegressionTargetEB_lowPt;
      Float_t       fEleRegressionTargetEB_highPt;
      Float_t       fEleRegressionTargetEE_lowPt;
      Float_t       fEleRegressionTargetEE_highPt;
      Float_t       fEleRegressionVarianceTargetEB_lowPt;
      Float_t       fEleRegressionVarianceTargetEB_highPt;
      Float_t       fEleRegressionVarianceTargetEE_lowPt;
      Float_t       fEleRegressionVarianceTargetEE_highPt;
      RealEleTree.InitTree();

      // Loading the weights
      TFile targetFile(targetFilename.c_str());
      TTree * targetee_lowPt_tree = (TTree*) targetFile.Get("targetee_lowPt_tree");
      TTree * targetee_highPt_tree = (TTree*) targetFile.Get("targetee_highPt_tree");
      TTree * targeteb_lowPt_tree = (TTree*) targetFile.Get("targeteb_lowPt_tree");
      TTree * targeteb_highPt_tree = (TTree*) targetFile.Get("targeteb_highPt_tree");
      TTree * targeteevar_lowPt_tree = (TTree*) targetFile.Get("targeteevar_lowPt_tree");
      TTree * targeteevar_highPt_tree = (TTree*) targetFile.Get("targeteevar_highPt_tree");
      TTree * targetebvar_lowPt_tree = (TTree*) targetFile.Get("targetebvar_lowPt_tree");
      TTree * targetebvar_highPt_tree = (TTree*) targetFile.Get("targetebvar_highPt_tree");

      targeteb_lowPt_tree->SetBranchAddress( "targeteb_lowPt", &fEleRegressionTargetEB_lowPt); 
      targeteb_highPt_tree->SetBranchAddress( "targeteb_highPt", &fEleRegressionTargetEB_highPt); 
      targetee_lowPt_tree->SetBranchAddress( "targetee_lowPt", &fEleRegressionTargetEE_lowPt); 
      targetee_highPt_tree->SetBranchAddress( "targetee_highPt", &fEleRegressionTargetEE_highPt); 
      targetebvar_lowPt_tree->SetBranchAddress( "targetebvar_lowPt", &fEleRegressionVarianceTargetEB_lowPt); 
      targetebvar_highPt_tree->SetBranchAddress( "targetebvar_highPt", &fEleRegressionVarianceTargetEB_highPt); 
      targeteevar_lowPt_tree->SetBranchAddress( "targeteevar_lowPt", &fEleRegressionVarianceTargetEE_lowPt); 
      targeteevar_highPt_tree->SetBranchAddress( "targeteevar_highPt", &fEleRegressionVarianceTargetEE_highPt); 

      for(UInt_t ientry=0; ientry < RealEleTree.tree_->GetEntries(); ientry++) {       	
	RealEleTree.tree_->GetEntry(ientry);
	targeteb_lowPt_tree->GetEntry(ientry);
	targeteb_highPt_tree->GetEntry(ientry);
	targetee_lowPt_tree->GetEntry(ientry);
	targetee_highPt_tree->GetEntry(ientry);
	targetebvar_lowPt_tree->GetEntry(ientry);
	targetebvar_highPt_tree->GetEntry(ientry);
	targeteevar_lowPt_tree->GetEntry(ientry);
	targeteevar_highPt_tree->GetEntry(ientry);

	if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

	//don't evaluate performance using training events
	if (RealEleTree.fEventNumber % 2 == 0) continue;
	if (RealEleTree.fElePt < 7) continue;

	//classify by eta and pt bins
	Int_t subdet = 0;
	if (fabs(RealEleTree.fEleSCEta) < 0.8) subdet = 0;
	else if (fabs(RealEleTree.fEleSCEta) < 1.479) subdet = 1;
	else subdet = 2;

	Int_t ptBin = 0;
	if (RealEleTree.fElePt > 10.0) ptBin = 1;
	if (RealEleTree.fElePt > 20.0) ptBin = 2;

	// Also separating by fBrem and EoP
	Int_t fbremBin = 0;
	if (RealEleTree.fEleNBrem > 0.2) fbremBin = 1;

	Int_t EoPBin = 0;
	if (RealEleTree.fEleEOverP > 0.8) EoPBin = 1;
	if (RealEleTree.fEleEOverP > 1.2) EoPBin = 2;

	Bool_t passCuts = kFALSE;
	if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
	if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
	if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
	if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
	if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
	if (Option == 5) passCuts = (subdet == 2 && ptBin == 1);    
	if (Option == 6) passCuts = (subdet == 0 && ptBin == 2);
	if (Option == 7) passCuts = (subdet == 1 && ptBin == 2);
	if (Option == 8) passCuts = (subdet == 2 && ptBin == 2);    
	if (Option == 10) passCuts = (ptBin == 0 );
	if (Option == -1) passCuts = kTRUE;

	// Now for "category" bins
	if (Option == 20) passCuts = (subdet < 2 && fbremBin == 0);			// Barrel, low brem
	if (Option == 21) passCuts = (subdet < 2 && fbremBin == 1 && EoPBin == 0);	// Barrel, high brem and low EoP	
	if (Option == 22) passCuts = (subdet < 2 && fbremBin == 1 && EoPBin == 1);	// Barrel, high brem and middle EoP	
	if (Option == 23) passCuts = (subdet < 2 && fbremBin == 1 && EoPBin == 2);	// Barrel, high brem and high EoP	

	if (Option == 30) passCuts = (subdet == 2 && fbremBin == 0);			// Endcap, low brem
	if (Option == 31) passCuts = (subdet == 2 && fbremBin == 1 && EoPBin == 0);	// Endcap, high brem and low EoP	
	if (Option == 32) passCuts = (subdet == 2 && fbremBin == 1 && EoPBin == 1);	// Endcap, high brem and middle EoP	
	if (Option == 33) passCuts = (subdet == 2 && fbremBin == 1 && EoPBin == 2);	// Endcap, high brem and high EoP	

	if (!passCuts) continue;    

	// Delta E over E from the energy regression training
	Double_t EoE_energy;
	if (RealEleTree.fElePt <= 15.0) {						// For pt < 15, use pt instead of SCRawEnergy
	  if (fabs(RealEleTree.fEleSCEta) < 1.479) {
	    EoE_energy = fabs((fEleRegressionTargetEB_lowPt * RealEleTree.fElePt * TMath::CosH(RealEleTree.fEleEta) - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1);
	  } else {
	    EoE_energy = fabs((fEleRegressionTargetEE_lowPt * RealEleTree.fElePt * TMath::CosH(RealEleTree.fEleEta) - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1);
	  }
	}
	
	else {
	  if (fabs(RealEleTree.fEleSCEta) < 1.479) {
	    EoE_energy = fabs((fEleRegressionTargetEB_highPt * RealEleTree.fSCRawEnergy - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1);
	  } else {
	    EoE_energy = fabs((fEleRegressionTargetEE_highPt * RealEleTree.fSCRawEnergy*(1+RealEleTree.fElePreShowerOverRaw) - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1);
	  }
	}

	// Delta E over E from the variance regression
	Double_t EoE_variance;
	if (RealEleTree.fElePt <= 15.0) {						// For pt < 15, use pt instead of SCRawEnergy
	  if (fabs(RealEleTree.fEleSCEta) < 1.479) {
	    EoE_variance = fabs(fEleRegressionVarianceTargetEB_lowPt / fEleRegressionTargetEB_lowPt);
	  } else {
	    EoE_variance = fabs(fEleRegressionVarianceTargetEB_highPt / fEleRegressionTargetEB_highPt);
	  }
	}
	
	else {
	  if (fabs(RealEleTree.fEleSCEta) < 1.479) {
	    EoE_variance = fabs(fEleRegressionVarianceTargetEE_lowPt / fEleRegressionTargetEE_lowPt);
	  } else {
	    EoE_variance = fabs(fEleRegressionVarianceTargetEE_highPt / fEleRegressionTargetEE_highPt);
	  }
	}

	EleRegressionVariance->Fill(EoE_variance, EoE_energy, RealEleTree.fWeight);
      } 
    }
  }
  TCanvas *cv = 0;
  TLegend *legend = 0;

  // Plotting histogram
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  EleRegressionVariance_hists[0]->Draw("colz");
  cv->SaveAs((outDirectoryName + "EnergyVariance2D_" + versionOptions[0] + "_binOption" + OptionStr + ".gif").c_str());
  cv->SaveAs((outDirectoryName + "EnergyVariance2D_" + versionOptions[0] + "_binOption" + OptionStr + ".pdf").c_str());
  cv->SaveAs((outDirectoryName + "EnergyVariance2D_" + versionOptions[0] + "_binOption" + OptionStr + ".root").c_str());

  cv = new TCanvas("cv", "cv", 800, 600);
  EleRegressionVariance_hists[1]->Draw("colz");
  cv->SaveAs((outDirectoryName + "EnergyVariance2D_" + versionOptions[1] + "_binOption" + OptionStr + ".gif").c_str());
  cv->SaveAs((outDirectoryName + "EnergyVariance2D_" + versionOptions[1] + "_binOption" + OptionStr + ".pdf").c_str());
  cv->SaveAs((outDirectoryName + "EnergyVariance2D_" + versionOptions[1] + "_binOption" + OptionStr + ".root").c_str());
} 

