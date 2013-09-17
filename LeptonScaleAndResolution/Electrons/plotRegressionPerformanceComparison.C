//================================================================================================
//			plotRegressionPerformance.C
// 
// Plots, for comparison, the performance of four different regressions for electron energy
// To be used with trainElectronEnergyRegression_option.C and applyElectronEnergyRegression_option.C
//
// USAGE
// plotRegressionPerformanceComparison(string applyingFilename, string targetFilename1, string targetFilename2, string targetFilename3, string targetFilename4, 
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

#include "CITCommon/CommonData/interface/ElectronTree.h"
#endif


//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
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
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void plotRegressionPerformanceComparison(string applyingFilename, string targetFilename1, string targetFilename2, string targetFilename3, string versionOption1, string versionOption2, string versionOption3, string outDirectoryName, Int_t Option) {  

  // Setting up graphic style
  gStyle->SetOptStat(0);

  string targetFilenames_array[] = {targetFilename1, targetFilename2, targetFilename3};
  string versionOptions_array[] = {versionOption1, versionOption2, versionOption3};
  vector<string> targetFilenames(targetFilenames_array, targetFilenames_array+3);
  vector<string> versionOptions(versionOptions_array, versionOptions_array+3);
  int nFiles = 3;

  // Preparing Option string to include in name
  string OptionStr; 
  ostringstream OptionConv;
  OptionConv << Option;
  OptionStr = OptionConv.str();

  vector<TH1F*> EleRegressionEnergyOverGenEnergy_hists;		// Vector to store all the histograms
  TH1F* EleEnergyOverGenEnergy;

  cout << "Applying weights to file " << applyingFilename.c_str() << endl;

  for (int i = 0; i < nFiles; i++) {
    string targetFilename = targetFilenames[i];
    string versionOption  = versionOptions[i];

    cout << "Using target file " << targetFilename.c_str() << " with option " << versionOption.c_str() << " in bin " << Option << endl; 

    // Preparing strings for name
    string binNameStr;
    if (Option == 0)  binNameStr = "Electron energy ratio, #eta < 0.8, 7 < p_{T} < 10";
    if (Option == 1)  binNameStr = "Electron energy ratio, 0.8 < #eta < 1.479, 7 < p_{T} < 10";
    if (Option == 2)  binNameStr = "Electron energy ratio, #eta > 1.479, 7 < p_{T} < 10";
    if (Option == 3)  binNameStr = "Electron energy ratio, #eta < 0.8, 10 < p_{T} < 20";
    if (Option == 4)  binNameStr = "Electron energy ratio, 0.8 < #eta < 1.479, 10 < p_{T} < 20";
    if (Option == 5)  binNameStr = "Electron energy ratio, #eta > 1.479, 10 < p_{T} < 20";
    if (Option == 6)  binNameStr = "Electron energy ratio, #eta < 0.8, p_{T} > 20";
    if (Option == 7)  binNameStr = "Electron energy ratio, 0.8 < #eta < 1.479, p_{T} > 20";
    if (Option == 8)  binNameStr = "Electron energy ratio, #eta > 1.479, p_{T} > 20";
    if (Option == 10) binNameStr = "Electron energy ratio, 7 < p_{T} < 10";
    if (Option == -1) binNameStr = "Electron energy ratio, p_{T} > 7";

    if (Option == 20) binNameStr = "Electron energy ratio, barrel, brem < 0.2";
    if (Option == 21) binNameStr = "Electron energy ratio, barrel, brem > 0.2, EoP < 0.8";
    if (Option == 22) binNameStr = "Electron energy ratio, barrel, brem > 0.2, 0.8 < EoP < 1.2";
    if (Option == 23) binNameStr = "Electron energy ratio, barrel, brem > 0.2, EoP > 1.2";

    if (Option == 30) binNameStr = "Electron energy ratio, endcap, brem < 0.2";
    if (Option == 31) binNameStr = "Electron energy ratio, endcap, brem > 0.2, EoP < 0.8";
    if (Option == 32) binNameStr = "Electron energy ratio, endcap, brem > 0.2, 0.8 < EoP < 1.2";
    if (Option == 33) binNameStr = "Electron energy ratio, endcap, brem > 0.2, EoP > 1.2";

    //--------------------------------------------------------------------------------------------------------------
    // Histograms
    //==============================================================================================================  
    EleEnergyOverGenEnergy = new TH1F("EleEnergyOverGenEnergy", (binNameStr + "; Electron Energy / Generated Energy ; Number of Events ").c_str(),  150, 0., 1.5);
    TH1F *EleRegressionEnergyOverGenEnergy = new TH1F(("EleRegressionEnergyOverGenEnergy_" + versionOption).c_str(), "; Electron Regression Energy / Generated Energy ; Number of Events ", 150, 0., 1.5);

    // Adding histogram to vector
    EleRegressionEnergyOverGenEnergy_hists.push_back(EleRegressionEnergyOverGenEnergy);
    
    //*****************************************************************************************
    //RealEleTree
    //*****************************************************************************************
    citana::ElectronTree RealEleTree;
    RealEleTree.LoadTree(applyingFilename.c_str());

    // If versionOption is V00 or V01, there aren't any different treatment for low- and high-pT bins
    if (versionOption == "V00" || versionOption == "V01" || versionOption == "V02") {
      //*****************************************************************************************
      //Load Regression
      //*****************************************************************************************
      Float_t       fEleRegressionTargetEB;
      Float_t       fEleRegressionTargetEE;
      RealEleTree.InitTree();

      // Loading the weights
      TFile targetFile(targetFilename.c_str());
      TTree * targetee_tree = (TTree*) targetFile.Get("targetee_tree");
      TTree * targeteb_tree = (TTree*) targetFile.Get("targeteb_tree");

      targeteb_tree->SetBranchAddress( "targeteb", &fEleRegressionTargetEB); 
      targetee_tree->SetBranchAddress( "targetee", &fEleRegressionTargetEE); 

      for(UInt_t ientry=0; ientry < RealEleTree.tree_->GetEntries(); ientry++) {       	
	RealEleTree.tree_->GetEntry(ientry);
	targeteb_tree->GetEntry(ientry);
	targetee_tree->GetEntry(ientry);

	if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

	//don't evaluate performance using training events
	if (RealEleTree.fEventNumber % 2 == 0) continue;
	if (RealEleTree.fElePt < 7) continue;
        if (!(RealEleTree.fGeneratedEnergyStatus3 >= RealEleTree.fGeneratedEnergyStatus1 && (RealEleTree.fGeneratedEnergyStatus3 - RealEleTree.fGeneratedEnergyStatus1)/RealEleTree.fGeneratedEnergyStatus3 < 0.01)) continue;

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

	// Regression energy
	Double_t Energy;
	if (fabs(RealEleTree.fEleSCEta) < 1.479) {
	  Energy = fEleRegressionTargetEB * RealEleTree.fSCRawEnergy;
	} else {
	  Energy = fEleRegressionTargetEE * ( RealEleTree.fSCRawEnergy * (1 + RealEleTree.fElePreShowerOverRaw ) );
	}

	EleEnergyOverGenEnergy->Fill(RealEleTree.fElePt * TMath::CosH(RealEleTree.fEleEta) / RealEleTree.fGeneratedEnergyStatus1 ,RealEleTree.fWeight);
	EleRegressionEnergyOverGenEnergy->Fill(Energy / RealEleTree.fGeneratedEnergyStatus1 ,RealEleTree.fWeight);
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
      RealEleTree.InitTree();

      // Loading the weights
      TFile targetFile(targetFilename.c_str());
      TTree * targetee_lowPt_tree = (TTree*) targetFile.Get("targetee_lowPt_tree");
      TTree * targetee_highPt_tree = (TTree*) targetFile.Get("targetee_highPt_tree");
      TTree * targeteb_lowPt_tree = (TTree*) targetFile.Get("targeteb_lowPt_tree");
      TTree * targeteb_highPt_tree = (TTree*) targetFile.Get("targeteb_highPt_tree");

      targeteb_lowPt_tree->SetBranchAddress( "targeteb_lowPt", &fEleRegressionTargetEB_lowPt); 
      targeteb_highPt_tree->SetBranchAddress( "targeteb_highPt", &fEleRegressionTargetEB_highPt); 
      targetee_lowPt_tree->SetBranchAddress( "targetee_lowPt", &fEleRegressionTargetEE_lowPt); 
      targetee_highPt_tree->SetBranchAddress( "targetee_highPt", &fEleRegressionTargetEE_highPt); 

      for(UInt_t ientry=0; ientry < RealEleTree.tree_->GetEntries(); ientry++) {       	
	RealEleTree.tree_->GetEntry(ientry);
	targeteb_lowPt_tree->GetEntry(ientry);
	targeteb_highPt_tree->GetEntry(ientry);
	targetee_lowPt_tree->GetEntry(ientry);
	targetee_highPt_tree->GetEntry(ientry);

	if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

	//don't evaluate performance using training events
	if (RealEleTree.fEventNumber % 2 == 0) continue;
	if (RealEleTree.fElePt < 7) continue;
        if (!(RealEleTree.fGeneratedEnergyStatus3 >= RealEleTree.fGeneratedEnergyStatus1 && (RealEleTree.fGeneratedEnergyStatus3 - RealEleTree.fGeneratedEnergyStatus1)/RealEleTree.fGeneratedEnergyStatus3 < 0.01)) continue;

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

	Double_t Energy;
	if (RealEleTree.fElePt <= 15.0) {						// For pt < 15, use pt instead of SCRawEnergy
	  if (fabs(RealEleTree.fEleSCEta) < 1.479) {
	    Energy = fEleRegressionTargetEB_lowPt * RealEleTree.fElePt * TMath::CosH(RealEleTree.fEleEta);
	  } else {
	    Energy = fEleRegressionTargetEE_lowPt * RealEleTree.fElePt * TMath::CosH(RealEleTree.fEleEta);
	  }
	}

	if (RealEleTree.fElePt > 15.0) {
	  if (fabs(RealEleTree.fEleSCEta) < 1.479) {
	    Energy = fEleRegressionTargetEB_highPt * RealEleTree.fSCRawEnergy;
	  } else {
	    Energy = fEleRegressionTargetEE_highPt * ( RealEleTree.fSCRawEnergy * (1 + RealEleTree.fElePreShowerOverRaw ) );
	  }
	}

	EleEnergyOverGenEnergy->Fill(RealEleTree.fElePt * TMath::CosH(RealEleTree.fEleEta) / RealEleTree.fGeneratedEnergyStatus1 ,RealEleTree.fWeight);
	EleRegressionEnergyOverGenEnergy->Fill(Energy / RealEleTree.fGeneratedEnergyStatus1 ,RealEleTree.fWeight);
      } 
    }
  }
  TCanvas *cv = 0;
  TLegend *legend = 0;

  // Plotting all for comparison
  // First, finding maximum Y range
  float max_yrange = 0;
  for (int j = 0; j < nFiles; j++) {
    float max = EleRegressionEnergyOverGenEnergy_hists[j]->GetMaximum();
    if (max > max_yrange) max_yrange = max;
  }
  float max0 = EleEnergyOverGenEnergy->GetMaximum();
  if (max0 > max_yrange) max_yrange = max0;

  // Settings for histograms
  EleEnergyOverGenEnergy->GetYaxis()->SetTitleOffset(1.5);
  EleEnergyOverGenEnergy->SetMaximum(max_yrange * 1.25);
  EleEnergyOverGenEnergy->SetLineColor(kBlack);
  EleRegressionEnergyOverGenEnergy_hists[0]->SetLineColor(kRed);
  EleRegressionEnergyOverGenEnergy_hists[1]->SetLineColor(kGreen);
  EleRegressionEnergyOverGenEnergy_hists[2]->SetLineColor(kBlue);
  EleEnergyOverGenEnergy->SetLineWidth(3);
  EleRegressionEnergyOverGenEnergy_hists[0]->SetLineWidth(2);
  EleRegressionEnergyOverGenEnergy_hists[1]->SetLineWidth(2);
  EleRegressionEnergyOverGenEnergy_hists[2]->SetLineWidth(2);

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(EleEnergyOverGenEnergy,"Standard", "L");
  EleEnergyOverGenEnergy->Draw("hist");
  for (int j = 0; j < nFiles; j++){
    // Prepare string for version name
    string versionName;
    if (versionOptions[j] == "V00") versionName = "ECAL only"; 
    if (versionOptions[j] == "V01") versionName = "ECAL + Trk Vars"; 
    if (versionOptions[j] == "V02") versionName = "ECAL + More Trk Vars"; 
    legend->AddEntry(EleRegressionEnergyOverGenEnergy_hists[j], ("Regression, " + versionName).c_str(), "L");
    EleRegressionEnergyOverGenEnergy_hists[j]->Draw("histsame");
  }
  legend->Draw();
 
  cv->SaveAs((outDirectoryName + "EnergyResponse_allversions_binOption" + OptionStr + ".gif").c_str());
//  cv->SaveAs((outDirectoryName + "EnergyResponse_allversions_binOption" + OptionStr + ".pdf").c_str());
//  cv->SaveAs((outDirectoryName + "EnergyResponse_allversions_binOption" + OptionStr + ".root").c_str());
} 

