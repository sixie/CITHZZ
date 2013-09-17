//================================================================================================
//			plotRegressionVarianceComparison.C
// 
// Plots, for comparison, the performance of four different regressions for electron energy
// To be used with trainElectronEnergyRegression_option.C and applyElectronEnergyRegression_option.C
//
// USAGE
// plotRegressionVarianceComparison(string applyingFilename, string targetFilename1, string targetFilename2, string targetFilename3, string targetFilename4, 
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
#include "TF1.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TProfile.h"
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
    hist->SetBinContent(b,(2.50663/0.14) * hist->GetBinContent(b) / norm);
    hist->SetBinError(b,(2.50663/0.14) * hist->GetBinError(b) / norm);
  } 

  norm = 0;
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  cout << "norm: " << norm << endl;

  return;
}

void ValidateRegressionVariance(string filename) {

  TFile *file = new TFile(filename.c_str(), "READ");
//   vector<TH1F*> EnergyVariance;
//   vector<TH1F*> EnergyResponse;

  for (int i=6; i < 7; ++i) {    
    TH1F *tmpEleEnergyVariance = (TH1F*)file->Get(Form("EleEnergyVariance_V00_Bin%d",i));
    TH1F *tmpEleEnergyResponse = (TH1F*)file->Get(Form("EleEnergyResponse_V00_Bin%d",i));
    TH1F *tmpEleRegressionEnergyNormalizedResponse = (TH1F*)file->Get(Form("EleRegressionEnergyNormalizedResponse_V00_Bin%d",i));
//     EnergyVariance.push_back(tmpEleEnergyVariance);
//     EnergyResponse.push_back(tmpEleEnergyResponse);
    NormalizeHist(tmpEleRegressionEnergyNormalizedResponse);

    tmpEleEnergyResponse->Fit("gaus", "","", -0.1, 0.1);

    cout << endl;
    cout << "Bin " << i << endl;
    cout << "Mean of Variance: " << tmpEleEnergyVariance->GetMean() << " +/- " << tmpEleEnergyVariance->GetMeanError() << endl;
//     cout << "RMS of Reponse: " << tmpEleEnergyResponse->GetRMS() << " +/- " << tmpEleEnergyResponse->GetRMSError() << endl;

//     TF1 *unitGauss = new TF1("unitGauss", "gaus", -5,5);
//     unitGauss->SetLineColor(kRed);
//     unitGauss->SetParameters(1,0,1);
//     cout << "Int: " << unitGauss->Integral(-10,10);

    TCanvas *cv = 0;
    TLegend *legend = 0;

    cv = new TCanvas("cv2", "cv2", 800, 600);
    legend = new TLegend(0.1,0.7,0.3,0.9);
    legend->SetTextSize(0.03);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(tmpEleRegressionEnergyNormalizedResponse, "Regression Response Pull", "L");
//     legend->AddEntry(unitGauss, "Unit Gaussian", "L");
    tmpEleRegressionEnergyNormalizedResponse->Draw("hist");
//     unitGauss->Draw("Lsame");
    legend->Draw();
    cv->SaveAs(Form("RegressionResponsePull_Bin%d.gif",i));
    

  }

 

}


//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void makeRegressionVarianceComparison(string applyingFilename, string targetFilename1, string targetFilename2, string versionOption1, string versionOption2, string outDirectoryName, Int_t Option) {  

  // Setting up graphic style
  gStyle->SetOptStat(0);

  string targetFilenames_array[] = {targetFilename1, targetFilename2};
  string versionOptions_array[] = {versionOption1, versionOption2};
  vector<string> targetFilenames(targetFilenames_array, targetFilenames_array+2);
  vector<string> versionOptions(versionOptions_array, versionOptions_array+2);
  int nFiles = 1;

  // Preparing Option string to include in name
  string OptionStr; 
  ostringstream OptionConv;
  OptionConv << Option;
  OptionStr = OptionConv.str();

  vector<TH2F*> EleEnergyError_RegressionVsTrue_hists;		// Vector to store all the histograms
  vector<TH2F*> EleEnergyErrorVsEt_hists;	    	        // Vector to store all the histograms
  vector<TH1F*> EleEnergyVariance_hists;		        // Vector to store all the histograms
  vector<TH1F*> EleEnergyResponse_hists;		        // Vector to store all the histograms
  vector<TH2F*> EleEnergyOverETrueVsEOverP_hists;		// Vector to store all the histograms
  vector<TH2F*> ElePOverETrueVsEOverP_hists;		        // Vector to store all the histograms

  vector<TH1F*> EleRegressionEnergyNormalizedResponse_hists;	// Vector to store all the histograms
  vector<TH2F*> EleEnergyResponseStandardVsRegression_hists;	// Vector to store all the histograms
  vector<TH2F*> EleEnergyResponseStandardVsRegression_Better_hists;	// Vector to store all the histograms
  vector<TH2F*> EleEnergyResponseStandardVsRegression_Worse_hists;	// Vector to store all the histograms
  vector<TH1F*> EleEnergyMigration_Better_hists;	// Vector to store all the histograms
  vector<TH1F*> EleEnergyMigration_Worse_hists;	// Vector to store all the histograms

  cout << "Applying weights to file " << applyingFilename.c_str() << endl;
  Double_t TotalElectrons = 0;
  Double_t RegressionImprovedElectrons = 0;
  Double_t RegressionDegradedElectrons = 0;


  for (int i = 0; i < nFiles; i++) {
    string targetFilename = targetFilenames[i];
    string versionOption  = versionOptions[i];

    cout << "Using target file " << targetFilename.c_str() << " with option " << versionOption.c_str() << " in bin " << Option << endl; 

    // Preparing strings for name
    string binNameStr;
    if (Option == 0)  binNameStr = "Energy regression variance profile, #eta < 0.8, 7 < p_{T} < 10";
    if (Option == 1)  binNameStr = "Energy regression variance profile, 0.8 < #eta < 1.479, 7 < p_{T} < 10";
    if (Option == 2)  binNameStr = "Energy regression variance profile, #eta > 1.479, 7 < p_{T} < 10";
    if (Option == 3)  binNameStr = "Energy regression variance profile, #eta < 0.8, 10 < p_{T} < 20";
    if (Option == 4)  binNameStr = "Energy regression variance profile, 0.8 < #eta < 1.479, 10 < p_{T} < 20";
    if (Option == 5)  binNameStr = "Energy regression variance profile, #eta > 1.479, 10 < p_{T} < 20";
    if (Option == 6)  binNameStr = "Energy regression variance profile, #eta < 0.8, p_{T} > 20";
    if (Option == 7)  binNameStr = "Energy regression variance profile, 0.8 < #eta < 1.479, p_{T} > 20";
    if (Option == 8)  binNameStr = "Energy regression variance profile, #eta > 1.479, p_{T} > 20";
    if (Option == 10) binNameStr = "Energy regression variance profile, 7 < p_{T} < 10";
    if (Option == 11) binNameStr = "Energy regression variance profile, 10 < p_{T} < 20";
    if (Option == 12) binNameStr = "Energy regression variance profile, p_{T} > 20";
    if (Option == -3) binNameStr = "Energy regression variance profile, Endcap";
    if (Option == -2) binNameStr = "Energy regression variance profile, Barrel";
    if (Option == -1) binNameStr = "Energy regression variance profile, p_{T} > 7";

    if (Option == 20) binNameStr = "Energy regression variance profile, barrel, brem < 0.2";
    if (Option == 21) binNameStr = "Energy regression variance profile, barrel, brem > 0.2, EoP < 0.8";
    if (Option == 22) binNameStr = "Energy regression variance profile, barrel, brem > 0.2, 0.8 < EoP < 1.2";
    if (Option == 23) binNameStr = "Energy regression variance profile, barrel, brem > 0.2, EoP > 1.2";

    if (Option == 30) binNameStr = "Energy regression variance profile, endcap, brem < 0.2";
    if (Option == 31) binNameStr = "Energy regression variance profile, endcap, brem > 0.2, EoP < 0.8";
    if (Option == 32) binNameStr = "Energy regression variance profile, endcap, brem > 0.2, 0.8 < EoP < 1.2";
    if (Option == 33) binNameStr = "Energy regression variance profile, endcap, brem > 0.2, EoP > 1.2";

    //--------------------------------------------------------------------------------------------------------------
    // Histograms
    //==============================================================================================================  
    TH2F *EleEnergyError_RegressionVsTrue = new TH2F(("EleEnergyError_RegressionVsTrue_" + versionOption).c_str(), (binNameStr + "; #Delta E / E (true) ; #Delta E/ E (regression)").c_str(), 100, 0., 0.1, 100, 0., 0.1);
    TH2F *EleEnergyErrorVsEt = new TH2F(("EleEnergyErrorVsEt_" + versionOption).c_str(), (binNameStr + "; E_{t} [GeV] ; #Delta E/ E (regression)").c_str(), 100, 0., 100, 100, 0., 0.1);
    TH1F *EleEnergyVariance = new TH1F(("EleEnergyVariance_" + versionOption).c_str(), (binNameStr + "; #Delta E / E ; Number of Events ").c_str(), 200, 0., 1.);
    TH1F *EleEnergyResponse = new TH1F(("EleEnergyResponse_" + versionOption).c_str(), (binNameStr + "; (E_{reco} - E_{true}) / E_{true} ; Number of Events ").c_str(), 200, -1., 1.);

    TH2F *EleEnergyOverETrueVsEOverP = new TH2F(("EleEnergyOverETrueVsEOverP_" + versionOption).c_str(), (binNameStr + "; E_{regression}/P ; E_{regression} / E_{true}").c_str(), 500, 0, 2.5, 500, 0., 2.5);
    TH2F *ElePOverETrueVsEOverP = new TH2F(("ElePOverETrueVsEOverP_" + versionOption).c_str(), (binNameStr + "; E_{regression}/P ; P / E_{true}").c_str(), 500, 0, 2.5, 500, 0., 2.5);

    TH1F *EleRegressionEnergyNormalizedResponse = new TH1F(("EleRegressionEnergyNormalizedResponse_" + versionOption).c_str(), (binNameStr + "; (E_{reco} - E_{true}) / #Delta E_{reco} ; Number of Events ").c_str(), 100, -7., 7.);
    TH2F *EleEnergyResponseStandardVsRegression = new TH2F(("EleEnergyResponseStandardVsRegression_" + versionOption).c_str(), (binNameStr + "; (E_{SC} - E_{true}) / E_{true}  ; (E_{regression} - E_{true}) / E_{true} ").c_str(), 500, -1, 1, 500, -1, 1);
    TH2F *EleEnergyResponseStandardVsRegression_Better = new TH2F(("EleEnergyResponseStandardVsRegression_Better_" + versionOption).c_str(), (binNameStr + "; (E_{SC} - E_{true}) / E_{true}  ; (E_{regression} - E_{true}) / E_{true} ").c_str(), 500, -1, 1, 500, -1, 1);
    TH2F *EleEnergyResponseStandardVsRegression_Worse = new TH2F(("EleEnergyResponseStandardVsRegression_Worse_" + versionOption).c_str(), (binNameStr + "; (E_{SC} - E_{true}) / E_{true}  ; (E_{regression} - E_{true}) / E_{true} ").c_str(), 500, -1, 1, 500, -1, 1);
    TH1F *EleEnergyMigration_Better = new TH1F(("EleEnergyMigration_Better_" + versionOption).c_str(), (binNameStr + "; |E_{corr SC} - E_{true}|/E_{true} - |E_{regression} - E_{true}|/E_{true} ; Number of Events ").c_str(), 100, 0.0, 0.5);
    TH1F *EleEnergyMigration_Worse = new TH1F(("EleEnergyMigration_Worse_" + versionOption).c_str(), (binNameStr + "; |E_{corr SC} - E_{true}|/E_{true} - |E_{regression} - E_{true}|/E_{true} ; Number of Events ").c_str(), 100, 0.0, 0.5);


   // Adding histogram to vector
    EleEnergyError_RegressionVsTrue_hists.push_back(EleEnergyError_RegressionVsTrue);
    EleEnergyErrorVsEt_hists.push_back(EleEnergyErrorVsEt);
    EleEnergyVariance_hists.push_back(EleEnergyVariance);
    EleEnergyResponse_hists.push_back(EleEnergyResponse);
    EleEnergyOverETrueVsEOverP_hists.push_back(EleEnergyOverETrueVsEOverP);
    ElePOverETrueVsEOverP_hists.push_back(ElePOverETrueVsEOverP);
    EleRegressionEnergyNormalizedResponse_hists.push_back(EleRegressionEnergyNormalizedResponse);
    EleEnergyResponseStandardVsRegression_hists.push_back(EleEnergyResponseStandardVsRegression);
    EleEnergyResponseStandardVsRegression_Better_hists.push_back(EleEnergyResponseStandardVsRegression_Better);
    EleEnergyResponseStandardVsRegression_Worse_hists.push_back(EleEnergyResponseStandardVsRegression_Worse);
    EleEnergyMigration_Better_hists.push_back(EleEnergyMigration_Better);
    EleEnergyMigration_Worse_hists.push_back(EleEnergyMigration_Worse);

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
        if (!(RealEleTree.fGeneratedEnergyStatus3 >= RealEleTree.fGeneratedEnergyStatus1 && (RealEleTree.fGeneratedEnergyStatus3 - RealEleTree.fGeneratedEnergyStatus1)/RealEleTree.fGeneratedEnergyStatus3 < 0.01)) continue;



	//classify by eta and pt bins
	Int_t subdet = 0;
	if (fabs(RealEleTree.fEleSCEta) < 0.8) subdet = 0;
	else if (fabs(RealEleTree.fEleSCEta) < 1.479) subdet = 1;
	else subdet = 2;

	Int_t ptBin = 0;
	if (RealEleTree.fElePt > 15.0) ptBin = 1;
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
	if (Option == 11) passCuts = (ptBin == 1 );
	if (Option == 12) passCuts = (ptBin == 2 );
	if (Option == -2) passCuts = (subdet < 2);
	if (Option == -3) passCuts = (subdet == 2);
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

	// Delta E over E truth
	Double_t EoE_energy;
        Double_t ERegression;
        Double_t DeltaERegression;
	if (fabs(RealEleTree.fEleSCEta) < 1.479) {
	  EoE_energy = fabs((fEleRegressionTargetEB * RealEleTree.fSCRawEnergy - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1);
          ERegression = fEleRegressionTargetEB * RealEleTree.fSCRawEnergy;
          DeltaERegression = fEleRegressionVarianceTargetEB * RealEleTree.fSCRawEnergy;
	} else {
	  EoE_energy = fabs((fEleRegressionTargetEE * RealEleTree.fSCRawEnergy*(1+RealEleTree.fElePreShowerOverRaw) - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1);
          ERegression = fEleRegressionTargetEE * RealEleTree.fSCRawEnergy*(1+RealEleTree.fElePreShowerOverRaw);
          DeltaERegression = fEleRegressionVarianceTargetEE* RealEleTree.fSCRawEnergy*(1+RealEleTree.fElePreShowerOverRaw);
	}

	// Delta E over E from the variance regression
	Double_t EoE_variance;
	if (fabs(RealEleTree.fEleSCEta) < 1.479) {
	  EoE_variance = (fEleRegressionVarianceTargetEB / fEleRegressionTargetEB);
	} else {
	  EoE_variance = (fEleRegressionVarianceTargetEE / fEleRegressionTargetEE);
	}

	EleEnergyError_RegressionVsTrue->Fill(EoE_energy, EoE_variance, RealEleTree.fWeight);
        EleEnergyErrorVsEt->Fill(ERegression/cosh(RealEleTree.fEleEta), EoE_variance, RealEleTree.fWeight);
	EleEnergyVariance->Fill(EoE_variance, RealEleTree.fWeight);
        EleEnergyResponse->Fill( (ERegression - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1, RealEleTree.fWeight); 
        EleEnergyOverETrueVsEOverP->Fill(ERegression/RealEleTree.fGsfTrackPIn, ERegression/ RealEleTree.fGeneratedEnergyStatus1, RealEleTree.fWeight);
        ElePOverETrueVsEOverP->Fill(ERegression/RealEleTree.fGsfTrackPIn, RealEleTree.fGsfTrackPIn/RealEleTree.fGeneratedEnergyStatus1, RealEleTree.fWeight);
        EleRegressionEnergyNormalizedResponse->Fill( (ERegression - RealEleTree.fGeneratedEnergyStatus1) / DeltaERegression,  RealEleTree.fWeight);
        EleEnergyResponseStandardVsRegression->Fill( (RealEleTree.fEleSCEt*cosh(RealEleTree.fEleSCEta) - RealEleTree.fGeneratedEnergyStatus1)/RealEleTree.fGeneratedEnergyStatus1,
                                                     (ERegression - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1, RealEleTree.fWeight);
        TotalElectrons++;
        if ( 
          fabs((RealEleTree.fEleSCEt*cosh(RealEleTree.fEleSCEta) - RealEleTree.fGeneratedEnergyStatus1)/RealEleTree.fGeneratedEnergyStatus1)
          > 
          fabs((ERegression - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1 )
          ) {
          EleEnergyResponseStandardVsRegression_Better->Fill( (RealEleTree.fEleSCEt*cosh(RealEleTree.fEleSCEta) - RealEleTree.fGeneratedEnergyStatus1)/RealEleTree.fGeneratedEnergyStatus1,
                                                              (ERegression - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1, RealEleTree.fWeight);
          EleEnergyMigration_Better->Fill( fabs(RealEleTree.fEleSCEt*cosh(RealEleTree.fEleSCEta) - RealEleTree.fGeneratedEnergyStatus1)/RealEleTree.fGeneratedEnergyStatus1 - 
                                           fabs(ERegression - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1, RealEleTree.fWeight);
          RegressionImprovedElectrons++;

        } else {
          EleEnergyResponseStandardVsRegression_Worse->Fill( (RealEleTree.fEleSCEt*cosh(RealEleTree.fEleSCEta) - RealEleTree.fGeneratedEnergyStatus1)/RealEleTree.fGeneratedEnergyStatus1,
                                                     (ERegression - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1, RealEleTree.fWeight);
          EleEnergyMigration_Worse->Fill( -1.0 * (fabs(RealEleTree.fEleSCEt*cosh(RealEleTree.fEleSCEta) - RealEleTree.fGeneratedEnergyStatus1)/RealEleTree.fGeneratedEnergyStatus1 - 
                                                  fabs(ERegression - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1), RealEleTree.fWeight);
          RegressionDegradedElectrons++;
        }


      }
    }
  

//     // Now if version is V10 or V11
//     if (versionOption == "V10" || versionOption == "V11") {
//       //*****************************************************************************************
//       //Load Regression
//       //*****************************************************************************************
//       Float_t       fEleRegressionTargetEB_lowPt;
//       Float_t       fEleRegressionTargetEB_highPt;
//       Float_t       fEleRegressionTargetEE_lowPt;
//       Float_t       fEleRegressionTargetEE_highPt;
//       Float_t       fEleRegressionVarianceTargetEB_lowPt;
//       Float_t       fEleRegressionVarianceTargetEB_highPt;
//       Float_t       fEleRegressionVarianceTargetEE_lowPt;
//       Float_t       fEleRegressionVarianceTargetEE_highPt;
//       RealEleTree.InitTree();

//       // Loading the weights
//       TFile targetFile(targetFilename.c_str());
//       TTree * targetee_lowPt_tree = (TTree*) targetFile.Get("targetee_lowPt_tree");
//       TTree * targetee_highPt_tree = (TTree*) targetFile.Get("targetee_highPt_tree");
//       TTree * targeteb_lowPt_tree = (TTree*) targetFile.Get("targeteb_lowPt_tree");
//       TTree * targeteb_highPt_tree = (TTree*) targetFile.Get("targeteb_highPt_tree");
//       TTree * targeteevar_lowPt_tree = (TTree*) targetFile.Get("targeteevar_lowPt_tree");
//       TTree * targeteevar_highPt_tree = (TTree*) targetFile.Get("targeteevar_highPt_tree");
//       TTree * targetebvar_lowPt_tree = (TTree*) targetFile.Get("targetebvar_lowPt_tree");
//       TTree * targetebvar_highPt_tree = (TTree*) targetFile.Get("targetebvar_highPt_tree");

//       targeteb_lowPt_tree->SetBranchAddress( "targeteb_lowPt", &fEleRegressionTargetEB_lowPt); 
//       targeteb_highPt_tree->SetBranchAddress( "targeteb_highPt", &fEleRegressionTargetEB_highPt); 
//       targetee_lowPt_tree->SetBranchAddress( "targetee_lowPt", &fEleRegressionTargetEE_lowPt); 
//       targetee_highPt_tree->SetBranchAddress( "targetee_highPt", &fEleRegressionTargetEE_highPt); 
//       targetebvar_lowPt_tree->SetBranchAddress( "targetebvar_lowPt", &fEleRegressionVarianceTargetEB_lowPt); 
//       targetebvar_highPt_tree->SetBranchAddress( "targetebvar_highPt", &fEleRegressionVarianceTargetEB_highPt); 
//       targeteevar_lowPt_tree->SetBranchAddress( "targeteevar_lowPt", &fEleRegressionVarianceTargetEE_lowPt); 
//       targeteevar_highPt_tree->SetBranchAddress( "targeteevar_highPt", &fEleRegressionVarianceTargetEE_highPt); 

//       for(UInt_t ientry=0; ientry < RealEleTree.tree_->GetEntries(); ientry++) {       	
// 	RealEleTree.tree_->GetEntry(ientry);
// 	targeteb_lowPt_tree->GetEntry(ientry);
// 	targeteb_highPt_tree->GetEntry(ientry);
// 	targetee_lowPt_tree->GetEntry(ientry);
// 	targetee_highPt_tree->GetEntry(ientry);
// 	targetebvar_lowPt_tree->GetEntry(ientry);
// 	targetebvar_highPt_tree->GetEntry(ientry);
// 	targeteevar_lowPt_tree->GetEntry(ientry);
// 	targeteevar_highPt_tree->GetEntry(ientry);

// 	if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

// 	//don't evaluate performance using training events
// 	if (RealEleTree.fEventNumber % 2 == 0) continue;
// 	if (RealEleTree.fElePt < 7) continue;
//         if (!(RealEleTree.fGeneratedEnergyStatus3 >= RealEleTree.fGeneratedEnergyStatus1 && (RealEleTree.fGeneratedEnergyStatus3 - RealEleTree.fGeneratedEnergyStatus1)/RealEleTree.fGeneratedEnergyStatus3 < 0.01)) continue;

// 	//classify by eta and pt bins
// 	Int_t subdet = 0;
// 	if (fabs(RealEleTree.fEleSCEta) < 0.8) subdet = 0;
// 	else if (fabs(RealEleTree.fEleSCEta) < 1.479) subdet = 1;
// 	else subdet = 2;

// 	Int_t ptBin = 0;
// 	if (RealEleTree.fElePt > 10.0) ptBin = 1;
// 	if (RealEleTree.fElePt > 20.0) ptBin = 2;

// 	// Also separating by fBrem and EoP
// 	Int_t fbremBin = 0;
// 	if (RealEleTree.fEleNBrem > 0.2) fbremBin = 1;

// 	Int_t EoPBin = 0;
// 	if (RealEleTree.fEleEOverP > 0.8) EoPBin = 1;
// 	if (RealEleTree.fEleEOverP > 1.2) EoPBin = 2;

// 	Bool_t passCuts = kFALSE;
// 	if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
// 	if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
// 	if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
// 	if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
// 	if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
// 	if (Option == 5) passCuts = (subdet == 2 && ptBin == 1);    
// 	if (Option == 6) passCuts = (subdet == 0 && ptBin == 2);
// 	if (Option == 7) passCuts = (subdet == 1 && ptBin == 2);
// 	if (Option == 8) passCuts = (subdet == 2 && ptBin == 2);    
// 	if (Option == 10) passCuts = (ptBin == 0 );
// 	if (Option == 11) passCuts = (ptBin == 1 );
// 	if (Option == 12) passCuts = (ptBin == 2 );
// 	if (Option == -2) passCuts = (subdet < 2);
// 	if (Option == -3) passCuts = (subdet == 2);
// 	if (Option == -1) passCuts = kTRUE;

// 	// Now for "category" bins
// 	if (Option == 20) passCuts = (subdet < 2 && fbremBin == 0);			// Barrel, low brem
// 	if (Option == 21) passCuts = (subdet < 2 && fbremBin == 1 && EoPBin == 0);	// Barrel, high brem and low EoP	
// 	if (Option == 22) passCuts = (subdet < 2 && fbremBin == 1 && EoPBin == 1);	// Barrel, high brem and middle EoP	
// 	if (Option == 23) passCuts = (subdet < 2 && fbremBin == 1 && EoPBin == 2);	// Barrel, high brem and high EoP	

// 	if (Option == 30) passCuts = (subdet == 2 && fbremBin == 0);			// Endcap, low brem
// 	if (Option == 31) passCuts = (subdet == 2 && fbremBin == 1 && EoPBin == 0);	// Endcap, high brem and low EoP	
// 	if (Option == 32) passCuts = (subdet == 2 && fbremBin == 1 && EoPBin == 1);	// Endcap, high brem and middle EoP	
// 	if (Option == 33) passCuts = (subdet == 2 && fbremBin == 1 && EoPBin == 2);	// Endcap, high brem and high EoP	

// 	if (!passCuts) continue;    

// 	// Delta E over E from the energy regression training
// 	Double_t EoE_energy;
//         Double_t ERegression;
//         Double_t DeltaERegression;
// 	if (RealEleTree.fElePt <= 15.0) {						// For pt < 15, use pt instead of SCRawEnergy
// 	  if (fabs(RealEleTree.fEleSCEta) < 1.479) {
// 	    EoE_energy = fabs((fEleRegressionTargetEB_lowPt * RealEleTree.fElePt * TMath::CosH(RealEleTree.fEleEta) - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1);
//             ERegression = fEleRegressionTargetEB_lowPt * RealEleTree.fElePt* TMath::CosH(RealEleTree.fEleEta);
//             DeltaERegression = fEleRegressionVarianceTargetEB_lowPt * RealEleTree.fElePt* TMath::CosH(RealEleTree.fEleEta);
// 	  } else {
// 	    EoE_energy = fabs((fEleRegressionTargetEE_lowPt * RealEleTree.fElePt * TMath::CosH(RealEleTree.fEleEta) - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1);
//             ERegression = fEleRegressionTargetEE_lowPt * RealEleTree.fElePt * TMath::CosH(RealEleTree.fEleEta);
//             DeltaERegression = fEleRegressionVarianceTargetEB_highPt * RealEleTree.fElePt * TMath::CosH(RealEleTree.fEleEta);
// 	  }
// 	}
	
// 	else {
// 	  if (fabs(RealEleTree.fEleSCEta) < 1.479) {
// 	    EoE_energy = fabs((fEleRegressionTargetEB_highPt * RealEleTree.fSCRawEnergy - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1);
//             ERegression = fEleRegressionTargetEB_highPt * RealEleTree.fSCRawEnergy;
//             DeltaERegression = fEleRegressionVarianceTargetEE_lowPt* RealEleTree.fSCRawEnergy;
//           } else {
// 	    EoE_energy = fabs((fEleRegressionTargetEE_highPt * RealEleTree.fSCRawEnergy*(1+RealEleTree.fElePreShowerOverRaw) - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1);
//             ERegression = fEleRegressionTargetEE_highPt * RealEleTree.fSCRawEnergy*(1+RealEleTree.fElePreShowerOverRaw);
//             DeltaERegression = fEleRegressionVarianceTargetEE_highPt * RealEleTree.fSCRawEnergy*(1+RealEleTree.fElePreShowerOverRaw);
// 	  }
// 	}

// 	// Delta E over E from the variance regression
// 	Double_t EoE_variance;
// 	if (RealEleTree.fElePt <= 15.0) {						// For pt < 15, use pt instead of SCRawEnergy
// 	  if (fabs(RealEleTree.fEleSCEta) < 1.479) {
// 	    EoE_variance = fabs(fEleRegressionVarianceTargetEB_lowPt / fEleRegressionTargetEB_lowPt);
// 	  } else {
// 	    EoE_variance = fabs(fEleRegressionVarianceTargetEB_highPt / fEleRegressionTargetEB_highPt);
// 	  }
// 	}
	
// 	else {
// 	  if (fabs(RealEleTree.fEleSCEta) < 1.479) {
// 	    EoE_variance = fabs(fEleRegressionVarianceTargetEE_lowPt / fEleRegressionTargetEE_lowPt);
// 	  } else {
// 	    EoE_variance = fabs(fEleRegressionVarianceTargetEE_highPt / fEleRegressionTargetEE_highPt);
// 	  }
// 	}

// 	EleEnergyError_RegressionVsTrue->Fill(EoE_energy, EoE_variance, RealEleTree.fWeight);
//         EleEnergyErrorVsEt->Fill(ERegression/cosh(RealEleTree.fEleEta), EoE_variance, RealEleTree.fWeight);
// 	EleEnergyVariance->Fill(EoE_variance, RealEleTree.fWeight);
//         EleEnergyResponse->Fill( (ERegression - RealEleTree.fGeneratedEnergyStatus1) / RealEleTree.fGeneratedEnergyStatus1); 
//         EleRegressionEnergyNormalizedResponse->Fill( (ERegression - RealEleTree.fGeneratedEnergyStatus1) / DeltaERegression,  RealEleTree.fWeight);
//       }
//     }
  

  } //loop over nfiles


  TCanvas *cv = 0;
  TLegend *legend = 0;


  //******************************************************************************************************
  //Save Histograms
  //******************************************************************************************************
  TFile *file = new TFile("RegressionVariancePlots.root", "UPDATE");
  file->cd();
  file->WriteTObject(EleEnergyResponse_hists[0], ("EleEnergyResponse_V00_Bin"+OptionStr).c_str(), "WriteDelete");
//   file->WriteTObject(EleEnergyResponse_hists[1], ("EleEnergyResponse_V01_Bin"+OptionStr).c_str(), "WriteDelete");
  file->WriteTObject(EleEnergyVariance_hists[0], ("EleEnergyVariance_V00_Bin"+OptionStr).c_str(), "WriteDelete");
//   file->WriteTObject(EleEnergyVariance_hists[1], ("EleEnergyVariance_V01_Bin"+OptionStr).c_str(), "WriteDelete");
  file->WriteTObject(EleRegressionEnergyNormalizedResponse_hists[0], ("EleRegressionEnergyNormalizedResponse_V00_Bin"+OptionStr).c_str(), "WriteDelete");
//   file->WriteTObject(EleRegressionEnergyNormalizedResponse_hists[1], ("EleRegressionEnergyNormalizedResponse_V01_Bin"+OptionStr).c_str(), "WriteDelete");

  file->WriteTObject(EleEnergyMigration_Better_hists[0], ("EleEnergyMigration_Better_V00_Bin"+OptionStr).c_str(), "WriteDelete");
  file->WriteTObject(EleEnergyMigration_Worse_hists[0], ("EleEnergyMigration_Worse_V00_Bin"+OptionStr).c_str(), "WriteDelete");


  //******************************************************************************************************
  //Plot variance comparison
  //******************************************************************************************************
  float max_yrange = 0;
  for (int j = 0; j < nFiles; j++) {
    float max = EleEnergyVariance_hists[j]->GetMaximum();
    if (max > max_yrange) max_yrange = max;
  }
  // Settings for histograms
  EleEnergyVariance_hists[0]->GetYaxis()->SetTitleOffset(1.5);
  EleEnergyVariance_hists[0]->SetMaximum(max_yrange * 1.25);

  EleEnergyVariance_hists[0]->SetLineColor(kRed);
//   EleEnergyVariance_hists[1]->SetLineColor(kGreen);
  EleEnergyVariance_hists[0]->SetLineWidth(2);
//   EleEnergyVariance_hists[1]->SetLineWidth(2);

  cv = new TCanvas("cv1", "cv1", 800, 600);
  legend = new TLegend(0.4,0.7,0.8,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  EleEnergyVariance_hists[0]->GetYaxis()->SetTitle("Number of Electrons");
  EleEnergyVariance_hists[0]->GetXaxis()->SetRangeUser(0,0.15);
  EleEnergyVariance_hists[0]->Draw("hist");
//   EleEnergyVariance_hists[1]->Draw("histsame");

  for (int j = 0; j < nFiles; j++) {
    // Prepare string for version name
    string versionName;
    if (versionOptions[j] == "V00") versionName = "no pt split, no track vars"; 
    if (versionOptions[j] == "V01") versionName = "no pt split, track vars"; 
    if (versionOptions[j] == "V10") versionName = "pt split, no track vars"; 
    if (versionOptions[j] == "V11") versionName = "pt split, track vars"; 
    if (versionOptions[j] == "V5")  versionName = "no pt split, track vars and track error";
    legend->AddEntry(EleEnergyVariance_hists[j], ("Regression, " + versionName).c_str(), "L");
  }
  legend->Draw();
 
  cv->SaveAs((outDirectoryName + "EnergyVariance_binOption" + OptionStr + ".gif").c_str());
  cv->SaveAs((outDirectoryName + "EnergyVariance_binOption" + OptionStr + ".pdf").c_str());
  cv->SaveAs((outDirectoryName + "EnergyVariance_binOption" + OptionStr + ".root").c_str());




  //******************************************************************************************************
  // Plotting histogram
  //******************************************************************************************************
  cv = new TCanvas("cv4", "cv4", 800, 600);
  legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);


  TProfile *EleEnergyError_RegressionVsTrue_Profile0 = EleEnergyError_RegressionVsTrue_hists[0]->ProfileX();
  EleEnergyError_RegressionVsTrue_Profile0->GetYaxis()->SetTitle("#Delta E / E (regression)");
//   EleEnergyError_RegressionVsTrue_Profile0->GetXaxis()->SetRangeUser(0,0.05);
  EleEnergyError_RegressionVsTrue_Profile0->GetYaxis()->SetTitleOffset(1.5);
  EleEnergyError_RegressionVsTrue_Profile0->Draw();
  cv->SaveAs((outDirectoryName + "EnergyVarianceProfile_" + versionOptions[0] + "_binOption" + OptionStr + ".gif").c_str());
  cv->SaveAs((outDirectoryName + "EnergyVarianceProfile_" + versionOptions[0] + "_binOption" + OptionStr + ".pdf").c_str());
  cv->SaveAs((outDirectoryName + "EnergyVarianceProfile_" + versionOptions[0] + "_binOption" + OptionStr + ".root").c_str());

//   cv = new TCanvas("cv5", "cv5", 800, 600);
//   TProfile *EleEnergyError_RegressionVsTrue_Profile1 = EleEnergyError_RegressionVsTrue_hists[1]->ProfileX();
//   EleEnergyError_RegressionVsTrue_Profile1->GetYaxis()->SetTitle("#Delta E / E (regression)");
// //   EleEnergyError_RegressionVsTrue_Profile1->GetXaxis()->SetRangeUser(0,0.05);
//   EleEnergyError_RegressionVsTrue_Profile1->GetYaxis()->SetTitleOffset(1.5);
//   EleEnergyError_RegressionVsTrue_Profile1->Draw();
//   cv->SaveAs((outDirectoryName + "EnergyVarianceProfile_" + versionOptions[1] + "_binOption" + OptionStr + ".gif").c_str());
//   cv->SaveAs((outDirectoryName + "EnergyVarianceProfile_" + versionOptions[1] + "_binOption" + OptionStr + ".pdf").c_str());
//   cv->SaveAs((outDirectoryName + "EnergyVarianceProfile_" + versionOptions[1] + "_binOption" + OptionStr + ".root").c_str());


  //******************************************************************************************************
  // Plot Regression Error Vs True Error
  //******************************************************************************************************
  cv = new TCanvas("cv2", "cv2", 800, 600);
  legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  EleEnergyError_RegressionVsTrue_hists[0]->GetYaxis()->SetTitleOffset(1.5);
  EleEnergyError_RegressionVsTrue_hists[0]->GetXaxis()->SetRangeUser(0,0.05);
  EleEnergyError_RegressionVsTrue_hists[0]->GetYaxis()->SetRangeUser(0,0.05);
  EleEnergyError_RegressionVsTrue_hists[0]->Draw("colz");
  cv->SaveAs((outDirectoryName + "EnergyVariance2D_" + versionOptions[0] + "_binOption" + OptionStr + ".gif").c_str());
  cv->SaveAs((outDirectoryName + "EnergyVariance2D_" + versionOptions[0] + "_binOption" + OptionStr + ".pdf").c_str());
  cv->SaveAs((outDirectoryName + "EnergyVariance2D_" + versionOptions[0] + "_binOption" + OptionStr + ".root").c_str());

//   cv = new TCanvas("cv3", "cv3", 800, 600);
//   EleEnergyError_RegressionVsTrue_hists[1]->GetYaxis()->SetTitleOffset(1.5);
//   EleEnergyError_RegressionVsTrue_hists[1]->GetXaxis()->SetRangeUser(0,0.05);
//   EleEnergyError_RegressionVsTrue_hists[1]->GetYaxis()->SetRangeUser(0,0.05);
//   EleEnergyError_RegressionVsTrue_hists[1]->Draw("colz");
//   cv->SaveAs((outDirectoryName + "EnergyVariance2D_" + versionOptions[1] + "_binOption" + OptionStr + ".gif").c_str());
//   cv->SaveAs((outDirectoryName + "EnergyVariance2D_" + versionOptions[1] + "_binOption" + OptionStr + ".pdf").c_str());
//   cv->SaveAs((outDirectoryName + "EnergyVariance2D_" + versionOptions[1] + "_binOption" + OptionStr + ".root").c_str());



  //******************************************************************************************************
  // Plotting histogram
  //******************************************************************************************************
  cv = new TCanvas("cv4", "cv4", 800, 600);
  legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);


  TProfile *EleEnergyErrorVsEt_Profile0 = EleEnergyErrorVsEt_hists[0]->ProfileX();
  EleEnergyErrorVsEt_Profile0->GetYaxis()->SetTitle("#Delta E / E (regression)");
  EleEnergyErrorVsEt_Profile0->GetXaxis()->SetRangeUser(7,50);
  EleEnergyErrorVsEt_Profile0->GetYaxis()->SetTitleOffset(1.5);
  EleEnergyErrorVsEt_Profile0->Draw();
  cv->SaveAs((outDirectoryName + "EnergyErrorVsEtProfile_" + versionOptions[0] + "_binOption" + OptionStr + ".gif").c_str());
  cv->SaveAs((outDirectoryName + "EnergyErrorVsEtProfile_" + versionOptions[0] + "_binOption" + OptionStr + ".pdf").c_str());
  cv->SaveAs((outDirectoryName + "EnergyErrorVsEtProfile_" + versionOptions[0] + "_binOption" + OptionStr + ".root").c_str());

//   cv = new TCanvas("cv5", "cv5", 800, 600);
//   TProfile *EleEnergyErrorVsEt_Profile1 = EleEnergyErrorVsEt_hists[1]->ProfileX();
//   EleEnergyErrorVsEt_Profile1->GetYaxis()->SetTitle("#Delta E / E (regression)");
//   EleEnergyErrorVsEt_Profile1->GetXaxis()->SetRangeUser(7,50);
//   EleEnergyErrorVsEt_Profile1->GetYaxis()->SetTitleOffset(1.5);
//   EleEnergyErrorVsEt_Profile1->Draw();
//   cv->SaveAs((outDirectoryName + "EnergyErrorVsEtProfile_" + versionOptions[1] + "_binOption" + OptionStr + ".gif").c_str());
//   cv->SaveAs((outDirectoryName + "EnergyErrorVsEtProfile_" + versionOptions[1] + "_binOption" + OptionStr + ".pdf").c_str());
//   cv->SaveAs((outDirectoryName + "EnergyErrorVsEtProfile_" + versionOptions[1] + "_binOption" + OptionStr + ".root").c_str());


  //******************************************************************************************************
  // Plot Regression Error Vs True Error
  //******************************************************************************************************
  cv = new TCanvas("cv2", "cv2", 800, 600);
  legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  EleEnergyOverETrueVsEOverP_hists[0]->GetYaxis()->SetTitleOffset(1.5);
  EleEnergyOverETrueVsEOverP_hists[0]->Draw("colz");
  cv->SaveAs((outDirectoryName + "EleEnergyOverETrueVsEOverP_" + versionOptions[0] + "_binOption" + OptionStr + ".gif").c_str());
  cv->SaveAs((outDirectoryName + "EleEnergyOverETrueVsEOverP_" + versionOptions[0] + "_binOption" + OptionStr + ".pdf").c_str());
  cv->SaveAs((outDirectoryName + "EleEnergyOverETrueVsEOverP_" + versionOptions[0] + "_binOption" + OptionStr + ".root").c_str());


  cv = new TCanvas("cv2", "cv2", 800, 600);
  legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  ElePOverETrueVsEOverP_hists[0]->GetYaxis()->SetTitleOffset(1.5);
  ElePOverETrueVsEOverP_hists[0]->Draw("colz");
  cv->SaveAs((outDirectoryName + "ElePOverETrueVsEOverP_" + versionOptions[0] + "_binOption" + OptionStr + ".gif").c_str());
  cv->SaveAs((outDirectoryName + "ElePOverETrueVsEOverP_" + versionOptions[0] + "_binOption" + OptionStr + ".pdf").c_str());
  cv->SaveAs((outDirectoryName + "ElePOverETrueVsEOverP_" + versionOptions[0] + "_binOption" + OptionStr + ".root").c_str());



  //******************************************************************************************************
  // Plot Regression Response Vs Standard Response
  //******************************************************************************************************
  cv = new TCanvas("cv2", "cv2", 800, 600);
  legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  EleEnergyResponseStandardVsRegression_hists[0]->GetYaxis()->SetTitleOffset(1.5);
//   EleEnergyResponseStandardVsRegression_hists[0]->GetXaxis()->SetRangeUser(0,0.05);
//   EleEnergyResponseStandardVsRegression_hists[0]->GetYaxis()->SetRangeUser(0,0.05);
  EleEnergyResponseStandardVsRegression_hists[0]->Draw("colz");
  cv->SaveAs((outDirectoryName + "EleEnergyResponseStandardVsRegression_" + versionOptions[0] + "_binOption" + OptionStr + ".gif").c_str());
  cv->SaveAs((outDirectoryName + "EleEnergyResponseStandardVsRegression_" + versionOptions[0] + "_binOption" + OptionStr + ".pdf").c_str());
  cv->SaveAs((outDirectoryName + "EleEnergyResponseStandardVsRegression_" + versionOptions[0] + "_binOption" + OptionStr + ".root").c_str());

  cv = new TCanvas("cv2", "cv2", 800, 600);
  legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  EleEnergyResponseStandardVsRegression_Better_hists[0]->GetYaxis()->SetTitleOffset(1.5);
//   EleEnergyResponseStandardVsRegression_Better_hists[0]->GetXaxis()->SetRangeUser(0,0.05);
//   EleEnergyResponseStandardVsRegression_Better_hists[0]->GetYaxis()->SetRangeUser(0,0.05);
  EleEnergyResponseStandardVsRegression_Better_hists[0]->Draw("colz");
  cv->SaveAs((outDirectoryName + "EleEnergyResponseStandardVsRegression_Better_" + versionOptions[0] + "_binOption" + OptionStr + ".gif").c_str());
  cv->SaveAs((outDirectoryName + "EleEnergyResponseStandardVsRegression_Better_" + versionOptions[0] + "_binOption" + OptionStr + ".pdf").c_str());
  cv->SaveAs((outDirectoryName + "EleEnergyResponseStandardVsRegression_Better_" + versionOptions[0] + "_binOption" + OptionStr + ".root").c_str());

  cv = new TCanvas("cv2", "cv2", 800, 600);
  legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  EleEnergyResponseStandardVsRegression_Worse_hists[0]->GetYaxis()->SetTitleOffset(1.5);
//   EleEnergyResponseStandardVsRegression_Worse_hists[0]->GetXaxis()->SetRangeUser(0,0.05);
//   EleEnergyResponseStandardVsRegression_Worse_hists[0]->GetYaxis()->SetRangeUser(0,0.05);
  EleEnergyResponseStandardVsRegression_Worse_hists[0]->Draw("colz");
  cv->SaveAs((outDirectoryName + "EleEnergyResponseStandardVsRegression_Worse_" + versionOptions[0] + "_binOption" + OptionStr + ".gif").c_str());
  cv->SaveAs((outDirectoryName + "EleEnergyResponseStandardVsRegression_Worse_" + versionOptions[0] + "_binOption" + OptionStr + ".pdf").c_str());
  cv->SaveAs((outDirectoryName + "EleEnergyResponseStandardVsRegression_Worse_" + versionOptions[0] + "_binOption" + OptionStr + ".root").c_str());

  cout << "Total Electrons: " << TotalElectrons << endl;
  cout << "Regression Improved Electrons: " << RegressionImprovedElectrons << endl;
  cout << "Regression Degraded Electrons: " << RegressionDegradedElectrons << endl;


  //******************************************************************************************************
  // Plot Energy Migration
  //******************************************************************************************************
  cv = new TCanvas("cv2", "cv2", 800, 600);
  legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  legend->AddEntry(EleEnergyMigration_Better_hists[0], "Better", "L");
  legend->AddEntry(EleEnergyMigration_Worse_hists[0], "Worse", "L");
  EleEnergyMigration_Better_hists[0]->SetLineColor(kBlue);
  EleEnergyMigration_Worse_hists[0]->SetLineColor(kRed);

  EleEnergyMigration_Better_hists[0]->GetYaxis()->SetTitleOffset(1.5);
//   EleEnergyMigration_Better_hists[0]->GetXaxis()->SetRangeUser(0,0.05);
//   EleEnergyMigration_Better_hists[0]->GetYaxis()->SetRangeUser(0,0.05);
  EleEnergyMigration_Better_hists[0]->Draw("hist");
  EleEnergyMigration_Worse_hists[0]->Draw("hist,same");
  legend->Draw();

  cv->SaveAs((outDirectoryName + "EleEnergyMigration_" + versionOptions[0] + "_binOption" + OptionStr + ".gif").c_str());
//   cv->SaveAs((outDirectoryName + "EleEnergyMigration_Better_" + versionOptions[0] + "_binOption" + OptionStr + ".pdf").c_str());
//   cv->SaveAs((outDirectoryName + "EleEnergyMigration_Better_" + versionOptions[0] + "_binOption" + OptionStr + ".root").c_str());




}


void plotRegressionVarianceComparison(string applyingFilename, string targetFilename1, string targetFilename2, string versionOption1, string versionOption2, string outDirectoryName, Int_t Option) {
 
  makeRegressionVarianceComparison(applyingFilename, targetFilename1, targetFilename2, versionOption1, versionOption2, outDirectoryName, Option);
}

void plotRegressionVarianceComparison() {
  ValidateRegressionVariance("RegressionVariancePlots.root");
}
