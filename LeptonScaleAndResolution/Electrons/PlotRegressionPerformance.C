//================================================================================================
//		PlotRegressionPerformance.C		(Example)
//	Plots pT over gen for an electron ntuple
//	USAGE:
//		plotRegressionPerformanceComparison(string applyingFilename, string outDirectoryName, Int_t Option)
//	
//	applyingFilename = ntuple file from which to read values
//	outDirectoryName = directory to which outputs the histograms (has to end in /)
//	Option = eta and pt bin to be plotted
//	0 = low eta, low pt
//	1 = medium eta, low pt
//	2 = endcap, low pt
//	3 = low eta, medium pt
//	4 = medium eta, medium pt
//	5 = endcap, medium pt
//	6 = low eta, high pt
//	7 = medium eta, high pt
//	8 = endcap, high pt
//	10= low pt
//	-1= all
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
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void PlotRegressionPerformance(string applyingFilename, string label, Int_t Option) {  
  // Plots, for comparison, the results of the regression of applyElectronEnergyRegression_option.C
  //
  // applyingFilename	= file to which the regression is going to be applied
  // Option		= specifies which eta and phi bins will be plotted

  string Label = label;
  if (label != "") Label = "_" + label;

  // Setting up graphic style
  gStyle->SetOptStat(0);

  // Preparing Option string to include in name
  string OptionStr; 
  ostringstream OptionConv;
  OptionConv << Option;
  OptionStr = OptionConv.str();

  // Preparing strings for name
  string binNameStr;
  if (Option == 0)  binNameStr = "Electron energy ratio, Eta < 0.8, 7 < Pt < 10";
  if (Option == 1)  binNameStr = "Electron energy ratio, 0.8 < Eta < 1.485, 7 < Pt < 10";
  if (Option == 2)  binNameStr = "Electron energy ratio, Eta > 1.485, 7 < Pt < 10";
  if (Option == 3)  binNameStr = "Electron energy ratio, Eta < 0.8, 10 < Pt < 20";
  if (Option == 4)  binNameStr = "Electron energy ratio, 0.8 < Eta < 1.485, 10 < Pt < 20";
  if (Option == 5)  binNameStr = "Electron energy ratio, Eta > 1.485, 10 < Pt < 20";
  if (Option == 6)  binNameStr = "Electron energy ratio, Eta < 0.8, Pt > 20";
  if (Option == 7)  binNameStr = "Electron energy ratio, 0.8 < Eta < 1.485, Pt > 20";
  if (Option == 8)  binNameStr = "Electron energy ratio, Eta > 1.485, Pt > 20";
  if (Option == 10) binNameStr = "Electron energy ratio, 7 < Pt < 10";
  if (Option == -1) binNameStr = "Electron energy ratio, Pt > 7";
  /*
     if (Option == 20) binNameStr = "Electron energy ratio, barrel, brem < 0.2";
     if (Option == 21) binNameStr = "Electron energy ratio, barrel, brem > 0.2, EoP < 0.8";
     if (Option == 22) binNameStr = "Electron energy ratio, barrel, brem > 0.2, 0.8 < EoP < 1.2";
     if (Option == 23) binNameStr = "Electron energy ratio, barrel, brem > 0.2, EoP > 1.2";
     if (Option == 24) binNameStr = "Electron energy ratio, barrel, brem > 0.2, EoP < 0.8 or EoP > 1.2";

     if (Option == 30) binNameStr = "Electron energy ratio, endcap, brem < 0.2";
     if (Option == 31) binNameStr = "Electron energy ratio, endcap, brem > 0.2, EoP < 0.8";
     if (Option == 32) binNameStr = "Electron energy ratio, endcap, brem > 0.2, 0.8 < EoP < 1.2";
     if (Option == 33) binNameStr = "Electron energy ratio, endcap, brem > 0.2, EoP > 1.2";
     if (Option == 34) binNameStr = "Electron energy ratio, endcap, brem > 0.2, EoP < 0.8 or EoP > 1.2";
   */
  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *EleEnergyOverGenEnergy = new TH1F("EleEnergyOverGenEnergy", (binNameStr + "; Electron Energy / Generated Energy ; Number of Events ").c_str(),  200, 0, 2);

  //*****************************************************************************************
  //RealEleTree
  //*****************************************************************************************
  citana::ElectronTree RealEleTree;
  RealEleTree.LoadTree(applyingFilename.c_str());
  RealEleTree.InitTree();

  for(UInt_t ientry=0; ientry < RealEleTree.tree_->GetEntries(); ientry++) {       	
    RealEleTree.tree_->GetEntry(ientry);

    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    if (RealEleTree.fElePt < 7) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(RealEleTree.fEleEta) < 0.8) subdet = 0;
    else if (fabs(RealEleTree.fEleEta) < 1.485) subdet = 1;
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

    /*
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
     */
    if (!passCuts) continue;    

    EleEnergyOverGenEnergy->Fill(RealEleTree.fElePt * TMath::CosH(RealEleTree.fEleEta) / RealEleTree.fGeneratedEnergyStatus1 ,RealEleTree.fWeight);
  } 

  TCanvas *cv = 0;
  TLegend *legend = 0;

  // Settings for histograms
  EleEnergyOverGenEnergy->GetYaxis()->SetTitleOffset(1.5);
  EleEnergyOverGenEnergy->SetLineColor(kBlack);

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.1,0.7,0.3,0.9);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(EleEnergyOverGenEnergy,"Standard", "L");
  EleEnergyOverGenEnergy->Draw("hist");
  legend->Draw();

  cv->SaveAs(("EnergyResponse" + Label + "binOption" + OptionStr + ".gif").c_str());

} 

