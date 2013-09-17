//================================================================================================
//					ValidateLL
// 
// Reads log-likelihood branches from fullsim and fastsim ntuples and plots them for comparison
//
// USAGE
//
// ValidateLL(string fullsimFilename, string fastsimFilename, string outputDirectory, Int_t norm_option = 1,  Int_t error_option = 2, const string Label = "ZZ")
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

#endif
//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
    if (b == 0 || (b == hist->GetXaxis()->GetNbins()+1)) 
      cout << "Bin " << b << "'s content: " << hist->GetBinContent(b) << endl;
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
void ValidateLL(string fullsimFilename, string fastsimFilename, string outputDirectory, Int_t norm_option = 1, Int_t error_option = 2,  const string Label = "ZZ") {

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
  TFile outRootFile((outputDirectory + "ValidateLL_AllHistograms" + label + ".root").c_str(), "recreate");

  TH1F* histLLscalar_fullsim_4e = new TH1F( "histLLscalar_fullsim_4e", "Log-likelihood, scalar, 4e;LLscalar;Number of Events", 40, -40, 0);
  TH1F* histLLscalar_fastsim_4e = new TH1F( "histLLscalar_fastsim_4e", "Log-likelihood, scalar, 4e;LLscalar;Number of Events", 40, -40, 0);
  TH1F* histLLscalar_fullsim_4m = new TH1F( "histLLscalar_fullsim_4m", "Log-likelihood, scalar, 4m;LLscalar;Number of Events", 40, -40, 0);
  TH1F* histLLscalar_fastsim_4m = new TH1F( "histLLscalar_fastsim_4m", "Log-likelihood, scalar, 4m;LLscalar;Number of Events", 40, -40, 0);
  TH1F* histLLscalar_fullsim_2e2m = new TH1F( "histLLscalar_fullsim_2e2m", "Log-likelihood, scalar, 2e2m;LLscalar;Number of Events", 40, -40, 0);
  TH1F* histLLscalar_fastsim_2e2m = new TH1F( "histLLscalar_fastsim_2e2m", "Log-likelihood, scalar, 2e2m;LLscalar;Number of Events", 40, -40, 0);

  TH1F* histLLpseudoscalar_fullsim_4e = new TH1F( "histLLpseudoscalar_fullsim_4e", "Log-likelihood, pseudoscalar, 4e;LLpseudoscalar;Number of Events", 40, -40, 0);
  TH1F* histLLpseudoscalar_fastsim_4e = new TH1F( "histLLpseudoscalar_fastsim_4e", "Log-likelihood, pseudoscalar, 4e;LLpseudoscalar;Number of Events", 40, -40, 0);
  TH1F* histLLpseudoscalar_fullsim_4m = new TH1F( "histLLpseudoscalar_fullsim_4m", "Log-likelihood, pseudoscalar, 4m;LLpseudoscalar;Number of Events", 40, -40, 0);
  TH1F* histLLpseudoscalar_fastsim_4m = new TH1F( "histLLpseudoscalar_fastsim_4m", "Log-likelihood, pseudoscalar, 4m;LLpseudoscalar;Number of Events", 40, -40, 0);
  TH1F* histLLpseudoscalar_fullsim_2e2m = new TH1F( "histLLpseudoscalar_fullsim_2e2m", "Log-likelihood, pseudoscalar, 2e2m;LLpseudoscalar;Number of Events", 40, -40, 0);
  TH1F* histLLpseudoscalar_fastsim_2e2m = new TH1F( "histLLpseudoscalar_fastsim_2e2m", "Log-likelihood, pseudoscalar, 2e2m;LLpseudoscalar;Number of Events", 40, -40, 0);

  TH1F* histLLratio_fullsim_4e = new TH1F( "histLLratio_fullsim_4e", "Log-likelihood-ratio (scalar/pseudoscalar), 4e;LLratio;Number of Events", 40, -20, 20);
  TH1F* histLLratio_fastsim_4e = new TH1F( "histLLratio_fastsim_4e", "Log-likelihood-ratio (scalar/pseudoscalar), 4e;LLratio;Number of Events", 40, -20, 20);
  TH1F* histLLratio_fullsim_4m = new TH1F( "histLLratio_fullsim_4m", "Log-likelihood-ratio (scalar/pseudoscalar), 4m;LLratio;Number of Events", 40, -20, 20);
  TH1F* histLLratio_fastsim_4m = new TH1F( "histLLratio_fastsim_4m", "Log-likelihood-ratio (scalar/pseudoscalar), 4m;LLratio;Number of Events", 40, -20, 20);
  TH1F* histLLratio_fullsim_2e2m = new TH1F( "histLLratio_fullsim_2e2m", "Log-likelihood-ratio (scalar/pseudoscalar), 2e2m;LLratio;Number of Events", 40, -20, 20);
  TH1F* histLLratio_fastsim_2e2m = new TH1F( "histLLratio_fastsim_2e2m", "Log-likelihood-ratio (scalar/pseudoscalar), 2e2m;LLratio;Number of Events", 40, -20, 20);

  //*****************************************************************************************
  // Loop over fullsim events
  //*****************************************************************************************
  for(UInt_t ientry=0; ientry < fullsimHZZEventTree.tree_->GetEntries(); ientry++) {       	
    fullsimHZZEventTree.tree_->GetEntry(ientry);
    if (ientry % 1000000 == 0) cout << "Event " << ientry << endl;

    // Setting up normalization options
    Float_t weight = fullsimHZZEventTree.fWeight; 

    // Filling Z masses
    if (fullsimHZZEventTree.fPassFullSelection) {

      // Filling log-likelihoods
      if (abs(fullsimHZZEventTree.fLep1Type) == 11 && abs(fullsimHZZEventTree.fLep3Type) == 11) {
	histLLscalar_fullsim_4e->Fill(fullsimHZZEventTree.fLogLikelihood_scalar, weight);
	histLLpseudoscalar_fullsim_4e->Fill(fullsimHZZEventTree.fLogLikelihood_pseudoscalar, weight);
	histLLratio_fullsim_4e->Fill(fullsimHZZEventTree.fLogLikelihood_scalar - fullsimHZZEventTree.fLogLikelihood_pseudoscalar, weight);
      } else if (abs(fullsimHZZEventTree.fLep1Type) == 13 && abs(fullsimHZZEventTree.fLep3Type) == 13) {
	histLLscalar_fullsim_4m->Fill(fullsimHZZEventTree.fLogLikelihood_scalar, weight);
	histLLpseudoscalar_fullsim_4m->Fill(fullsimHZZEventTree.fLogLikelihood_pseudoscalar, weight);
	histLLratio_fullsim_4m->Fill(fullsimHZZEventTree.fLogLikelihood_scalar - fullsimHZZEventTree.fLogLikelihood_pseudoscalar, weight);
      } else {
	histLLscalar_fullsim_2e2m->Fill(fullsimHZZEventTree.fLogLikelihood_scalar, weight);
	histLLpseudoscalar_fullsim_2e2m->Fill(fullsimHZZEventTree.fLogLikelihood_pseudoscalar, weight);
	histLLratio_fullsim_2e2m->Fill(fullsimHZZEventTree.fLogLikelihood_scalar - fullsimHZZEventTree.fLogLikelihood_pseudoscalar, weight);
      }
    }
  }

  //*****************************************************************************************
  // Loop over fastsim events
  //*****************************************************************************************
  for(UInt_t ientry=0; ientry < fastsimHZZEventTree.tree_->GetEntries(); ientry++) {       	
    fastsimHZZEventTree.tree_->GetEntry(ientry);
    if (ientry % 1000000 == 0) cout << "Event " << ientry << endl;

    // Setting up normalization options
    Float_t weight = fastsimHZZEventTree.fWeight; 

    // Filling Z masses
    if (fastsimHZZEventTree.fPassFullSelection) {

      // Filling log-likelihoods
      if (abs(fastsimHZZEventTree.fLep1Type) == 11 && abs(fastsimHZZEventTree.fLep3Type) == 11) {
	histLLscalar_fastsim_4e->Fill(fastsimHZZEventTree.fLogLikelihood_scalar, weight);
	histLLpseudoscalar_fastsim_4e->Fill(fastsimHZZEventTree.fLogLikelihood_pseudoscalar, weight);
	histLLratio_fastsim_4e->Fill(fastsimHZZEventTree.fLogLikelihood_scalar - fastsimHZZEventTree.fLogLikelihood_pseudoscalar, weight);
      } else if (abs(fastsimHZZEventTree.fLep1Type) == 13 && abs(fastsimHZZEventTree.fLep3Type) == 13) {
	histLLscalar_fastsim_4m->Fill(fastsimHZZEventTree.fLogLikelihood_scalar, weight);
	histLLpseudoscalar_fastsim_4m->Fill(fastsimHZZEventTree.fLogLikelihood_pseudoscalar, weight);
	histLLratio_fastsim_4m->Fill(fastsimHZZEventTree.fLogLikelihood_scalar - fastsimHZZEventTree.fLogLikelihood_pseudoscalar, weight);
      } else {
	histLLscalar_fastsim_2e2m->Fill(fastsimHZZEventTree.fLogLikelihood_scalar, weight);
	histLLpseudoscalar_fastsim_2e2m->Fill(fastsimHZZEventTree.fLogLikelihood_pseudoscalar, weight);
	histLLratio_fastsim_2e2m->Fill(fastsimHZZEventTree.fLogLikelihood_scalar - fastsimHZZEventTree.fLogLikelihood_pseudoscalar, weight);
      }
    }
  }

  //*****************************************************************************************
  // Now, normalizing if option one was selected
  if (norm_option) {
    NormalizeHist(histLLscalar_fullsim_4e);
    NormalizeHist(histLLscalar_fastsim_4e);
    NormalizeHist(histLLscalar_fullsim_4m);
    NormalizeHist(histLLscalar_fastsim_4m);
    NormalizeHist(histLLscalar_fullsim_2e2m);
    NormalizeHist(histLLscalar_fastsim_2e2m);
    
    NormalizeHist(histLLpseudoscalar_fullsim_4e);
    NormalizeHist(histLLpseudoscalar_fastsim_4e);
    NormalizeHist(histLLpseudoscalar_fullsim_4m);
    NormalizeHist(histLLpseudoscalar_fastsim_4m);
    NormalizeHist(histLLpseudoscalar_fullsim_2e2m);
    NormalizeHist(histLLpseudoscalar_fastsim_2e2m);

    NormalizeHist(histLLratio_fullsim_4e);
    NormalizeHist(histLLratio_fastsim_4e);
    NormalizeHist(histLLratio_fullsim_4m);
    NormalizeHist(histLLratio_fastsim_4m);
    NormalizeHist(histLLratio_fullsim_2e2m);
    NormalizeHist(histLLratio_fastsim_2e2m);
  }

  //*****************************************************************************************
  // Plot
  //*****************************************************************************************

  // Setting label offset
  histLLscalar_fullsim_4e->GetYaxis()->SetTitleOffset(1.7);
  histLLscalar_fastsim_4e->GetYaxis()->SetTitleOffset(1.7);
  histLLscalar_fullsim_4m->GetYaxis()->SetTitleOffset(1.7);
  histLLscalar_fastsim_4m->GetYaxis()->SetTitleOffset(1.7);
  histLLscalar_fullsim_2e2m->GetYaxis()->SetTitleOffset(1.7);
  histLLscalar_fastsim_2e2m->GetYaxis()->SetTitleOffset(1.7);

  histLLpseudoscalar_fullsim_4e->GetYaxis()->SetTitleOffset(1.7);
  histLLpseudoscalar_fastsim_4e->GetYaxis()->SetTitleOffset(1.7);
  histLLpseudoscalar_fullsim_4m->GetYaxis()->SetTitleOffset(1.7);
  histLLpseudoscalar_fastsim_4m->GetYaxis()->SetTitleOffset(1.7);
  histLLpseudoscalar_fullsim_2e2m->GetYaxis()->SetTitleOffset(1.7);
  histLLpseudoscalar_fastsim_2e2m->GetYaxis()->SetTitleOffset(1.7);

  histLLratio_fullsim_4e->GetYaxis()->SetTitleOffset(1.7);
  histLLratio_fastsim_4e->GetYaxis()->SetTitleOffset(1.7);
  histLLratio_fullsim_4m->GetYaxis()->SetTitleOffset(1.7);
  histLLratio_fastsim_4m->GetYaxis()->SetTitleOffset(1.7);
  histLLratio_fullsim_2e2m->GetYaxis()->SetTitleOffset(1.7);
  histLLratio_fastsim_2e2m->GetYaxis()->SetTitleOffset(1.7);

  // Saving all histograms to file
  outRootFile.Write();

  // Plotting all the histograms
  PlotTwo(histLLscalar_fullsim_4e, histLLscalar_fastsim_4e, (outputDirectory + "ValidateLL_LLscalar_4e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLLscalar_fullsim_4m, histLLscalar_fastsim_4m, (outputDirectory + "ValidateLL_LLscalar_4m" + label + ".gif").c_str(), error_option);
  PlotTwo(histLLscalar_fullsim_2e2m, histLLscalar_fastsim_2e2m, (outputDirectory + "ValidateLL_LLscalar_2e2m" + label + ".gif").c_str(), error_option);

  PlotTwo(histLLpseudoscalar_fullsim_4e, histLLpseudoscalar_fastsim_4e, (outputDirectory + "ValidateLL_LLpseudoscalar_4e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLLpseudoscalar_fullsim_4m, histLLpseudoscalar_fastsim_4m, (outputDirectory + "ValidateLL_LLpseudoscalar_4m" + label + ".gif").c_str(), error_option);
  PlotTwo(histLLpseudoscalar_fullsim_2e2m, histLLpseudoscalar_fastsim_2e2m, (outputDirectory + "ValidateLL_LLpseudoscalar_2e2m" + label + ".gif").c_str(), error_option);

  PlotTwo(histLLratio_fullsim_4e, histLLratio_fastsim_4e, (outputDirectory + "ValidateLL_LLratio_4e" + label + ".gif").c_str(), error_option);
  PlotTwo(histLLratio_fullsim_4m, histLLratio_fastsim_4m, (outputDirectory + "ValidateLL_LLratio_4m" + label + ".gif").c_str(), error_option);
  PlotTwo(histLLratio_fullsim_2e2m, histLLratio_fastsim_2e2m, (outputDirectory + "ValidateLL_LLratio_2e2m" + label + ".gif").c_str(), error_option);
}
