

//======================================================== //
// For ECAL with CMS Detector at LHC                       //
// Roofit Macro for Unbinned fit to Z peak                 //
//======================================================== //

#ifndef __CINT__
#include<stdio.h>
#include<string>
#include<sstream> 
#include<iostream>
#include<fstream>
#endif

#include "RooGlobalFunc.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TAxis.h"
#include "CITCommon/CommonData/interface/ZeeEventTree.h"
#include "CITHZZ/CommonCode/CommonDefs.hh"
#include "tdrstyle.C"

using namespace RooFit;

void drawPlot(string workspaceFilename, string inputFilename, string outFilename,  
                      Int_t EnergyType,
                      Int_t CategoryBin,
                      double minMass, double maxMass, 
                      double mean_bw, double gamma_bw, double cutoff_cb, double power_cb,
                      const char *plotOpt, const int nbins);

void PlotZMassScaleAndResolutionFit(string workspaceFilename, string inputFilename, string outFilename, Int_t EnergyType, Int_t CategoryBin) {

  // Define Fit Inputs and Call Fit
  double minMass = 70;
  double maxMass = 110;
  double mean_bw = 91.1876;
  double gamma_bw = 2.4952;
  double cutoff_cb = 1.0;
//double power_cb = 1.40;		// Use to fix some fits
  double power_cb = 2.45;
  const char *plotOpt = "NEU";
  const int nbins = 40;

  // Call the fitting program and output a workspace with a root file
  // of the model and data as well as a pdf of the fit
  drawPlot(workspaceFilename, inputFilename, outFilename, EnergyType, CategoryBin, minMass,  maxMass,  mean_bw,  gamma_bw,  cutoff_cb, power_cb, plotOpt, nbins);

}
//______________________________________________________________


void drawPlot(string workspaceFilename, string inputFilename, string outFilename, 
                      Int_t EnergyType,
                      Int_t CategoryBin,
                      double minMass, double maxMass, 
                      double mean_bw, double gamma_bw, double cutoff_cb, double power_cb, 
                      const char* plotOpt, const int nbins) {


  TFile *workspaceFile = new TFile(workspaceFilename.c_str(), "read");
  RooWorkspace* w = (RooWorkspace*)workspaceFile->Get("ZeeMassScaleAndResolutionFit");
 
  


  //Create Data Set
  RooRealVar *mass = (RooRealVar*)w->var("mass");

  // Reading everything from root tree instead
  citana::ZeeEventTree *zeeTree = new citana::ZeeEventTree();
  zeeTree->LoadTree(inputFilename.c_str(), citana::ZeeEventTree::kCITZeeEvent);
  
  RooArgSet zMassArgSet(*mass);
  RooDataSet* data = new RooDataSet("data", "ntuple parameters", zMassArgSet);

  for (int i = 0; i < zeeTree->tree_->GetEntries(); i++) {
    zeeTree->tree_->GetEntry(i);

    //*************************************************************************
    //Electron Selection
    //*************************************************************************
    if (!(zeeTree->fEle1PassHZZICHEP2012 == 1 && zeeTree->fEle2PassHZZICHEP2012 == 1)) continue;

    //*************************************************************************
    //Compute electron four vector;
    //*************************************************************************
    double ele1pt = zeeTree->fEle1Pt;
    double ele2pt = zeeTree->fEle2Pt;
    if (EnergyType == 1) {
      ele1pt = zeeTree->fEle1EnergyRegressionWithTrkVarTwoPtBins / TMath::CosH(zeeTree->fEle1Eta);
      ele2pt = zeeTree->fEle2EnergyRegressionWithTrkVarTwoPtBins / TMath::CosH(zeeTree->fEle2Eta);
    }
    else if (EnergyType == 2) {
      ele1pt = zeeTree->fEle1EnergyRegressionWithTrkVar / TMath::CosH(zeeTree->fEle1Eta);
      ele2pt = zeeTree->fEle2EnergyRegressionWithTrkVar / TMath::CosH(zeeTree->fEle2Eta);
    }
    else if (EnergyType == 3) {
      ele1pt = zeeTree->fEle1EnergyRegressionNoTrkVarTwoPtBins / TMath::CosH(zeeTree->fEle1Eta);
      ele2pt = zeeTree->fEle2EnergyRegressionNoTrkVarTwoPtBins / TMath::CosH(zeeTree->fEle2Eta);
    }
    else if (EnergyType == 4) {
      ele1pt = zeeTree->fEle1EnergyRegressionNoTrkVar / TMath::CosH(zeeTree->fEle1Eta);
      ele2pt = zeeTree->fEle2EnergyRegressionNoTrkVar / TMath::CosH(zeeTree->fEle2Eta);
    }
    TLorentzVector ele1FourVector;
    ele1FourVector.SetPtEtaPhiM(ele1pt, zeeTree->fEle1Eta, zeeTree->fEle1Phi, ELECTRONMASS);
    TLorentzVector ele2FourVector;
    ele2FourVector.SetPtEtaPhiM(ele2pt, zeeTree->fEle2Eta, zeeTree->fEle2Phi, ELECTRONMASS);
    
    //*************************************************************************
    //pt and eta cuts on electron
    //*************************************************************************
    if (! (ele1pt > 7 && ele2pt > 7 
           && fabs( zeeTree->fEle1Eta) < 2.5 
           && fabs( zeeTree->fEle2Eta) < 2.5 )) continue;

    //*************************************************************************
    //pt bins and eta bins
    //*************************************************************************
    Int_t Ele1PtBin = -1;
    Int_t Ele1EtaBin = -1;
    Int_t Ele2PtBin = -1;
    Int_t Ele2EtaBin = -1;
    if (ele1pt > 10 && ele1pt < 20) Ele1PtBin = 0;
    else if (ele1pt < 30) Ele1PtBin = 1;
    else if (ele1pt < 40) Ele1PtBin = 2;
    else Ele1PtBin = 3;
    if (ele2pt > 10 && ele2pt < 20) Ele2PtBin = 0;
    else if (ele2pt < 30) Ele2PtBin = 1;
    else if (ele2pt < 40) Ele2PtBin = 2;
    else Ele2PtBin = 3;
    if (fabs(zeeTree->fEle1SCEta) < 1.0) Ele1EtaBin = 0;
    else if (fabs(zeeTree->fEle1SCEta) < 1.479) Ele1EtaBin = 1;
    else Ele1EtaBin = 2;
    if (fabs(zeeTree->fEle2SCEta) < 1.0) Ele2EtaBin = 0;
    else if (fabs(zeeTree->fEle2SCEta) < 1.479) Ele2EtaBin = 1;
    else Ele2EtaBin = 2;

    if (CategoryBin == 0) { 
      if (!(Ele1EtaBin == 0 && Ele2EtaBin == 0)) continue; 
    }
    else if (CategoryBin == 1) {
      if (!(Ele1EtaBin == 1 && Ele2EtaBin == 1)) continue; 
    }
    else if (CategoryBin == 2) {
      if (!(Ele1EtaBin == 2 && Ele2EtaBin == 2)) continue;
    }
    
    //*************************************************************************
    // restrict range of mass
    //*************************************************************************
    double zMass = (ele1FourVector+ele2FourVector).M();
    if (zMass < minMass || zMass > maxMass) continue;

    //*************************************************************************
    //set mass variable
    //*************************************************************************
    zMassArgSet.setRealValue("mass", zMass);

    data->add(zMassArgSet);
  }

  cout << "dataset size: " << data->numEntries() << endl;


  RooRealVar *cbBias = (RooRealVar*)w->var("#Deltam_{CB}");
  RooRealVar *cbSigma = (RooRealVar*)w->var("sigma_{CB}");
  RooRealVar *cbCut   = (RooRealVar*)w->var("a_{CB}");
  RooRealVar *cbPower = (RooRealVar*)w->var("n_{CB}");

//   // Now if it's a restricted fit, fix values of cbCut and cbPower to MC values.
//   if (isRestricted) {
//     cbCut.setConstant(kTRUE);
//     cbPower.setConstant(kTRUE);
//   }

  // Mass model for signal electrons p.d.f.
  RooAddPdf *model = (RooAddPdf*)w->pdf("model");


  TCanvas* c = new TCanvas("c","c", 0,0,800,600);

  //========================== Plotting  ============================
  //Create a frame
  RooPlot* plot = mass->frame(Range(minMass,maxMass),Bins(nbins));
  // Add data and model to canvas
  plot->SetTitle("");
  plot->GetYaxis()->SetTitleOffset(1.4);
  plot->GetXaxis()->SetTitle("m_{ee} (GeV/c^{2})");

  data->plotOn(plot);
  model->plotOn(plot);
//   model->paramOn(plot, Format(plotOpt, AutoPrecision(1)), Parameters(RooArgSet(*cbBias, *cbSigma, *cbCut, *cbPower)), Layout(0.60,0.90,0.90));
  model->paramOn(plot, Format(plotOpt, AutoPrecision(1)), Parameters(RooArgSet(*cbBias, *cbSigma, *cbCut, *cbPower)), Layout(0.12,0.38,0.60));
  plot->getAttText()->SetTextSize(.025);
  plot->Draw();

  // Print Fit Values
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(.04);
  tex->SetTextFont(2);
  tex->DrawLatex(0.195,0.775, "Run 2012A/B");
  tex->Draw();
//   tex->SetTextSize(0.022);
//   tex->DrawLatex(0.195, 0.75, "Z #rightarrow ee^{+}");
//   tex->SetTextSize(0.024);
//   tex->DrawLatex(0.645, 0.59, Form("BW Mean = %.2f GeV/c^{2}", bwMean.getVal()));
//   tex->DrawLatex(0.645, 0.54, Form("BW #sigma = %.2f GeV/c^{2}", bwWidth.getVal()));
  c->Update();
  c->SaveAs((outFilename + ".gif").c_str());


}
