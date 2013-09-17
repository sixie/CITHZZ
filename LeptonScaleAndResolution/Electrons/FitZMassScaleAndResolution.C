

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
#include "CITCommon/CommonData/interface/ZeeEventTree.h"
#include "CITHZZ/CommonCode/CommonDefs.hh"
#include "tdrstyle.C"
 
using namespace RooFit;

RooWorkspace* makefit(string inputFilename, string outFilename,  
                      Int_t EnergyType,
                      Int_t CategoryBin,
                      double minMass, double maxMass, 
                      double mean_bw, double gamma_bw, double cutoff_cb, double power_cb,
                      const char *plotOpt, const int nbins, bool isRestricted);

void FitZMassScaleAndResolution(string inputFilename, string outFilename, Int_t EnergyType, Int_t CategoryBin, bool isRestricted, string mcWorkspaceFilename = "") {

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

  // Now if restricted fit is desired, read values from workspace
  if (isRestricted) {
    TFile mcWorkspaceFile(mcWorkspaceFilename.c_str());
    RooWorkspace* mcWorkspace = (RooWorkspace*) mcWorkspaceFile.Get("ZeeMassScaleAndResolutionFit");
    cutoff_cb = mcWorkspace->var("a_{CB}")->getVal();
    power_cb = mcWorkspace->var("n_{CB}")->getVal();
  } 

  // Call the fitting program and output a workspace with a root file
  // of the model and data as well as a pdf of the fit
  RooWorkspace *w =  makefit(inputFilename, outFilename, EnergyType, CategoryBin, minMass,  maxMass,  mean_bw,  gamma_bw,  cutoff_cb, power_cb, plotOpt, nbins, isRestricted);
  w->writeToFile((outFilename + ".root").c_str());
  delete w;

}
//______________________________________________________________


RooWorkspace* makefit(string inputFilename, string outFilename, 
                      Int_t EnergyType,
                      Int_t CategoryBin,
                      double minMass, double maxMass, 
                      double mean_bw, double gamma_bw, double cutoff_cb, double power_cb, 
                      const char* plotOpt, const int nbins, bool isRestricted) {

//   gROOT->ProcessLine(".L tdrstyle.C");
//   setTDRStyle();
//   gStyle->SetPadRightMargin(0.05);

  //Create Data Set
  RooRealVar mass("mass","m(EE)",minMass,maxMass,"GeV/c^{2}");
  mass.setRange(80,100);

  // Reading everything from root tree instead
  citana::ZeeEventTree *zeeTree = new citana::ZeeEventTree();
  zeeTree->LoadTree(inputFilename.c_str(), citana::ZeeEventTree::kCITZeeEvent);
  
  RooArgSet zMassArgSet(mass);
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
      ele1pt = zeeTree->fEle1EnergyRegressionV0 / TMath::CosH(zeeTree->fEle1Eta);
      ele2pt = zeeTree->fEle2EnergyRegressionV0 / TMath::CosH(zeeTree->fEle2Eta);
    }
    else if (EnergyType == 2) {
      ele1pt = zeeTree->fEle1EnergyRegressionV1 / TMath::CosH(zeeTree->fEle1Eta);
      ele2pt = zeeTree->fEle2EnergyRegressionV1 / TMath::CosH(zeeTree->fEle2Eta);
    }
    else if (EnergyType == 3) {
      ele1pt = zeeTree->fEle1EnergyRegressionV2 / TMath::CosH(zeeTree->fEle1Eta);
      ele2pt = zeeTree->fEle2EnergyRegressionV2 / TMath::CosH(zeeTree->fEle2Eta);
    }
    TLorentzVector ele1FourVector;
    ele1FourVector.SetPtEtaPhiM(ele1pt, zeeTree->fEle1Eta, zeeTree->fEle1Phi, ELECTRONMASS);
    TLorentzVector ele2FourVector;
    ele2FourVector.SetPtEtaPhiM(ele2pt, zeeTree->fEle2Eta, zeeTree->fEle2Phi, ELECTRONMASS);
    
    //*************************************************************************
    //pt and eta cuts on electron
    //*************************************************************************
    if (! (ele1pt > 15 && ele2pt > 15
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

//   // Closing file
//   treeFile->Close();
  //====================== Parameters===========================

  //Crystal Ball parameters
//   RooRealVar cbBias ("#Deltam_{CB}", "CB Bias", -.01, -10, 10, "GeV/c^{2}");
//   RooRealVar cbSigma("sigma_{CB}", "CB Width", 1.7, 0.8, 5.0, "GeV/c^{2}");
//   RooRealVar cbCut  ("a_{CB}","CB Cut", 1.05, 1.0, 3.0);
//   RooRealVar cbPower("n_{CB}","CB Order", 2.45, 0.1, 20.0);
  RooRealVar cbBias ("#Deltam_{CB}", "CB Bias", -.01, -10, 10, "GeV/c^{2}");
  RooRealVar cbSigma("sigma_{CB}", "CB Width", 1.5, 0.8, 5.0, "GeV/c^{2}");
  RooRealVar cbCut  ("a_{CB}","CB Cut", 1.0, 1.0, 3.0);
  RooRealVar cbPower("n_{CB}","CB Order", 2.5, 0.1, 20.0);
  cbCut.setVal(cutoff_cb);
  cbPower.setVal(power_cb);

  // Just checking
  //cbCut.Print();
  //cbPower.Print();

  // Now if it's a restricted fit, fix values of cbCut and cbPower to MC values.
  if (isRestricted) {
    cbCut.setConstant(kTRUE);
    cbPower.setConstant(kTRUE);
  }

  //Breit_Wigner parameters
  RooRealVar bwMean("m_{Z}","BW Mean", 91.1876, "GeV/c^{2}");
  bwMean.setVal(mean_bw);
  RooRealVar bwWidth("#Gamma_{Z}", "BW Width", 2.4952, "GeV/c^{2}");
  bwWidth.setVal(gamma_bw);

  // Fix the Breit-Wigner parameters to PDG values
  bwMean.setConstant(kTRUE);
  bwWidth.setConstant(kTRUE);

  // Exponential Background parameters
  RooRealVar expRate("#lambda_{exp}", "Exponential Rate", -0.064, -1, 1);
  RooRealVar c0("c_{0}", "c0", 1., 0., 50.);

  //Number of Signal and Background events
  RooRealVar nsig("N_{S}", "# signal events", 524, 0.1, 10000000000.);
  RooRealVar nbkg("N_{B}", "# background events", 43, 1., 10000000.);

  //============================ P.D.F.s=============================

  // Mass signal for two decay electrons p.d.f.
  RooBreitWigner bw("bw", "bw", mass, bwMean, bwWidth);
  RooCBShape  cball("cball", "Crystal Ball", mass, cbBias, cbSigma, cbCut, cbPower);
  RooFFTConvPdf BWxCB("BWxCB", "bw X crystal ball", mass, bw, cball);

  // Mass background p.d.f.
  RooExponential bg("bg", "exp. background", mass, expRate);

  // Mass model for signal electrons p.d.f.
  RooAddPdf model("model", "signal", RooArgList(BWxCB), RooArgList(nsig));

  TStopwatch t ;
  t.Start() ;
  model.fitTo(*data,FitOptions("mh"),Optimize(0),Timer(1));
  t.Print() ;

  TCanvas* c = new TCanvas("c","Unbinned Invariant Mass Fit", 0,0,800,600);

  //========================== Plotting  ============================
  //Create a frame
  RooPlot* plot = mass.frame(Range(minMass,maxMass),Bins(nbins));
  // Add data and model to canvas
  data->plotOn(plot);
  model.plotOn(plot);
  model.paramOn(plot, Format(plotOpt, AutoPrecision(1)), Parameters(RooArgSet(cbBias, cbSigma, cbCut, cbPower, bwMean, bwWidth, expRate, nsig, nbkg)), Layout(0.66,0.63));
  plot->getAttText()->SetTextSize(.03);
  plot->Draw();

  // Print Fit Values
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(.04);
  tex->SetTextFont(2);
  tex->DrawLatex(0.195,0.875, "CMS ECAL, 2012");
  tex->Draw();
  tex->SetTextSize(0.022);
  tex->DrawLatex(0.195, 0.75, "Z #rightarrow ee^{+}");
  tex->SetTextSize(0.024);
  tex->DrawLatex(0.645, 0.59, Form("BW Mean = %.2f GeV/c^{2}", bwMean.getVal()));
  tex->DrawLatex(0.645, 0.54, Form("BW #sigma = %.2f GeV/c^{2}", bwWidth.getVal()));
  c->Update();
  c->SaveAs((outFilename + ".gif").c_str());

  // tablefile << Form(Outfile + "& $ %f $ & $ %f $ & $ %f $\\ \hline",cbBias.getVal(), cbSigma.getVal(), cbCut.getVal());
  // Output workspace with model and data

  RooWorkspace *w = new RooWorkspace("ZeeMassScaleAndResolutionFit");
  w->import(model);
//  w->import(*data);		// Don't save data together, it's not necessary
  return w;
}
