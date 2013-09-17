

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
#include "RooDataHist.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"

#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "CITCommon/CommonData/interface/ZeeTreeHZZ4l.h"
#include "CITHZZ/CommonCode/CommonDefs.hh"
#include "CITHZZ/LeptonScaleAndResolution/Electrons/tdrstyle.C"
 
using namespace RooFit;

void makefit(string inputFilename, string outFilename,  
	     Int_t CategoryBin, Int_t R9Bin,
	     double minMass, double maxMass, 
	     double mean_bw, double gamma_bw, double cutoff_cb, double power_cb,
	     const char *plotOpt, const int nbins, Int_t isMC);

void FitInclusiveZ(string inputFilename, string outFilename, Int_t CategoryBin, Int_t R9Bin, Int_t isMC) {

  // Define Fit Inputs and Call Fit
  double minMass = 75;
  double maxMass = 105;
  double mean_bw = 91.1876;
  double gamma_bw = 2.4952;
  double cutoff_cb = 1.0;
//double power_cb = 1.40;		// Use to fix some fits
  double power_cb = 2.45;
  const char *plotOpt = "NEU";
  const int nbins = 40;

  // Call the fitting program and output a workspace with a root file
  // of the model and data as well as a pdf of the fit
  makefit(inputFilename, outFilename, CategoryBin, R9Bin, minMass,  maxMass,  mean_bw,  gamma_bw,  cutoff_cb, power_cb, plotOpt, nbins, isMC);

}
//______________________________________________________________


void makefit(string inputFilename, string outFilename, 
	     Int_t CategoryBin, Int_t R9Bin,
	     double minMass, double maxMass, 
	     double mean_bw, double gamma_bw, double cutoff_cb, double power_cb, 
	     const char* plotOpt, const int nbins, Int_t isMC) {

//   gROOT->ProcessLine(".L tdrstyle.C");
//   setTDRStyle();
//   gStyle->SetPadRightMargin(0.05);

  //Create Data Set
  RooRealVar mass("zmass","m(e^{+}e^{-})",minMass,maxMass,"GeV/c^{2}");
  //  mass.setRange(80,100);

  // Reading everything from root tree instead
  TFile *tfile = TFile::Open(inputFilename.c_str());
  TTree *ttree = (TTree*)tfile->Get("zeetree/probe_tree");
  citana::ZeeTreeHZZ4l *zeeTree = new citana::ZeeTreeHZZ4l(ttree);
  
  RooArgSet zMassArgSet(mass);
  RooDataSet* data = new RooDataSet("data", "ntuple parameters", zMassArgSet);

  for (int i = 0; i < zeeTree->fChain->GetEntries(); i++) {
    if(i%100000==0) cout << "Processing Event " << i << endl;
    zeeTree->fChain->GetEntry(i);

    //*************************************************************************
    //Electron Selection
    //*************************************************************************
    // already passed for this tree

    //*************************************************************************
    //Compute electron four vector;
    //*************************************************************************
    double ele1pt = zeeTree->l1pt;
    double ele2pt = zeeTree->l2pt;

    TLorentzVector ele1FourVector;
    ele1FourVector.SetPtEtaPhiM(zeeTree->l1pt, zeeTree->l1eta, zeeTree->l1phi, ELECTRONMASS);
    TLorentzVector ele2FourVector;
    ele2FourVector.SetPtEtaPhiM(zeeTree->l2pt, zeeTree->l2eta, zeeTree->l2phi, ELECTRONMASS);

    
    //*************************************************************************
    //pt and eta cuts on electron
    //*************************************************************************
    if (! (ele1pt > 7 && ele2pt > 7
           && fabs( zeeTree->l1eta) < 2.5 
           && fabs( zeeTree->l2eta) < 2.5 )) continue;

    //*************************************************************************
    //pt bins and eta bins
    //*************************************************************************
    Int_t Ele1PtBin = -1;
    Int_t Ele1EtaBin = -1;
    Int_t Ele2PtBin = -1;
    Int_t Ele2EtaBin = -1;
    Int_t Ele1R9Bin = -1;
    Int_t Ele2R9Bin = -1;
    if (ele1pt > 10 && ele1pt < 20) Ele1PtBin = 0;
    else if (ele1pt < 30) Ele1PtBin = 1;
    else if (ele1pt < 40) Ele1PtBin = 2;
    else Ele1PtBin = 3;
    if (ele2pt > 10 && ele2pt < 20) Ele2PtBin = 0;
    else if (ele2pt < 30) Ele2PtBin = 1;
    else if (ele2pt < 40) Ele2PtBin = 2;
    else Ele2PtBin = 3;
    if (fabs(zeeTree->l1sceta) < 1.0) Ele1EtaBin = 0;
    else if (fabs(zeeTree->l1sceta) < 1.479) Ele1EtaBin = 1;
    else if (fabs(zeeTree->l1sceta) < 2.0) Ele1EtaBin = 2;
    else Ele1EtaBin = 3;
    if (fabs(zeeTree->l2sceta) < 1.0) Ele2EtaBin = 0;
    else if (fabs(zeeTree->l2sceta) < 1.479) Ele2EtaBin = 1;
    else if (fabs(zeeTree->l2sceta) < 2.0) Ele2EtaBin = 2;
    else Ele2EtaBin = 3;
    if (zeeTree->l1r9 > 0.94) Ele1R9Bin = 0;
    else Ele1R9Bin = 1;
    if (zeeTree->l2r9 > 0.94) Ele2R9Bin = 0;
    else Ele2R9Bin = 1;

    if (CategoryBin == 0) { 
      if (!(Ele1EtaBin == 0 && Ele2EtaBin == 0)) continue; 
    }
    else if (CategoryBin == 1) {
      if (!(Ele1EtaBin == 1 && Ele2EtaBin == 1)) continue; 
    }
    else if (CategoryBin == 2) {
      if (!(Ele1EtaBin == 2 && Ele2EtaBin == 2)) continue;
    }
    else if (CategoryBin == 3) {
      if (!(Ele1EtaBin == 3 && Ele2EtaBin == 3)) continue;
    }
    
    if (R9Bin == 0) { 
      if (!(Ele1R9Bin == 0 && Ele2R9Bin == 0)) continue; 
    }
    else if (R9Bin == 1) {
      if (!(Ele1R9Bin == 1 && Ele2R9Bin == 1)) continue; 
    }

    
    //*************************************************************************
    // restrict range of mass
    //*************************************************************************
    double zMass = (ele1FourVector+ele2FourVector).M();
    if (zMass < minMass || zMass > maxMass) continue;

    //*************************************************************************
    //set mass variable
    //*************************************************************************
    zMassArgSet.setRealValue("zmass", zMass);    

    data->add(zMassArgSet);
  }

  // do binned fit to gain time...
  mass.setBins(nbins);
  RooDataHist *bdata = new RooDataHist("data_binned","data_binned", zMassArgSet, *data);

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
  RooRealVar cbSigma("#sigma_{CB}", "CB Width", 1.5, 0.8, 5.0, "GeV/c^{2}");
  RooRealVar cbCut  ("a_{CB}","CB Cut", 1.0, 1.0, 3.0);
  RooRealVar cbPower("n_{CB}","CB Order", 2.0);
  cbCut.setVal(cutoff_cb);
  cbPower.setVal(power_cb);

  // Just checking
  //cbCut.Print();
  //cbPower.Print();

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
  RooFitResult *fitres = model.fitTo(*bdata,Hesse(1),Minos(1),Timer(1),Save(1));
  fitres->SetName("fitres");
  t.Print() ;

  TCanvas* c = new TCanvas("c","Unbinned Invariant Mass Fit", 0,0,800,600);

  //========================== Plotting  ============================
  //Create a frame
  RooPlot* plot = mass.frame(Range(minMass,maxMass),Bins(nbins));
  // Add data and model to canvas
  int col = (isMC ? kAzure+4 : kGreen+1);
  data->plotOn(plot);
  model.plotOn(plot,LineColor(col));
  data->plotOn(plot);
  model.paramOn(plot, Format(plotOpt, AutoPrecision(1)), Parameters(RooArgSet(cbBias, cbSigma, cbCut, cbPower, bwMean, bwWidth, expRate, nsig, nbkg)), Layout(0.15,0.45,0.80));
  plot->getAttText()->SetTextSize(.03);
  plot->SetTitle("");
  plot->Draw();

  // Print Fit Values
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(.1);
  tex->SetTextFont(132);
  //  tex->Draw();
  tex->SetTextSize(0.057);
  if(isMC) tex->DrawLatex(0.65, 0.75, "Z #rightarrow e^{+}e^{-} MC");
  else tex->DrawLatex(0.65, 0.75, "Z #rightarrow e^{+}e^{-} data");
  tex->SetTextSize(0.030);
  tex->DrawLatex(0.645, 0.65, Form("BW Mean = %.2f GeV/c^{2}", bwMean.getVal()));
  tex->DrawLatex(0.645, 0.60, Form("BW #sigma = %.2f GeV/c^{2}", bwWidth.getVal()));
  c->Update();
  c->SaveAs((outFilename + ".pdf").c_str());
  c->SaveAs((outFilename + ".png").c_str());

  // tablefile << Form(Outfile + "& $ %f $ & $ %f $ & $ %f $\\ \hline",cbBias.getVal(), cbSigma.getVal(), cbCut.getVal());
  // Output workspace with model and data

  RooWorkspace *w = new RooWorkspace("ZeeMassScaleAndResolutionFit");
  w->import(model);
  w->import(*bdata);
  w->writeToFile((outFilename + ".root").c_str());  

  TFile *tfileo = TFile::Open((outFilename + ".root").c_str(),"update");
  fitres->Write();
  tfileo->Close();

}
