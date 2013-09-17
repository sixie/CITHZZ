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

#include "TCanvas.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"

#include "tdrstyle.C"

using namespace RooFit;

RooWorkspace* makefit(char *filename, char* FitTitle, char* Outfile, double minMass, double maxMass, double mean_bw, double gamma_bw, double cutoff_cb, const char *plotOpt, const int nbins);

void fitZraw() {
  // ======================================================= //
  //                  Setting Program Mode                   //
  //                                                         //
  //     Auto:  Runs preset text files with file names       //
  //            title and output file names                  //
  //                                                         //
  //     User:  Allows user input of individual file,        //
  //            title and output file names                  //
  // ======================================================= //
  string mode;
  cout << "Choose Mode (Auto/User):"<<endl;
  getline(cin, mode);

  // Declaring tablefile
  ofstream tablefile;

  // For Auto Mode
  if (mode == "Auto") {
    // Line Counters
    string set;
    cout<< "Choose Data Set (tight, mva, mvat, NULL)"<<endl;
    getline(cin, set);
    Int_t nfile = 0;
    Int_t ntitle = 0;
    Int_t noutput = 0;
    // Lines
    string fline;
    string tline;
    string oline;
    string table;
    // Arrays of Lines
    string files[200];
    string titles[200];
    string outputs[200];
    // Open Files
    string datafiles = "fitfiles/"+ set +"datafiles.txt";
    string titlefiles = "fitfiles/"+ set+"titlefiles.txt";
    string outputfiles = "fitfiles/"+ set+"outputfiles.txt";
    string tablefiles = "fitfiles/"+ set+ "tablefiles.txt";
    ifstream dfiles(datafiles.c_str());
    ifstream tfiles(titlefiles.c_str());
    ifstream ofiles(outputfiles.c_str());
    ifstream Tfiles(tablefiles.c_str());

    // Check that files are open
    if (!dfiles.is_open()){
      cout<<"ERROR: Given data file could not be opened. File does not exist"<<endl;
      exit(1);
    }
    if (!tfiles.is_open()){
      cout<<"ERROR: Given title file could not be opened. File does not exist"<<endl;
      exit(1);
    }
    if (!ofiles.is_open()){
      cout<<"ERROR: Given output file could not be opened. File does not exist"<<endl;
      exit(1);
    }
    if (!Tfiles.is_open()){
      cout<<"ERROR: Given Table file could not be opened. File does not exist"<<endl;
      exit(1);
    }

    // Read Files to Arrays and count number of lines
    while (true){
      getline(dfiles, fline);
      if (dfiles.eof()) {
	break;
      }
      files[nfile] = fline;
      nfile++;
    }
    while (true){
      getline(tfiles, tline);
      if (tfiles.eof()) {
	break;
      }
      titles[ntitle] = tline;
      ntitle++;
    }
    while (true){
      getline(ofiles, oline);
      if (ofiles.eof()) {
	break;
      }
      outputs[noutput] = oline;
      noutput++;
    }
    getline(Tfiles, table);
    cout<<ntitle<<endl;
    cout<<nfile<<endl;
    cout<<noutput<<endl;
    // Require all files to have same length
    if (nfile == ntitle && ntitle == noutput) {
      // Define Fit Inputs
      double minMass = 60;
      double maxMass = 120;
      double mean_bw = 91.1876;
      double gamma_bw = 2.4952;
      double cutoff_cb = 1.0;
      const char *plotOpt = "NEU";
      const int nbins = 40;
      string autofile;
      string autoTitle;
      string autooutput;
      char* filename;
      char* FitTitle;
      char* Outfile;
      tablefile.open(table.c_str());
      //Loop through each of the files and call the fitting program
      for (Int_t j=0; j<nfile;j++) {
	autofile = files[j];
	autoTitle = titles[j];
	autooutput = outputs[j];
	filename = (char*) autofile.c_str();
	FitTitle = (char*) autoTitle.c_str();
	Outfile = (char*) autooutput.c_str();
	// Call the fitting program and output a workspace with a root file
	// of the model and data as well as a pdf of the fit
	RooWorkspace *w = makefit(filename, FitTitle, Outfile, minMass,  maxMass,  mean_bw,  gamma_bw,  cutoff_cb,  plotOpt, nbins);
	w->writeToFile(Form("%s.root", Outfile));
	delete w;
      }
      tablefile.close();
    } else {
      cout<<"ERROR: File Sizes Disagree"<<endl;
    }
  }

  // For Interactive User mode:
  if (mode == "User") { 
    // Choose the File Name, Title and Output File for Fit
    string file;
    cout << "Enter File Name:"<<endl;
    getline(cin, file);
    string Title;
    cout << "Enter Fit Title"<<endl;
    getline(cin, Title);
    string output;
    cout << "Enter Output File Name"<<endl;
    getline(cin, output);
    string table;
    cout << "Enter Table File Name"<<endl;
    getline(cin, table);
    char* filename = (char*) file.c_str();
    char* FitTitle = (char*) Title.c_str();
    char* Outfile = (char*) output.c_str();
    tablefile.open(table.c_str());
    // Define Fit Inputs and Call Fit
    double minMass = 60;
    double maxMass = 120;
    double mean_bw = 91.1876;
    double gamma_bw = 2.4952;
    double cutoff_cb = 1.0;
    const char *plotOpt = "NEU";
    const int nbins = 40;
    // Call the fitting program and output a workspace with a root file
    // of the model and data as well as a pdf of the fit
    RooWorkspace *w =  makefit(filename, FitTitle, Outfile, minMass,  maxMass,  mean_bw,  gamma_bw,  cutoff_cb,  plotOpt, nbins);
    w->writeToFile(Form("%s.root", Outfile));
    delete w;
  }  
  tablefile.close();
}
//______________________________________________________________


RooWorkspace* makefit(char* filename, char* FitTitle, char* Outfile, double minMass, double maxMass, double mean_bw, double gamma_bw, double cutoff_cb, const char* plotOpt, const int nbins) {
  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();
  gStyle->SetPadRightMargin(0.05);

  //Create Data Set
  RooRealVar mass("mass","m(EE)",minMass,maxMass,"GeV/c^{2}");
 // RooDataSet *data = RooDataSet::read(filename,RooArgSet(mass));

  // Reading everything from root tree instead
  TFile* treeFile = new TFile(filename);
  TTree* zeeTree = (TTree*) treeFile->Get("ZeeEvent");

  float zMass;
  zeeTree->SetBranchAddress("mass", &zMass);

  RooArgSet zMassArgSet(mass);
  RooDataSet* data = new RooDataSet("data", "ntuple parameters", zMassArgSet);

  int nEntries = zeeTree->GetEntries();
  for (int i = 0; i < nEntries; i++) {
    zeeTree->GetEntry(i);

    // Cutting manually the values outside of range (ugly but works...)
    if (zMass < minMass || zMass > maxMass) continue;

    zMassArgSet.setRealValue("mass", zMass);

    data->add(zMassArgSet);
  }

  // Closing file
  treeFile->Close();
  //====================== Parameters===========================

  //Crystal Ball parameters
  RooRealVar cbBias ("#Deltam_{CB}", "CB Bias", -.01, -10, 10, "GeV/c^{2}");
  RooRealVar cbSigma("sigma_{CB}", "CB Width", 1.7, 0.02, 5.0, "GeV/c^{2}");
  RooRealVar cbCut  ("a_{CB}","CB Cut", 1.05, 0.1, 3.0);
  RooRealVar cbPower("n_{CB}","CB Order", 2.45, 0.1, 20.0);
  cbCut.setVal(cutoff_cb);

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
  tex->DrawLatex(0.195, 0.81, FitTitle);
  tex->DrawLatex(0.195, 0.75, "Z #rightarrow ee^{+}");
  tex->SetTextSize(0.024);
  tex->DrawLatex(0.645, 0.59, Form("BW Mean = %.2f GeV/c^{2}", bwMean.getVal()));
  tex->DrawLatex(0.645, 0.54, Form("BW #sigma = %.2f GeV/c^{2}", bwWidth.getVal()));
  c->Update();
  c->Print(Outfile);
  // tablefile << Form(Outfile + "& $ %f $ & $ %f $ & $ %f $\\ \hline",cbBias.getVal(), cbSigma.getVal(), cbCut.getVal());
  // Output workspace with model and data
  RooWorkspace *w = new RooWorkspace("zfit");
  w->import(model);
  w->import(*data);
  return w;
}
