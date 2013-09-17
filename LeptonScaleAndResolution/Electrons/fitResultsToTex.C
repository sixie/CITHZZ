#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TGraphAsymmErrors.h>      // graphs
#include <TH2F.h>                   // 2D histograms
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <string>		    // string class for handling filename
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <sstream>
#include <fstream>                  // functions for file I/O

#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#endif

void fitResultsToTex (string dataNameModel, string mcNameModel, string outName) {
  // Reads fit results from .root files containing RooWorkspaces and creates a .tex file with compairing tables
  // Each variable in the result is composed of three different tables

  // Parameters to be read from the input files
  float deltaMC_eta0[5], deltaData_eta0[5], sigmaMC_eta0[5], sigmaData_eta0[5];
  float deltaMC_eta1[5], deltaData_eta1[5], sigmaMC_eta1[5], sigmaData_eta1[5];
  float deltaMC_eta2[5], deltaData_eta2[5], sigmaMC_eta2[5], sigmaData_eta2[5];

  float deltaMCError_eta0[5], deltaDataError_eta0[5], sigmaMCError_eta0[5], sigmaDataError_eta0[5];
  float deltaMCError_eta1[5], deltaDataError_eta1[5], sigmaMCError_eta1[5], sigmaDataError_eta1[5];
  float deltaMCError_eta2[5], deltaDataError_eta2[5], sigmaMCError_eta2[5], sigmaDataError_eta2[5];

  // Opening files
  for (int i = 0; i < 5; i++) {

    string label;
    stringstream labelStream;
    labelStream << "_EnergyType" << i;
    label = labelStream.str();
    TFile dataFile_eta0((dataNameModel + label + "_CategoryBin0.root").c_str());
    TFile dataFile_eta1((dataNameModel + label + "_CategoryBin1.root").c_str());
    TFile dataFile_eta2((dataNameModel + label + "_CategoryBin2.root").c_str());
    TFile mcFile_eta0((mcNameModel + label + "_CategoryBin0.root").c_str());
    TFile mcFile_eta1((mcNameModel + label + "_CategoryBin1.root").c_str());
    TFile mcFile_eta2((mcNameModel + label + "_CategoryBin2.root").c_str());

    RooWorkspace* dataW_eta0 = (RooWorkspace*) dataFile_eta0.Get("ZeeMassScaleAndResolutionFit");
    RooWorkspace* dataW_eta1 = (RooWorkspace*) dataFile_eta1.Get("ZeeMassScaleAndResolutionFit");
    RooWorkspace* dataW_eta2 = (RooWorkspace*) dataFile_eta2.Get("ZeeMassScaleAndResolutionFit");
    RooWorkspace* mcW_eta0 = (RooWorkspace*) mcFile_eta0.Get("ZeeMassScaleAndResolutionFit");
    RooWorkspace* mcW_eta1 = (RooWorkspace*) mcFile_eta1.Get("ZeeMassScaleAndResolutionFit");
    RooWorkspace* mcW_eta2 = (RooWorkspace*) mcFile_eta2.Get("ZeeMassScaleAndResolutionFit");

    assert(dataW_eta0);
    assert(dataW_eta1);
    assert(dataW_eta2);
    assert(mcW_eta0);
    assert(mcW_eta1);
    assert(mcW_eta2);

    deltaData_eta0[i] = dataW_eta0->var("#Deltam_{CB}")->getVal(); deltaDataError_eta0[i] = dataW_eta0->var("#Deltam_{CB}")->getError();
    deltaMC_eta0[i] = mcW_eta0->var("#Deltam_{CB}")->getVal(); deltaMCError_eta0[i] = mcW_eta0->var("#Deltam_{CB}")->getError();
    sigmaData_eta0[i] = dataW_eta0->var("sigma_{CB}")->getVal(); sigmaDataError_eta0[i] = dataW_eta0->var("sigma_{CB}")->getError();
    sigmaMC_eta0[i] = mcW_eta0->var("sigma_{CB}")->getVal(); sigmaMCError_eta0[i] = mcW_eta0->var("sigma_{CB}")->getError();

    deltaData_eta1[i] = dataW_eta1->var("#Deltam_{CB}")->getVal(); deltaDataError_eta1[i] = dataW_eta1->var("#Deltam_{CB}")->getError();
    deltaMC_eta1[i] = mcW_eta1->var("#Deltam_{CB}")->getVal(); deltaMCError_eta1[i] = mcW_eta1->var("#Deltam_{CB}")->getError();
    sigmaData_eta1[i] = dataW_eta1->var("sigma_{CB}")->getVal(); sigmaDataError_eta1[i] = dataW_eta1->var("sigma_{CB}")->getError();
    sigmaMC_eta1[i] = mcW_eta1->var("sigma_{CB}")->getVal(); sigmaMCError_eta1[i] = mcW_eta1->var("sigma_{CB}")->getError();

    deltaData_eta2[i] = dataW_eta2->var("#Deltam_{CB}")->getVal(); deltaDataError_eta2[i] = dataW_eta2->var("#Deltam_{CB}")->getError();
    deltaMC_eta2[i] = mcW_eta2->var("#Deltam_{CB}")->getVal(); deltaMCError_eta2[i] = mcW_eta2->var("#Deltam_{CB}")->getError();
    sigmaData_eta2[i] = dataW_eta2->var("sigma_{CB}")->getVal(); sigmaDataError_eta2[i] = dataW_eta2->var("sigma_{CB}")->getError();
    sigmaMC_eta2[i] = mcW_eta2->var("sigma_{CB}")->getVal(); sigmaMCError_eta2[i] = mcW_eta2->var("sigma_{CB}")->getError();
  }

  // Opening TeX file
  ofstream texfile;
  texfile.open(outName.c_str());
  assert(texfile.is_open());

  texfile << "\\documentclass{article}" << endl;
  texfile << "\\begin{document}" << endl; 


  //***********************************************************************************
  // Compare different regressions
  //***********************************************************************************

  texfile << "\\begin{table}[!ht]" << endl;
  texfile << "\\begin{center} " << endl;
  texfile << "$ 0.0 < |\\eta| < 1.0 $ \\\\" << endl;

  texfile << "\\begin{tabular}{|c|c|c|c|c|}" << endl;
  texfile << "\\hline" << endl;

  texfile << "Standard/Regression type   &   $\\Delta m$ (MC)   &   $\\Delta m$ (Data)  &  $\\sigma$ (MC)  &   $\\sigma$ (Data)  \\\\  ";
  texfile << "\\hline" << endl;

  // Now writing the values for the standard and the regressions
  for (int i = 0; i < 5; i++) {
    // Preparing string for plotting
    string versionName;
    if (i == 0) versionName = "Standard";
    else if (i == 1) versionName = "No trk var";
    else if (i == 2) versionName = "With trk var";
    else if (i == 3) versionName = "No trk var, two pt bins";
    else if (i == 4) versionName = "With trk var, two pt bins";

    // Preparing line to put in latex table

    string tableLine = Form((versionName + " & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ \\\\").c_str(),
		    deltaMC_eta0[i], deltaMCError_eta0[i], deltaData_eta0[i], deltaDataError_eta0[i],
		    sigmaMC_eta0[i], sigmaMCError_eta0[i], sigmaData_eta0[i], sigmaDataError_eta0[i]);

    texfile << tableLine;
    texfile << "\\hline" << endl;
  }

  texfile << "\\end{tabular}" << endl;
  texfile << "\\end{center}" << endl;
  texfile << "\\end{table}" << endl;

  texfile << endl;

  texfile << "\\begin{table}[!ht]" << endl;
  texfile << "\\begin{center} " << endl;
  texfile << "$ 1.0 < |\\eta| < 1.479 $ \\\\";

  texfile << "\\begin{tabular}{|c|c|c|c|c|}" << endl;
  texfile << "\\hline" << endl;

  texfile << "Standard/Regression type   &   $\\Delta m$ (MC)   &   $\\Delta m$ (Data)  &  $\\sigma$ (MC)  &   $\\sigma$ (Data)  \\\\  ";
  texfile << "\\hline" << endl;

  // Now writing the values for the standard and the regressions
  for (int i = 0; i < 5; i++) {
    // Preparing string for plotting
    string versionName;
    if (i == 0) versionName = "Standard";
    else if (i == 1) versionName = "No trk var";
    else if (i == 2) versionName = "With trk var";
    else if (i == 3) versionName = "No trk var, two pt bins";
    else if (i == 4) versionName = "With trk var, two pt bins";

    // Preparing line to put in latex table

    string tableLine = Form((versionName + " & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ \\\\").c_str(),
		    deltaMC_eta1[i], deltaMCError_eta1[i], deltaData_eta1[i], deltaDataError_eta1[i],
		    sigmaMC_eta1[i], sigmaMCError_eta1[i], sigmaData_eta1[i], sigmaDataError_eta1[i]);

    texfile << tableLine;
    texfile << "\\hline" << endl;
  }

  texfile << "\\end{tabular}" << endl;
  texfile << "\\end{center}" << endl;
  texfile << "\\end{table}" << endl;

  texfile << endl;

  texfile << "\\begin{table}[!ht]" << endl;
  texfile << "\\begin{center} " << endl;
  texfile << "$ |\\eta| > 1.479 $ \\\\";

  texfile << "\\begin{tabular}{|c|c|c|c|c|}" << endl;
  texfile << "\\hline" << endl;

  texfile << "Standard/Regression type   &   $\\Delta m$ (MC)   &   $\\Delta m$ (Data)  &  $\\sigma$ (MC)  &   $\\sigma$ (Data)  \\\\  ";
  texfile << "\\hline" << endl;

  // Now writing the values for the standard and the regressions
  for (int i = 0; i < 5; i++) {
    // Preparing string for plotting
    string versionName;
    if (i == 0) versionName = "Standard";
    else if (i == 1) versionName = "No trk var";
    else if (i == 2) versionName = "With trk var";
    else if (i == 3) versionName = "No trk var, two pt bins";
    else if (i == 4) versionName = "With trk var, two pt bins";

    // Preparing line to put in latex table

    string tableLine = Form((versionName + " & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ \\\\").c_str(),
		    deltaMC_eta2[i], deltaMCError_eta2[i], deltaData_eta2[i], deltaDataError_eta2[i],
		    sigmaMC_eta2[i], sigmaMCError_eta2[i], sigmaData_eta2[i], sigmaDataError_eta2[i]);

    texfile << tableLine;
    texfile << "\\hline" << endl;
  }

  texfile << "\\end{tabular}" << endl;
  texfile << "\\end{center}" << endl;
  texfile << "\\end{table}" << endl;



  //***********************************************************************************
  // Group by regression type 
  //***********************************************************************************

  for (int t = 0; t < 5; t++) {

    texfile << "\\begin{table}[!ht]" << endl;
    texfile << "\\begin{center} " << endl;

    if (t== 0) {
      texfile << " Standard Electron Momentum \\\\" << endl;
    } else if (t==1) {
      texfile << " Regression without trk var \\\\" << endl;
    }else if (t==2) {
      texfile << " Regression with trk var \\\\" << endl;
    }else if (t==3) {
      texfile << " Regression without trk var, two pt bins \\\\" << endl;
    }else if (t==4) {
      texfile << " Regression with trk var, two pt bins \\\\" << endl;
    } 

    texfile << "\\begin{tabular}{|c|c|c|c|c|}" << endl;
    texfile << "\\hline" << endl;
    texfile << "Bin   &   $\\Delta m$ (MC)   &   $\\Delta m$ (Data)  &  $\\sigma$ (MC)  &   $\\sigma$ (Data)  \\\\  ";
    texfile << "\\hline" << endl;

    for (int i = 0; i < 3; i++) {
      string binLabel;
      if (i == 0) binLabel = "$0.0 \\le |\\eta| < 1.0$";
      else if (i == 1) binLabel = "$1.0 \\le |\\eta| < 1.479$";
      else if (i == 2) binLabel = "$1.479 \\le |\\eta| < 2.5$";

      string tableLine;
      if (i==0) {
        tableLine = Form((binLabel + " & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ \\\\").c_str(),
                         deltaMC_eta0[t], deltaMCError_eta0[t], deltaData_eta0[t], deltaDataError_eta0[t],
                         sigmaMC_eta0[t], sigmaMCError_eta0[t], sigmaData_eta0[t], sigmaDataError_eta0[t]);
      } else if (i==1) {
        tableLine = Form((binLabel + " & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ \\\\").c_str(),
                         deltaMC_eta1[t], deltaMCError_eta1[t], deltaData_eta1[t], deltaDataError_eta1[t],
                         sigmaMC_eta1[t], sigmaMCError_eta1[t], sigmaData_eta1[t], sigmaDataError_eta1[t]);

      } else if (i==2) {
        tableLine = Form((binLabel + " & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ & $%4.3f \\pm %4.3f$ \\\\").c_str(),
                         deltaMC_eta2[t], deltaMCError_eta2[t], deltaData_eta2[t], deltaDataError_eta2[t],
                         sigmaMC_eta2[t], sigmaMCError_eta2[t], sigmaData_eta2[t], sigmaDataError_eta2[t]);

      }
      texfile << tableLine;
      texfile << "\\hline" << endl;
    }

    texfile << "\\end{tabular}" << endl;
    texfile << "\\caption{caption}" << endl;
    texfile << "\\label{tab:label}" << endl;


    texfile << "\\end{center}" << endl;
    texfile << "\\end{table}" << endl;

    texfile << endl;
    texfile << endl;
  }


  //***********************************************************************************
  // Compute Energy Scale Corrections 
  //***********************************************************************************
  for (int t = 0; t < 5; t++) {

    texfile << "\\begin{table}[!ht]" << endl;
    texfile << "\\begin{center} " << endl;

    if (t== 0) {
      texfile << " Standard Electron Momentum \\\\" << endl;
    } else if (t==1) {
      texfile << " Regression without trk var \\\\" << endl;
    }else if (t==2) {
      texfile << " Regression with trk var \\\\" << endl;
    }else if (t==3) {
      texfile << " Regression without trk var, two pt bins \\\\" << endl;
    }else if (t==4) {
      texfile << " Regression with trk var, two pt bins \\\\" << endl;
    } 

    texfile << "\\begin{tabular}{|c|c|c|}" << endl;
    texfile << "\\hline" << endl;
    texfile << "Bin   &   $\\Delta E/E$ &   $\\sigma_{E}/E$ \\\\  ";
    texfile << "\\hline" << endl;

    for (int i = 0; i < 3; i++) {
      string binLabel;
      if (i == 0) binLabel = "$0.0 \\le |\\eta| < 1.0$";
      else if (i == 1) binLabel = "$1.0 \\le |\\eta| < 1.479$";
      else if (i == 2) binLabel = "$1.479 \\le |\\eta| < 2.5$";

      double deltaEOverE = 0;
      double deltaEOverEError = 0;
      double sigmaEOverE = 0;
      double sigmaEOverEError = 0;
      string tableLine;      
      if (i==0) {
        deltaEOverE = ( deltaData_eta0[t] - deltaMC_eta0[t] ) / ( 91.2 + deltaMC_eta0[t]);
        deltaEOverEError = TMath::Sqrt( pow(deltaDataError_eta0[t],2)+pow(deltaMCError_eta0[t],2)) / (91.2 + deltaMC_eta0[t]);
        sigmaEOverE = TMath::Sqrt(2) * TMath::Sqrt( pow(sigmaData_eta0[t],2) - pow(sigmaMC_eta0[t],2) ) / (91.2 + deltaMC_eta0[t]);
        sigmaEOverEError = TMath::Sqrt(2)*TMath::Sqrt( pow(sigmaData_eta0[t]*sigmaDataError_eta0[t]/TMath::Sqrt( pow(sigmaData_eta0[t],2) - pow(sigmaMC_eta0[t],2)),2) +
                                                       pow(sigmaMC_eta0[t]*sigmaMCError_eta0[t]/TMath::Sqrt( pow(sigmaData_eta0[t],2) - pow(sigmaMC_eta0[t],2)),2) ) / (91.2 + deltaMC_eta0[t]);
      } else if (i==1) {
        deltaEOverE = ( deltaData_eta1[t] - deltaMC_eta1[t] ) / ( 91.2 + deltaMC_eta1[t]);
        deltaEOverEError = TMath::Sqrt( pow(deltaDataError_eta1[t],2)+pow(deltaMCError_eta1[t],2)) / (91.2 + deltaMC_eta1[t]);
        sigmaEOverE = TMath::Sqrt(2) * TMath::Sqrt( pow(sigmaData_eta1[t],2) - pow(sigmaMC_eta1[t],2) ) / (91.2 + deltaMC_eta1[t]);
        sigmaEOverEError = TMath::Sqrt(2)*TMath::Sqrt( pow(sigmaData_eta1[t]*sigmaDataError_eta1[t]/TMath::Sqrt( pow(sigmaData_eta1[t],2) - pow(sigmaMC_eta1[t],2)),2) +
                                                       pow(sigmaMC_eta1[t]*sigmaMCError_eta1[t]/TMath::Sqrt( pow(sigmaData_eta1[t],2) - pow(sigmaMC_eta1[t],2)),2) ) / (91.2 + deltaMC_eta1[t]);
      } else if (i==2) {
        deltaEOverE = ( deltaData_eta2[t] - deltaMC_eta2[t] ) / ( 91.2 + deltaMC_eta2[t]);
        deltaEOverEError = TMath::Sqrt( pow(deltaDataError_eta2[t],2)+pow(deltaMCError_eta2[t],2)) / (91.2 + deltaMC_eta2[t]);
        sigmaEOverE = TMath::Sqrt(2) * TMath::Sqrt( pow(sigmaData_eta2[t],2) - pow(sigmaMC_eta2[t],2) ) / (91.2 + deltaMC_eta2[t]);
        sigmaEOverEError = TMath::Sqrt(2)*TMath::Sqrt( pow(sigmaData_eta2[t]*sigmaDataError_eta2[t]/TMath::Sqrt( pow(sigmaData_eta2[t],2) - pow(sigmaMC_eta2[t],2)),2) +
                                                       pow(sigmaMC_eta2[t]*sigmaMCError_eta2[t]/TMath::Sqrt( pow(sigmaData_eta2[t],2) - pow(sigmaMC_eta2[t],2)),2) ) / (91.2 + deltaMC_eta2[t]);
      }

      if (i==0) {
        tableLine = Form((binLabel + " & $%4.5f \\pm %4.5f$ & $%4.4f \\pm %4.4f$ \\\\").c_str(),
                         deltaEOverE, deltaEOverEError, sigmaEOverE, sigmaEOverEError);
      } else {
        tableLine = Form((binLabel + " & $%4.4f \\pm %4.4f$ & $%4.4f \\pm %4.4f$ \\\\").c_str(),
                         deltaEOverE, deltaEOverEError, sigmaEOverE, sigmaEOverEError);
      }
      texfile << tableLine;
      texfile << "\\hline" << endl;
    }

    texfile << "\\end{tabular}" << endl;
    texfile << "\\caption{caption}" << endl;
    texfile << "\\label{tab:label}" << endl;


    texfile << "\\end{center}" << endl;
    texfile << "\\end{table}" << endl;

    texfile << endl;
    texfile << endl;
  }



  //***********************************************************************************
  // Compute Energy Scale Corrections 
  //***********************************************************************************






  texfile << "\\end{document}" << endl;
  texfile.close();
  cout << outName << " created" << endl;
} 
