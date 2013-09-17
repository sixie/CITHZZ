//***************************************************
//2012 Data
//root -l CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/compareRegressionZMasses.C+'("/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_HCP2012.root","DataHCP2012",0,kTRUE,kFALSE,kFALSE)'
//root -l CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/compareRegressionZMasses.C+'("/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Summer12DY53X.root","Summer1253X",0,kTRUE,kFALSE,kTRUE)'
//root -l CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/compareRegressionZMasses.C+'("DataHCP2012",0)'
//root -l CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/compareRegressionZMasses.C+'("Summer1253X",0)'
//***************************************************

//***************************************************
//2011 Data
//root -l CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/compareRegressionZMasses.C+'("/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Run2011.root","Data2011",0,kTRUE,kTRUE,kFALSE)'
//root -l CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/compareRegressionZMasses.C+'("/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Fall11DY42X.root","Fall11DY42X",0,kTRUE,kTRUE,kTRUE)'
//root -l CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/compareRegressionZMasses.C+'("Data2011",0)'
//root -l CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/compareRegressionZMasses.C+'("Fall11DY42X",0)'
//***************************************************




/* ====================================================================

   compareRegressionZMasses

   Compares the performance of electron energy regression in Z->ee ntuples

   USAGE
   compareRegressionZeePerformance (string zeeFilename, string label)

   zeeFilename		file with the ntuples produced by iggsAna/HZZ4l/LeptonScaleAndResolution/MakeZeeEventNtuples.C
   label		filename to which the plots will be saved

   Saves files 	<label>_zMass.gif, <label>_zMass.root
   		<label>_ele1Pt.gif, <label>_ele1Pt.root
   		<label>_ele2Pt.gif, <label>_ele2Pt.root
   _____________________________________________________________________
 */


// Compares the Z masses obtained from electron energies that used different electrons
#include "TSystem.h"
#include "TStyle.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include <string>
#include <cmath>
#include <iostream>
#include "TLatex.h"

// mass definitions
#include "CITHZZ/CommonCode/CommonDefs.hh"
#include "HiggsAna/Utils/LeptonScaleCorrections.hh"

// 
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
#include "RooDataHist.h"

using namespace RooFit;


//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}



//*************************************************************************************************
//Computes Eff Sigma
//*************************************************************************************************


Double_t effSigma(TH1 *hist )
{

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  // Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    cout << "effsigma: Too few entries " << total << endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;
  
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }   
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
  
  return widmin;
  
}


///-----------------------------------------------------------------------------
Double_t effSigma(RooAbsPdf *pdf, RooRealVar *obs, Int_t nbins)
{
  TH1 *hist = pdf->createHistogram(obs->GetName(), nbins);
  hist->Scale(nbins);

  return effSigma( hist);
}



void compareRegressionZMasses (string zeeFilename, string label, Int_t bin, Bool_t doCorr, Bool_t is7TeV, Bool_t isMC) {

  gStyle->SetOptStat(0);

  TRandom3 *random = new TRandom3(9999);

  // Opening ntuple file and reading the tree
  TFile zeeFile(zeeFilename.c_str());
  TTree* zeeTree = (TTree*) zeeFile.Get("ZeeEvent");

  // Reading branches
  UInt_t run;
  float zMass;			// Standard Z mass
  float weight;

  float ele1Eta;		// Electron 1 parameters
  float ele1Phi;
  float ele1Pt;
  float ele1Energy;
  float ele1GenPt;
  float ele1R9;
  float ele1EnergyRegressionV0;
  float ele1EnergyRegressionV1;
  float ele1EnergyRegressionV2;
  Bool_t ele1PassHZZICHEP2012; 

  float ele2Eta;		// Electron 2 parameters
  float ele2Phi;
  float ele2Pt;
  float ele2Energy;
  float ele2GenPt;
  float ele2R9;
  float ele2EnergyRegressionV0;
  float ele2EnergyRegressionV1;
  float ele2EnergyRegressionV2;
  Bool_t ele2PassHZZICHEP2012; 

  // for debugging 
  float eventNumber;

  zeeTree->SetBranchAddress("run", &run);
  zeeTree->SetBranchAddress("mass", &zMass);
  zeeTree->SetBranchAddress("weight", &weight);

  zeeTree->SetBranchAddress("Ele1Eta", &ele1Eta);
  zeeTree->SetBranchAddress("Ele1Phi", &ele1Phi);
  zeeTree->SetBranchAddress("Ele1Pt", &ele1Pt);
  zeeTree->SetBranchAddress("Ele1Energy", &ele1Energy);
  zeeTree->SetBranchAddress("Ele1GenPt", &ele1GenPt);
  zeeTree->SetBranchAddress("Ele1R9", &ele1R9);
  zeeTree->SetBranchAddress("Ele1EnergyRegressionV0", &ele1EnergyRegressionV0);
  zeeTree->SetBranchAddress("Ele1EnergyRegressionV1", &ele1EnergyRegressionV1);
  zeeTree->SetBranchAddress("Ele1EnergyRegressionV2", &ele1EnergyRegressionV2);
  zeeTree->SetBranchAddress("Ele1PassHZZICHEP2012", &ele1PassHZZICHEP2012);

  zeeTree->SetBranchAddress("Ele2Eta", &ele2Eta);
  zeeTree->SetBranchAddress("Ele2Phi", &ele2Phi);
  zeeTree->SetBranchAddress("Ele2Pt", &ele2Pt);
  zeeTree->SetBranchAddress("Ele2Energy", &ele2Energy);
  zeeTree->SetBranchAddress("Ele2GenPt", &ele2GenPt);
  zeeTree->SetBranchAddress("Ele2R9", &ele2R9);
  zeeTree->SetBranchAddress("Ele2EnergyRegressionV0", &ele2EnergyRegressionV0);
  zeeTree->SetBranchAddress("Ele2EnergyRegressionV1", &ele2EnergyRegressionV1);
  zeeTree->SetBranchAddress("Ele2EnergyRegressionV2", &ele2EnergyRegressionV2);
  zeeTree->SetBranchAddress("Ele2PassHZZICHEP2012", &ele2PassHZZICHEP2012);

  zeeTree->SetBranchAddress("event", &eventNumber);

  // Setting up the histgrams
  TH1F* zMass_hist = new TH1F("zMass_hist", "Z mass regression comparison ; Z mass ; Number of Events", 200, 0., 200.);
  TH1F* zMassRegressionV0_hist = new TH1F("zMassRegressionV0_hist", "Z mass regression comparison ; Z mass ; Number of Events", 200, 0., 200.);
  TH1F* zMassRegressionV1_hist = new TH1F("zMassRegressionV1_hist", "Z mass regression comparison ; Z mass ; Number of Events", 200, 0., 200.);
  TH1F* zMassRegressionV2_hist = new TH1F("zMassRegressionV2_hist", "Z mass regression comparison ; Z mass ; Number of Events", 200, 0., 200.);

  TH1F* zMass_hist_finebins = new TH1F("zMass_hist_finebins", "Z mass regression comparison ; Z mass ; Number of Events", 50000, 0., 200.);
  TH1F* zMassRegressionV0_hist_finebins = new TH1F("zMassRegressionV0_hist_finebins", "Z mass regression comparison ; Z mass ; Number of Events", 50000, 0., 200.);
  TH1F* zMassRegressionV1_hist_finebins = new TH1F("zMassRegressionV1_hist_finebins", "Z mass regression comparison ; Z mass ; Number of Events", 50000, 0., 200.);
  TH1F* zMassRegressionV2_hist_finebins = new TH1F("zMassRegressionV2_hist_finebins", "Z mass regression comparison ; Z mass ; Number of Events", 50000, 0., 200.);

  // For comparing regression in electron Pt
  TH1F* ele1Pt_hist = new TH1F("ele1Pt_hist", "Electron 1 Pt comparison ; Electron Pt / Generated electron Pt; Number of Events", 100, 0, 2);
  TH1F* ele1PtRegressionV0_hist = new TH1F("ele1PtRegressionV0_hist", "Electron 1 Pt comparison ; Electron Pt / Generated electron Pt; Number of Events", 100, 0, 2);
  TH1F* ele1PtRegressionV1_hist = new TH1F("ele1PtRegressionV1_hist", "Electron 1 Pt comparison ; Electron Pt / Generated electron Pt; Number of Events", 100, 0, 2);
  TH1F* ele1PtRegressionV2_hist = new TH1F("ele1PtRegressionV2_hist", "Electron 1 Pt comparison ; Electron Pt / Generated electron Pt; Number of Events", 100, 0, 2);

  TH1F* ele2Pt_hist = new TH1F("ele2Pt_hist", "Electron 2 Pt comparison ; Electron Pt / Generated electron Pt; Number of Events", 100, 0, 2);
  TH1F* ele2PtRegressionV0_hist = new TH1F("ele2PtRegressionV0_hist", "Electron 2 Pt comparison ; Electron Pt / Generated electron Pt; Number of Events", 100, 0, 2);
  TH1F* ele2PtRegressionV1_hist = new TH1F("ele2PtRegressionV1_hist", "Electron 2 Pt comparison ; Electron Pt / Generated electron Pt; Number of Events", 100, 0, 2);
  TH1F* ele2PtRegressionV2_hist = new TH1F("ele2PtRegressionV2_hist", "Electron 2 Pt comparison ; Electron Pt / Generated electron Pt; Number of Events", 100, 0, 2);

  int nEntries = zeeTree->GetEntries();
  // Setting up TLorentzVectors to hold electron data
  TLorentzVector* ele1Std = new TLorentzVector();
  TLorentzVector* ele1V0 = new TLorentzVector();
  TLorentzVector* ele1V1 = new TLorentzVector();
  TLorentzVector* ele1V2 = new TLorentzVector();

  TLorentzVector* ele2Std = new TLorentzVector();
  TLorentzVector* ele2V0 = new TLorentzVector();
  TLorentzVector* ele2V1 = new TLorentzVector();
  TLorentzVector* ele2V2 = new TLorentzVector();

  cout << "Total Events: " << nEntries << endl;
  for (int i = 0; i < nEntries; i++) {

    if (i%100000 == 0) cout << "Event: " << i << endl;
    zeeTree->GetEntry(i);

    // Cutting out events where electrons have too low pt (pt < 7)
    if (fabs(ele1Pt) < 20 || fabs(ele2Pt) < 20) continue;
    if (!(ele1PassHZZICHEP2012 && ele2PassHZZICHEP2012)) continue;


    //Select kinematic bin
    if (bin == 0) {
      if (!(fabs(ele1Eta) < 1.479 && fabs(ele2Eta) < 1.479 )) continue; 
    }
    if (bin == 1) {
      if (!( (fabs(ele1Eta) >= 1.479 && fabs(ele2Eta) < 1.479 )
            || 
             (fabs(ele1Eta) < 1.479 && fabs(ele2Eta) >= 1.479 )
            )
        ) continue; 
    }
    if (bin == 2) {
      if (!(fabs(ele1Eta) >= 1.479 && fabs(ele2Eta) >= 1.479 )) continue; 
    }

    double Ele1EnergyCorr = 0;
    double Ele1RegressionV0EnergyCorr = 0;
    double Ele1RegressionV1EnergyCorr = 0;
    double Ele1RegressionV2EnergyCorr = 0;
    double Ele2EnergyCorr = 0;
    double Ele2RegressionV0EnergyCorr = 0;
    double Ele2RegressionV1EnergyCorr = 0;
    double Ele2RegressionV2EnergyCorr = 0;


    if (is7TeV) {
      Ele1EnergyCorr = correctedElectronEnergy( ele1Energy, ele1Eta, ele1R9, run, 0, "2011", isMC, random );
      Ele1RegressionV0EnergyCorr = correctedElectronEnergy( ele1EnergyRegressionV0, ele1Eta, ele1R9, run, 0, "2011", isMC, random );
      Ele1RegressionV1EnergyCorr = correctedElectronEnergy( ele1EnergyRegressionV1, ele1Eta, ele1R9, run, 0, "2011", isMC, random );
      Ele1RegressionV2EnergyCorr = correctedElectronEnergy( ele1EnergyRegressionV2, ele1Eta, ele1R9, run, 0, "2011", isMC, random );
      Ele2EnergyCorr = correctedElectronEnergy( ele2Energy, ele2Eta, ele2R9, run, 0, "2011", isMC, random );
      Ele2RegressionV0EnergyCorr = correctedElectronEnergy( ele2EnergyRegressionV0, ele2Eta, ele2R9, run, 0, "2011", isMC, random );
      Ele2RegressionV1EnergyCorr = correctedElectronEnergy( ele2EnergyRegressionV1, ele2Eta, ele2R9, run, 0, "2011", isMC, random );
      Ele2RegressionV2EnergyCorr = correctedElectronEnergy( ele2EnergyRegressionV2, ele2Eta, ele2R9, run, 0, "2011", isMC, random );
    } else {
      Ele1EnergyCorr = correctedElectronEnergy( ele1Energy, ele1Eta, ele1R9, run, 0, "HCP2012", isMC, random );
      Ele1RegressionV0EnergyCorr = correctedElectronEnergy( ele1EnergyRegressionV0, ele1Eta, ele1R9, run, 1, "HCP2012", isMC, random );
      Ele1RegressionV1EnergyCorr = correctedElectronEnergy( ele1EnergyRegressionV1, ele1Eta, ele1R9, run, 1, "HCP2012", isMC, random );
      Ele1RegressionV2EnergyCorr = correctedElectronEnergy( ele1EnergyRegressionV2, ele1Eta, ele1R9, run, 1, "HCP2012", isMC, random );
      Ele2EnergyCorr = correctedElectronEnergy( ele2Energy, ele2Eta, ele2R9, run, 0, "HCP2012", isMC, random );
      Ele2RegressionV0EnergyCorr = correctedElectronEnergy( ele2EnergyRegressionV0, ele2Eta, ele2R9, run, 1, "HCP2012", isMC, random );
      Ele2RegressionV1EnergyCorr = correctedElectronEnergy( ele2EnergyRegressionV1, ele2Eta, ele2R9, run, 1, "HCP2012", isMC, random );
      Ele2RegressionV2EnergyCorr = correctedElectronEnergy( ele2EnergyRegressionV2, ele2Eta, ele2R9, run, 1, "HCP2012", isMC, random );
    }

//     cout << ele1Energy << " " << Ele1EnergyCorr << " : " << ele1EnergyRegressionV0 << " " 
//          << Ele1RegressionV0EnergyCorr << " " 
//          << Ele1RegressionV1EnergyCorr << " " 
//          << Ele1RegressionV2EnergyCorr << " " 
//          <<  endl;

    // Now for each regression fill in the electron using TLorentzVectors
    if (doCorr) {
      ele1Std->SetPtEtaPhiM( Ele1EnergyCorr/ TMath::CosH(ele1Eta),  ele1Eta, ele1Phi, ELECTRONMASS);
      ele1V0->SetPtEtaPhiM( Ele1RegressionV0EnergyCorr / TMath::CosH(ele1Eta),  ele1Eta, ele1Phi, ELECTRONMASS);
      ele1V1->SetPtEtaPhiM( Ele1RegressionV1EnergyCorr / TMath::CosH(ele1Eta), ele1Eta, ele1Phi, ELECTRONMASS);
      ele1V2->SetPtEtaPhiM( Ele1RegressionV2EnergyCorr / TMath::CosH(ele1Eta), ele1Eta, ele1Phi, ELECTRONMASS);
      ele2Std->SetPtEtaPhiM( Ele2EnergyCorr / TMath::CosH(ele2Eta), ele2Eta, ele2Phi, ELECTRONMASS);
      ele2V0->SetPtEtaPhiM( Ele2RegressionV0EnergyCorr / TMath::CosH(ele2Eta), ele2Eta, ele2Phi, ELECTRONMASS);
      ele2V1->SetPtEtaPhiM( Ele2RegressionV1EnergyCorr / TMath::CosH(ele2Eta), ele2Eta, ele2Phi, ELECTRONMASS);
      ele2V2->SetPtEtaPhiM( Ele2RegressionV2EnergyCorr / TMath::CosH(ele2Eta), ele2Eta, ele2Phi, ELECTRONMASS);

    } else {
      ele1Std->SetPtEtaPhiM( ele1Energy/ TMath::CosH(ele1Eta),  ele1Eta, ele1Phi, ELECTRONMASS);
      ele1V0->SetPtEtaPhiM( ele1EnergyRegressionV0/ TMath::CosH(ele1Eta),  ele1Eta, ele1Phi, ELECTRONMASS);
      ele1V1->SetPtEtaPhiM( ele1EnergyRegressionV1 / TMath::CosH(ele1Eta), ele1Eta, ele1Phi, ELECTRONMASS);
      ele1V2->SetPtEtaPhiM( ele1EnergyRegressionV2 / TMath::CosH(ele1Eta), ele1Eta, ele1Phi, ELECTRONMASS);
      ele2Std->SetPtEtaPhiM( ele2Energy / TMath::CosH(ele2Eta), ele2Eta, ele2Phi, ELECTRONMASS);
      ele2V0->SetPtEtaPhiM( ele2EnergyRegressionV0 / TMath::CosH(ele2Eta), ele2Eta, ele2Phi, ELECTRONMASS);
      ele2V1->SetPtEtaPhiM( ele2EnergyRegressionV1 / TMath::CosH(ele2Eta), ele2Eta, ele2Phi, ELECTRONMASS);
      ele2V2->SetPtEtaPhiM( ele2EnergyRegressionV2 / TMath::CosH(ele2Eta), ele2Eta, ele2Phi, ELECTRONMASS);

    }

    // Now filling up the histograms for Z mass
    zMass_hist->Fill( ((*ele1Std) + (*ele2Std)).M(), weight);
    zMassRegressionV0_hist->Fill( ((*ele1V0) + (*ele2V0)).M(), weight);
    zMassRegressionV1_hist->Fill( ((*ele1V1) + (*ele2V1)).M(), weight);
    zMassRegressionV2_hist->Fill( ((*ele1V2) + (*ele2V2)).M(), weight);
    zMass_hist_finebins->Fill( ((*ele1Std) + (*ele2Std)).M(), weight);
    zMassRegressionV0_hist_finebins->Fill( ((*ele1V0) + (*ele2V0)).M(), weight);
    zMassRegressionV1_hist_finebins->Fill( ((*ele1V1) + (*ele2V1)).M(), weight);
    zMassRegressionV2_hist_finebins->Fill( ((*ele1V2) + (*ele2V2)).M(), weight);

  }

  //normalize histograms
  NormalizeHist(zMass_hist);
  NormalizeHist(zMassRegressionV0_hist);
  NormalizeHist(zMassRegressionV1_hist);
  NormalizeHist(zMassRegressionV2_hist);

  // Now plotting for comparison
  // First, finding maximum Y range
  float maxYRange = zMass_hist->GetMaximum();
  maxYRange = max(maxYRange, (float) zMassRegressionV0_hist->GetMaximum());
  maxYRange = max(maxYRange, (float) zMassRegressionV1_hist->GetMaximum());
  maxYRange = max(maxYRange, (float) zMassRegressionV2_hist->GetMaximum());

  // Graphic settings
  zMass_hist->SetMaximum(maxYRange * 1.25);
  zMass_hist->SetLineColor(kBlack);
  zMass_hist->SetLineWidth(3);
  zMassRegressionV0_hist->SetLineColor(kRed);
  zMassRegressionV0_hist->SetLineWidth(2);
  zMassRegressionV1_hist->SetLineColor(kGreen);
  zMassRegressionV1_hist->SetLineWidth(2);
  zMassRegressionV2_hist->SetLineColor(kBlue);
  zMassRegressionV2_hist->SetLineWidth(2);


  TCanvas* cv = new TCanvas("cv", "cv", 800, 600); 
  TLegend* legend = new TLegend(0.25,0.7,0.5,0.9); 
  legend->SetTextSize(0.03); 
  legend->SetBorderSize(0); 
  legend->SetFillStyle(0); 
  legend->AddEntry(zMass_hist,"No Regression", "L");
  legend->AddEntry(zMassRegressionV0_hist, "Regression", 	 "L");
  legend->AddEntry(zMassRegressionV1_hist, "Regression w/trk vars",	"L");
  legend->AddEntry(zMassRegressionV2_hist, "Regression w/more trk vars",	"L");

  zMass_hist->SetTitle("");
  zMass_hist->GetYaxis()->SetTitle("Fraction of Events");
  zMass_hist->GetXaxis()->SetRangeUser(60,120);
  zMass_hist->Draw("hist");
  zMassRegressionV0_hist->Draw("histsame");
  zMassRegressionV1_hist->Draw("histsame");
  zMassRegressionV2_hist->Draw("histsame");
  legend->Draw();
  cv->SaveAs(Form("ZeeMass_%s_bin%i.gif",label.c_str(),bin));


  
  TFile *file = new TFile("ZeeMass.root", "UPDATE");
  file->WriteTObject(zMass_hist, Form("ZeeMassStandard_%s_bin%i",label.c_str(),bin), "WriteDelete");
  file->WriteTObject(zMassRegressionV0_hist, Form("ZeeMassRegressionV0_%s_bin%i",label.c_str(),bin), "WriteDelete");
  file->WriteTObject(zMassRegressionV1_hist, Form("ZeeMassRegressionV1_%s_bin%i",label.c_str(),bin), "WriteDelete");
  file->WriteTObject(zMassRegressionV2_hist, Form("ZeeMassRegressionV2_%s_bin%i",label.c_str(),bin), "WriteDelete");
  file->WriteTObject(zMass_hist_finebins, Form("ZeeMassStandard_finebins_%s_bin%i",label.c_str(),bin), "WriteDelete");
  file->WriteTObject(zMassRegressionV0_hist_finebins, Form("ZeeMassRegressionV0_finebins_%s_bin%i",label.c_str(),bin), "WriteDelete");
  file->WriteTObject(zMassRegressionV1_hist_finebins, Form("ZeeMassRegressionV1_finebins_%s_bin%i",label.c_str(),bin), "WriteDelete");
  file->WriteTObject(zMassRegressionV2_hist_finebins, Form("ZeeMassRegressionV2_finebins_%s_bin%i",label.c_str(),bin), "WriteDelete");
  file->Close();
  delete file;


}

void compareRegressionZMasses ( string label, Int_t bin) {

  TFile *file = new TFile("ZeeMass.root", "UPDATE");

  TH1F *massHistStd = 0;
  TH1F *massHistRegressionV0 = 0;
  TH1F *massHistRegressionV1 = 0;
  TH1F *massHistRegressionV2 = 0;
  TH1F *massHistStd_finebins = 0;
  TH1F *massHistRegressionV0_finebins = 0;
  TH1F *massHistRegressionV1_finebins = 0;
  TH1F *massHistRegressionV2_finebins = 0;

  massHistStd = (TH1F*)file->Get(Form("ZeeMassStandard_%s_bin%i",label.c_str(),bin));
  massHistRegressionV0 = (TH1F*)file->Get(Form("ZeeMassRegressionV0_%s_bin%i",label.c_str(),bin));
  massHistRegressionV1 = (TH1F*)file->Get(Form("ZeeMassRegressionV1_%s_bin%i",label.c_str(),bin));
  massHistRegressionV2 = (TH1F*)file->Get(Form("ZeeMassRegressionV2_%s_bin%i",label.c_str(),bin));
  massHistStd_finebins = (TH1F*)file->Get(Form("ZeeMassStandard_finebins_%s_bin%i",label.c_str(),bin));
  massHistRegressionV0_finebins = (TH1F*)file->Get(Form("ZeeMassRegressionV0_finebins_%s_bin%i",label.c_str(),bin));
  massHistRegressionV1_finebins = (TH1F*)file->Get(Form("ZeeMassRegressionV1_finebins_%s_bin%i",label.c_str(),bin));
  massHistRegressionV2_finebins = (TH1F*)file->Get(Form("ZeeMassRegressionV2_finebins_%s_bin%i",label.c_str(),bin));


  assert( massHistStd );
  assert( massHistRegressionV0 );
  assert( massHistRegressionV1 );
  assert( massHistRegressionV2 );

  //Create Data Set
  RooRealVar mass("mass","m(EE)",60,120,"GeV/c^{2}");
  
  RooArgSet zMassArgSet(mass);

  //********************************************************************************
  //Setup
  //********************************************************************************
  double minMass = 70;
  double maxMass = 110;
  double mean_bw = 91.1876;
  double gamma_bw = 2.4952;
  double cutoff_cb = 1.0;
  double power_cb = 2.45;
  const char *plotOpt = "NEU";
  const int nbins = 40;

  //********************************************************************************
  //Data
  //********************************************************************************
  cout << "here01\n";
  RooDataHist* dataStd = new RooDataHist("dataStd", "zee mass", zMassArgSet,Import(*massHistStd, kFALSE) );
  RooDataHist* dataRegressionV0 = new RooDataHist("dataRegressionV0", "zee mass", zMassArgSet,Import(*massHistRegressionV0, kFALSE) );
  RooDataHist* dataRegressionV1 = new RooDataHist("dataRegressionV1", "zee mass", zMassArgSet,Import(*massHistRegressionV1, kFALSE) );
  RooDataHist* dataRegressionV2 = new RooDataHist("dataRegressionV2", "zee mass", zMassArgSet,Import(*massHistRegressionV2, kFALSE) );
  cout << "here02\n";
  if (massHistStd) {cout << "done\n";}

  //********************************************************************************
  //Fit For No Regression
  //********************************************************************************
  RooRealVar cbBias_std ("#Deltam_{CB}_std", "CB Bias", -.01, -10, 10, "GeV/c^{2}");
  RooRealVar cbSigma_std("sigma_{CB}_std", "CB Width", 1.5, 0.8, 5.0, "GeV/c^{2}");
  RooRealVar cbCut_std  ("a_{CB}_std","CB Cut", 1.0, 1.0, 3.0);
  RooRealVar cbPower_std("n_{CB}_std","CB Order", 2.5, 0.1, 20.0);
  cbCut_std.setVal(cutoff_cb);
  cbPower_std.setVal(power_cb);

  //Breit_Wigner parameters
  RooRealVar bwMean_std("m_{Z}","BW Mean", 91.1876, "GeV/c^{2}");
  bwMean_std.setVal(mean_bw);
  RooRealVar bwWidth_std("#Gamma_{Z}", "BW Width", 2.4952, "GeV/c^{2}");
  bwWidth_std.setVal(gamma_bw);

  // Fix the Breit-Wigner parameters to PDG values
  bwMean_std.setConstant(kTRUE);
  bwWidth_std.setConstant(kTRUE);

  // Exponential Background parameters
  RooRealVar expRate_std("#lambda_{exp}_std", "Exponential Rate", -0.064, -1, 1);
  RooRealVar c0_std("c_{0}_std", "c0", 1., 0., 50.);

  //Number of Signal and Background events
  RooRealVar nsig_std("N_{S}_std", "# signal events", 524, 0.1, 10000000000.);
  RooRealVar nbkg_std("N_{B}_std", "# background events", 43, 1., 10000000.);

  //============================ P.D.F.s=============================
  // Mass signal for two decay electrons p.d.f.
  RooBreitWigner bw_std("bw_std", "bw_std", mass, bwMean_std, bwWidth_std);
  RooCBShape  cball_std("cball_std", "Crystal Ball", mass, cbBias_std, cbSigma_std, cbCut_std, cbPower_std);
  RooFFTConvPdf BWxCB_std("BWxCB_std", "bw X crystal ball", mass, bw_std, cball_std);

  // Mass model for signal electrons p.d.f.
  RooAddPdf model_std("model_std", "signal", RooArgList(BWxCB_std), RooArgList(nsig_std));

  mass.setRange(60,120);
  model_std.fitTo(*dataStd,FitOptions("mh"),Optimize(0),Timer(1));




  //********************************************************************************
  //Fit For RegressionV0
  //********************************************************************************
  RooRealVar cbBias_regressionV0 ("#Deltam_{CB}_regressionV0", "CB Bias", -.01, -10, 10, "GeV/c^{2}");
  RooRealVar cbSigma_regressionV0("sigma_{CB}_regressionV0", "CB Width", 1.5, 0.8, 5.0, "GeV/c^{2}");
  RooRealVar cbCut_regressionV0  ("a_{CB}_regressionV0","CB Cut", 1.0, 1.0, 3.0);
  RooRealVar cbPower_regressionV0("n_{CB}_regressionV0","CB Order", 2.5, 0.1, 20.0);
  cbCut_regressionV0.setVal(cutoff_cb);
  cbPower_regressionV0.setVal(power_cb);

  //Breit_Wigner parameters
  RooRealVar bwMean_regressionV0("m_{Z}","BW Mean", 91.1876, "GeV/c^{2}");
  bwMean_regressionV0.setVal(mean_bw);
  RooRealVar bwWidth_regressionV0("#Gamma_{Z}", "BW Width", 2.4952, "GeV/c^{2}");
  bwWidth_regressionV0.setVal(gamma_bw);

  // Fix the Breit-Wigner parameters to PDG values
  bwMean_regressionV0.setConstant(kTRUE);
  bwWidth_regressionV0.setConstant(kTRUE);

  // Exponential Background parameters
  RooRealVar expRate_regressionV0("#lambda_{exp}_regressionV0", "Exponential Rate", -0.064, -1, 1);
  RooRealVar c0_regressionV0("c_{0}_regressionV0", "c0", 1., 0., 50.);

  //Number of Signal and Background events
  RooRealVar nsig_regressionV0("N_{S}_regressionV0", "# signal events", 524, 0.1, 10000000000.);
  RooRealVar nbkg_regressionV0("N_{B}_regressionV0", "# background events", 43, 1., 10000000.);

  //============================ P.D.F.s=============================
  // Mass signal for two decay electrons p.d.f.
  RooBreitWigner bw_regressionV0("bw_regressionV0", "bw_regressionV0", mass, bwMean_regressionV0, bwWidth_regressionV0);
  RooCBShape  cball_regressionV0("cball_regressionV0", "Crystal Ball", mass, cbBias_regressionV0, cbSigma_regressionV0, cbCut_regressionV0, cbPower_regressionV0);
  RooFFTConvPdf BWxCB_regressionV0("BWxCB_regressionV0", "bw X crystal ball", mass, bw_regressionV0, cball_regressionV0);

  // Mass model for signal electrons p.d.f.
  RooAddPdf model_regressionV0("model_regressionV0", "signal", RooArgList(BWxCB_regressionV0), RooArgList(nsig_regressionV0));

  mass.setRange(60,120);
  model_regressionV0.fitTo(*dataRegressionV0,FitOptions("mh"),Optimize(0),Timer(1));



  //********************************************************************************
  //Fit For RegressionV1
  //********************************************************************************
  RooRealVar cbBias_regressionV1 ("#Deltam_{CB}_regressionV1", "CB Bias", -.01, -10, 10, "GeV/c^{2}");
  RooRealVar cbSigma_regressionV1("sigma_{CB}_regressionV1", "CB Width", 1.5, 0.8, 5.0, "GeV/c^{2}");
  RooRealVar cbCut_regressionV1  ("a_{CB}_regressionV1","CB Cut", 1.0, 1.0, 3.0);
  RooRealVar cbPower_regressionV1("n_{CB}_regressionV1","CB Order", 2.5, 0.1, 20.0);
  cbCut_regressionV1.setVal(cutoff_cb);
  cbPower_regressionV1.setVal(power_cb);

  //Breit_Wigner parameters
  RooRealVar bwMean_regressionV1("m_{Z}","BW Mean", 91.1876, "GeV/c^{2}");
  bwMean_regressionV1.setVal(mean_bw);
  RooRealVar bwWidth_regressionV1("#Gamma_{Z}", "BW Width", 2.4952, "GeV/c^{2}");
  bwWidth_regressionV1.setVal(gamma_bw);

  // Fix the Breit-Wigner parameters to PDG values
  bwMean_regressionV1.setConstant(kTRUE);
  bwWidth_regressionV1.setConstant(kTRUE);

  // Exponential Background parameters
  RooRealVar expRate_regressionV1("#lambda_{exp}_regressionV1", "Exponential Rate", -0.064, -1, 1);
  RooRealVar c0_regressionV1("c_{0}_regressionV1", "c0", 1., 0., 50.);

  //Number of Signal and Background events
  RooRealVar nsig_regressionV1("N_{S}_regressionV1", "# signal events", 524, 0.1, 10000000000.);
  RooRealVar nbkg_regressionV1("N_{B}_regressionV1", "# background events", 43, 1., 10000000.);

  //============================ P.D.F.s=============================
  // Mass signal for two decay electrons p.d.f.
  RooBreitWigner bw_regressionV1("bw_regressionV1", "bw_regressionV1", mass, bwMean_regressionV1, bwWidth_regressionV1);
  RooCBShape  cball_regressionV1("cball_regressionV1", "Crystal Ball", mass, cbBias_regressionV1, cbSigma_regressionV1, cbCut_regressionV1, cbPower_regressionV1);
  RooFFTConvPdf BWxCB_regressionV1("BWxCB_regressionV1", "bw X crystal ball", mass, bw_regressionV1, cball_regressionV1);

  // Mass model for signal electrons p.d.f.
  RooAddPdf model_regressionV1("model_regressionV1", "signal", RooArgList(BWxCB_regressionV1), RooArgList(nsig_regressionV1));

  mass.setRange(60,120);
  model_regressionV1.fitTo(*dataRegressionV1,FitOptions("mh"),Optimize(0),Timer(1));


  //********************************************************************************
  //Fit For RegressionV2
  //********************************************************************************
  RooRealVar cbBias_regressionV2 ("#Deltam_{CB}_regressionV2", "CB Bias", -.01, -10, 10, "GeV/c^{2}");
  RooRealVar cbSigma_regressionV2("sigma_{CB}_regressionV2", "CB Width", 1.5, 0.8, 5.0, "GeV/c^{2}");
  RooRealVar cbCut_regressionV2  ("a_{CB}_regressionV2","CB Cut", 1.0, 1.0, 3.0);
  RooRealVar cbPower_regressionV2("n_{CB}_regressionV2","CB Order", 2.5, 0.1, 20.0);
  cbCut_regressionV2.setVal(cutoff_cb);
  cbPower_regressionV2.setVal(power_cb);

  //Breit_Wigner parameters
  RooRealVar bwMean_regressionV2("m_{Z}","BW Mean", 91.1876, "GeV/c^{2}");
  bwMean_regressionV2.setVal(mean_bw);
  RooRealVar bwWidth_regressionV2("#Gamma_{Z}", "BW Width", 2.4952, "GeV/c^{2}");
  bwWidth_regressionV2.setVal(gamma_bw);

  // Fix the Breit-Wigner parameters to PDG values
  bwMean_regressionV2.setConstant(kTRUE);
  bwWidth_regressionV2.setConstant(kTRUE);

  // Exponential Background parameters
  RooRealVar expRate_regressionV2("#lambda_{exp}_regressionV2", "Exponential Rate", -0.064, -1, 1);
  RooRealVar c0_regressionV2("c_{0}_regressionV2", "c0", 1., 0., 50.);

  //Number of Signal and Background events
  RooRealVar nsig_regressionV2("N_{S}_regressionV2", "# signal events", 524, 0.1, 10000000000.);
  RooRealVar nbkg_regressionV2("N_{B}_regressionV2", "# background events", 43, 1., 10000000.);

  //============================ P.D.F.s=============================
  // Mass signal for two decay electrons p.d.f.
  RooBreitWigner bw_regressionV2("bw_regressionV2", "bw_regressionV2", mass, bwMean_regressionV2, bwWidth_regressionV2);
  RooCBShape  cball_regressionV2("cball_regressionV2", "Crystal Ball", mass, cbBias_regressionV2, cbSigma_regressionV2, cbCut_regressionV2, cbPower_regressionV2);
  RooFFTConvPdf BWxCB_regressionV2("BWxCB_regressionV2", "bw X crystal ball", mass, bw_regressionV2, cball_regressionV2);

  // Mass model for signal electrons p.d.f.
  RooAddPdf model_regressionV2("model_regressionV2", "signal", RooArgList(BWxCB_regressionV2), RooArgList(nsig_regressionV2));

  mass.setRange(60,120);
  model_regressionV2.fitTo(*dataRegressionV2,FitOptions("mh"),Optimize(0),Timer(1));




  //********************************************************************************
  //Compute Effective Sigma
  //********************************************************************************
  float sigmaeff_std = effSigma(massHistStd_finebins );  
  float sigmaeff_regressionV0 = effSigma(massHistRegressionV0_finebins );  
  float sigmaeff_regressionV1 = effSigma(massHistRegressionV1_finebins );  
  float sigmaeff_regressionV2 = effSigma(massHistRegressionV2_finebins );  

//   RooAbsReal *cdf = model_std.createCdf(mass);
//   float testmass = 91.1876; 
//   float center = testmass-10.0;
//   float minwidth = 999.0;
//   float mlmin = 0.0;
//   float mhmin = 0.0;
//   float step=0.01;
//   int Nstep = int(15/step+0.1);
  
//   int kkk = 0; 
//   for (int i=0; i<Nstep; ++i) {
//     cout << i << endl;
//     float mlow = center+i*step;
//     mass.setVal(mlow);
//     float cdflo = cdf->getVal();
//     for (int j=i+1; j<Nstep; ++j) {
//       float mhigh = center+j*step;
//       mass.setVal(mhigh);
//       float cdfhi = cdf->getVal();
//       if ( (cdfhi-cdflo)>0.683 ) {
//         if ( (mhigh-mlow)<minwidth) {
//           minwidth = mhigh-mlow;
//           mlmin = mlow;
//           mhmin = mhigh;
//         }
//         break;
//       }
//     }
//   }
//   sigmaeff_std = minwidth/2.0;
  
//   cdf = model_regressionV0.createCdf(mass);
//   testmass = 91.1876; 
//   center = testmass-10.0;
//   minwidth = 999.0;
//   mlmin = 0.0;
//   mhmin = 0.0;
//   step=0.01;
//   Nstep = int(15/step+0.1);
  
//   kkk = 0; 
//   for (int i=0; i<Nstep; ++i) {
//     cout << i << endl;
//     float mlow = center+i*step;
//     mass.setVal(mlow);
//     float cdflo = cdf->getVal();
//     for (int j=i+1; j<Nstep; ++j) {
//       float mhigh = center+j*step;
//       mass.setVal(mhigh);
//       float cdfhi = cdf->getVal();
//       if ( (cdfhi-cdflo)>0.683 ) {
//         if ( (mhigh-mlow)<minwidth) {
//           minwidth = mhigh-mlow;
//           mlmin = mlow;
//           mhmin = mhigh;
//         }
//         break;
//       }
//     }
//   }
//   sigmaeff_regressionV0 = minwidth/2.0;

 
  cout << sigmaeff_std << endl;
  cout << sigmaeff_regressionV0 << endl;
  cout << sigmaeff_regressionV1 << endl;
  cout << sigmaeff_regressionV2 << endl;


  //********************************************************************************
  //Plot
  //********************************************************************************


  TCanvas* c = new TCanvas("c","Unbinned Invariant Mass Fit", 0,0,800,600);
  gStyle->SetOptTitle(0);
  //========================== Plotting  ============================
  //Create a frame
  RooPlot* plot = mass.frame(Range(75,105),Bins(30));
  // Add data and model to canvas
  dataStd->plotOn(plot, MarkerColor(kBlack));
  model_std.plotOn(plot, LineColor(kBlack));
  dataRegressionV0->plotOn(plot, MarkerColor(kRed));
  model_regressionV0.plotOn(plot, LineColor(kRed));
  dataRegressionV1->plotOn(plot, MarkerColor(kGreen+2));
  model_regressionV1.plotOn(plot, LineColor(kGreen+2));
  dataRegressionV2->plotOn(plot, MarkerColor(kBlue));
  model_regressionV2.plotOn(plot, LineColor(kBlue));
  plot->Draw();

  // Print Fit Values
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(.04);
  tex->SetTextFont(2);
  tex->Draw();
//   tex->SetTextSize(0.022);
  tex->SetTextSize(0.024);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.745, 0.64, "No Regression");
  tex->DrawLatex(0.745, 0.59, Form("CB Bias = %.2f GeV/c^{2}", cbBias_std.getVal()));
//   tex->DrawLatex(0.745, 0.54, Form("#sigma = %.2f GeV/c^{2}", cbSigma_std.getVal()));
  tex->DrawLatex(0.745, 0.54, Form("#sigma Effective = %.2f GeV/c^{2}", sigmaeff_std));

  tex->SetTextSize(0.024);
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.745, 0.84, "RegressionV0");
  tex->DrawLatex(0.745, 0.79, Form("CB Bias = %.2f GeV/c^{2}", cbBias_regressionV0.getVal()));
//   tex->DrawLatex(0.745, 0.74, Form("#sigma = %.2f GeV/c^{2}", cbSigma_regressionV0.getVal()));
   tex->DrawLatex(0.745, 0.74, Form("#sigma Effective = %.2f GeV/c^{2}", sigmaeff_regressionV0));


  tex->SetTextSize(0.024);
  tex->SetTextColor(kGreen+2);
  tex->DrawLatex(0.245, 0.64, "RegressionV1");
  tex->DrawLatex(0.245, 0.59, Form("CB Bias = %.2f GeV/c^{2}", cbBias_regressionV1.getVal()));
//   tex->DrawLatex(0.245, 0.54, Form("#sigma = %.2f GeV/c^{2}", cbSigma_regressionV1.getVal()));
   tex->DrawLatex(0.245, 0.54, Form("#sigma Effective = %.2f GeV/c^{2}", sigmaeff_regressionV1));


  tex->SetTextSize(0.024);
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.245, 0.84, "RegressionV2");
  tex->DrawLatex(0.245, 0.79, Form("CB Bias = %.2f GeV/c^{2}", cbBias_regressionV2.getVal()));
//   tex->DrawLatex(0.245, 0.74, Form("#sigma = %.2f GeV/c^{2}", cbSigma_regressionV2.getVal()));
  tex->DrawLatex(0.245, 0.74, Form("#sigma Effective = %.2f GeV/c^{2}", sigmaeff_regressionV2));


  c->Update();
  c->SaveAs(Form("ZeeMass_%s_bin%i.gif",label.c_str(),bin));



  //********************************************************************************
  //Plot just the mass spectra in fine bins
  //********************************************************************************
  TCanvas* cv = new TCanvas("c","Unbinned Invariant Mass Fit", 0,0,800,600);
  gStyle->SetOptTitle(0);
  
  massHistStd_finebins->Rebin(40);
  massHistRegressionV0_finebins->Rebin(40);
  massHistRegressionV1_finebins->Rebin(40);
  massHistRegressionV2_finebins->Rebin(40);

  massHistStd_finebins->SetLineColor(kBlack);
  massHistRegressionV0_finebins->SetLineColor(kRed);
  massHistRegressionV1_finebins->SetLineColor(kGreen+2);
  massHistRegressionV2_finebins->SetLineColor(kBlue);

  massHistStd_finebins->GetXaxis()->SetRangeUser(80,100);
  massHistRegressionV0_finebins->GetXaxis()->SetRangeUser(80,100);
  massHistRegressionV1_finebins->GetXaxis()->SetRangeUser(80,100);
  massHistRegressionV2_finebins->GetXaxis()->SetRangeUser(80,100);

  double max = 1.2 * massHistRegressionV2_finebins->GetMaximum();

  massHistStd_finebins->SetMaximum(max);
  massHistStd_finebins->Draw("hist");
  massHistRegressionV0_finebins->Draw("hist,same");
  massHistRegressionV1_finebins->Draw("hist,same");
  massHistRegressionV2_finebins->Draw("hist,same");

  // Print Fit Values
  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(.04);
  tex->SetTextFont(2);
  tex->Draw();
//   tex->SetTextSize(0.022);
  tex->SetTextSize(0.024);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.245, 0.84, "No Regression");
  tex->DrawLatex(0.245, 0.79, Form("#sigma Effective = %.2f GeV/c^{2}", sigmaeff_std));

  tex->SetTextSize(0.024);
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.245, 0.74, "RegressionV0");
   tex->DrawLatex(0.245, 0.69, Form("#sigma Effective = %.2f GeV/c^{2}", sigmaeff_regressionV0));

  tex->SetTextSize(0.024);
  tex->SetTextColor(kGreen+2);
  tex->DrawLatex(0.245, 0.64, "RegressionV1");
   tex->DrawLatex(0.245, 0.59, Form("#sigma Effective = %.2f GeV/c^{2}", sigmaeff_regressionV1));


  tex->SetTextSize(0.024);
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.245, 0.54, "RegressionV2");
  tex->DrawLatex(0.245, 0.49, Form("#sigma Effective = %.2f GeV/c^{2}", sigmaeff_regressionV2));



  cv->SaveAs(Form("ZeeMassHistogram_%s_bin%i.gif",label.c_str(),bin));


}

