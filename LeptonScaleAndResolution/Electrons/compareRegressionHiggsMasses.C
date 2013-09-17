//***************************************************
//2012 Data
//root -l CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/compareRegressionHiggsMasses.C+'("hzzTree_id1125.noregression.root","hzzTree_id1125.regressionV0.root","hzzTree_id1125.regressionV1.root","hzzTree_id1125.regressionV2.root", "4e",1)'
//root -l CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/compareRegressionHiggsMasses.C+'("hzzTree_id1125.noregression.root","hzzTree_id1125.regressionV0.root","hzzTree_id1125.regressionV1.root","hzzTree_id1125.regressionV2.root", "2e2mu",2)'
//root -l CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/compareRegressionHiggsMasses.C+'("hzzTree_id1125.noregression.root","hzzTree_id1125.regressionV0.root","hzzTree_id1125.regressionV1.root","hzzTree_id1125.regressionV2.root", "2mu2e",3)'
//***************************************************

//***************************************************
//2011 Data
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



void compareRegressionHiggsMasses ( string filenameNoRegression,
                                    string filenameRegressionV0,  
                                    string filenameRegressionV1,  
                                    string filenameRegressionV2,
                                    string label,
                                    Int_t channel ) {

  gStyle->SetOptStat(0);

  TFile *fileNoRegression = new TFile(filenameNoRegression.c_str(), "READ");
  TFile *fileRegressionV0 = new TFile(filenameRegressionV0.c_str(), "READ");
  TFile *fileRegressionV1 = new TFile(filenameRegressionV1.c_str(), "READ");
  TFile *fileRegressionV2 = new TFile(filenameRegressionV2.c_str(), "READ");
  
  TTree *treeNoRegression = (TTree*)fileNoRegression->Get("zz4lTree/probe_tree");
  TTree *treeRegressionV0 = (TTree*)fileRegressionV0->Get("zz4lTree/probe_tree");
  TTree *treeRegressionV1 = (TTree*)fileRegressionV1->Get("zz4lTree/probe_tree");
  TTree *treeRegressionV2 = (TTree*)fileRegressionV2->Get("zz4lTree/probe_tree");

  assert(treeNoRegression);
  assert(treeRegressionV0);
  assert(treeRegressionV1);
  assert(treeRegressionV2);

  TH1F *massNoRegression = new TH1F( "massNoRegression", ";m_{4l} [GeV/c^{2}]; Fraction of Events", 160, 100, 150);
  TH1F *massRegressionV0 = new TH1F( "massRegressionV0", ";m_{4l} [GeV/c^{2}]; Fraction of Events", 160, 100, 150);
  TH1F *massRegressionV1 = new TH1F( "massRegressionV1", ";m_{4l} [GeV/c^{2}]; Fraction of Events", 160, 100, 150);
  TH1F *massRegressionV2 = new TH1F( "massRegressionV2", ";m_{4l} [GeV/c^{2}]; Fraction of Events", 160, 100, 150);
  TH1F *massNoRegression_finebins = new TH1F( "massNoRegression_finebins", ";m_{4l} [GeV/c^{2}]; Fraction of Events", 50000, 0, 200);
  TH1F *massRegressionV0_finebins = new TH1F( "massRegressionV0_finebins", ";m_{4l} [GeV/c^{2}]; Fraction of Events", 50000, 0, 200);
  TH1F *massRegressionV1_finebins = new TH1F( "massRegressionV1_finebins", ";m_{4l} [GeV/c^{2}]; Fraction of Events", 50000, 0, 200);
  TH1F *massRegressionV2_finebins = new TH1F( "massRegressionV2_finebins", ";m_{4l} [GeV/c^{2}]; Fraction of Events", 50000, 0, 200);

  //********************************************************************************
  //Fill Histogram
  //********************************************************************************
  treeNoRegression->Draw("mass>>massNoRegression_finebins", Form("channel==%d",channel));
  treeRegressionV0->Draw("mass>>massRegressionV0_finebins", Form("channel==%d",channel));
  treeRegressionV1->Draw("mass>>massRegressionV1_finebins", Form("channel==%d",channel));
  treeRegressionV2->Draw("mass>>massRegressionV2_finebins", Form("channel==%d",channel));

  treeNoRegression->Draw("mass>>massNoRegression", Form("channel==%d",channel));
  treeRegressionV0->Draw("mass>>massRegressionV0", Form("channel==%d",channel));
  treeRegressionV1->Draw("mass>>massRegressionV1", Form("channel==%d",channel));
  treeRegressionV2->Draw("mass>>massRegressionV2", Form("channel==%d",channel));

  NormalizeHist(massNoRegression);
  NormalizeHist(massRegressionV0);
  NormalizeHist(massRegressionV1);
  NormalizeHist(massRegressionV2);

  //********************************************************************************
  //Compute Effective Sigma
  //********************************************************************************
  float sigmaeff_std = effSigma(massNoRegression_finebins );  
  float sigmaeff_regressionV0 = effSigma(massRegressionV0_finebins );  
  float sigmaeff_regressionV1 = effSigma(massRegressionV1_finebins );  
  float sigmaeff_regressionV2 = effSigma(massRegressionV2_finebins );  
 
  cout << sigmaeff_std << endl;
  cout << sigmaeff_regressionV0 << endl;
  cout << sigmaeff_regressionV1 << endl;
  cout << sigmaeff_regressionV2 << endl;



  //********************************************************************************
  //Plot just the mass spectra in fine bins
  //********************************************************************************
  TCanvas* cv = new TCanvas("c","Unbinned Invariant Mass Fit", 0,0,800,600);
  gStyle->SetOptTitle(0);
  
  massNoRegression->SetLineColor(kBlack);
  massRegressionV0->SetLineColor(kRed);
  massRegressionV1->SetLineColor(kGreen+2);
  massRegressionV2->SetLineColor(kBlue);

  massNoRegression->GetXaxis()->SetRangeUser(110,140);
  massRegressionV0->GetXaxis()->SetRangeUser(110,140);
  massRegressionV1->GetXaxis()->SetRangeUser(110,140);
  massRegressionV2->GetXaxis()->SetRangeUser(110,140);

  double max = 1.2 * massRegressionV2->GetMaximum();

  massNoRegression->SetMaximum(max);
  massNoRegression->Draw("hist");
  massRegressionV0->Draw("hist,same");
  massRegressionV1->Draw("hist,same");
  massRegressionV2->Draw("hist,same");


  // Print Fit Values
  TLatex *tex = 0;

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



  cv->SaveAs(Form("HiggsMassHistogram_%s.gif",label.c_str()));


}

