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
#include "TAxis.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"

#include "hzzStyle.C"

using namespace RooFit;

void plotZResolution() {

  double etabins[5] = {0,1,1.479,2.0,2.5};
  TH2D *hScale = new TH2D("hScale","",4,etabins,2,-0.5,1.5);
  TH2D *hReso = new TH2D("hReso","",4,etabins,2,-0.5,1.5);

  for(int ieta=0; ieta<4; ++ieta) {
    for(int ir9=0; ir9<2; ++ir9) {

      cout << "Analyzing eta bin: " << ieta << "  and r9 bin: " << ir9 << endl;

    stringstream mcfile, datafile;
    mcfile << "mc42X_EtaBin" << ieta << "_R9Bin" << ir9 << ".root";
    datafile << "data2011_EtaBin" << ieta << "_R9Bin" << ir9 << ".root";
    
    TFile *tmcfile = TFile::Open(mcfile.str().c_str());
    RooFitResult *mcfr = (RooFitResult*)tmcfile->Get("fitres");
    float mcDM = ((RooRealVar*)(mcfr->floatParsFinal().find("#Deltam_{CB}")))->getVal();
    float mcDM_err = ((RooRealVar*)(mcfr->floatParsFinal().find("#Deltam_{CB}")))->getError();
    float mcS = ((RooRealVar*)(mcfr->floatParsFinal().find("#sigma_{CB}")))->getVal();
    float mcS_err = ((RooRealVar*)(mcfr->floatParsFinal().find("#sigma_{CB}")))->getError();

    TFile *tdatafile = TFile::Open(datafile.str().c_str());
    RooFitResult *datafr = (RooFitResult*)tdatafile->Get("fitres");
    float dataDM = ((RooRealVar*)(datafr->floatParsFinal().find("#Deltam_{CB}")))->getVal();
    float dataDM_err = ((RooRealVar*)(datafr->floatParsFinal().find("#Deltam_{CB}")))->getError();
    float dataS = ((RooRealVar*)(datafr->floatParsFinal().find("#sigma_{CB}")))->getVal();
    float dataS_err = ((RooRealVar*)(datafr->floatParsFinal().find("#sigma_{CB}")))->getError();
    
    float rM = (dataDM-mcDM)/(dataDM + 91.19);
    float rM_err = rM * sqrt(dataDM_err*dataDM_err + mcDM_err*mcDM_err);
    float rS = (dataS-mcS)/(dataS);
    cout << "rS = " << rS << endl;
    float rS_err = rS * sqrt(dataS_err*dataS_err + mcS_err*mcS_err);

    hScale->SetBinContent(ieta+1,ir9+1,rM);
    hScale->SetBinError(ieta+1,ir9+1,rM_err);
    hReso->SetBinContent(ieta+1,ir9+1,rS);
    hReso->SetBinError(ieta+1,ir9+1,rS_err);

    }
  }

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.3f");
  gStyle->SetTextFont(62);
  gStyle->SetTextSize(100);

  TCanvas *c1 = new TCanvas("c1","",600,600);

  hScale->GetXaxis()->SetTitle("electron |#eta|");
  hScale->GetYaxis()->SetTitle("R9 bin");
  hScale->SetTitle("#Delta m/m (data - sim.)");
  hScale->GetYaxis()->SetTitleOffset(2.0);
  hScale->GetYaxis()->SetNdivisions(4,kFALSE);  
  hScale->GetYaxis()->SetTitle("R_{9} class");
  hScale->SetMarkerSize(2);
  hScale->Draw("texte col");
  c1->SaveAs("delta-scale.pdf");
  c1->SaveAs("delta-scale.png");

  TCanvas *c2 = new TCanvas("c2","",600,600);
  hReso->GetXaxis()->SetTitle("electron |#eta|");
  hReso->GetYaxis()->SetTitle("R9 bin");
  hReso->SetTitle("#Delta #sigma/#sigma (data - sim.)");
  hReso->GetYaxis()->SetTitleOffset(2.0);
  hReso->GetYaxis()->SetNdivisions(4,kFALSE);  
  hReso->GetYaxis()->SetTitle("R_{9} class");
  hReso->SetMarkerSize(2);
  hReso->Draw("texte col");
  c2->SaveAs("delta-resolution.pdf");
  c2->SaveAs("delta-resolution.png");


}

void overlayDataMC(bool do7TeV) {

  TStyle *hzzstyle = getStyle("ZZ");
  hzzstyle->cd();

  double etabins[5] = {0,1,1.479,2.0,2.5};
  TH2D *hScale = new TH2D("hScale","",4,etabins,2,-0.5,1.5);
  TH2D *hReso = new TH2D("hReso","",4,etabins,2,-0.5,1.5);
  
  for(int ieta=0; ieta<4; ++ieta) {
    for(int ir9=0; ir9<2; ++ir9) {
      
      cout << "Analyzing eta bin: " << ieta << "  and r9 bin: " << ir9 << endl;
      
      stringstream mcfile, datafile;
      mcfile << (do7TeV ? "mc42x_EtaBin" : "mc_EtaBin") << ieta << "_R9Bin" << ir9 << ".root";
      datafile << (do7TeV ? "data2011_Jan22_EtaBin" : "data2012_Jan22_EtaBin") << ieta << "_R9Bin" << ir9 << ".root";


      TFile *tmcfile = TFile::Open(mcfile.str().c_str());
      RooWorkspace *wmc = (RooWorkspace*)tmcfile->Get("ZeeMassScaleAndResolutionFit");

      RooRealVar *zmass = (RooRealVar*)wmc->var("zmass");
      RooPlot* plot = zmass->frame(Range(75,105),Bins(40));
      plot->SetTitle("");

      RooAbsPdf* mcpdf = (RooAbsPdf*)wmc->pdf("model");
      RooDataHist *mch = (RooDataHist*)wmc->data("data_binned");
      float mcn = mch->sum(true);

      TFile *tdatafile = TFile::Open(datafile.str().c_str());
      RooWorkspace *wdata = (RooWorkspace*)tdatafile->Get("ZeeMassScaleAndResolutionFit");
      RooAbsPdf* datapdf = (RooAbsPdf*)wdata->pdf("model");
      RooDataHist *datah = (RooDataHist*)wdata->data("data_binned");
      float datan = datah->sum(true);

      datah->plotOn(plot,MarkerStyle(kFullCircle),MarkerColor(kRed+2),LineColor(kRed+2),Rescale(1./datan));
      datapdf->plotOn(plot,LineColor(kRed+2),Normalization(1./datan));
      datah->plotOn(plot,MarkerStyle(kFullCircle),MarkerColor(kRed+2),LineColor(kRed+2),Rescale(1./datan));

      mch->plotOn(plot,MarkerStyle(kFullSquare),MarkerColor(kBlue+2),LineColor(kBlue+2),Rescale(1./mcn));
      mcpdf->plotOn(plot,LineColor(kBlue+2),Normalization(1./mcn));
      mch->plotOn(plot,MarkerStyle(kFullSquare),MarkerColor(kBlue+2),LineColor(kBlue+2),Rescale(1./mcn));

      // parameters
      RooFitResult *mcfr = (RooFitResult*)tmcfile->Get("fitres");
      float mcDM = ((RooRealVar*)(mcfr->floatParsFinal().find("#Deltam_{CB}")))->getVal();
      float mcDM_err = ((RooRealVar*)(mcfr->floatParsFinal().find("#Deltam_{CB}")))->getError();
      float mcS = ((RooRealVar*)(mcfr->floatParsFinal().find("#sigma_{CB}")))->getVal();
      float mcS_err = ((RooRealVar*)(mcfr->floatParsFinal().find("#sigma_{CB}")))->getError();
      
      RooFitResult *datafr = (RooFitResult*)tdatafile->Get("fitres");
      float dataDM = ((RooRealVar*)(datafr->floatParsFinal().find("#Deltam_{CB}")))->getVal();
      float dataDM_err = ((RooRealVar*)(datafr->floatParsFinal().find("#Deltam_{CB}")))->getError();
      float dataS = ((RooRealVar*)(datafr->floatParsFinal().find("#sigma_{CB}")))->getVal();
      float dataS_err = ((RooRealVar*)(datafr->floatParsFinal().find("#sigma_{CB}")))->getError();
      
      float rM = (dataDM-mcDM)/(dataDM + 91.19);
      float rM_err = sqrt(pow(dataDM_err/dataDM,2)+pow(mcDM_err/mcDM,2))/91.19;
      float rS = (dataS-mcS)/dataS;
      float rS_err = sqrt(pow(dataS_err/dataS,2)+pow(mcS_err/mcS,2));
      
      stringstream fss;
      fss << "overlayZee" << (do7TeV ? "_2011" : "_2012") << "_EtaBin" << ieta << "_R9Bin" << ir9;

      TCanvas *c = new TCanvas("c","",600,600);
      float maxy=plot->GetMaximum()/datan;
      plot->GetYaxis()->SetRangeUser(0.0,1.5*maxy);
      plot->Draw();

      TLegend* legend = new TLegend(0.18, 0.77, 0.30, 0.87);
      
      legend->SetBorderSize(     0);
      legend->SetFillColor (  4000);
      legend->SetTextAlign (    12);
      legend->SetTextFont  (   132);
      legend->SetTextSize  (  0.04);

      TH1F *dummydata = new TH1F("dummydata","",1,0,1);
      TH1F *dummymc = new TH1F("dummymc","",1,0,1);
      dummydata->SetMarkerStyle(kFullCircle);
      dummydata->SetMarkerColor(kRed+2);
      dummymc->SetMarkerStyle(kFullSquare);
      dummymc->SetMarkerColor(kBlue+2);

      legend->AddEntry(dummydata,"data","pe");
      legend->AddEntry(dummymc,"simulation","pe");
      legend->Draw();

      TLatex* CP = new TLatex(75,1.52*maxy, do7TeV ? "CMS Preliminary                               #sqrt{s} = 7 TeV, L = 5.05 fb^{-1}" : "CMS Preliminary                               #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
      CP->SetTextSize(0.030);
      CP->SetTextFont  (   132);
      CP->SetTextSize  (   0.035);
      CP->Draw();

      TLatex *tex = new TLatex();
      tex->SetNDC();
      tex->SetTextSize(.1);
      tex->SetTextFont(132);
      tex->SetTextSize(0.037);
      tex->DrawLatex(0.55, 0.87, Form("#Delta M/M = %.3f #pm %.3f",rM,rM_err));
      tex->DrawLatex(0.55, 0.77, Form("#Delta #sigma/#sigma = %.2f #pm %.2f",rS,rS_err));

      TLine *line = new TLine(75,1.1*maxy,105,1.1*maxy);
      line->SetLineWidth(2.);
      line->Draw();

      c->Update();

      c->SaveAs((fss.str()+".pdf").c_str());
      c->SaveAs((fss.str()+".png").c_str());
    }
  }

}
