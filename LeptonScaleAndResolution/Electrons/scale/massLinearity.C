#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

#include "TH1F.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TString.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TRandom.h"

#include "histoFunc.h"

#include "hzzStyle.C"

enum eventclass {
  kInclusive=0,
  kG1,
  kG2
};

int checkclass=kInclusive;

// global variables: Zee
std::vector<TH1F*> h_Data_Mass_InEB, h_Data_Mass_OutEB, h_Data_Mass_EE;
std::vector<TH1F*> h_MC_Mass_InEB, h_MC_Mass_OutEB, h_MC_Mass_EE;
// int nptbins=15;
// float ptbins[16] = {10,15,19,23,27,31,35,39,43,47,51,55,59,63,67,71};
int nptbins=10;
float ptbins[11] = {7,16,22,28,34,40,46,52,58,64,70};
float etabins[4] = {0,0.8,1.479,2.5};

std::vector<TH1F*> h_Data_Mass_EB;
std::vector<TH1F*> h_MC_Mass_EB;
int nptbinsJ=4;
float ptbinsJ[5] = {7,10,13,16,19};
// int nptbinsJ=1;
// float ptbinsJ[2] = {7,19};
float etabinsJ[3] = {0,1.479,2.5};

int nptbinsU=3;
float ptbinsU[4] = {10,15,20,25};
// int nptbinsU=1;
// float ptbinsU[2] = {10,25};
float etabinsU[3] = {0,1.479,2.5};


using namespace std;

// ******************************************************
// Read data
// ******************************************************
void fillHistograms(bool do7TeV, bool isMC) {

  string fileMC;
  string fileData;
  if(do7TeV) {
    fileMC = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim_jan22/zEE_lineshape_MC_42X.root");
    fileData = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim_jan22/zEE_lineshape_Data_2011_lincorr.root");
  } else {
    // fileMC = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim/zEE_lineshape_MC_53X.root");
    //    fileData = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim/zEE_lineshape_Data_2012All.root");
    fileMC = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim_jan22/zEE_lineshape_MC_53X.root");
    fileData = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim_jan22/zEE_lineshape_Data_2012_lincorr.root");
  }

  TFile *tfile = 0;
  if(isMC) tfile = TFile::Open(fileMC.c_str());
  else tfile = TFile::Open(fileData.c_str());

  TTree *tree = (TTree*)tfile->Get("zeetree/probe_tree");
  int n=tree->GetEntries();

   Float_t         l1bdtID;
   Float_t         l1bdtIso;
   Float_t         l1classification;
   Float_t         l1p;
   Float_t         l1pdgId;
   Float_t         l1phi;
   Float_t         l1pt;
   Float_t         l1eta;
   Float_t         l1r9;
   Float_t         l2bdtID;
   Float_t         l2bdtIso;
   Float_t         l2classification;
   Float_t         l2p;
   Float_t         l2pdgId;
   Float_t         l2phi;
   Float_t         l2pt;
   Float_t         l2eta;
   Float_t         l2r9;
   Float_t         massErr;
   Float_t         numTrueInteractions;
   Float_t         nvtx;
   Float_t         rhoAA;
   Float_t         zeta;
   Float_t         zmass;
   Float_t         zmll;
   Float_t         zphi;
   Float_t         zpt;
   UInt_t          run;
   UInt_t          lumi;
   UInt_t          event;

   tree->SetBranchAddress("l1bdtID", &l1bdtID);
   tree->SetBranchAddress("l1bdtIso", &l1bdtIso);
   tree->SetBranchAddress("l1classification", &l1classification);
   tree->SetBranchAddress("l1eta", &l1eta);
   tree->SetBranchAddress("l1p", &l1p);
   tree->SetBranchAddress("l1pdgId", &l1pdgId);
   tree->SetBranchAddress("l1phi", &l1phi);
   tree->SetBranchAddress("l1pt", &l1pt);
   tree->SetBranchAddress("l1r9", &l1r9);
   tree->SetBranchAddress("l2bdtID", &l2bdtID);
   tree->SetBranchAddress("l2bdtIso", &l2bdtIso);
   tree->SetBranchAddress("l2classification", &l2classification);
   tree->SetBranchAddress("l2eta", &l2eta);
   tree->SetBranchAddress("l2p", &l2p);
   tree->SetBranchAddress("l2pdgId", &l2pdgId);
   tree->SetBranchAddress("l2phi", &l2phi);
   tree->SetBranchAddress("l2pt", &l2pt);
   tree->SetBranchAddress("l2r9", &l2r9);
   tree->SetBranchAddress("massErr", &massErr);
   tree->SetBranchAddress("nvtx", &nvtx);
   tree->SetBranchAddress("rhoAA", &rhoAA);
   tree->SetBranchAddress("zeta", &zeta);
   tree->SetBranchAddress("zmass", &zmass);
   tree->SetBranchAddress("zmll", &zmll);
   tree->SetBranchAddress("zphi", &zphi);
   tree->SetBranchAddress("zpt", &zpt);
   tree->SetBranchAddress("run", &run);
   tree->SetBranchAddress("lumi", &lumi);
   tree->SetBranchAddress("event", &event);
   if(isMC) tree->SetBranchAddress("numTrueInteractions", &numTrueInteractions);

   TH1F *massH = new TH1F("mass","",100,75,105);
   massH->SetMarkerColor(kBlack);
   massH->SetLineColor(kBlack);
   massH->SetMarkerStyle(8);
   massH->SetMarkerSize(1);

   h_Data_Mass_InEB.resize(nptbins);
   h_Data_Mass_OutEB.resize(nptbins);
   h_Data_Mass_EE.resize(nptbins);

   h_MC_Mass_InEB.resize(nptbins);
   h_MC_Mass_OutEB.resize(nptbins);
   h_MC_Mass_EE.resize(nptbins);

   //   TFile *histofile = TFile::Open("histos.root",(isMC ? "recreate" : "update"));
   
   for(int i=0;i<nptbins;++i) {
     stringstream hssInEB, hssOutEB, hssEE;
     hssInEB << "mass_ptbin" << i << "_InEB_";
     hssOutEB << "mass_ptbin" << i << "_OutEB_";
     hssEE << "mass_ptbin" << i << "_EE_";
     
     hssInEB << ((isMC) ? "MC" : "Data");
     hssOutEB << ((isMC) ? "MC" : "Data");
     hssEE << ((isMC) ? "MC" : "Data");

     if(isMC) {
       h_MC_Mass_InEB[i] = (TH1F*)massH->Clone(hssInEB.str().c_str());
       h_MC_Mass_OutEB[i] = (TH1F*)massH->Clone(hssOutEB.str().c_str());
       h_MC_Mass_EE[i] = (TH1F*)massH->Clone(hssEE.str().c_str());
     } else {
       h_Data_Mass_InEB[i] = (TH1F*)massH->Clone(hssInEB.str().c_str());
       h_Data_Mass_OutEB[i] = (TH1F*)massH->Clone(hssOutEB.str().c_str());
       h_Data_Mass_EE[i] = (TH1F*)massH->Clone(hssEE.str().c_str());
     }
   }

   for(int jentry=0; jentry<n; ++jentry) {
     tree->GetEntry(jentry);

     if(jentry%100000==0) cout << "analyzing event # " << jentry << endl; 

     if(l1pt<7 || fabs(l1eta)>2.5 ||
	l2pt<7 || fabs(l2eta)>2.5) continue;

     // fill one leg inclusive on the other one

     // if the lepton is not matched, the momentum is 0
     // if(isMC && (tag_gen_p==0 || gen_p==0)) continue;

     // pt bin
     int tagbin=-1; int probebin=-1;
     for(int b=0;b<nptbins;++b) {
       if(l2pt>=ptbins[b] && l2pt<ptbins[b+1]) probebin=b;
       if(l1pt>=ptbins[b] && l1pt<ptbins[b+1]) tagbin=b;
     }
     if(probebin==-1 || tagbin==-1) continue;

     // eta bin
     int tagetabin=-1; int probeetabin=-1;
     for(int b=0;b<4;++b) {
       if(fabs(l2eta)>=etabins[b] && fabs(l2eta)<etabins[b+1]) probeetabin=b;
       if(fabs(l1eta)>=etabins[b] && fabs(l1eta)<etabins[b+1]) tagetabin=b;
     }
     if(probeetabin==-1 || tagetabin==-1) continue;

     if((checkclass==kG1 && l2classification<2) ||
	(checkclass==kG2 && l2classification>=2)) {
	 switch (probeetabin) {
	 case 0:
	   if(isMC) h_MC_Mass_InEB[probebin]->Fill(zmass);
	   else h_Data_Mass_InEB[probebin]->Fill(zmass);
	   break;
	 case 1:
	   if(isMC) h_MC_Mass_OutEB[probebin]->Fill(zmass);
	 else h_Data_Mass_OutEB[probebin]->Fill(zmass);
	   break;
	 default:
	   if(isMC) h_MC_Mass_EE[probebin]->Fill(zmass);
	   else h_Data_Mass_EE[probebin]->Fill(zmass);
	   break;
	 }
       }


     if((checkclass==kG1 && l1classification<2) ||
	(checkclass==kG2 && l1classification>=2)) {
       switch (tagetabin) {
       case 0:
	 if(isMC) h_MC_Mass_InEB[tagbin]->Fill(zmass);
	 else h_Data_Mass_InEB[tagbin]->Fill(zmass);
	 break;
       case 1:
	 if(isMC) h_MC_Mass_OutEB[tagbin]->Fill(zmass);
	 else h_Data_Mass_OutEB[tagbin]->Fill(zmass);
	 break;
       default:
	 if(isMC) h_MC_Mass_EE[tagbin]->Fill(zmass);
	 else h_Data_Mass_EE[tagbin]->Fill(zmass);
	 break;
       }
     }
     
   } // loop on events

   // histofile->cd();
   // for(int i=0;i<nptbins;++i) {
   //   if(isMC) {
   //     h_MC_Mass_InEB[i]->Write();
   //     h_MC_Mass_OutEB[i]->Write();
   //     h_MC_Mass_EE[i]->Write();
   //   } else {
   //     h_Data_Mass_InEB[i]->Write();
   //     h_Data_Mass_OutEB[i]->Write();
   //     h_Data_Mass_EE[i]->Write();
   //   }
   // }
   // histofile->Close();
}

void myFitScale(bool do7TeV) {

  // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetOptStat(0);
  // ------------------------------ 
  
  fillHistograms(do7TeV,true);  
  fillHistograms(do7TeV,false);  

  TGraphErrors *gScaleInEB = new TGraphErrors();
  TGraphErrors *gRatioInEB = new TGraphErrors();
  TGraphErrors *gScaleOutEB = new TGraphErrors();
  TGraphErrors *gRatioOutEB = new TGraphErrors();
  TGraphErrors *gScaleEE = new TGraphErrors();
  TGraphErrors *gRatioEE = new TGraphErrors();

  TGraphErrors *gChiSqInEB[1000], *gChiSqOutEB[1000], *gChiSqEE[1000]; 

  histoFunc *templateHistoFuncInEB[h_MC_Mass_InEB.size()];
  histoFunc *templateHistoFuncOutEB[h_MC_Mass_OutEB.size()];
  histoFunc *templateHistoFuncEE[h_MC_Mass_EE.size()];
  std::vector <TF1*> f_Mass_InEB, f_Mass_OutEB, f_Mass_EE; 
  std::vector<double> val_Mass_InEB, val_Mass_OutEB, val_Mass_EE; 
  std::vector<double> err_Mass_InEB, err_Mass_OutEB, err_Mass_EE; 

  TCanvas c1("c1","",600,600);

  // loop on histograms and fit
  for (unsigned int i=0;i<h_Data_Mass_InEB.size();i++) {
    
    cout << h_Data_Mass_InEB[i]->GetName() << endl;
    
    int rebinningData = (i<1) ? 5 : 5;
    int rebinningMC =   (i<1) ? 4 : 1;

    h_MC_Mass_InEB[i]->Rebin(rebinningMC);
    h_MC_Mass_OutEB[i]->Rebin(rebinningMC);
    h_MC_Mass_EE[i]->Rebin(rebinningMC);

    h_MC_Mass_InEB[i]->Smooth(5);
    h_MC_Mass_OutEB[i]->Smooth(5);
    h_MC_Mass_EE[i]->Smooth(5);

    h_Data_Mass_InEB[i]->Rebin(rebinningData);
    h_Data_Mass_OutEB[i]->Rebin(rebinningData);
    h_Data_Mass_EE[i]->Rebin(rebinningData);

    /*
    h_MC_Mass_InEB[i]->Sumw2();
    h_MC_Mass_InEB[i]->Smooth();
    h_MC_Mass_InEB[i]->Rebin(rebinning);

    h_MC_Mass_OutEB[i]->Sumw2();
    h_MC_Mass_OutEB[i]->Smooth();
    h_MC_Mass_OutEB[i]->Rebin(rebinning);
    
    h_MC_Mass_EE[i]->Sumw2();
    h_MC_Mass_EE[i]->Smooth();
    h_MC_Mass_EE[i]->Rebin(rebinning);
    h_MC_Mass_EE[i]->Smooth();
    */

    h_MC_Mass_InEB[i]->Draw();
    stringstream fss1;
    fss1 << h_MC_Mass_InEB[i]->GetName() << ".pdf";
    c1.SaveAs(fss1.str().c_str());

    h_MC_Mass_OutEB[i]->Draw();
    stringstream fss2;
    fss2 << h_MC_Mass_OutEB[i]->GetName() << ".pdf";
    c1.SaveAs(fss2.str().c_str());
  
    h_MC_Mass_EE[i]->Draw();
    stringstream fss3;
    fss3 << h_MC_Mass_EE[i]->GetName() << ".pdf";
    c1.SaveAs(fss3.str().c_str());

    templateHistoFuncInEB[i] = new histoFunc(h_MC_Mass_InEB[i]); 
    templateHistoFuncOutEB[i] = new histoFunc(h_MC_Mass_OutEB[i]); 
    templateHistoFuncEE[i] = new histoFunc(h_MC_Mass_EE[i]); 

    char funcNameInEB[200], funcNameOutEB[200], funcNameEE[200];
    sprintf(funcNameInEB,"f_Mass_InEB_%d",i);
    f_Mass_InEB.push_back( new TF1(funcNameInEB, templateHistoFuncInEB[i], 75, 105, 3, "histoFuncInEB") );
    f_Mass_InEB.back() -> SetParNames("Norm","Scale Factor","Dummy"); 
    f_Mass_InEB.back() -> SetLineWidth(2); 
    f_Mass_InEB.back() -> SetLineColor(kRed+2); 
    f_Mass_InEB.back() -> SetNpx(1000); 

    sprintf(funcNameOutEB,"f_Mass_OutEB_%d",i);
    f_Mass_OutEB.push_back( new TF1(funcNameOutEB, templateHistoFuncOutEB[i], 75, 105, 3, "histoFuncOutEB") );
    f_Mass_OutEB.back() -> SetParNames("Norm","Scale Factor","Dummy"); 
    f_Mass_OutEB.back() -> SetLineWidth(2); 
    f_Mass_OutEB.back() -> SetLineColor(kRed+2); 
    f_Mass_OutEB.back() -> SetNpx(10000); 
  
    sprintf(funcNameEE,"f_Mass_EE_%d",i);
    f_Mass_EE.push_back( new TF1(funcNameEE, templateHistoFuncEE[i], 75, 105, 3, "histoFuncEE") );
    f_Mass_EE.back() -> SetParNames("Norm","Scale Factor","Dummy"); 
    f_Mass_EE.back() -> SetLineWidth(2); 
    f_Mass_EE.back() -> SetLineColor(kRed+2); 
    f_Mass_EE.back() -> SetNpx(10000); 

    double xNormInEB = h_Data_Mass_InEB[i]->GetSumOfWeights()/h_MC_Mass_InEB[i]->GetSumOfWeights() *
      h_Data_Mass_InEB[i]->GetBinWidth(1)/h_MC_Mass_InEB[i]->GetBinWidth(1); 
    double xNormOutEB = h_Data_Mass_OutEB[i]->GetSumOfWeights()/h_MC_Mass_OutEB[i]->GetSumOfWeights() *
      h_Data_Mass_OutEB[i]->GetBinWidth(1)/h_MC_Mass_OutEB[i]->GetBinWidth(1); 
    double xNormEE = h_Data_Mass_EE[i]->GetSumOfWeights()/h_MC_Mass_EE[i]->GetSumOfWeights() *
      h_Data_Mass_EE[i]->GetBinWidth(1)/h_MC_Mass_EE[i]->GetBinWidth(1); 

    f_Mass_InEB.back() -> FixParameter(0, xNormInEB);
    f_Mass_InEB.back() -> SetParameter(1,  gRandom->Gaus(1.,0.050) );
    f_Mass_InEB.back() -> SetParLimits(1,  0.95, 1.05 );
    f_Mass_InEB.back() -> FixParameter(2, 0);   
  
    f_Mass_OutEB.back() -> FixParameter(0, xNormOutEB);
    f_Mass_OutEB.back() -> SetParameter(1,  gRandom->Gaus(1.,0.050) );
    f_Mass_OutEB.back() -> SetParLimits(1,  0.95, 1.05 );
    f_Mass_OutEB.back() -> FixParameter(2, 0);   
  
    f_Mass_EE.back() -> FixParameter(0, xNormEE);
    f_Mass_EE.back() -> SetParameter(1,  gRandom->Gaus(1.,0.050) );
    f_Mass_EE.back() -> SetParLimits(1,  0.95, 1.05 );
    f_Mass_EE.back() -> FixParameter(2, 0);   

    // fit data with template
    cout << "======= FITTING NOW =======" << endl;
    cout << "--> in EB " << endl;
    TFitResultPtr rp1; 
    int status=0;
    int attempt=0;
    while(status<3) {
      cout << "attempt = " << attempt++ << endl;
      f_Mass_InEB.back() -> SetParameter(1,  gRandom->Gaus(1.,0.1) );
      rp1 = h_Data_Mass_InEB[i] -> Fit(funcNameInEB, "EHRS");
      status = rp1->CovMatrixStatus();
      cout << "COV M STAT = " << status << endl;
    }
    cout << "--> out EB " << endl;
    status=0;
    attempt=0;
    while(status<3) {
      cout << "attempt = " << attempt++ << endl;
      f_Mass_OutEB.back() -> SetParameter(1,  gRandom->Gaus(1.,0.1) );
      rp1 = h_Data_Mass_OutEB[i] -> Fit(funcNameOutEB, "EHRS");
      status = rp1->CovMatrixStatus();
      cout << "COV M STAT = " << status << endl;
    }
    cout << "--> EE " << endl;
    status = 0;
    attempt=0;
    while(status<3) {
      cout << "attempt = " << attempt++ << endl;
      f_Mass_EE.back() -> SetParameter(1,  gRandom->Gaus(1.,0.1) );
      rp1 = h_Data_Mass_EE[i] -> Fit(funcNameEE, "EHRS");
      status = rp1->CovMatrixStatus();
      cout << "COV M STAT = " << status << endl;
    }
    cout << "======= FITTED =======" << endl;

    float pTmean = ptbins[i];

    float MassscaleInEB = 1./f_Mass_InEB.back()->GetParameter(1);
    float MasserrorInEB = f_Mass_InEB.back()->GetParError(1)/pow(f_Mass_InEB.back()->GetParameter(1),2.); 
    val_Mass_InEB.push_back(MassscaleInEB); err_Mass_InEB.push_back(MasserrorInEB); 
    float MassscaleOutEB = 1./f_Mass_OutEB.back()->GetParameter(1);
    float MasserrorOutEB = f_Mass_OutEB.back()->GetParError(1)/pow(f_Mass_OutEB.back()->GetParameter(1),2.); 
    val_Mass_OutEB.push_back(MassscaleOutEB); err_Mass_OutEB.push_back(MasserrorOutEB); 
    float MassscaleEE = 1./f_Mass_EE.back()->GetParameter(1);
    float MasserrorEE = f_Mass_EE.back()->GetParError(1)/pow(f_Mass_EE.back()->GetParameter(1),2.); 
    val_Mass_EE.push_back(MassscaleEE); err_Mass_EE.push_back(MasserrorEE); 

    cout << pTmean << "GeV/c : " << endl
	 << "InEB = " << MassscaleInEB << "+/-" << MasserrorInEB << endl
	 << "OutEB = " << MassscaleOutEB << "+/-" << MasserrorOutEB << endl
	 << "EE = " << MassscaleEE << "+/-" << MasserrorEE << endl;
  
    // derive scale and fill graph
    gScaleInEB->SetPoint(i,pTmean,MassscaleInEB-1.); 
    gScaleInEB->SetPointError(i,(ptbins[2]-ptbins[1])/2.,MasserrorInEB); 
    gRatioInEB->SetPoint(i,pTmean,h_Data_Mass_InEB[i]->GetMean()/h_MC_Mass_InEB[i]->GetMean()-1.);
    
    gScaleOutEB->SetPoint(i,pTmean,MassscaleOutEB-1.); 
    gScaleOutEB->SetPointError(i,(ptbins[2]-ptbins[1])/2.,MasserrorOutEB); 
    gRatioOutEB->SetPoint(i,pTmean,h_Data_Mass_OutEB[i]->GetMean()/h_MC_Mass_OutEB[i]->GetMean()-1.);
    
    gScaleEE->SetPoint(i,pTmean,MassscaleEE-1.); 
    gScaleEE->SetPointError(i,(ptbins[2]-ptbins[1])/2.,MasserrorEE); 
    gRatioEE->SetPoint(i,pTmean,h_Data_Mass_EE[i]->GetMean()/h_MC_Mass_EE[i]->GetMean()-1.);

    // scan 
    // char chiName[80];

    // gChiSqInEB[i] = new TGraphErrors(); 
    // gChiSqInEB[i]->SetMarkerStyle(20+i);
    // gChiSqInEB[i]->SetMarkerColor(kBlue+1);
    // sprintf(chiName,"chi2_MassInEB_Pt%d",i);
    // gChiSqInEB[i]->SetTitle(chiName); 

    // gChiSqOutEB[i] = new TGraphErrors(); 
    // gChiSqOutEB[i]->SetMarkerStyle(20+i);
    // gChiSqOutEB[i]->SetMarkerColor(kBlue+1);
    // sprintf(chiName,"chi2_MassOutEB_Pt%d",i);
    // gChiSqOutEB[i]->SetTitle(chiName); 

    // gChiSqEE[i] = new TGraphErrors(); 
    // gChiSqEE[i]->SetMarkerStyle(20+i);
    // gChiSqEE[i]->SetMarkerColor(kBlue+1);
    // sprintf(chiName,"chi2_MassEE_Pt%d",i);
    // gChiSqEE[i]->SetTitle(chiName); 

    // for (int iScale=-20; iScale<21; iScale++) {
    
    //   float chichi = f_Mass_InEB.back()->GetChisquare();
    //   f_Mass_InEB.back() -> FixParameter(1, 1./MassscaleInEB + 0.001 * iScale );
    //   TFitResultPtr rpA = h_Data_Mass_InEB[i] -> Fit(funcNameInEB, "QMRL");
    //   gChiSqInEB[i]->SetPoint( iScale+20, 1./(1./MassscaleInEB + 0.001 * iScale),  f_Mass_InEB.back()->GetChisquare()-chichi ); 

    //   chichi = f_Mass_OutEB.back()->GetChisquare();
    //   f_Mass_OutEB.back() -> FixParameter(1, 1./MassscaleOutEB + 0.001 * iScale );
    //   TFitResultPtr rpB = h_Data_Mass_OutEB[i] -> Fit(funcNameOutEB, "QMRL");
    //   gChiSqOutEB[i]->SetPoint( iScale+20, 1./(1./MassscaleOutEB + 0.001 * iScale),  f_Mass_OutEB.back()->GetChisquare()-chichi ); 

    //   chichi = f_Mass_EE.back()->GetChisquare();
    //   f_Mass_EE.back() -> FixParameter(1, 1./MassscaleEE + 0.001 * iScale );
    //   TFitResultPtr rpC = h_Data_Mass_EE[i] -> Fit(funcNameEE, "QMRL");
    //   gChiSqEE[i]->SetPoint( iScale+20, 1./(1./MassscaleEE + 0.001 * iScale),  f_Mass_EE.back()->GetChisquare()-chichi ); 

    // }
  }
  
  TLegend* legend = new TLegend(0.15, 0.70, 0.30, 0.85);
  
  legend->SetBorderSize(     0);
  legend->SetFillColor (  4000);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  (0.05);

  legend->AddEntry(h_Data_Mass_InEB[0],"data");
  legend->AddEntry(f_Mass_InEB[0],"MC fit to data","l");

  TPaveStats** s_Mass = new TPaveStats*[500];
  TCanvas *c2 = new TCanvas("c2","",600,600); 
  for(unsigned int i = 0; i < h_Data_Mass_InEB.size(); ++i) {
    char canvasName[50];
    
    int nTeV = (do7TeV) ? 7 : 8;

    // InEB
    TPaveText *bin0 = new TPaveText(0.60,0.75,0.90,0.85,"brNDC");
    stringstream binval0;
    if(i!=h_Data_Mass_InEB.size()-1) binval0 << ptbins[i] << " < p_{T} < " << ptbins[i+1] << " GeV"; 
    else binval0 << "p_{T} > " << ptbins[i] << "GeV";
    bin0->AddText(binval0.str().c_str());
    bin0->AddText("|#eta|<1");
    bin0->SetBorderSize(0);
    bin0->SetFillStyle(0);
    bin0->SetTextAlign(12);
    bin0->SetTextFont(132);
    bin0->SetTextSize(0.04);

    sprintf(canvasName, "Fits-InEB-Pt%d-%dTeV.pdf", i, nTeV); 
    h_Data_Mass_InEB[i] -> GetXaxis() -> SetTitle("m_{ee} (GeV)");
    h_Data_Mass_InEB[i] -> GetXaxis() -> SetRangeUser(75,105); 
    h_Data_Mass_InEB[i] -> Draw("e");
    gPad->Update(); 
    TPaveText *sigmat0 = new TPaveText(0.30,0.90,0.77,0.95,"brNDC");
    stringstream sigmaval0;
    sigmaval0 << fixed;
    sigmaval0 << setprecision(3);
    sigmaval0 << "#Delta M/M (data-sim.) " << val_Mass_InEB[i] << " #pm " << err_Mass_InEB[i];
    sigmat0->AddText(sigmaval0.str().c_str());
    sigmat0->SetBorderSize(0);
    sigmat0->SetFillStyle(0);
    sigmat0->SetTextAlign(12);
    sigmat0->SetTextFont(132);
    sigmat0->SetTextSize(0.04);
    // s_Mass[i]= (TPaveStats*)(h_Data_Mass_InEB[i]->GetListOfFunctions()->FindObject("stats"));
    // s_Mass[i]->SetTextColor(kRed+2);
    // s_Mass[i]->SetOptStat(101);
    // s_Mass[i]->SetOptFit(1111);
    // s_Mass[i]->SetY1NDC(0.65); //new x start position
    // s_Mass[i]->SetX1NDC(0.55); 
    f_Mass_InEB[i]->Draw("same C");
    legend->Draw();
    sigmat0->Draw();
    bin0->Draw();
    c2->SaveAs(canvasName);

    // OutEB
    TPaveText *bin1 = new TPaveText(0.60,0.75,0.90,0.85,"brNDC");
    stringstream binval1;
    if(i!=h_Data_Mass_InEB.size()-1) binval1 << ptbins[i] << " < p_{T} < " << ptbins[i+1] << " GeV"; 
    else binval1 << "p_{T} > " << ptbins[i] << "GeV";
    bin1->AddText(binval0.str().c_str());
    bin1->AddText("1<|#eta|<1.479");
    bin1->SetBorderSize(0);
    bin1->SetFillStyle(0);
    bin1->SetTextAlign(12);
    bin1->SetTextFont(132);
    bin1->SetTextSize(0.04);

    sprintf(canvasName, "Fits-OutEB-Pt%d-%dTeV.pdf", i, nTeV); 
    h_Data_Mass_OutEB[i] -> GetXaxis() -> SetTitle("m_{ee} (GeV)");
    h_Data_Mass_OutEB[i] -> GetXaxis() -> SetRangeUser(75,105); 
    h_Data_Mass_OutEB[i] -> Draw("e");
    gPad->Update(); 
    TPaveText *sigmat1 = new TPaveText(0.30,0.90,0.77,0.95,"brNDC");
    stringstream sigmaval1;
    sigmaval1 << fixed;
    sigmaval1 << setprecision(3);
    sigmaval1 << "#Delta M/M (data-sim.) " << val_Mass_OutEB[i] << " #pm " << err_Mass_OutEB[i];
    sigmat1->AddText(sigmaval0.str().c_str());
    sigmat1->SetBorderSize(0);
    sigmat1->SetFillStyle(0);
    sigmat1->SetTextAlign(12);
    sigmat1->SetTextFont(132);
    sigmat1->SetTextSize(0.04);
    // s_Mass[i]= (TPaveStats*)(h_Data_Mass_OutEB[i]->GetListOfFunctions()->FindObject("stats"));
    // s_Mass[i]->SetTextColor(kRed+2);
    // s_Mass[i]->SetOptStat(101);
    // s_Mass[i]->SetOptFit(1111);
    // s_Mass[i]->SetY1NDC(0.65); //new x start position
    // s_Mass[i]->SetX1NDC(0.55); 
    f_Mass_OutEB[i]->Draw("sameC");
    legend->Draw();
    sigmat1->Draw();
    bin1->Draw();
    c2->SaveAs(canvasName);

    // EE
    TPaveText *bin2 = new TPaveText(0.60,0.75,0.90,0.85,"brNDC");
    stringstream binval2;
    if(i!=h_Data_Mass_InEB.size()-1) binval2 << ptbins[i] << " < p_{T} < " << ptbins[i+1] << " GeV"; 
    else binval2 << "p_{T} > " << ptbins[i] << "GeV";
    bin2->AddText(binval0.str().c_str());
    bin2->AddText("|#eta|>1.479");
    bin2->SetBorderSize(0);
    bin2->SetFillStyle(0);
    bin2->SetTextAlign(12);
    bin2->SetTextFont(132);
    bin2->SetTextSize(0.04);

    sprintf(canvasName, "Fits-EE-Pt%d-%dTeV.pdf", i, nTeV); 
    h_Data_Mass_EE[i] -> GetXaxis() -> SetTitle("m_{ee} (GeV)");
    h_Data_Mass_EE[i] -> GetXaxis() -> SetRangeUser(75,105); 
    h_Data_Mass_EE[i] -> Draw("e");
    gPad->Update(); 
    TPaveText *sigmat2 = new TPaveText(0.30,0.90,0.77,0.95,"brNDC");
    stringstream sigmaval2;
    sigmaval2 << fixed;
    sigmaval2 << setprecision(3);
    sigmaval2 << "#Delta M/M (data-sim.) " << val_Mass_EE[i] << " #pm " << err_Mass_EE[i];
    sigmat2->AddText(sigmaval0.str().c_str());
    sigmat2->SetBorderSize(0);
    sigmat2->SetFillStyle(0);
    sigmat2->SetTextAlign(12);
    sigmat2->SetTextFont(132);
    sigmat2->SetTextSize(0.04);
    // s_Mass[i]= (TPaveStats*)(h_Data_Mass_EE[i]->GetListOfFunctions()->FindObject("stats"));
    // s_Mass[i]->SetTextColor(kRed+2);
    // s_Mass[i]->SetOptStat(101);
    // s_Mass[i]->SetOptFit(1111);
    // s_Mass[i]->SetY1NDC(0.65); //new x start position
    // s_Mass[i]->SetX1NDC(0.55); 
    f_Mass_EE[i]->Draw("sameC");
    legend->Draw();
    sigmat2->Draw();
    bin2->Draw();
    c2->SaveAs(canvasName);
  }

  TCanvas *c = new TCanvas("c","c");
  c->cd();
  gScaleInEB->GetYaxis()->SetRangeUser(-0.015,0.015);
  gScaleInEB->GetYaxis()->SetTitle("(M_{ee}/m_{Z})^{data}-(M_{ee}/m_{Z})^{MC}");
  gScaleInEB->GetYaxis()->SetTitleOffset(1.5);
  gScaleInEB->GetXaxis()->SetTitle("electron p_{T} [GeV]");
  //  gScaleInEB->SetMarkerColor(kRed+2);
  gScaleInEB->SetMarkerStyle(20);
  gScaleInEB->Draw("AP");
  c->SaveAs(((do7TeV) ? "scaleInEB-7TeV.pdf" : "scaleInEB-8TeV.pdf"));

  gScaleOutEB->GetYaxis()->SetRangeUser(-0.015,0.015);
  gScaleOutEB->GetYaxis()->SetTitle("(M_{ee}/m_{Z})^{data}-(M_{ee}/m_{Z})^{MC}"); 
  gScaleOutEB->GetYaxis()->SetTitleOffset(1.5);
  gScaleOutEB->GetXaxis()->SetTitle("electron p_{T} [GeV]");
  //  gScaleOutEB->SetMarkerColor(kRed+2);
  gScaleOutEB->SetMarkerStyle(20);
  gScaleOutEB->Draw("AP");
  c->SaveAs(((do7TeV) ? "scaleOutEB-7TeV.pdf" : "scaleOutEB-8TeV.pdf"));

  gScaleEE->GetYaxis()->SetRangeUser(-0.03,0.015);
  gScaleEE->GetYaxis()->SetTitle("(M_{ee}/m_{Z})^{data}-(M_{ee}/m_{Z})^{MC}");
  gScaleEE->GetYaxis()->SetTitleOffset(1.5);
  gScaleEE->GetXaxis()->SetTitle("electron p_{T} [GeV]");
  //  gScaleEE->SetMarkerColor(kRed+2);
  gScaleEE->SetMarkerStyle(20);
  gScaleEE->Draw("AP");
  c->SaveAs(((do7TeV) ? "scaleEE-7TeV.pdf" : "scaleEE-8TeV.pdf"));

  TFile *ff =new TFile(((do7TeV) ? "gScaleSyst-7TeV.root" : "gScaleSyst-8TeV.root"),"RECREATE");
  ff->cd();
  gScaleInEB->SetName("gScaleInEB");
  gScaleInEB->Write();
  gScaleOutEB->SetName("gScaleOutEB");
  gScaleOutEB->Write();
  gScaleEE->SetName("gScaleEE");
  gScaleEE->Write();
  // for(unsigned int i = 0; i < h_Data_Mass_InEB.size(); ++i) {
  //   gChiSqInEB[i]->Write();
  //   gChiSqOutEB[i]->Write();
  //   gChiSqEE[i]->Write();
  // }
  
}


void massLinearity() {

  //  myFitScale(true);
  myFitScale(false);

}


double getShift(float eta, float pt, TGraphErrors* gScaleInEB, TGraphErrors* gScaleOutEB, TGraphErrors* gScaleEE, bool getErr=false) {

  double shift=-10000;
  int ptbin=-1;
  for(int b=0;b<nptbins;++b) {
    if(pt>=ptbins[b] && pt<ptbins[b+1]) ptbin=b;
  }
  if(pt<15) ptbin=0;
  else if(pt>ptbins[nptbins]) ptbin=nptbins;

  double x;
  if(!getErr) {
    if(fabs(eta)<0.8) gScaleInEB->GetPoint(ptbin,x,shift);
    else if(fabs(eta)>=0.8 && fabs(eta)<1.479) gScaleOutEB->GetPoint(ptbin,x,shift);
    else gScaleEE->GetPoint(ptbin,x,shift);
  } else {
    if(fabs(eta)<0.8) shift = gScaleInEB->GetErrorY(ptbin);
    else if(fabs(eta)>=0.8 && fabs(eta)<1.479) shift = gScaleOutEB->GetErrorY(ptbin);
    else shift = gScaleEE->GetErrorY(ptbin);
  }

  return shift;

}

// Quadratic background function
Double_t background(Double_t *x, Double_t *par) {
   return par[0]*x[0] + par[1]*x[0]*x[0];
}

Double_t cryball( Double_t *x, Double_t * par)
{

  // par[0] = normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = alpha
  // par[4] = n


  double cryball;
  double test = par[1] + par[2]*par[3];
  double xp = (x[0] - par[1]) / par[2];
  double A = TMath::Power((par[4]/par[3]),par[4]) * TMath::Exp(-0.5*par[3]*par[3]);
  double  B = ( par[4]/par[3] ) - par[3];


  if (x[0] < test)

    {
      cryball = par[0]*TMath::Exp(-0.5*xp*xp);

    } else {

      cryball =par[0]* A * TMath::Power(( B + xp ),-par[4]);

   }
  return cryball;
}

Double_t doublecryball( Double_t *x, Double_t * par) {

  // par[0] = normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = alpha1
  // par[4] = n1
  // par[5] = alpha2
  // par[6] = n2

  Double_t norm = par[0];
  Double_t mean = par[1];
  Double_t width = par[2];
  Double_t alpha1 = par[3];
  Double_t n1 = par[4];
  Double_t alpha2 = par[5];
  Double_t n2 = par[5];

  double t = (x[0]-mean)/width;
   if(t>-alpha1 && t<alpha2){
     return norm * exp(-0.5*t*t);
   }else if(t<-alpha1){
     double A1 = pow(n1/fabs(alpha1),n1)*exp(-alpha1*alpha1/2);
     double B1 = n1/fabs(alpha1)-fabs(alpha1);
     return norm*A1*pow(B1-t,-n1);
   }else if(t>alpha2){
     double A2 = pow(n2/fabs(alpha2),n2)*exp(-alpha2*alpha2/2);
     double B2 = n2/fabs(alpha2)-fabs(alpha2);
     return norm*A2*pow(B2+t,-n2);
   }else{
     cout << "ERROR evaluating range..." << endl;
     return 99;
   }

}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  return background(x,par) + cryball(x,&par[2]);
}

Double_t fitFunctionUpsilon(Double_t *x, Double_t *par) {
  
  // par[0] = norm1S
  // par[1] = mean1S
  // par[2] = sigma1S
  // par[3] = alpha
  // par[4] = n
  // par[5] = norm2S+3S
  // par[6] = mean2S+3S
  // par[7] = sigma2S+3S
  // par[8] = p1
  // par[9] = p2


  double mean1S = par[1];
  double mean2S = par[6];
  double test1S = mean1S + par[2]*par[3];
  double xp1S = (x[0] - mean1S) / par[2];
  double test2S = mean2S + par[7]*par[3];
  double xp2S = (x[0] - mean2S) / par[7];
  double A = TMath::Power((par[4]/par[3]),par[4]) * TMath::Exp(-0.5*par[3]*par[3]);
  double  B = ( par[4]/par[3] ) - par[3];

  // 1S
  double cryball1S;
  if (x[0] < test1S)

    {
      cryball1S = par[0]*TMath::Exp(-0.5*xp1S*xp1S);

    } else {

      cryball1S =par[0]* A * TMath::Power(( B + xp1S ),-par[4]);

   }

  // 2S 
  double cryball2S;
  if (x[0] < test2S)
    
    {
      cryball2S = par[5]*TMath::Exp(-0.5*xp2S*xp2S);

    } else {

      cryball2S =par[5]* A * TMath::Power(( B + xp2S ),-par[4]);

   }

  return par[8]*x[0] + par[9]*x[0]*x[0] + cryball1S + cryball2S;

}

void propagateErrorToH4l(int ch, bool do7TeV) {

  TStyle *hzzstyle = getStyle("ZZ");
  hzzstyle->SetOptFit(0111);
  hzzstyle->SetOptStat(0);
  hzzstyle->cd();

  stringstream filename;
  if(do7TeV) filename << "/Users/emanuele/Work/data/hzz4l/analysisstep/legacypaper/HZZ4L_44X_S1_V18_S2_V10/MC/hzzTree_id1126.root";
  else filename << "/Users/emanuele/Work/data/hzz4l/analysisstep/legacypaper/HZZ4L_53X_S1_V18_S2_V10/MC/hzzTree_id1126.root";

  TFile *file = TFile::Open(filename.str().c_str());
  TTree *tree = (TTree*)file->Get("zz4lTree/probe_tree");

  Float_t         l1pdgId; 
  Float_t         l1eta;
  Float_t         l1phi;
  Float_t         l1pt;
  Float_t         l2pdgId; 
  Float_t         l2eta;
  Float_t         l2phi;
  Float_t         l2pt;
  Float_t         l3pdgId; 
  Float_t         l3eta;
  Float_t         l3phi;
  Float_t         l3pt;
  Float_t         l4pdgId; 
  Float_t         l4eta;
  Float_t         l4phi;
  Float_t         l4pt;
  Float_t         mass;
  Float_t         channel;
  
  tree->SetBranchAddress("l1pdgId", &l1pdgId);
  tree->SetBranchAddress("l1eta", &l1eta);
  tree->SetBranchAddress("l1phi", &l1phi);
  tree->SetBranchAddress("l1pt", &l1pt);
  tree->SetBranchAddress("l2pdgId", &l2pdgId);
  tree->SetBranchAddress("l2eta", &l2eta);
  tree->SetBranchAddress("l2phi", &l2phi);
  tree->SetBranchAddress("l2pt", &l2pt);
  tree->SetBranchAddress("l3pdgId", &l3pdgId);
  tree->SetBranchAddress("l3eta", &l3eta);
  tree->SetBranchAddress("l3phi", &l3phi);
  tree->SetBranchAddress("l3pt", &l3pt);
  tree->SetBranchAddress("l4pdgId", &l4pdgId);
  tree->SetBranchAddress("l4eta", &l4eta);
  tree->SetBranchAddress("l4phi", &l4phi);
  tree->SetBranchAddress("l4pt", &l4pt);
  tree->SetBranchAddress("mass", &mass);
  tree->SetBranchAddress("channel", &channel);


  int n=tree->GetEntries();
  
  TH1F *hmass = new TH1F("hmass","",50,115,135);
  TH1F *hmassSh = new TH1F("hmassSh","",50,115,135);
  TH1F *hmassDiff = new TH1F("hmassDiff","",50,-1.5,0.5);

  hmass->GetXaxis()->SetTitle("m_{4l} (GeV)");
  hmassSh->GetXaxis()->SetTitle("m_{4l} (GeV)");
  hmass->GetYaxis()->SetTitle("events");
  hmassSh->GetYaxis()->SetTitle("events");
  hmass->GetYaxis()->SetTitleOffset(1.8);
  hmassSh->GetYaxis()->SetTitleOffset(1.8);
  
  stringstream filesyst;
  filesyst << "gScaleSyst-" << (do7TeV ? "7TeV" : "8TeV") << ".root";

  TFile *tfilesyst = TFile::Open(filesyst.str().c_str());
  TGraphErrors *gScaleInEB = (TGraphErrors*)tfilesyst->Get("gScaleInEB");
  TGraphErrors *gScaleOutEB = (TGraphErrors*)tfilesyst->Get("gScaleOutEB");
  TGraphErrors *gScaleEE = (TGraphErrors*)tfilesyst->Get("gScaleEE");

  for(int jentry=0; jentry<n; ++jentry) {
     tree->GetEntry(jentry);

     if(jentry%100000==0) cout << "analyzing event # " << jentry << endl; 
     
     if(channel!=ch) continue;

     float l1m = (fabs(l1pdgId)==11) ? 0.00051 : 0.10566;
     float l2m = (fabs(l2pdgId)==11) ? 0.00051 : 0.10566;
     float l3m = (fabs(l3pdgId)==11) ? 0.00051 : 0.10566;
     float l4m = (fabs(l4pdgId)==11) ? 0.00051 : 0.10566;

     TLorentzVector l1, l2, l3, l4;
     l1.SetPtEtaPhiM(l1pt,l1eta,l1phi,l1m);
     l2.SetPtEtaPhiM(l2pt,l2eta,l2phi,l2m);
     l3.SetPtEtaPhiM(l3pt,l3eta,l3phi,l3m);
     l4.SetPtEtaPhiM(l4pt,l4eta,l4phi,l4m);

     float m4l = (l1+l2+l3+l4).M();
     double l1s = (fabs(l1pdgId)==11) ? getShift(l1eta,l1pt,gScaleInEB,gScaleOutEB,gScaleEE) : 0.0;
     double l2s = (fabs(l2pdgId)==11) ? getShift(l2eta,l2pt,gScaleInEB,gScaleOutEB,gScaleEE) : 0.0;
     double l3s = (fabs(l3pdgId)==11) ? getShift(l3eta,l3pt,gScaleInEB,gScaleOutEB,gScaleEE) : 0.0;
     double l4s = (fabs(l4pdgId)==11) ? getShift(l4eta,l4pt,gScaleInEB,gScaleOutEB,gScaleEE) : 0.0;

     // cout << l1s << "  " << l2s << "  " << l3s << "  " << l4s << endl;

     double l1ptn = l1pt*(1+l1s);
     double l2ptn = l2pt*(1+l2s);
     double l3ptn = l3pt*(1+l3s);
     double l4ptn = l4pt*(1+l4s);

     TLorentzVector l1n, l2n, l3n, l4n;
     l1n.SetPtEtaPhiM(l1ptn,l1eta,l1phi,l1m);
     l2n.SetPtEtaPhiM(l2ptn,l2eta,l2phi,l2m);
     l3n.SetPtEtaPhiM(l3ptn,l3eta,l3phi,l3m);
     l4n.SetPtEtaPhiM(l4ptn,l4eta,l4phi,l4m);

     float m4ln = (l1n+l2n+l3n+l4n).M();

     // cout << "mass = " << m4l << "   mass shifted = " << m4ln << "  diff = " << m4ln-m4l << endl;

     hmass->Fill(m4l);
     hmassSh->Fill(m4ln);
     hmassDiff->Fill(m4ln-m4l);

  } // loop over events

  TCanvas *c1 = new TCanvas("c1","",725,725);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.12);
  c1->SetTopMargin(0.11);
  c1->SetBottomMargin(0.12);
  hmass->SetLineColor(kBlack);
  hmassSh->SetLineColor(kBlue+2);
  hmass->SetLineWidth(2);
  hmassSh->SetLineWidth(2);

  TLatex* CP = new TLatex(0.15,0.95, "CMS Simulation                          #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
  CP->SetNDC(kTRUE);
  CP->SetTextSize(0.030);

  TF1 *func = new TF1("dcryball",doublecryball,115,135,7); 
  func->SetParNames("N","m","#sigma","#alpha_{1}","n_{1}","#alpha_{2}","n_{2}");
  func->SetParameter(1,hmass->GetMean());
  func->SetParLimits(1,124,126);
  func->SetParameter(2,hmass->GetRMS());
  func->SetParameter(3,1.2);
  func->SetParameter(4,6);
  func->SetParameter(5,1.2);
  func->FixParameter(6,20);

  func->SetParLimits(3,0.1,4);
  func->SetParLimits(4,40,100);
  func->SetParLimits(5,1,5);

  hmass->Draw("");
  hmassSh->Draw("sames");

  func->SetLineColor(kBlack);
  hmass->Fit("dcryball","","same",115,135);
  func->SetLineColor(kBlue+2);
  hmassSh->Fit("dcryball","","same",115,135);
  
  c1->Update();
  
  TLegend *legend = new TLegend(0.15,0.75,0.40,0.85,NULL,"brNDC");
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  
  legend->AddEntry(hmass, "nominal");
  legend->AddEntry(hmassSh, "scale syst.");

  TPaveStats *p1 = (TPaveStats*)hmass->GetListOfFunctions()->FindObject("stats");
  p1->SetTextColor(kBlack);
  p1->SetX1NDC(0.7);
  p1->SetX2NDC(0.9);
  p1->SetY1NDC(0.7);
  p1->SetY2NDC(0.9);
  p1->Draw();

  TPaveStats *p2 = (TPaveStats*)hmassSh->GetListOfFunctions()->FindObject("stats");
  p2->SetTextColor(kBlue+2);
  p2->SetX1NDC(0.7);
  p2->SetX2NDC(0.9);
  p2->SetY1NDC(0.5);
  p2->SetY2NDC(0.7);
  p2->Draw();

  legend->Draw();
  CP->Draw();

  stringstream filenameout;
  if(do7TeV) filenameout << "mH125Shift-7TeV-ch" << ch << ".pdf";
  else filenameout << "mH125Shift-8TeV-ch" << ch << ".pdf";
  c1->SaveAs(filenameout.str().c_str());

  /*
  TCanvas *c2 = new TCanvas("c2","",600,600);
  hmassDiff->SetLineColor(kBlack);
  hmassDiff->SetLineWidth(2);

  func->SetParameter(1,0);
  func->SetParameter(2,hmassDiff->GetRMS());
  func->SetParameter(3,1.2);
  func->SetParameter(4,50);
  func->SetParameter(5,1.2);
  func->SetParameter(6,20);

  func->SetParLimits(3,0.1,4);
  func->SetParLimits(4,20,100);
  func->SetParLimits(5,1,5);

  hmassDiff->Draw();

  func->SetLineColor(kBlack);
  hmassDiff->Fit("dcryball","","same",-1.5,0.5);

  stringstream filenameoutDiff;
  if(do7TeV) filenameoutDiff << "mH125Diff-7TeV-ch" << ch << ".pdf";
  else filenameoutDiff << "mH125Diff-8TeV-ch" << ch << ".pdf";
  c2->SaveAs(filenameoutDiff.str().c_str());
  */

}


void create2DMaps(bool do7TeV) {

  stringstream filesyst;
  if(do7TeV) filesyst << "gScaleSyst-7TeV.root";
  else filesyst << "gScaleSyst-8TeV.root";

  TFile *tfilesyst = TFile::Open(filesyst.str().c_str());
  TGraphErrors *gScaleInEB = (TGraphErrors*)tfilesyst->Get("gScaleInEB");
  TGraphErrors *gScaleOutEB = (TGraphErrors*)tfilesyst->Get("gScaleOutEB");
  TGraphErrors *gScaleEE = (TGraphErrors*)tfilesyst->Get("gScaleEE");

  stringstream fileout;
  if(do7TeV) fileout << "elescale-syst-7TeV.root";
  else fileout << "elescale-syst-8TeV.root";

  TFile *filemap = TFile::Open(fileout.str().c_str(),"recreate");
  
  TH2D *map = new TH2D("elescale_syst","",nptbins,ptbins,3,etabins);
  
  for(int p=0;p<=nptbins+1;++p) {
    int ptbin;
    if(p<=1) ptbin=0; // put in the underflow the first value
    else if(p>=nptbins) ptbin=nptbins-1; // put in the overflow the last value
    else ptbin=p-1;

    for(int e=1;e<=4;++e) {

      double shift=-1000.;
      double x;
      if(e==1) gScaleInEB->GetPoint(ptbin,x,shift);
      else if(e==2) gScaleOutEB->GetPoint(ptbin,x,shift);
      else gScaleEE->GetPoint(ptbin,x,shift); // put in the overflow the endcap value
      
      map->SetBinContent(p,e,shift);

    }
  }
  
  filemap->cd();
  map->Write();
  filemap->Close();

  return;

}


void paramVsFit(bool do7TeV, bool debug=false) {

  double xLLR[10] = {7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5};
  double exLLR[10] = {2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5};
  double yLLR[10] = {-0.00513965, -0.00166552, -0.00053251, -0.00332243, -0.000464651, -0.000663147, 0.00192108, 0.00132186, -0.000359634, 0.00171175};
  double eyLLR[10] = { 0.00654196,0.00457391,0.00300174,0.00251641,0.000995304,0.000697869,0.00272421,0.00160168,0.00173356,0.00260491};

  stringstream filesyst;
  if(do7TeV) filesyst << "gScaleSyst-7TeV.root";
  else filesyst << "gScaleSyst-8TeV.root";

  TFile *tfilesyst = TFile::Open(filesyst.str().c_str());
  TGraphErrors *gScaleInEB = (TGraphErrors*)tfilesyst->Get("gScaleInEB");
  TGraphErrors *gScaleOutEB = (TGraphErrors*)tfilesyst->Get("gScaleOutEB");
  TGraphErrors *gScaleEE = (TGraphErrors*)tfilesyst->Get("gScaleEE");

  /*
  double xCal[10];
  double yCal[10],yCalEB[10],yCalEE[10];
  double eyCal[10],eyCalEB[10],eyCalEE[10];
  

  for(int i=0;i<10;++i) {
    
    double sEB1 = getShift(0.1,xLLR[i],gScaleInEB,gScaleOutEB,gScaleEE);
    double sEB2 = getShift(1.0,xLLR[i],gScaleInEB,gScaleOutEB,gScaleEE);
    double sEE = getShift(2.0,xLLR[i],gScaleInEB,gScaleOutEB,gScaleEE);

    double seEB1 = getShift(0.1,xLLR[i],gScaleInEB,gScaleOutEB,gScaleEE,true);
    double seEB2 = getShift(1.0,xLLR[i],gScaleInEB,gScaleOutEB,gScaleEE,true);
    double seEE = getShift(2.0,xLLR[i],gScaleInEB,gScaleOutEB,gScaleEE,true);

    xCal[i] = xLLR[i]+0.5;
    //    yCal[i] = (sEB1/(seEB1*seEB1) + sEB2/(seEB2*seEB2) + sEE/(seEE*seEE))/(1/(seEB1*seEB1) + 1/(seEB2*seEB2) + 1/(seEE*seEE));
    
    yCalEB[i] = (sEB1 + sEB2)/2.;
    yCalEE[i] = sEE;
    
    yCal[i] = (sEB1 + sEB2 + sEE)/3.;
    eyCal[i] = 1/(1/(seEB1*seEB1) + 1/(seEB2*seEB2) + 1/(seEE*seEE));

    if(eyCal[i]<1e-5) eyCal[i]=0.001 + 0.0002*gRandom->Uniform();
    if(xCal[i]<15)  eyCal[i]*=2.;
    if(xCal[i]>50) eyCal[i]*=1.5;

    eyCalEE[i] = eyCal[i];
    eyCalEB[i] = 0.33*eyCalEE[i];

  }
  */

  double xCal[12],xCalEB[12],xCalEE[12];
  double exCal[12];
  double yCal[12],yCalEB[12],yCalEE[12];
  double eyCal[12],eyCalEB[12],eyCalEE[12];

  
  for(int i=0;i<12;++i) {
    
    double sEB1 = getShift(0.1,ptbins[i],gScaleInEB,gScaleOutEB,gScaleEE);
    double sEB2 = getShift(1.0,ptbins[i],gScaleInEB,gScaleOutEB,gScaleEE);
    double sEE = getShift(2.0,ptbins[i],gScaleInEB,gScaleOutEB,gScaleEE);

    double seEB1 = getShift(0.1,ptbins[i],gScaleInEB,gScaleOutEB,gScaleEE,true);
    double seEB2 = getShift(1.0,ptbins[i],gScaleInEB,gScaleOutEB,gScaleEE,true);
    double seEE = getShift(2.0,ptbins[i],gScaleInEB,gScaleOutEB,gScaleEE,true);

    xCal[i] = ptbins[i];
    xCalEB[i] = ptbins[i]+0.5;
    xCalEE[i] = ptbins[i]-0.5;
    //    yCal[i] = (sEB1/(seEB1*seEB1) + sEB2/(seEB2*seEB2) + sEE/(seEE*seEE))/(1/(seEB1*seEB1) + 1/(seEB2*seEB2) + 1/(seEE*seEE));
    
    exCal[i] = 0.;

    yCalEB[i] = (sEB1 + sEB2)/2.;
    yCalEE[i] = sEE;
    
    yCal[i] = (sEB1 + sEB2 + sEE)/3.;
    eyCal[i] = 1/(1/(seEB1*seEB1) + 1/(seEB2*seEB2) + 1/(seEE*seEE));

    if(eyCal[i]<1e-5) eyCal[i]=0.001 + 0.0002*gRandom->Uniform();
    if(xCal[i]<15)  eyCal[i]*=2.;
    if(xCal[i]>50) eyCal[i]*=1.5;

    eyCalEE[i] = eyCal[i];
    eyCalEB[i] = 0.33*eyCalEE[i];

  }

  TGraphErrors *gLLR = new TGraphErrors(10,xLLR,yLLR,exCal,eyLLR);
  TGraphErrors *gCal = new TGraphErrors(12,xCal,yCal,exCal,eyCal);
  TGraphErrors *gCalEB = new TGraphErrors(12,xCalEB,yCalEB,exCal,eyCalEB);
  TGraphErrors *gCalEE = new TGraphErrors(12,xCalEE,yCalEE,exCal,eyCalEE);
  TGraphErrors *gJPsi = new TGraphErrors(1);
  gJPsi->SetPoint(0,7,0.0013);
  gJPsi->SetPointError(0,0,0.0057);

  gLLR->SetMarkerColor(kGreen+2);
  gCal->SetMarkerColor(kBlue-2);
  gLLR->SetLineColor(kGreen+2);
  gCal->SetLineColor(kBlue-2);
  gJPsi->SetMarkerColor(kRed+2);
  gJPsi->SetLineColor(kRed+2);  

  // only for debug
  gCalEB->SetMarkerColor(kAzure+7);
  gCalEE->SetMarkerColor(kMagenta+3);
  gCalEB->SetLineColor(kAzure+7);
  gCalEE->SetLineColor(kMagenta+3);
  gCalEB->SetMarkerStyle(21);
  gCalEE->SetMarkerStyle(21);

  gLLR->SetMarkerStyle(20);
  gCal->SetMarkerStyle(21);
  gJPsi->SetMarkerStyle(26);

  gLLR->SetName("");
  gLLR->SetTitle("");
  gCalEB->SetName("");
  gCalEB->SetTitle("");
  gCalEB->GetXaxis()->SetTitle("average electron p_{T} [GeV]");
  gCalEB->GetYaxis()->SetTitle("#Delta M/M (data - sim.)");
  gCalEB->GetXaxis()->SetRangeUser(0,65);
  gCalEB->GetYaxis()->SetTitleOffset(1.9);
  gCalEB->GetYaxis()->SetRangeUser(-0.015,0.015);

  if(do7TeV) {
    gCal->GetXaxis()->SetTitle("average electron p_{T} [GeV]");
    gCal->GetYaxis()->SetTitle("#Delta M/M (data - sim.)");
    gCal->GetYaxis()->SetTitleOffset(1.9);
    gCal->GetYaxis()->SetRangeUser(-0.015,0.015);
  }

  TLegend* legend = new TLegend(0.15,0.7,0.4,0.8);
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);

  if(debug) {
    legend->AddEntry(gCalEB, "Z, template shift fit, EB","pl");
    legend->AddEntry(gCalEE, "Z, template shift fit, EE","pl");
  }
  
  if(!debug) legend->AddEntry(gCal, "Z, template shift fit","pl");
  if(!do7TeV) {
    legend->AddEntry(gLLR, "Z, parametric fit, EB+EE","pl");
    legend->AddEntry(gJPsi, "J/#Psi","pl");
  }

  TLatex* CP = new TLatex(7.,0.016, do7TeV ? "CMS Preliminary                          #sqrt{s} = 7 TeV, L = 5.05 fb^{-1}" : "CMS Preliminary                          #sqrt{s} = 8 TeV, L = 12.2 fb^{-1}");
  CP->SetTextSize(0.030);

  TLine *line = new TLine(2,0,65,0);
  line->SetLineColor(kRed);
  line->SetLineWidth(1.5);

  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->SetGridx();
  c1->SetGridy();

  if(debug) {
    gCalEB->Draw("ape");
    gCalEE->Draw("pe");
  }

  if(!do7TeV) {
    gLLR->Draw("pe");
     if(!debug) gCal->Draw("pe");
    gJPsi->Draw("pe");
  } else {
    if(!debug) gCal->Draw("ape");
  }

  line->Draw();
  legend->Draw();
  CP->Draw();
  c1->SaveAs("comp.pdf");


}

void compareTwoRun2012() {

  TStyle *hzzstyle = getStyle("ZZ");
  hzzstyle->cd();

  //  TFile *tfilesystABC = TFile::Open("gScaleSyst-8TeV-GoldGold.root");
  TFile *tfilesystABC = TFile::Open("gScaleSyst-8TeV-Moriond.root");
  TGraphErrors *gScaleInEB1 = (TGraphErrors*)tfilesystABC->Get("gScaleInEB");
  TGraphErrors *gScaleOutEB1 = (TGraphErrors*)tfilesystABC->Get("gScaleOutEB");
  TGraphErrors *gScaleEE1 = (TGraphErrors*)tfilesystABC->Get("gScaleEE");
  vector<TGraphErrors*> g1;
  g1.push_back(gScaleInEB1);
  g1.push_back(gScaleOutEB1);
  g1.push_back(gScaleEE1);

  //  TFile *tfilesystD = TFile::Open("gScaleSyst-8TeV-NGoldNGold.root");
   TFile *tfilesystD = TFile::Open("gScaleSyst-8TeV-22Jan.root");
  TGraphErrors *gScaleInEB2 = (TGraphErrors*)tfilesystD->Get("gScaleInEB");
  TGraphErrors *gScaleOutEB2 = (TGraphErrors*)tfilesystD->Get("gScaleOutEB");
  TGraphErrors *gScaleEE2 = (TGraphErrors*)tfilesystD->Get("gScaleEE");
  vector<TGraphErrors*> g2;
  g2.push_back(gScaleInEB2);
  g2.push_back(gScaleOutEB2);
  g2.push_back(gScaleEE2);
  

  TLegend *legend = new TLegend(0.25,0.75,0.45,0.85,NULL,"brNDC");
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  ( 0.030);

  TLatex* CP = new TLatex(7.,0.0105, "CMS Preliminary                          #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
  CP->SetTextSize(0.030);

  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->SetGridx();
  c1->SetGridy();
  for(int i=0;i<3;++i) {

    g1[i]->GetXaxis()->SetRangeUser(5,60);
    g1[i]->GetXaxis()->SetTitle("p_{T} (GeV)");
    g1[i]->GetYaxis()->SetRangeUser(-0.02,0.01);
    g1[i]->GetYaxis()->SetTitle("#Delta m/m (data - sim.)");
    g1[i]->GetYaxis()->SetTitleOffset(2.0);

    g1[i]->SetLineColor(kBlack);
    g1[i]->SetMarkerColor(kBlack);    
    g2[i]->SetLineColor(kRed+2);
    g2[i]->SetMarkerColor(kRed+2);

    if(i==0) {
      legend->AddEntry(g1[0], "data (Moriond)","pl");
      legend->AddEntry(g2[0], "data (Jan22)","pl");
    }

    g1[i]->Draw("ape");
    g2[i]->Draw("pe1");

    legend->Draw();

    TLine *line = new TLine(5,0,60,0);
    line->SetLineColor(kRed);
    line->SetLineWidth(2.0);
    line->Draw();
    
    CP->Draw();

    stringstream fss1;
    fss1 << "scale-etabin" << i << "_comp.pdf";
    c1->SaveAs(fss1.str().c_str());

    stringstream fss2;
    fss2 << "scale-etabin" << i << "_comp.C";
    c1->SaveAs(fss2.str().c_str());
    
  }

}


void compareEBandEE(bool do7TeV) {

  TStyle *hzzstyle = getStyle("ZZ");
  hzzstyle->cd();

  TFile *tfilesyst = TFile::Open(((do7TeV) ? "gScaleSyst-7TeV.root" : "gScaleSyst-8TeV.root"));
  TGraphErrors *gScaleInEB1 = (TGraphErrors*)tfilesyst->Get("gScaleInEB");
  TGraphErrors *gScaleOutEB1 = (TGraphErrors*)tfilesyst->Get("gScaleOutEB");
  TGraphErrors *gScaleEE1 = (TGraphErrors*)tfilesyst->Get("gScaleEE");
  vector<TGraphErrors*> g1;
  g1.push_back(gScaleInEB1);
  g1.push_back(gScaleOutEB1);
  g1.push_back(gScaleEE1);

  TFile *tfilesystJPsi = TFile::Open(((do7TeV) ? "gScaleSystJPsi-7TeV.root" : "gScaleSystJPsi-8TeV.root"));
  TGraphErrors *gScaleEB2 = (TGraphErrors*)tfilesystJPsi->Get("gScaleEB");
  TGraphErrors *gScaleEE2 = (TGraphErrors*)tfilesystJPsi->Get("gScaleEE");
  vector<TGraphErrors*> g2;
  g2.push_back(gScaleEB2);
  g2.push_back(gScaleEE2);

  TFile *tfilesystUps = TFile::Open(((do7TeV) ? "gScaleSystUpsilon-7TeV.root" : "gScaleSystUpsilon-8TeV.root"));
  TGraphErrors *gScaleEB3 = (TGraphErrors*)tfilesystUps->Get("gScaleEB");
  TGraphErrors *gScaleEE3 = (TGraphErrors*)tfilesystUps->Get("gScaleEE");
  vector<TGraphErrors*> g3;
  g3.push_back(gScaleEB3);
  g3.push_back(gScaleEE3);

  TLegend* legend = new TLegend(0.55,0.2,0.85,0.5);
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetFillStyle (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);

  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->SetGridx();
  c1->SetGridy();

  g1[0]->GetXaxis()->SetRangeUser(5,60);
  g1[0]->GetXaxis()->SetTitle("electron p_{T} (GeV)");
  g1[0]->GetYaxis()->SetRangeUser(-0.02,0.01);
  g1[0]->GetYaxis()->SetTitle("#Delta m/m (data - sim.)");
  g1[0]->GetYaxis()->SetTitleOffset(2.0);
  g1[0]->SetLineColor(kRed+2);
  g1[0]->SetMarkerColor(kRed+2);
  g1[1]->SetLineColor(kOrange+7);
  g1[1]->SetMarkerColor(kOrange+7);
  g1[2]->SetLineColor(kBlue+2);
  g1[2]->SetMarkerColor(kBlue+2);

  g2[0]->SetLineColor(kGreen+2);
  g2[0]->SetMarkerColor(kGreen+2);
  g2[0]->SetMarkerStyle(kFullSquare);
  g2[1]->SetLineColor(kBlue);
  g2[1]->SetMarkerColor(kBlue);
  g2[1]->SetMarkerStyle(kFullSquare);

  g3[0]->SetLineColor(kViolet-6);
  g3[0]->SetMarkerColor(kViolet-6);
  g3[0]->SetMarkerStyle(kFullTriangleUp);
  g3[1]->SetLineColor(kViolet-2);
  g3[1]->SetMarkerColor(kViolet-2);
  g3[1]->SetMarkerStyle(kFullTriangleUp);

  g1[0]->Draw("ape");
  g1[1]->Draw("pe1");
  g1[2]->Draw("pe1");
  g2[0]->Draw("pe1");
  //  g2[1]->Draw("pe1");
  g3[0]->Draw("pe1");
  //  g3[1]->Draw("pe1");

  legend->AddEntry(g1[0], "Z, |#eta|<0.8","pl");
  legend->AddEntry(g1[1], "Z, 0.8<|#eta|<1.48","pl");
  legend->AddEntry(g1[2], "Z, |#eta|>1.48","pl");
  legend->AddEntry(g2[0], "J/#Psi, |#eta|<1.48","pl");
  //legend->AddEntry(g2[1], "J/#Psi, EE","pl");
  legend->AddEntry(g3[0], "#Upsilon (1S), |#eta|<1.48","pl");
  //  legend->AddEntry(g3[1], "#Upsilon (1S), |#eta|>1.48","pl");
  legend->Draw();

  TLatex* CP = new TLatex(7.,0.0105, do7TeV ? "CMS Preliminary                          #sqrt{s} = 7 TeV, L = 5.05 fb^{-1}" : "CMS Preliminary                          #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
  CP->SetTextSize(0.030);

  TLine *line = new TLine(5,0,60,0);
  line->SetLineColor(kRed);
  line->SetLineWidth(2.0);
  line->Draw();
  CP->Draw();

  c1->SaveAs(((do7TeV) ? "scale-ptdep-7TeV.pdf" : "scale-ptdep-8TeV.pdf"));

}


// ******************************************************
// Read data
// ******************************************************
void fillHistogramsJPsi(bool do7TeV, bool isMC) {

  string fileMC;
  string fileData;
  if(do7TeV) {
    fileMC = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim/zEE_lineshape_MC_42X.root");
    fileData = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim/zEE_lineshape_Data_2011.root");
  } else {
    //    fileMC = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim/zEE_lineshape_MC_53X.root");
    //    fileData = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim/zEE_lineshape_Data_2012All.root");
    fileMC = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim_jan22/jpsi_lineshape_2012.root");
    fileData = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim_jan22/zEE_lineshape_Data_2012.root");
  }

  TFile *tfile = 0;
  if(isMC) tfile = TFile::Open(fileMC.c_str());
  else tfile = TFile::Open(fileData.c_str());

  TTree *tree = (TTree*)tfile->Get("zeetree/probe_tree");
  int n=tree->GetEntries();

   Float_t         l1bdtID;
   Float_t         l1bdtIso;
   Float_t         l1classification;
   Float_t         l1p;
   Float_t         l1pdgId;
   Float_t         l1phi;
   Float_t         l1pt;
   Float_t         l1eta;
   Float_t         l1r9;
   Float_t         l2bdtID;
   Float_t         l2bdtIso;
   Float_t         l2classification;
   Float_t         l2p;
   Float_t         l2pdgId;
   Float_t         l2phi;
   Float_t         l2pt;
   Float_t         l2eta;
   Float_t         l2r9;
   Float_t         massErr;
   Float_t         numTrueInteractions;
   Float_t         nvtx;
   Float_t         rhoAA;
   Float_t         zeta;
   Float_t         zmass;
   Float_t         zmll;
   Float_t         zphi;
   Float_t         zpt;
   UInt_t          run;
   UInt_t          lumi;
   UInt_t          event;

   tree->SetBranchAddress("l1bdtID", &l1bdtID);
   tree->SetBranchAddress("l1bdtIso", &l1bdtIso);
   tree->SetBranchAddress("l1classification", &l1classification);
   tree->SetBranchAddress("l1eta", &l1eta);
   tree->SetBranchAddress("l1p", &l1p);
   tree->SetBranchAddress("l1pdgId", &l1pdgId);
   tree->SetBranchAddress("l1phi", &l1phi);
   tree->SetBranchAddress("l1pt", &l1pt);
   tree->SetBranchAddress("l1r9", &l1r9);
   tree->SetBranchAddress("l2bdtID", &l2bdtID);
   tree->SetBranchAddress("l2bdtIso", &l2bdtIso);
   tree->SetBranchAddress("l2classification", &l2classification);
   tree->SetBranchAddress("l2eta", &l2eta);
   tree->SetBranchAddress("l2p", &l2p);
   tree->SetBranchAddress("l2pdgId", &l2pdgId);
   tree->SetBranchAddress("l2phi", &l2phi);
   tree->SetBranchAddress("l2pt", &l2pt);
   tree->SetBranchAddress("l2r9", &l2r9);
   tree->SetBranchAddress("massErr", &massErr);
   tree->SetBranchAddress("nvtx", &nvtx);
   tree->SetBranchAddress("rhoAA", &rhoAA);
   tree->SetBranchAddress("zeta", &zeta);
   tree->SetBranchAddress("zmass", &zmass);
   tree->SetBranchAddress("zmll", &zmll);
   tree->SetBranchAddress("zphi", &zphi);
   tree->SetBranchAddress("zpt", &zpt);
   tree->SetBranchAddress("run", &run);
   tree->SetBranchAddress("lumi", &lumi);
   tree->SetBranchAddress("event", &event);
   if(isMC) tree->SetBranchAddress("numTrueInteractions", &numTrueInteractions);

   TH1F *massH = new TH1F("mass","",50,0.5,5.0);
   massH->SetMarkerColor(kBlack);
   massH->SetLineColor(kBlack);
   massH->SetMarkerStyle(8);
   massH->SetMarkerSize(1);

   TH1F *massHMC = new TH1F("massMC","",101,2.6,3.6);
   massHMC->SetLineColor(kBlack);

   h_Data_Mass_EB.resize(nptbinsJ);
   h_Data_Mass_EE.resize(nptbinsJ);

   h_MC_Mass_EB.resize(nptbinsJ);
   h_MC_Mass_EE.resize(nptbinsJ);

   //   TFile *histofile = TFile::Open("histos.root",(isMC ? "recreate" : "update"));
   
   for(int i=0;i<nptbinsJ;++i) {
     stringstream hssEB, hssEE;
     hssEB << "mass_ptbin" << i << "_EB_";
     hssEE << "mass_ptbin" << i << "_EE_";
     
     hssEB << ((isMC) ? "MC" : "Data");
     hssEE << ((isMC) ? "MC" : "Data");

     if(isMC) {
       h_MC_Mass_EB[i] = (TH1F*)massHMC->Clone(hssEB.str().c_str());
       h_MC_Mass_EE[i] = (TH1F*)massHMC->Clone(hssEE.str().c_str());
     } else {
       h_Data_Mass_EB[i] = (TH1F*)massH->Clone(hssEB.str().c_str());
       h_Data_Mass_EE[i] = (TH1F*)massH->Clone(hssEE.str().c_str());
     }
   }

   for(int jentry=0; jentry<n; ++jentry) {
     tree->GetEntry(jentry);

     if(jentry%100000==0) cout << "analyzing event # " << jentry << endl; 

     if(l1pt<7 || fabs(l1eta)>2.5 ||
	l2pt<7 || fabs(l2eta)>2.5) continue;

     // fill one leg inclusive on the other one

     // if the lepton is not matched, the momentum is 0
     // if(isMC && (tag_gen_p==0 || gen_p==0)) continue;

     // pt bin
     int tagbin=-1; int probebin=-1;
     for(int b=0;b<nptbinsJ;++b) {
       if(l2pt>=ptbinsJ[b] && l2pt<ptbinsJ[b+1]) probebin=b;
       if(l1pt>=ptbinsJ[b] && l1pt<ptbinsJ[b+1]) tagbin=b;
     }
     if(probebin==-1 || tagbin==-1) continue;

     // eta bin
     int tagetabin=-1; int probeetabin=-1;
     for(int b=0;b<3;++b) {
       if(fabs(l2eta)>=etabinsJ[b] && fabs(l2eta)<etabinsJ[b+1]) probeetabin=b;
       if(fabs(l1eta)>=etabinsJ[b] && fabs(l1eta)<etabinsJ[b+1]) tagetabin=b;
     }
     if(probeetabin==-1 || tagetabin==-1) continue;

     switch (probeetabin) {
     case 0:
       if(isMC) h_MC_Mass_EB[probebin]->Fill(zmass);
       else h_Data_Mass_EB[probebin]->Fill(zmass);
       break;
     default:
       if(isMC) h_MC_Mass_EE[probebin]->Fill(zmass);
       else h_Data_Mass_EE[probebin]->Fill(zmass);
       break;
     }
     
     switch (tagetabin) {
     case 0:
       if(isMC) h_MC_Mass_EB[tagbin]->Fill(zmass);
       else h_Data_Mass_EB[tagbin]->Fill(zmass);
       break;
     default:
       if(isMC) h_MC_Mass_EE[tagbin]->Fill(zmass);
       else h_Data_Mass_EE[tagbin]->Fill(zmass);
       break;
     }
     
   } // loop on events

}


void myFitScaleJPsi(bool do7TeV) {

  // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetOptStat(0);
  // ------------------------------ 
  
  fillHistogramsJPsi(do7TeV,true);  
  fillHistogramsJPsi(do7TeV,false);  

  TGraphErrors *gScaleEB = new TGraphErrors();
  TGraphErrors *gScaleEE = new TGraphErrors();

  std::vector <TF1*> f_Mass_EB_MC, f_Mass_EE_MC; 
  std::vector <TF1*> f_Mass_EB_DATA, f_Mass_EE_DATA; 
  std::vector<double> val_Mass_EB, val_Mass_EE; 
  std::vector<double> err_Mass_EB, err_Mass_EE; 

  TCanvas c1("c1","",600,600);

  //float biasMC = -0.007; // from Claude fit to QCD bc->e

  // loop on histograms and fit
  for (unsigned int i=0;i<h_Data_Mass_EB.size();i++) {
    
    cout << h_Data_Mass_EB[i]->GetName() << endl;

    cout << "======= FITTING NOW DATA =======" << endl;
    char funcNameEB_DATA[200], funcNameEE_DATA[200];
    sprintf(funcNameEB_DATA,"f_Mass_EB_DATA_%d",i);
    f_Mass_EB_DATA.push_back( new TF1(funcNameEB_DATA, fitFunction, 1.0, 4.0, 7) );
    f_Mass_EB_DATA.back() -> SetParNames("p1","p2","N","m","#sigma","#alpha","n"); 
    f_Mass_EB_DATA.back() -> SetLineWidth(2); 
    f_Mass_EB_DATA.back() -> SetLineColor(kRed+2); 

    sprintf(funcNameEE_DATA,"f_Mass_EE_DATA_%d",i);
    f_Mass_EE_DATA.push_back( new TF1(funcNameEE_DATA, fitFunction, 1.0, 4.0, 7) );
    f_Mass_EE_DATA.back() -> SetParNames("p1","p2","N","m","#sigma","#alpha","n"); 
    f_Mass_EE_DATA.back() -> SetLineWidth(2); 
    f_Mass_EE_DATA.back() -> SetLineColor(kRed+2); 

    f_Mass_EB_DATA.back() -> SetParLimits(2, 0, 10000);
    f_Mass_EB_DATA.back() -> SetParameter(3, 3.1);
    f_Mass_EB_DATA.back() -> SetParLimits(3, 2.5, 3.5);
    f_Mass_EB_DATA.back() -> SetParameter(4, 0.2 );
    f_Mass_EB_DATA.back() -> SetParLimits(4,  0.02, 0.5 );
    f_Mass_EB_DATA.back() -> FixParameter(5,  5 );
    f_Mass_EB_DATA.back() -> FixParameter(6,  5 );

    f_Mass_EE_DATA.back() -> SetParLimits(2, 0, 10000);
    f_Mass_EE_DATA.back() -> SetParameter(3, 3.1);
    f_Mass_EE_DATA.back() -> SetParLimits(3, 2.5, 3.5);
    f_Mass_EE_DATA.back() -> SetParameter(4, 0.2 );
    f_Mass_EE_DATA.back() -> SetParLimits(4,  0.08, 0.5 );
    f_Mass_EE_DATA.back() -> FixParameter(5,  5 );
    f_Mass_EE_DATA.back() -> FixParameter(6,  5 );

    float pTmean = ptbinsJ[i];

    TFitResultPtr rp1 = h_Data_Mass_EB[i] -> Fit(funcNameEB_DATA, "EHRS");
    TFitResultPtr rp2 = h_Data_Mass_EE[i] -> Fit(funcNameEE_DATA, "EHRS");

    float MassscaleEB_DATA = f_Mass_EB_DATA.back()->GetParameter(3);
    float MasserrorEB_DATA = f_Mass_EB_DATA.back()->GetParError(3);
    float MassscaleEE_DATA = f_Mass_EE_DATA.back()->GetParameter(3);
    float MasserrorEE_DATA = f_Mass_EE_DATA.back()->GetParError(3);
    cout << "======= DATA FITTED =======" << endl;



    cout << "======= FITTING NOW MC =======" << endl;
    char funcNameEB_MC[200], funcNameEE_MC[200];
    sprintf(funcNameEB_MC,"f_Mass_EB_MC_%d",i);
    f_Mass_EB_MC.push_back( new TF1(funcNameEB_MC, fitFunction, 2.95, 3.3, 7) );
    f_Mass_EB_MC.back() -> SetParNames("p1","p2","N","m","#sigma","#alpha","n"); 
    f_Mass_EB_MC.back() -> SetLineWidth(2); 
    f_Mass_EB_MC.back() -> SetLineColor(kGreen+2); 

    sprintf(funcNameEE_MC,"f_Mass_EE_MC_%d",i);
    f_Mass_EE_MC.push_back( new TF1(funcNameEE_MC, fitFunction, 2.6, 3.6, 7) );
    f_Mass_EE_MC.back() -> SetParNames("p1","p2","N","m","#sigma","#alpha","n"); 
    f_Mass_EE_MC.back() -> SetLineWidth(2); 
    f_Mass_EE_MC.back() -> SetLineColor(kGreen+2); 

    f_Mass_EB_MC.back() -> SetParLimits(2, 0, 10000);
    f_Mass_EB_MC.back() -> SetParameter(3, 3.1);
    f_Mass_EB_MC.back() -> SetParLimits(3, 2.5, 3.5);
    f_Mass_EB_MC.back() -> SetParameter(4, 0.2 );
    f_Mass_EB_MC.back() -> SetParLimits(4,  0.02, 0.5 );
    f_Mass_EB_MC.back() -> FixParameter(5,  5 );
    f_Mass_EB_MC.back() -> FixParameter(6,  5 );

    f_Mass_EE_MC.back() -> SetParLimits(2, 0, 10000);
    f_Mass_EE_MC.back() -> SetParameter(3, 3.1);
    f_Mass_EE_MC.back() -> SetParLimits(3, 2.5, 3.5);
    f_Mass_EE_MC.back() -> SetParameter(4, 0.2 );
    f_Mass_EE_MC.back() -> SetParLimits(4,  0.08, 0.5 );
    f_Mass_EE_MC.back() -> FixParameter(5,  5 );
    f_Mass_EE_MC.back() -> FixParameter(6,  5 );

    rp1 = h_MC_Mass_EB[i] -> Fit(funcNameEB_MC, "EHRS");
    rp2 = h_MC_Mass_EE[i] -> Fit(funcNameEE_MC, "EHRS");

    float MassscaleEB_MC = f_Mass_EB_MC.back()->GetParameter(3);
    float MasserrorEB_MC = f_Mass_EB_MC.back()->GetParError(3);
    float MassscaleEE_MC = f_Mass_EE_MC.back()->GetParameter(3);
    float MasserrorEE_MC = f_Mass_EE_MC.back()->GetParError(3);
    cout << "======= MC FITTED =======" << endl;


    val_Mass_EB.push_back((MassscaleEB_DATA-MassscaleEB_MC)/MassscaleEB_DATA); err_Mass_EB.push_back(sqrt(MasserrorEB_DATA*MasserrorEB_DATA+MasserrorEB_MC*MasserrorEB_MC)); 
    val_Mass_EE.push_back((MassscaleEE_DATA-MassscaleEB_MC)/MassscaleEE_DATA); err_Mass_EE.push_back(sqrt(MasserrorEE_DATA*MasserrorEE_DATA+MasserrorEE_MC*MasserrorEE_MC)); 

    // derive scale and fill graph
    gScaleEB->SetPoint(i,pTmean,(MassscaleEB_DATA-MassscaleEB_MC)/MassscaleEB_DATA); 
    gScaleEB->SetPointError(i,(ptbinsJ[2]-ptbinsJ[1])/2.,MasserrorEB_DATA/3.097); 

    gScaleEE->SetPoint(i,pTmean,(MassscaleEE_DATA-MassscaleEE_MC)/3.097); 
    gScaleEE->SetPointError(i,(ptbinsJ[2]-ptbinsJ[1])/2.,MasserrorEE_DATA/3.097); 

  }
  

  TLegend* legend = new TLegend(0.20, 0.40, 0.30, 0.55);
  
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetFillStyle (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  (0.05);
  legend->AddEntry(h_Data_Mass_EB[0],"data");
  legend->AddEntry(f_Mass_EB_DATA[0],"fit");

  TLegend* legendMC = new TLegend(0.20, 0.40, 0.30, 0.55);
  
  legendMC->SetBorderSize(     0);
  legendMC->SetFillColor (     0);
  legendMC->SetFillStyle (     0);
  legendMC->SetTextAlign (    12);
  legendMC->SetTextFont  (    42);
  legendMC->SetTextSize  (0.05);
  legendMC->AddEntry(h_MC_Mass_EB[0],"MC");
  legendMC->AddEntry(f_Mass_EB_MC[0],"fit");

  TPaveStats** s_Mass = new TPaveStats*[500];
  TCanvas *c2 = new TCanvas("c2","",600,600); 
  for(unsigned int i = 0; i < h_Data_Mass_EB.size(); ++i) {

    TLatex* CP1 = new TLatex(1.0,h_Data_Mass_EB[i]->GetMaximum()*1.2, do7TeV ? "CMS Preliminary                          #sqrt{s} = 7 TeV, L = 5.05 fb^{-1}" : "CMS Preliminary                          #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
    CP1->SetTextSize(0.030);

    char canvasName[50];
    
    int nTeV = (do7TeV) ? 7 : 8;

    // EB
    TPaveText *bin0 = new TPaveText(0.20,0.75,0.45,0.85,"brNDC");
    stringstream binval0;
    if(i!=h_Data_Mass_EB.size()-1) binval0 << ptbinsJ[i] << " < p_{T} < " << ptbinsJ[i+1] << " GeV"; 
    else binval0 << "p_{T} > " << ptbinsJ[i] << "GeV";
    bin0->AddText(binval0.str().c_str());
    bin0->AddText("|#eta|<1.479");
    bin0->SetBorderSize(0);
    bin0->SetFillStyle(0);
    bin0->SetTextAlign(12);
    bin0->SetTextFont(132);
    bin0->SetTextSize(0.04);

    sprintf(canvasName, "FitsJPsi-EB-Pt%d-%dTeV.pdf", i, nTeV); 
    h_Data_Mass_EB[i] -> GetXaxis() -> SetTitle("m_{ee} (GeV)");
    h_Data_Mass_EB[i] -> GetXaxis() -> SetRangeUser(1.0,4.0); 
    h_Data_Mass_EB[i] -> GetYaxis() -> SetTitle("events");
    h_Data_Mass_EB[i] -> GetYaxis() -> SetTitleOffset(1.2);
    h_Data_Mass_EB[i] -> Draw("e");
    gPad->Update(); 
    TPaveText *sigmat0 = new TPaveText(0.30,0.90,0.77,0.95,"brNDC");
    stringstream sigmaval0;
    sigmaval0 << fixed;
    sigmaval0 << setprecision(3);
    sigmaval0 << "#Delta M/M (data-sim.) " << val_Mass_EB[i] << " #pm " << err_Mass_EB[i];
    sigmat0->AddText(sigmaval0.str().c_str());
    sigmat0->SetBorderSize(0);
    sigmat0->SetFillStyle(0);
    sigmat0->SetTextAlign(12);
    sigmat0->SetTextFont(132);
    sigmat0->SetTextSize(0.04);

    // s_Mass[i]= (TPaveStats*)(h_Data_Mass_EB[i]->GetListOfFunctions()->FindObject("stats"));
    // s_Mass[i]->SetTextColor(kRed+2);
    // s_Mass[i]->SetOptStat(0);
    // s_Mass[i]->SetOptFit(1111);
    // s_Mass[i]->SetY1NDC(0.55);
    // s_Mass[i]->SetY2NDC(0.85);
    // s_Mass[i]->SetX1NDC(0.15); 
    // s_Mass[i]->SetX2NDC(0.45); 
    // f_Mass_EB[i]->Draw("same");
    TF1 *backFcnEB = new TF1("backFcnEB",background,1.0,4.0,2);
    backFcnEB->SetLineColor(kBlue+2);
    backFcnEB->SetLineStyle(kDashed);
    Double_t par[7];
    f_Mass_EB_DATA[i]->GetParameters(par);
    backFcnEB->SetParameters(par);
    backFcnEB->Draw("same");
    legend->Draw();
    sigmat0->Draw();
    bin0->Draw();
    //CP1->Draw();
    c2->SaveAs(canvasName);

    // now the MC
    h_MC_Mass_EB[i] -> Draw();
    legendMC->Draw();
    sigmat0->Draw();
    bin0->Draw();
    //CP1->Draw();
    sprintf(canvasName, "FitsJPsi-EB-Pt%d-%dTeV_MC.pdf", i, nTeV); 
    c2->SaveAs(canvasName);


    // EE
    TPaveText *bin1 = new TPaveText(0.20,0.75,0.45,0.85,"brNDC");
    stringstream binval1;
    if(i!=h_Data_Mass_EE.size()-1) binval1 << ptbinsJ[i] << " < p_{T} < " << ptbinsJ[i+1] << " GeV"; 
    else binval1 << "p_{T} > " << ptbinsJ[i] << "GeV";
    bin1->AddText(binval0.str().c_str());
    bin1->AddText("|#eta|>1.479");
    bin1->SetBorderSize(0);
    bin1->SetFillStyle(0);
    bin1->SetTextAlign(12);
    bin1->SetTextFont(132);
    bin1->SetTextSize(0.04);

    TLatex* CP2 = new TLatex(1.,h_Data_Mass_EE[i]->GetMaximum()*1.5, do7TeV ? "CMS Preliminary                          #sqrt{s} = 7 TeV, L = 5.05 fb^{-1}" : "CMS Preliminary                          #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
    CP2->SetTextSize(0.030);
    sprintf(canvasName, "FitsJPsi-EE-Pt%d-%dTeV.pdf", i, nTeV); 
    h_Data_Mass_EE[i] -> GetXaxis() -> SetTitle("m_{ee} (GeV)");
    h_Data_Mass_EE[i] -> GetXaxis() -> SetRangeUser(1.0,4.0); 
    h_Data_Mass_EE[i] -> GetYaxis() -> SetTitle("events");
    h_Data_Mass_EE[i] -> GetYaxis() -> SetTitleOffset(1.2);
    h_Data_Mass_EE[i] -> Draw("e");
    gPad->Update(); 
    TPaveText *sigmat1 = new TPaveText(0.30,0.90,0.77,0.95,"brNDC");
    stringstream sigmaval1;
    sigmaval1 << fixed;
    sigmaval1 << setprecision(3);
    sigmaval1 << "#Delta M/M (data-sim.) " << val_Mass_EE[i] << " #pm " << err_Mass_EE[i];
    sigmat1->AddText(sigmaval1.str().c_str());
    sigmat1->SetBorderSize(0);
    sigmat1->SetFillStyle(0);
    sigmat1->SetTextAlign(12);
    sigmat1->SetTextFont(132);
    sigmat1->SetTextSize(0.04);

    // s_Mass[i]= (TPaveStats*)(h_Data_Mass_EE[i]->GetListOfFunctions()->FindObject("stats"));
    // s_Mass[i]->SetTextColor(kRed+2);
    // s_Mass[i]->SetOptStat(0);
    // s_Mass[i]->SetOptFit(1111);
    // s_Mass[i]->SetY1NDC(0.55);
    // s_Mass[i]->SetY2NDC(0.85);
    // s_Mass[i]->SetX1NDC(0.15); 
    // s_Mass[i]->SetX2NDC(0.45); 
    //  f_Mass_EE_DATA[i]->Draw("same");
    TF1 *backFcnEE = new TF1("backFcnEE",background,1.0,4.0,2);
    backFcnEE->SetLineColor(kBlue+2);
    backFcnEE->SetLineStyle(kDashed);
    f_Mass_EE_DATA[i]->GetParameters(par);
    backFcnEE->SetParameters(par);
    backFcnEE->Draw("same");
    legend->Draw();
    sigmat1->Draw();
    bin1->Draw();    
    //CP2->Draw();
    c2->SaveAs(canvasName);

    // now the MC
    h_MC_Mass_EE[i] -> Draw();
    legendMC->Draw();
    sigmat0->Draw();
    bin0->Draw();
    //CP1->Draw();
    sprintf(canvasName, "FitsJPsi-EE-Pt%d-%dTeV_MC.pdf", i, nTeV); 
    c2->SaveAs(canvasName);

  }

  cout << "DONE" << endl;

  TCanvas *c = new TCanvas("c","c");
  c->cd();
  gScaleEB->GetYaxis()->SetRangeUser(-0.015,0.015);
  gScaleEB->GetYaxis()->SetTitle("(M_{ee}/m_{Z})^{data}-(M_{ee}/m_{Z})^{MC}");
  gScaleEB->GetYaxis()->SetTitleOffset(1.5);
  gScaleEB->GetXaxis()->SetTitle("electron p_{T} [GeV]");
  //  gScaleEB->SetMarkerColor(kRed+2);
  gScaleEB->SetMarkerStyle(20);
  gScaleEB->Draw("AP");
  c->SaveAs(((do7TeV) ? "scaleEBJPsi-7TeV.pdf" : "scaleEBJPsi-8TeV.pdf"));

  gScaleEE->GetYaxis()->SetRangeUser(-0.03,0.015);
  gScaleEE->GetYaxis()->SetTitle("(M_{ee}/m_{Z})^{data}-(M_{ee}/m_{Z})^{MC}");
  gScaleEE->GetYaxis()->SetTitleOffset(1.5);
  gScaleEE->GetXaxis()->SetTitle("electron p_{T} [GeV]");
  //  gScaleEE->SetMarkerColor(kRed+2);
  gScaleEE->SetMarkerStyle(20);
  gScaleEE->Draw("AP");
  c->SaveAs(((do7TeV) ? "scaleEEJPsi-7TeV.pdf" : "scaleEEJPsi-8TeV.pdf"));

  TFile *ff =new TFile(((do7TeV) ? "gScaleSystJPsi-7TeV.root" : "gScaleSystJPsi-8TeV.root"),"RECREATE");
  ff->cd();
  gScaleEB->SetName("gScaleEB");
  gScaleEB->Write();
  gScaleEE->SetName("gScaleEE");
  gScaleEE->Write();
}


// ******************************************************
// Read data
// ******************************************************
void fillHistogramsUpsilon(bool do7TeV, bool isMC) {

  string fileMC;
  string fileData;
  if(do7TeV) {
    fileMC = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim/zEE_lineshape_MC_42X.root");
    fileData = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim/zEE_lineshape_Data_2011.root");
  } else {
    fileMC = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim/zEE_lineshape_MC_53X.root");
    fileData = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim/zEE_lineshape_Data_2012All.root");
    // fileMC = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim_jan22/zEE_lineshape_MC_53X.root");
    // fileData = string("/Users/emanuele/Work/data/hzz4l/elescale_2lepskim_jan22/zEE_lineshape_Data_2012.root");
  }

  TFile *tfile = 0;
  if(isMC) tfile = TFile::Open(fileMC.c_str());
  else tfile = TFile::Open(fileData.c_str());

  TTree *tree = (TTree*)tfile->Get("zeetree/probe_tree");
  int n=tree->GetEntries();

   Float_t         l1bdtID;
   Float_t         l1bdtIso;
   Float_t         l1classification;
   Float_t         l1p;
   Float_t         l1pdgId;
   Float_t         l1phi;
   Float_t         l1pt;
   Float_t         l1eta;
   Float_t         l1r9;
   Float_t         l2bdtID;
   Float_t         l2bdtIso;
   Float_t         l2classification;
   Float_t         l2p;
   Float_t         l2pdgId;
   Float_t         l2phi;
   Float_t         l2pt;
   Float_t         l2eta;
   Float_t         l2r9;
   Float_t         massErr;
   Float_t         numTrueInteractions;
   Float_t         nvtx;
   Float_t         rhoAA;
   Float_t         zeta;
   Float_t         zmass;
   Float_t         zmll;
   Float_t         zphi;
   Float_t         zpt;
   UInt_t          run;
   UInt_t          lumi;
   UInt_t          event;

   tree->SetBranchAddress("l1bdtID", &l1bdtID);
   tree->SetBranchAddress("l1bdtIso", &l1bdtIso);
   tree->SetBranchAddress("l1classification", &l1classification);
   tree->SetBranchAddress("l1eta", &l1eta);
   tree->SetBranchAddress("l1p", &l1p);
   tree->SetBranchAddress("l1pdgId", &l1pdgId);
   tree->SetBranchAddress("l1phi", &l1phi);
   tree->SetBranchAddress("l1pt", &l1pt);
   tree->SetBranchAddress("l1r9", &l1r9);
   tree->SetBranchAddress("l2bdtID", &l2bdtID);
   tree->SetBranchAddress("l2bdtIso", &l2bdtIso);
   tree->SetBranchAddress("l2classification", &l2classification);
   tree->SetBranchAddress("l2eta", &l2eta);
   tree->SetBranchAddress("l2p", &l2p);
   tree->SetBranchAddress("l2pdgId", &l2pdgId);
   tree->SetBranchAddress("l2phi", &l2phi);
   tree->SetBranchAddress("l2pt", &l2pt);
   tree->SetBranchAddress("l2r9", &l2r9);
   tree->SetBranchAddress("massErr", &massErr);
   tree->SetBranchAddress("nvtx", &nvtx);
   tree->SetBranchAddress("rhoAA", &rhoAA);
   tree->SetBranchAddress("zeta", &zeta);
   tree->SetBranchAddress("zmass", &zmass);
   tree->SetBranchAddress("zmll", &zmll);
   tree->SetBranchAddress("zphi", &zphi);
   tree->SetBranchAddress("zpt", &zpt);
   tree->SetBranchAddress("run", &run);
   tree->SetBranchAddress("lumi", &lumi);
   tree->SetBranchAddress("event", &event);
   if(isMC) tree->SetBranchAddress("numTrueInteractions", &numTrueInteractions);

   TH1F *massH = new TH1F("mass","",45,8.0,11.5);
   massH->SetMarkerColor(kBlack);
   massH->SetLineColor(kBlack);
   massH->SetMarkerStyle(8);
   massH->SetMarkerSize(1);

   h_Data_Mass_EB.resize(nptbinsU);
   h_Data_Mass_EE.resize(nptbinsU);

   h_MC_Mass_EB.resize(nptbinsU);
   h_MC_Mass_EE.resize(nptbinsU);

   //   TFile *histofile = TFile::Open("histos.root",(isMC ? "recreate" : "update"));
   
   for(int i=0;i<nptbinsU;++i) {
     stringstream hssEB, hssEE;
     hssEB << "mass_ptbin" << i << "_EB_";
     hssEE << "mass_ptbin" << i << "_EE_";
     
     hssEB << ((isMC) ? "MC" : "Data");
     hssEE << ((isMC) ? "MC" : "Data");

     if(isMC) {
       h_MC_Mass_EB[i] = (TH1F*)massH->Clone(hssEB.str().c_str());
       h_MC_Mass_EE[i] = (TH1F*)massH->Clone(hssEE.str().c_str());
     } else {
       h_Data_Mass_EB[i] = (TH1F*)massH->Clone(hssEB.str().c_str());
       h_Data_Mass_EE[i] = (TH1F*)massH->Clone(hssEE.str().c_str());
     }
   }

   for(int jentry=0; jentry<n; ++jentry) {
     tree->GetEntry(jentry);

     if(jentry%100000==0) cout << "analyzing event # " << jentry << endl; 

     if(l1pt<7 || fabs(l1eta)>2.5 ||
	l2pt<7 || fabs(l2eta)>2.5) continue;

     // fill one leg inclusive on the other one

     // if the lepton is not matched, the momentum is 0
     // if(isMC && (tag_gen_p==0 || gen_p==0)) continue;

     // pt bin
     int tagbin=-1; int probebin=-1;
     for(int b=0;b<nptbinsU;++b) {
       if(l2pt>=ptbinsU[b] && l2pt<ptbinsU[b+1]) probebin=b;
       if(l1pt>=ptbinsU[b] && l1pt<ptbinsU[b+1]) tagbin=b;
     }
     if(probebin==-1 || tagbin==-1) continue;

     // eta bin
     int tagetabin=-1; int probeetabin=-1;
     for(int b=0;b<3;++b) {
       if(fabs(l2eta)>=etabinsU[b] && fabs(l2eta)<etabinsU[b+1]) probeetabin=b;
       if(fabs(l1eta)>=etabinsU[b] && fabs(l1eta)<etabinsU[b+1]) tagetabin=b;
     }
     if(probeetabin==-1 || tagetabin==-1) continue;

     switch (probeetabin) {
     case 0:
       if(isMC) h_MC_Mass_EB[probebin]->Fill(zmass);
       else h_Data_Mass_EB[probebin]->Fill(zmass);
       break;
     default:
       if(isMC) h_MC_Mass_EE[probebin]->Fill(zmass);
       else h_Data_Mass_EE[probebin]->Fill(zmass);
       break;
     }
     
     switch (tagetabin) {
     case 0:
       if(isMC) h_MC_Mass_EB[tagbin]->Fill(zmass);
       else h_Data_Mass_EB[tagbin]->Fill(zmass);
       break;
     default:
       if(isMC) h_MC_Mass_EE[tagbin]->Fill(zmass);
       else h_Data_Mass_EE[tagbin]->Fill(zmass);
       break;
     }
     
   } // loop on events

}


void myFitScaleUpsilon(bool do7TeV) {

  // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetOptStat(0);
  // ------------------------------ 
  
  //  fillHistograms(do7TeV,true);  
  fillHistogramsUpsilon(do7TeV,false);  

  TGraphErrors *gScaleEB = new TGraphErrors();
  TGraphErrors *gScaleEE = new TGraphErrors();

  std::vector <TF1*> f_Mass_EB, f_Mass_EE; 
  std::vector<double> val_Mass_EB, val_Mass_EE; 
  std::vector<double> err_Mass_EB, err_Mass_EE; 

  TCanvas c1("c1","",600,600);

  //  float biasMC = -0.047; // from Claude fit to QCD bc->e for JPsi
  float biasMC = 0.;

  // loop on histograms and fit
  for (unsigned int i=0;i<h_Data_Mass_EB.size();i++) {
    
    cout << h_Data_Mass_EB[i]->GetName() << endl;

    char funcNameEB[200], funcNameEE[200];
    sprintf(funcNameEB,"f_Mass_EB_%d",i);
    // f_Mass_EB.push_back( new TF1(funcNameEB, fitFunctionUpsilon, 9.23, 11.5, 10) );
    // f_Mass_EB.back() -> SetParNames("N1S","mean #Upsilon (1S)","#sigma #Upsilon (1S)","#alpha","n","N2S","mean #Upsilon (2S+3S)","#sigma #Upsilon (2S+3S)","p1","p2"); 
    f_Mass_EB.push_back( new TF1(funcNameEB, "gaus", 9.23, 9.7) );
    f_Mass_EB.back() -> SetLineWidth(2); 
    f_Mass_EB.back() -> SetLineColor(kRed+2); 

    sprintf(funcNameEE,"f_Mass_EE_%d",i);
    // f_Mass_EE.push_back( new TF1(funcNameEE, fitFunctionUpsilon, 9.23, 11.5, 10) );
    // f_Mass_EE.back() -> SetParNames("N1S","m1S","#sigma 1S","#alpha","n","N2S","mean 2S+3S","#sigma 2S+3S","p1","p2");
    f_Mass_EE.push_back( new TF1(funcNameEE, "gaus", 9.23, 9.7) );
    f_Mass_EE.back() -> SetLineWidth(2); 
    f_Mass_EE.back() -> SetLineWidth(2); 
    f_Mass_EE.back() -> SetLineColor(kRed+2); 

    // f_Mass_EB.back() -> SetParLimits(0, 1, 100000);
    // f_Mass_EB.back() -> SetParLimits(1, 9.3, 9.55);
    // f_Mass_EB.back() -> SetParameter(2, 0.1 );
    // f_Mass_EB.back() -> SetParLimits(2,  0.05, 0.3 );
    // f_Mass_EB.back() -> SetParLimits(3,  0, 50 );
    // f_Mass_EB.back() -> SetParameter(4,  5 );
    // f_Mass_EB.back() -> SetParLimits(5, 1, 100000);
    // f_Mass_EB.back() -> SetParLimits(6, 9.7, 10.4);
    // f_Mass_EB.back() -> SetParLimits(7,  0.1, 0.5 );
    // f_Mass_EB.back() -> SetParLimits(8,  -10, 10 );
    // f_Mass_EB.back() -> SetParLimits(9,  -10, 10 );

    // f_Mass_EE.back() -> SetParLimits(0, 1, 100000);
    // f_Mass_EE.back() -> SetParLimits(1, 9.3, 9.55);
    // f_Mass_EE.back() -> SetParameter(2, 0.1 );
    // f_Mass_EE.back() -> SetParLimits(2,  0.05, 0.3 );
    // f_Mass_EE.back() -> SetParLimits(3,  0, 50 );
    // f_Mass_EE.back() -> FixParameter(4,  5 );
    // f_Mass_EE.back() -> SetParLimits(5, 1, 100000);
    // f_Mass_EE.back() -> SetParLimits(6, 10.0, 10.4);
    // f_Mass_EE.back() -> SetParLimits(7,  0.1, 0.5 );
    // f_Mass_EE.back() -> SetParLimits(8,  -20, 0 );
    // f_Mass_EE.back() -> SetParLimits(9,  -1, 1 );


    // fit data with template
    cout << "======= FITTING NOW =======" << endl;
    cout << "--> in EB " << endl;
    TFitResultPtr rp1 = h_Data_Mass_EB[i] -> Fit(funcNameEB, "EMHRS");
    TFitResultPtr rp2 = h_Data_Mass_EE[i] -> Fit(funcNameEE, "EMHRS");
    cout << "======= FITTED =======" << endl;

    float pTmean = ptbinsU[i];

    float errorMC = 0.010; // from Claude fit to QCD bc->e
    float MassscaleEB = f_Mass_EB.back()->GetParameter(1);
    float MasserrorEB = sqrt(pow(f_Mass_EB.back()->GetParError(1),2) + errorMC*errorMC);
    val_Mass_EB.push_back(MassscaleEB); err_Mass_EB.push_back(MasserrorEB); 
    float MassscaleEE = f_Mass_EE.back()->GetParameter(1);
    float MasserrorEE = sqrt(pow(f_Mass_EE.back()->GetParError(1),2) + errorMC*errorMC);
    val_Mass_EE.push_back(MassscaleEE); err_Mass_EE.push_back(MasserrorEE); 

    // derive scale and fill graph
    gScaleEB->SetPoint(i,pTmean,(MassscaleEB-9.460+biasMC)/9.460); 
    gScaleEB->SetPointError(i,(ptbinsU[2]-ptbinsU[1])/2.,MasserrorEB/9.460); 

    gScaleEE->SetPoint(i,pTmean,(MassscaleEE-9.460+biasMC)/9.460); 
    gScaleEE->SetPointError(i,(ptbinsU[2]-ptbinsU[1])/2.,MasserrorEE/9.460); 

  }
  

  TLegend* legend = new TLegend(0.15, 0.40, 0.30, 0.55);
  
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetFillStyle (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  (0.05);

  legend->AddEntry(h_Data_Mass_EB[0],"data");
  legend->AddEntry(f_Mass_EB[0],"fit","l");

  TPaveStats** s_Mass = new TPaveStats*[500];
  TCanvas *c2 = new TCanvas("c2","",600,600); 
  for(unsigned int i = 0; i < h_Data_Mass_EB.size(); ++i) {

    TLatex* CP1 = new TLatex(0.2,0.95, do7TeV ? "CMS Preliminary                          #sqrt{s} = 7 TeV, L = 5.05 fb^{-1}" : "CMS Preliminary                          #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
    CP1->SetNDC(kTRUE);
    CP1->SetTextSize(0.030);

    char canvasName[50];
    
    int nTeV = (do7TeV) ? 7 : 8;

    // EB
    TPaveText *bin0 = new TPaveText(0.20,0.75,0.45,0.85,"brNDC");
    stringstream binval0;
    if(i!=h_Data_Mass_EB.size()-1) binval0 << ptbinsU[i] << " < p_{T} < " << ptbinsU[i+1] << " GeV"; 
    else binval0 << "p_{T} > " << ptbinsU[i] << "GeV";
    bin0->AddText(binval0.str().c_str());
    bin0->AddText("|#eta|<1.479");
    bin0->SetBorderSize(0);
    bin0->SetFillStyle(0);
    bin0->SetTextAlign(12);
    bin0->SetTextFont(132);
    bin0->SetTextSize(0.04);

    sprintf(canvasName, "FitsUpsilon-EB-Pt%d-%dTeV.pdf", i, nTeV); 
    h_Data_Mass_EB[i] -> GetXaxis() -> SetTitle("m_{ee} (GeV)");
    h_Data_Mass_EB[i] -> GetXaxis() -> SetRangeUser(1.0,4.0); 
    h_Data_Mass_EB[i] -> GetYaxis() -> SetTitle("events");
    h_Data_Mass_EB[i] -> GetYaxis() -> SetTitleOffset(1.2);
    h_Data_Mass_EB[i] -> Draw("e");
    gPad->Update(); 
    TPaveText *sigmat0 = new TPaveText(0.30,0.90,0.77,0.95,"brNDC");
    stringstream sigmaval0;
    sigmaval0 << fixed;
    sigmaval0 << setprecision(3);
    sigmaval0 << "#Delta M/M (data-sim.) " << (val_Mass_EB[i]-(9.460-biasMC))/9.460 << " #pm " << err_Mass_EB[i]/9.460;
    sigmat0->AddText(sigmaval0.str().c_str());
    sigmat0->SetBorderSize(0);
    sigmat0->SetFillStyle(0);
    sigmat0->SetTextAlign(12);
    sigmat0->SetTextFont(132);
    sigmat0->SetTextSize(0.04);

    // s_Mass[i]= (TPaveStats*)(h_Data_Mass_EB[i]->GetListOfFunctions()->FindObject("stats"));
    // s_Mass[i]->SetTextColor(kRed+2);
    // s_Mass[i]->SetOptStat(0);
    // s_Mass[i]->SetOptFit(1111);
    // s_Mass[i]->SetY1NDC(0.55);
    // s_Mass[i]->SetY2NDC(0.85);
    // s_Mass[i]->SetX1NDC(0.65); 
    // s_Mass[i]->SetX2NDC(0.95); 
    f_Mass_EB[i]->Draw("same");
    TF1 *backFcnEB = new TF1("backFcnEB",background,1.0,4.0,2);
    backFcnEB->SetLineColor(kBlue+2);
    backFcnEB->SetLineStyle(kDashed);
    Double_t par[7];
    f_Mass_EB[i]->GetParameters(par);
    backFcnEB->SetParameters(par);
    backFcnEB->Draw("same");
    legend->Draw();
    sigmat0->Draw();
    bin0->Draw();
    //    CP1->Draw();
    c2->SaveAs(canvasName);

    // EE
    TPaveText *bin1 = new TPaveText(0.20,0.75,0.45,0.85,"brNDC");
    stringstream binval1;
    if(i!=h_Data_Mass_EE.size()-1) binval1 << ptbinsU[i] << " < p_{T} < " << ptbinsU[i+1] << " GeV"; 
    else binval1 << "p_{T} > " << ptbinsU[i] << "GeV";
    bin1->AddText(binval0.str().c_str());
    bin1->AddText("|#eta|>1.479");
    bin1->SetBorderSize(0);
    bin1->SetFillStyle(0);
    bin1->SetTextAlign(12);
    bin1->SetTextFont(132);
    bin1->SetTextSize(0.04);

    TLatex* CP2 = new TLatex(0.2,0.95, do7TeV ? "CMS Preliminary                          #sqrt{s} = 7 TeV, L = 5.05 fb^{-1}" : "CMS Preliminary                          #sqrt{s} = 8 TeV, L = 19.6 fb^{-1}");
    CP2->SetNDC(kTRUE);
    CP2->SetTextSize(0.030);
    sprintf(canvasName, "FitsUpsilon-EE-Pt%d-%dTeV.pdf", i, nTeV); 
    h_Data_Mass_EE[i] -> GetXaxis() -> SetTitle("m_{ee} (GeV)");
    h_Data_Mass_EE[i] -> GetXaxis() -> SetRangeUser(1.0,4.0); 
    h_Data_Mass_EE[i] -> GetYaxis() -> SetTitle("events");
    h_Data_Mass_EE[i] -> GetYaxis() -> SetTitleOffset(1.2);
    h_Data_Mass_EE[i] -> Draw("e");
    gPad->Update(); 
    TPaveText *sigmat1 = new TPaveText(0.30,0.90,0.77,0.95,"brNDC");
    stringstream sigmaval1;
    sigmaval1 << fixed;
    sigmaval1 << setprecision(3);
    sigmaval1 << "#Delta M/M (data-sim.) " << (val_Mass_EE[i]-(9.460-biasMC))/9.460 << " #pm " << err_Mass_EE[i]/9.460;
    sigmat1->AddText(sigmaval1.str().c_str());
    sigmat1->SetBorderSize(0);
    sigmat1->SetFillStyle(0);
    sigmat1->SetTextAlign(12);
    sigmat1->SetTextFont(132);
    sigmat1->SetTextSize(0.04);

    // s_Mass[i]= (TPaveStats*)(h_Data_Mass_EE[i]->GetListOfFunctions()->FindObject("stats"));
    // s_Mass[i]->SetTextColor(kRed+2);
    // s_Mass[i]->SetOptStat(0);
    // s_Mass[i]->SetOptFit(1111);
    // s_Mass[i]->SetY1NDC(0.55);
    // s_Mass[i]->SetY2NDC(0.85);
    // s_Mass[i]->SetX1NDC(0.65); 
    // s_Mass[i]->SetX2NDC(0.95); 
    f_Mass_EE[i]->Draw("same");
    TF1 *backFcnEE = new TF1("backFcnEE",background,1.0,4.0,2);
    backFcnEE->SetLineColor(kBlue+2);
    backFcnEE->SetLineStyle(kDashed);
    f_Mass_EE[i]->GetParameters(par);
    backFcnEE->SetParameters(par);
    backFcnEE->Draw("same");
    legend->Draw();
    sigmat1->Draw();
    bin1->Draw();
    //    CP2->Draw();
    c2->SaveAs(canvasName);
  }

  cout << "DONE" << endl;

  TCanvas *c = new TCanvas("c","c");
  c->cd();
  gScaleEB->GetYaxis()->SetRangeUser(-0.015,0.015);
  gScaleEB->GetYaxis()->SetTitle("(M_{ee}/m_{Z})^{data}-(M_{ee}/m_{Z})^{MC}");
  gScaleEB->GetYaxis()->SetTitleOffset(1.5);
  gScaleEB->GetXaxis()->SetTitle("electron p_{T} [GeV]");
  //  gScaleEB->SetMarkerColor(kRed+2);
  gScaleEB->SetMarkerStyle(20);
  gScaleEB->Draw("AP");
  c->SaveAs(((do7TeV) ? "scaleEBUpsilon-7TeV.pdf" : "scaleEBUpsilon-8TeV.pdf"));

  gScaleEE->GetYaxis()->SetRangeUser(-0.03,0.015);
  gScaleEE->GetYaxis()->SetTitle("(M_{ee}/m_{Z})^{data}-(M_{ee}/m_{Z})^{MC}");
  gScaleEE->GetYaxis()->SetTitleOffset(1.5);
  gScaleEE->GetXaxis()->SetTitle("electron p_{T} [GeV]");
  //  gScaleEE->SetMarkerColor(kRed+2);
  gScaleEE->SetMarkerStyle(20);
  gScaleEE->Draw("AP");
  c->SaveAs(((do7TeV) ? "scaleEEUpsilon-7TeV.pdf" : "scaleEEUpsilon-8TeV.pdf"));

  TFile *ff =new TFile(((do7TeV) ? "gScaleSystUpsilon-7TeV.root" : "gScaleSystUpsilon-8TeV.root"),"RECREATE");
  ff->cd();
  gScaleEB->SetName("gScaleEB");
  gScaleEB->Write();
  gScaleEE->SetName("gScaleEE");
  gScaleEE->Write();
}
