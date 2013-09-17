#include "rootheader.h"
#include "plotfigures.cc"

void plotZeeSmear(){
  
  
  TFile *f1 = new TFile("makeZeeWsShapeSmear_dm1_regver1_eleID2_escale2_smear0.root","read");

  TFile *f2 = new TFile("makeZeeWsShapeSmear_dm2_regver1_eleID2_escale0_smear0.root","read");
  TFile *f3 = new TFile("makeZeeWsShapeSmear_dm2_regver1_eleID2_escale0_smear1.root","read");
  
  gStyle->SetOptStat(0);
    
  string smearcatName[10] ={
    "e_{1}:|#eta|<1,r_{9}>0.94",
    "e_{1}:|#eta|<1,r_{9}<0.94",
    "e_{1}:1<|#eta|<1.5,r_{9}>0.94",
    "e_{1}:1<|#eta|<1.5,r_{9}<0.94",
    "e_{1}:|#eta|<1,r_{9}>0.94",
    "e_{1}:|#eta|<1,r_{9}>0.94",
    "e_{1}:|#eta|<1,r_{9}>0.94",
    "e_{1}:|#eta|<1,r_{9}<0.94",
    "e_{1}:|#eta|<1,r_{9}<0.94",
    "e_{1}:1<|#eta|<1.5,r_{9}>0.94"
    
  };
  string smearcatName2[10] ={
    "e_{2}:|#eta|<1,r_{9}>0.94",
    "e_{2}:|#eta|<1,r_{9}<0.94",
    "e_{2}:1<|#eta|<1.5,r_{9}>0.94",
    "e_{2}:1<|#eta|<1.5,r_{9}<0.94",
    "e_{2}:|#eta|<1,r_{9}<0.94",
    "e_{2}:1<|#eta|<1.5,r_{9}>0.94",
    "e_{2}:1<|#eta|<1.5,r_{9}<0.94",
    "e_{2}:1<|#eta|<1.5,r_{9}>0.94",
    "e_{2}:1<|#eta|<1.5,r_{9}<0.94",
    "e_{2}:1<|#eta|<1.5,r_{9}<0.94"
  };
  


      
  string smearcatNameEE[10] ={
    "e_{1}:1.5<|#eta|<2,r_{9}>0.94",
    "e_{1}:1.5<|#eta|<2,r_{9}<0.94",
    "e_{1}:2<|#eta|<2.5,r_{9}>0.94",
    "e_{1}:2<|#eta|<2.5,r_{9}<0.94",
    "e_{1}:1.5<|#eta|<2,r_{9}>0.94",
    "e_{1}:1.5<|#eta|<2,r_{9}>0.94",
    "e_{1}:1.5<|#eta|<2,r_{9}>0.94",
    "e_{1}:1.5<|#eta|<2,r_{9}<0.94",
    "e_{1}:1.5<|#eta|<2,r_{9}<0.94",
    "e_{1}:2<|#eta|<2.5,r_{9}>0.94"
    
  };
  string smearcatNameEE2[10] ={
    "e_{2}:1.5<|#eta|<2,r_{9}>0.94",
    "e_{2}:1.5<|#eta|<2,r_{9}<0.94",
    "e_{2}:2<|#eta|<2.5,r_{9}>0.94",
    "e_{2}:2<|#eta|<2.5,r_{9}<0.94",
    "e_{2}:1.5<|#eta|<2,r_{9}<0.94",
    "e_{2}:2<|#eta|<2.5,r_{9}>0.94",
    "e_{2}:2<|#eta|<2.5,r_{9}<0.94",
    "e_{2}:2<|#eta|<2.5,r_{9}>0.94",
    "e_{2}:2<|#eta|<2.5,r_{9}<0.94",
    "e_{2}:2<|#eta|<2.5,r_{9}<0.94"
  };
  

  
  float xmin = 75; 
  float xmax = 105; 


  for(int j=0; j<10; j++){
    
    TString histname = TString(Form("th1f_ebeb_smearmcat%d",j));
    TH1F *h1 = (TH1F*)f1->Get(histname);
    TH1F *h2 = (TH1F*)f2->Get(histname);
    TH1F *h3 = (TH1F*)f3->Get(histname);
    float x0 = h1->GetXaxis()->GetXmin();
    
    h1->Rebin(2);
    h2->Rebin(2);
    h3->Rebin(2);

    float binwidth = h1->GetBinWidth(1);
    int bin1 = int((xmin-x0)/binwidth +0.1)+1;
    int bin2 = int((xmax-x0)/binwidth +0.1);
    float data = h1->Integral(bin1,bin2);
    float mc = h2->Integral(bin1,bin2);
    float mc2 = h3->Integral(bin1,bin2);
    h2->Scale(data/mc);
    h3->Scale(data/mc2);
    
    string ytitle = string(Form("Events /%2.2f GeV",binwidth));
    
    continue; 

    plot_threeHist1FGeneral_dataAndtwoMC(h1,h2,h3,"m_{ee} (GeV)",ytitle.c_str(),xmin,xmax,0,0,1.4,"Data","Simulation","Smeared simulation","plots",histname,
					 0.5,0.77,0.94,0.92, 510,
					 0.2,0.85, "CMS 2012 Z#rightarrow ee",
					 0.2,0.8,smearcatName[j].c_str(),0.2,0.75,smearcatName2[j].c_str()
					 );
    //return; 
    
  }
  
  xmin = 75; 
  xmax = 105; 
  
  
  for(int j=0; j<10; j++){
    TString histname = TString(Form("th1f_eeee_smearmcat%d",j));
    TH1F *h1 = (TH1F*)f1->Get(histname);
    TH1F *h2 = (TH1F*)f2->Get(histname);
    TH1F *h3 = (TH1F*)f3->Get(histname);
    float x0 = h1->GetXaxis()->GetXmin();
    
    h1->Rebin(4);
    h2->Rebin(4);
    h3->Rebin(4);
    
    float binwidth = h1->GetBinWidth(1);
    int bin1 = int((xmin-x0)/binwidth +0.1)+1;
    int bin2 = int((xmax-x0)/binwidth +0.1);
    float data = h1->Integral(bin1,bin2);
    float mc = h2->Integral(bin1,bin2);
    float mc2 = h3->Integral(bin1,bin2);
    h2->Scale(data/mc);
    h3->Scale(data/mc2);
    
    string ytitle = string(Form("Events /%2.2f GeV",binwidth));
    continue; 
    plot_threeHist1FGeneral_dataAndtwoMC(h1,h2,h3,"m_{ee} (GeV)",ytitle.c_str(),xmin,xmax,0,0,1.4,"Data","Simulation","Smeared simulation","plots",histname,
					 0.5,0.77,0.94,0.92, 510,
					 0.2,0.85, "CMS 2012 Z#rightarrow ee",
					 0.2,0.8,smearcatNameEE[j].c_str(),0.2,0.75,smearcatNameEE2[j].c_str()
					 );
    //return; 
    
  }
  
  string etaname[5] = {"|#eta^{max}|<1.5,r^{min}_{9} > 0.94",
			"|#eta^{max}|<1.5,r^{min}_{9} < 0.94",
		       "|#eta^{max}|>1.5,r^{min}_{9} > 0.94",
			"|#eta^{max}|>1.5,r^{min}_{9} < 0.94",
  };
  
  for(int j=0; j<4; j++){
    TString histname = TString(Form("th1f_mpair_cat%d",j));
    TH1F *h1 = (TH1F*)f1->Get(histname);
    TH1F *h2 = (TH1F*)f2->Get(histname);
    TH1F *h3 = (TH1F*)f3->Get(histname);
    float x0 = h1->GetXaxis()->GetXmin();
    
    h1->Rebin(2);
    h2->Rebin(2);
    h3->Rebin(2);
    
    float binwidth = h1->GetBinWidth(1);
    int bin1 = int((xmin-x0)/binwidth +0.1)+1;
    int bin2 = int((xmax-x0)/binwidth +0.1);
    float data = h1->Integral(bin1,bin2);
    float mc = h2->Integral(bin1,bin2);
    float mc2 = h3->Integral(bin1,bin2);
    h2->Scale(data/mc);
    h3->Scale(data/mc2);
    
    string ytitle = string(Form("Events /%2.2f GeV",binwidth));
    
    plot_threeHist1FGeneral_dataAndtwoMC(h1,h2,h3,"m_{ee} (GeV)",ytitle.c_str(),xmin,xmax,0,0,1.4,"Data","Simulation","Smeared simulation","plots",histname,
					 0.5,0.77,0.94,0.92, 510,
					 0.2,0.85, "CMS 2012 Z#rightarrow ee",
					 0.2,0.8,etaname[j].c_str()
					 );
    //return; 
    
    //ratio
    TH1F *hhtmp = (TH1F*)h1->Clone("test");
    for(int b=1; b<= h1->GetNbinsX(); b++){
      float d = h1->GetBinContent(b);
      float m = h3->GetBinContent(b);
      float dE = h1->GetBinError(b);
      float mE = h3->GetBinError(b);
      
      float r = 0; 
      float rE = 0 ; 
      if(d>0 && m>0){
	r = d/m; 
	rE = r * sqrt( dE* dE / (d*d) + mE* mE / (m*m)); 
      }
      hhtmp->SetBinContent(b,r);
      hhtmp->SetBinError(b,rE);
    }
    //print1Fhistogram(TH1F *hhtmp,const char *xtitle,char *ytitle,int setgrid,float xmin, float xmax, int logy, float ymin, float ymax, char *dirName, const char *gifName , float text_x, float text_y, char *texName){

    TString histname2 = histname + TString("ratio");
    print1Fhistogram(hhtmp,"m_{ee} (GeV)","Data/Smeared simulation",0,xmin,xmax,0,0.5,1.5,"plots",histname2,  0.2,0.8,etaname[j].c_str());
    
    //return; 
    
  }
  
  
  
  
  
}
