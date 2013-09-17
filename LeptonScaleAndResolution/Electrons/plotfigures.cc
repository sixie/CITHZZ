
void setTCanvasNice(TCanvas *can0){
  
  
  can0->SetFillColor(0);
  can0->SetBorderMode(0);
  can0->SetBorderSize(2);
  can0->SetTickx(1);
  can0->SetTicky(1);
  can0->SetLeftMargin(0.17);
  can0->SetRightMargin(0.11);
  can0->SetTopMargin(0.06);
  can0->SetBottomMargin(0.13);
  can0->SetFrameFillStyle(0);
  can0->SetFrameBorderMode(0);
  can0->SetFrameFillStyle(0);
  can0->SetFrameBorderMode(0);
}

void setTCanvasNicev1(TCanvas *can0){
  
  
  can0->SetFillColor(0);
  can0->SetBorderMode(0);
  can0->SetBorderSize(2);
  can0->SetTickx(1);
  can0->SetTicky(1);
  can0->SetLeftMargin(0.17);
  can0->SetRightMargin(0.05);
  can0->SetTopMargin(0.05);
  can0->SetBottomMargin(0.13);
  can0->SetFrameFillStyle(0);
  can0->SetFrameBorderMode(0);
  can0->SetFrameFillStyle(0);
  can0->SetFrameBorderMode(0);
}


void setTLegendNice(TLegend *leg){
  leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetShadowColor(0);

  
}


void setTH1BinErorZero(TH1 *h1){
   int nbin1 = h1->GetNbinsX();
   for(int b=1; b<=nbin1; b++){
     h1->SetBinError(b,0);
   }
}



void plot_threeHist1FGeneral_dataAndtwoMC(TH1F *hhtmp1, TH1F *hhtmp2,  TH1F *hhtmp3,const  char *xtitle, const char *ytitle,float xmin, float xmax, int logy, float ymin, float ykmax, const char *leg1Name, const char *leg2Name, char *leg3Name,char *dirName, const char *gifName , float legx1=0.55,float legy1=0.75,float legx2 = 0.9,float legy2 =0.9, int ndivx = 510,  float text_x = -1, float text_y = -1, const char *texName= "", float text_x2 = -1,float text_y2 = -1, const char *texName2= "CMS Simulation", float text_x3 = -1,float text_y3 = -1, const char *texName3= ""){
  
  
  TH1F *hh[3] = {hhtmp1, hhtmp2,hhtmp3};
  
  TCanvas *can0 = new TCanvas("can0","c000",200,10,550,500);
  can0->SetLogy(logy);

  if(xmax<0.01){
    setTCanvasNice(can0);
  }else{
    setTCanvasNicev1(can0);
  }
  
  
  float ymax = hhtmp1->GetMaximum() > hhtmp2->GetMaximum() ? hhtmp1->GetMaximum() : hhtmp2->GetMaximum();
  ymax *= ykmax; 

  hh[1]->GetYaxis()->SetRangeUser(ymin,ymax);
  hh[1]->GetXaxis()->SetRangeUser(xmin,xmax);
  hh[1]->GetXaxis()->SetTitle(xtitle);
  hh[1]->GetYaxis()->SetTitle(ytitle);
  hh[1]->SetLineWidth(2);
  hh[2]->SetLineWidth(2);
  //hh[0]->SetMarkerSize(1.5);
  hh[1]->GetXaxis()->SetNdivisions(ndivx);
  hh[1]->GetYaxis()->SetNdivisions(512);
  
  
  //hh[1]->SetLineColor(kRed);
  // hh[1]->SetMarkerColor(kRed);
  
  
  //hh[2]->SetLineColor(kBlue);
  //hh[2]->SetMarkerColor(kBlue);
  
  //hh[0]->Draw("e"); // data
    
  hh[1]->SetLineStyle(2);
  
  hh[1]->Draw("hist"); //MC
  hh[2]->Draw("histsame");
  hh[0]->Draw("esame"); // data
  
  
  //TLegend *leg1 = new TLegend(0.55,0.75,0.9,0.9,NULL,"brNDC");
  TLegend *leg1 = new TLegend(legx1,legy1,legx2,legy2,NULL,"brNDC");
  
  
  string sName = gifName; 

  if( sName.substr(0,13) == "npizRec_nclus"){
    leg1 = new TLegend(0.43,0.23,0.88,0.48,NULL,"brNDC");
  }
  if( sName.substr(0,11) == "dndeta_data"){
    leg1 = new TLegend(0.33,0.21,0.78,0.47,NULL,"brNDC");
  }
  if( sName.substr(0,15) == "dndeta_datacorr"){
    leg1 = new TLegend(0.30,0.69,0.75,0.94,NULL,"brNDC");
  }
  if( sName.substr(0,10) == "mpair_ebeb"){
    leg1 = new TLegend(0.6,0.8,0.94,0.92,NULL,"brNDC");
  }

  
  setTLegendNice(leg1); 
  
  
  leg1->AddEntry(hh[0],leg1Name,"p");
  leg1->AddEntry(hh[1],leg2Name,"l");
  leg1->AddEntry(hh[2],leg3Name,"l");
  leg1->Draw();
  


  
  if( text_x >0 && text_y >0){
    TLatex *   tex = new TLatex(text_x, text_y, texName);
    tex->SetNDC();
    tex->SetTextSize(0.038);
    tex->Draw();
  }
  
  
  
  if( text_x2 >0 && text_y2 >0){
    TLatex *   tex = new TLatex(text_x2, text_y2, texName2);
    tex->SetNDC();
    tex->SetTextSize(0.038);
    ///tex->SetLineWidth(2);
    tex->Draw();
  }
  
  if( text_x3 >0 && text_y3 >0){
    TLatex *   tex = new TLatex(text_x3, text_y3, texName3);
    tex->SetNDC();
    tex->SetTextSize(0.038);
    tex->Draw();
  }
  
  
  string histName = string(Form("%s/%s.pdf",dirName,gifName));
  can0->Print(histName.c_str());
  histName = string(Form("%s/%s.C",dirName,gifName));
  can0->Print(histName.c_str());
  
}


void print1Fhistogram(TH1F *hhtmp,const char *xtitle,const char *ytitle,int setgrid,float xmin, float xmax, int logy, float ymin, float ymax, char *dirName, const char *gifName , float text_x, float text_y, const char *texName){

  
   TCanvas *can0 = new TCanvas("can0","c000",200,10,550,500);
   setTCanvasNice(can0);
   
   can0->SetLogy(logy);
   
   if(setgrid>=1){
     can0->SetGridx(setgrid);
     can0->SetGridy(setgrid);
   }
   
  hhtmp->GetXaxis()->SetTitle(xtitle);
  hhtmp->GetXaxis()->SetRangeUser(xmin,xmax);
  hhtmp->GetYaxis()->SetRangeUser(ymin,ymax);
  hhtmp->GetYaxis()->SetTitle(ytitle);
  hhtmp->Draw();
  
  if( text_x >0 && text_y >0){
    TLatex *   tex = new TLatex(text_x, text_y, texName);
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->Draw();
  }
  
  string histName = string(Form("%s/%s.pdf",dirName,gifName));
  can0->Print(histName.c_str());
  histName = string(Form("%s/%s.C",dirName,gifName));
  can0->Print(histName.c_str());
  
}



void plot_twoHist1FGeneral(TH1 *hhtmp1, TH1 *hhtmp2, const char *xtitle, const char *ytitle,float xmin, float xmax, int logy, float ymin, float ykmax, char *leg1Name, char *leg2Name, char *dirName, const char *gifName , float text_x= -1, float text_y = -1, const char *texName = "", float text_x1 = -1, float text_y1 = -1, const char *texName1 = "", float text_x2 =-1,float text_y2 = -1,const char*texName2 = ""){
  
  
  gStyle->SetOptStat(0);
  

  TH1 *hh[6] = {hhtmp1, hhtmp2};
    
  TCanvas *can0 = new TCanvas("can0","c000",200,10,550,500);
  setTCanvasNice(can0);
  can0->SetLogy(logy);
  
  //can0->SetGridx(1);
  //can0->SetGridy(1);
  
  float ymax = hhtmp1->GetMaximum() > hhtmp2->GetMaximum() ? hhtmp1->GetMaximum() : hhtmp2->GetMaximum();
  ymax *= ykmax; 

  hh[1]->GetYaxis()->SetRangeUser(ymin,ymax);
  hh[1]->GetXaxis()->SetRangeUser(xmin,xmax);
  hh[1]->GetXaxis()->SetTitle(xtitle);
  hh[1]->GetYaxis()->SetTitle(ytitle);
  hh[1]->SetLineWidth(2);
  hh[0]->SetLineWidth(2);
  hh[1]->SetLineStyle(1);
   hh[1]->SetLineColor(kRed);
   hh[1]->SetMarkerColor(kRed);
   hh[0]->SetLineColor(kBlue);
   hh[0]->SetMarkerColor(kBlue);
   hh[0]->SetMarkerSize(1);
   //hh[0]->SetMarkerStyle(4);
    
  
  hh[1]->Draw("hist");
  //hh[1]->Draw("");
  
  hh[0]->Draw("same");
  
  
  TLegend *leg1 = new TLegend(0.52,0.8,0.87,0.92,NULL,"brNDC");
  
  
  setTLegendNice(leg1); 
  
  
  leg1->AddEntry(hh[0],leg1Name,"p");
  leg1->AddEntry(hh[1],leg2Name,"l");
  leg1->Draw();
  
  
//    string etaname[6] = {"|#eta^{max}|<1.5,r^{min}_{9} > 0.94",
// 			"|#eta^{max}|<1.5,r^{min}_{9} < 0.94",
// 			"|#eta^{max}|>1.5,r^{min}_{9} > 0.94",
// 			"|#eta^{max}|>1.5,r^{min}_{9} < 0.94",
// 			"VBF category",
// 			"Combined categories",
//   };
 
  
//   if( diphtcat>=0 && diphtcat <=5){
//     TLatex l; 
//     l.SetNDC();
//     l.SetTextSize(0.04);
//     l.DrawLatex(0.2,0.8,etaname[diphtcat].c_str());
//   }
  
  
  if( text_x >0 && text_y >0){
    TLatex *   tex = new TLatex(text_x, text_y, texName);
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    tex->Draw();
  }
   
  if( text_x1 >0 && text_y1 >0){
    TLatex *   tex = new TLatex(text_x1, text_y1, texName1);
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    tex->Draw();
  }

  if( text_x2 >0 && text_y2 >0){
    TLatex *   tex = new TLatex(text_x2, text_y2, texName2);
    tex->SetNDC();
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  
  char *histName = new char[100];
  sprintf(histName,"%s/%s.gif",dirName,gifName);
  can0->Print(histName);
  sprintf(histName,"%s/%s.pdf",dirName,gifName);
  can0->Print(histName);
  sprintf(histName,"%s/%s.C",dirName,gifName);
  can0->Print(histName);
  
  
}



///sum over the cut and overflow
void sumupOverFlow(TH1 *hhtmp, float xend ){
  
  float sum = 0; 
  int nbins = hhtmp->GetNbinsX();
  
  int nbin_end = int( (xend - hhtmp->GetXaxis()->GetXmin())/hhtmp->GetBinWidth(1));
  for(int bin = nbin_end; bin<= nbins+1; bin++){ /// bins +1 for overflowsbin
    sum += hhtmp->GetBinContent(bin);
  }
  hhtmp->SetBinContent(nbin_end,sum);
  
}

/// Input, others, z+jets, z+gamma, data , signal 
void plot_stack_two(TH1F *hhtmp1, TH1F *hhtmp2, TH1F *hhtmp3, int flagUnit, char *xtitle ,
		    char *leg1Name,char *leg2Name,char *leg3Name,
		    int logy, float ymin, float ymax, int xstart, int xend, char *dirName, const char *gifName){
  
  TCanvas *can00s = new TCanvas("can00s", "c000",205,47,550,500);
   gStyle->SetOptStat(0);
   
   setTCanvasNicev1(can00s);
   can00s->SetLogy(logy);
   can00s->SetTopMargin(0.06991526);
   char *filename = new char[100];
   
   THStack *hs = new THStack();
   hs->SetName("hs");
   hs->SetTitle("test stacked histograms");
   hs->SetMinimum(ymin);
   hs->SetMaximum(ymax);
   
   
   
   int nbins = hhtmp1->GetNbinsX();
   float xmin = hhtmp1->GetXaxis()->GetXmin(); 
   float xmax = hhtmp1->GetXaxis()->GetXmax(); 
   
  

   //int icol[10]= {3,2,920,70,150,200};

  //  hhtmp1->SetFillColor(icol[0]);
//    hhtmp2->SetFillColor(icol[1]);
//    hhtmp3->SetFillColor(icol[2]);
//    hhtmp4->SetFillColor(icol[3]);

   
   hhtmp1->SetFillColor(kBlue-5);
   //hhtmp3->SetFillColor(kGreen-5);
   //   hhtmp4->SetFillColor(kWhitle
   

   //hhtmp5->SetFillColor(icol[4]);

   // sumupOverFlow(hhtmp1,xend);
   //   sumupOverFlow(hhtmp2,xend);

  //  sumupOverFlow(hhtmp3,xend);
//    sumupOverFlow(hhtmp4,xend);
//    sumupOverFlow(hhtmp5,xend);

   hs->Add(hhtmp1 ); 
   hs->Add(hhtmp2 ); 
//    hs->Add(hhtmp3 ); 
//    hs->Add(hhtmp4 ); 
   // hs->Add(hhtmp5 ); 
   
   
   
   hs->Draw();
   hs->GetXaxis()->SetTitle(xtitle);
   hs->GetXaxis()->SetRangeUser(xstart+1,xend-1);
    
   
   if(flagUnit==1){
     sprintf(filename,"Events / %2.2f GeV/c ",  hhtmp1->GetBinWidth(1));
   }else if(flagUnit==2) {
     sprintf(filename,"Events / %2.2f GeV/c^{2}", hhtmp1->GetBinWidth(1));
   }else if(flagUnit==3) {
     sprintf(filename,"Events / %2.2f GeV", hhtmp1->GetBinWidth(1));
   }

   hs->GetYaxis()->SetTitle(filename);
   
   
  //  //signal point
//    TH1F *hh_signal = (TH1F*)hhtmp6->Clone();
//    sumupOverFlow(hh_signal,xend);
//    hh_signal->SetLineWidth(2);
//    hh_signal->Draw("histsame");
//    hh_signal ->SetLineColor(kRed);
   
   hhtmp3->Draw("esame");
   
   
   // cout<<"dddddd."<<endl;

   //TLegend *leg = new TLegend(0.45,0.65,0.9,0.88,NULL,"brNDC");
   TLegend *leg = new TLegend(0.5,0.68,0.9,0.88,NULL,"brNDC");
   setTLegendNice(leg);
   
   //leg->AddEntry(hhtmp6,leg6Name,"p");
   leg->AddEntry(hhtmp3,leg3Name,"p");
   leg->AddEntry(hhtmp2,leg2Name,"f");
   leg->AddEntry(hhtmp1,leg1Name,"f");
   
   
   cout<<"ddddde."<<endl;
   leg->Draw();
   
   string histName = string(Form("%s/%s.pdf",dirName,gifName));
   can00s->Print(histName.c_str());
   histName = string(Form("%s/%s.C",dirName,gifName));
   can00s->Print(histName.c_str());
   
   
}



/// Input, others, z+jets, z+gamma, data , signal 
void plot_stack_five(TH1F *hhtmp1, TH1F *hhtmp2, TH1F *hhtmp3, TH1F *hhtmp4,TH1F *hhtmp5,TH1F *hhtmp6,int flagUnit, char *xtitle ,
		     char *leg1Name,char *leg2Name,char *leg3Name,char *leg4Name,char *leg5Name,char *leg6Name,
		     int logy, float ymin, float ymax, int xstart, int xend, char *dirName, const char *gifName, float lumiData){
  
  TCanvas *can00s = new TCanvas("can00s", "c000",205,47,550,500);
   gStyle->SetOptStat(0);
   
   setTCanvasNicev1(can00s);
   can00s->SetLogy(logy);
   can00s->SetTopMargin(0.06991526);
   
   
   THStack *hs = new THStack();
   hs->SetName("hs");
   hs->SetTitle("test stacked histograms");
   hs->SetMinimum(ymin);
   hs->SetMaximum(ymax);
   
   
   
   int nbins = hhtmp1->GetNbinsX();
   float xmin = hhtmp1->GetXaxis()->GetXmin(); 
   float xmax = hhtmp1->GetXaxis()->GetXmax(); 
   
   char *filename = new char[100];

   //int icol[10]= {3,2,920,70,150,200};

  //  hhtmp1->SetFillColor(icol[0]);
//    hhtmp2->SetFillColor(icol[1]);
//    hhtmp3->SetFillColor(icol[2]);
//    hhtmp4->SetFillColor(icol[3]);

   
   hhtmp1->SetFillColor(kCyan-6);
   hhtmp2->SetFillColor(kBlue-8);
   hhtmp3->SetFillColor(kGreen-5);
   //   hhtmp4->SetFillColor(kWhitle
   

   //hhtmp5->SetFillColor(icol[4]);

   sumupOverFlow(hhtmp1,xend);
   sumupOverFlow(hhtmp2,xend);
   sumupOverFlow(hhtmp3,xend);
   sumupOverFlow(hhtmp4,xend);
   sumupOverFlow(hhtmp5,xend);
   

   hs->Add(hhtmp1 ); 
   hs->Add(hhtmp2 ); 
   hs->Add(hhtmp3 ); 
   hs->Add(hhtmp4 ); 
   // hs->Add(hhtmp5 ); 
   
   
   
   hs->Draw();
   hs->GetXaxis()->SetTitle(xtitle);
   hs->GetXaxis()->SetRangeUser(xstart+1,xend-1);
    
   
   if(flagUnit==1){
     sprintf(filename,"Events / %2.2f GeV/c ",  hhtmp1->GetBinWidth(1));
   }else if(flagUnit==2) {
     sprintf(filename,"Events / %2.2f GeV/c^{2}", hhtmp1->GetBinWidth(1));
   }else if(flagUnit==3) {
     sprintf(filename,"Events / %3.3f", hhtmp1->GetBinWidth(1));
   }

   hs->GetYaxis()->SetTitle(filename);
   
   
   //signal point
   TH1F *hh_signal = (TH1F*)hhtmp6->Clone();
   sumupOverFlow(hh_signal,xend);
   hh_signal->SetLineWidth(2);
   hh_signal->Draw("histsame");
   hh_signal ->SetLineColor(kRed);
   
   hhtmp5->Draw("esame");
   
   
   // cout<<"dddddd."<<endl;

   //TLegend *leg = new TLegend(0.45,0.65,0.9,0.88,NULL,"brNDC");
   TLegend *leg = new TLegend(0.5,0.58,0.9,0.88,NULL,"brNDC");
   setTLegendNice(leg);
   
   //leg->AddEntry(hhtmp6,leg6Name,"p");
   leg->AddEntry(hhtmp5,leg5Name,"p");
   leg->AddEntry(hhtmp4,leg4Name,"f");
   leg->AddEntry(hhtmp3,leg3Name,"f");
   leg->AddEntry(hhtmp2,leg2Name,"f");
   leg->AddEntry(hhtmp1,leg1Name,"f");
   leg->AddEntry(hh_signal,leg6Name,"l");
   
   
   
   cout<<"ddddde."<<endl;
   leg->Draw();
   
   sprintf(filename,"CMS Preliminary 7 TeV %3.0f pb^{-1}",lumiData);
   TLatex *   tex = new TLatex(0.32,0.945,filename);
   tex->SetNDC();
   tex->SetTextSize(0.05);
   tex->SetLineWidth(2);
   tex->Draw();
   can00s->Modified();
   can00s->cd();
   can00s->SetSelected(can00s);
     string histName = string(Form("%s/%s.pdf",dirName,gifName));
  can00s->Print(histName.c_str());
  histName = string(Form("%s/%s.C",dirName,gifName));
  can00s->Print(histName.c_str());
      
}
