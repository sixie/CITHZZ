#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <RooDataSet.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveText.h>

struct fitbinning {
  enum ptBins {kPtLow = 0,
               kPtHigh = 1}; 
  
  enum etaBins {kCentralEB = 0,
                kOuterEB,
                kEE};
};

void drawSPlot(int ptbin, int etabin, const char *varname, 
               const char *axistitle, const char *extracut, int nbins, float min, float max, int logy=0) {

  // Load results file...
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0) ;
  gStyle->SetOptTitle(0) ;

  char cut[200];
  if(ptbin==fitbinning::kPtLow) sprintf(cut,"pt>=7 && pt<12");
  if(ptbin==fitbinning::kPtHigh) sprintf(cut,"pt>=12");
  if(etabin==fitbinning::kCentralEB) sprintf(cut,"%s && abs(eta)<0.8",cut);
  if(etabin==fitbinning::kOuterEB) sprintf(cut,"%s && abs(eta)>=0.8 && abs(eta)<1.479",cut);
  if(etabin==fitbinning::kEE) sprintf(cut,"%s && abs(eta)>=1.479",cut);

  char finalcut[500];
  sprintf(finalcut,"%s && %s",cut,extracut);

  // MC
  TFile *fileMC = TFile::Open("root://eoscms//eos/cms/store/group/phys_egamma/emanuele/elereg/eletrees/V1/ZjetsMad_Summer12.root");
  TFile *fileMC2 = TFile::Open("root://eoscms//eos/cms/store/group/phys_egamma/emanuele/elereg/eletrees/V1/RelValZEE_g4emtest-START50_V8_special.root");
  TTree *treeMC = (TTree*)fileMC->Get("eleIDdir/T1");
  TTree *treeMC2 = (TTree*)fileMC2->Get("eleIDdir/T1");

  // data with sWeights
  char filesmall[200];
  sprintf(filesmall,"sPlots/electrons_ptbin%d_etabin%d.root",ptbin,etabin);
  TFile *file = TFile::Open(filesmall);
  TTree *tree = (TTree*)file->Get("eleIDdir/T1");
  char filesweight[200];
  sprintf(filesweight,"sPlots/sweights_ptbin%d_etabin%d.root",ptbin,etabin);
  tree->AddFriend("sweights=dataset",filesweight);

  TH1F *signalMC = new TH1F(varname, "", nbins, min, max);
  treeMC->Project(varname,varname,finalcut);
  TH1F *signalMC2 = new TH1F((string(varname)+string("2")).c_str(), "", nbins, min, max);
  treeMC2->Project((string(varname)+string("2")).c_str(),varname,finalcut);

  char splotname[200];
  sprintf(splotname,"%s_sPlotSig",varname);
  TH1F *signalsPlot = new TH1F(splotname,"", nbins, min, max);
  
  char finalcutweight[600];
  sprintf(finalcutweight,"(%s)*N_sig_sw",finalcut);
  tree->Project(splotname,varname,finalcutweight);
  signalsPlot->Sumw2();

  // sPlot has the correct normalization. Normalize the MC to this area
  float integral = signalsPlot->Integral();
  signalMC->Scale(integral/signalMC->Integral());
  signalMC2->Scale(integral/signalMC2->Integral());

  float maxY = TMath::Max(signalsPlot->GetMaximum(),signalMC->GetMaximum());
  maxY = TMath::Max(maxY,float(signalMC2->GetMaximum()));

  TPaveText *text = new TPaveText(0.65,0.50,0.90,0.70,"brNDC");
  if(ptbin==0) text->AddText("p_{T}<10 GeV");
  else text->AddText("p_{T}>10 GeV");
  if(etabin==0) text->AddText("|#eta|<0.8");
  else if(etabin==1) text->AddText("0.8<|#eta|<1.479");
  else text->AddText("1.479<|#eta|<2.5");
  text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextAlign(12);
  text->SetTextFont(132);
  text->SetTextSize(0.04);

  TLegend* leg = new TLegend(0.65,0.70,0.90,0.80);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(signalMC,"Z(ee) MC (Summer12)","fl");
  leg->AddEntry(signalMC2,"Z(ee) MC (G4 emtest)","fl");
  leg->AddEntry(signalsPlot,"s-weighted data","pl");

  TCanvas c1;
  if(logy) c1.SetLogy(1);
  signalMC->SetLineColor(kOrange+7);
  signalMC->SetFillColor(kOrange+7);
  signalMC->SetLineWidth(1.5);
  signalMC->SetFillStyle(3005);

  signalMC2->SetLineColor(kAzure-6);
  signalMC2->SetFillColor(kAzure-6);
  signalMC2->SetLineWidth(1.5);
  signalMC2->SetFillStyle(3004);

  signalsPlot->SetMarkerStyle(8);
  signalsPlot->SetMarkerSize(1);
  signalsPlot->GetXaxis()->SetTitle(axistitle);
  signalsPlot->GetYaxis()->SetTitle("weighted events");

  signalsPlot->SetMaximum(maxY+sqrt(maxY));

  signalsPlot->Draw("pe1");
  signalMC->Draw("same hist");
  signalMC2->Draw("same hist");
  leg->Draw();
  text->Draw();

  char figfilename[500];
  sprintf(figfilename,"%s_ptbin%d_etabin%d_sPlot.png",varname,ptbin,etabin);
  c1.SaveAs(figfilename);
  sprintf(figfilename,"%s_ptbin%d_etabin%d_sPlot.pdf",varname,ptbin,etabin);
  c1.SaveAs(figfilename);

}
  

void makeSPlots() {

  gSystem->Load("libRooFit");
  
  // barrel
  for(int ieta=0;ieta<2;ieta++) {
    drawSPlot(1,ieta,"EoP","E/P_{in}","vertices<10",50,0.0,3.0);
    drawSPlot(1,ieta,"HoE","H/E","vertices<10",50,0.0,0.1,1);
    drawSPlot(1,ieta,"deta","#Delta #eta_{in}","vertices<10",50,-0.02,0.02);
    drawSPlot(1,ieta,"dphi","#Delta #phi_{in}","vertices<10",50,-0.1,0.1);
    drawSPlot(1,ieta,"see","#sigma_{i#eta i#eta}","vertices<10",50,0.0,0.02);
    //     drawSPlot(1,ieta,"sep","#sigma_{i#eta i#phi}","vertices<10",50,0.0,0.05);
    //     drawSPlot(1,ieta,"spp","#sigma_{i#phi i#phi}","vertices<10",50,0.0,0.05);
    drawSPlot(1,ieta,"etawidth","#eta width","vertices<10",50,0.0,0.02);
    drawSPlot(1,ieta,"phiwidth","#phi width","vertices<10",50,0.0,0.2);
    drawSPlot(1,ieta,"missHits","miss pix hits","vertices<10",50,0.0,10.0);
    drawSPlot(1,ieta,"fbrem","f_{brem}","vertices<10",50,-0.2,1.0);
    drawSPlot(1,ieta,"R9","R9","vertices<10",50,0.0,1.2); 

    drawSPlot(0,ieta,"EoP","E/P_{in}","vertices<10",20,0.0,3.0);
    drawSPlot(0,ieta,"HoE","H/E","vertices<10",30,0.0,0.1,1);
    drawSPlot(0,ieta,"deta","#Delta #eta_{in}","vertices<10",30,-0.02,0.02);
    drawSPlot(0,ieta,"dphi","#Delta #phi_{in}","vertices<10",15,-0.1,0.1);
    drawSPlot(0,ieta,"see","#sigma_{i#eta i#eta}","vertices<10",30,0.0,0.05);
    //     drawSPlot(0,ieta,"sep","#sigma_{i#eta i#phi}","vertices<10",30,0.0,0.05);
    //     drawSPlot(0,ieta,"spp","#sigma_{i#phi i#phi}","vertices<10",30,0.0,0.05);
    drawSPlot(0,ieta,"etawidth","#eta width","vertices<10",30,0.0,0.06);
    drawSPlot(0,ieta,"phiwidth","#phi width","vertices<10",30,0.0,0.2);
    drawSPlot(0,ieta,"missHits","miss pix hits","vertices<10",10,0.0,10.0);
    drawSPlot(0,ieta,"fbrem","f_{brem}","vertices<10",20,-0.2,1.0);
    drawSPlot(0,ieta,"R9","R9","vertices<10",10,0.0,1.2); 
  }

  // endcap
  drawSPlot(1,2,"EoP","E/P_{in}","vertices<10",50,0.0,5.0);
  drawSPlot(1,2,"HoE","H/E","vertices<10",50,0.0,0.1,1);
  drawSPlot(1,2,"deta","#Delta #eta_{in}","vertices<10",50,-0.04,0.04);
  drawSPlot(1,2,"dphi","#Delta #phi_{in}","vertices<10",50,-0.1,0.1);
  drawSPlot(1,2,"see","#sigma_{i#eta i#eta}","vertices<10",50,0.0,0.05);
  //   drawSPlot(1,2,"sep","#sigma_{i#eta i#phi}","vertices<10",50,0.0,0.05);
  //   drawSPlot(1,2,"spp","#sigma_{i#phi i#phi}","vertices<10",50,0.03,0.2);
  drawSPlot(1,2,"etawidth","#eta width","vertices<10",50,0.0,0.06);
  drawSPlot(1,2,"phiwidth","#phi width","vertices<10",50,0.0,0.2);
  drawSPlot(1,2,"missHits","miss pix hits","vertices<10",50,0.,10.);
  drawSPlot(1,2,"fbrem","f_{brem}","vertices<10",50,-0.2,1.0);
  drawSPlot(1,2,"R9","R9","vertices<10",50,0.0,1.2);

  drawSPlot(0,2,"EoP","E/P_{in}","vertices<10",20,0.0,5.0);
  drawSPlot(0,2,"HoE","H/E","vertices<10",25,0.0,0.1,1);
  drawSPlot(0,2,"deta","#Delta #eta_{in}","vertices<10",25,-0.04,0.04);
  drawSPlot(0,2,"dphi","#Delta #phi_{in}","vertices<10",15,-0.1,0.1);
  drawSPlot(0,2,"see","#sigma_{i#eta i#eta}","vertices<10",20,0.0,0.1);
  //   drawSPlot(0,2,"sep","#sigma_{i#eta i#phi}","vertices<10",25,0.0,0.05);
  //   drawSPlot(0,2,"spp","#sigma_{i#phi i#phi}","vertices<10",25,0.03,0.2);
  drawSPlot(0,2,"etawidth","#eta width","vertices<10",15,0.0,0.1);
  drawSPlot(0,2,"phiwidth","#phi width","vertices<10",15,0.0,0.2);
  drawSPlot(0,2,"missHits","miss pix hits","vertices<10",10,0.,10.);
  drawSPlot(0,2,"fbrem","f_{brem}","vertices<10",20,-0.2,1.0);
  drawSPlot(0,2,"R9","R9","vertices<10",10,0.0,1.2);

}

