
#include "../Reducer/rootheader.h"
#include "../Selector/roofitheader.h"
int debug_ = -2; 

using namespace RooFit ;

#include "../Reducer/RecoAnalyzer.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#endif
using namespace TMVA;



#include "TMinuit.h"

#include "/home/yongy/backup/plot_stack.cc"

//#include "th1fmorph.C"

// float fitrangeLow = 91.1876 - 2*3.2; 
// float fitrangeHigh = 91.1986 + 2*3.2; 

double fitrangeLow = 75;
double fitrangeHigh = 105;
double binwidth; 

int transferfun = 0;

TH1F *hh_mcinput;
TH1F *hh_datainput;
TH1F *hhdata2;

TChain *fChainMC; 
TChain *fChainData; 

int NMC; 
int NData;

int etcut = 25;
int nsccat = 8;

map<string,TH1F*> map_th1f; 

float smearStep = 0.002;


int fitdet;


  Float_t         elescr9[2];
   Float_t         elesceta[2];
   Float_t         eleen[2];
   Float_t         eleen0[2];
   Float_t         eleenscraw[2];
   //Float_t         mpair;
   Float_t         mpair0;
   Float_t         eleeta[2];
   Int_t           eleieta[2];
   Int_t           eleiphi[2];
 //   Int_t           runNumber;
//    Int_t           evtNumber;
//    Int_t           lumiBlock;
   Float_t         eleetrue[2];
   Float_t         eleetruenofsr[2];
   Float_t         eleetrueplusfsrdrbc003[2];
   Float_t         eleetrueplusfsrdrbc004[2];
   Float_t         eleetrueplusfsrdrbc005[2];
   Float_t         eleetrueplusfsrdrbc006[2];
   Float_t         eleetrueplusfsrdrbc007[2];
   Float_t         eleetrueplusfsrdrbc008[2];
   Float_t         eleetrueplusfsrdrbc009[2];
   Float_t         eleetrueplusfsrdrbc01[2];
   Int_t           eletrueelematched[2];
   Int_t           eletrueelematchednofsr[2];
   Float_t         weight;




//int phtcorr = 360;
//int phtcorr = 358;
//int phtcorr = 94;

#include "../Reducer/testSelection.h"
#include "../Reducer/preselectionCut.cc"
#include "../Reducer/physUtils.cc"
#include "../Reducer/utils.cc"
#include "../Reducer/energyScaleCorrectionnew.cc"

string fitmethod;
bool printHist; 

#include "..//Selector/roodatasetth1.cc"

ofstream txtout;

int testOnlyCat;

string smearMethod;
string testphtcat;


vector<int> vsc1cat;
vector<int> vsc2cat;
vector<float> vele1eta;
vector<float> vele1r9;
vector<float> vele2r9;
vector<float> vele1etrue;
vector<float> vele2etrue;
vector<int> vele1ieta;
vector<int> vele1iphi;
vector<int> vele2ieta;
vector<int> vele2iphi;

vector<float> vele2eta;
vector<float> vele1en;
vector<float> vele2en;
vector<float> vmpair;
vector<float> vweight;
vector<int> vindpair;


// int scCategoryFour(float sc_eta, float sc_r9){
  
//   if( fabs(sc_eta)<=1.48){
//     return sc_r9 >0.94  ? 0 : 1; 
//   }else{
//     return sc_r9 >0.94  ? 2 : 3; 
//   }
  
// }


// int scCategoryEight(float sc_eta, float sc_r9){
  
//   if( fabs(sc_eta)<=1.0){
//     return sc_r9 >0.94  ? 0 : 1; 
//   }else if( fabs(sc_eta)<=1.48){
//     return sc_r9 >0.94 ? 2: 3;
//   }else if( fabs(sc_eta)<=2.0){
//     return sc_r9 > 0.94 ? 4:5; 
//   }
//   else{
//     return sc_r9 >0.94  ? 6 : 7; 
//   }
  
// }


int scEnergy4Bin(float en){
  if( en>=25 && en<=45){
    return 0; 
  }else if( en<=55){
    return 1; 
  }else if( en<=70){
    return 2; 
  }else if( en<=150){
    return 3;
  }
  else return -1;
  
}

int scCategoryEightEn4Bin(float en, float sc_eta, float sc_r9){
  
  int enbin = scEnergy4Bin(en);
  if(enbin<0) return -1;
  
  
  if( fabs(sc_eta)<=1.0){
    return sc_r9 >0.94  ? 2*enbin : 2*enbin+1; 
  }else if( fabs(sc_eta)<=1.48){
    return sc_r9 >0.94  ? 8+2*enbin : 9+2*enbin; 
  }else if( fabs(sc_eta)<=2.0){
    return sc_r9 >0.94  ? 16+2*enbin : 17+2*enbin; 
  }else{
    return sc_r9 >0.94  ? 24+2*enbin : 25+2*enbin; 
  }
  
  
  
}


void exChangeTwoNumber(int &tmp1, int &tmp2){
  int tmp = tmp1; 
  tmp1 = tmp2;
  tmp2 = tmp;
}
    

// int scCategoryTen(float sc_eta, float sc_r9, int sc_ieta, int sc_iphi){
  
//   if( fabs(sc_eta)<1.0){
//     if( !isAtEcalBarrelModuleCracksv1(sc_ieta,sc_iphi)){
//       return sc_r9 >0.94  ? 0 : 1; 
//     }else{
//       return sc_r9 >0.94  ? 2 : 3; 
//     }
//   }else if( fabs(sc_eta)<1.48){
//     return sc_r9 >0.94 ? 4: 5;
//   }else if( fabs(sc_eta)<2.0){
//     return sc_r9 > 0.94 ? 6:7; 
//   }
//   else{
//     return sc_r9 >0.94  ? 8:9; 
//   }
  
// }

map<int, map<int,int> > map_indscpairBarrel; 

map<int, map<int,int> > map_indscpairEndcap; 



map<int,int> map_indsc1;
map<int,int> map_indsc2;


void fillMapIndex(){
  map_indscpairBarrel[0][0] = 0;
  map_indscpairBarrel[1][1] = 1;
  map_indscpairBarrel[2][2] = 2;
  map_indscpairBarrel[3][3] = 3;
  
  map_indscpairBarrel[0][1] = 4;
  map_indscpairBarrel[0][2] = 5;
  map_indscpairBarrel[0][3] = 6;
  map_indscpairBarrel[1][2] = 7;
  map_indscpairBarrel[1][3] = 8;
  map_indscpairBarrel[2][3] = 9;
  

  ///Endcap 4,5,6,7
  map_indscpairEndcap[4][4] = 0;
  map_indscpairEndcap[5][5] = 1;
  map_indscpairEndcap[6][6] = 2;
  map_indscpairEndcap[7][7] = 3;
  map_indscpairEndcap[4][5] = 4;
  map_indscpairEndcap[4][6] = 5;
  map_indscpairEndcap[4][7] = 6;
  map_indscpairEndcap[5][6] = 7;
  map_indscpairEndcap[5][7] = 8;
  map_indscpairEndcap[6][7] = 9;
  

  //get indscpair 
  map_indsc1[0] = 0;
  map_indsc2[0] = 0;
  map_indsc1[1] = 1;
  map_indsc2[1] = 1;
  map_indsc1[2] = 2;
  map_indsc2[2] = 2;
  map_indsc1[3] = 3;
  map_indsc2[3] = 3;
  map_indsc1[4] = 0;
  map_indsc2[4] = 1;
  map_indsc1[5] = 0;
  map_indsc2[5] = 2;
  map_indsc1[6] = 0;
  map_indsc2[6] = 3;
  map_indsc1[7] = 1;
  map_indsc2[7] = 2;
  map_indsc1[8] = 1;
  map_indsc2[8] = 3;
  map_indsc1[9] = 2;
  map_indsc2[9] = 3;
  
  
}


////16 category,  4enx2etax2r9
void fillMapIndex2(){
  
  int n = 0; 
  for(int j=0;j<16;j++){
    for(int k=j;k<16;k++){
      map_indscpairBarrel[j][k] = n;

      cout<<" map_indscpairBarrel " << j<<" "<<k <<" "<< n <<endl; 

      n ++; 
    }
  }
  
}








double getInterpolatedvaluefromTH1(TH1 *hh1, float val){
  
  double total = hh1->Integral();
  if(total !=1){
    hh1->Scale(1.0/total);
  }
  
  
  int bin = hh1->FindFixBin(val);
  binwidth =  hh1->GetBinWidth(bin);
  int nbins = hh1->GetNbinsX();
  
  
  if( val < hh1->GetXaxis()->GetXmin()) return 0; 
  if( val > hh1->GetXaxis()->GetXmax()) return 0; 
  
  double pdf = 0;
  if( bin==1) {
    double y1 = hh1->GetBinContent(1);
    double y2 = hh1->GetBinContent(2);
    double x1 = hh1->GetBinCenter(1);
    double x2 = hh1->GetBinCenter(2);
    pdf = (y1-y2)/(x1-x2) * ( val- x1) + y1; 

    binwidth =  0.5 *( hh1->GetBinWidth(1) +  hh1->GetBinWidth(2));
    
  }else if( bin== nbins) {
    double y1 = hh1->GetBinContent(nbins-1);
    double y2 = hh1->GetBinContent(nbins);
    double x1 = hh1->GetBinCenter(nbins-1);
    double x2 = hh1->GetBinCenter(nbins);
    pdf = (y1-y2)/(x1-x2) * ( val- x1) + y1; 

    binwidth =  0.5 *( hh1->GetBinWidth(nbins-1) +  hh1->GetBinWidth(nbins));
  }
  else{
    if(val <=  hh1->GetBinCenter(bin)){
      double y1 = hh1->GetBinContent(bin-1);
      double y2 = hh1->GetBinContent(bin);
      double x1 = hh1->GetBinCenter(bin-1);
      double x2 = hh1->GetBinCenter(bin);
      pdf = (y1-y2)/(x1-x2) * ( val- x1) + y1; 
      
      binwidth =  0.5 *( hh1->GetBinWidth(bin) +  hh1->GetBinWidth(bin-1));

    }else{
      double y1 = hh1->GetBinContent(bin);
      double y2 = hh1->GetBinContent(bin+1);
      double x1 = hh1->GetBinCenter(bin);
      double x2 = hh1->GetBinCenter(bin+1);
      pdf = (y1-y2)/(x1-x2) * ( val- x1) + y1; 

      binwidth =  0.5 *( hh1->GetBinWidth(bin) +  hh1->GetBinWidth(bin+1));
    }
  }
  pdf /= binwidth;
  
  return pdf;
  
  
}



int Ntrial = 100; ///12 secs for testcat0 with pregenerated gaussian randoms


int Nrand = 2000000 *200;

vector<float> vrandgaus; 
void generateGaussRandom(){

  TRandom3 *gr = new TRandom3(12345);
  cout<<"generating gaus rand for use "<< Nrand<<endl; 
  for(int j=0;j<Nrand;j++){
    if(j%100000000==0) cout<<" j " << j <<endl; 
    float tmp = gr->Gaus();
    vrandgaus.push_back(tmp);
  }
  
}


double function2(double par[], bool drawHist = true)
{ 
  
  double logL = 0; 
  string sname; 
  
  //int nscpair = 1;
  int nscpair = 10;
    
  for(int indpair = 0; indpair < nscpair ; indpair++){
    sname = string(Form("mpairmcsmeared_indscpair%d",indpair));
    TString tname = TString("th1f_") +  TString(Form("mpairmcsmeared_indscpair%d",indpair));
    th1f_map[sname] = new TH1F(tname,tname, int( (fitrangeHigh-fitrangeLow)/binwidth+0.1),fitrangeLow,fitrangeHigh);
    th1f_map[sname]->Sumw2();
    
  }
  
  NMC = int(vsc1cat.size());

  //cout<<" NMC " << NMC <<endl; 
  TDatime *t1 = new TDatime();
  t1->Print();
  
  int nn = 0; 
  for(int j=0; j< NMC; j++){
    
    //if(j%500000 ==0) cout<<" j" << j <<endl;
    
    //     int sc1cat = vsc1cat[j];
    //     int sc2cat = vsc2cat[j];
    //     if(sc2cat < sc1cat){
    //       exChangeTwoNumber(sc1cat,sc2cat);
    //     }
    ///int indpair = map_indscpairBarrel[sc1cat][sc2cat]; 
    int indpair = vindpair[j];
    
    //test one category
    if(testOnlyCat >=0 && indpair != testOnlyCat) continue;
    
    int sc1cat = vsc1cat[j];
    int sc2cat = vsc2cat[j];
    

    float e1true = vele1etrue[j]; 
    float e2true = vele2etrue[j]; 
    
    eleen[0] = vele1en[j];
    eleen[1] = vele2en[j];
    eleeta[0] = vele1eta[j];
    eleeta[1] = vele2eta[j];
    eleieta[0] = vele1ieta[j];
    eleieta[1] = vele2ieta[j];
    eleiphi[0] = vele1iphi[j];
    eleiphi[1] = vele2iphi[j];

    float e1new;
    float e2new;
    
    sname = string(Form("mpairmcsmeared_indscpair%d",indpair));
    
    weight = vweight[j];
    mpair = vmpair[j];
    if(smearMethod=="corrSmear"){
      e1new = e1true *  ( ( par[sc1cat] * ( eleen[0]/e1true -1) + 1) + par[sc1cat+4]);
      e2new = e2true *  ( ( par[sc2cat] * ( eleen[1]/e2true -1) + 1) + par[sc2cat+4]);
      float s1 = e1new / eleen[0];
      float s2 = e2new / eleen[1];
      float mpair_new = mpair *sqrt( s1 * s2 );
      if( mpair_new < fitrangeLow || mpair_new > fitrangeHigh) continue; 
      if( e1new * sin(2*atan(exp(-eleeta[0]))) < etcut) continue; 
      if( e2new * sin(2*atan(exp(-eleeta[1]))) < etcut) continue; 
      th1f_map[sname]->Fill(mpair_new,weight);
      
    }else  if(smearMethod=="uncorrSmear"){
      
      for(int n=0;n<Ntrial;n++){
	if( nn < Nrand){
	  float rand = vrandgaus[nn];
	  e1new = eleen[0] + rand * par[sc1cat] * eleen[0];
	  rand = vrandgaus[nn+1];
	  e2new = eleen[1] + rand * par[sc2cat] * eleen[1];
	  nn++;
	}else{
	  cout<<"warning Nrand reached " << nn <<endl;
	  e1new = getGaussian(eleen[0],par[sc1cat]*eleen[0],0,10E3,eleieta[0],eleiphi[0],n);
	  e2new = getGaussian(eleen[1],par[sc2cat]*eleen[1],0,10E3,eleieta[1],eleiphi[1],n);
	}
	
	//scale
	e1new *= (1+ par[sc1cat+4]);
	e2new *= (1+ par[sc2cat+4]);
	
	
	float s1 = e1new / eleen[0];
	float s2 = e2new / eleen[1];
	float mpair_new = mpair *sqrt( s1 * s2 );
	if( mpair_new < fitrangeLow || mpair_new > fitrangeHigh) continue; 
	if( e1new * sin(2*atan(exp(-eleeta[0]))) < etcut) continue; 
	if( e2new * sin(2*atan(exp(-eleeta[1]))) < etcut) continue; 
	th1f_map[sname]->Fill(mpair_new,weight);
      }
      
    }else{
      cout<<"smearMethod  " << smearMethod.c_str()<<" NA " <<endl; 
      exit(1);
    }
    
    //cout<<"mapir_new " << par[sc1cat]<<" "<< par[sc2cat] <<" "<< s1 <<" "<< s2<<" "<< mpair_new <<" "<<mpair<<endl; 
    
  }
  
  
  for(int indpair = 0; indpair < nscpair ; indpair++){
    
    //test one category
    if(testOnlyCat >=0 && indpair != testOnlyCat) continue;
    
    sname = string(Form("mpairmcsmeared_indscpair%d",indpair));
    float mc =  th1f_map[sname] ->Integral();
    string snamed = string(Form("mpairdata_indscpair%d",indpair));
    float data =  th1f_map[snamed] ->Integral();
    if(mc<=0 || data <=0){
      cout<<"emtpy categories ? " << indpair <<" "<<mc<<" "<<data <<endl; 
    }
    
    th1f_map[sname] ->Scale(data/mc);
        
    if(fitmethod=="lh"){
      
    }else if( fitmethod== "lhpoisson"){
      ///      logL = 0; 
      int nbins =  th1f_map[sname]->GetNbinsX();
      for(int b=1; b<= nbins; b++){
	double n = th1f_map[snamed]->GetBinContent(b);
	double m = th1f_map[sname]->GetBinContent(b);
	double p = 0; 
	if( m <=0) continue; 
	p = n * log(m) -m ; 
	for(int k=n; k>=1; k--){
	  p -= log(k);
	}
	logL += p;
      }
    }
    
    if(drawHist){
      //draw
      TCanvas *can0 = new TCanvas("can0","c000",200,10,550,500);
      
      float ymax = th1f_map[sname]->GetMaximum();
      th1f_map[snamed]->GetYaxis()->SetRangeUser(1,1.1*ymax);
      th1f_map[snamed]->Draw("e");
      th1f_map[sname]->Draw("histsame");  
      
      double p = th1f_map[snamed]->Chi2Test(th1f_map[sname],"UW P");
      cout<<"chi2test p " << p <<endl; 

      double chi2 = th1f_map[snamed]->Chi2Test(th1f_map[sname],"UW CHI2");
      double chi2ndf = th1f_map[snamed]->Chi2Test(th1f_map[sname],"UW CHI2/NDF");
      
      TString chi2res = TString(Form("#chi^{2}/ndf = %4.2f/%d",chi2,int(chi2/chi2ndf+0.1)));
      TLatex *   tex = new TLatex(0.62,0.8,chi2res);
      tex->SetNDC();
      tex->SetTextSize(0.05);
      tex->Draw();
      TString gifname = TString(Form("th1f_mpairs_data_indpair%d_%4.4f.pdf",indpair,par[indpair]));
      can0->Print(gifname);
      gifname = TString(Form("th1f_mpairs_data_indpair%d_%4.4f.C",indpair,par[indpair]));
      can0->Print(gifname);
    }
    
  }
  for(int indpair = 0; indpair < nscpair ; indpair++){
    ///delete histograms
    sname = string(Form("mpairmcsmeared_indscpair%d",indpair));
    th1f_map[sname] ->Delete();
  }
  
  
  double f = -logL; 
  
  cout<<" par"; 
  for(int j=0; j< 8; j++){
    cout<<" "<< par[j]; 
  }
  cout<<" f2: "<< f <<endl; 
  
  //TDatime *t2 = new TDatime();
  //t2->Print();
  
  
  return f;
  
  
  
}


void function2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{ 
  f = function2(par,false);
  return ;
}





//void testZeeMCSmear(char *var, char *det,float xmin, float xmax, float xmin_plot, float xmax_plot, string fitmethod = "chisq",int test_transerfun=0, int nbins_fit=200, int nbins_plot=200,float ainit =1.0, bool doFit=true){

void testZeeMCSmearWithScale(string datareco = "jan16 or nov30",int phtcorr = 231, int test_fitdet= 1, string test_smearmethod = "corrSmear or uncorrSmear",string test_fitmethod = "lh or lhpoisson",int fixscale = 0,int fitstra = 1, int testpaircat=0, double test_fitlow = 70, double test_fithigh = 110, double test_binwidth = 1.0){
  
  
  

  int phtcorr_test =  phtcorr; 
  if( phtcorr_test == 486){
    phtcorr = 485;
  }


  fitdet = test_fitdet; 
  
  smearMethod = test_smearmethod ;
  
  testOnlyCat = testpaircat; 
  
  fitrangeLow = test_fitlow; 
  fitrangeHigh = test_fithigh; 
  binwidth = test_binwidth;
  
  
  fitmethod = test_fitmethod; 
  
  
  gStyle->SetOptStat(0);
  gStyle->SetNdivisions(508,"X");
  gStyle->SetNdivisions(512,"Y");
    
  
  int mc = 1; 

  fillMapIndex();
  

  TString filename; 
  
  
    
  fChainMC = new TChain("Analysis");
  fChainData = new TChain("Analysis");
  
  
  
  string dataset; 
  int iflag = 2; 
  if( datareco == "nov30"){
    iflag = 1; 
    dataset = "DoubleElectronRun2011AB30Nov2011v1AOD";
  }else if( datareco == "jan16"){
    iflag = 2; 
    dataset = "DoubleElectronRun2011AB16Jan2012v1AOD"; 
  }else{
    cout<<"dataset ? " << datareco.c_str()<<endl; 
    return; 
  }
  

  ///MC

  filename = TString(Form("/home/raid2/yangyong/data/CMSSW/v1/CMSSW_4_2_8/src/zSelector/dielectrontree/makeDiElectronTree.v3.DYJetsToLL_TuneZ2_M50_7TeVmadgraphtauolaFall11PU_S6_START42_V14Bv1AODSIMallstat.etcut20.corr%d.eleid1.datapu6.mcpu1.r1to92.scale0.root",phtcorr));
  
  cout<<filename<<endl; 
  fChainMC->Add(filename);
  
  //fChainMC->SetBranchAddress("isRealData", &isRealData);
  fChainMC->SetBranchAddress("elescr9", elescr9);
  fChainMC->SetBranchAddress("elesceta", elesceta);
  fChainMC->SetBranchAddress("eleen", eleen);
  
  //fChainMC->SetBranchAddress("eleen0", eleen0);
  fChainMC->SetBranchAddress("mpair", &mpair);
//   fChainMC->SetBranchAddress("mpair0", &mpair0);
  fChainMC->SetBranchAddress("eleeta", eleeta);
   fChainMC->SetBranchAddress("eleieta", eleieta);
  fChainMC->SetBranchAddress("eleiphi", eleiphi);
   fChainMC->SetBranchAddress("runNumber", &runNumber);
   fChainMC->SetBranchAddress("evtNumber", &evtNumber);
   fChainMC->SetBranchAddress("lumiBlock", &lumiBlock);
  fChainMC->SetBranchAddress("eleetrue", eleetrue);
//   fChainMC->SetBranchAddress("eleetruenofsr", eleetruenofsr);
//   fChainMC->SetBranchAddress("eleetrueplusfsrdrbc003", eleetrueplusfsrdrbc003);
//   fChainMC->SetBranchAddress("eleetrueplusfsrdrbc004", eleetrueplusfsrdrbc004);
   fChainMC->SetBranchAddress("eleetrueplusfsrdrbc005", eleetrueplusfsrdrbc005);
//   fChainMC->SetBranchAddress("eleetrueplusfsrdrbc006", eleetrueplusfsrdrbc006);
//   fChainMC->SetBranchAddress("eleetrueplusfsrdrbc007", eleetrueplusfsrdrbc007);
//   fChainMC->SetBranchAddress("eleetrueplusfsrdrbc008", eleetrueplusfsrdrbc008);
//   fChainMC->SetBranchAddress("eleetrueplusfsrdrbc009", eleetrueplusfsrdrbc009);
   fChainMC->SetBranchAddress("eleetrueplusfsrdrbc01", eleetrueplusfsrdrbc01);
//   fChainMC->SetBranchAddress("eletrueelematched", eletrueelematched);
//   fChainMC->SetBranchAddress("eletrueelematchednofsr", eletrueelematchednofsr);
   fChainMC->SetBranchAddress("weight", &weight);
   
  //Data
   filename = TString(Form("/home/raid2/yangyong/data/CMSSW/v1/CMSSW_4_2_8/src/zSelector/dielectrontree/makeDiElectronTree.v3.%s.etcut20.corr%d.eleid1.datapu0.mcpu0.r1to131.scale0.root",dataset.c_str(),phtcorr));
   
   cout<<filename<<endl; 
  fChainData->Add(filename);
  
  
  fChainData->SetBranchAddress("isRealData", &isRealData);
  fChainData->SetBranchAddress("elescr9", elescr9);
  fChainData->SetBranchAddress("elesceta", elesceta);
  fChainData->SetBranchAddress("eleen", eleen);
  fChainData->SetBranchAddress("eleen0", eleen0);
  fChainData->SetBranchAddress("mpair", &mpair);
  fChainData->SetBranchAddress("mpair0", &mpair0);
  fChainData->SetBranchAddress("eleeta", eleeta);
  fChainData->SetBranchAddress("eleieta", eleieta);
  fChainData->SetBranchAddress("eleiphi", eleiphi);
  fChainData->SetBranchAddress("runNumber", &runNumber);
  fChainData->SetBranchAddress("evtNumber", &evtNumber);
  fChainData->SetBranchAddress("lumiBlock", &lumiBlock);
  
  
  NMC = fChainMC->GetEntries();
  NData = fChainData->GetEntries();
  cout<<"nMC " << NMC <<" "<< NData <<endl; 
    
  
  filename = TString(Form("testZeeMCSmearWithScale.%s.%s.%s.testpair%d.fitrange%dto%d.binwidth%2.2f.fitdet%d.corr%d.fitstra%d.fixscale%d.etcut%d.root",datareco.c_str(),test_smearmethod.c_str(),test_fitmethod.c_str(),testpaircat,int(fitrangeLow+0.1),int(fitrangeHigh+0.1),binwidth,fitdet,phtcorr_test,fitstra,fixscale,etcut));
  TFile *fnew = new TFile(filename,"recreate");
  
  //filename = TString(Form("testZeeMCSmearWithScale.%s.%s.%s.testpair%d.fitrange%dto%d.binwidth%2.2f.fitdet%d.corr%d.fitstra%d.fixscale%d.etcut%d.txt",datareco.c_str(),test_smearmethod.c_str(),test_fitmethod.c_str(),testpaircat,int(fitrangeLow+0.1),int(fitrangeHigh+0.1),binwidth,fitdet,phtcorr,fitstra,fixscale,etcut));
  //txtout.open(filename,ios::out);
  

  TTree *mytree = new TTree("Analysis","");
  
  float mcfitpar[8] = {0};
  float mcfitparErr[8] = {0};
  int mcfitStatus;
  float fAmin;
  mytree->Branch("mcfitpar",mcfitpar,"mcfitpar[8]/F");
  mytree->Branch("mcfitparErr",mcfitparErr,"mcfitparErr[8]/F");
  mytree->Branch("mcfitStatus",&mcfitStatus,"mcfitStatus/I");
  mytree->Branch("fAmin",&fAmin,"fAmin/F");
  
  rv_mass = new RooRealVar("rv_mass","mass",100,0,1E6);
  rv_weight = new RooRealVar("rv_weight","weight",1.0,0,1E6);
  
  for(int j=0;j<20;j++){
    string sname = string(Form("mpairmc_indscpair%d",j));
    makeRootDataSetAndTH1F(sname,int( (fitrangeHigh-fitrangeLow)/binwidth+0.1),fitrangeLow,fitrangeHigh);
    sname = string(Form("mpairdata_indscpair%d",j));
    makeRootDataSetAndTH1F(sname,int( (fitrangeHigh-fitrangeLow)/binwidth+0.1),fitrangeLow,fitrangeHigh);
  }
  
  int sc1cat; 
  int sc2cat; 
  
  
  cout<<" nMC " <<NMC <<endl; 
  for(int n=0; n< NMC; n++){
    fChainMC->GetEntry(n);
    if(phtcorr==94){
      if(evtNumber%2==0) continue;
    }
    
    ///same as 485, for test purpose only
    if(phtcorr_test==486){
      if(evtNumber%2==0) continue;
    }
    
    

    if(n%500000==0) cout<<"n " << n <<endl;
    
    if(fitdet==1){
      if( fabs(elesceta[0])>1.48 || fabs(elesceta[1])>1.48) continue; 
    }else if( fitdet==2){
      if( fabs(elesceta[0])<1.48 || fabs(elesceta[1])<1.48) continue; 
    }
    
    
    ///scale r9 in MC
    for(int j=0; j<2; j++){
      if(fabs(elesceta[j])<1.48) elescr9[j] *= 1.004;
      else elescr9[j] *= 1.006;
    }
    
    
    ///get pair index
    sc1cat = scCategoryEight(elesceta[0],elescr9[0]);
    sc2cat = scCategoryEight(elesceta[1],elescr9[1]);
    if(sc2cat < sc1cat){
      exChangeTwoNumber(sc1cat,sc2cat);
    }
    int indpair; 
    if(fitdet==1){
      indpair = map_indscpairBarrel[sc1cat][sc2cat]; 
    }else{
      indpair = map_indscpairEndcap[sc1cat][sc2cat]; 
    }
    
    
    //test one category
    if(testOnlyCat >=0 && indpair != testOnlyCat) continue;
    float e1true;
    if( fabs(elesceta[0])<1.48){
      e1true = eleetrueplusfsrdrbc005[0];
    }else{
      e1true = eleetrueplusfsrdrbc01[0];
    }
    float e2true;
    if( fabs(elesceta[1])<1.48){
      e2true = eleetrueplusfsrdrbc005[1];
    }else{
      e2true = eleetrueplusfsrdrbc01[1];
    }
    if(e1true <1 || e2true<1) continue; 
    
    vindpair.push_back(indpair);
    sc1cat = scCategoryEight(elesceta[0],elescr9[0]);
    sc2cat = scCategoryEight(elesceta[1],elescr9[1]);

    
    if(fitdet==2){ //Endcap, starting with 4,5,6,7
      sc1cat -=4; 
      sc2cat -=4; 
    }
    
    vsc1cat.push_back(sc1cat);
    vsc2cat.push_back(sc2cat);
    vele1eta.push_back(eleeta[0]);
    vele1ieta.push_back(eleieta[0]);
    vele1iphi.push_back(eleiphi[0]);
    vele2eta.push_back(eleeta[1]);
    vele2ieta.push_back(eleieta[1]);
    vele2iphi.push_back(eleiphi[1]);
    vele1en.push_back(eleen[0]);
    vele2en.push_back(eleen[1]);
    vmpair.push_back(mpair);
    vweight.push_back(weight);
    
    
    vele1etrue.push_back(e1true);
    vele2etrue.push_back(e2true);
    
    
    if( eleen[0]* sin(2*atan(exp(-eleeta[0]))) < etcut) continue; 
    if( eleen[1]* sin(2*atan(exp(-eleeta[1]))) < etcut) continue; 
    
    string sname = string(Form("mpairmc_indscpair%d",indpair));
    fillRootDataSetAndTH1F(sname,mpair,weight);
    
  }
  
  
  setMap_cat_runbin();
  loadRunbyRunScaleCorrectionNoR9(iflag,phtcorr);
  loadScaleCorrectionR9(iflag,phtcorr);
  
  
  for(int n=0; n< NData; n++){
    fChainData->GetEntry(n);
    if(n%500000==0) cout<<"n " << n <<endl;

    if(fitdet==1){
      if( fabs(elesceta[0])>1.48 || fabs(elesceta[1])>1.48) continue; 
    }else if( fitdet==2){
      if( fabs(elesceta[0])<1.48 || fabs(elesceta[1])<1.48) continue; 
    }
    
    float scorr = 1.0; 
    for(int j=0; j<2; j++){
      float dEoE = energyScaleRunNoR9(elesceta[j]);
      float dEoE1 = energyScaleR9(elesceta[j],elescr9[j]);
      eleen[j] *= 1.0/(1+dEoE) * 1.0/(1+dEoE1);
      scorr *= 1.0/(1+dEoE) * 1.0/(1+dEoE1);
    }
    mpair *= sqrt(scorr);
    if( eleen[0]*sin(2*atan(exp(-eleeta[0]))) < etcut) continue; 
    if( eleen[1]*sin(2*atan(exp(-eleeta[1]))) < etcut) continue; 
    
    sc1cat = scCategoryEight(elesceta[0],elescr9[0]);
    sc2cat = scCategoryEight(elesceta[1],elescr9[1]);
    if(sc2cat < sc1cat){
      exChangeTwoNumber(sc1cat,sc2cat);
    }
    int indpair; 
    if(fitdet==1){
      indpair = map_indscpairBarrel[sc1cat][sc2cat]; 
    }else{
      indpair = map_indscpairEndcap[sc1cat][sc2cat]; 
    }
    
    
    string sname = string(Form("mpairdata_indscpair%d",indpair));
    fillRootDataSetAndTH1F(sname,mpair,1);
  }
  
  
  for(int j=0;j<10;j++){
    string sname = string(Form("mpairmc_indscpair%d",j));
    string snamed = string(Form("mpairdata_indscpair%d",j));
    cout<<"mpair catpair "     <<  th1f_map[sname]->GetEntries()<<" "<<  th1f_map[snamed]->GetEntries()<<endl; 
  }
  
  generateGaussRandom();
  
  //return; 

  TMinuit *minuit; 
  int npar = 8; 
  minuit  = new TMinuit(npar); 
  
  //minuit->SetFCN(function);
  //minuit->SetFCN(function1);
  minuit->SetFCN(function2);
  
  
  //settings
  Double_t arglist[1];
  Int_t ierflg = 0;
  //double STEPMN = 0.01;
  double STEPMN = 0.0001;
  
  // 1 for Chi square
  // 0.5 for negative log likelihood
  if(fitmethod== "lh" || fitmethod == "lhpoisson"){
    minuit->SetErrorDef(0.5);
  }else{
    minuit->SetErrorDef(1);
  }
    double fitpar[10];
    double fitparErr[10];
    
  double smearcat[4] = {1.1,1.1,1.1,1.1};
  double smearcatMax[4] = {3,3,3,3};

  double deltaEcat[4] = {0,0,0,0};
  
  
  if(smearMethod=="uncorrSmear"){
    smearcatMax[0] = 0.05;
    smearcatMax[1] = 0.05;
    smearcatMax[2] = 0.05;
    smearcatMax[3] = 0.05;

    smearcat[0] = 0.007;
    smearcat[1] = 0.01;
    smearcat[2] = 0.015;
    smearcat[3] = 0.02;
    
    if(fitdet==2){
      smearcat[0] = 0.02;
      smearcat[1] = 0.02;
      smearcat[2] = 0.02;
      smearcat[3] = 0.02;
    }
    
  }
  

  for(int j=0; j<4; j++){
    TString parname = TString (Form("smearcat%d",j));
    if(smearMethod=="corrSmear"){
      minuit->mnparm(j, parname, smearcat[j], STEPMN, 0.5,smearcatMax[j],ierflg);
    }else  if(smearMethod=="uncorrSmear"){
      minuit->mnparm(j, parname, smearcat[j], STEPMN, 0,smearcatMax[j],ierflg);
    }

    fitpar[j] = smearcat[j]; //initialized values
  }
  
  for(int j=4; j<npar; j++){
    TString parname = TString (Form("scalecat%d",j-4));
    minuit->mnparm(j, parname, deltaEcat[j-4], STEPMN, -0.03,0.03,ierflg);
    fitpar[j] = deltaEcat[j-4]; //initialized values
  }
  

  if(fixscale==1){
    for(int j=4; j<npar; j++){
      minuit->FixParameter(j);
    }
  }
  


  if(testpaircat==1){
    minuit->FixParameter(0);
    minuit->FixParameter(2);
    minuit->FixParameter(3);

    minuit->FixParameter(4);
    minuit->FixParameter(6);
    minuit->FixParameter(7);

  }else if( testpaircat==0){
    minuit->FixParameter(1);
    minuit->FixParameter(2);
    minuit->FixParameter(3);

    minuit->FixParameter(5);
    minuit->FixParameter(6);
    minuit->FixParameter(7);

  } else if( testpaircat==2){
    minuit->FixParameter(0);
    minuit->FixParameter(1);
    minuit->FixParameter(3);

    minuit->FixParameter(4);
    minuit->FixParameter(5);
    minuit->FixParameter(7);

  } else if( testpaircat==3){
    minuit->FixParameter(0);
    minuit->FixParameter(1);
    minuit->FixParameter(2);

    minuit->FixParameter(4);
    minuit->FixParameter(5);
    minuit->FixParameter(6);
    
  }  
  
  arglist[0] = 0.0001;
  minuit->mnexcm("SET EPS",arglist,1,ierflg);
  
  //minuit->mnsimp();
  ////arglist[0] = 1; 
  ///minuit->mnexcm("SET STR",arglist,1,ierflg);
  //minuit->mnexcm("MIGRAD", arglist ,1,ierflg);
  
  arglist[0] = fitstra; 
  minuit->mnexcm("SET STR",arglist,1,ierflg);

  bool dofit = true; 
  
  if( dofit){
    minuit->Migrad();
    if (!minuit->fCstatu.Contains("CONVERGED")) {
      mcfitStatus = 1; //first try not converged.
      minuit->Migrad();
    }else{
      mcfitStatus = 0;
      
    }
  }
    
  double ftest[1000]; 

  double atest[1000];
  
  //minuit->GetParameter(0,fitpar[0],fitparErr[0]);

  if (!minuit->fCstatu.Contains("CONVERGED")) {
    printf("No convergence at fitting, routine Fit \n");
    printf("Minuit return string: %s \n",minuit->fCstatu.Data());
    mcfitStatus = 2;  //2nd try not converged.
  }else{
    mcfitStatus = -1;
  }
  fAmin = minuit->fAmin;

  for(int j=0; j<npar; j++){
    minuit->GetParameter(j,fitpar[j],fitparErr[j]);
    mcfitpar[j] = fitpar[j];
    mcfitparErr[j] = fitparErr[j];
  }
  
  
  double par[10];
  
  for(int j=0;j<npar;j++){
    par[j] = fitpar[j];
  }
  
  
  bool printScan = false; 
  //bool printScan = true; 
  par[0] = 0.0076; 
  
  if(printScan && testpaircat>=0){
    double  fmin = function2(par,true);
    
    
    float x1 = fitpar[testpaircat] - 0.1;
    float x2 = fitpar[testpaircat] + 0.1;
    float teststep = 0.01;
    if(smearMethod=="uncorrSmear"){
      x1 =  fitpar[testpaircat] - 0.005 >0 ? fitpar[testpaircat] - 0.005: 0;
      x2 =  fitpar[testpaircat] +0.005;
      teststep = 0.001;
    }
    
    double dmin1 = 1;
    double dmin2 = 1;
    double fitparErrLow[10];
    double fitparErrHigh[10];
  
    int ntest = 0; 
    for(float x = x1; x <= x2; x += teststep){
      
      cout<<"x " << x <<endl; 
      
      
      par[testpaircat]  = x ;
    
      //double  f = function2(par,false);
      double  f = function2(par,true);
    
      if(fabs(f-fmin-1)< dmin1 && x < fitpar[testpaircat]){
	dmin1 = fabs(f-fmin-0.5);
	fitparErrLow[testpaircat] = fitpar[testpaircat] -x;
      }
      if(fabs(f-fmin-1)< dmin2 && x > fitpar[testpaircat]){
	dmin2 = fabs(f-fmin-0.5);
	fitparErrHigh[testpaircat] = x-fitpar[testpaircat];
      }
      atest[ntest] = x; 
      ftest[ntest] = f; 
      ntest ++; 
    }
    
    cout<<"fitpar[testpaircat]" << fitpar[testpaircat]<<" +/- "<<fitparErr[testpaircat]<<" - "<< fitparErrLow[testpaircat]<<" + "<< fitparErrHigh[testpaircat]<<endl; 
  
    TCanvas *can0 = new TCanvas("can0","c000",200,10,550,500);
    setTCanvasNice(can0);
    TGraph *gr = new TGraph(ntest,atest,ftest);
    gr->Draw("ap");
    //TF1 *ff = new TF1("ff","[0]+[1]*x+[2]*x*x",0,0.0);
    //gr->Fit(ff,"r");
  
    can0->Print("testconverge1.pdf");
    can0->Print("testconverge1.C");
    
  }
  
  mytree->Fill();
  mytree->Write();
  fnew->Write();
  fnew->Close();


}
