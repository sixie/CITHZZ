#include "rootheader.h"
TChain *fChainData;
TChain *fChainMC;

#include "zeevariables.h"
#include "setbranches.cc"

#include "TMinuit.h"

#include "utils.cc"
#include "energyScaleCorrection.cc"


vector<int> v_sc1cat;
vector<int> v_sc2cat;
vector<int> v_mcat;
vector<float> v_mee;

vector<float> v_e1;
vector<float> v_e2;
vector<float> v_wt;

vector<float> v_rand1;
vector<float> v_rand2;
vector<float> v_e1eta;
vector<float> v_e2eta;


TH1F *th1f_Mee_data[10];
TH1F *th1f_Mee_mc[10];
TH1F *th1f_Mee_mc0[10];

double par_fit[10];

  
int test_cat; 

TRandom3 *grand; 


int nbins = 80;
// float xmin = 70;
// float xmax = 110;

int Ntrial = 50; 


double fitrangeLow = 75;
double fitrangeHigh = 105;
double binwidth; 

long int Nrand; 
int eleID;


void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  
 
  
  int Nmc =  int(v_e1.size());
    
  
  for(int j=0; j<10; j++){
    TString histname = TString(Form("th1f_Mee_mc_%d",j));
    th1f_Mee_mc[j] =new TH1F(histname,histname,nbins,fitrangeLow,fitrangeHigh);
  }
  
  int nrand = 0; 
  float rand1; 
  float rand2; 

  for(int j=0; j< Nmc; j++){
    float e1 = v_e1[j];
    float e2 = v_e2[j];
    int sc1cat = v_sc1cat[j];
    int sc2cat = v_sc2cat[j];
    int mcat = v_mcat[j];

    if( test_cat >=0 && mcat != test_cat) continue; 
    
    float mee = v_mee[j];

  
    float wt = v_wt[j];

    float e1eta = v_e1eta[j];
    float e2eta = v_e2eta[j];
    
    ///cout<<"sccat1 " << sc1cat <<" "<< sc2cat <<endl; 
    

    for(int n=0; n<Ntrial; n++){
      if( nrand < Nrand ){
	rand1 = v_rand1[nrand];
	rand2 = v_rand2[nrand];
      }else{
	cout<<"warning Nrand reached " << nrand <<endl; //this should not happen..
	rand1  = gRandom->Gaus();
	rand2  = gRandom->Gaus();
      }
      
      nrand ++;
      
      float e1s = e1 * (1+ par[sc1cat+4]) * ( 1+ rand1 * par[sc1cat]);
      float e2s = e2 * (1+ par[sc2cat+4]) * ( 1+ rand2 * par[sc2cat]);
      
      
      float et1 = e1s * sin(2*atan(exp(-e1eta)));
      float et2 = e2s * sin(2*atan(exp(-e2eta)));
      
      
      if(eleID==1){
	if(et1<15 || et2 <15) continue; 
      }
      if(eleID==2){
	if(et1<20 || et2 <20) continue; 
      }
      
      float mees = sqrt( (e1s*e2s)/(e1*e2)) * mee; 
      th1f_Mee_mc[mcat]->Fill(mees,wt);
    }
    
  }

  //cout<<"intergral "<<th1f_Mee_mc[0]->GetEntries()<<" " <<th1f_Mee_mc[0]->Integral()<<endl; 
  //cout<<"rms "<< th1f_Mee_mc[0]->GetRMS()<<" "<< th1f_Mee_mc[0]->GetMean()<<endl; 
  
  
  double logL = 0; 
  for(int icat =0; icat<10; icat++){
    double tmp1 = th1f_Mee_data[icat]->Integral();
    double tmp2 = th1f_Mee_mc[icat]->Integral();
    th1f_Mee_mc[icat]->Scale(tmp1/tmp2);
    
    if( test_cat >=0 && icat != test_cat) continue; 
    
    for(int b=1; b<= nbins; b++){
      double data = th1f_Mee_data[icat]->GetBinContent(b);
      double mc = th1f_Mee_mc[icat]->GetBinContent(b);

      //      logL += data - mc; 
      // if(data >0 && mc > 0){   ///likelhood ratio
      //	logL += -data *log(data/mc); 
      // }

      if( mc <=0) continue; 
      double p = data * log(mc) -mc ;  ///pure likehilhood
      for(int k=data; k>=1; k--){ 
	p -= log(k);
      }
      logL += p;
      
      
    }
    
  }
  
  for(int j=0;j<10;j++){
    th1f_Mee_mc[j]->Delete();
  }
  
  f= -2*logL;
  
  cout<<"f " << f << " "<< par[0]<<" "<<par[1]<<" "<<  par[2]<<" "<<par[3]<<endl; 
  
}




void function(double par[]){
  int npar; 
  double *gin; 
  //double *par;
  double f; 
  int iflag; 
  //fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  fcn(npar,gin,f,par,iflag);

  cout<<"f: " << f <<endl; 
  
}


void testEnergySmearZee(int barrelorEndcap=1, int energyRegrVer = 0, int test_eleID = 2, int testcat=-1 ,int fixscale = 1, 
			double test_fitlow = 70, double test_fithigh = 110, double test_binwidth = 0.5, int test_Ntrial=50){
  
  Ntrial = test_Ntrial; 
  eleID = test_eleID; 
  test_cat = testcat; 
  fitrangeLow = test_fitlow; 
  fitrangeHigh = test_fithigh; 
  binwidth = test_binwidth;
  
  TString inputdatafile = "/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_HCP2012.root";
  TString inputmcfile = "/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Summer12DY53X.root";

  
  TString filename = TString(Form("testEnergySmearZee_barend%d_regver%d_eleID%d_fitrange%dto%d_binwidth%2.2f_Ntrial%d.root",barrelorEndcap,energyRegrVer,eleID,
				  int(fitrangeLow+0.1),int(fitrangeHigh+0.1),binwidth, Ntrial));
  
  TFile *fnew = new TFile(filename,"recreate");
  

  nbins = int( (fitrangeHigh-fitrangeLow)/binwidth+0.1); 

  for(int j=0; j<10; j++){
    TString histname = TString(Form("th1f_Mee_data_%d",j));
    th1f_Mee_data[j] =new TH1F(histname,histname,nbins,fitrangeLow,fitrangeHigh);
    histname = TString(Form("th1f_Mee_mc0_%d",j));
    th1f_Mee_mc0[j] =new TH1F(histname,histname,nbins,fitrangeLow,fitrangeHigh);
  }
  
  fChainData = new TChain("ZeeEvent");
  fChainMC = new TChain("ZeeEvent");

  //data
  fChainData->Add(inputdatafile);
  fChainMC->Add(inputmcfile);
  
  setbranchaddress();
  
  
  int totalEntriesMC = fChainMC->GetEntries();
  int totalEntriesData = fChainData->GetEntries();
  

  cout<<" totalEntriesMC /data " << totalEntriesMC << " "<< totalEntriesData <<endl; 

  //totalEntriesData = 1E6;
  
  load_energyScaleCorrection_eta_RunByRun();
  load_energyScaleCorrection_etaR9bin();
  

  cout<< " fill data histograms " <<endl; 
  for(int entry=0; entry< totalEntriesData; entry++){
    fChainData->GetEntry(entry);
    if(entry%100000==0) cout<<"entry " << entry <<endl; 
    
    if(barrelorEndcap==1 && (fabs(Ele1SCEta)>1.49 || fabs(Ele2SCEta)>1.49)) continue; 
    if(barrelorEndcap==2 && (fabs(Ele1SCEta)<1.49 || fabs(Ele2SCEta)<1.49)) continue; 
    
    ///
    
    
    if(eleID==1 || eleID==2){
      if( Ele1PassMediumSimpleCuts ==0 || Ele2PassMediumSimpleCuts == 0) continue; 
    }
    
    //here missing run-by-run corrections. TBD
    

    int sccat1 = scCategoryEight(Ele1SCEta,Ele1R9);
    int sccat2 = scCategoryEight(Ele2SCEta,Ele2R9);
    
    if( barrelorEndcap==2){
      sccat1 -=4; 
      sccat2 -=4; 
    }
    
    int mcat = MeesmearCategory(sccat1,sccat2);

    if(testcat >=0 &&  testcat != mcat) continue; 
    
    
    float e1 = Ele1Energy;
    float e2 = Ele2Energy;

    if(energyRegrVer==0){
      e1 = Ele1EnergyRegressionV0;
      e2 = Ele2EnergyRegressionV0;
    } else if(energyRegrVer==1){
      e1 = Ele1EnergyRegressionV1;
      e2 = Ele2EnergyRegressionV1;
    } else if(energyRegrVer==2){
      e1 = Ele1EnergyRegressionV2;
      e2 = Ele2EnergyRegressionV2;
    }
    
    float corr1 = energyScaleCorrection_eta_RunByRun(Ele1SCEta,run);
    e1 *= ( 1- corr1);
    float corr2 = energyScaleCorrection_eta_RunByRun(Ele2SCEta,run);
    e2 *= ( 1- corr2);
    corr1 = energyScaleCorrection_etaR9bin(Ele1SCEta,Ele1R9);
    e1 *= (1-corr1);
    corr2 = energyScaleCorrection_etaR9bin(Ele2SCEta,Ele2R9);
    e2 *= (1-corr2);
    
    

    float et1 = e1 * sin(2*atan(exp(-Ele1Eta)));
    float et2 = e2 * sin(2*atan(exp(-Ele2Eta)));
    
    if(eleID==1){
      if(et1<15 || et2 <15) continue; 
    }
    if(eleID==2){
      if(et1<20 || et2 <20) continue; 
    }
    
    float mee = calcZmass(e1,Ele1Eta,Ele1Phi,e2,Ele2Eta,Ele2Phi);
    th1f_Mee_data[mcat]->Fill(mee);
    
  }
  
  
  //totalEntriesMC = 1E6;
  
  int nMC = 0; 

  ///MC
  for(int entry=0; entry< totalEntriesMC; entry++){
    fChainMC->GetEntry(entry);
    if(entry%100000==0) cout<<"entry " << entry <<endl; 
    
    if( event%2==0) continue; ///used for training.
    

    if(barrelorEndcap==1 && (fabs(Ele1SCEta)>1.49 || fabs(Ele2SCEta)>1.49)) continue; 
    if(barrelorEndcap==2 && (fabs(Ele1SCEta)<1.49 || fabs(Ele2SCEta)<1.49)) continue; 
    ///
    
    
    if(eleID==1 || eleID==2){
      if( Ele1PassMediumSimpleCuts ==0 || Ele2PassMediumSimpleCuts == 0) continue; 
    }
    
    
    int sccat1 = scCategoryEight(Ele1SCEta,Ele1R9);
    int sccat2 = scCategoryEight(Ele2SCEta,Ele2R9);
    if( barrelorEndcap==2){
      sccat1 -=4; 
      sccat2 -=4; 
    }
    
    int mcat = MeesmearCategory(sccat1,sccat2);

    if(testcat >=0 &&  testcat != mcat) continue; 
    
    nMC ++; 
    
    float e1 = Ele1Energy;
    float e2 = Ele2Energy;
    if(energyRegrVer==0){
      e1 = Ele1EnergyRegressionV0;
      e2 = Ele2EnergyRegressionV0;
    } else if(energyRegrVer==1){
      e1 = Ele1EnergyRegressionV1;
      e2 = Ele2EnergyRegressionV1;
    } else if(energyRegrVer==2){
      e1 = Ele1EnergyRegressionV2;
      e2 = Ele2EnergyRegressionV2;
    }
    
    float mee = calcZmass(e1,Ele1Eta,Ele1Phi,e2,Ele2Eta,Ele2Phi);
    v_e1.push_back(e1);
    v_e2.push_back(e2);
    v_sc1cat.push_back(sccat1);
    v_sc2cat.push_back(sccat2);
    v_mee.push_back(mee);
    v_mcat.push_back(mcat);
    v_wt.push_back(weight);
    v_e1eta.push_back(Ele1Eta);
    v_e2eta.push_back(Ele2Eta);


    float et1 = e1 * sin(2*atan(exp(-Ele1Eta)));
    float et2 = e2 * sin(2*atan(exp(-Ele2Eta)));

    
    if(eleID==1){
      if(et1<15 || et2 <15) continue; 
    }
    if(eleID==2){
      if(et1<20 || et2 <20) continue; 
    }
    th1f_Mee_mc0[mcat]->Fill(mee,weight);
    
  }
  
  cout<<"mc "<< th1f_Mee_mc0[0]->GetEntries()<<" " << th1f_Mee_mc0[0]->Integral()<<endl;
  
  Nrand = Ntrial * nMC;
  cout<<"filling random " << Ntrial * nMC <<endl; 
  
  grand = new TRandom3(12345);
  for(int j=0; j< Nrand; j++){
    v_rand1.push_back( grand->Gaus() );
    v_rand2.push_back( grand->Gaus() );
  }
  
 //  par_fit[0] = 0.005;
//   for(int j=0;j<10;j++){
//     function(par_fit);
//     par_fit[0] += 0.0005;
//   }
  

  TTree *mytree = new TTree("Analysis","");
  float mcfitpar[8] = {0};
  float mcfitparErr[8] = {0};
  int mcfitStatus;
  float fAmin;
  mytree->Branch("mcfitpar",mcfitpar,"mcfitpar[8]/F");
  mytree->Branch("mcfitparErr",mcfitparErr,"mcfitparErr[8]/F");
  mytree->Branch("mcfitStatus",&mcfitStatus,"mcfitStatus/I");
  mytree->Branch("fAmin",&fAmin,"fAmin/F");
  
   //return; 

  TMinuit *minuit; 
  int npar = 8; 
  minuit  = new TMinuit(npar); 

  minuit->SetFCN(fcn);
  
  //settings
  Double_t arglist[1];
  Int_t ierflg = 0;
  //double STEPMN = 0.01;
  double STEPMN = 0.0001;

  double fitpar[10];
  double fitparErr[10];
  
  double smearcat[4] = {0.02,0.02,0.02,0.02};
  double smearcatMax[4] = {0.1,0.1,0.1,0.1};
  
  double deltaEcat[4] = {0,0,0,0};
  
  smearcatMax[0] = 0.05;
  smearcatMax[1] = 0.05;
  smearcatMax[2] = 0.05;
  smearcatMax[3] = 0.05;
  
  smearcat[0] = 0.005;
  smearcat[1] = 0.01;
  smearcat[2] = 0.015;
  smearcat[3] = 0.02;
  
  if(barrelorEndcap==2){
    smearcat[0] = 0.02;
    smearcat[1] = 0.02;
    smearcat[2] = 0.02;
    smearcat[3] = 0.02;
  }
  
  for(int j=0; j<4; j++){
    TString parname = TString (Form("smearcat%d",j));
    minuit->mnparm(j, parname, smearcat[j], STEPMN, 0,smearcatMax[j],ierflg);
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

  int testpaircat = testcat; 
  
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
  
  arglist[0] = 1; 
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
    //minuit->mnmnos();
  }
  
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
  
  
  
  
  mytree->Fill();
  mytree->Write();
  fnew->Write();
  fnew->Close();

  
  
}
