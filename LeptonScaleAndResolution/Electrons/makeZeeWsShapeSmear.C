#include "rootheader.h"
#include "roofitheader.h"


TChain *fChain; 

#include "zeevariables.h"
#include "setbranchesv1.cc"

#include "roodatasetth1.cc"

#include "utils.cc"
#include "energyScaleCorrection.cc"

int Ntrial = 100;


///for validation of the smearing 
void makeZeeWsShapeSmear(int dataOrMC =1, int energyRegrVer =1 ,int eleID = 2, int applyEScale =0, int smear=0, int sys=0){
  
  TString inputdatafile = "/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_HCP2012.root";
  TString inputmcfile = "/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Summer12DY53X.root"; ///event%2==0 for training. 

  if(dataOrMC==1 && smear>=1){
    cout<<"data !"<<endl; 
    return; 
  }
  
  fChain = new TChain("ZeeEvent");
  setbranchaddress();
  
  TString filename = TString(Form("makeZeeWsShapeSmear_dm%d_regver%d_eleID%d_escale%d_smear%d_sys%d.root",dataOrMC,energyRegrVer,eleID,applyEScale,smear,sys));
  TFile *fnew = new TFile(filename,"recreate");
  
  
  if(dataOrMC ==1){
    fChain->Add(inputdatafile);
  }else if(dataOrMC==2){
    fChain->Add(inputmcfile);
  }else{
    cout<<"dataOrMC 1/2 for data/mc " <<  endl; 
    return; 
  }
  
  cout<<"data chain " << fChain->GetNtrees()<<endl; 
  
  int totalEntries = fChain->GetEntries();
  
  
//   vector<string> mpair_var;
//   mpair_var.push_back("mpair_ebeb");
//   mpair_var.push_back("mpair_ebeb_highr9");
//   mpair_var.push_back("mpair_ebeb_lowr9");
//   mpair_var.push_back("mpair_ebebc3_highr9");
//   mpair_var.push_back("mpair_ebebc3_lowr9");
//   mpair_var.push_back("mpair_ebebo3_highr9");
//   mpair_var.push_back("mpair_ebebo3_lowr9");
//   mpair_var.push_back("mpair_ebebc3");
//   mpair_var.push_back("mpair_ebebo3");
  
//   mpair_var.push_back("mpair_eeee");
//   mpair_var.push_back("mpair_eeee_highr9");
//   mpair_var.push_back("mpair_eeee_lowr9");
//   mpair_var.push_back("mpair_eeeec_highr9");
//   mpair_var.push_back("mpair_eeeec_lowr9");
//   mpair_var.push_back("mpair_eeeeo_highr9");
//   mpair_var.push_back("mpair_eeeeo_lowr9");
  
//   mpair_var.push_back("mpair_eeeec");
//   mpair_var.push_back("mpair_eeeeo");

//   //for histogry
//   mpair_var.push_back("mpair_ebebpc3");
//   mpair_var.push_back("mpair_ebebpo3");
//   mpair_var.push_back("mpair_ebebmc3");
//   mpair_var.push_back("mpair_ebebmo3");
//   mpair_var.push_back("mpair_eeeepc");
//   mpair_var.push_back("mpair_eeeepo");
//   mpair_var.push_back("mpair_eeeemc");
//   mpair_var.push_back("mpair_eeeemo");
  
//   rv_mass = new RooRealVar("rv_mass","mass",100,0,1000);
//   rv_weight = new RooRealVar("rv_weight","weight",1.0,0,1E6);
//   for(int j=0; j< int( mpair_var.size()); j++){
//     string mpairname = mpair_var[j];
//     TString rname = TString( Form("rhs_%s",mpairname.c_str()) );
//     makeRootDataSetAndTH1F(mpairname,120,60,120);
//   }
  
  
  ///10 smearing categories
  for(int j=0;j<10;j++){
    string histname = string(Form("ebeb_smearmcat%d",j));
    makeTH1F(histname,240,60,120);
    histname = string(Form("eeee_smearmcat%d",j));
    makeTH1F(histname,240,60,120);
  }
  for(int j=0;j<4;j++){
    string histname = string(Form("mpair_cat%d",j));
    makeTH1F(histname,240,60,120);
  }
  
  
//   map<string,int> map_mpaircategory; 
  
  cout<<"totalEntries  " << totalEntries <<endl; 
  
  //totalEntries = 1E6; 
  
  if(applyEScale>=1){
    load_energyScaleCorrection_eta_RunByRun();

    if(applyEScale==2){
      load_energyScaleCorrection_etaR9bin();
    }
  }

  float smearcat[8]  = {
    0.0095,
    0.0116,
    0.0180,
    0.0207,
    0.0289,
    0.0296,
    0.0326,
    0.0340
  };

  float smearsys = 0.001;
  
  for(int j=0; j<10; j++){
    smearcat[j] += sys * smearsys; 
    if( smearcat[j]<0) smearcat[j] = 0; 
  }
  

  TRandom3 *grnd = new TRandom3(12345);
    
  
  for(int entry=0; entry< totalEntries; entry++){
    fChain->GetEntry(entry);
    if(entry%100000==0) cout<<"entry " << entry <<endl; 
    
    if(dataOrMC==2 && event%2==0) continue; 
    
    
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
    
    if(applyEScale>=1 && dataOrMC ==1 ){


      float corr1 = energyScaleCorrection_eta_RunByRun(Ele1SCEta,run);
      e1 *= ( 1- corr1);
      float corr2 = energyScaleCorrection_eta_RunByRun(Ele2SCEta,run);
      e2 *= ( 1- corr2);

      if(applyEScale==2){
	corr1 = energyScaleCorrection_etaR9bin(Ele1SCEta,Ele1R9);
	e1 *= (1-corr1);
	corr2 = energyScaleCorrection_etaR9bin(Ele2SCEta,Ele2R9);
	e2 *= (1-corr2);
      }
      
    }
    
    
    float et1 = e1 * sin(2*atan(exp(-Ele1Eta)));
    float et2 = e2 * sin(2*atan(exp(-Ele2Eta)));

    if(eleID==2) {
      if( Ele1PassMediumSimpleCuts ==0 || Ele2PassMediumSimpleCuts == 0) continue; 
    }
    
    
    float mee = calcZmass(e1,Ele1Eta,Ele1Phi,e2,Ele2Eta,Ele2Phi);
    bool bothEB =  fabs(Ele1SCEta) < 1.49 && fabs(Ele2SCEta) < 1.49;
    bool bothEE =  fabs(Ele1SCEta) > 1.49 && fabs(Ele2SCEta) > 1.49;

    int sccat1 = scCategoryEight(Ele1SCEta,Ele1R9);
    int sccat2 = scCategoryEight(Ele2SCEta,Ele2R9);

    int sc1cat = sccat1;
    int sc2cat = sccat2;

    int mcat = -1; 
    if(bothEB){
      mcat = MeesmearCategory(sccat1,sccat2);
    }
    if(bothEE){
      sccat1 -=4;
      sccat2 -=4;
      mcat = MeesmearCategory(sccat1,sccat2);
    }
    
    int mcat4 = pairscCategoryFour(Ele1SCEta,Ele1R9,Ele2SCEta,Ele2R9);
    
    if(smear==0){
      if(eleID==2) {
	if(et1<20 || et2 <20) continue; 
      }	
      if(mcat>=0 && bothEB){
	string histname = string(Form("ebeb_smearmcat%d",mcat));
	fillTH1F(histname,mee,weight);
      }
      if(mcat>=0 && bothEE){
	string histname = string(Form("eeee_smearmcat%d",mcat));
	fillTH1F(histname,mee,weight);
      }
      string histname = string(Form("mpair_cat%d",mcat4));
      fillTH1F(histname,mee,weight);
      
    }else{
      
      for(int n=0; n< Ntrial; n++){
	float e1s = grnd->Gaus(e1,e1*smearcat[sc1cat]);
	float e2s = grnd->Gaus(e2,e2*smearcat[sc2cat]);
	float et1s = e1s * sin(2*atan(exp(-Ele1Eta)));
	float et2s = e2s * sin(2*atan(exp(-Ele2Eta)));
	if(eleID==2) {
	  if(et1s<20 || et2s <20) continue; 
	}
	float mees = sqrt( (e1s*e2s)/(e1*e2)) * mee;
	if(mcat>=0 && bothEB){
	  string histname = string(Form("ebeb_smearmcat%d",mcat));
	  fillTH1F(histname,mees,weight);
	}
	if(mcat>=0 && bothEE){
	  string histname = string(Form("eeee_smearmcat%d",mcat));
	  fillTH1F(histname,mees,weight);
	}
	string histname = string(Form("mpair_cat%d",mcat4));
	fillTH1F(histname,mees,weight);
      }
      
    }
    
  }
  
  

//   fnew->cd();
//   //RooContainer_Save();
//   RooWorkspace *w = new RooWorkspace("zeeShape","workspace") ;
//   for (std::map<string,RooDataSet*>::iterator it_data = rds_map.begin()
//          ;it_data != rds_map.end();it_data++)   {
//     w->import(*it_data->second);
//   }
  
  
//   w->Write();

  fnew->Write();
  fnew->Close();



}
