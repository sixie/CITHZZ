#include "rootheader.h"
#include "roofitheader.h"


TChain *fChain; 

#include "zeevariables.h"
#include "setbranchesv1.cc"

#include "roodatasetth1.cc"

#include "utils.cc"
#include "energyScaleCorrection.cc"



void makeZeeWsShape(int dataOrMC =1, int energyRegrVer =1 ,int eleID = 2, int applyEScale =0){
  
  TString inputdatafile = "/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_HCP2012.root";
  TString inputmcfile = "/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Summer12DY53X.root";
  
  
  fChain = new TChain("ZeeEvent");
  setbranchaddress();
  
  TString filename = TString(Form("makeZeeWsShape_dm%d_regver%d_eleID%d_escale%d.root",dataOrMC,energyRegrVer,eleID,applyEScale));
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
  
  
  vector<string> mpair_var;
  mpair_var.push_back("mpair_ebeb");
  mpair_var.push_back("mpair_ebeb_highr9");
  mpair_var.push_back("mpair_ebeb_lowr9");
  mpair_var.push_back("mpair_ebebc3_highr9");
  mpair_var.push_back("mpair_ebebc3_lowr9");
  mpair_var.push_back("mpair_ebebo3_highr9");
  mpair_var.push_back("mpair_ebebo3_lowr9");
  mpair_var.push_back("mpair_ebebc3");
  mpair_var.push_back("mpair_ebebo3");
  
  mpair_var.push_back("mpair_eeee");
  mpair_var.push_back("mpair_eeee_highr9");
  mpair_var.push_back("mpair_eeee_lowr9");
  mpair_var.push_back("mpair_eeeec_highr9");
  mpair_var.push_back("mpair_eeeec_lowr9");
  mpair_var.push_back("mpair_eeeeo_highr9");
  mpair_var.push_back("mpair_eeeeo_lowr9");
  
  mpair_var.push_back("mpair_eeeec");
  mpair_var.push_back("mpair_eeeeo");

  //for histogry
  mpair_var.push_back("mpair_ebebpc3");
  mpair_var.push_back("mpair_ebebpo3");
  mpair_var.push_back("mpair_ebebmc3");
  mpair_var.push_back("mpair_ebebmo3");
  mpair_var.push_back("mpair_eeeepc");
  mpair_var.push_back("mpair_eeeepo");
  mpair_var.push_back("mpair_eeeemc");
  mpair_var.push_back("mpair_eeeemo");
  
  rv_mass = new RooRealVar("rv_mass","mass",100,0,1000);
  rv_weight = new RooRealVar("rv_weight","weight",1.0,0,1E6);
  for(int j=0; j< int( mpair_var.size()); j++){
    string mpairname = mpair_var[j];
    TString rname = TString( Form("rhs_%s",mpairname.c_str()) );
    makeRootDataSetAndTH1F(mpairname,120,60,120);
  }

  map<string,int> map_mpaircategory; 
  
  cout<<"totalEntries  " << totalEntries <<endl; 
  
  //totalEntries = 1E6; 
  
  if(applyEScale>=1){
    load_energyScaleCorrection_eta_RunByRun();

    if(applyEScale==2){
      load_energyScaleCorrection_etaR9bin();
    }
  }
  

  for(int entry=0; entry< totalEntries; entry++){
    fChain->GetEntry(entry);
    if(entry%100000==0) cout<<"entry " << entry <<endl; 
    
    
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
      if(et1<20 || et2 <20) continue; 
      if( Ele1PassMediumSimpleCuts ==0 || Ele2PassMediumSimpleCuts == 0) continue; 
    }
    
    
    float mee = calcZmass(e1,Ele1Eta,Ele1Phi,e2,Ele2Eta,Ele2Phi);
    
    
    bool bothEB =  fabs(Ele1SCEta) < 1.49 && fabs(Ele2SCEta) < 1.49;
    bool bothEBCentral =  fabs(Ele1SCEta) < 1. && fabs(Ele2SCEta) < 1.;
    bool bothEBOuter = bothEB && fabs(Ele1SCEta) > 1. && fabs(Ele2SCEta) > 1.;
    bool bothEBMinus = bothEB && Ele1SCEta <0 && Ele2SCEta <0; 
    bool bothEBPlus = bothEB && Ele1SCEta >0 && Ele2SCEta >0; 
    
    bool bothEE =  fabs(Ele1SCEta) > 1.49 && fabs(Ele2SCEta) > 1.49;
    bool bothEEHighEta = bothEE &&  fabs(Ele1SCEta) > 2 && fabs(Ele2SCEta) > 2;
    bool bothEELowEta = bothEE &&  fabs(Ele1SCEta) < 2 && fabs(Ele2SCEta) < 2;
    
    bool bothEEMinus = bothEE && Ele1SCEta <0 && Ele2SCEta < 0; 
    bool bothEEPlus = bothEE && Ele1SCEta >0 && Ele2SCEta > 0;
    
    bool bothHighR9 = Ele1R9 > 0.94 && Ele2R9 > 0.94; 
    bool bothLowR9 = Ele1R9 < 0.94 && Ele2R9 < 0.94; 
    
    
    for(int j=0; j< int( mpair_var.size()); j++){
      string mpairname = mpair_var[j];
      map_mpaircategory[mpairname] = -1; 
    }
    
    map_mpaircategory["mpair_ebeb"] = bothEB ;
    map_mpaircategory["mpair_ebeb_highr9"] = bothEB && bothHighR9;
    map_mpaircategory["mpair_ebeb_lowr9"] = bothEB && bothLowR9;
    map_mpaircategory["mpair_ebebc3_highr9"] = bothEBCentral && bothHighR9; 
    map_mpaircategory["mpair_ebebc3_lowr9"] = bothEBCentral && bothLowR9; 
    map_mpaircategory["mpair_ebebo3_highr9"] = bothEBOuter && bothHighR9; 
    map_mpaircategory["mpair_ebebo3_lowr9"] = bothEBOuter && bothLowR9; 
    
    map_mpaircategory["mpair_eeee"] = bothEE ;
    map_mpaircategory["mpair_eeee_highr9"] = bothEE && bothHighR9;
    map_mpaircategory["mpair_eeee_lowr9"] = bothEE && bothLowR9;
    map_mpaircategory["mpair_eeeec_highr9"] = bothEELowEta && bothHighR9; 
    map_mpaircategory["mpair_eeeec_lowr9"] = bothEELowEta && bothLowR9; 
    map_mpaircategory["mpair_eeeeo_highr9"] = bothEEHighEta && bothHighR9; 
    map_mpaircategory["mpair_eeeeo_lowr9"] = bothEEHighEta && bothLowR9; 
        
    //history
    map_mpaircategory["mpair_ebebpc3"] = bothEBCentral && bothEBPlus; 
    map_mpaircategory["mpair_ebebmc3"] = bothEBCentral && bothEBMinus; 
    map_mpaircategory["mpair_ebebpo3"] = bothEBOuter && bothEBPlus; 
    map_mpaircategory["mpair_ebebmo3"] = bothEBOuter && bothEBMinus; 
    map_mpaircategory["mpair_eeeepc"] = bothEELowEta && bothEEPlus; 
    map_mpaircategory["mpair_eeeemc"] = bothEELowEta && bothEEMinus; 
    map_mpaircategory["mpair_eeeepo"] = bothEEHighEta && bothEEPlus; 
    map_mpaircategory["mpair_eeeemo"] = bothEEHighEta && bothEEMinus; 
    
    map_mpaircategory["mpair_ebebc3"] = bothEBCentral; 
    map_mpaircategory["mpair_ebebo3"] = bothEBOuter; 
    map_mpaircategory["mpair_eeeec"] = bothEELowEta;
    map_mpaircategory["mpair_eeeeo"] = bothEEHighEta;
    
    
    double evtweight = weight; 
    if(dataOrMC==1){
      evtweight = run;
    }
    
    for(int j=0; j< int( mpair_var.size()); j++){
      string mpairname = mpair_var[j];
      if( map_mpaircategory[mpairname] ==-1){
	cout<<"wrong category ? " << mpairname.c_str()<<endl; 
	return; 
      }
      if( map_mpaircategory[mpairname]==1){
	if(dataOrMC==1) fillRootDataSetAndTH1F(mpairname,mee,evtweight,true);
	else  fillRootDataSetAndTH1F(mpairname,mee,evtweight);
      }
    }
  }
  
  

  fnew->cd();
  //RooContainer_Save();
  RooWorkspace *w = new RooWorkspace("zeeShape","workspace") ;
  for (std::map<string,RooDataSet*>::iterator it_data = rds_map.begin()
         ;it_data != rds_map.end();it_data++)   {
    w->import(*it_data->second);
  }
  
  
  w->Write();

  fnew->Write();
  fnew->Close();



}
