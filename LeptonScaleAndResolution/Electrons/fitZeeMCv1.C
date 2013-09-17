#include "rootheader.h"
#include "roofitheader.h"

using namespace RooFit ;

#include "usefullcoderoofit.cc"
#include "fitZToMuMuGammaMassUnbinnedwtv3.C"


void fitZeeMCv1( int energyRegrVer =1 ,int eleID = 2, int mMin=70, int mMax = 110){
  
  
  TString workingDIR = "/afs/cern.ch/work/y/yangyong/CMSSW/CMSSW_5_2_4/src/CITHZZ/LeptonScaleAndResolution/Electrons";
  
  TString filename = workingDIR + TString(Form("/zeeWS/makeZeeWsShape_dm2_regver%d_eleID%d_escale0.root",energyRegrVer,eleID));
  cout<<filename<<endl;
  
  TFile *f1 =new TFile(filename,"read");
  RooWorkspace* w = (RooWorkspace*) f1->Get("zeeShape") ;
  
  vector<string> mpair_var;
  

//   mpair_var.push_back("mpair_ebeb");
//   mpair_var.push_back("mpair_eeee");
    

  mpair_var.push_back("mpair_ebebc3_highr9");
  mpair_var.push_back("mpair_ebebc3_lowr9");
  mpair_var.push_back("mpair_ebebo3_highr9");
  mpair_var.push_back("mpair_ebebo3_lowr9");
  
  //mpair_var.push_back("mpair_eeee_highr9");
  //mpair_var.push_back("mpair_eeee_lowr9");
  
  
  mpair_var.push_back("mpair_eeeec_highr9");
  mpair_var.push_back("mpair_eeeec_lowr9");
  mpair_var.push_back("mpair_eeeeo_highr9");
  mpair_var.push_back("mpair_eeeeo_lowr9");
  
  
//     mpair_var.push_back("mpair_ebebc3_highr9"); //both electron inside |etasc|<0.9
//    mpair_var.push_back("mpair_ebebc3_lowr9"); //both electron inside |etasc|<0.9
//   mpair_var.push_back("mpair_ebebo3_highr9"); //both electron inside |etasc|<0.9
//   mpair_var.push_back("mpair_ebebo3_lowr9"); //both electron inside |etasc|<0.9
  
//   mpair_var.push_back("mpair_ebebc3_highr9_eleebminus_posebplus"); ///central 0.9 all r9 
//    mpair_var.push_back("mpair_ebebc3_highr9_posebminus_eleebplus"); ///central 0.9 all r9 
//    mpair_var.push_back("mpair_ebebc3_highr9_ebminus");
//    mpair_var.push_back("mpair_ebebc3_highr9_ebplus"); 
   

  
   
   
  ///for history correction
  mpair_var.push_back("mpair_ebebc3");
  mpair_var.push_back("mpair_ebebo3");

  //mpair_var.push_back("mpair_ebebmc3");
  //mpair_var.push_back("mpair_ebebpc3");
  
  //mpair_var.push_back("mpair_ebebmo3");
  //mpair_var.push_back("mpair_ebebpo3");
  //mpair_var.push_back("mpair_eeeemc");
  //mpair_var.push_back("mpair_eeeemo");
  //mpair_var.push_back("mpair_eeeepc");
  //mpair_var.push_back("mpair_eeeepo");
  mpair_var.push_back("mpair_eeeec");
  mpair_var.push_back("mpair_eeeeo");

  
  
  RooRealVar *rv_mass = (RooRealVar*)w->var("rv_mass");
  RooRealVar *rv_weight = new RooRealVar("rv_weight","weight",1.0,0,1E6);

  float a_fit; 
  float n_fit; 
  float dm_fit; 
  float dm_fitErr;
  float sigcb_fit; 
  float sigcb_fitErr;

  filename = TString(Form("fitZeeMCv1_regver%d_eleID%d_fitres",energyRegrVer,eleID));
  ofstream txtout(filename,ios::out);
  
  ///bool IsData = 0; 

  bool docdfrms = false; 
  //bool docdfrms = true;
  

  int fitstatus;
  

  map<string,string> map_mpair_textName; 

  //map_mpair_textName["mpair_ebebc3_highr9"] = "|#eta|<1.0,r_{9} >0.94 ";
  //map_mpair_textName["mpair_ebebc3_lowr9"] = "|#eta|<1.0,r_{9} <0.94 ";
  //map_mpair_textName["mpair_ebebo3_highr9"] = "1.0<|#eta|<1.5,r_{9} >0.94 ";
  //map_mpair_textName["mpair_ebebo3_lowr9"] = "1.0<|#eta|<1.5,r_{9} <0.94 ";

  //map_mpair_textName["mpair_ebeb"] = "ECAL Barrel";
  //map_mpair_textName["mpair_eeee"] = "ECAL Endcap";
  

//   string corrName = "Factorized correction";
//   if( corr >0){
//     corrName = "Regression correction";
//   }
  
//     
  string corrName = string(Form("Regression correction V%d",energyRegrVer));
  
  for(int j=0; j< int( mpair_var.size()); j++){
    string mpairname = mpair_var[j];
    rv_mass->setRange(0,1000);
    
    TString rname = TString( Form("rds_%s",mpairname.c_str()) );
    RooDataSet *rd0 = (RooDataSet*)w->data(rname);
    if(rd0==NULL){
      cout<<" no dataset " << rname <<endl; 
      return; 
    }
    TString name = rname + TString("massrange");
    RooDataSet *rd = cwdset(rd0, rv_mass,rv_weight, name, mMin,mMax,1);

    string gifname = mpairname + string( Form("_regver%d_eleID%d",energyRegrVer,eleID));
    
    fitZToMuMuGammaMassUnbinnedv3(dm_fit,dm_fitErr,a_fit,n_fit,sigcb_fit,sigcb_fitErr,rd,rv_mass,mMin,mMax,1,5,docdfrms,false,".",gifname.c_str());
    
    txtout<<mpairname.c_str()<<" "<<a_fit<<" "<<n_fit<<" "<<dm_fit<<" "<< dm_fitErr<<" "<<sigcb_fit<<" "<<sigcb_fitErr<<" "<<endl; 
    
    //return; 
    
    ///if(j==2) return; 
    
    
  }
  
  
}
