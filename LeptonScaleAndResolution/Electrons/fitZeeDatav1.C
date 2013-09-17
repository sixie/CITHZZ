#include "rootheader.h"
#include "roofitheader.h"

using namespace RooFit ;

#include "usefullcoderoofit.cc"
#include "fitZToMuMuGammaMassUnbinnedwtv3data.C"


void fitZeeDatav1( int energyRegrVer =1 ,int eleID = 2,int fixtoMC=1, int scaled=1, int doHistory= 0, int mMin=70, int mMax = 110){
  
  
  TString workingDIR = "/afs/cern.ch/work/y/yangyong/CMSSW/CMSSW_5_2_4/src/CITHZZ/LeptonScaleAndResolution/Electrons";
  
  TString filename = workingDIR + TString(Form("/zeeWS/makeZeeWsShape_dm1_regver%d_eleID%d_escale%d.root",energyRegrVer,eleID,scaled));
  cout<<filename<<endl;
  
  TFile *f1 =new TFile(filename,"read");
  RooWorkspace* w = (RooWorkspace*) f1->Get("zeeShape") ;
  
  vector<string> mpair_var;
  //mpair_var.push_back("mpair");
  //mpair_var.push_back("mpair_ebeb");
  
  
  //bool doHistory = true; 
  
  map<string,string> map_name; 

  if(doHistory){
    mpair_var.push_back("mpair_ebebpc3");
    mpair_var.push_back("mpair_ebebpo3");
    mpair_var.push_back("mpair_ebebmc3");
    mpair_var.push_back("mpair_ebebmo3");
    mpair_var.push_back("mpair_eeeepc");
    mpair_var.push_back("mpair_eeeepo");
    mpair_var.push_back("mpair_eeeemc");
    mpair_var.push_back("mpair_eeeemo");
    

  }else{

    //mpair_var.push_back("mpair_ebeb");
    //mpair_var.push_back("mpair_eeee");
    
    
    mpair_var.push_back("mpair_ebebc3_highr9"); //both electron inside |etasc|<1
    mpair_var.push_back("mpair_ebebc3_lowr9"); //both electron inside |etasc|<1
    mpair_var.push_back("mpair_ebebo3_highr9"); //both electron inside |etasc|<1
    mpair_var.push_back("mpair_ebebo3_lowr9"); //both electron inside |etasc|<1
    mpair_var.push_back("mpair_eeeec_highr9");
    mpair_var.push_back("mpair_eeeec_lowr9");
    mpair_var.push_back("mpair_eeeeo_highr9");
    mpair_var.push_back("mpair_eeeeo_lowr9");
  }
  
  
  
  RooRealVar *rv_mass = (RooRealVar*)w->var("rv_mass");
  RooRealVar *rv_weight = new RooRealVar("rv_weight","weight",1.0,0,1E6);
  
  
  
  map<string, float> acb_MC_cut1; 
  map<string, float> ncb_MC_cut1; 
  map<string, float> dcb_MC_cut1; 
  map<string, float> sigcb_MC_cut1; 
  map<string, float> sigcbErr_MC_cut1; 
  
  map<string, float> dErrcb_MC_cut1; 

  map<string, float> acb_MC_cut2; 
  map<string, float> ncb_MC_cut2; 
  map<string, float> dcb_MC_cut2; 
  map<string, float> dErrcb_MC_cut2; 
  
  
  ///datapu1(EPS jul5)  mcpile1


  //test training options
  //format: name, acb,ncb,dcb,dErrcb,sigcb,sigcbErr
  ///mpair_ebeb_highr9 1.5158 2.26039 0.205366 0.0165984 0.855609 0.0212227
    

  string stmp;
  float acb; 
  float ncb; 
  float dcb; 
  float dErrcb; 
  float sigcb; 
  float sigcbErr; 
  
  
  TString input =  workingDIR + TString(Form("/fitres/fitZeeMCv1_regver%d_eleID%d_fitres",energyRegrVer,eleID));
  
  cout<<input<<endl; 
  
  int fitstatus; 
  
  TString checkinput = TString("ls ") + input;
  if( gSystem->Exec(checkinput)!=0){
    return; 
  }
    
  ifstream txtin(input,ios::in);
  while( txtin.good() ){
    txtin>> stmp >> acb >> ncb >> dcb >> dErrcb >> sigcb >> sigcbErr ; 
      
    acb_MC_cut1[stmp] = acb; 
    ncb_MC_cut1[stmp] = ncb; 
    dcb_MC_cut1[stmp] = dcb; 
    dErrcb_MC_cut1[stmp] = dErrcb; 
    sigcb_MC_cut1[stmp] = sigcb; 
    sigcbErr_MC_cut1[stmp] = sigcbErr; 
      
    if( txtin.eof()) break; 
  }
    
  
  cout<<"MC parameters loaded " <<endl; 
    
  
  
  
  
  float a_fit; 
  float n_fit; 
  float dm_fit; 
  float sigcb_fit; 
    
  
  float a_fitErr; 
  float n_fitErr; 
  float dm_fitErr; 
  float sigcb_fitErr; 
  
  
  float aMC; 
  float nMC; 
  float dMC; 
  float dErrMC; 
  float sigMC; 
  float sigMCErr; 
  
  float weightscale = 1.0;
  
  filename   = TString (Form("fitZeeDatav1.regver%d.eleid%d.fixtoMC%d.scale%d.doHistory%d.txt",energyRegrVer,eleID,fixtoMC,scaled,doHistory));
  ofstream txtout(filename,ios::out);
  
  int runMin = 132440;
  int runMax = 999999;
  
  
//   if(runRange==-1){ //evtweight set to be  1; 
//     runMin = 0; 
//     runMax = 10; 
//   }
//   else if(runRange==0){
//   }
  

//   else if(runRange == 1){
//     runMin = 160431; 
//     runMax = 167913; 
//   }
//   else if(runRange == 2){
//     runMin = 170000;
//     runMax = 172619; 
//   }
//   else if(runRange == 3){
//     runMin = 172620;
//     runMax = 173692; 
//   }
//   else if(runRange == 4){
//     runMin = 175860;
//     runMax = 177139; 
//   }
//   else if(runRange == 5){
//     runMin = 177140;
//     runMax = 178421; 
//   }
//   else if(runRange == 6){
//     runMin = 178424;
//     runMax = 180252; 
//   }
  
  
  bool plotcdfrms = false; 
  //bool plotcdfrms = true; 
    
  

  float mZ = 91.1876;
  

  //2011
//   int runMinAll[10] = {0,     160431,170000,172620,175860,177140,178424};
//   int runMaxAll[10] = {999999,167913,172619,173692,177139,178421,180252};
  
  //2012
  int runMinAll[100] = {0,     190645,190782,191043,193556,194151,194533,195114,195916,198116,199804,200049,200152,200491,200532,201657,202306};
  int runMaxAll[100] = {999999,190781,191042,193555,194150,194532,195113,195915,198115,199803,200048,200151,200490,200531,201656,202305,203002};
  
  
  int runRangeStart = 0; 
  int runRangeEnd = 0; 

  if( doHistory){
    runRangeStart = 1; 
    runRangeEnd = 16; 
  }
  //runRangeEnd = 1;
  
  for(int runRange = runRangeStart; runRange <= runRangeEnd; runRange++){
    
    runMin = runMinAll[runRange];
    runMax = runMaxAll[runRange];
     
     

    for(int j=0; j< int( mpair_var.size()); j++){
      string mpairname = mpair_var[j];
      TString rname = TString( Form("rds_%s",mpairname.c_str()) );
      RooDataSet *rdinput = (RooDataSet*)w->data(rname);
      ///string gifname = mpairname + string(Form("%s_etcut%d_corr%d_fixtoMC%d_runRange%d",sample,etcut,corr,fixtoMC,runRange));
      string gifname = mpairname + string( Form("fixtoMC%d_regver%d_eleID%d_runRange%d",fixtoMC,energyRegrVer,eleID,runRange));
      
      RooDataSet *rd = cwdsetv2(rdinput,rv_mass,rv_weight,rname,mMin,mMax,runMin,runMax,weightscale,1);
    
      string sname = mpairname; 

      if(sname.find("eeeep")!= string::npos 
	 || sname.find("eeeem")!= string::npos 
	 || sname.find("ebebp") != string::npos
	 || sname.find("ebebm") != string::npos
	 ){
	sname.erase(10,1);
      }
      
      
      aMC = acb_MC_cut1[sname];
      nMC = ncb_MC_cut1[sname];
      dMC = dcb_MC_cut1[sname];
      dErrMC = dErrcb_MC_cut1[sname];
      sigMC = sigcb_MC_cut1[sname];
      sigMCErr = sigcbErr_MC_cut1[sname];
    
      if(aMC ==0 || nMC ==0 || dMC ==0 ){
	cout<<"bad MC " << aMC <<" "<< nMC<<" "<< dMC <<" "<<mpairname.c_str()<<endl; 
	return; 
      }

      
      if(fixtoMC==0){
      	fitZToMuMuGammaMassUnbinnedv3data(dm_fit,dm_fitErr,a_fit,n_fit,sigcb_fit,sigcb_fitErr,rd,rv_mass,mMin,mMax,aMC,nMC,plotcdfrms,false,".",gifname.c_str());
	rname = TString( Form("rds_%s_fixtoMC%d",mpairname.c_str(),0));
	float scale =  ( mZ + dMC )/(mZ + dm_fit);
	txtout<<" scaleBeforeShift " << mpairname.c_str()<<" "<<  dm_fit<<" +/- "<< dm_fitErr<< " mc "<< dMC <<" +/- " << dErrMC <<endl; 
	float scaleShift = (dm_fit - dMC) / ( mZ + dm_fit ); 
	float scaleShiftErr = sqrt(dm_fitErr*dm_fitErr + dErrMC* dErrMC ) / ( mZ + dm_fit);  
	
	RooDataSet *rdv1 =  cwdsetv2(rdinput,rv_mass,rv_weight,rname,mMin,mMax,runMin,runMax,weightscale,scale);
	//	gifname = mpairname + string(Form("_%s_etcut%d_corr%d_fixtoMC%d_scaletoMC_runRange%d",sample,etcut,corr,fixtoMC,runRange));
	string gifname = mpairname + string( Form("fixtoMC%d_scaletoMC_regver%d_eleID%d_runRange%d",fixtoMC,energyRegrVer,eleID,runRange));
	
	fitZToMuMuGammaMassUnbinnedv3data(dm_fit,dm_fitErr,a_fit,n_fit,sigcb_fit,sigcb_fitErr,rdv1,rv_mass,mMin,mMax,aMC,nMC,plotcdfrms,false,".",gifname.c_str());

	txtout<< mpairname<<" dm "<< dm_fit<<" +/- "<< dm_fitErr <<" reso "<<sigcb_fit <<" +/- "<< sigcb_fitErr <<endl; 
	float smearing = 0; 
	float smearingStatErr = 0;
	if( sigcb_fit > sigMC){
	  smearing = sqrt(2)/ (mZ + dMC ) * sqrt( sigcb_fit * sigcb_fit - sigMC * sigMC);
	  smearingStatErr = sqrt(2)/ (mZ + dMC )  * 
	    sqrt( pow( sigcb_fit * sigcb_fitErr/sqrt( sigcb_fit * sigcb_fit - sigMC * sigMC),2) +
		  pow( sigMC * sigMCErr /sqrt( sigcb_fit * sigcb_fit - sigMC * sigMC),2) );
	  
	}
	txtout<<"scale_smearing " << mpairname.c_str()<<" "<< scaleShift <<" smearing "<< smearing*100 <<" +/- " << smearingStatErr*100<<" %" <<endl;
	
      }else{
	
	plotcdfrms = false;
	
	fitZToMuMuGammaMassUnbinnedv3data(dm_fit,dm_fitErr,a_fit,n_fit,sigcb_fit,sigcb_fitErr,rd,rv_mass,mMin,mMax,aMC,nMC,plotcdfrms,true,".",gifname.c_str());
		
	///scale data to MC
	rname = TString( Form("rds_%s_scaletoMC",mpairname.c_str()));
	float scale =  ( mZ + dMC )/(mZ + dm_fit);
      
	txtout<<" scaleBeforeShift " << mpairname.c_str()<<" "<<  dm_fit<<" +/- "<< dm_fitErr<<" "<<" mc "<< dMC <<" +/- " << dErrMC <<" sigfit "<<  sigcb_fit <<" +/- " << sigcb_fitErr <<endl; 
      	float scaleShift = (dm_fit - dMC) / ( mZ + dm_fit ); 
	float scaleShiftErr = sqrt(dm_fitErr*dm_fitErr + dErrMC* dErrMC ) / ( mZ + dm_fit);  
	
	TString textout = TString(Form("scaleForText %s & %d-%d & $%3.2f \\pm %3.2f$ \\\\ \\hline", mpairname.c_str(),runMin,runMax,scaleShift*100,scaleShiftErr*100));
	txtout<<textout<<endl; 
	
	textout = TString(Form("scaleForRead %s %dto%d %5.4f %5.4f ", mpairname.c_str(),runMin,runMax,scaleShift,scaleShiftErr));
	txtout<<textout<<endl; 
	
	
	///return; 
		
	RooDataSet *rdv1 =  cwdsetv2(rdinput,rv_mass,rv_weight,rname,mMin,mMax,runMin,runMax,weightscale,scale);
	gifname = mpairname + string( Form("fixtoMC%d_scaletoMC_regver%d_eleID%d_runRange%d",fixtoMC,energyRegrVer,eleID,runRange));
	
	fitZToMuMuGammaMassUnbinnedv3data(dm_fit,dm_fitErr,a_fit,n_fit,sigcb_fit,sigcb_fitErr,rdv1,rv_mass,mMin,mMax,aMC,nMC,plotcdfrms,true,".",gifname.c_str());
	
	///print out (d-m)/d and sqrt(2)/(m) *sqrt(Sd*Sd -Sm*Sm);
	
	txtout<< mpairname<<" dm "<< dm_fit<<" +/- "<< dm_fitErr <<" reso "<<sigcb_fit <<" +/- "<< sigcb_fitErr <<" mc "<< sigMC <<"+/-"<< sigMCErr<< endl; 
	
	
	float smearing = 0; 
	float smearingStatErr = 0;
	if( sigcb_fit > sigMC){
	  smearing = sqrt(2)/ (mZ + dMC ) * sqrt( sigcb_fit * sigcb_fit - sigMC * sigMC);
	  smearingStatErr = sqrt(2)/ (mZ + dMC )  * 
	    sqrt( pow( sigcb_fit * sigcb_fitErr/sqrt( sigcb_fit * sigcb_fit - sigMC * sigMC),2) +
		  pow( sigMC * sigMCErr /sqrt( sigcb_fit * sigcb_fit - sigMC * sigMC),2) );
	  
	}
	txtout<<"scale_smearing " << mpairname.c_str()<<" "<< scaleShift <<" smearing "<< smearing*100 <<" +/- " << smearingStatErr*100<<" %" <<endl;
	
	textout = TString( Form("gettext %s $ %3.2f\\pm%3.2f $ & $ %3.2f\\pm%3.2f$ & $%3.2f\\pm%3.2f$ \\\\ \\hline",mpairname.c_str(), sigcb_fit,sigcb_fitErr,sigMC,sigMCErr,smearing*100,100*smearingStatErr));
	txtout<<textout<<endl;
	
	
      }
      
    }

  
  }
  
  
}
