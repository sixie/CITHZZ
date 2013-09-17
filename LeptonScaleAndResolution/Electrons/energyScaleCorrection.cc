


int etabinForRunByRunCorrection(float etasc){
  if( fabs(etasc)<1.49){
    if( etasc>0) return etasc<1 ? 0:1;
    else return fabs(etasc)<1 ? 2:3;
  }else{
    if( etasc>0) return etasc<2 ? 4:5;
    else return fabs(etasc)<2 ? 6:7;
  }
  
}


float corr_etaR9bin[10];

///in different rapidity as a function of runs
void load_energyScaleCorrection_etaR9bin(){
  
  
  
  TString workingDIR = "/afs/cern.ch/work/y/yangyong/CMSSW/CMSSW_5_2_4/src/CITHZZ/LeptonScaleAndResolution/Electrons";
  
  TString inputfile = workingDIR + "/fitresdata/fitZeeDatav1.regver1.eleid2.fixtoMC1.scale1.doHistory0.txt.forread";
  ifstream txtin(inputfile,ios::in);
  
  
  string stmp;
  float corr;
  float corrErr; 
  int k = 0; 
  while(txtin.good()){
    txtin>>stmp>>stmp>>stmp>>corr>>corrErr; 
    if(txtin.eof()){
      break; 
    }
    corr_etaR9bin[k] = corr; 
    k++;
  }
  for(int j=0;j<8;j++){
    cout<<"corrbin "<< corr_etaR9bin[j]<<endl; 
  }
  
}

float energyScaleCorrection_etaR9bin(float etasc,float r9sc){
  int cat = scCategoryEight(etasc, r9sc);
  return corr_etaR9bin[cat];
}


map<int, vector<float> > map_correction_run; 
vector< pair<int,int> > run_period; 


///in different rapidity as a function of runs
void load_energyScaleCorrection_eta_RunByRun(){
  
  
  TString workingDIR = "/afs/cern.ch/work/y/yangyong/CMSSW/CMSSW_5_2_4/src/CITHZZ/LeptonScaleAndResolution/Electrons";
  
  TString inputfile = workingDIR + "/fitresdata/fitZeeDatav1.regver1.eleid2.fixtoMC1.scale0.doHistory1.txt.forread";
  ifstream txtin(inputfile,ios::in);
  
  string stmp;
  float corr;
  float corrErr; 
  int k = 0; 
  while(txtin.good()){
    txtin>>stmp>>stmp>>stmp>>corr>>corrErr; 
    
    if(txtin.eof()){
      break; 
    }
    
    string run1s = stmp.substr(0,6);
    string run2s = stmp.substr(8,6);

    int run1 = atoi(run1s.c_str());
    int run2 = atoi(run2s.c_str());
    
    int kk = k%8;
    map_correction_run[kk].push_back(corr);
    
    if(kk==0){
      run_period.push_back( make_pair(run1,run2));
    }
    
    k++;
    
  }
  
  cout<<run_period.size()<<" run periods corrections load " <<endl; 
  vector< pair<int,int> >::const_iterator it =  run_period.begin();
  for(; it !=  run_period.end(); it++){
    int run1 = it->first; 
    int run2 = it->second;
    cout<<"runtorun " << run1 <<" "<< run2 <<endl; 
  }
  
}


float energyScaleCorrection_eta_RunByRun(float etasc,int runnumb){
  
  int kk = etabinForRunByRunCorrection(etasc);

  vector< pair<int,int> >::const_iterator it =  run_period.begin();
  
  int nn = -1; 
  for(; it !=  run_period.end(); it++){
    int run1 = it->first; 
    int run2 = it->second;
    if( run1<= runnumb  && runnumb <= run2){
      nn = int( it - run_period.begin());
      break; 
    }
  }
  if(nn <0){
    cout<<"energyScaleCorrection_eta_RunByRu NA!! " << etasc <<" "<< runnumb <<endl; 
    exit(1);
  }

  if( nn < int(map_correction_run[kk].size())){
    return map_correction_run[kk][nn];
  }else{
    cout<<"wrong! map_corr " << map_correction_run[kk].size() <<" period " << nn <<endl; 
    exit(1);
  }
  
  return 0;
  
}


