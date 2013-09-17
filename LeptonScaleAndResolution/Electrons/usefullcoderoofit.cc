

/////simply copy and scale weight , get only fraction of the data 
RooDataSet *cwdset(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, TString name, double vmin, double vmax,Double_t weightscale, int flag) {
    
  
  RooDataSet *outdata = new RooDataSet(name,"",RooArgList(*mvar,*wvar),wvar->GetName());
  for (Int_t ient=0; ient<indata->numEntries(); ++ient) {

    if(flag==1){
      if( ient%2==0) continue; 
    }
    else if(flag==2){
      if( ient%2==1) continue; 
    }
    else if(flag==3){
      if( ient%4!=0) continue; 
    }
    else if(flag==4){
      if( ient%4==0) continue; 
    }
    else if(flag==5){
      if( ient%10!=0) continue; 
    }
    else if(flag==6){
      if( ient%10==0) continue; 
    }
    else if(flag==7){
      if( ient%50!=0) continue; 
    }
    else if(flag==8){
      if( ient%50==0) continue; 
    }
    else if(flag==9){
      if( ient%100!=0) continue; 
    }
    else if(flag==10){
      if( ient%100==0) continue; 
    }
    else if(flag==11){
      if( ient%200!=0) continue; 
    }
    else if(flag==12){
      if( ient%200==0) continue; 
    }

    const RooArgSet *ent = indata->get(ient);
    double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
    if( val >vmin && val < vmax){
      mvar->setVal(val);
      outdata->add(*mvar,weightscale*indata->weight());
      //cout<<"val " << val <<" "<< indata ->weight()<<endl; 
    }
  }
  
  return outdata;
    
}



/////simply copy and scale weight , get only fraction of the data 
RooDataSet *cwdset(RooDataSet *indata, RooFormulaVar *mvar1, RooRealVar *mvar, RooRealVar *wvar, TString name, double vmin, double vmax,Double_t weightscale, int flag) {
  
  
  RooDataSet *outdata = new RooDataSet(name,"",RooArgList(*mvar1,*wvar),wvar->GetName());
  for (Int_t ient=0; ient<indata->numEntries(); ++ient) {

    if(flag==1){
      if( ient%2==0) continue; 
    }
    else if(flag==2){
      if( ient%2==1) continue; 
    }
    if(flag==3){
      if( ient%4!=0) continue; 
    }
    else if(flag==4){
      if( ient%4==0) continue; 
    }
    if(flag==5){
      if( ient%10!=0) continue; 
    }
    else if(flag==6){
      if( ient%10==0) continue; 
    }
    if(flag==7){
      if( ient%50!=0) continue; 
    }
    else if(flag==8){
      if( ient%50==0) continue; 
    }
    if(flag==9){
      if( ient%100!=0) continue; 
    }
    else if(flag==10){
      if( ient%100==0) continue; 
    }
    
    const RooArgSet *ent = indata->get(ient);
    double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
    if( val >vmin && val < vmax){
      mvar->setVal(val);
      outdata->add(*mvar,weightscale*indata->weight());
      //cout<<"val " << val <<" "<< indata ->weight()<<endl; 
    }
  }
  
  return outdata;
  
}

/////simply copy and scale weight 
RooDataSet *cwdset(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, TString name, double vmin, double vmax,Double_t weightscale) {
  
  RooDataSet *outdata = new RooDataSet(name,"",RooArgList(*mvar,*wvar),wvar->GetName());
  for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
    const RooArgSet *ent = indata->get(ient);
    double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
    if( val >vmin && val < vmax){
      mvar->setVal(val);
      outdata->add(*mvar,weightscale*indata->weight());
      //cout<<"val " << val <<" "<< indata ->weight()<<endl; 
    }
  }
  
  return outdata;
    
}


// /////simply copy and scale weight 
// RooDataSet *cwdsetv1(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, TString name, double vmin, double vmax,Double_t weightscale,  float scale,float smear) {
  
//   RooDataSet *outdata = new RooDataSet(name,"",RooArgList(*mvar,*wvar),wvar->GetName());
//   for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
//     const RooArgSet *ent = indata->get(ient);
//     double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
//     float tmp = getGaussian(1+scale,smear,0,1E3) * getGaussian(1+scale,smear,0,1E3);
//     val *= sqrt(tmp);
//     if( val >vmin && val < vmax){
//       mvar->setVal(val);
//       outdata->add(*mvar,weightscale*indata->weight());
//       //cout<<"val " << val <<" "<< indata ->weight()<<endl; 
//     }
//   }
  
//   return outdata;
  
// }

TH1F *convertRooDataSetToTH1F(RooDataSet *indata,RooRealVar *mvar, TString name,  double vmin, double vmax,Double_t weightscale){
  TH1F *outdata = new TH1F(name,name,60,vmin,vmax);
  for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
    const RooArgSet *ent = indata->get(ient);
    double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
    if( val >vmin && val < vmax){
      //mvar->setVal(val);
      //outdata->add(*mvar,weightscale*indata->weight());
      outdata->Fill(val,weightscale*indata->weight());
      //cout<<"val " << val <<" "<< indata ->weight()<<endl; 
    }
  }
  
  return outdata;
    
}


/////simply copy and scale weight 
RooDataSet *cwdsetv1(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, TString name, double vmin, double vmax,Double_t weightscale,  double scaleMass) {
  
  RooDataSet *outdata = new RooDataSet(name,"",RooArgList(*mvar,*wvar),wvar->GetName());
  for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
    const RooArgSet *ent = indata->get(ient);
    double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
    val *= scaleMass;
    if( val >vmin && val < vmax){
      mvar->setVal(val);
      outdata->add(*mvar,weightscale*indata->weight());
      //cout<<"val " << val <<" "<< indata ->weight()<<endl; 
    }
  }
  
  return outdata;
  
}




vector<float> getweightds(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, double vmin, double vmax){
  
  ///RooDataSet *outdata = new RooDataSet(name,"",RooArgList(*mvar,*wvar),wvar->GetName());
  vector<float> wts; 
  for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
    const RooArgSet *ent = indata->get(ient);
    double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
    
    if(val >vmin && val < vmax){
      wts.push_back( indata->weight() );
    }
  }
  return wts; 
    
}





/////simply copy , choose weight ( for data , choose runNumber)
RooDataSet *cwdsetv2(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, TString name, double vmin, double vmax,double wmin, double wmax,Double_t weightscale, double scaleMass) {
  
  RooDataSet *outdata = new RooDataSet(name,"",RooArgList(*mvar,*wvar),wvar->GetName());
  for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
    const RooArgSet *ent = indata->get(ient);
    double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
    val *= scaleMass;
    if(indata->weight() >= wmin && indata->weight() <= wmax &&  val >vmin && val < vmax){
      mvar->setVal(val);
      outdata->add(*mvar,weightscale*1);
      //cout<<"val " << val <<" "<< indata ->weight()<<endl; 
    }
  }
  
  return outdata;
  
}


/// addpend and scalw weight 
void appendcwd(RooDataSet *outdata, RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar,  double vmin, double vmax,Double_t weightscale) {
  
  
  
    for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
      const RooArgSet *ent = indata->get(ient);
      double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
      if(val >vmin && val < vmax){
	mvar->setVal(val);

	//cout<<"check " << val <<" " << indata->weight()<<endl;  
	
	outdata->add(*mvar,weightscale*indata->weight());
      }
      //      if( ient < 10) cout<<"checke " << mvar->getVal()<<" "<< indata->weight() <<endl; 
      
    }
    ///return outdata;
}


// /// addpend and scalw weight 
// void appendcwdv1(RooDataSet *outdata, RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar,  double vmin, double vmax,Double_t weightscale, float scale,float smear) {
  
//     for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
//       const RooArgSet *ent = indata->get(ient);
//       double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
//       float tmp = getGaussian(1+scale,smear,0,1E3) * getGaussian(1+scale,smear,0,1E3);
//       val *= sqrt(tmp);
      
//       if(val >vmin && val < vmax){
// 	mvar->setVal(val);
	
// 	//cout<<"check " << val <<" " << indata->weight()<<endl;  
	
// 	outdata->add(*mvar,weightscale*indata->weight());
//       }
//       //      if( ient < 10) cout<<"checke " << mvar->getVal()<<" "<< indata->weight() <<endl; 
      
//     }
//     ///return outdata;
// }


/// addpend and scalw weight 
void printRooDataSetToFile(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, const char *filename) {
  
  ofstream testme(filename,ios::out);
  
  for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
    const RooArgSet *ent = indata->get(ient);
    double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
    testme<<  val<<" "<< indata->weight() <<endl; 
    
  }
    ///return outdata;
}



