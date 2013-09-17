
int scCategoryEight(float sc_eta, float sc_r9){
  
  if( fabs(sc_eta)<1.0){
    return sc_r9 >0.94  ? 0 : 1; 
  }else if( fabs(sc_eta)<1.49){
    return sc_r9 >0.94 ? 2: 3;
  }else if( fabs(sc_eta)<2.0){
    return sc_r9 > 0.94 ? 4:5; 
  }
  else{
    return sc_r9 >0.94  ? 6 : 7; 
  }
  
}


int scCategoryFour(float sc_eta, float sc_r9){
  
  if( fabs(sc_eta)<1.49){
    return sc_r9 >0.94  ? 0 : 1; 
  }else{
    return sc_r9 >0.94  ? 2 : 3; 
  }
  
}


int pairscCategoryFour(float sc1_eta, float sc1_r9, float sc2_eta, float sc2_r9){
  
  int sc1cat = scCategoryFour( sc1_eta, sc1_r9);
  int sc2cat = scCategoryFour( sc2_eta, sc2_r9);
  
  if( sc1cat ==0 && sc2cat==0) return 0;  ///barerl highr9
  else if( (sc1cat==0 &&sc2cat==1) || (sc1cat==1 && sc2cat==0) || (sc1cat==1 && sc2cat ==1)) return 1;  //barrel low/highr9 or both lowr9 
  else if( (sc1cat==0 &&sc2cat==2) || (sc1cat==2 && sc2cat==0) || (sc1cat==2 && sc2cat ==2)) return 2; ///highr9 barrel/endcap or both endcap 
  else return 3; 
  
}


///4 cats --> 10 Mee categories
int MeesmearCategory(int sccat1, int sccat2){
  
  if(sccat1>sccat2){
    swap(sccat1,sccat2);    
  }
  
  if( sccat1== sccat2){
    return sccat1; ///0,1,2,3
  }else if( sccat1==0 ){
    return sccat2+3; //4,5,6
  }else if( sccat1==1){
    return sccat2+5; //7,8
  }else if( sccat1==2 && sccat2==3){
    return 9; 
  }else{
    cout<<"wrong MeesmearCategory!" << sccat1 <<" "<< sccat2<<endl; 
    exit(1);
  }
  
  
}

float calcZmass(float e1,float eta1,float phi1,float e2,float eta2,float phi2){
  float me = 0.511*0.001;
  float p1 = e1>me? sqrt(e1*e1-me*me): e1;
  TLorentzVector v1(p1 /cosh(eta1) * cos(phi1), p1 /cosh(eta1) * sin(phi1), p1 * tanh(eta1), e1);
  float p2 = e2>me? sqrt(e2*e2-me*me): e2;
  TLorentzVector v2(p2 /cosh(eta2) * cos(phi2), p2 /cosh(eta2) * sin(phi2), p2 * tanh(eta2), e2);
  
  return (v1+v2).M();
    
}
