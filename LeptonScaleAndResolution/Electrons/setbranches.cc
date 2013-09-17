void setbranchaddress(){

 fChainData->SetBranchAddress("weight", &weight);
  fChainData->SetBranchAddress("run", &run);
  fChainData->SetBranchAddress("event", &event);
  fChainData->SetBranchAddress("rho", &rho);
  fChainData->SetBranchAddress("mass", &mass);
  fChainData->SetBranchAddress("Ele1Eta", &Ele1Eta);
  fChainData->SetBranchAddress("Ele1Pt", &Ele1Pt);
  fChainData->SetBranchAddress("Ele1Phi", &Ele1Phi);
  fChainData->SetBranchAddress("Ele1SCEta", &Ele1SCEta);
  fChainData->SetBranchAddress("Ele1SCPhi", &Ele1SCPhi);
  fChainData->SetBranchAddress("Ele1Energy", &Ele1Energy);
  fChainData->SetBranchAddress("Ele1EnergyRegression", &Ele1EnergyRegression);
  fChainData->SetBranchAddress("Ele1EnergyRegressionV0", &Ele1EnergyRegressionV0);
  fChainData->SetBranchAddress("Ele1EnergyRegressionV1", &Ele1EnergyRegressionV1);
  fChainData->SetBranchAddress("Ele1EnergyRegressionV2", &Ele1EnergyRegressionV2);
  fChainData->SetBranchAddress("Ele1EnergyRegressionErrorV0", &Ele1EnergyRegressionErrorV0);
  fChainData->SetBranchAddress("Ele1EnergyRegressionErrorV1", &Ele1EnergyRegressionErrorV1);
  fChainData->SetBranchAddress("Ele1EnergyRegressionErrorV2", &Ele1EnergyRegressionErrorV2);
  fChainData->SetBranchAddress("Ele1HZZICHEP2012IDMVA", &Ele1HZZICHEP2012IDMVA);
  fChainData->SetBranchAddress("Ele1PFIso04", &Ele1PFIso04);
  fChainData->SetBranchAddress("Ele1R9", &Ele1R9);
  fChainData->SetBranchAddress("Ele1PassLooseSimpleCuts", &Ele1PassLooseSimpleCuts);
  fChainData->SetBranchAddress("Ele1PassMediumSimpleCuts", &Ele1PassMediumSimpleCuts);
  fChainData->SetBranchAddress("Ele1PassTightSimpleCuts", &Ele1PassTightSimpleCuts);
  fChainData->SetBranchAddress("Ele1PassHZZICHEP2012", &Ele1PassHZZICHEP2012);
  fChainData->SetBranchAddress("Ele2Pt", &Ele2Pt);
  fChainData->SetBranchAddress("Ele2Eta", &Ele2Eta);
  fChainData->SetBranchAddress("Ele2Phi", &Ele2Phi);
  fChainData->SetBranchAddress("Ele2SCEt", &Ele2SCEt);
  fChainData->SetBranchAddress("Ele2SCEta", &Ele2SCEta);
  fChainData->SetBranchAddress("Ele2SCPhi", &Ele2SCPhi);
  fChainData->SetBranchAddress("Ele2Energy", &Ele2Energy);
  fChainData->SetBranchAddress("Ele2EnergyRegression", &Ele2EnergyRegression);
  fChainData->SetBranchAddress("Ele2EnergyRegressionV0", &Ele2EnergyRegressionV0);
  fChainData->SetBranchAddress("Ele2EnergyRegressionV1", &Ele2EnergyRegressionV1);
  fChainData->SetBranchAddress("Ele2EnergyRegressionV2", &Ele2EnergyRegressionV2);
  fChainData->SetBranchAddress("Ele2EnergyRegressionErrorV0", &Ele2EnergyRegressionErrorV0);
  fChainData->SetBranchAddress("Ele2EnergyRegressionErrorV1", &Ele2EnergyRegressionErrorV1);
  fChainData->SetBranchAddress("Ele2EnergyRegressionErrorV2", &Ele2EnergyRegressionErrorV2);
  fChainData->SetBranchAddress("Ele2HZZICHEP2012IDMVA", &Ele2HZZICHEP2012IDMVA);
  fChainData->SetBranchAddress("Ele2PFIso04", &Ele2PFIso04);
  fChainData->SetBranchAddress("Ele2R9", &Ele2R9);
  fChainData->SetBranchAddress("Ele2PassLooseSimpleCuts", &Ele2PassLooseSimpleCuts);
  fChainData->SetBranchAddress("Ele2PassMediumSimpleCuts", &Ele2PassMediumSimpleCuts);
  fChainData->SetBranchAddress("Ele2PassTightSimpleCuts", &Ele2PassTightSimpleCuts);
  fChainData->SetBranchAddress("Ele2PassHZZICHEP2012", &Ele2PassHZZICHEP2012);

  
  
  fChainMC->SetBranchAddress("event", &event);
  fChainMC->SetBranchAddress("weight", &weight);
  fChainMC->SetBranchAddress("run", &run);
  fChainMC->SetBranchAddress("rho", &rho);
  fChainMC->SetBranchAddress("mass", &mass);
  fChainMC->SetBranchAddress("Ele1Eta", &Ele1Eta);
  fChainMC->SetBranchAddress("Ele1Phi", &Ele1Phi);
  fChainMC->SetBranchAddress("Ele1Pt", &Ele1Pt);
  fChainMC->SetBranchAddress("Ele1SCEta", &Ele1SCEta);
  fChainMC->SetBranchAddress("Ele1SCPhi", &Ele1SCPhi);
  fChainMC->SetBranchAddress("Ele1Energy", &Ele1Energy);
  fChainMC->SetBranchAddress("Ele1EnergyRegression", &Ele1EnergyRegression);
  fChainMC->SetBranchAddress("Ele1EnergyRegressionV0", &Ele1EnergyRegressionV0);
  fChainMC->SetBranchAddress("Ele1EnergyRegressionV1", &Ele1EnergyRegressionV1);
  fChainMC->SetBranchAddress("Ele1EnergyRegressionV2", &Ele1EnergyRegressionV2);
  fChainMC->SetBranchAddress("Ele1EnergyRegressionErrorV0", &Ele1EnergyRegressionErrorV0);
  fChainMC->SetBranchAddress("Ele1EnergyRegressionErrorV1", &Ele1EnergyRegressionErrorV1);
  fChainMC->SetBranchAddress("Ele1EnergyRegressionErrorV2", &Ele1EnergyRegressionErrorV2);
  fChainMC->SetBranchAddress("Ele1HZZICHEP2012IDMVA", &Ele1HZZICHEP2012IDMVA);
  fChainMC->SetBranchAddress("Ele1PFIso04", &Ele1PFIso04);
  fChainMC->SetBranchAddress("Ele1R9", &Ele1R9);
  fChainMC->SetBranchAddress("Ele1PassLooseSimpleCuts", &Ele1PassLooseSimpleCuts);
  fChainMC->SetBranchAddress("Ele1PassMediumSimpleCuts", &Ele1PassMediumSimpleCuts);
  fChainMC->SetBranchAddress("Ele1PassTightSimpleCuts", &Ele1PassTightSimpleCuts);
  fChainMC->SetBranchAddress("Ele1PassHZZICHEP2012", &Ele1PassHZZICHEP2012);
  fChainMC->SetBranchAddress("Ele2Pt", &Ele2Pt);
  fChainMC->SetBranchAddress("Ele2Eta", &Ele2Eta);
  fChainMC->SetBranchAddress("Ele2Phi", &Ele2Phi);
  fChainMC->SetBranchAddress("Ele2SCEt", &Ele2SCEt);
  fChainMC->SetBranchAddress("Ele2SCEta", &Ele2SCEta);
  fChainMC->SetBranchAddress("Ele2SCPhi", &Ele2SCPhi);
  fChainMC->SetBranchAddress("Ele2Energy", &Ele2Energy);
  fChainMC->SetBranchAddress("Ele2EnergyRegression", &Ele2EnergyRegression);
  fChainMC->SetBranchAddress("Ele2EnergyRegressionV0", &Ele2EnergyRegressionV0);
  fChainMC->SetBranchAddress("Ele2EnergyRegressionV1", &Ele2EnergyRegressionV1);
  fChainMC->SetBranchAddress("Ele2EnergyRegressionV2", &Ele2EnergyRegressionV2);
  fChainMC->SetBranchAddress("Ele2EnergyRegressionErrorV0", &Ele2EnergyRegressionErrorV0);
  fChainMC->SetBranchAddress("Ele2EnergyRegressionErrorV1", &Ele2EnergyRegressionErrorV1);
  fChainMC->SetBranchAddress("Ele2EnergyRegressionErrorV2", &Ele2EnergyRegressionErrorV2);
  fChainMC->SetBranchAddress("Ele2HZZICHEP2012IDMVA", &Ele2HZZICHEP2012IDMVA);
  fChainMC->SetBranchAddress("Ele2PFIso04", &Ele2PFIso04);
  fChainMC->SetBranchAddress("Ele2R9", &Ele2R9);
  fChainMC->SetBranchAddress("Ele2PassLooseSimpleCuts", &Ele2PassLooseSimpleCuts);
  fChainMC->SetBranchAddress("Ele2PassMediumSimpleCuts", &Ele2PassMediumSimpleCuts);
  fChainMC->SetBranchAddress("Ele2PassTightSimpleCuts", &Ele2PassTightSimpleCuts);
  fChainMC->SetBranchAddress("Ele2PassHZZICHEP2012", &Ele2PassHZZICHEP2012);

}
