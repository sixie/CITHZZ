//================================================================================================
//			trainElectronEnergyRegression.C
//
// Given an appropriate electron ntuple, uses the classes GBRTrees etc. in order to train a regression
// to recover the generated electron energy
//
// USAGE
// trainElectronEnergyRegression(char* trainingFile, char* outWeightFile, char* optionChar)
//
// trainingFile		= ntuple file on which to perform the training
// outWeightFile	= output file, to which it will save the weights in a .root file
// optionChar		= denotes which version of the training is being performed
//
// V00	no pT split	no tracker variables
// V01	no pT split	includes tracker variables
// V10	pT split	no tracker variables
// V11	pT split	includes tracker variables
//________________________________________________________________________________________________

//Need to check out CondFormats/EgammaObjects/interface/GBRForest.h

#include "TFile.h"
#include "TTree.h"
#include "GBRTrainer.h"
#include "GBRForest.h"
#include "GBRApply.h"
#include "TCut.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TChain.h"
#include "TObjArray.h"
#include "TTreeFormula.h" 
#include "TCanvas.h"
#include "Cintex/Cintex.h"

int GetTotalEvents(TChain *chain) {

  int numevents = 0;

  TObjArray *files = chain->GetListOfFiles();
  for (int i=0; i<files->GetEntries(); ++i) {    
    TFile *file = TFile::Open(files->At(i)->GetTitle(),"READ");

    TDirectory *fwkdir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
    TH1D *hevents = (TH1D*)fwkdir->Get("hDAllEvents");

    numevents += hevents->GetEntries();

    file->Close();

  }

  return numevents;

}

void trainElectronEnergyRegression(char* trainingFile, char* outWeightFile, char* optionChar, Bool_t restrictNEvents = false) {
  
  // Setting up training option
  std::string optionStr(optionChar);

  // ******** If option is V00, V01, V02, etc. ********* //
  if (optionStr == "V00" || optionStr == "V01" || optionStr == "V02") {

    GBRTrainer *traineb = new GBRTrainer;
    GBRTrainer *trainebvar = new GBRTrainer;
    GBRTrainer *trainee = new GBRTrainer;
    GBRTrainer *traineevar = new GBRTrainer;

    TTree *intree = 0;

    cout << "Training on file " << trainingFile << " with version " << optionChar << endl;
//     TChain *chainele = new TChain("Electrons");
    TChain *chainele = new TChain("T1");
    chainele->Add(trainingFile);
    chainele->LoadTree(0);    
    chainele->SetCacheSize(64*1024*1024);
    chainele->SetCacheLearnEntries();
    intree = chainele;

    traineb->AddTree(chainele);
    trainebvar->AddTree(chainele);
    trainee->AddTree(chainele);
    traineevar->AddTree(chainele);

    TCut traincut = "spp == spp && sep == sep && GeneratedEnergyStatus3 >= GeneratedEnergyStatus1 && (GeneratedEnergyStatus3 - GeneratedEnergyStatus1)/GeneratedEnergyStatus3 < 0.01 && pt>5";

    TCut evtcut;
    TCut evtcutvar;

    //Do this if we only want to train for energy 
    //   evtcut = "evt_num%2==0 ";
    //   evtcutvar = "evt_num%2==0 ";

    //Do this if we also want to train for the per electron energy variance
    if (!restrictNEvents) {
      evtcut = "event%2==0 && (event/2)%4!=3";
      evtcutvar = "event%2==0 && (event/2)%4==3";
    } else {
      evtcut = "event%2==0 && (event/2)%4!=3 && event < 30000000";
      evtcutvar = "event%2==0 && (event/2)%4==3 && event < 30000000";
    }

    traineb->SetTrainingCut(std::string(traincut && evtcut && "abs(scEta)<1.479"));
    trainee->SetTrainingCut(std::string(traincut && evtcut && "abs(scEta)>=1.479"));
    //turn this off for now
    trainebvar->SetTrainingCut(std::string(traincut && evtcutvar && "IsEB"));
    traineevar->SetTrainingCut(std::string(traincut && evtcutvar && "IsEE"));

    const double maxsig = 3.0;
    const double shrinkage = 0.1;

    traineb->SetMinEvents(200);
    traineb->SetShrinkage(shrinkage);  
    traineb->SetMinCutSignificance(maxsig);

    trainebvar->SetMinEvents(200);
    trainebvar->SetShrinkage(shrinkage);  
    trainebvar->SetMinCutSignificance(maxsig);  

    trainee->SetMinEvents(200);
    trainee->SetShrinkage(shrinkage);  
    trainee->SetMinCutSignificance(maxsig);  

    traineevar->SetMinEvents(200);
    traineevar->SetShrinkage(shrinkage);  
    traineevar->SetMinCutSignificance(maxsig);    

    traineb->SetTargetVar("GeneratedEnergyStatus1/SCRawEnergy");
    trainebvar->SetTargetVar("1.253*abs( targeteb - GeneratedEnergyStatus1/SCRawEnergy) ");
    trainee->SetTargetVar("GeneratedEnergyStatus1/(SCRawEnergy*(1+PreShowerOverRaw))");
    traineevar->SetTargetVar("1.253*abs( targetee - GeneratedEnergyStatus1/(SCRawEnergy*(1+PreShowerOverRaw)))");

    std::vector<std::string> *varsf = new std::vector<std::string>;
    varsf->push_back("SCRawEnergy");
    varsf->push_back("scEta");
    varsf->push_back("scPhi");
    varsf->push_back("R9");  
    varsf->push_back("E5x5Seed/SCRawEnergy");  
    varsf->push_back("etawidth");
    varsf->push_back("phiwidth");  
    varsf->push_back("NClusters");
    varsf->push_back("HoE");
    varsf->push_back("rho"); 
    varsf->push_back("vertices");  

    varsf->push_back("EtaSeed-scEta");
    varsf->push_back("atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi))");
    varsf->push_back("ESeed/SCRawEnergy");
    varsf->push_back("E3x3Seed/ESeed");
    varsf->push_back("E5x5Seed/ESeed");
    varsf->push_back("see");   
    varsf->push_back("spp");   
    varsf->push_back("sep");
    varsf->push_back("EMaxSeed/ESeed");
    varsf->push_back("E2ndSeed/ESeed");
    varsf->push_back("ETopSeed/ESeed");
    varsf->push_back("EBottomSeed/ESeed");
    varsf->push_back("ELeftSeed/ESeed");
    varsf->push_back("ERightSeed/ESeed");
    varsf->push_back("E2x5MaxSeed/ESeed");
    varsf->push_back("E2x5TopSeed/ESeed");
    varsf->push_back("E2x5BottomSeed/ESeed");
    varsf->push_back("E2x5LeftSeed/ESeed");
    varsf->push_back("E2x5RightSeed/ESeed");
    varsf->push_back("ecaldriven");

    // Also training on tracker variables if option 1 is selected
    if (optionStr == "V01") {
      varsf->push_back("pmodegsf");
      varsf->push_back("fbrem");
      varsf->push_back("Charge");
      varsf->push_back("EoP");
      varsf->push_back("perror/pmodegsf");
      varsf->push_back("ecalenergyerror/SCRawEnergy");
      varsf->push_back("EleClassification");
    }
    if (optionStr == "V02") {
      varsf->push_back("pmodegsf");
      varsf->push_back("fbrem");
      varsf->push_back("Charge");
      varsf->push_back("EoP");
      varsf->push_back("perror/pmodegsf");
      varsf->push_back("ecalenergyerror/SCRawEnergy");
      varsf->push_back("EleClassification");
      varsf->push_back("deta");
      varsf->push_back("dphi");
      varsf->push_back("detacalo");
      varsf->push_back("dphicalo");
      varsf->push_back("gsfchi2");
      varsf->push_back("kflayers");
      varsf->push_back("EEleoPout");
    }
    


    std::vector<std::string> *varseb = new std::vector<std::string>(*varsf);
    std::vector<std::string> *varsee = new std::vector<std::string>(*varsf);

    varseb->push_back("IEtaSeed");
    varseb->push_back("IPhiSeed");
    varseb->push_back("IEtaSeed%5");
    varseb->push_back("IPhiSeed%2");       
    varseb->push_back("(abs(IEtaSeed)<=25)*(IEtaSeed%25) + (abs(IEtaSeed)>25)*((IEtaSeed-25*abs(IEtaSeed)/IEtaSeed)%20)");
    varseb->push_back("IPhiSeed%20"); 
    varseb->push_back("EtaCrySeed");
    varseb->push_back("PhiCrySeed");

    varsee->push_back("PreShowerOverRaw");


    for (int i=0; i<varseb->size(); ++i) {
      cout << "var " << i << " = " << varseb->at(i) << endl;
      traineb->AddInputVar(varseb->at(i));
      trainebvar->AddInputVar(varseb->at(i));
    }

    for (int i=0; i<varsee->size(); ++i) {
      cout << "var " << i << " = " << varsee->at(i) << endl;
      trainee->AddInputVar(varsee->at(i));
      traineevar->AddInputVar(varsee->at(i));
    }

    ROOT::Cintex::Cintex::Enable();   

    //  TFile *ftmp = new TFile("tmpfile.root","RECREATE");    
    GBRApply gbrapply;


    //Train Barrel Energy Regression
    intree->LoadTree(0);  
    const GBRForest *foresteb = traineb->TrainForest(200);
    delete traineb;

    //Apply Barrel Energy Regression
    intree->LoadTree(0);  
    gbrapply.ApplyAsFriend(intree, foresteb, *varseb, "targeteb");


    //Train Barrel Variance Regression
    intree->LoadTree(0);
    const GBRForest *forestebvar = trainebvar->TrainForest(200);
    delete trainebvar;

    //Train Endcap Energy Regression
    intree->LoadTree(0);
    const GBRForest *forestee = trainee->TrainForest(200);
    delete trainee;

    //Apply Endcap Energy Regression
    intree->LoadTree(0);
    gbrapply.ApplyAsFriend(intree, forestee, *varsee, "targetee");

    //Train Endcap Variance Regression
    intree->LoadTree(0);  
    const GBRForest *foresteevar = traineevar->TrainForest(200);
    delete traineevar;  

    TString fname;
    fname = outWeightFile;

    TFile *fout = new TFile(fname,"RECREATE");  
    cout << "Saving weights to file " << outWeightFile << endl;

    fout->WriteObject(foresteb,"EBCorrection");
    fout->WriteObject(forestebvar,"EBUncertainty");
    fout->WriteObject(forestee,"EECorrection");
    fout->WriteObject(foresteevar,"EEUncertainty");

    fout->WriteObject(varseb, "varlisteb");
    fout->WriteObject(varsee, "varlistee");

    //  ftmp->Close();  
    //  fout->Close();
  }  




}
