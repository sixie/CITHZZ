

//======================================================== //
// For ECAL with CMS Detector at LHC                       //
// Roofit Macro for Unbinned fit to Z peak                 //
//======================================================== //

#ifndef __CINT__
#include<stdio.h>
#include<string>
#include<sstream> 
#include<iostream>
#include<fstream>
#endif

#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
// #include "TROOT.h"
// #include "TStopwatch.h"
// #include "TStyle.h"
// #include "TLatex.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "CITCommon/CommonData/interface/ZeeEventTree.h"
#include "CITHZZ/CommonCode/CommonDefs.hh"
#include "tdrstyle.C"
#include "CITCommon/Utils/ElectronEnergyCorrection.hh"


//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}


void DoValidateZMassScaleAndResolution(string datafilename, string mcfilename,
                                       string Label, 
                                       Int_t EnergyType,
                                       Bool_t applyCorrections
  ) {
  
  //*************************************************************************
  //setup
  //*************************************************************************
  string label = Label;
  if (label != "") label = "_"+Label;
  Double_t minMass = 70;
  Double_t maxMass = 110;
  TRandom3 *random = new TRandom3();
  
  //*************************************************************************
  // Histograms
  //*************************************************************************
  const UInt_t ncategory = 6;
  TH1F *ZeeMass_Data[ncategory];
  TH1F *ZeeMass_MC[ncategory];
  for (uint i=0; i < ncategory; ++i) {
    ZeeMass_Data[i] = new TH1F( Form("ZeeMass_Data_Category%d",i), " ; m_{ee} [GeV/c^{2}] ; Fraction of Events", 40, 70, 110);
    ZeeMass_MC[i] = new TH1F( Form("ZeeMass_MC_Category%d",i), " ; m_{ee} [GeV/c^{2}] ; Fraction of Events", 40, 70, 110);
  }



  //*************************************************************************
  //Load input files
  //*************************************************************************
  
  citana::ZeeEventTree *dataTree = new citana::ZeeEventTree();
  dataTree->LoadTree(datafilename.c_str(), citana::ZeeEventTree::kCITZeeEvent);
  citana::ZeeEventTree *mcTree = new citana::ZeeEventTree();
  mcTree->LoadTree(mcfilename.c_str(), citana::ZeeEventTree::kCITZeeEvent);
  

  //*************************************************************************
  //Loop over data events
  //*************************************************************************
  for (int i = 0; i < dataTree->tree_->GetEntries(); i++) {
    if (i % 100000 == 0) cout << "Event " << i << endl;
    dataTree->tree_->GetEntry(i);
    
    //*************************************************************************
    //Electron Selection
    //*************************************************************************
    if (!(dataTree->fEle1PassHZZICHEP2012 == 1 && dataTree->fEle2PassHZZICHEP2012 == 1)) continue;

    //*************************************************************************
    //Compute electron four vector;
    //*************************************************************************
    double ele1pt = dataTree->fEle1Pt;
    double ele2pt = dataTree->fEle2Pt;
    if (EnergyType == 1) {
      ele1pt = dataTree->fEle1EnergyRegressionV0 / TMath::CosH(dataTree->fEle1Eta);
      ele2pt = dataTree->fEle2EnergyRegressionV0 / TMath::CosH(dataTree->fEle2Eta);
    }
    else if (EnergyType == 2) {
      ele1pt = dataTree->fEle1EnergyRegressionV1 / TMath::CosH(dataTree->fEle1Eta);
      ele2pt = dataTree->fEle2EnergyRegressionV1 / TMath::CosH(dataTree->fEle2Eta);
    }
    else if (EnergyType == 3) {
      ele1pt = dataTree->fEle1EnergyRegressionV2 / TMath::CosH(dataTree->fEle1Eta);
      ele2pt = dataTree->fEle2EnergyRegressionV2 / TMath::CosH(dataTree->fEle2Eta);
    }
    //*************************************************************************
    //perform energy scale and resolution correction
    //*************************************************************************
    Double_t ele1ptCorr = ele1pt;
    Double_t ele2ptCorr = ele2pt;

    if (applyCorrections) {
      Double_t ele1pCorr = ElectronEnergyScaleAndResolutionCorrection( random, ele1pt*TMath::CosH(dataTree->fEle1Eta), 
                                                                       dataTree->fEle1SCEta, kFALSE, EnergyType, "Prompt" );
      ele1ptCorr = ele1pCorr / TMath::CosH(dataTree->fEle1Eta);
      Double_t ele2pCorr = ElectronEnergyScaleAndResolutionCorrection( random, ele2pt*TMath::CosH(dataTree->fEle2Eta), 
                                                                     dataTree->fEle2SCEta, kFALSE, EnergyType, "Prompt" );
      ele2ptCorr = ele2pCorr / TMath::CosH(dataTree->fEle2Eta);
    }

    //*************************************************************************
    //set four vectors
    //*************************************************************************
    TLorentzVector ele1FourVector;
    ele1FourVector.SetPtEtaPhiM(ele1ptCorr, dataTree->fEle1Eta, dataTree->fEle1Phi, ELECTRONMASS);
    TLorentzVector ele2FourVector;
    ele2FourVector.SetPtEtaPhiM(ele2ptCorr, dataTree->fEle2Eta, dataTree->fEle2Phi, ELECTRONMASS);
    

    //*************************************************************************
    //pt and eta cuts on electron
    //*************************************************************************
    if (! (ele1ptCorr > 20 && ele2ptCorr > 20
           && fabs( dataTree->fEle1Eta) < 2.5 
           && fabs( dataTree->fEle2Eta) < 2.5 )) continue;


    //*************************************************************************
    //pt bins and eta bins
    //*************************************************************************

    Int_t Ele1EtaBin = -1;
    Int_t Ele2EtaBin = -1;
    if (fabs(dataTree->fEle1SCEta) < 1.0) Ele1EtaBin = 0;
    else if (fabs(dataTree->fEle1SCEta) < 1.479) Ele1EtaBin = 1;
    else Ele1EtaBin = 2;
    if (fabs(dataTree->fEle2SCEta) < 1.0) Ele2EtaBin = 0;
    else if (fabs(dataTree->fEle2SCEta) < 1.479) Ele2EtaBin = 1;
    else Ele2EtaBin = 2;

    UInt_t CategoryBin = -1;
    if (Ele1EtaBin == 0 && Ele2EtaBin == 0) CategoryBin = 0;
    else if ((Ele1EtaBin == 0 && Ele2EtaBin == 1)
             || (Ele1EtaBin == 1 && Ele2EtaBin == 0)
      ) CategoryBin = 1;
    else if ((Ele1EtaBin == 0 && Ele2EtaBin == 2)
             || (Ele1EtaBin == 2 && Ele2EtaBin == 0)
      ) CategoryBin = 2;
    else if (Ele1EtaBin == 1 && Ele2EtaBin == 1)
      CategoryBin = 3;
    else if ((Ele1EtaBin == 1 && Ele2EtaBin == 2)
             || (Ele1EtaBin == 2 && Ele2EtaBin == 1)
      ) CategoryBin = 4;
    else if (Ele1EtaBin == 2 && Ele2EtaBin == 2)
      CategoryBin = 5;
    if (CategoryBin == -1) assert(false);

    //*************************************************************************
    // restrict range of mass
    //*************************************************************************
    double zMass = (ele1FourVector+ele2FourVector).M();
    if (zMass < minMass || zMass > maxMass) continue;

    //*************************************************************************
    //fill histogram    
    //*************************************************************************
    ZeeMass_Data[CategoryBin]->Fill(zMass);
  }



  //*************************************************************************
  //Loop over mc events
  //*************************************************************************
  for (int i = 0; i < mcTree->tree_->GetEntries(); i++) {
    if (i % 100000 == 0) cout << "Event " << i << endl;
    mcTree->tree_->GetEntry(i);

    //*************************************************************************
    //Electron Selection
    //*************************************************************************
    if (!(mcTree->fEle1PassHZZICHEP2012 == 1 && mcTree->fEle2PassHZZICHEP2012 == 1)) continue;

    //*************************************************************************
    //Compute electron four vector;
    //*************************************************************************
    double ele1pt = mcTree->fEle1Pt;
    double ele2pt = mcTree->fEle2Pt;
    if (EnergyType == 1) {
      ele1pt = mcTree->fEle1EnergyRegressionV0 / TMath::CosH(mcTree->fEle1Eta);
      ele2pt = mcTree->fEle2EnergyRegressionV0 / TMath::CosH(mcTree->fEle2Eta);
    }
    else if (EnergyType == 2) {
      ele1pt = mcTree->fEle1EnergyRegressionV1 / TMath::CosH(mcTree->fEle1Eta);
      ele2pt = mcTree->fEle2EnergyRegressionV1 / TMath::CosH(mcTree->fEle2Eta);
    }
    else if (EnergyType == 3) {
      ele1pt = mcTree->fEle1EnergyRegressionV2 / TMath::CosH(mcTree->fEle1Eta);
      ele2pt = mcTree->fEle2EnergyRegressionV2 / TMath::CosH(mcTree->fEle2Eta);
    }

    //*************************************************************************
    //perform energy scale and resolution correction
    //*************************************************************************
    Double_t ele1ptCorr = ele1pt;
    Double_t ele2ptCorr = ele2pt;

    if (applyCorrections) {
      Double_t ele1pCorr = ElectronEnergyScaleAndResolutionCorrection( random, ele1pt*TMath::CosH(mcTree->fEle1Eta), 
                                                                       mcTree->fEle1SCEta, kTRUE, EnergyType, "Summer12" );
      ele1ptCorr = ele1pCorr / TMath::CosH(mcTree->fEle1Eta);
      Double_t ele2pCorr = ElectronEnergyScaleAndResolutionCorrection( random, ele2pt*TMath::CosH(mcTree->fEle2Eta), 
                                                                       mcTree->fEle2SCEta, kTRUE, EnergyType, "Summer12");
      ele2ptCorr = ele2pCorr / TMath::CosH(mcTree->fEle2Eta);
    }

    //*************************************************************************
    //set four vectors
    //*************************************************************************
    TLorentzVector ele1FourVector;
    ele1FourVector.SetPtEtaPhiM(ele1ptCorr, mcTree->fEle1Eta, mcTree->fEle1Phi, ELECTRONMASS);
    TLorentzVector ele2FourVector;
    ele2FourVector.SetPtEtaPhiM(ele2ptCorr, mcTree->fEle2Eta, mcTree->fEle2Phi, ELECTRONMASS);
    

    //*************************************************************************
    //pt and eta cuts on electron
    //*************************************************************************
    if (! (ele1ptCorr > 20 && ele2ptCorr > 20
           && fabs( mcTree->fEle1Eta) < 2.5 
           && fabs( mcTree->fEle2Eta) < 2.5 )) continue;


    //*************************************************************************
    //pt bins and eta bins
    //*************************************************************************

    Int_t Ele1EtaBin = -1;
    Int_t Ele2EtaBin = -1;
    if (fabs(mcTree->fEle1SCEta) < 1.0) Ele1EtaBin = 0;
    else if (fabs(mcTree->fEle1SCEta) < 1.479) Ele1EtaBin = 1;
    else Ele1EtaBin = 2;
    if (fabs(mcTree->fEle2SCEta) < 1.0) Ele2EtaBin = 0;
    else if (fabs(mcTree->fEle2SCEta) < 1.479) Ele2EtaBin = 1;
    else Ele2EtaBin = 2;

    UInt_t CategoryBin = -1;
    if (Ele1EtaBin == 0 && Ele2EtaBin == 0) CategoryBin = 0;
    else if ((Ele1EtaBin == 0 && Ele2EtaBin == 1)
             || (Ele1EtaBin == 1 && Ele2EtaBin == 0)
      ) CategoryBin = 1;
    else if ((Ele1EtaBin == 0 && Ele2EtaBin == 2)
             || (Ele1EtaBin == 2 && Ele2EtaBin == 0)
      ) CategoryBin = 2;
    else if (Ele1EtaBin == 1 && Ele2EtaBin == 1)
      CategoryBin = 3;
    else if ((Ele1EtaBin == 1 && Ele2EtaBin == 2)
             || (Ele1EtaBin == 2 && Ele2EtaBin == 1)
      ) CategoryBin = 4;
    else if (Ele1EtaBin == 2 && Ele2EtaBin == 2)
      CategoryBin = 5;
    if (CategoryBin == -1) assert(false);

    //*************************************************************************
    // restrict range of mass
    //*************************************************************************
    double zMass = (ele1FourVector+ele2FourVector).M();
    if (zMass < minMass || zMass > maxMass) continue;

    //*************************************************************************
    //fill histogram    
    //*************************************************************************
    ZeeMass_MC[CategoryBin]->Fill(zMass);
  }


  //*************************************************************************
  //save histograms to file 
  //*************************************************************************
  TFile *outputfile = new TFile(Form("ZeeMassScaleAndResolutionValidationPlots_EnergyType%d_ApplyCorr%d.root",EnergyType, applyCorrections), "UPDATE");
  for (uint i=0; i < ncategory; ++i) {
    outputfile->WriteTObject(ZeeMass_Data[i], ZeeMass_Data[i]->GetName(), "WriteDelete");
    outputfile->WriteTObject(ZeeMass_MC[i], ZeeMass_MC[i]->GetName(), "WriteDelete");
  }
  outputfile->Close();
  delete outputfile;

}



void CompareZMassScaleAndResolution(string datafilename, 
                                    string Label
  ) {
  
  //*************************************************************************
  //setup
  //*************************************************************************
  string label = Label;
  if (label != "") label = "_"+Label;
  Double_t minMass = 70;
  Double_t maxMass = 110;
  TRandom3 *random = new TRandom3();
  
  //*************************************************************************
  // Histograms
  //*************************************************************************
  const UInt_t ncategory = 6;
  TH1F *ZeeMass[ncategory];
  TH1F *ZeeMassRegressionV0[ncategory];
  TH1F *ZeeMassRegressionV1[ncategory];
  TH1F *ZeeMassRegressionV2[ncategory];
  for (uint i=0; i < ncategory; ++i) {
    ZeeMass[i] = new TH1F( Form("ZeeMass_%s_Category%d",Label.c_str(),i), " ; m_{ee} [GeV/c^{2}] ; Fraction of Events", 40, 70, 110);
    ZeeMassRegressionV0[i] = new TH1F( Form("ZeeMassRegressionV0_%s_Category%d",Label.c_str(),i), " ; m_{ee} [GeV/c^{2}] ; Fraction of Events", 40, 70, 110);
    ZeeMassRegressionV1[i] = new TH1F( Form("ZeeMassRegressionV1_%s_Category%d",Label.c_str(),i), " ; m_{ee} [GeV/c^{2}] ; Fraction of Events", 40, 70, 110);
    ZeeMassRegressionV2[i] = new TH1F( Form("ZeeMassRegressionV2_%s_Category%d",Label.c_str(),i), " ; m_{ee} [GeV/c^{2}] ; Fraction of Events", 40, 70, 110);
  }



  //*************************************************************************
  //Load input files
  //*************************************************************************
  
  citana::ZeeEventTree *dataTree = new citana::ZeeEventTree();
  dataTree->LoadTree(datafilename.c_str(), citana::ZeeEventTree::kCITZeeEvent);

  //*************************************************************************
  //Loop over data events
  //*************************************************************************
  for (int i = 0; i < dataTree->tree_->GetEntries(); i++) {
    if (i % 100000 == 0) cout << "Event " << i << endl;
    dataTree->tree_->GetEntry(i);
    
    //*************************************************************************
    //Electron Selection
    //*************************************************************************
    if (!(dataTree->fEle1PassHZZICHEP2012 == 1 && dataTree->fEle2PassHZZICHEP2012 == 1)) continue;

    //*************************************************************************
    //Compute electron four vector;
    //*************************************************************************

    double ele1pt = dataTree->fEle1Pt;
    double ele2pt = dataTree->fEle2Pt;
    double ele1pt_RegressionV0 = dataTree->fEle1EnergyRegressionV0 / TMath::CosH(dataTree->fEle1Eta);
    double ele2pt_RegressionV0 = dataTree->fEle2EnergyRegressionV0 / TMath::CosH(dataTree->fEle2Eta);
    double ele1pt_RegressionV1 = dataTree->fEle1EnergyRegressionV1 / TMath::CosH(dataTree->fEle1Eta);
    double ele2pt_RegressionV1 = dataTree->fEle2EnergyRegressionV1 / TMath::CosH(dataTree->fEle2Eta);
    double ele1pt_RegressionV2 = dataTree->fEle1EnergyRegressionV2 / TMath::CosH(dataTree->fEle1Eta);
    double ele2pt_RegressionV2 = dataTree->fEle2EnergyRegressionV2 / TMath::CosH(dataTree->fEle2Eta);

    //*************************************************************************
    //set four vectors
    //*************************************************************************
    TLorentzVector ele1FourVector;
    ele1FourVector.SetPtEtaPhiM(ele1pt, dataTree->fEle1Eta, dataTree->fEle1Phi, ELECTRONMASS);
    TLorentzVector ele2FourVector;
    ele2FourVector.SetPtEtaPhiM(ele2pt, dataTree->fEle2Eta, dataTree->fEle2Phi, ELECTRONMASS);

    TLorentzVector ele1FourVectorRegressionV0;
    ele1FourVectorRegressionV0.SetPtEtaPhiM(ele1pt_RegressionV0, dataTree->fEle1Eta, dataTree->fEle1Phi, ELECTRONMASS);
    TLorentzVector ele2FourVectorRegressionV0;
    ele2FourVectorRegressionV0.SetPtEtaPhiM(ele2pt_RegressionV0, dataTree->fEle2Eta, dataTree->fEle2Phi, ELECTRONMASS);
    TLorentzVector ele1FourVectorRegressionV1;
    ele1FourVectorRegressionV1.SetPtEtaPhiM(ele1pt_RegressionV1, dataTree->fEle1Eta, dataTree->fEle1Phi, ELECTRONMASS);
    TLorentzVector ele2FourVectorRegressionV1;
    ele2FourVectorRegressionV1.SetPtEtaPhiM(ele2pt_RegressionV1, dataTree->fEle2Eta, dataTree->fEle2Phi, ELECTRONMASS);
    TLorentzVector ele1FourVectorRegressionV2;
    ele1FourVectorRegressionV2.SetPtEtaPhiM(ele1pt_RegressionV2, dataTree->fEle1Eta, dataTree->fEle1Phi, ELECTRONMASS);
    TLorentzVector ele2FourVectorRegressionV2;
    ele2FourVectorRegressionV2.SetPtEtaPhiM(ele2pt_RegressionV2, dataTree->fEle2Eta, dataTree->fEle2Phi, ELECTRONMASS);
    

    //*************************************************************************
    //pt and eta cuts on electron
    //*************************************************************************
    if (! (ele1pt > 20 && ele2pt > 20
           && fabs( dataTree->fEle1Eta) < 2.5 
           && fabs( dataTree->fEle2Eta) < 2.5 )) continue;


    //*************************************************************************
    //pt bins and eta bins
    //*************************************************************************

    Int_t Ele1EtaBin = -1;
    Int_t Ele2EtaBin = -1;
    if (fabs(dataTree->fEle1SCEta) < 1.0) Ele1EtaBin = 0;
    else if (fabs(dataTree->fEle1SCEta) < 1.479) Ele1EtaBin = 1;
    else Ele1EtaBin = 2;
    if (fabs(dataTree->fEle2SCEta) < 1.0) Ele2EtaBin = 0;
    else if (fabs(dataTree->fEle2SCEta) < 1.479) Ele2EtaBin = 1;
    else Ele2EtaBin = 2;

    UInt_t CategoryBin = -1;
    if (Ele1EtaBin == 0 && Ele2EtaBin == 0) CategoryBin = 0;
    else if ((Ele1EtaBin == 0 && Ele2EtaBin == 1)
             || (Ele1EtaBin == 1 && Ele2EtaBin == 0)
      ) CategoryBin = 1;
    else if ((Ele1EtaBin == 0 && Ele2EtaBin == 2)
             || (Ele1EtaBin == 2 && Ele2EtaBin == 0)
      ) CategoryBin = 2;
    else if (Ele1EtaBin == 1 && Ele2EtaBin == 1)
      CategoryBin = 3;
    else if ((Ele1EtaBin == 1 && Ele2EtaBin == 2)
             || (Ele1EtaBin == 2 && Ele2EtaBin == 1)
      ) CategoryBin = 4;
    else if (Ele1EtaBin == 2 && Ele2EtaBin == 2)
      CategoryBin = 5;
    if (CategoryBin == -1) assert(false);

    //*************************************************************************
    // restrict range of mass
    //*************************************************************************
    double zMass = (ele1FourVector+ele2FourVector).M();
    double zMassRegressionV0 = (ele1FourVectorRegressionV0+ele2FourVectorRegressionV0).M();
    double zMassRegressionV1 = (ele1FourVectorRegressionV1+ele2FourVectorRegressionV1).M();
    double zMassRegressionV2 = (ele1FourVectorRegressionV2+ele2FourVectorRegressionV2).M();
    
    if (!(zMass < minMass || zMass > maxMass)) {
      ZeeMass[CategoryBin]->Fill(zMass);
      ZeeMassRegressionV0[CategoryBin]->Fill(zMassRegressionV0);
      ZeeMassRegressionV1[CategoryBin]->Fill(zMassRegressionV1);
      ZeeMassRegressionV2[CategoryBin]->Fill(zMassRegressionV2);
    }

    //*************************************************************************
    //fill histogram    
    //*************************************************************************
  }


  //*************************************************************************
  //save histograms to file 
  //*************************************************************************
  TFile *outputfile = new TFile("ZeeMassComparisonPlots.root", "UPDATE");
  for (uint i=0; i < ncategory; ++i) {
    outputfile->WriteTObject(ZeeMass[i], ZeeMass[i]->GetName(), "WriteDelete");
    outputfile->WriteTObject(ZeeMassRegressionV0[i], ZeeMassRegressionV0[i]->GetName(), "WriteDelete");
    outputfile->WriteTObject(ZeeMassRegressionV1[i], ZeeMassRegressionV1[i]->GetName(), "WriteDelete");
    outputfile->WriteTObject(ZeeMassRegressionV2[i], ZeeMassRegressionV2[i]->GetName(), "WriteDelete");
  }
  outputfile->Close();
  delete outputfile;

}



void MakePlots( Int_t EnergyType, Bool_t applyCorrections ) {

  //*************************************************************************
  //Load Plots
  //*************************************************************************
  const UInt_t ncategory = 6;
  TH1F *ZeeMass_Data[ncategory];
  TH1F *ZeeMassRegressionV0_Data[ncategory];
  TH1F *ZeeMassRegressionV1_Data[ncategory];
  TH1F *ZeeMassRegressionV2_Data[ncategory];
  TH1F *ZeeMass_MC[ncategory];


//   TFile *file = new TFile(Form("ZeeMassScaleAndResolutionValidationPlots_EnergyType%d_ApplyCorr%d.root",EnergyType, applyCorrections), "READ");
//   for (uint i=0; i < ncategory; ++i) {
//     ZeeMass_Data[i] = (TH1F*)file->Get(Form("ZeeMass_Data_Category%d",i));
//     ZeeMass_MC[i] = (TH1F*)file->Get(Form("ZeeMass_MC_Category%d",i));
//   }

  TFile *file = new TFile("ZeeMassComparisonPlots.root", "READ");
  for (uint i=0; i < ncategory; ++i) {
    ZeeMass_Data[i] = (TH1F*)file->Get(Form("ZeeMass_Data_Category%d",i));
    ZeeMassRegressionV0_Data[i] = (TH1F*)file->Get(Form("ZeeMassRegressionV0_Data_Category%d",i));
    ZeeMassRegressionV1_Data[i] = (TH1F*)file->Get(Form("ZeeMassRegressionV1_Data_Category%d",i));
    ZeeMassRegressionV2_Data[i] = (TH1F*)file->Get(Form("ZeeMassRegressionV2_Data_Category%d",i));
    ZeeMass_MC[i] = (TH1F*)file->Get(Form("ZeeMass_MC_Category%d",i));
  }


  //*************************************************************************
  //plot
  //*************************************************************************
  TCanvas* cv;
  TLegend* legend;

  
//   //*************************************************************************
//   //Data vs MC
//   //*************************************************************************
//   for (uint i=0; i < ncategory; ++i) {
    
//     cv = new TCanvas("cv","cv", 800,600);
//     legend = new TLegend( 0.25, 0.7, 0.5, 0.9);
//     legend->SetTextSize(0.04);
//     legend->SetBorderSize(0);
//     legend->SetFillStyle(0);
//     legend->AddEntry(ZeeMass_Data[i], "Data", "L");
//     legend->AddEntry(ZeeMass_MC[i], "Monte Carlo", "L");

//     NormalizeHist(ZeeMass_Data[i]);
//     NormalizeHist(ZeeMass_MC[i]);

//     ZeeMass_Data[i]->SetTitle("");
//     ZeeMass_Data[i]->SetMarkerSize(0.75);
//     ZeeMass_Data[i]->SetMarkerColor(kBlue);
//     ZeeMass_Data[i]->GetYaxis()->SetTitleOffset(1.2);
//     ZeeMass_Data[i]->GetXaxis()->SetTitleOffset(1.05);
//     ZeeMass_Data[i]->Draw("e1");
//     ZeeMass_Data[i]->SetLineColor(kBlue);
//     ZeeMass_Data[i]->SetLineWidth(2);
//     ZeeMass_MC[i]->Draw("hist,same");
//     ZeeMass_MC[i]->SetLineColor(kRed);
//     ZeeMass_MC[i]->SetLineWidth(2);
// //     ZeeMass_Data[i]->SetMaximum(1.2*max(ZeeMass_Data[i]->GetMaximum(),ZeeMass_MC[i]->GetMaximum()));
//     ZeeMass_Data[i]->SetMaximum(0.18);
//     legend->Draw();  
//     if (applyCorrections) {
//       cv->SaveAs(Form("ZeeMassScaleAndResolutionValidationPlots_EnergyType%d_AfterCorr_CategoryBin%d.gif",EnergyType,i));
//     } else {
//       cv->SaveAs(Form("ZeeMassScaleAndResolutionValidationPlots_EnergyType%d_BeforeCorr_CategoryBin%d.gif",EnergyType,i));
//     }

//   }


  gStyle->SetOptStat(0);

  //*************************************************************************
  //Standard Vs Regression
  //*************************************************************************
  for (uint i=0; i < ncategory; ++i) {
    
    cv = new TCanvas("cv","cv", 800,600);
    legend = new TLegend( 0.25, 0.7, 0.5, 0.9);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(ZeeMass_Data[i], "Standard", "L");
    legend->AddEntry(ZeeMassRegressionV0_Data[i], "Regression", "L");

    NormalizeHist(ZeeMass_Data[i]);
    NormalizeHist(ZeeMassRegressionV0_Data[i]);

    ZeeMass_Data[i]->SetTitle("");
    ZeeMass_Data[i]->SetMarkerSize(0.75);
    ZeeMass_Data[i]->SetMarkerColor(kBlue);
    ZeeMass_Data[i]->GetYaxis()->SetTitleOffset(1.2);
    ZeeMass_Data[i]->GetXaxis()->SetTitleOffset(1.05);
    ZeeMass_Data[i]->Draw("hist");
    ZeeMass_Data[i]->SetLineColor(kBlue);
    ZeeMass_Data[i]->SetLineWidth(2);
    ZeeMassRegressionV0_Data[i]->Draw("hist,same");
    ZeeMassRegressionV0_Data[i]->SetLineColor(kRed);
    ZeeMassRegressionV0_Data[i]->SetLineWidth(2);
    ZeeMass_Data[i]->SetMaximum(0.18);
    legend->Draw();  
    cv->SaveAs(Form("ZeeMassComparison_CategoryBin%d.gif",i));

  }





}

void MakeComparisonPlots( string Label ) {

  //*************************************************************************
  //Load Plots
  //*************************************************************************
  const UInt_t ncategory = 6;
  TH1F *ZeeMass[ncategory];
  TH1F *ZeeMassRegressionV0[ncategory];
  TH1F *ZeeMassRegressionV1[ncategory];
  TH1F *ZeeMassRegressionV2[ncategory];

  TFile *file = new TFile("ZeeMassComparisonPlots.root", "READ");
  for (uint i=0; i < ncategory; ++i) {
    ZeeMass[i] = (TH1F*)file->Get(Form("ZeeMass_%s_Category%d",Label.c_str(),i));
    ZeeMassRegressionV0[i] = (TH1F*)file->Get(Form("ZeeMassRegressionV0_%s_Category%d",Label.c_str(),i));
    ZeeMassRegressionV1[i] = (TH1F*)file->Get(Form("ZeeMassRegressionV1_%s_Category%d",Label.c_str(),i));
    ZeeMassRegressionV2[i] = (TH1F*)file->Get(Form("ZeeMassRegressionV2_%s_Category%d",Label.c_str(),i));
  }


  //*************************************************************************
  //plot
  //*************************************************************************
  TCanvas* cv;
  TLegend* legend;

  gStyle->SetOptStat(0);

  //*************************************************************************
  //Standard Vs Regression
  //*************************************************************************
  for (uint i=0; i < ncategory; ++i) {
    
    cv = new TCanvas("cv","cv", 800,600);
    legend = new TLegend( 0.25, 0.7, 0.5, 0.9);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(ZeeMass[i], "Standard", "L");
    legend->AddEntry(ZeeMassRegressionV0[i], "Regression", "L");

    NormalizeHist(ZeeMass[i]);
    NormalizeHist(ZeeMassRegressionV0[i]);

    ZeeMass[i]->SetTitle("");
    ZeeMass[i]->SetMarkerSize(0.75);
    ZeeMass[i]->SetMarkerColor(kBlue);
    ZeeMass[i]->GetYaxis()->SetTitleOffset(1.2);
    ZeeMass[i]->GetXaxis()->SetTitleOffset(1.05);
    ZeeMass[i]->Draw("hist");
    ZeeMass[i]->SetLineColor(kBlue);
    ZeeMass[i]->SetLineWidth(2);
    ZeeMassRegressionV0[i]->Draw("hist,same");
    ZeeMassRegressionV0[i]->SetLineColor(kRed);
    ZeeMassRegressionV0[i]->SetLineWidth(2);
    ZeeMass[i]->SetMaximum(0.18);
    legend->Draw();  
    cv->SaveAs(Form("ZeeMassComparison_CategoryBin%d.gif",i));

  }





}



void ValidateZMassScaleAndResolution() {


  //****************************************************
  //Standard Vs Regression1
  //****************************************************

//   DoValidateZMassScaleAndResolution("/afs/cern.ch/user/s/sixie/work/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Data2012.root", 
//                                     "/afs/cern.ch/user/s/sixie/work/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Summer12DY.root",
//                                     "0",
//                                     0, kFALSE);
 
//   MakePlots(0);



//   DoValidateZMassScaleAndResolution("/afs/cern.ch/user/s/sixie/work/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Data2012.root", 
//                                     "/afs/cern.ch/user/s/sixie/work/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Summer12DY.root",
//                                     "0",
//                                     0, kFALSE);
  
//   MakePlots(0);



  //****************************************************
  //regression type 2
  //****************************************************

//   DoValidateZMassScaleAndResolution("/afs/cern.ch/user/s/sixie/work/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Data2012.root", 
//                                      "/afs/cern.ch/user/s/sixie/work/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Summer12DY.root",
//                                      "Regression2",
//                                     2, kTRUE);

//   MakePlots(2, kTRUE);


//   DoValidateZMassScaleAndResolution("/afs/cern.ch/user/s/sixie/work/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Data2012.root", 
//                                     "/afs/cern.ch/user/s/sixie/work/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Summer12DY.root",
//                                     "Regression2",
//                                     2, kFALSE);

//   MakePlots(2, kFALSE);




  //****************************************************
  //Compare Standard Vs Regression
  //****************************************************
  CompareZMassScaleAndResolution("/afs/cern.ch/user/s/sixie/work/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Summer12DY53X.root", "Summer12DY53X");
  MakeComparisonPlots("Summer12DY53X");
  
//   CompareZMassScaleAndResolution("/afs/cern.ch/user/s/sixie/work/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_HCP2012.root", "HCP2012");
//   MakeComparisonPlots("HCP2012");

}
