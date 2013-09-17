//================================================================================================
//                      applyElectronEnergyRegression.C
// 
// Applies the weights created by trainElectronEnergyRegression.C to the desired ntuple file, and
// generates a .root file with "targets" i.e. the ratios to obtain the regression energy, as well
// as .gif plots for diagnostic
//
// USAGE
// applyElectronEnergyRegression(string applyingFile, string weightFilename, string outFileDirectory, string option)
//
// applyingFile		= ntuple file to which apply the regression
// weightFilename	= file outputted by trainElectronEnergyRegression.C containing the weights
// outFileDirectory	= directory to which output the .gif and .root files (must end in /)
// option		= describes which version ("V00" etc.) corresponds to weightFilename (see trainElectronEnergyRegression_option.C for details)
//________________________________________________________________________________________________

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
#include "TCanvas.h"
#include <string.h>


void applyElectronEnergyRegression(string applyingFile, string weightFilename, string outFileDirectory, string option) {
  ROOT::Cintex::Cintex::Enable();   

  if (option == "V00" || option == "V01" || option == "V02" ) {
    //  TFile * weight_file = new TFile("photon_2011_no_tk_vars_energy_reg_weights.root");
    TFile * weight_file = new TFile(weightFilename.c_str());
    const GBRForest * forest_eb_correction = (GBRForest *)weight_file->Get("EBCorrection");
    const GBRForest * forest_eb_uncertainty = (GBRForest *)weight_file->Get("EBUncertainty");
    const GBRForest * forest_ee_correction = (GBRForest *)weight_file->Get("EECorrection");
    const GBRForest * forest_ee_uncertainty = (GBRForest *)weight_file->Get("EEUncertainty");
    assert(forest_eb_correction);
    assert(forest_eb_uncertainty);
    assert(forest_ee_correction);
    assert(forest_ee_uncertainty);

    TTree *intree = 0;

    TChain *chainele = new TChain("Electrons");
    chainele->Add(applyingFile.c_str());
    chainele->LoadTree(0);    
    chainele->SetCacheSize(64*1024*1024);
    chainele->SetCacheLearnEntries();
    intree = chainele;

    GBRApply gbrapply;

    std::vector<string> *varseb = (std::vector<string> *)weight_file->Get("varlisteb");
    std::vector<string> *varsee = (std::vector<string> *)weight_file->Get("varlistee");

    TFile * outFile = new TFile((outFileDirectory + "ElectronEnergyRegression_targets_" + option + ".root").c_str(), "RECREATE");

    intree->LoadTree(0);
    TTree *intree_eb = gbrapply.ApplyAsFriend(intree, forest_eb_correction, *varseb, "targeteb");

    intree->LoadTree(0);
    TTree *intree_ebvar = gbrapply.ApplyAsFriend(intree, forest_eb_uncertainty, *varseb, "targetebvar");

    intree->LoadTree(0);
    TTree *intree_ee = gbrapply.ApplyAsFriend(intree, forest_ee_correction, *varsee, "targetee");

    intree->LoadTree(0);
    TTree *intree_eevar = gbrapply.ApplyAsFriend(intree, forest_ee_uncertainty, *varsee, "targeteevar");

    TCanvas* c1 = new TCanvas("c1", "targeteb", 1);
    intree_eb->Draw("targeteb");
    c1->SaveAs((outFileDirectory + "targeteb_" + option + ".gif").c_str());

    c1 = new TCanvas("c1", "targetee", 1);
    intree_ee->Draw("targetee");
    c1->SaveAs((outFileDirectory + "targetee_" + option + ".gif").c_str());

    // Writing trees
    intree_eb->SetName("targeteb_tree");
    intree_eb->Write();
    intree_ee->SetName("targetee_tree");
    intree_ee->Write();

    intree_ebvar->SetName("targetebvar_tree");
    intree_ebvar->Write();
    intree_eevar->SetName("targeteevar_tree");
    intree_eevar->Write();
  }

  // Now if the option is V10 or V11, we have to consider separately the two pT bins as well
  if (option == "V10" || option == "V11") {
    //  TFile * weight_file = new TFile("photon_2011_no_tk_vars_energy_reg_weights.root");
    TFile * weight_file = new TFile(weightFilename.c_str());
    const GBRForest * forest_eb_correction_lowPt = (GBRForest *)weight_file->Get("EBCorrection_lowPt");
    const GBRForest * forest_eb_correction_highPt = (GBRForest *)weight_file->Get("EBCorrection_highPt");
    const GBRForest * forest_eb_uncertainty_lowPt = (GBRForest *)weight_file->Get("EBUncertainty_lowPt");
    const GBRForest * forest_eb_uncertainty_highPt = (GBRForest *)weight_file->Get("EBUncertainty_highPt");
    const GBRForest * forest_ee_correction_lowPt = (GBRForest *)weight_file->Get("EECorrection_lowPt");
    const GBRForest * forest_ee_correction_highPt = (GBRForest *)weight_file->Get("EECorrection_highPt");
    const GBRForest * forest_ee_uncertainty_lowPt = (GBRForest *)weight_file->Get("EEUncertainty_lowPt");
    const GBRForest * forest_ee_uncertainty_highPt = (GBRForest *)weight_file->Get("EEUncertainty_highPt");

    assert(forest_eb_correction_lowPt);
    assert(forest_eb_correction_highPt);
    assert(forest_eb_uncertainty_lowPt);
    assert(forest_eb_uncertainty_highPt);
    assert(forest_ee_correction_lowPt);
    assert(forest_ee_correction_highPt);
    assert(forest_ee_uncertainty_lowPt);
    assert(forest_ee_uncertainty_highPt);


    TTree *intree = 0;

    TChain *chainele = new TChain("Electrons");
    chainele->Add(applyingFile.c_str());
    chainele->LoadTree(0);    
    chainele->SetCacheSize(64*1024*1024);
    chainele->SetCacheLearnEntries();
    intree = chainele;

    GBRApply gbrapply;

    std::vector<string> *varseb = (std::vector<string> *)weight_file->Get("varlisteb");
    std::vector<string> *varsee = (std::vector<string> *)weight_file->Get("varlistee");

    TFile * outFile = new TFile((outFileDirectory + "ElectronEnergyRegression_targets_" + option + ".root").c_str(), "RECREATE");

    intree->LoadTree(0);
    TTree *intree_eb_lowPt = gbrapply.ApplyAsFriend(intree, forest_eb_correction_lowPt, *varseb, "targeteb_lowPt");

    intree->LoadTree(0);
    TTree *intree_eb_highPt = gbrapply.ApplyAsFriend(intree, forest_eb_correction_highPt, *varseb, "targeteb_highPt");

    intree->LoadTree(0);
    TTree *intree_ebvar_lowPt = gbrapply.ApplyAsFriend(intree, forest_eb_uncertainty_lowPt, *varseb, "targetebvar_lowPt");

    intree->LoadTree(0);
    TTree *intree_ebvar_highPt = gbrapply.ApplyAsFriend(intree, forest_eb_uncertainty_highPt, *varseb, "targetebvar_highPt");

    intree->LoadTree(0);
    TTree *intree_ee_lowPt = gbrapply.ApplyAsFriend(intree, forest_ee_correction_lowPt, *varsee, "targetee_lowPt");

    intree->LoadTree(0);
    TTree *intree_ee_highPt = gbrapply.ApplyAsFriend(intree, forest_ee_correction_highPt, *varsee, "targetee_highPt");

    intree->LoadTree(0);
    TTree *intree_eevar_lowPt = gbrapply.ApplyAsFriend(intree, forest_ee_uncertainty_lowPt, *varsee, "targeteevar_lowPt");

    intree->LoadTree(0);
    TTree *intree_eevar_highPt = gbrapply.ApplyAsFriend(intree, forest_ee_uncertainty_highPt, *varsee, "targeteevar_highPt");

    TCanvas* c1 = new TCanvas("c1", "targeteb_lowPt", 1);
    intree_eb_lowPt->Draw("targeteb_lowPt");
    c1->SaveAs((outFileDirectory + "targeteb_lowPt_" + option + ".gif").c_str());

    c1 = new TCanvas("c1", "targeteb_highPt", 1);
    intree_eb_highPt->Draw("targeteb_highPt");
    c1->SaveAs((outFileDirectory + "targeteb_highPt_" + option + ".gif").c_str());

    c1 = new TCanvas("c1", "targetee_lowPt", 1);
    intree_ee_lowPt->Draw("targetee_lowPt");
    c1->SaveAs((outFileDirectory + "targetee_lowPt_" + option + ".gif").c_str());

    c1 = new TCanvas("c1", "targetee_highPt", 1);
    intree_ee_highPt->Draw("targetee_highPt");
    c1->SaveAs((outFileDirectory + "targetee_highPt_" + option + ".gif").c_str());

    // Not cloning the tree for now
    //TTree * intree_clone = intree->CloneTree();

    // Writing trees
    intree_eb_lowPt->SetName("targeteb_lowPt_tree");
    intree_eb_lowPt->Write();
    intree_eb_highPt->SetName("targeteb_highPt_tree");
    intree_eb_highPt->Write();
    intree_ee_lowPt->SetName("targetee_lowPt_tree");
    intree_ee_lowPt->Write();
    intree_ee_highPt->SetName("targetee_highPt_tree");
    intree_ee_highPt->Write();

    intree_ebvar_lowPt->SetName("targetebvar_lowPt_tree");
    intree_ebvar_lowPt->Write();
    intree_ebvar_highPt->SetName("targetebvar_highPt_tree");
    intree_ebvar_highPt->Write();
    intree_eevar_lowPt->SetName("targeteevar_lowPt_tree");
    intree_eevar_lowPt->Write();
    intree_eevar_highPt->SetName("targeteevar_highPt_tree");
    intree_eevar_highPt->Write();
  }
}
