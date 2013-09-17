#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TBranch.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <utility>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "scales2012.h"
#include "findAndReplace.h"

#include <RooAddPdf.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h> 
#include <RooWorkspace.h> 
#include <RooLandau.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h> 
#include <RooFFTConvPdf.h>

#include "CITHZZ/AnalysisStep/macro/FakeRateCalculator.h"
#include "CITHZZ/AnalysisStep/macro/CardTemplate.h"
#include "CITHZZ/AnalysisStep/macro/FitMaker.h"
#include "CITHZZ/AnalysisStep/macro/YieldMaker.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"

using namespace RooFit;

//gStyle->SetOptStat(0);			// Apparently doesn't work..

struct HiggsMassPointInfo {

  static float lumi;
  static float melacut;
  static float z1min;
  static float z2min;
  static float massLow;
  static float massHigh;
  static float massLowBkgFit;
  static float massHighBkgFit;
  static int   nBinsMass2D;
  static int   nBinsMELA2D;
  static bool  doSS;
  static bool  doShapeAnalysis;
  static bool  do1D;
  static std::string treeFolder;

  int mass;
  int massLowSigFit;
  int massHighSigFit;
  int gghid;
  int vbfid;

  HiggsMassPointInfo(int m, int mLow, int mHigh, int gid, int vid):
	  mass(m),
	  massLowSigFit(mLow),
	  massHighSigFit(mHigh),
	  gghid(gid),
	  vbfid(vid)
  {}
};

void analysisEngine(HiggsMassPointInfo hinfo, char* filename) {

  stringstream mass_str_ss;
  stringstream gghid_ss;
  stringstream vbfid_ss;

  mass_str_ss << hinfo.mass;
  gghid_ss << hinfo.gghid;
  vbfid_ss << hinfo.vbfid;

  std::string mass_str = mass_str_ss.str();
  std::string gghid = gghid_ss.str();
  std::string vbfid = vbfid_ss.str();

  // Initializing all necessary parameters
  std::cout << "Analyzing " << mass_str << " GeV mass point ... " << std::endl;
  std::string base_folder = hinfo.treeFolder;
  init();
  // Reading from different files for data and MC
  std::string data_folder = base_folder + "data/";
  std::string mc_folder = base_folder + "mc/";
  FakeRateCalculator FR(data_folder+"hzzTree.root", true, 40, 120, 0.0, 0.0, true);

  std::string ggh_rootfile = "hzzTree_id";
  ggh_rootfile += gghid;
  ggh_rootfile += ".root";

  std::string vbf_rootfile = "hzzTree_id";
  vbf_rootfile += vbfid;
  vbf_rootfile += ".root";

  DataYieldMaker ymaker_data;
  ZXYieldMaker   ymaker_zxss;
  ZZYieldMaker   ymaker_qqzz;
  ZZYieldMaker   ymaker_ggzz;
  ZZYieldMaker   ymaker_ghzz;
  ZZYieldMaker   ymaker_qhzz;

  // Initiating all weights (function in interface/scales.h or interface/scales2012.h)
  init();

  ymaker_data.fill(data_folder+"hzzTree.root");
  ymaker_zxss.fill(data_folder+"hzzTree.root"       , 1.0, FR, hinfo.doSS);
  ymaker_qqzz.fill(mc_folder+"hzzTree_id102.root" , xsecweights[102]        *hinfo.lumi, 0.0, false);
  ymaker_qqzz.fill(mc_folder+"hzzTree_id103.root" , xsecweights[103]        *hinfo.lumi, 0.0, false);
  ymaker_qqzz.fill(mc_folder+"hzzTree_id104.root" , xsecweights[104]        *hinfo.lumi, 0.0, false);
  ymaker_qqzz.fill(mc_folder+"hzzTree_id105.root" , xsecweights[105]        *hinfo.lumi, 0.0, false);
  ymaker_qqzz.fill(mc_folder+"hzzTree_id106.root" , xsecweights[106]        *hinfo.lumi, 0.0, false);
  ymaker_qqzz.fill(mc_folder+"hzzTree_id107.root" , xsecweights[107]        *hinfo.lumi, 0.0, false);
  ymaker_ggzz.fill(mc_folder+"hzzTree_id101.root" , xsecweights[101]        *hinfo.lumi, 0.0, false);
  ymaker_ggzz.fill(mc_folder+"hzzTree_id100.root" , xsecweights[100]        *hinfo.lumi, 0.0, false);
  ymaker_ghzz.fill(mc_folder+ggh_rootfile         , xsecweights[hinfo.gghid]*hinfo.lumi, 0.0, true );
  ymaker_qhzz.fill(mc_folder+vbf_rootfile         , xsecweights[hinfo.vbfid]*hinfo.lumi, 0.0, true );

  // File for the histograms
  TFile* outFile = new TFile(filename, "recreate");

  // Plotting the desired histograms
  TH1F data_hist_mass("data_hist_mass", "Mass, all, data", 20, 100, 160);
  TH1F data_hist_mass_4m("data_hist_mass_4m", "Mass, 4m, data", 20, 100, 160);
  TH1F data_hist_mass_4e("data_hist_mass_4e", "Mass, 4e, data", 20, 100, 160);
  TH1F data_hist_mass_2e2m("data_hist_mass_2e2m", "Mass, 2e2m, data", 20, 100, 160);
  TH1F data_hist_z1mass("data_hist_z1mass", "Z1 mass, all, data", 50, 0, 150);
  TH1F data_hist_z1mass_2m("data_hist_z1mass_2m", "Z1 mass, 2m, data", 50, 0, 150);
  TH1F data_hist_z1mass_2e("data_hist_z1mass_2e", "Z1 mass, 2e, data", 50, 0, 150);
  TH1F data_hist_z2mass("data_hist_z2mass", "Z2 mass, all, data", 50, 0, 150);
  TH1F data_hist_z2mass_2m("data_hist_z2mass_2m", "Z2 mass, 2m, data", 50, 0, 150);
  TH1F data_hist_z2mass_2e("data_hist_z2mass_2e", "Z2 mass, 2e, data", 50, 0, 150);
  TH1F data_hist_mela("data_hist_mela", "MELA, all, data", 50, 0, 1);
  TH1F data_hist_mela_4m("data_hist_mela_4m", "MELA, 4m, data", 50, 0, 1);
  TH1F data_hist_mela_4e("data_hist_mela_4e", "MELA, 4e, data", 50, 0, 1);
  TH1F data_hist_mela_2e2m("data_hist_mela_2e2m", "MELA, 2e2m, data", 50, 0, 1);
  TH2F data_hist_2D("data_hist_2D", "Z1 vs. Z2 masses, all, data", 50, 0, 150, 50, 0, 150);
  TH2F data_hist_2D_4m("data_hist_2D_4m", "Z1 vs. Z2 masses, 4m, data", 50, 0, 150, 50, 0, 150);
  TH2F data_hist_2D_4e("data_hist_2D_4e", "Z1 vs. Z2 masses, 4e, data", 50, 0, 150, 50, 0, 150);
  TH2F data_hist_2D_2e2m("data_hist_2D_2e2m", "Z1 vs. Z2 masses, 2e2m, data", 50, 0, 150, 50, 0, 150);
  TH2F data_hist_2D_2m2e("data_hist_2D_2m2e", "Z1 vs. Z2 masses, 2m2e, data", 50, 0, 150, 50, 0, 150);

  TH1F bkg_zxss_hist_mass("bkg_zxss_hist_mass", "Mass, all, ZX background", 20, 100, 160);
  TH1F bkg_zxss_hist_mass_4m("bkg_zxss_hist_mass_4m", "Mass, 4m, ZX background", 20, 100, 160);
  TH1F bkg_zxss_hist_mass_4e("bkg_zxss_hist_mass_4e", "Mass, 4e, ZX background", 20, 100, 160);
  TH1F bkg_zxss_hist_mass_2e2m("bkg_zxss_hist_mass_2e2m", "Mass, 2e2m, ZX background", 20, 100, 160);
  TH1F bkg_zxss_hist_z1mass("bkg_zxss_hist_z1mass", "Z1 mass, all, ZX background", 50, 0, 150);
  TH1F bkg_zxss_hist_z1mass_2m("bkg_zxss_hist_z1mass_2m", "Z1 mass, 2m, ZX background", 50, 0, 150);
  TH1F bkg_zxss_hist_z1mass_2e("bkg_zxss_hist_z1mass_2e", "Z1 mass, 2e, ZX background", 50, 0, 150);
  TH1F bkg_zxss_hist_z2mass("bkg_zxss_hist_z2mass", "Z2 mass, all, ZX background", 50, 0, 150);
  TH1F bkg_zxss_hist_z2mass_2m("bkg_zxss_hist_z2mass_2m", "Z2 mass, 2m, ZX background", 50, 0, 150);
  TH1F bkg_zxss_hist_z2mass_2e("bkg_zxss_hist_z2mass_2e", "Z2 mass, 2e, ZX background", 50, 0, 150);
  TH1F bkg_zxss_hist_mela("bkg_zxss_hist_mela", "MELA, all, ZX background", 50, 0, 1);
  TH1F bkg_zxss_hist_mela_4m("bkg_zxss_hist_mela_4m", "MELA, 4m, ZX background", 50, 0, 1);
  TH1F bkg_zxss_hist_mela_4e("bkg_zxss_hist_mela_4e", "MELA, 4e, ZX background", 50, 0, 1);
  TH1F bkg_zxss_hist_mela_2e2m("bkg_zxss_hist_mela_2e2m", "MELA, 2e2m, ZX background", 50, 0, 1);
  TH2F bkg_zxss_hist_2D("bkg_zxss_hist_2D", "Z1 vs. Z2 masses, all, ZX background", 50, 0, 150, 50, 0, 150);
  TH2F bkg_zxss_hist_2D_4m("bkg_zxss_hist_2D_4m", "Z1 vs. Z2 masses, 4m, ZX background", 50, 0, 150, 50, 0, 150);
  TH2F bkg_zxss_hist_2D_4e("bkg_zxss_hist_2D_4e", "Z1 vs. Z2 masses, 4e, ZX background", 50, 0, 150, 50, 0, 150);
  TH2F bkg_zxss_hist_2D_2e2m("bkg_zxss_hist_2D_2e2m", "Z1 vs. Z2 masses, 2e2m, ZX background", 50, 0, 150, 50, 0, 150);
  TH2F bkg_zxss_hist_2D_2m2e("bkg_zxss_hist_2D_2m2e", "Z1 vs. Z2 masses, 2m2e, ZX background", 50, 0, 150, 50, 0, 150);

  TH1F bkg_qqzz_hist_mass("bkg_qqzz_hist_mass", "Mass, all, qqZZ background", 20, 100, 160);
  TH1F bkg_qqzz_hist_mass_4m("bkg_qqzz_hist_mass_4m", "Mass, 4m, qqZZ background", 20, 100, 160);
  TH1F bkg_qqzz_hist_mass_4e("bkg_qqzz_hist_mass_4e", "Mass, 4e, qqZZ background", 20, 100, 160);
  TH1F bkg_qqzz_hist_mass_2e2m("bkg_qqzz_hist_mass_2e2m", "Mass, 2e2m, qqZZ background", 20, 100, 160);
  TH1F bkg_qqzz_hist_z1mass("bkg_qqzz_hist_z1mass", "Z1 mass, all, qqZZ background", 50, 0, 150);
  TH1F bkg_qqzz_hist_z1mass_2m("bkg_qqzz_hist_z1mass_2m", "Z1 mass, 2m, qqZZ background", 50, 0, 150);
  TH1F bkg_qqzz_hist_z1mass_2e("bkg_qqzz_hist_z1mass_2e", "Z1 mass, 2e, qqZZ background", 50, 0, 150);
  TH1F bkg_qqzz_hist_z2mass("bkg_qqzz_hist_z2mass", "Z2 mass, all, qqZZ background", 50, 0, 150);
  TH1F bkg_qqzz_hist_z2mass_2m("bkg_qqzz_hist_z2mass_2m", "Z2 mass, 2m, qqZZ background", 50, 0, 150);
  TH1F bkg_qqzz_hist_z2mass_2e("bkg_qqzz_hist_z2mass_2e", "Z2 mass, 2e, qqZZ background", 50, 0, 150);
  TH1F bkg_qqzz_hist_mela("bkg_qqzz_hist_mela", "MELA, all, qqZZ background", 50, 0, 1);
  TH1F bkg_qqzz_hist_mela_4m("bkg_qqzz_hist_mela_4m", "MELA, 4m, qqZZ background", 50, 0, 1);
  TH1F bkg_qqzz_hist_mela_4e("bkg_qqzz_hist_mela_4e", "MELA, 4e, qqZZ background", 50, 0, 1);
  TH1F bkg_qqzz_hist_mela_2e2m("bkg_qqzz_hist_mela_2e2m", "MELA, 2e2m, qqZZ background", 50, 0, 1);
  TH2F bkg_qqzz_hist_2D("bkg_qqzz_hist_2D", "Z1 vs. Z2 masses, all, qqZZ background", 50, 0, 150, 50, 0, 150);
  TH2F bkg_qqzz_hist_2D_4m("bkg_qqzz_hist_2D_4m", "Z1 vs. Z2 masses, 4m, qqZZ background", 50, 0, 150, 50, 0, 150);
  TH2F bkg_qqzz_hist_2D_4e("bkg_qqzz_hist_2D_4e", "Z1 vs. Z2 masses, 4e, qqZZ background", 50, 0, 150, 50, 0, 150);
  TH2F bkg_qqzz_hist_2D_2e2m("bkg_qqzz_hist_2D_2e2m", "Z1 vs. Z2 masses, 2e2m, qqZZ background", 50, 0, 150, 50, 0, 150);
  TH2F bkg_qqzz_hist_2D_2m2e("bkg_qqzz_hist_2D_2m2e", "Z1 vs. Z2 masses, 2m2e, qqZZ background", 50, 0, 150, 50, 0, 150);

  TH1F bkg_ggzz_hist_mass("bkg_ggzz_hist_mass", "Mass, all, ggZZ background", 20, 100, 160);
  TH1F bkg_ggzz_hist_mass_4m("bkg_ggzz_hist_mass_4m", "Mass, 4m, ggZZ background", 20, 100, 160);
  TH1F bkg_ggzz_hist_mass_4e("bkg_ggzz_hist_mass_4e", "Mass, 4e, ggZZ background", 20, 100, 160);
  TH1F bkg_ggzz_hist_mass_2e2m("bkg_ggzz_hist_mass_2e2m", "Mass, 2e2m, ggZZ background", 20, 100, 160);
  TH1F bkg_ggzz_hist_z1mass("bkg_ggzz_hist_z1mass", "Z1 mass, all, ggZZ background", 50, 0, 150);
  TH1F bkg_ggzz_hist_z1mass_2m("bkg_ggzz_hist_z1mass_2m", "Z1 mass, 2m, ggZZ background", 50, 0, 150);
  TH1F bkg_ggzz_hist_z1mass_2e("bkg_ggzz_hist_z1mass_2e", "Z1 mass, 2e, ggZZ background", 50, 0, 150);
  TH1F bkg_ggzz_hist_z2mass("bkg_ggzz_hist_z2mass", "Z2 mass, all, ggZZ background", 50, 0, 150);
  TH1F bkg_ggzz_hist_z2mass_2m("bkg_ggzz_hist_z2mass_2m", "Z2 mass, 2m, ggZZ background", 50, 0, 150);
  TH1F bkg_ggzz_hist_z2mass_2e("bkg_ggzz_hist_z2mass_2e", "Z2 mass, 2e, ggZZ background", 50, 0, 150);
  TH1F bkg_ggzz_hist_mela("bkg_ggzz_hist_mela", "MELA, all, ggZZ background", 50, 0, 1);
  TH1F bkg_ggzz_hist_mela_4m("bkg_ggzz_hist_mela_4m", "MELA, 4m, ggZZ background", 50, 0, 1);
  TH1F bkg_ggzz_hist_mela_4e("bkg_ggzz_hist_mela_4e", "MELA, 4e, ggZZ background", 50, 0, 1);
  TH1F bkg_ggzz_hist_mela_2e2m("bkg_ggzz_hist_mela_2e2m", "MELA, 2e2m, ggZZ background", 50, 0, 1);
  TH2F bkg_ggzz_hist_2D("bkg_ggzz_hist_2D", "Z1 vs. Z2 masses, all, ggZZ background", 50, 0, 150, 50, 0, 150);
  TH2F bkg_ggzz_hist_2D_4m("bkg_ggzz_hist_2D_4m", "Z1 vs. Z2 masses, 4m, ggZZ background", 50, 0, 150, 50, 0, 150);
  TH2F bkg_ggzz_hist_2D_4e("bkg_ggzz_hist_2D_4e", "Z1 vs. Z2 masses, 4e, ggZZ background", 50, 0, 150, 50, 0, 150);
  TH2F bkg_ggzz_hist_2D_2e2m("bkg_ggzz_hist_2D_2e2m", "Z1 vs. Z2 masses, 2e2m, ggZZ background", 50, 0, 150, 50, 0, 150);
  TH2F bkg_ggzz_hist_2D_2m2e("bkg_ggzz_hist_2D_2m2e", "Z1 vs. Z2 masses, 2m2e, ggZZ background", 50, 0, 150, 50, 0, 150);

  TH1F sig_ghzz_hist_mass("sig_ghzz_hist_mass", "Mass, all, gHZZ signal", 20, 100, 160);
  TH1F sig_ghzz_hist_mass_4m("sig_ghzz_hist_mass_4m", "Mass, 4m, gHZZ signal", 20, 100, 160);
  TH1F sig_ghzz_hist_mass_4e("sig_ghzz_hist_mass_4e", "Mass, 4e, gHZZ signal", 20, 100, 160);
  TH1F sig_ghzz_hist_mass_2e2m("sig_ghzz_hist_mass_2e2m", "Mass, 2e2m, gHZZ signal", 20, 100, 160);
  TH1F sig_ghzz_hist_z1mass("sig_ghzz_hist_z1mass", "Z1 mass, all, gHZZ signal", 50, 0, 150);
  TH1F sig_ghzz_hist_z1mass_2m("sig_ghzz_hist_z1mass_2m", "Z1 mass, 2m, gHZZ signal", 50, 0, 150);
  TH1F sig_ghzz_hist_z1mass_2e("sig_ghzz_hist_z1mass_2e", "Z1 mass, 2e, gHZZ signal", 50, 0, 150);
  TH1F sig_ghzz_hist_z2mass("sig_ghzz_hist_z2mass", "Z2 mass, all, gHZZ signal", 50, 0, 150);
  TH1F sig_ghzz_hist_z2mass_2m("sig_ghzz_hist_z2mass_2m", "Z2 mass, 2m, gHZZ signal", 50, 0, 150);
  TH1F sig_ghzz_hist_z2mass_2e("sig_ghzz_hist_z2mass_2e", "Z2 mass, 2e, gHZZ signal", 50, 0, 150);
  TH1F sig_ghzz_hist_mela("sig_ghzz_hist_mela", "MELA, all, gHZZ signal", 50, 0, 1);
  TH1F sig_ghzz_hist_mela_4m("sig_ghzz_hist_mela_4m", "MELA, 4m, gHZZ signal", 50, 0, 1);
  TH1F sig_ghzz_hist_mela_4e("sig_ghzz_hist_mela_4e", "MELA, 4e, gHZZ signal", 50, 0, 1);
  TH1F sig_ghzz_hist_mela_2e2m("sig_ghzz_hist_mela_2e2m", "MELA, 2e2m, gHZZ signal", 50, 0, 1);
  TH2F sig_ghzz_hist_2D("sig_ghzz_hist_2D", "Z1 vs. Z2 masses, all, gHZZ signal", 50, 0, 150, 50, 0, 150);
  TH2F sig_ghzz_hist_2D_4m("sig_ghzz_hist_2D_4m", "Z1 vs. Z2 masses, 4m, gHZZ signal", 50, 0, 150, 50, 0, 150);
  TH2F sig_ghzz_hist_2D_4e("sig_ghzz_hist_2D_4e", "Z1 vs. Z2 masses, 4e, gHZZ signal", 50, 0, 150, 50, 0, 150);
  TH2F sig_ghzz_hist_2D_2e2m("sig_ghzz_hist_2D_2e2m", "Z1 vs. Z2 masses, 2e2m, gHZZ signal", 50, 0, 150, 50, 0, 150);
  TH2F sig_ghzz_hist_2D_2m2e("sig_ghzz_hist_2D_2m2e", "Z1 vs. Z2 masses, 2m2e, gHZZ signal", 50, 0, 150, 50, 0, 150);

  TH1F sig_qhzz_hist_mass("sig_qhzz_hist_mass", "Mass, all, qHZZ signal", 20, 100, 160);
  TH1F sig_qhzz_hist_mass_4m("sig_qhzz_hist_mass_4m", "Mass, 4m, qHZZ signal", 20, 100, 160);
  TH1F sig_qhzz_hist_mass_4e("sig_qhzz_hist_mass_4e", "Mass, 4e, qHZZ signal", 20, 100, 160);
  TH1F sig_qhzz_hist_mass_2e2m("sig_qhzz_hist_mass_2e2m", "Mass, 2e2m, qHZZ signal", 20, 100, 160);
  TH1F sig_qhzz_hist_z1mass("sig_qhzz_hist_z1mass", "Z1 mass, all, qHZZ signal", 50, 0, 150);
  TH1F sig_qhzz_hist_z1mass_2m("sig_qhzz_hist_z1mass_2m", "Z1 mass, 2m, qHZZ signal", 50, 0, 150);
  TH1F sig_qhzz_hist_z1mass_2e("sig_qhzz_hist_z1mass_2e", "Z1 mass, 2e, qHZZ signal", 50, 0, 150);
  TH1F sig_qhzz_hist_z2mass("sig_qhzz_hist_z2mass", "Z2 mass, all, qHZZ signal", 50, 0, 150);
  TH1F sig_qhzz_hist_z2mass_2m("sig_qhzz_hist_z2mass_2m", "Z2 mass, 2m, qHZZ signal", 50, 0, 150);
  TH1F sig_qhzz_hist_z2mass_2e("sig_qhzz_hist_z2mass_2e", "Z2 mass, 2e, qHZZ signal", 50, 0, 150);
  TH1F sig_qhzz_hist_mela("sig_qhzz_hist_mela", "MELA, all, qHZZ signal", 50, 0, 1);
  TH1F sig_qhzz_hist_mela_4m("sig_qhzz_hist_mela_4m", "MELA, 4m, qHZZ signal", 50, 0, 1);
  TH1F sig_qhzz_hist_mela_4e("sig_qhzz_hist_mela_4e", "MELA, 4e, qHZZ signal", 50, 0, 1);
  TH1F sig_qhzz_hist_mela_2e2m("sig_qhzz_hist_mela_2e2m", "MELA, 2e2m, qHZZ signal", 50, 0, 1);
  TH2F sig_qhzz_hist_2D("sig_qhzz_hist_2D", "Z1 vs. Z2 masses, all, qHZZ signal", 50, 0, 150, 50, 0, 150);
  TH2F sig_qhzz_hist_2D_4m("sig_qhzz_hist_2D_4m", "Z1 vs. Z2 masses, 4m, qHZZ signal", 50, 0, 150, 50, 0, 150);
  TH2F sig_qhzz_hist_2D_4e("sig_qhzz_hist_2D_4e", "Z1 vs. Z2 masses, 4e, qHZZ signal", 50, 0, 150, 50, 0, 150);
  TH2F sig_qhzz_hist_2D_2e2m("sig_qhzz_hist_2D_2e2m", "Z1 vs. Z2 masses, 2e2m, qHZZ signal", 50, 0, 150, 50, 0, 150);
  TH2F sig_qhzz_hist_2D_2m2e("sig_qhzz_hist_2D_2m2e", "Z1 vs. Z2 masses, 2m2e, qHZZ signal", 50, 0, 150, 50, 0, 150);

  // Using function defined in YieldMaker.h to fill histograms according to channels and cuts
  ymaker_data.get1DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mass);
  ymaker_data.get1DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mass_4m);
  ymaker_data.get1DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mass);
  ymaker_data.get1DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mass_4e);
  ymaker_data.get1DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mass);
  ymaker_data.get1DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mass_2e2m);
  ymaker_data.get1DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mass);
  ymaker_data.get1DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mass_2e2m);

  ymaker_zxss.get1DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mass);
  ymaker_zxss.get1DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mass_4m);
  ymaker_zxss.get1DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mass);
  ymaker_zxss.get1DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mass_4e);
  ymaker_zxss.get1DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mass);
  ymaker_zxss.get1DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mass_2e2m);
  ymaker_zxss.get1DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mass);
  ymaker_zxss.get1DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mass_2e2m);

  ymaker_qqzz.get1DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mass);
  ymaker_qqzz.get1DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mass_4m);
  ymaker_qqzz.get1DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mass);
  ymaker_qqzz.get1DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mass_4e);
  ymaker_qqzz.get1DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mass);
  ymaker_qqzz.get1DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mass_2e2m);
  ymaker_qqzz.get1DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mass);
  ymaker_qqzz.get1DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mass_2e2m);

  ymaker_ggzz.get1DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mass);
  ymaker_ggzz.get1DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mass_4m);
  ymaker_ggzz.get1DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mass);
  ymaker_ggzz.get1DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mass_4e);
  ymaker_ggzz.get1DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mass);
  ymaker_ggzz.get1DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mass_2e2m);
  ymaker_ggzz.get1DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mass);
  ymaker_ggzz.get1DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mass_2e2m);

  ymaker_ghzz.get1DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mass);
  ymaker_ghzz.get1DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mass_4m);
  ymaker_ghzz.get1DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mass);
  ymaker_ghzz.get1DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mass_4e);
  ymaker_ghzz.get1DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mass);
  ymaker_ghzz.get1DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mass_2e2m);
  ymaker_ghzz.get1DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mass);
  ymaker_ghzz.get1DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mass_2e2m);

  ymaker_qhzz.get1DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mass);
  ymaker_qhzz.get1DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mass_4m);
  ymaker_qhzz.get1DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mass);
  ymaker_qhzz.get1DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mass_4e);
  ymaker_qhzz.get1DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mass);
  ymaker_qhzz.get1DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mass_2e2m);
  ymaker_qhzz.get1DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mass);
  ymaker_qhzz.get1DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mass_2e2m);

  ymaker_data.get1DHist_z1mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z1mass);
  ymaker_data.get1DHist_z1mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z1mass_2m);
  ymaker_data.get1DHist_z1mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z1mass);
  ymaker_data.get1DHist_z1mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z1mass_2e);
  ymaker_data.get1DHist_z1mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z1mass);
  ymaker_data.get1DHist_z1mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z1mass_2e);
  ymaker_data.get1DHist_z1mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z1mass);
  ymaker_data.get1DHist_z1mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z1mass_2m);

  ymaker_zxss.get1DHist_z1mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z1mass);
  ymaker_zxss.get1DHist_z1mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z1mass_2m);
  ymaker_zxss.get1DHist_z1mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z1mass);
  ymaker_zxss.get1DHist_z1mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z1mass_2e);
  ymaker_zxss.get1DHist_z1mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z1mass);
  ymaker_zxss.get1DHist_z1mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z1mass_2e);
  ymaker_zxss.get1DHist_z1mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z1mass);
  ymaker_zxss.get1DHist_z1mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z1mass_2m);

  ymaker_qqzz.get1DHist_z1mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z1mass);
  ymaker_qqzz.get1DHist_z1mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z1mass_2m);
  ymaker_qqzz.get1DHist_z1mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z1mass);
  ymaker_qqzz.get1DHist_z1mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z1mass_2e);
  ymaker_qqzz.get1DHist_z1mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z1mass);
  ymaker_qqzz.get1DHist_z1mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z1mass_2e);
  ymaker_qqzz.get1DHist_z1mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z1mass);
  ymaker_qqzz.get1DHist_z1mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z1mass_2m);

  ymaker_ggzz.get1DHist_z1mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass);
  ymaker_ggzz.get1DHist_z1mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass_2m);
  ymaker_ggzz.get1DHist_z1mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass);
  ymaker_ggzz.get1DHist_z1mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass_2e);
  ymaker_ggzz.get1DHist_z1mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass);
  ymaker_ggzz.get1DHist_z1mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass_2e);
  ymaker_ggzz.get1DHist_z1mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass);
  ymaker_ggzz.get1DHist_z1mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass_2m);

  ymaker_ggzz.get1DHist_z1mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass);
  ymaker_ggzz.get1DHist_z1mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass_2m);
  ymaker_ggzz.get1DHist_z1mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass);
  ymaker_ggzz.get1DHist_z1mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass_2e);
  ymaker_ggzz.get1DHist_z1mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass);
  ymaker_ggzz.get1DHist_z1mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass_2e);
  ymaker_ggzz.get1DHist_z1mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass);
  ymaker_ggzz.get1DHist_z1mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z1mass_2m);

  ymaker_ghzz.get1DHist_z1mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z1mass);
  ymaker_ghzz.get1DHist_z1mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z1mass_2m);
  ymaker_ghzz.get1DHist_z1mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z1mass);
  ymaker_ghzz.get1DHist_z1mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z1mass_2e);
  ymaker_ghzz.get1DHist_z1mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z1mass);
  ymaker_ghzz.get1DHist_z1mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z1mass_2e);
  ymaker_ghzz.get1DHist_z1mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z1mass);
  ymaker_ghzz.get1DHist_z1mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z1mass_2m);

  ymaker_qhzz.get1DHist_z1mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z1mass);
  ymaker_qhzz.get1DHist_z1mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z1mass_2m);
  ymaker_qhzz.get1DHist_z1mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z1mass);
  ymaker_qhzz.get1DHist_z1mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z1mass_2e);
  ymaker_qhzz.get1DHist_z1mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z1mass);
  ymaker_qhzz.get1DHist_z1mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z1mass_2e);
  ymaker_qhzz.get1DHist_z1mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z1mass);
  ymaker_qhzz.get1DHist_z1mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z1mass_2m);

  ymaker_data.get1DHist_z2mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z2mass);
  ymaker_data.get1DHist_z2mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z2mass_2m);
  ymaker_data.get1DHist_z2mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z2mass);
  ymaker_data.get1DHist_z2mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z2mass_2e);
  ymaker_data.get1DHist_z2mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z2mass);
  ymaker_data.get1DHist_z2mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z2mass_2m);
  ymaker_data.get1DHist_z2mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z2mass);
  ymaker_data.get1DHist_z2mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_z2mass_2e);

  ymaker_zxss.get1DHist_z2mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z2mass);
  ymaker_zxss.get1DHist_z2mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z2mass_2m);
  ymaker_zxss.get1DHist_z2mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z2mass);
  ymaker_zxss.get1DHist_z2mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z2mass_2e);
  ymaker_zxss.get1DHist_z2mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z2mass);
  ymaker_zxss.get1DHist_z2mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z2mass_2m);
  ymaker_zxss.get1DHist_z2mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z2mass);
  ymaker_zxss.get1DHist_z2mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_z2mass_2e);

  ymaker_qqzz.get1DHist_z2mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z2mass);
  ymaker_qqzz.get1DHist_z2mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z2mass_2m);
  ymaker_qqzz.get1DHist_z2mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z2mass);
  ymaker_qqzz.get1DHist_z2mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z2mass_2e);
  ymaker_qqzz.get1DHist_z2mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z2mass);
  ymaker_qqzz.get1DHist_z2mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z2mass_2m);
  ymaker_qqzz.get1DHist_z2mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z2mass);
  ymaker_qqzz.get1DHist_z2mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_z2mass_2e);

  ymaker_ggzz.get1DHist_z2mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass);
  ymaker_ggzz.get1DHist_z2mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass_2m);
  ymaker_ggzz.get1DHist_z2mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass);
  ymaker_ggzz.get1DHist_z2mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass_2e);
  ymaker_ggzz.get1DHist_z2mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass);
  ymaker_ggzz.get1DHist_z2mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass_2m);
  ymaker_ggzz.get1DHist_z2mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass);
  ymaker_ggzz.get1DHist_z2mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass_2e);

  ymaker_ggzz.get1DHist_z2mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass);
  ymaker_ggzz.get1DHist_z2mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass_2m);
  ymaker_ggzz.get1DHist_z2mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass);
  ymaker_ggzz.get1DHist_z2mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass_2e);
  ymaker_ggzz.get1DHist_z2mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass);
  ymaker_ggzz.get1DHist_z2mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass_2m);
  ymaker_ggzz.get1DHist_z2mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass);
  ymaker_ggzz.get1DHist_z2mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_z2mass_2e);

  ymaker_ghzz.get1DHist_z2mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z2mass);
  ymaker_ghzz.get1DHist_z2mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z2mass_2m);
  ymaker_ghzz.get1DHist_z2mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z2mass);
  ymaker_ghzz.get1DHist_z2mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z2mass_2e);
  ymaker_ghzz.get1DHist_z2mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z2mass);
  ymaker_ghzz.get1DHist_z2mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z2mass_2m);
  ymaker_ghzz.get1DHist_z2mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z2mass);
  ymaker_ghzz.get1DHist_z2mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_z2mass_2e);

  ymaker_qhzz.get1DHist_z2mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z2mass);
  ymaker_qhzz.get1DHist_z2mass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z2mass_2m);
  ymaker_qhzz.get1DHist_z2mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z2mass);
  ymaker_qhzz.get1DHist_z2mass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z2mass_2e);
  ymaker_qhzz.get1DHist_z2mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z2mass);
  ymaker_qhzz.get1DHist_z2mass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z2mass_2m);
  ymaker_qhzz.get1DHist_z2mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z2mass);
  ymaker_qhzz.get1DHist_z2mass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_z2mass_2e);

  ymaker_data.get1DHist_mela(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mela);
  ymaker_data.get1DHist_mela(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mela_4m);
  ymaker_data.get1DHist_mela(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mela);
  ymaker_data.get1DHist_mela(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mela_4e);
  ymaker_data.get1DHist_mela(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mela);
  ymaker_data.get1DHist_mela(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mela_2e2m);
  ymaker_data.get1DHist_mela(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mela);
  ymaker_data.get1DHist_mela(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_mela_2e2m);

  ymaker_zxss.get1DHist_mela(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mela);
  ymaker_zxss.get1DHist_mela(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mela_4m);
  ymaker_zxss.get1DHist_mela(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mela);
  ymaker_zxss.get1DHist_mela(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mela_4e);
  ymaker_zxss.get1DHist_mela(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mela);
  ymaker_zxss.get1DHist_mela(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mela_2e2m);
  ymaker_zxss.get1DHist_mela(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mela);
  ymaker_zxss.get1DHist_mela(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_mela_2e2m);

  ymaker_qqzz.get1DHist_mela(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mela);
  ymaker_qqzz.get1DHist_mela(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mela_4m);
  ymaker_qqzz.get1DHist_mela(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mela);
  ymaker_qqzz.get1DHist_mela(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mela_4e);
  ymaker_qqzz.get1DHist_mela(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mela);
  ymaker_qqzz.get1DHist_mela(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mela_2e2m);
  ymaker_qqzz.get1DHist_mela(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mela);
  ymaker_qqzz.get1DHist_mela(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_mela_2e2m);

  ymaker_ggzz.get1DHist_mela(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mela);
  ymaker_ggzz.get1DHist_mela(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mela_4m);
  ymaker_ggzz.get1DHist_mela(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mela);
  ymaker_ggzz.get1DHist_mela(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mela_4e);
  ymaker_ggzz.get1DHist_mela(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mela);
  ymaker_ggzz.get1DHist_mela(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mela_2e2m);
  ymaker_ggzz.get1DHist_mela(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mela);
  ymaker_ggzz.get1DHist_mela(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_mela_2e2m);

  ymaker_ghzz.get1DHist_mela(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mela);
  ymaker_ghzz.get1DHist_mela(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mela_4m);
  ymaker_ghzz.get1DHist_mela(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mela);
  ymaker_ghzz.get1DHist_mela(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mela_4e);
  ymaker_ghzz.get1DHist_mela(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mela);
  ymaker_ghzz.get1DHist_mela(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mela_2e2m);
  ymaker_ghzz.get1DHist_mela(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mela);
  ymaker_ghzz.get1DHist_mela(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_mela_2e2m);

  ymaker_qhzz.get1DHist_mela(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mela);
  ymaker_qhzz.get1DHist_mela(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mela_4m);
  ymaker_qhzz.get1DHist_mela(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mela);
  ymaker_qhzz.get1DHist_mela(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mela_4e);
  ymaker_qhzz.get1DHist_mela(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mela);
  ymaker_qhzz.get1DHist_mela(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mela_2e2m);
  ymaker_qhzz.get1DHist_mela(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mela);
  ymaker_qhzz.get1DHist_mela(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_mela_2e2m);

  // Now getting the 2D histograms
  ymaker_data.get2DHist_zmass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_2D);
  ymaker_data.get2DHist_zmass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_2D_4m);
  ymaker_data.get2DHist_zmass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_2D);
  ymaker_data.get2DHist_zmass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_2D_4e);
  ymaker_data.get2DHist_zmass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_2D);
  ymaker_data.get2DHist_zmass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_2D_2e2m);
  ymaker_data.get2DHist_zmass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_2D);
  ymaker_data.get2DHist_zmass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &data_hist_2D_2m2e);

  ymaker_zxss.get2DHist_zmass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_2D);
  ymaker_zxss.get2DHist_zmass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_2D_4m);
  ymaker_zxss.get2DHist_zmass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_2D);
  ymaker_zxss.get2DHist_zmass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_2D_4e);
  ymaker_zxss.get2DHist_zmass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_2D);
  ymaker_zxss.get2DHist_zmass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_2D_2e2m);
  ymaker_zxss.get2DHist_zmass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_2D);
  ymaker_zxss.get2DHist_zmass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_zxss_hist_2D_2m2e);

  ymaker_qqzz.get2DHist_zmass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_2D);
  ymaker_qqzz.get2DHist_zmass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_2D_4m);
  ymaker_qqzz.get2DHist_zmass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_2D);
  ymaker_qqzz.get2DHist_zmass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_2D_4e);
  ymaker_qqzz.get2DHist_zmass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_2D);
  ymaker_qqzz.get2DHist_zmass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_2D_2e2m);
  ymaker_qqzz.get2DHist_zmass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_2D);
  ymaker_qqzz.get2DHist_zmass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_qqzz_hist_2D_2m2e);

  ymaker_ggzz.get2DHist_zmass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D);
  ymaker_ggzz.get2DHist_zmass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D_4m);
  ymaker_ggzz.get2DHist_zmass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D);
  ymaker_ggzz.get2DHist_zmass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D_4e);
  ymaker_ggzz.get2DHist_zmass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D);
  ymaker_ggzz.get2DHist_zmass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D_2e2m);
  ymaker_ggzz.get2DHist_zmass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D);
  ymaker_ggzz.get2DHist_zmass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D_2m2e);

  ymaker_ggzz.get2DHist_zmass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D);
  ymaker_ggzz.get2DHist_zmass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D_4m);
  ymaker_ggzz.get2DHist_zmass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D);
  ymaker_ggzz.get2DHist_zmass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D_4e);
  ymaker_ggzz.get2DHist_zmass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D);
  ymaker_ggzz.get2DHist_zmass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D_2e2m);
  ymaker_ggzz.get2DHist_zmass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D);
  ymaker_ggzz.get2DHist_zmass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &bkg_ggzz_hist_2D_2m2e);

  ymaker_ghzz.get2DHist_zmass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_2D);
  ymaker_ghzz.get2DHist_zmass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_2D_4m);
  ymaker_ghzz.get2DHist_zmass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_2D);
  ymaker_ghzz.get2DHist_zmass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_2D_4e);
  ymaker_ghzz.get2DHist_zmass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_2D);
  ymaker_ghzz.get2DHist_zmass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_2D_2e2m);
  ymaker_ghzz.get2DHist_zmass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_2D);
  ymaker_ghzz.get2DHist_zmass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_ghzz_hist_2D_2m2e);

  ymaker_qhzz.get2DHist_zmass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_2D);
  ymaker_qhzz.get2DHist_zmass(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_2D_4m);
  ymaker_qhzz.get2DHist_zmass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_2D);
  ymaker_qhzz.get2DHist_zmass(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_2D_4e);
  ymaker_qhzz.get2DHist_zmass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_2D);
  ymaker_qhzz.get2DHist_zmass(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_2D_2e2m);
  ymaker_qhzz.get2DHist_zmass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_2D);
  ymaker_qhzz.get2DHist_zmass(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &sig_qhzz_hist_2D_2m2e);

  // Saving all histograms to file
  outFile->Write();			// Also, who thought this was a good idea??
}

float       HiggsMassPointInfo::lumi = 5.26;		// Using 2012 value
float       HiggsMassPointInfo::z1min = 40.;
float       HiggsMassPointInfo::z2min = 12.;
float       HiggsMassPointInfo::massLow = 118.;
float       HiggsMassPointInfo::massHigh = 135.;
float       HiggsMassPointInfo::massLowBkgFit = 100.;
float       HiggsMassPointInfo::massHighBkgFit = 600.;
float       HiggsMassPointInfo::melacut = -1.0;
int         HiggsMassPointInfo::nBinsMass2D = 100;
int         HiggsMassPointInfo::nBinsMELA2D = 30;
bool        HiggsMassPointInfo::doShapeAnalysis = true;
bool        HiggsMassPointInfo::do1D = true;
bool        HiggsMassPointInfo::doSS = true;
std::string HiggsMassPointInfo::treeFolder = "/mnt/hadoop/user/sixie/ntuples/hzz4l/step2/ichep2012/";

void doHZZHistograms_2012_cut(char* filename) {

  // HiggsMassPointInfo h126(mass, mass_low, mass_high, ggh_id, vbh_id)
  HiggsMassPointInfo h126(126, 100, 135, 1126, 2126);

  analysisEngine(h126, filename);
}
