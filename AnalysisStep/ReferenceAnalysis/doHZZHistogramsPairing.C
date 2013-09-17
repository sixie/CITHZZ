//=======================================================================================================
//					doHZZHistogramsPairing
// 
// Reads ntuples with the kinematic parameters of four leptons, pairs them according to various criteria 
// and plots the resulting Z1 and Z2 masses
//
// USAGE
//
// doHZZHistogramsPairing(char* directoryName, bool isCut = 0) {
// 
// INPUT
// directoryName		Name of the directory to which the plotted histograms are saved.
//				Directory must have the following subdirectories:
//				pairing0/  pairing1/  ... pairing7/    allpairs/
//
// isCut			Whether the program should look only into the signal region (118 - 135 GeV)
//
//
// OUTPUT
// 
// Saves all individual histograms to the directory "directoryName" and .gif superimposed plots in its
// subdirectories (specified above). The pairing are as follows, where (l1, l2) is assigned to Z1, and (l3, l4) is assigned to Z2:
//
//	0          keep ordering of leptons
//	1          (l1, l2) has smaller delta R than (l3, l4)
//	2          (l1, l2) has larger 3-momenta dot product than (l3, l4)
//	3          (l1, l2) has reconstructed mass closer to the Z pole than (l3, l4)
//	4          (l1, l2) has larger sum of absolute values of transverse momentum than (l3, l4)
//	5          pairing is such that Z1 and Z2 are most back-to-back, i.e. have the biggest delta phi
//	6          pairing is such that the mass difference between the two reconstructed bosons is reduced
//	7          pairing is such that, in the center of mass frame, the momenta of the Z's is minimized
//
// NOTE: DOES NOT CONSIDER ZX BACKGROUND
//________________________________________________________________________________________________________

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TBranch.h>
#include <iostream>
#include <iomanip>
#include <utility>
#include <string>
#include <sstream>
#include <iostream>
#include "FakeRateCalculator.h"
#include "YieldMaker.h"
#include "doHZZHistogramsPlot.C"
#include "doHZZHistogramsAllPairsPlot.C"

	struct MassPointInfoPairing {

	  float lumi;
	  float melacut;
	  float z1min, z1max;
	  float z2min, z2max;
	  float massLow, massHigh;
	  bool  do7TeV;
	  bool  doSS;
	  std::string treeFolder;

	  FakeRateCalculator FR;

	  DataYieldMaker ymaker_data;
	  ZXYieldMaker   ymaker_zxss;
	  ZZYieldMaker   ymaker_qqzz;
	  ZZYieldMaker   ymaker_ggzz;
	  ZZYieldMaker   ymaker_ghzz;
	  ZZYieldMaker   ymaker_qhzz;

	  void makeHistogramsPairing(float mass, int igghid, int ivbfid, std::string filename, int pairOption) {
	    // mass      = Higgs mass
	    // igghid    = ID for the gHZZ signal MC
	    // ivbfid    = ID for the qHZZ signal MC

	    // Prepare strings to read files
	    stringstream mass_str_ss;
	    stringstream gghid_ss;
	    stringstream vbfid_ss;
	    mass_str_ss << mass;
	    gghid_ss << igghid;
	    vbfid_ss << ivbfid;
	    std::string mass_str = mass_str_ss.str();
	    std::string gghid = gghid_ss.str();
	    std::string vbfid = vbfid_ss.str();

	    //  std::cout << "Analyzing " << mass_str << " GeV mass point ... " << std::endl;

	    std::string ggh_rootfile = "hzzTree_id";
	    ggh_rootfile += gghid;
	    ggh_rootfile += ".root";

	    std::string vbf_rootfile = "hzzTree_id";
	    vbf_rootfile += vbfid;
	    vbf_rootfile += ".root";

	    // Opening file to hold all the histograms
	    TFile outFile(filename.c_str(), "recreate");

	    // Creating the desired histograms
	    TH1F data_hist_mass("data_hist_mass", "Mass, all, data", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F data_hist_mass_4m("data_hist_mass_4m", "Mass, 4m, data", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F data_hist_mass_4e("data_hist_mass_4e", "Mass, 4e, data", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F data_hist_mass_2e2m("data_hist_mass_2e2m", "Mass, 2e2m, data", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F data_hist_z1mass("data_hist_z1mass", "Z1 mass, all, data", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F data_hist_z1mass_2m("data_hist_z1mass_2m", "Z1 mass, 2m, data", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F data_hist_z1mass_2e("data_hist_z1mass_2e", "Z1 mass, 2e, data", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F data_hist_z2mass("data_hist_z2mass", "Z2 mass, all, data", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F data_hist_z2mass_2m("data_hist_z2mass_2m", "Z2 mass, 2m, data", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F data_hist_z2mass_2e("data_hist_z2mass_2e", "Z2 mass, 2e, data", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F data_hist_mela("data_hist_mela", "MELA, all, data", 50, 0, 1);
	    TH1F data_hist_mela_4m("data_hist_mela_4m", "MELA, 4m, data", 50, 0, 1);
	    TH1F data_hist_mela_4e("data_hist_mela_4e", "MELA, 4e, data", 50, 0, 1);
	    TH1F data_hist_mela_2e2m("data_hist_mela_2e2m", "MELA, 2e2m, data", 50, 0, 1);
	    TH2F data_hist_2D("data_hist_2D", "Z1 vs. Z2 masses, all, data", (int) ( (z2max- z2min)/3 ), z2min, z2max, 50, 0, 150);
	    TH2F data_hist_2D_4m("data_hist_2D_4m", "Z1 vs. Z2 masses, 4m, data", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F data_hist_2D_4e("data_hist_2D_4e", "Z1 vs. Z2 masses, 4e, data", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F data_hist_2D_2e2m("data_hist_2D_2e2m", "Z1 vs. Z2 masses, 2e2m, data", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F data_hist_2D_2m2e("data_hist_2D_2m2e", "Z1 vs. Z2 masses, 2m2e, data", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);

	    TH1F bkg_zxss_hist_mass("bkg_zxss_hist_mass", "Mass, all, ZX background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_zxss_hist_mass_4m("bkg_zxss_hist_mass_4m", "Mass, 4m, ZX background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_zxss_hist_mass_4e("bkg_zxss_hist_mass_4e", "Mass, 4e, ZX background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_zxss_hist_mass_2e2m("bkg_zxss_hist_mass_2e2m", "Mass, 2e2m, ZX background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_zxss_hist_z1mass("bkg_zxss_hist_z1mass", "Z1 mass, all, ZX background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_zxss_hist_z1mass_2m("bkg_zxss_hist_z1mass_2m", "Z1 mass, 2m, ZX background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_zxss_hist_z1mass_2e("bkg_zxss_hist_z1mass_2e", "Z1 mass, 2e, ZX background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_zxss_hist_z2mass("bkg_zxss_hist_z2mass", "Z2 mass, all, ZX background", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F bkg_zxss_hist_z2mass_2m("bkg_zxss_hist_z2mass_2m", "Z2 mass, 2m, ZX background", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F bkg_zxss_hist_z2mass_2e("bkg_zxss_hist_z2mass_2e", "Z2 mass, 2e, ZX background", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F bkg_zxss_hist_mela("bkg_zxss_hist_mela", "MELA, all, ZX background", 50, 0, 1);
	    TH1F bkg_zxss_hist_mela_4m("bkg_zxss_hist_mela_4m", "MELA, 4m, ZX background", 50, 0, 1);
	    TH1F bkg_zxss_hist_mela_4e("bkg_zxss_hist_mela_4e", "MELA, 4e, ZX background", 50, 0, 1);
	    TH1F bkg_zxss_hist_mela_2e2m("bkg_zxss_hist_mela_2e2m", "MELA, 2e2m, ZX background", 50, 0, 1);
	    TH2F bkg_zxss_hist_2D("bkg_zxss_hist_2D", "Z1 vs. Z2 masses, all, ZX background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F bkg_zxss_hist_2D_4m("bkg_zxss_hist_2D_4m", "Z1 vs. Z2 masses, 4m, ZX background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F bkg_zxss_hist_2D_4e("bkg_zxss_hist_2D_4e", "Z1 vs. Z2 masses, 4e, ZX background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F bkg_zxss_hist_2D_2e2m("bkg_zxss_hist_2D_2e2m", "Z1 vs. Z2 masses, 2e2m, ZX background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F bkg_zxss_hist_2D_2m2e("bkg_zxss_hist_2D_2m2e", "Z1 vs. Z2 masses, 2m2e, ZX background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);

	    TH1F bkg_qqzz_hist_mass("bkg_qqzz_hist_mass", "Mass, all, qqZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_qqzz_hist_mass_4m("bkg_qqzz_hist_mass_4m", "Mass, 4m, qqZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_qqzz_hist_mass_4e("bkg_qqzz_hist_mass_4e", "Mass, 4e, qqZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_qqzz_hist_mass_2e2m("bkg_qqzz_hist_mass_2e2m", "Mass, 2e2m, qqZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_qqzz_hist_z1mass("bkg_qqzz_hist_z1mass", "Z1 mass, all, qqZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_qqzz_hist_z1mass_2m("bkg_qqzz_hist_z1mass_2m", "Z1 mass, 2m, qqZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_qqzz_hist_z1mass_2e("bkg_qqzz_hist_z1mass_2e", "Z1 mass, 2e, qqZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_qqzz_hist_z2mass("bkg_qqzz_hist_z2mass", "Z2 mass, all, qqZZ background", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F bkg_qqzz_hist_z2mass_2m("bkg_qqzz_hist_z2mass_2m", "Z2 mass, 2m, qqZZ background", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F bkg_qqzz_hist_z2mass_2e("bkg_qqzz_hist_z2mass_2e", "Z2 mass, 2e, qqZZ background", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F bkg_qqzz_hist_mela("bkg_qqzz_hist_mela", "MELA, all, qqZZ background", 50, 0, 1);
	    TH1F bkg_qqzz_hist_mela_4m("bkg_qqzz_hist_mela_4m", "MELA, 4m, qqZZ background", 50, 0, 1);
	    TH1F bkg_qqzz_hist_mela_4e("bkg_qqzz_hist_mela_4e", "MELA, 4e, qqZZ background", 50, 0, 1);
	    TH1F bkg_qqzz_hist_mela_2e2m("bkg_qqzz_hist_mela_2e2m", "MELA, 2e2m, qqZZ background", 50, 0, 1);
	    TH2F bkg_qqzz_hist_2D("bkg_qqzz_hist_2D", "Z1 vs. Z2 masses, all, qqZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F bkg_qqzz_hist_2D_4m("bkg_qqzz_hist_2D_4m", "Z1 vs. Z2 masses, 4m, qqZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F bkg_qqzz_hist_2D_4e("bkg_qqzz_hist_2D_4e", "Z1 vs. Z2 masses, 4e, qqZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F bkg_qqzz_hist_2D_2e2m("bkg_qqzz_hist_2D_2e2m", "Z1 vs. Z2 masses, 2e2m, qqZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F bkg_qqzz_hist_2D_2m2e("bkg_qqzz_hist_2D_2m2e", "Z1 vs. Z2 masses, 2m2e, qqZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);

	    TH1F bkg_ggzz_hist_mass("bkg_ggzz_hist_mass", "Mass, all, ggZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_ggzz_hist_mass_4m("bkg_ggzz_hist_mass_4m", "Mass, 4m, ggZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_ggzz_hist_mass_4e("bkg_ggzz_hist_mass_4e", "Mass, 4e, ggZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_ggzz_hist_mass_2e2m("bkg_ggzz_hist_mass_2e2m", "Mass, 2e2m, ggZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_ggzz_hist_z1mass("bkg_ggzz_hist_z1mass", "Z1 mass, all, ggZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_ggzz_hist_z1mass_2m("bkg_ggzz_hist_z1mass_2m", "Z1 mass, 2m, ggZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_ggzz_hist_z1mass_2e("bkg_ggzz_hist_z1mass_2e", "Z1 mass, 2e, ggZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_ggzz_hist_z2mass("bkg_ggzz_hist_z2mass", "Z2 mass, all, ggZZ background", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F bkg_ggzz_hist_z2mass_2m("bkg_ggzz_hist_z2mass_2m", "Z2 mass, 2m, ggZZ background", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F bkg_ggzz_hist_z2mass_2e("bkg_ggzz_hist_z2mass_2e", "Z2 mass, 2e, ggZZ background", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F bkg_ggzz_hist_mela("bkg_ggzz_hist_mela", "MELA, all, ggZZ background", 50, 0, 1);
	    TH1F bkg_ggzz_hist_mela_4m("bkg_ggzz_hist_mela_4m", "MELA, 4m, ggZZ background", 50, 0, 1);
	    TH1F bkg_ggzz_hist_mela_4e("bkg_ggzz_hist_mela_4e", "MELA, 4e, ggZZ background", 50, 0, 1);
	    TH1F bkg_ggzz_hist_mela_2e2m("bkg_ggzz_hist_mela_2e2m", "MELA, 2e2m, ggZZ background", 50, 0, 1);
	    TH2F bkg_ggzz_hist_2D("bkg_ggzz_hist_2D", "Z1 vs. Z2 masses, all, ggZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F bkg_ggzz_hist_2D_4m("bkg_ggzz_hist_2D_4m", "Z1 vs. Z2 masses, 4m, ggZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F bkg_ggzz_hist_2D_4e("bkg_ggzz_hist_2D_4e", "Z1 vs. Z2 masses, 4e, ggZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F bkg_ggzz_hist_2D_2e2m("bkg_ggzz_hist_2D_2e2m", "Z1 vs. Z2 masses, 2e2m, ggZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F bkg_ggzz_hist_2D_2m2e("bkg_ggzz_hist_2D_2m2e", "Z1 vs. Z2 masses, 2m2e, ggZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);

	    TH1F sig_ghzz_hist_mass("sig_ghzz_hist_mass", "Mass, all, gHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_ghzz_hist_mass_4m("sig_ghzz_hist_mass_4m", "Mass, 4m, gHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_ghzz_hist_mass_4e("sig_ghzz_hist_mass_4e", "Mass, 4e, gHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_ghzz_hist_mass_2e2m("sig_ghzz_hist_mass_2e2m", "Mass, 2e2m, gHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_ghzz_hist_z1mass("sig_ghzz_hist_z1mass", "Z1 mass, all, gHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F sig_ghzz_hist_z1mass_2m("sig_ghzz_hist_z1mass_2m", "Z1 mass, 2m, gHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F sig_ghzz_hist_z1mass_2e("sig_ghzz_hist_z1mass_2e", "Z1 mass, 2e, gHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F sig_ghzz_hist_z2mass("sig_ghzz_hist_z2mass", "Z2 mass, all, gHZZ signal", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F sig_ghzz_hist_z2mass_2m("sig_ghzz_hist_z2mass_2m", "Z2 mass, 2m, gHZZ signal", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F sig_ghzz_hist_z2mass_2e("sig_ghzz_hist_z2mass_2e", "Z2 mass, 2e, gHZZ signal", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F sig_ghzz_hist_mela("sig_ghzz_hist_mela", "MELA, all, gHZZ signal", 50, 0, 1);
	    TH1F sig_ghzz_hist_mela_4m("sig_ghzz_hist_mela_4m", "MELA, 4m, gHZZ signal", 50, 0, 1);
	    TH1F sig_ghzz_hist_mela_4e("sig_ghzz_hist_mela_4e", "MELA, 4e, gHZZ signal", 50, 0, 1);
	    TH1F sig_ghzz_hist_mela_2e2m("sig_ghzz_hist_mela_2e2m", "MELA, 2e2m, gHZZ signal", 50, 0, 1);
	    TH2F sig_ghzz_hist_2D("sig_ghzz_hist_2D", "Z1 vs. Z2 masses, all, gHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F sig_ghzz_hist_2D_4m("sig_ghzz_hist_2D_4m", "Z1 vs. Z2 masses, 4m, gHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F sig_ghzz_hist_2D_4e("sig_ghzz_hist_2D_4e", "Z1 vs. Z2 masses, 4e, gHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F sig_ghzz_hist_2D_2e2m("sig_ghzz_hist_2D_2e2m", "Z1 vs. Z2 masses, 2e2m, gHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F sig_ghzz_hist_2D_2m2e("sig_ghzz_hist_2D_2m2e", "Z1 vs. Z2 masses, 2m2e, gHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);

	    TH1F sig_qhzz_hist_mass("sig_qhzz_hist_mass", "Mass, all, qHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_qhzz_hist_mass_4m("sig_qhzz_hist_mass_4m", "Mass, 4m, qHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_qhzz_hist_mass_4e("sig_qhzz_hist_mass_4e", "Mass, 4e, qHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_qhzz_hist_mass_2e2m("sig_qhzz_hist_mass_2e2m", "Mass, 2e2m, qHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_qhzz_hist_z1mass("sig_qhzz_hist_z1mass", "Z1 mass, all, qHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F sig_qhzz_hist_z1mass_2m("sig_qhzz_hist_z1mass_2m", "Z1 mass, 2m, qHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F sig_qhzz_hist_z1mass_2e("sig_qhzz_hist_z1mass_2e", "Z1 mass, 2e, qHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F sig_qhzz_hist_z2mass("sig_qhzz_hist_z2mass", "Z2 mass, all, qHZZ signal", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F sig_qhzz_hist_z2mass_2m("sig_qhzz_hist_z2mass_2m", "Z2 mass, 2m, qHZZ signal", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F sig_qhzz_hist_z2mass_2e("sig_qhzz_hist_z2mass_2e", "Z2 mass, 2e, qHZZ signal", (int) ( (z2max- z2min)/3 ), z2min, z2max);
	    TH1F sig_qhzz_hist_mela("sig_qhzz_hist_mela", "MELA, all, qHZZ signal", 50, 0, 1);
	    TH1F sig_qhzz_hist_mela_4m("sig_qhzz_hist_mela_4m", "MELA, 4m, qHZZ signal", 50, 0, 1);
	    TH1F sig_qhzz_hist_mela_4e("sig_qhzz_hist_mela_4e", "MELA, 4e, qHZZ signal", 50, 0, 1);
	    TH1F sig_qhzz_hist_mela_2e2m("sig_qhzz_hist_mela_2e2m", "MELA, 2e2m, qHZZ signal", 50, 0, 1);
	    TH2F sig_qhzz_hist_2D("sig_qhzz_hist_2D", "Z1 vs. Z2 masses, all, qHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F sig_qhzz_hist_2D_4m("sig_qhzz_hist_2D_4m", "Z1 vs. Z2 masses, 4m, qHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F sig_qhzz_hist_2D_4e("sig_qhzz_hist_2D_4e", "Z1 vs. Z2 masses, 4e, qHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F sig_qhzz_hist_2D_2e2m("sig_qhzz_hist_2D_2e2m", "Z1 vs. Z2 masses, 2e2m, qHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);
	    TH2F sig_qhzz_hist_2D_2m2e("sig_qhzz_hist_2D_2m2e", "Z1 vs. Z2 masses, 2m2e, qHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max, (z2max- z2min), z2min, z2max);

	    // Using function defined in YieldMaker.h to fill histograms according to channels and cuts
	    ymaker_data.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass);
	    ymaker_data.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass_4m);
	    ymaker_data.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass);
	    ymaker_data.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass_4e);
	    ymaker_data.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass);
	    ymaker_data.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass_2e2m);
	    ymaker_data.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass);
	    ymaker_data.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass_2e2m);

	    ymaker_zxss.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass);
	    ymaker_zxss.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass_4m);
	    ymaker_zxss.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass);
	    ymaker_zxss.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass_4e);
	    ymaker_zxss.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass);
	    ymaker_zxss.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass_2e2m);
	    ymaker_zxss.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass);
	    ymaker_zxss.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass_2e2m);

	    ymaker_qqzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass);
	    ymaker_qqzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass_4m);
	    ymaker_qqzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass);
	    ymaker_qqzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass_4e);
	    ymaker_qqzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass);
	    ymaker_qqzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass_2e2m);
	    ymaker_qqzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass);
	    ymaker_qqzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass_2e2m);

	    ymaker_ggzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass);
	    ymaker_ggzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass_4m);
	    ymaker_ggzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass);
	    ymaker_ggzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass_4e);
	    ymaker_ggzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass);
	    ymaker_ggzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass_2e2m);
	    ymaker_ggzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass);
	    ymaker_ggzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass_2e2m);

	    ymaker_ghzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass);
	    ymaker_ghzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass_4m);
	    ymaker_ghzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass);
	    ymaker_ghzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass_4e);
	    ymaker_ghzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass);
	    ymaker_ghzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass_2e2m);
	    ymaker_ghzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass);
	    ymaker_ghzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass_2e2m);

	    ymaker_qhzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass);
	    ymaker_qhzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass_4m);
	    ymaker_qhzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass);
	    ymaker_qhzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass_4e);
	    ymaker_qhzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass);
	    ymaker_qhzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass_2e2m);
	    ymaker_qhzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass);
	    ymaker_qhzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass_2e2m);

	    ymaker_data.get1DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_z1mass, &data_hist_z2mass, pairOption);
	    ymaker_data.get1DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_z1mass_2m, &data_hist_z2mass_2m, pairOption);
	    ymaker_data.get1DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_z1mass, &data_hist_z2mass, pairOption);
	    ymaker_data.get1DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_z1mass_2e, &data_hist_z2mass_2e, pairOption);
	    ymaker_data.get1DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_z1mass, &data_hist_z2mass, pairOption);
	    ymaker_data.get1DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_z1mass_2e, &data_hist_z2mass_2e, pairOption);
	    ymaker_data.get1DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_z1mass, &data_hist_z2mass, pairOption);
	    ymaker_data.get1DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_z1mass_2m, &data_hist_z2mass_2m, pairOption);

	    ymaker_zxss.get1DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_z1mass, &bkg_zxss_hist_z2mass, pairOption);
	    ymaker_zxss.get1DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_z1mass_2m, &bkg_zxss_hist_z2mass_2m, pairOption);
	    ymaker_zxss.get1DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_z1mass, &bkg_zxss_hist_z2mass, pairOption);
	    ymaker_zxss.get1DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_z1mass_2e, &bkg_zxss_hist_z2mass_2e, pairOption);
	    ymaker_zxss.get1DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_z1mass, &bkg_zxss_hist_z2mass, pairOption);
	    ymaker_zxss.get1DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_z1mass_2e, &bkg_zxss_hist_z2mass_2e, pairOption);
	    ymaker_zxss.get1DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_z1mass, &bkg_zxss_hist_z2mass, pairOption);
	    ymaker_zxss.get1DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_z1mass_2m, &bkg_zxss_hist_z2mass_2m, pairOption);

	    ymaker_ggzz.get1DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_z1mass, &bkg_ggzz_hist_z2mass, pairOption);
	    ymaker_ggzz.get1DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_z1mass_2m, &bkg_ggzz_hist_z2mass_2m, pairOption);
	    ymaker_ggzz.get1DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_z1mass, &bkg_ggzz_hist_z2mass, pairOption);
	    ymaker_ggzz.get1DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_z1mass_2e, &bkg_ggzz_hist_z2mass_2e, pairOption);
	    ymaker_ggzz.get1DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_z1mass, &bkg_ggzz_hist_z2mass, pairOption);
	    ymaker_ggzz.get1DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_z1mass_2e, &bkg_ggzz_hist_z2mass_2e, pairOption);
	    ymaker_ggzz.get1DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_z1mass, &bkg_ggzz_hist_z2mass, pairOption);
	    ymaker_ggzz.get1DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_z1mass_2m, &bkg_ggzz_hist_z2mass_2m, pairOption);

	    ymaker_qqzz.get1DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_z1mass, &bkg_qqzz_hist_z2mass, pairOption);
	    ymaker_qqzz.get1DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_z1mass_2m, &bkg_qqzz_hist_z2mass_2m, pairOption);
	    ymaker_qqzz.get1DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_z1mass, &bkg_qqzz_hist_z2mass, pairOption);
	    ymaker_qqzz.get1DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_z1mass_2e, &bkg_qqzz_hist_z2mass_2e, pairOption);
	    ymaker_qqzz.get1DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_z1mass, &bkg_qqzz_hist_z2mass, pairOption);
	    ymaker_qqzz.get1DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_z1mass_2e, &bkg_qqzz_hist_z2mass_2e, pairOption);
	    ymaker_qqzz.get1DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_z1mass, &bkg_qqzz_hist_z2mass, pairOption);
	    ymaker_qqzz.get1DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_z1mass_2m, &bkg_qqzz_hist_z2mass_2m, pairOption);

	    ymaker_qhzz.get1DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_z1mass, &sig_qhzz_hist_z2mass, pairOption);
	    ymaker_qhzz.get1DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_z1mass_2m, &sig_qhzz_hist_z2mass_2m, pairOption);
	    ymaker_qhzz.get1DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_z1mass, &sig_qhzz_hist_z2mass, pairOption);
	    ymaker_qhzz.get1DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_z1mass_2e, &sig_qhzz_hist_z2mass_2e, pairOption);
	    ymaker_qhzz.get1DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_z1mass, &sig_qhzz_hist_z2mass, pairOption);
	    ymaker_qhzz.get1DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_z1mass_2e, &sig_qhzz_hist_z2mass_2e, pairOption);
	    ymaker_qhzz.get1DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_z1mass, &sig_qhzz_hist_z2mass, pairOption);
	    ymaker_qhzz.get1DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_z1mass_2m, &sig_qhzz_hist_z2mass_2m, pairOption);

	    ymaker_ghzz.get1DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_z1mass, &sig_ghzz_hist_z2mass, pairOption);
	    ymaker_ghzz.get1DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_z1mass_2m, &sig_ghzz_hist_z2mass_2m, pairOption);
	    ymaker_ghzz.get1DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_z1mass, &sig_ghzz_hist_z2mass, pairOption);
	    ymaker_ghzz.get1DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_z1mass_2e, &sig_ghzz_hist_z2mass_2e, pairOption);
	    ymaker_ghzz.get1DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_z1mass, &sig_ghzz_hist_z2mass, pairOption);
	    ymaker_ghzz.get1DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_z1mass_2e, &sig_ghzz_hist_z2mass_2e, pairOption);
	    ymaker_ghzz.get1DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_z1mass, &sig_ghzz_hist_z2mass, pairOption);
	    ymaker_ghzz.get1DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_z1mass_2m, &sig_ghzz_hist_z2mass_2m, pairOption);

	    ymaker_data.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela);
	    ymaker_data.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela_4m);
	    ymaker_data.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela);
	    ymaker_data.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela_4e);
	    ymaker_data.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela);
	    ymaker_data.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela_2e2m);
	    ymaker_data.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela);
	    ymaker_data.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela_2e2m);

	    ymaker_zxss.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela);
	    ymaker_zxss.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela_4m);
	    ymaker_zxss.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela);
	    ymaker_zxss.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela_4e);
	    ymaker_zxss.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela);
	    ymaker_zxss.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela_2e2m);
	    ymaker_zxss.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela);
	    ymaker_zxss.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela_2e2m);

	    ymaker_qqzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela);
	    ymaker_qqzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela_4m);
	    ymaker_qqzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela);
	    ymaker_qqzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela_4e);
	    ymaker_qqzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela);
	    ymaker_qqzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela_2e2m);
	    ymaker_qqzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela);
	    ymaker_qqzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela_2e2m);

	    ymaker_ggzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela);
	    ymaker_ggzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela_4m);
	    ymaker_ggzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela);
	    ymaker_ggzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela_4e);
	    ymaker_ggzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela);
	    ymaker_ggzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela_2e2m);
	    ymaker_ggzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela);
	    ymaker_ggzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela_2e2m);

	    ymaker_ghzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela);
	    ymaker_ghzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela_4m);
	    ymaker_ghzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela);
	    ymaker_ghzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela_4e);
	    ymaker_ghzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela);
	    ymaker_ghzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela_2e2m);
	    ymaker_ghzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela);
	    ymaker_ghzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela_2e2m);

	    ymaker_qhzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela);
	    ymaker_qhzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela_4m);
	    ymaker_qhzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela);
	    ymaker_qhzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela_4e);
	    ymaker_qhzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela);
	    ymaker_qhzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela_2e2m);
	    ymaker_qhzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela);
	    ymaker_qhzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela_2e2m);

	    // Now getting the 2D histograms
	    ymaker_data.get2DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_2D, pairOption);
	    ymaker_data.get2DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_2D_4m, pairOption);
	    ymaker_data.get2DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_2D, pairOption);
	    ymaker_data.get2DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_2D_4e, pairOption);
	    ymaker_data.get2DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_2D, pairOption);
	    ymaker_data.get2DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_2D_2e2m, pairOption);
	    ymaker_data.get2DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_2D, pairOption);
	    ymaker_data.get2DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_2D_2m2e, pairOption);

	    ymaker_zxss.get2DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_2D, pairOption);
	    ymaker_zxss.get2DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_2D_4m, pairOption);
	    ymaker_zxss.get2DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_2D, pairOption);
	    ymaker_zxss.get2DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_2D_4e, pairOption);
	    ymaker_zxss.get2DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_2D, pairOption);
	    ymaker_zxss.get2DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_2D_2e2m, pairOption);
	    ymaker_zxss.get2DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_2D, pairOption);
	    ymaker_zxss.get2DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_2D_2m2e, pairOption);

	    ymaker_qqzz.get2DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_2D, pairOption);
	    ymaker_qqzz.get2DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_2D_4m, pairOption);
	    ymaker_qqzz.get2DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_2D, pairOption);
	    ymaker_qqzz.get2DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_2D_4e, pairOption);
	    ymaker_qqzz.get2DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_2D, pairOption);
	    ymaker_qqzz.get2DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_2D_2e2m, pairOption);
	    ymaker_qqzz.get2DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_2D, pairOption);
	    ymaker_qqzz.get2DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_2D_2m2e, pairOption);

	    ymaker_ggzz.get2DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D, pairOption);
	    ymaker_ggzz.get2DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D_4m, pairOption);
	    ymaker_ggzz.get2DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D, pairOption);
	    ymaker_ggzz.get2DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D_4e, pairOption);
	    ymaker_ggzz.get2DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D, pairOption);
	    ymaker_ggzz.get2DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D_2e2m, pairOption);
	    ymaker_ggzz.get2DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D, pairOption);
	    ymaker_ggzz.get2DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D_2m2e, pairOption);

	    ymaker_ggzz.get2DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D, pairOption);
	    ymaker_ggzz.get2DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D_4m, pairOption);
	    ymaker_ggzz.get2DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D, pairOption);
	    ymaker_ggzz.get2DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D_4e, pairOption);
	    ymaker_ggzz.get2DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D, pairOption);
	    ymaker_ggzz.get2DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D_2e2m, pairOption);
	    ymaker_ggzz.get2DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D, pairOption);
	    ymaker_ggzz.get2DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_2D_2m2e, pairOption);

	    ymaker_ghzz.get2DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_2D, pairOption);
	    ymaker_ghzz.get2DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_2D_4m, pairOption);
	    ymaker_ghzz.get2DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_2D, pairOption);
	    ymaker_ghzz.get2DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_2D_4e, pairOption);
	    ymaker_ghzz.get2DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_2D, pairOption);
	    ymaker_ghzz.get2DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_2D_2e2m, pairOption);
	    ymaker_ghzz.get2DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_2D, pairOption);
	    ymaker_ghzz.get2DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_2D_2m2e, pairOption);

	    ymaker_qhzz.get2DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_2D, pairOption);
	    ymaker_qhzz.get2DHist_zmass_pairing(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_2D_4m, pairOption);
	    ymaker_qhzz.get2DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_2D, pairOption);
	    ymaker_qhzz.get2DHist_zmass_pairing(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_2D_4e, pairOption);
	    ymaker_qhzz.get2DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_2D, pairOption);
	    ymaker_qhzz.get2DHist_zmass_pairing(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_2D_2e2m, pairOption);
	    ymaker_qhzz.get2DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_2D, pairOption);
	    ymaker_qhzz.get2DHist_zmass_pairing(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_2D_2m2e, pairOption);

	    // Saving all the histograms to file
	    cout << "Saving all histograms to file " << filename << endl;
	    outFile.Write();
	  }



	  void makeHistogramsAllPairs(float mass, int igghid, int ivbfid, std::string filename) {
	    // mass      = Higgs mass
	    // igghid    = ID for the gHZZ signal MC
	    // ivbfid    = ID for the qHZZ signal MC

	    // Prepare strings to read files
	    stringstream mass_str_ss;
	    stringstream gghid_ss;
	    stringstream vbfid_ss;
	    mass_str_ss << mass;
	    gghid_ss << igghid;
	    vbfid_ss << ivbfid;
	    std::string mass_str = mass_str_ss.str();
	    std::string gghid = gghid_ss.str();
	    std::string vbfid = vbfid_ss.str();

	    //  std::cout << "Analyzing " << mass_str << " GeV mass point ... " << std::endl;

	    // Opening file to hold all the histograms
	    TFile outFile(filename.c_str(), "recreate");

	    // Creating the desired histograms
	    TH1F data_hist_mass("data_hist_mass", "Mass, all, data", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F data_hist_mass_4m("data_hist_mass_4m", "Mass, 4m, data", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F data_hist_mass_4e("data_hist_mass_4e", "Mass, 4e, data", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F data_hist_mass_2e2m("data_hist_mass_2e2m", "Mass, 2e2m, data", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F data_hist_allzmass("data_hist_allzmass", "Z1 mass, all, data", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F data_hist_allzmass_2m("data_hist_allzmass_2m", "Z1 mass, 2m, data", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F data_hist_allzmass_2e("data_hist_allzmass_2e", "Z1 mass, 2e, data", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F data_hist_mela("data_hist_mela", "MELA, all, data", 50, 0, 1);
	    TH1F data_hist_mela_4m("data_hist_mela_4m", "MELA, 4m, data", 50, 0, 1);
	    TH1F data_hist_mela_4e("data_hist_mela_4e", "MELA, 4e, data", 50, 0, 1);
	    TH1F data_hist_mela_2e2m("data_hist_mela_2e2m", "MELA, 2e2m, data", 50, 0, 1);

	    TH1F bkg_zxss_hist_mass("bkg_zxss_hist_mass", "Mass, all, ZX background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_zxss_hist_mass_4m("bkg_zxss_hist_mass_4m", "Mass, 4m, ZX background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_zxss_hist_mass_4e("bkg_zxss_hist_mass_4e", "Mass, 4e, ZX background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_zxss_hist_mass_2e2m("bkg_zxss_hist_mass_2e2m", "Mass, 2e2m, ZX background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_zxss_hist_allzmass("bkg_zxss_hist_allzmass", "Z1 mass, all, ZX background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_zxss_hist_allzmass_2m("bkg_zxss_hist_allzmass_2m", "Z1 mass, 2m, ZX background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_zxss_hist_allzmass_2e("bkg_zxss_hist_allzmass_2e", "Z1 mass, 2e, ZX background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_zxss_hist_mela("bkg_zxss_hist_mela", "MELA, all, ZX background", 50, 0, 1);
	    TH1F bkg_zxss_hist_mela_4m("bkg_zxss_hist_mela_4m", "MELA, 4m, ZX background", 50, 0, 1);
	    TH1F bkg_zxss_hist_mela_4e("bkg_zxss_hist_mela_4e", "MELA, 4e, ZX background", 50, 0, 1);
	    TH1F bkg_zxss_hist_mela_2e2m("bkg_zxss_hist_mela_2e2m", "MELA, 2e2m, ZX background", 50, 0, 1);

	    TH1F bkg_qqzz_hist_mass("bkg_qqzz_hist_mass", "Mass, all, qqZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_qqzz_hist_mass_4m("bkg_qqzz_hist_mass_4m", "Mass, 4m, qqZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_qqzz_hist_mass_4e("bkg_qqzz_hist_mass_4e", "Mass, 4e, qqZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_qqzz_hist_mass_2e2m("bkg_qqzz_hist_mass_2e2m", "Mass, 2e2m, qqZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_qqzz_hist_allzmass("bkg_qqzz_hist_allzmass", "Z1 mass, all, qqZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_qqzz_hist_allzmass_2m("bkg_qqzz_hist_allzmass_2m", "Z1 mass, 2m, qqZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_qqzz_hist_allzmass_2e("bkg_qqzz_hist_allzmass_2e", "Z1 mass, 2e, qqZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_qqzz_hist_mela("bkg_qqzz_hist_mela", "MELA, all, qqZZ background", 50, 0, 1);
	    TH1F bkg_qqzz_hist_mela_4m("bkg_qqzz_hist_mela_4m", "MELA, 4m, qqZZ background", 50, 0, 1);
	    TH1F bkg_qqzz_hist_mela_4e("bkg_qqzz_hist_mela_4e", "MELA, 4e, qqZZ background", 50, 0, 1);
	    TH1F bkg_qqzz_hist_mela_2e2m("bkg_qqzz_hist_mela_2e2m", "MELA, 2e2m, qqZZ background", 50, 0, 1);

	    TH1F bkg_ggzz_hist_mass("bkg_ggzz_hist_mass", "Mass, all, ggZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_ggzz_hist_mass_4m("bkg_ggzz_hist_mass_4m", "Mass, 4m, ggZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_ggzz_hist_mass_4e("bkg_ggzz_hist_mass_4e", "Mass, 4e, ggZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_ggzz_hist_mass_2e2m("bkg_ggzz_hist_mass_2e2m", "Mass, 2e2m, ggZZ background", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F bkg_ggzz_hist_allzmass("bkg_ggzz_hist_allzmass", "Z1 mass, all, ggZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_ggzz_hist_allzmass_2m("bkg_ggzz_hist_allzmass_2m", "Z1 mass, 2m, ggZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_ggzz_hist_allzmass_2e("bkg_ggzz_hist_allzmass_2e", "Z1 mass, 2e, ggZZ background", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F bkg_ggzz_hist_mela("bkg_ggzz_hist_mela", "MELA, all, ggZZ background", 50, 0, 1);
	    TH1F bkg_ggzz_hist_mela_4m("bkg_ggzz_hist_mela_4m", "MELA, 4m, ggZZ background", 50, 0, 1);
	    TH1F bkg_ggzz_hist_mela_4e("bkg_ggzz_hist_mela_4e", "MELA, 4e, ggZZ background", 50, 0, 1);
	    TH1F bkg_ggzz_hist_mela_2e2m("bkg_ggzz_hist_mela_2e2m", "MELA, 2e2m, ggZZ background", 50, 0, 1);

	    TH1F sig_ghzz_hist_mass("sig_ghzz_hist_mass", "Mass, all, gHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_ghzz_hist_mass_4m("sig_ghzz_hist_mass_4m", "Mass, 4m, gHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_ghzz_hist_mass_4e("sig_ghzz_hist_mass_4e", "Mass, 4e, gHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_ghzz_hist_mass_2e2m("sig_ghzz_hist_mass_2e2m", "Mass, 2e2m, gHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_ghzz_hist_allzmass("sig_ghzz_hist_allzmass", "Z1 mass, all, gHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F sig_ghzz_hist_allzmass_2m("sig_ghzz_hist_allzmass_2m", "Z1 mass, 2m, gHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F sig_ghzz_hist_allzmass_2e("sig_ghzz_hist_allzmass_2e", "Z1 mass, 2e, gHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F sig_ghzz_hist_mela("sig_ghzz_hist_mela", "MELA, all, gHZZ signal", 50, 0, 1);
	    TH1F sig_ghzz_hist_mela_4m("sig_ghzz_hist_mela_4m", "MELA, 4m, gHZZ signal", 50, 0, 1);
	    TH1F sig_ghzz_hist_mela_4e("sig_ghzz_hist_mela_4e", "MELA, 4e, gHZZ signal", 50, 0, 1);
	    TH1F sig_ghzz_hist_mela_2e2m("sig_ghzz_hist_mela_2e2m", "MELA, 2e2m, gHZZ signal", 50, 0, 1);

	    TH1F sig_qhzz_hist_mass("sig_qhzz_hist_mass", "Mass, all, qHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_qhzz_hist_mass_4m("sig_qhzz_hist_mass_4m", "Mass, 4m, qHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_qhzz_hist_mass_4e("sig_qhzz_hist_mass_4e", "Mass, 4e, qHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_qhzz_hist_mass_2e2m("sig_qhzz_hist_mass_2e2m", "Mass, 2e2m, qHZZ signal", (int) ( (massHigh - massLow)/3 ), massLow, massHigh);
	    TH1F sig_qhzz_hist_allzmass("sig_qhzz_hist_allzmass", "Z1 mass, all, qHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F sig_qhzz_hist_allzmass_2m("sig_qhzz_hist_allzmass_2m", "Z1 mass, 2m, qHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F sig_qhzz_hist_allzmass_2e("sig_qhzz_hist_allzmass_2e", "Z1 mass, 2e, qHZZ signal", (int) ( (z1max- z1min)/3 ), z1min, z1max);
	    TH1F sig_qhzz_hist_mela("sig_qhzz_hist_mela", "MELA, all, qHZZ signal", 50, 0, 1);
	    TH1F sig_qhzz_hist_mela_4m("sig_qhzz_hist_mela_4m", "MELA, 4m, qHZZ signal", 50, 0, 1);
	    TH1F sig_qhzz_hist_mela_4e("sig_qhzz_hist_mela_4e", "MELA, 4e, qHZZ signal", 50, 0, 1);
	    TH1F sig_qhzz_hist_mela_2e2m("sig_qhzz_hist_mela_2e2m", "MELA, 2e2m, qHZZ signal", 50, 0, 1);

	    // Using function defined in YieldMaker.h to fill histograms according to channels and cuts
	    ymaker_data.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass);
	    ymaker_data.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass_4m);
	    ymaker_data.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass);
	    ymaker_data.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass_4e);
	    ymaker_data.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass);
	    ymaker_data.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass_2e2m);
	    ymaker_data.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass);
	    ymaker_data.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mass_2e2m);

	    ymaker_zxss.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass);
	    ymaker_zxss.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass_4m);
	    ymaker_zxss.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass);
	    ymaker_zxss.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass_4e);
	    ymaker_zxss.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass);
	    ymaker_zxss.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass_2e2m);
	    ymaker_zxss.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass);
	    ymaker_zxss.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mass_2e2m);

	    ymaker_qqzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass);
	    ymaker_qqzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass_4m);
	    ymaker_qqzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass);
	    ymaker_qqzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass_4e);
	    ymaker_qqzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass);
	    ymaker_qqzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass_2e2m);
	    ymaker_qqzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass);
	    ymaker_qqzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mass_2e2m);

	    ymaker_ggzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass);
	    ymaker_ggzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass_4m);
	    ymaker_ggzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass);
	    ymaker_ggzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass_4e);
	    ymaker_ggzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass);
	    ymaker_ggzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass_2e2m);
	    ymaker_ggzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass);
	    ymaker_ggzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mass_2e2m);

	    ymaker_ghzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass);
	    ymaker_ghzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass_4m);
	    ymaker_ghzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass);
	    ymaker_ghzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass_4e);
	    ymaker_ghzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass);
	    ymaker_ghzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass_2e2m);
	    ymaker_ghzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass);
	    ymaker_ghzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mass_2e2m);

	    ymaker_qhzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass);
	    ymaker_qhzz.get1DHist_mass(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass_4m);
	    ymaker_qhzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass);
	    ymaker_qhzz.get1DHist_mass(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass_4e);
	    ymaker_qhzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass);
	    ymaker_qhzz.get1DHist_mass(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass_2e2m);
	    ymaker_qhzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass);
	    ymaker_qhzz.get1DHist_mass(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mass_2e2m);

	    ymaker_data.get1DHist_all_zmasses(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_allzmass);
	    ymaker_data.get1DHist_all_zmasses(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_allzmass_2m);
	    ymaker_data.get1DHist_all_zmasses(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_allzmass);
	    ymaker_data.get1DHist_all_zmasses(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_allzmass_2e);
	    ymaker_data.get1DHist_all_zmasses(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_allzmass);
	    ymaker_data.get1DHist_all_zmasses(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_allzmass_2e);
	    ymaker_data.get1DHist_all_zmasses(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_allzmass);
	    ymaker_data.get1DHist_all_zmasses(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_allzmass_2m);

	    ymaker_zxss.get1DHist_all_zmasses(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_allzmass);
	    ymaker_zxss.get1DHist_all_zmasses(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_allzmass_2m);
	    ymaker_zxss.get1DHist_all_zmasses(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_allzmass);
	    ymaker_zxss.get1DHist_all_zmasses(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_allzmass_2e);
	    ymaker_zxss.get1DHist_all_zmasses(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_allzmass);
	    ymaker_zxss.get1DHist_all_zmasses(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_allzmass_2e);
	    ymaker_zxss.get1DHist_all_zmasses(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_allzmass);
	    ymaker_zxss.get1DHist_all_zmasses(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_allzmass_2m);

	    ymaker_qqzz.get1DHist_all_zmasses(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_allzmass);
	    ymaker_qqzz.get1DHist_all_zmasses(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_allzmass_2m);
	    ymaker_qqzz.get1DHist_all_zmasses(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_allzmass);
	    ymaker_qqzz.get1DHist_all_zmasses(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_allzmass_2e);
	    ymaker_qqzz.get1DHist_all_zmasses(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_allzmass);
	    ymaker_qqzz.get1DHist_all_zmasses(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_allzmass_2e);
	    ymaker_qqzz.get1DHist_all_zmasses(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_allzmass);
	    ymaker_qqzz.get1DHist_all_zmasses(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_allzmass_2m);

	    ymaker_ggzz.get1DHist_all_zmasses(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass);
	    ymaker_ggzz.get1DHist_all_zmasses(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass_2m);
	    ymaker_ggzz.get1DHist_all_zmasses(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass);
	    ymaker_ggzz.get1DHist_all_zmasses(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass_2e);
	    ymaker_ggzz.get1DHist_all_zmasses(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass);
	    ymaker_ggzz.get1DHist_all_zmasses(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass_2e);
	    ymaker_ggzz.get1DHist_all_zmasses(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass);
	    ymaker_ggzz.get1DHist_all_zmasses(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass_2m);

	    ymaker_ggzz.get1DHist_all_zmasses(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass);
	    ymaker_ggzz.get1DHist_all_zmasses(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass_2m);
	    ymaker_ggzz.get1DHist_all_zmasses(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass);
	    ymaker_ggzz.get1DHist_all_zmasses(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass_2e);
	    ymaker_ggzz.get1DHist_all_zmasses(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass);
	    ymaker_ggzz.get1DHist_all_zmasses(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass_2e);
	    ymaker_ggzz.get1DHist_all_zmasses(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass);
	    ymaker_ggzz.get1DHist_all_zmasses(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_allzmass_2m);

	    ymaker_ghzz.get1DHist_all_zmasses(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_allzmass);
	    ymaker_ghzz.get1DHist_all_zmasses(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_allzmass_2m);
	    ymaker_ghzz.get1DHist_all_zmasses(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_allzmass);
	    ymaker_ghzz.get1DHist_all_zmasses(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_allzmass_2e);
	    ymaker_ghzz.get1DHist_all_zmasses(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_allzmass);
	    ymaker_ghzz.get1DHist_all_zmasses(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_allzmass_2e);
	    ymaker_ghzz.get1DHist_all_zmasses(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_allzmass);
	    ymaker_ghzz.get1DHist_all_zmasses(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_allzmass_2m);

	    ymaker_qhzz.get1DHist_all_zmasses(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_allzmass);
	    ymaker_qhzz.get1DHist_all_zmasses(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_allzmass_2m);
	    ymaker_qhzz.get1DHist_all_zmasses(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_allzmass);
	    ymaker_qhzz.get1DHist_all_zmasses(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_allzmass_2e);
	    ymaker_qhzz.get1DHist_all_zmasses(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_allzmass);
	    ymaker_qhzz.get1DHist_all_zmasses(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_allzmass_2e);
	    ymaker_qhzz.get1DHist_all_zmasses(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_allzmass);
	    ymaker_qhzz.get1DHist_all_zmasses(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_allzmass_2m);

	    ymaker_data.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela);
	    ymaker_data.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela_4m);
	    ymaker_data.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela);
	    ymaker_data.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela_4e);
	    ymaker_data.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela);
	    ymaker_data.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela_2e2m);
	    ymaker_data.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela);
	    ymaker_data.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &data_hist_mela_2e2m);

	    ymaker_zxss.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela);
	    ymaker_zxss.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela_4m);
	    ymaker_zxss.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela);
	    ymaker_zxss.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela_4e);
	    ymaker_zxss.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela);
	    ymaker_zxss.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela_2e2m);
	    ymaker_zxss.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela);
	    ymaker_zxss.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_zxss_hist_mela_2e2m);

	    ymaker_qqzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela);
	    ymaker_qqzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela_4m);
	    ymaker_qqzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela);
	    ymaker_qqzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela_4e);
	    ymaker_qqzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela);
	    ymaker_qqzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela_2e2m);
	    ymaker_qqzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela);
	    ymaker_qqzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_qqzz_hist_mela_2e2m);

	    ymaker_ggzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela);
	    ymaker_ggzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela_4m);
	    ymaker_ggzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela);
	    ymaker_ggzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela_4e);
	    ymaker_ggzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela);
	    ymaker_ggzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela_2e2m);
	    ymaker_ggzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela);
	    ymaker_ggzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &bkg_ggzz_hist_mela_2e2m);

	    ymaker_ghzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela);
	    ymaker_ghzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela_4m);
	    ymaker_ghzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela);
	    ymaker_ghzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela_4e);
	    ymaker_ghzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela);
	    ymaker_ghzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela_2e2m);
	    ymaker_ghzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela);
	    ymaker_ghzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_ghzz_hist_mela_2e2m);

	    ymaker_qhzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela);
	    ymaker_qhzz.get1DHist_mela(0, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela_4m);
	    ymaker_qhzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela);
	    ymaker_qhzz.get1DHist_mela(1, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela_4e);
	    ymaker_qhzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela);
	    ymaker_qhzz.get1DHist_mela(2, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela_2e2m);
	    ymaker_qhzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela);
	    ymaker_qhzz.get1DHist_mela(3, z1min, z1max, z2min, z2max, massLow, massHigh, melacut, &sig_qhzz_hist_mela_2e2m);

	    // Saving all the histograms to file
	    cout << "Saving all histograms to file " << filename << endl;
	    outFile.Write();
	  }
	};


void doHZZHistogramsPairing(char* directoryName, bool isCut = 0) {

  //==============================================================================================================//
  //				READING DATA FROM FILES AND SAVING HISTOGRAMS					//	
  //______________________________________________________________________________________________________________//

  // Creating MassPointInfoPairing for 7 TeV  with desired parameters
  MassPointInfoPairing hmpi7;
  hmpi7.lumi = 5.05;
  hmpi7.z1min = 0.;
  hmpi7.z1max = 120.;
  hmpi7.z2min = 0.;
  hmpi7.z2max = 120.;

  // If isCut is selected, then do a cut on the mass to the signal region
  if (isCut) {
    hmpi7.massLow  = 118.;
    hmpi7.massHigh = 135;
  } else {
    hmpi7.massLow  = 100.;
    hmpi7.massHigh = 220.;
  }

  hmpi7.melacut = -1.0;
  hmpi7.doSS = true;
  hmpi7.do7TeV = true;
  hmpi7.treeFolder = "/mnt/hadoop/user/sixie/ntuples/hzz4l/step2/HZZ4L_42X_S1_V04_S2_V01_adish/";

  // Initialize common definitions (from macro/scales2.h)
  init(hmpi7.do7TeV);

  FakeRateCalculator FR_7TeV(hmpi7.treeFolder+"hzzTree.root", hmpi7.do7TeV, 40, 120, 0.0, 0.0, true);
  hmpi7.FR = FR_7TeV;


  // Filling data, background and signal YieldMakers
  hmpi7.ymaker_data.fill(hmpi7.treeFolder+"hzzTree.root");
  hmpi7.ymaker_zxss.fill(hmpi7.treeFolder+"hzzTree.root"       , 1.0, hmpi7.FR, hmpi7.doSS);
  hmpi7.ymaker_qqzz.fill(hmpi7.treeFolder+"hzzTree_id121.root" , getBkgXsec(121)*hmpi7.lumi/evt_7TeV[121], 0.0, false);
  hmpi7.ymaker_qqzz.fill(hmpi7.treeFolder+"hzzTree_id122.root" , getBkgXsec(122)*hmpi7.lumi/evt_7TeV[122], 0.0, false);
  hmpi7.ymaker_qqzz.fill(hmpi7.treeFolder+"hzzTree_id123.root" , getBkgXsec(123)*hmpi7.lumi/evt_7TeV[123], 0.0, false);
  hmpi7.ymaker_qqzz.fill(hmpi7.treeFolder+"hzzTree_id124.root" , getBkgXsec(124)*hmpi7.lumi/evt_7TeV[124], 0.0, false);
  hmpi7.ymaker_qqzz.fill(hmpi7.treeFolder+"hzzTree_id125.root" , getBkgXsec(125)*hmpi7.lumi/evt_7TeV[125], 0.0, false);
  hmpi7.ymaker_qqzz.fill(hmpi7.treeFolder+"hzzTree_id126.root" , getBkgXsec(126)*hmpi7.lumi/evt_7TeV[126], 0.0, false);
  hmpi7.ymaker_ggzz.fill(hmpi7.treeFolder+"hzzTree_id101.root" , getBkgXsec(101)*hmpi7.lumi/evt_7TeV[101], 0.0, false);
  hmpi7.ymaker_ggzz.fill(hmpi7.treeFolder+"hzzTree_id100.root" , getBkgXsec(100)*hmpi7.lumi/evt_7TeV[100], 0.0, false);
  hmpi7.ymaker_ghzz.fill(hmpi7.treeFolder+"hzzTree_id201.root" , getXsecggH(120)*hmpi7.lumi/evt_7TeV[201], 0.0, true );
  hmpi7.ymaker_qhzz.fill(hmpi7.treeFolder+"hzzTree_id251.root" , getXsecVBF(120)*hmpi7.lumi/evt_7TeV[251], 0.0, true );

  // Generating only histograms for 120 datapoint (signal from here will not be used)
  cout << "Making histograms for 7 TeV data..." << endl;
  std::string fileString_7TeV = std::string(directoryName) + "AllHistograms_7TeV";
  hmpi7.makeHistogramsAllPairs(120, 201, 251, fileString_7TeV + "_allpairs.root");
  hmpi7.makeHistogramsPairing(120, 201, 251, fileString_7TeV + "_pairing0.root", 0);
  hmpi7.makeHistogramsPairing(120, 201, 251, fileString_7TeV + "_pairing1.root", 1);
  hmpi7.makeHistogramsPairing(120, 201, 251, fileString_7TeV + "_pairing2.root", 2);
  hmpi7.makeHistogramsPairing(120, 201, 251, fileString_7TeV + "_pairing3.root", 3);
  hmpi7.makeHistogramsPairing(120, 201, 251, fileString_7TeV + "_pairing4.root", 4);
  hmpi7.makeHistogramsPairing(120, 201, 251, fileString_7TeV + "_pairing5.root", 5);
  hmpi7.makeHistogramsPairing(120, 201, 251, fileString_7TeV + "_pairing6.root", 6);
  hmpi7.makeHistogramsPairing(120, 201, 251, fileString_7TeV + "_pairing7.root", 7);

  // Doing the same for the 8 TeV
  MassPointInfoPairing hmpi8;
  hmpi8.lumi = 5.26; 
  hmpi8.z1min = 0.;
  hmpi8.z1max = 120.;
  hmpi8.z2min = 0.;
  hmpi8.z2max = 120.;

  // If isCut is selected, then do a cut on the mass to the signal region
  if (isCut) {
    hmpi8.massLow  = 118.;
    hmpi8.massHigh = 135;
  } else {
    hmpi8.massLow  = 100.;
    hmpi8.massHigh = 220.;
  }

  hmpi8.melacut = -1.0;
  hmpi8.doSS = true;    
  hmpi8.do7TeV = false;
  hmpi8.treeFolder = "/mnt/hadoop/user/sixie/ntuples/hzz4l/step2/ichep2012/";

  init(hmpi8.do7TeV);

  FakeRateCalculator FR_8TeV(hmpi8.treeFolder+"data/hzzTree.root", hmpi8.do7TeV, 40, 120, 0.0, 0.0, true);
  hmpi8.FR = FR_8TeV;

  hmpi8.ymaker_data.fill(hmpi8.treeFolder+"data/hzzTree.root");
  hmpi8.ymaker_zxss.fill(hmpi8.treeFolder+"data/hzzTree.root"       , 1.0, hmpi8.FR, hmpi8.doSS);
  hmpi8.ymaker_qqzz.fill(hmpi8.treeFolder+"mc/hzzTree_id102.root" , getBkgXsec(102)*hmpi8.lumi/evt_8TeV[102], 0.0, false);
  hmpi8.ymaker_qqzz.fill(hmpi8.treeFolder+"mc/hzzTree_id103.root" , getBkgXsec(103)*hmpi8.lumi/evt_8TeV[103], 0.0, false);
  hmpi8.ymaker_qqzz.fill(hmpi8.treeFolder+"mc/hzzTree_id104.root" , getBkgXsec(104)*hmpi8.lumi/evt_8TeV[104], 0.0, false);
  hmpi8.ymaker_qqzz.fill(hmpi8.treeFolder+"mc/hzzTree_id105.root" , getBkgXsec(105)*hmpi8.lumi/evt_8TeV[105], 0.0, false);
  hmpi8.ymaker_qqzz.fill(hmpi8.treeFolder+"mc/hzzTree_id106.root" , getBkgXsec(106)*hmpi8.lumi/evt_8TeV[106], 0.0, false);
  hmpi8.ymaker_qqzz.fill(hmpi8.treeFolder+"mc/hzzTree_id107.root" , getBkgXsec(107)*hmpi8.lumi/evt_8TeV[107], 0.0, false);
  hmpi8.ymaker_ggzz.fill(hmpi8.treeFolder+"mc/hzzTree_id101.root" , getBkgXsec(101)*hmpi8.lumi/evt_8TeV[101], 0.0, false);
  hmpi8.ymaker_ggzz.fill(hmpi8.treeFolder+"mc/hzzTree_id100.root" , getBkgXsec(100)*hmpi8.lumi/evt_8TeV[100], 0.0, false);
  hmpi8.ymaker_ghzz.fill(hmpi8.treeFolder+"mc/hzzTree_id1126.root", getXsecggH(126)*hmpi8.lumi/evt_8TeV[1126], 0.0, true );
  hmpi8.ymaker_qhzz.fill(hmpi8.treeFolder+"mc/hzzTree_id2126.root", getXsecVBF(126)*hmpi8.lumi/evt_8TeV[2126], 0.0, true );

  // Analyzing the 126 GeV point for 2012 data
  cout << "Making histograms for 8 TeV data..." << endl;
  std::string fileString_8TeV = std::string(directoryName) + "AllHistograms_8TeV";
  hmpi8.makeHistogramsAllPairs(126, 1126, 2126, fileString_8TeV + "_allpairs.root");
  hmpi8.makeHistogramsPairing(126, 1126, 2126, fileString_8TeV + "_pairing0.root", 0);
  hmpi8.makeHistogramsPairing(126, 1126, 2126, fileString_8TeV + "_pairing1.root", 1);
  hmpi8.makeHistogramsPairing(126, 1126, 2126, fileString_8TeV + "_pairing2.root", 2);
  hmpi8.makeHistogramsPairing(126, 1126, 2126, fileString_8TeV + "_pairing3.root", 3);
  hmpi8.makeHistogramsPairing(126, 1126, 2126, fileString_8TeV + "_pairing4.root", 4);
  hmpi8.makeHistogramsPairing(126, 1126, 2126, fileString_8TeV + "_pairing5.root", 5);
  hmpi8.makeHistogramsPairing(126, 1126, 2126, fileString_8TeV + "_pairing6.root", 6);
  hmpi8.makeHistogramsPairing(126, 1126, 2126, fileString_8TeV + "_pairing7.root", 7);


  //==============================================================================================================//
  //					PLOTTING HISTOGRAMS							//	
  //______________________________________________________________________________________________________________//

  cout << "Plotting histograms together, and saving it to directory " << directoryName <<endl;
  doHZZHistogramsAllPairsPlot(fileString_7TeV + "_allpairs.root", fileString_8TeV + "_allpairs.root", std::string(directoryName) + "allpairs/", 0, hmpi7.lumi, hmpi8.lumi);
  doHZZHistogramsPlot(fileString_7TeV + "_pairing0.root", fileString_8TeV + "_pairing0.root", std::string(directoryName) + "pairing0/", 0, hmpi7.lumi, hmpi8.lumi);
  doHZZHistogramsPlot(fileString_7TeV + "_pairing1.root", fileString_8TeV + "_pairing1.root", std::string(directoryName) + "pairing1/", 0, hmpi7.lumi, hmpi8.lumi);
  doHZZHistogramsPlot(fileString_7TeV + "_pairing2.root", fileString_8TeV + "_pairing2.root", std::string(directoryName) + "pairing2/", 0, hmpi7.lumi, hmpi8.lumi);
  doHZZHistogramsPlot(fileString_7TeV + "_pairing3.root", fileString_8TeV + "_pairing3.root", std::string(directoryName) + "pairing3/", 0, hmpi7.lumi, hmpi8.lumi);
  doHZZHistogramsPlot(fileString_7TeV + "_pairing4.root", fileString_8TeV + "_pairing4.root", std::string(directoryName) + "pairing4/", 0, hmpi7.lumi, hmpi8.lumi);
  doHZZHistogramsPlot(fileString_7TeV + "_pairing5.root", fileString_8TeV + "_pairing5.root", std::string(directoryName) + "pairing5/", 0, hmpi7.lumi, hmpi8.lumi);
  doHZZHistogramsPlot(fileString_7TeV + "_pairing6.root", fileString_8TeV + "_pairing6.root", std::string(directoryName) + "pairing6/", 0, hmpi7.lumi, hmpi8.lumi);
  doHZZHistogramsPlot(fileString_7TeV + "_pairing7.root", fileString_8TeV + "_pairing7.root", std::string(directoryName) + "pairing7/", 0, hmpi7.lumi, hmpi8.lumi);
}
