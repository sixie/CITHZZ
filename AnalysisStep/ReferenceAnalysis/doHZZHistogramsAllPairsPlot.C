#include <sstream>
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "THStack.h"
#include "TMath.h"
#include "TFile.h"

void doHZZHistogramsAllPairsPlot (std::string filename2011, std::string filename2012, std::string outDirectory, bool plotZXbackground = 1, float lumi_2011 = 5.05, float lumi_2012 = 5.26) {

  // Setting up graphic style
  gStyle->SetOptStat(0);

  // Retrieving histograms from file
  TFile* file_2011 = new TFile(filename2011.c_str());
  TFile* file_2012 = new TFile(filename2012.c_str());

  TH1F* data_mass_2011 = (TH1F*) file_2011->Get("data_hist_mass");
  TH1F* data_mass_4m_2011 = (TH1F*) file_2011->Get("data_hist_mass_4m");
  TH1F* data_mass_4e_2011 = (TH1F*) file_2011->Get("data_hist_mass_4e");
  TH1F* data_mass_2e2m_2011 = (TH1F*) file_2011->Get("data_hist_mass_2e2m");
  TH1F* data_allzmass_2011 = (TH1F*) file_2011->Get("data_hist_allzmass");
  TH1F* data_allzmass_2m_2011 = (TH1F*) file_2011->Get("data_hist_allzmass_2m");
  TH1F* data_allzmass_2e_2011 = (TH1F*) file_2011->Get("data_hist_allzmass_2e");
  TH1F* data_mela_2011 = (TH1F*) file_2011->Get("data_hist_mela");
  TH1F* data_mela_4m_2011 = (TH1F*) file_2011->Get("data_hist_mela_4m");
  TH1F* data_mela_4e_2011 = (TH1F*) file_2011->Get("data_hist_mela_4e");
  TH1F* data_mela_2e2m_2011 = (TH1F*) file_2011->Get("data_hist_mela_2e2m");
  TH1F* data_mass_2012 = (TH1F*) file_2012->Get("data_hist_mass");
  TH1F* data_mass_4m_2012 = (TH1F*) file_2012->Get("data_hist_mass_4m");
  TH1F* data_mass_4e_2012 = (TH1F*) file_2012->Get("data_hist_mass_4e");
  TH1F* data_mass_2e2m_2012 = (TH1F*) file_2012->Get("data_hist_mass_2e2m");
  TH1F* data_allzmass_2012 = (TH1F*) file_2012->Get("data_hist_allzmass");
  TH1F* data_allzmass_2m_2012 = (TH1F*) file_2012->Get("data_hist_allzmass_2m");
  TH1F* data_allzmass_2e_2012 = (TH1F*) file_2012->Get("data_hist_allzmass_2e");
  TH1F* data_mela_2012 = (TH1F*) file_2012->Get("data_hist_mela");
  TH1F* data_mela_4m_2012 = (TH1F*) file_2012->Get("data_hist_mela_4m");
  TH1F* data_mela_4e_2012 = (TH1F*) file_2012->Get("data_hist_mela_4e");
  TH1F* data_mela_2e2m_2012 = (TH1F*) file_2012->Get("data_hist_mela_2e2m");

  TH1F *bkg_zxss_mass_2011, *bkg_zxss_mass_4m_2011, *bkg_zxss_mass_4e_2011, *bkg_zxss_mass_2e2m_2011;
  TH1F *bkg_zxss_allzmass_2011, *bkg_zxss_allzmass_2m_2011, *bkg_zxss_allzmass_2e_2011;
  TH1F *bkg_zxss_mela_2011, *bkg_zxss_mela_4m_2011, *bkg_zxss_mela_4e_2011, *bkg_zxss_mela_2e2m_2011;
  TH1F *bkg_zxss_mass_2012, *bkg_zxss_mass_4m_2012, *bkg_zxss_mass_4e_2012, *bkg_zxss_mass_2e2m_2012;
  TH1F *bkg_zxss_allzmass_2012, *bkg_zxss_allzmass_2m_2012, *bkg_zxss_allzmass_2e_2012;
  TH1F *bkg_zxss_mela_2012, *bkg_zxss_mela_4m_2012, *bkg_zxss_mela_4e_2012, *bkg_zxss_mela_2e2m_2012;

  if (plotZXbackground) {
  bkg_zxss_mass_2011 = (TH1F*) file_2011->Get("bkg_zxss_hist_mass");
  bkg_zxss_mass_4m_2011 = (TH1F*) file_2011->Get("bkg_zxss_hist_mass_4m");
  bkg_zxss_mass_4e_2011 = (TH1F*) file_2011->Get("bkg_zxss_hist_mass_4e");
  bkg_zxss_mass_2e2m_2011 = (TH1F*) file_2011->Get("bkg_zxss_hist_mass_2e2m");
  bkg_zxss_allzmass_2011 = (TH1F*) file_2011->Get("bkg_zxss_hist_allzmass");
  bkg_zxss_allzmass_2m_2011 = (TH1F*) file_2011->Get("bkg_zxss_hist_allzmass_2m");
  bkg_zxss_allzmass_2e_2011 = (TH1F*) file_2011->Get("bkg_zxss_hist_allzmass_2e");
  bkg_zxss_mela_2011 = (TH1F*) file_2011->Get("bkg_zxss_hist_mela");
  bkg_zxss_mela_4m_2011 = (TH1F*) file_2011->Get("bkg_zxss_hist_mela_4m");
  bkg_zxss_mela_4e_2011 = (TH1F*) file_2011->Get("bkg_zxss_hist_mela_4e");
  bkg_zxss_mela_2e2m_2011 = (TH1F*) file_2011->Get("bkg_zxss_hist_mela_2e2m");
  bkg_zxss_mass_2012 = (TH1F*) file_2012->Get("bkg_zxss_hist_mass");
  bkg_zxss_mass_4m_2012 = (TH1F*) file_2012->Get("bkg_zxss_hist_mass_4m");
  bkg_zxss_mass_4e_2012 = (TH1F*) file_2012->Get("bkg_zxss_hist_mass_4e");
  bkg_zxss_mass_2e2m_2012 = (TH1F*) file_2012->Get("bkg_zxss_hist_mass_2e2m");
  bkg_zxss_allzmass_2012 = (TH1F*) file_2012->Get("bkg_zxss_hist_allzmass");
  bkg_zxss_allzmass_2m_2012 = (TH1F*) file_2012->Get("bkg_zxss_hist_allzmass_2m");
  bkg_zxss_allzmass_2e_2012 = (TH1F*) file_2012->Get("bkg_zxss_hist_allzmass_2e");
  bkg_zxss_mela_2012 = (TH1F*) file_2012->Get("bkg_zxss_hist_mela");
  bkg_zxss_mela_4m_2012 = (TH1F*) file_2012->Get("bkg_zxss_hist_mela_4m");
  bkg_zxss_mela_4e_2012 = (TH1F*) file_2012->Get("bkg_zxss_hist_mela_4e");
  bkg_zxss_mela_2e2m_2012 = (TH1F*) file_2012->Get("bkg_zxss_hist_mela_2e2m");
  }

  TH1F* bkg_ggzz_mass_2011 = (TH1F*) file_2011->Get("bkg_ggzz_hist_mass");
  TH1F* bkg_ggzz_mass_4m_2011 = (TH1F*) file_2011->Get("bkg_ggzz_hist_mass_4m");
  TH1F* bkg_ggzz_mass_4e_2011 = (TH1F*) file_2011->Get("bkg_ggzz_hist_mass_4e");
  TH1F* bkg_ggzz_mass_2e2m_2011 = (TH1F*) file_2011->Get("bkg_ggzz_hist_mass_2e2m");
  TH1F* bkg_ggzz_allzmass_2011 = (TH1F*) file_2011->Get("bkg_ggzz_hist_allzmass");
  TH1F* bkg_ggzz_allzmass_2m_2011 = (TH1F*) file_2011->Get("bkg_ggzz_hist_allzmass_2m");
  TH1F* bkg_ggzz_allzmass_2e_2011 = (TH1F*) file_2011->Get("bkg_ggzz_hist_allzmass_2e");
  TH1F* bkg_ggzz_mela_2011 = (TH1F*) file_2011->Get("bkg_ggzz_hist_mela");
  TH1F* bkg_ggzz_mela_4m_2011 = (TH1F*) file_2011->Get("bkg_ggzz_hist_mela_4m");
  TH1F* bkg_ggzz_mela_4e_2011 = (TH1F*) file_2011->Get("bkg_ggzz_hist_mela_4e");
  TH1F* bkg_ggzz_mela_2e2m_2011 = (TH1F*) file_2011->Get("bkg_ggzz_hist_mela_2e2m");
  TH1F* bkg_ggzz_mass_2012 = (TH1F*) file_2012->Get("bkg_ggzz_hist_mass");
  TH1F* bkg_ggzz_mass_4m_2012 = (TH1F*) file_2012->Get("bkg_ggzz_hist_mass_4m");
  TH1F* bkg_ggzz_mass_4e_2012 = (TH1F*) file_2012->Get("bkg_ggzz_hist_mass_4e");
  TH1F* bkg_ggzz_mass_2e2m_2012 = (TH1F*) file_2012->Get("bkg_ggzz_hist_mass_2e2m");
  TH1F* bkg_ggzz_allzmass_2012 = (TH1F*) file_2012->Get("bkg_ggzz_hist_allzmass");
  TH1F* bkg_ggzz_allzmass_2m_2012 = (TH1F*) file_2012->Get("bkg_ggzz_hist_allzmass_2m");
  TH1F* bkg_ggzz_allzmass_2e_2012 = (TH1F*) file_2012->Get("bkg_ggzz_hist_allzmass_2e");
  TH1F* bkg_ggzz_mela_2012 = (TH1F*) file_2012->Get("bkg_ggzz_hist_mela");
  TH1F* bkg_ggzz_mela_4m_2012 = (TH1F*) file_2012->Get("bkg_ggzz_hist_mela_4m");
  TH1F* bkg_ggzz_mela_4e_2012 = (TH1F*) file_2012->Get("bkg_ggzz_hist_mela_4e");
  TH1F* bkg_ggzz_mela_2e2m_2012 = (TH1F*) file_2012->Get("bkg_ggzz_hist_mela_2e2m");

  TH1F* bkg_qqzz_mass_2011 = (TH1F*) file_2011->Get("bkg_qqzz_hist_mass");
  TH1F* bkg_qqzz_mass_4m_2011 = (TH1F*) file_2011->Get("bkg_qqzz_hist_mass_4m");
  TH1F* bkg_qqzz_mass_4e_2011 = (TH1F*) file_2011->Get("bkg_qqzz_hist_mass_4e");
  TH1F* bkg_qqzz_mass_2e2m_2011 = (TH1F*) file_2011->Get("bkg_qqzz_hist_mass_2e2m");
  TH1F* bkg_qqzz_allzmass_2011 = (TH1F*) file_2011->Get("bkg_qqzz_hist_allzmass");
  TH1F* bkg_qqzz_allzmass_2m_2011 = (TH1F*) file_2011->Get("bkg_qqzz_hist_allzmass_2m");
  TH1F* bkg_qqzz_allzmass_2e_2011 = (TH1F*) file_2011->Get("bkg_qqzz_hist_allzmass_2e");
  TH1F* bkg_qqzz_mela_2011 = (TH1F*) file_2011->Get("bkg_qqzz_hist_mela");
  TH1F* bkg_qqzz_mela_4m_2011 = (TH1F*) file_2011->Get("bkg_qqzz_hist_mela_4m");
  TH1F* bkg_qqzz_mela_4e_2011 = (TH1F*) file_2011->Get("bkg_qqzz_hist_mela_4e");
  TH1F* bkg_qqzz_mela_2e2m_2011 = (TH1F*) file_2011->Get("bkg_qqzz_hist_mela_2e2m");
  TH1F* bkg_qqzz_mass_2012 = (TH1F*) file_2012->Get("bkg_qqzz_hist_mass");
  TH1F* bkg_qqzz_mass_4m_2012 = (TH1F*) file_2012->Get("bkg_qqzz_hist_mass_4m");
  TH1F* bkg_qqzz_mass_4e_2012 = (TH1F*) file_2012->Get("bkg_qqzz_hist_mass_4e");
  TH1F* bkg_qqzz_mass_2e2m_2012 = (TH1F*) file_2012->Get("bkg_qqzz_hist_mass_2e2m");
  TH1F* bkg_qqzz_allzmass_2012 = (TH1F*) file_2012->Get("bkg_qqzz_hist_allzmass");
  TH1F* bkg_qqzz_allzmass_2m_2012 = (TH1F*) file_2012->Get("bkg_qqzz_hist_allzmass_2m");
  TH1F* bkg_qqzz_allzmass_2e_2012 = (TH1F*) file_2012->Get("bkg_qqzz_hist_allzmass_2e");
  TH1F* bkg_qqzz_mela_2012 = (TH1F*) file_2012->Get("bkg_qqzz_hist_mela");
  TH1F* bkg_qqzz_mela_4m_2012 = (TH1F*) file_2012->Get("bkg_qqzz_hist_mela_4m");
  TH1F* bkg_qqzz_mela_4e_2012 = (TH1F*) file_2012->Get("bkg_qqzz_hist_mela_4e");
  TH1F* bkg_qqzz_mela_2e2m_2012 = (TH1F*) file_2012->Get("bkg_qqzz_hist_mela_2e2m");

  TH1F* sig_qhzz_mass_2012 = (TH1F*) file_2012->Get("sig_qhzz_hist_mass");
  TH1F* sig_qhzz_mass_4m_2012 = (TH1F*) file_2012->Get("sig_qhzz_hist_mass_4m");
  TH1F* sig_qhzz_mass_4e_2012 = (TH1F*) file_2012->Get("sig_qhzz_hist_mass_4e");
  TH1F* sig_qhzz_mass_2e2m_2012 = (TH1F*) file_2012->Get("sig_qhzz_hist_mass_2e2m");
  TH1F* sig_qhzz_allzmass_2012 = (TH1F*) file_2012->Get("sig_qhzz_hist_allzmass");
  TH1F* sig_qhzz_allzmass_2m_2012 = (TH1F*) file_2012->Get("sig_qhzz_hist_allzmass_2m");
  TH1F* sig_qhzz_allzmass_2e_2012 = (TH1F*) file_2012->Get("sig_qhzz_hist_allzmass_2e");
  TH1F* sig_qhzz_mela_2012 = (TH1F*) file_2012->Get("sig_qhzz_hist_mela");
  TH1F* sig_qhzz_mela_4m_2012 = (TH1F*) file_2012->Get("sig_qhzz_hist_mela_4m");
  TH1F* sig_qhzz_mela_4e_2012 = (TH1F*) file_2012->Get("sig_qhzz_hist_mela_4e");
  TH1F* sig_qhzz_mela_2e2m_2012 = (TH1F*) file_2012->Get("sig_qhzz_hist_mela_2e2m");

  TH1F* sig_ghzz_mass_2012 = (TH1F*) file_2012->Get("sig_ghzz_hist_mass");
  TH1F* sig_ghzz_mass_4m_2012 = (TH1F*) file_2012->Get("sig_ghzz_hist_mass_4m");
  TH1F* sig_ghzz_mass_4e_2012 = (TH1F*) file_2012->Get("sig_ghzz_hist_mass_4e");
  TH1F* sig_ghzz_mass_2e2m_2012 = (TH1F*) file_2012->Get("sig_ghzz_hist_mass_2e2m");
  TH1F* sig_ghzz_allzmass_2012 = (TH1F*) file_2012->Get("sig_ghzz_hist_allzmass");
  TH1F* sig_ghzz_allzmass_2m_2012 = (TH1F*) file_2012->Get("sig_ghzz_hist_allzmass_2m");
  TH1F* sig_ghzz_allzmass_2e_2012 = (TH1F*) file_2012->Get("sig_ghzz_hist_allzmass_2e");
  TH1F* sig_ghzz_mela_2012 = (TH1F*) file_2012->Get("sig_ghzz_hist_mela");
  TH1F* sig_ghzz_mela_4m_2012 = (TH1F*) file_2012->Get("sig_ghzz_hist_mela_4m");
  TH1F* sig_ghzz_mela_4e_2012 = (TH1F*) file_2012->Get("sig_ghzz_hist_mela_4e");
  TH1F* sig_ghzz_mela_2e2m_2012 = (TH1F*) file_2012->Get("sig_ghzz_hist_mela_2e2m");

  // Now adding up all 2011 and 2012 data
  TH1F* data_mass = (TH1F*) data_mass_2011->Clone("data_mass"); data_mass->Add(data_mass_2012);
  TH1F* data_mass_4m = (TH1F*) data_mass_4m_2011->Clone("data_mass_4m"); data_mass_4m->Add(data_mass_4m_2012);
  TH1F* data_mass_4e = (TH1F*) data_mass_4e_2011->Clone("data_mass_4e"); data_mass_4e->Add(data_mass_4e_2012);
  TH1F* data_mass_2e2m = (TH1F*) data_mass_2e2m_2011->Clone("data_mass_2e2m"); data_mass_2e2m->Add(data_mass_2e2m_2012);
  TH1F* data_allzmass = (TH1F*) data_allzmass_2011->Clone("data_allzmass"); data_allzmass->Add(data_allzmass_2012);
  TH1F* data_allzmass_2m = (TH1F*) data_allzmass_2m_2011->Clone("data_allzmass_2m"); data_allzmass_2m->Add(data_allzmass_2m_2012);
  TH1F* data_allzmass_2e = (TH1F*) data_allzmass_2e_2011->Clone("data_allzmass_2e"); data_allzmass_2e->Add(data_allzmass_2e_2012);
  TH1F* data_mela = (TH1F*) data_mela_2011->Clone("data_mela"); data_mela->Add(data_mela_2012);
  TH1F* data_mela_4m = (TH1F*) data_mela_4m_2011->Clone("data_mela_4m"); data_mela_4m->Add(data_mela_4m_2012);
  TH1F* data_mela_4e = (TH1F*) data_mela_4e_2011->Clone("data_mela_4e"); data_mela_4e->Add(data_mela_4e_2012);
  TH1F* data_mela_2e2m = (TH1F*) data_mela_2e2m_2011->Clone("data_mela_2e2m"); data_mela_2e2m->Add(data_mela_2e2m_2012);

  TH1F *bkg_zxss_mass, *bkg_zxss_mass_4m, *bkg_zxss_mass_4e, *bkg_zxss_mass_2e2m;
  TH1F *bkg_zxss_allzmass, *bkg_zxss_allzmass_2m, *bkg_zxss_allzmass_2e;
  TH1F *bkg_zxss_mela, *bkg_zxss_mela_4m, *bkg_zxss_mela_4e, *bkg_zxss_mela_2e2m;

  if (plotZXbackground) {
  bkg_zxss_mass = (TH1F*) bkg_zxss_mass_2011->Clone("bkg_zxss_mass"); bkg_zxss_mass->Add(bkg_zxss_mass_2012);
  bkg_zxss_mass_4m = (TH1F*) bkg_zxss_mass_4m_2011->Clone("bkg_zxss_mass_4m"); bkg_zxss_mass_4m->Add(bkg_zxss_mass_4m_2012);
  bkg_zxss_mass_4e = (TH1F*) bkg_zxss_mass_4e_2011->Clone("bkg_zxss_mass_4e"); bkg_zxss_mass_4e->Add(bkg_zxss_mass_4e_2012);
  bkg_zxss_mass_2e2m = (TH1F*) bkg_zxss_mass_2e2m_2011->Clone("bkg_zxss_mass_2e2m"); bkg_zxss_mass_2e2m->Add(bkg_zxss_mass_2e2m_2012);
  bkg_zxss_allzmass = (TH1F*) bkg_zxss_allzmass_2011->Clone("bkg_zxss_allzmass"); bkg_zxss_allzmass->Add(bkg_zxss_allzmass_2012);
  bkg_zxss_allzmass_2m = (TH1F*) bkg_zxss_allzmass_2m_2011->Clone("bkg_zxss_allzmass_2m"); bkg_zxss_allzmass_2m->Add(bkg_zxss_allzmass_2m_2012);
  bkg_zxss_allzmass_2e = (TH1F*) bkg_zxss_allzmass_2e_2011->Clone("bkg_zxss_allzmass_2e"); bkg_zxss_allzmass_2e->Add(bkg_zxss_allzmass_2e_2012);
  bkg_zxss_mela = (TH1F*) bkg_zxss_mela_2011->Clone("bkg_zxss_mela"); bkg_zxss_mela->Add(bkg_zxss_mela_2012);
  bkg_zxss_mela_4m = (TH1F*) bkg_zxss_mela_4m_2011->Clone("bkg_zxss_mela_4m"); bkg_zxss_mela_4m->Add(bkg_zxss_mela_4m_2012);
  bkg_zxss_mela_4e = (TH1F*) bkg_zxss_mela_4e_2011->Clone("bkg_zxss_mela_4e"); bkg_zxss_mela_4e->Add(bkg_zxss_mela_4e_2012);
  bkg_zxss_mela_2e2m = (TH1F*) bkg_zxss_mela_2e2m_2011->Clone("bkg_zxss_mela_2e2m"); bkg_zxss_mela_2e2m->Add(bkg_zxss_mela_2e2m_2012);
  }

  TH1F* bkg_ggzz_mass = (TH1F*) bkg_ggzz_mass_2011->Clone("bkg_ggzz_mass"); bkg_ggzz_mass->Add(bkg_ggzz_mass_2012);
  TH1F* bkg_ggzz_mass_4m = (TH1F*) bkg_ggzz_mass_4m_2011->Clone("bkg_ggzz_mass_4m"); bkg_ggzz_mass_4m->Add(bkg_ggzz_mass_4m_2012);
  TH1F* bkg_ggzz_mass_4e = (TH1F*) bkg_ggzz_mass_4e_2011->Clone("bkg_ggzz_mass_4e"); bkg_ggzz_mass_4e->Add(bkg_ggzz_mass_4e_2012);
  TH1F* bkg_ggzz_mass_2e2m = (TH1F*) bkg_ggzz_mass_2e2m_2011->Clone("bkg_ggzz_mass_2e2m"); bkg_ggzz_mass_2e2m->Add(bkg_ggzz_mass_2e2m_2012);
  TH1F* bkg_ggzz_allzmass = (TH1F*) bkg_ggzz_allzmass_2011->Clone("bkg_ggzz_allzmass"); bkg_ggzz_allzmass->Add(bkg_ggzz_allzmass_2012);
  TH1F* bkg_ggzz_allzmass_2m = (TH1F*) bkg_ggzz_allzmass_2m_2011->Clone("bkg_ggzz_allzmass_2m"); bkg_ggzz_allzmass_2m->Add(bkg_ggzz_allzmass_2m_2012);
  TH1F* bkg_ggzz_allzmass_2e = (TH1F*) bkg_ggzz_allzmass_2e_2011->Clone("bkg_ggzz_allzmass_2e"); bkg_ggzz_allzmass_2e->Add(bkg_ggzz_allzmass_2e_2012);
  TH1F* bkg_ggzz_mela = (TH1F*) bkg_ggzz_mela_2011->Clone("bkg_ggzz_mela"); bkg_ggzz_mela->Add(bkg_ggzz_mela_2012);
  TH1F* bkg_ggzz_mela_4m = (TH1F*) bkg_ggzz_mela_4m_2011->Clone("bkg_ggzz_mela_4m"); bkg_ggzz_mela_4m->Add(bkg_ggzz_mela_4m_2012);
  TH1F* bkg_ggzz_mela_4e = (TH1F*) bkg_ggzz_mela_4e_2011->Clone("bkg_ggzz_mela_4e"); bkg_ggzz_mela_4e->Add(bkg_ggzz_mela_4e_2012);
  TH1F* bkg_ggzz_mela_2e2m = (TH1F*) bkg_ggzz_mela_2e2m_2011->Clone("bkg_ggzz_mela_2e2m"); bkg_ggzz_mela_2e2m->Add(bkg_ggzz_mela_2e2m_2012);

  TH1F* bkg_qqzz_mass = (TH1F*) bkg_qqzz_mass_2011->Clone("bkg_qqzz_mass"); bkg_qqzz_mass->Add(bkg_qqzz_mass_2012);
  TH1F* bkg_qqzz_mass_4m = (TH1F*) bkg_qqzz_mass_4m_2011->Clone("bkg_qqzz_mass_4m"); bkg_qqzz_mass_4m->Add(bkg_qqzz_mass_4m_2012);
  TH1F* bkg_qqzz_mass_4e = (TH1F*) bkg_qqzz_mass_4e_2011->Clone("bkg_qqzz_mass_4e"); bkg_qqzz_mass_4e->Add(bkg_qqzz_mass_4e_2012);
  TH1F* bkg_qqzz_mass_2e2m = (TH1F*) bkg_qqzz_mass_2e2m_2011->Clone("bkg_qqzz_mass_2e2m"); bkg_qqzz_mass_2e2m->Add(bkg_qqzz_mass_2e2m_2012);
  TH1F* bkg_qqzz_allzmass = (TH1F*) bkg_qqzz_allzmass_2011->Clone("bkg_qqzz_allzmass"); bkg_qqzz_allzmass->Add(bkg_qqzz_allzmass_2012);
  TH1F* bkg_qqzz_allzmass_2m = (TH1F*) bkg_qqzz_allzmass_2m_2011->Clone("bkg_qqzz_allzmass_2m"); bkg_qqzz_allzmass_2m->Add(bkg_qqzz_allzmass_2m_2012);
  TH1F* bkg_qqzz_allzmass_2e = (TH1F*) bkg_qqzz_allzmass_2e_2011->Clone("bkg_qqzz_allzmass_2e"); bkg_qqzz_allzmass_2e->Add(bkg_qqzz_allzmass_2e_2012);
  TH1F* bkg_qqzz_mela = (TH1F*) bkg_qqzz_mela_2011->Clone("bkg_qqzz_mela"); bkg_qqzz_mela->Add(bkg_qqzz_mela_2012);
  TH1F* bkg_qqzz_mela_4m = (TH1F*) bkg_qqzz_mela_4m_2011->Clone("bkg_qqzz_mela_4m"); bkg_qqzz_mela_4m->Add(bkg_qqzz_mela_4m_2012);
  TH1F* bkg_qqzz_mela_4e = (TH1F*) bkg_qqzz_mela_4e_2011->Clone("bkg_qqzz_mela_4e"); bkg_qqzz_mela_4e->Add(bkg_qqzz_mela_4e_2012);
  TH1F* bkg_qqzz_mela_2e2m = (TH1F*) bkg_qqzz_mela_2e2m_2011->Clone("bkg_qqzz_mela_2e2m"); bkg_qqzz_mela_2e2m->Add(bkg_qqzz_mela_2e2m_2012);

  // Using the luminosity parameters to provide an estimate for the 2011 signal mc
  TH1F* sig_ghzz_mass = (TH1F*) sig_ghzz_mass_2012->Clone("sig_ghzz_mass"); sig_ghzz_mass->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_ghzz_mass_4m = (TH1F*) sig_ghzz_mass_4m_2012->Clone("sig_ghzz_mass_4m"); sig_ghzz_mass_4m->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_ghzz_mass_4e = (TH1F*) sig_ghzz_mass_4e_2012->Clone("sig_ghzz_mass_4e"); sig_ghzz_mass_4e->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_ghzz_mass_2e2m = (TH1F*) sig_ghzz_mass_2e2m_2012->Clone("sig_ghzz_mass_2e2m"); sig_ghzz_mass_2e2m->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_ghzz_allzmass = (TH1F*) sig_ghzz_allzmass_2012->Clone("sig_ghzz_allzmass"); sig_ghzz_allzmass->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_ghzz_allzmass_2m = (TH1F*) sig_ghzz_allzmass_2m_2012->Clone("sig_ghzz_allzmass_2m"); sig_ghzz_allzmass_2m->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_ghzz_allzmass_2e = (TH1F*) sig_ghzz_allzmass_2e_2012->Clone("sig_ghzz_allzmass_2e"); sig_ghzz_allzmass_2e->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_ghzz_mela = (TH1F*) sig_ghzz_mela_2012->Clone("sig_ghzz_mela"); sig_ghzz_mela->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_ghzz_mela_4m = (TH1F*) sig_ghzz_mela_4m_2012->Clone("sig_ghzz_mela_4m"); sig_ghzz_mela_4m->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_ghzz_mela_4e = (TH1F*) sig_ghzz_mela_4e_2012->Clone("sig_ghzz_mela_4e"); sig_ghzz_mela_4e->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_ghzz_mela_2e2m = (TH1F*) sig_ghzz_mela_2e2m_2012->Clone("sig_ghzz_mela_2e2m"); sig_ghzz_mela_2e2m->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  
  TH1F* sig_qhzz_mass = (TH1F*) sig_qhzz_mass_2012->Clone("sig_qhzz_mass"); sig_qhzz_mass->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_qhzz_mass_4m = (TH1F*) sig_qhzz_mass_4m_2012->Clone("sig_qhzz_mass_4m"); sig_qhzz_mass_4m->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_qhzz_mass_4e = (TH1F*) sig_qhzz_mass_4e_2012->Clone("sig_qhzz_mass_4e"); sig_qhzz_mass_4e->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_qhzz_mass_2e2m = (TH1F*) sig_qhzz_mass_2e2m_2012->Clone("sig_qhzz_mass_2e2m"); sig_qhzz_mass_2e2m->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_qhzz_allzmass = (TH1F*) sig_qhzz_allzmass_2012->Clone("sig_qhzz_allzmass"); sig_qhzz_allzmass->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_qhzz_allzmass_2m = (TH1F*) sig_qhzz_allzmass_2m_2012->Clone("sig_qhzz_allzmass_2m"); sig_qhzz_allzmass_2m->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_qhzz_allzmass_2e = (TH1F*) sig_qhzz_allzmass_2e_2012->Clone("sig_qhzz_allzmass_2e"); sig_qhzz_allzmass_2e->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_qhzz_mela = (TH1F*) sig_qhzz_mela_2012->Clone("sig_qhzz_mela"); sig_qhzz_mela->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_qhzz_mela_4m = (TH1F*) sig_qhzz_mela_4m_2012->Clone("sig_qhzz_mela_4m"); sig_qhzz_mela_4m->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_qhzz_mela_4e = (TH1F*) sig_qhzz_mela_4e_2012->Clone("sig_qhzz_mela_4e"); sig_qhzz_mela_4e->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH1F* sig_qhzz_mela_2e2m = (TH1F*) sig_qhzz_mela_2e2m_2012->Clone("sig_qhzz_mela_2e2m"); sig_qhzz_mela_2e2m->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  
  // Saving all histograms individually to separate files
  TCanvas *cv = 0;
  TLegend *legend = 0;

  cv = new TCanvas("cv", "Canvas", 800, 600);
//   data_mass->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_data_mass.gif").c_str());
//   data_mass_4m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_data_mass_4m.gif").c_str());
//   data_mass_4e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_data_mass_4e.gif").c_str());
//   data_mass_2e2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_data_mass_2e2m.gif").c_str());
//   data_allzmass->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_data_allzmass.gif").c_str());
//   data_allzmass_2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_data_allzmass_2m.gif").c_str());
//   data_allzmass_2e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_data_allzmass_2e.gif").c_str());
//   data_mela->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_data_mela.gif").c_str());
//   data_mela_4m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_data_mela_4m.gif").c_str());
//   data_mela_4e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_data_mela_4e.gif").c_str());
//   data_mela_2e2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_data_mela_2e2m.gif").c_str());

//   if (plotZXbackground) {
//   bkg_zxss_mass->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_zxss_mass.gif").c_str());
//   bkg_zxss_mass_4m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_zxss_mass_4m.gif").c_str());
//   bkg_zxss_mass_4e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_zxss_mass_4e.gif").c_str());
//   bkg_zxss_mass_2e2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_zxss_mass_2e2m.gif").c_str());
//   bkg_zxss_allzmass->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_zxss_allzmass.gif").c_str());
//   bkg_zxss_allzmass_2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_zxss_allzmass_2m.gif").c_str());
//   bkg_zxss_allzmass_2e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_zxss_allzmass_2e.gif").c_str());
//   bkg_zxss_mela->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_zxss_mela.gif").c_str());
//   bkg_zxss_mela_4m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_zxss_mela_4m.gif").c_str());
//   bkg_zxss_mela_4e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_zxss_mela_4e.gif").c_str());
//   bkg_zxss_mela_2e2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_zxss_mela_2e2m.gif").c_str());
//   }

//   bkg_ggzz_mass->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_ggzz_mass.gif").c_str());
//   bkg_ggzz_mass_4m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_ggzz_mass_4m.gif").c_str());
//   bkg_ggzz_mass_4e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_ggzz_mass_4e.gif").c_str());
//   bkg_ggzz_mass_2e2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_ggzz_mass_2e2m.gif").c_str());
//   bkg_ggzz_allzmass->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_ggzz_allzmass.gif").c_str());
//   bkg_ggzz_allzmass_2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_ggzz_allzmass_2m.gif").c_str());
//   bkg_ggzz_allzmass_2e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_ggzz_allzmass_2e.gif").c_str());
//   bkg_ggzz_mela->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_ggzz_mela.gif").c_str());
//   bkg_ggzz_mela_4m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_ggzz_mela_4m.gif").c_str());
//   bkg_ggzz_mela_4e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_ggzz_mela_4e.gif").c_str());
//   bkg_ggzz_mela_2e2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_ggzz_mela_2e2m.gif").c_str());

//   bkg_qqzz_mass->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_qqzz_mass.gif").c_str());
//   bkg_qqzz_mass_4m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_qqzz_mass_4m.gif").c_str());
//   bkg_qqzz_mass_4e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_qqzz_mass_4e.gif").c_str());
//   bkg_qqzz_mass_2e2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_qqzz_mass_2e2m.gif").c_str());
//   bkg_qqzz_allzmass->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_qqzz_allzmass.gif").c_str());
//   bkg_qqzz_allzmass_2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_qqzz_allzmass_2m.gif").c_str());
//   bkg_qqzz_allzmass_2e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_qqzz_allzmass_2e.gif").c_str());
//   bkg_qqzz_mela->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_qqzz_mela.gif").c_str());
//   bkg_qqzz_mela_4m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_qqzz_mela_4m.gif").c_str());
//   bkg_qqzz_mela_4e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_qqzz_mela_4e.gif").c_str());
//   bkg_qqzz_mela_2e2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_bkg_qqzz_mela_2e2m.gif").c_str());

//   sig_qhzz_mass->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_qhzz_mass.gif").c_str());
//   sig_qhzz_mass_4m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_qhzz_mass_4m.gif").c_str());
//   sig_qhzz_mass_4e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_qhzz_mass_4e.gif").c_str());
//   sig_qhzz_mass_2e2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_qhzz_mass_2e2m.gif").c_str());
//   sig_qhzz_allzmass->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_qhzz_allzmass.gif").c_str());
//   sig_qhzz_allzmass_2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_qhzz_allzmass_2m.gif").c_str());
//   sig_qhzz_allzmass_2e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_qhzz_allzmass_2e.gif").c_str());
//   sig_qhzz_mela->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_qhzz_mela.gif").c_str());
//   sig_qhzz_mela_4m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_qhzz_mela_4m.gif").c_str());
//   sig_qhzz_mela_4e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_qhzz_mela_4e.gif").c_str());
//   sig_qhzz_mela_2e2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_qhzz_mela_2e2m.gif").c_str());

//   sig_ghzz_mass->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_ghzz_mass.gif").c_str());
//   sig_ghzz_mass_4m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_ghzz_mass_4m.gif").c_str());
//   sig_ghzz_mass_4e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_ghzz_mass_4e.gif").c_str());
//   sig_ghzz_mass_2e2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_ghzz_mass_2e2m.gif").c_str());
//   sig_ghzz_allzmass->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_ghzz_allzmass.gif").c_str());
//   sig_ghzz_allzmass_2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_ghzz_allzmass_2m.gif").c_str());
//   sig_ghzz_allzmass_2e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_ghzz_allzmass_2e.gif").c_str());
//   sig_ghzz_mela->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_ghzz_mela.gif").c_str());
//   sig_ghzz_mela_4m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_ghzz_mela_4m.gif").c_str());
//   sig_ghzz_mela_4e->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_ghzz_mela_4e.gif").c_str());
//   sig_ghzz_mela_2e2m->Draw("hist"); cv->SaveAs((outDirectory + "doHZZAnalysis_sig_ghzz_mela_2e2m.gif").c_str());

  // Now making stacks to plot everything on the same canvas
  // Background
  if (plotZXbackground) {
  bkg_zxss_mass->SetFillColor(12);
  bkg_zxss_mass_4m->SetFillColor(12);
  bkg_zxss_mass_4e->SetFillColor(12);
  bkg_zxss_mass_2e2m->SetFillColor(12);
  bkg_zxss_allzmass->SetFillColor(12);
  bkg_zxss_allzmass_2m->SetFillColor(12);
  bkg_zxss_allzmass_2e->SetFillColor(12);
  bkg_zxss_mela->SetFillColor(12);
  bkg_zxss_mela_4m->SetFillColor(12);
  bkg_zxss_mela_4e->SetFillColor(12);
  bkg_zxss_mela_2e2m->SetFillColor(12);
  }

  bkg_ggzz_mass->SetFillColor(14);
  bkg_ggzz_mass_4m->SetFillColor(14);
  bkg_ggzz_mass_4e->SetFillColor(14);
  bkg_ggzz_mass_2e2m->SetFillColor(14);
  bkg_ggzz_allzmass->SetFillColor(14);
  bkg_ggzz_allzmass_2m->SetFillColor(14);
  bkg_ggzz_allzmass_2e->SetFillColor(14);
  bkg_ggzz_mela->SetFillColor(14);
  bkg_ggzz_mela_4m->SetFillColor(14);
  bkg_ggzz_mela_4e->SetFillColor(14);
  bkg_ggzz_mela_2e2m->SetFillColor(14);

  bkg_qqzz_mass->SetFillColor(16);
  bkg_qqzz_mass_4m->SetFillColor(16);
  bkg_qqzz_mass_4e->SetFillColor(16);
  bkg_qqzz_mass_2e2m->SetFillColor(16);
  bkg_qqzz_allzmass->SetFillColor(16);
  bkg_qqzz_allzmass_2m->SetFillColor(16);
  bkg_qqzz_allzmass_2e->SetFillColor(16);
  bkg_qqzz_mela->SetFillColor(16);
  bkg_qqzz_mela_4m->SetFillColor(16);
  bkg_qqzz_mela_4e->SetFillColor(16);
  bkg_qqzz_mela_2e2m->SetFillColor(16);

  THStack stack_mass("stack_mass", "Mass");
  if (plotZXbackground) stack_mass.Add(bkg_zxss_mass);
  stack_mass.Add(bkg_qqzz_mass);
  stack_mass.Add(bkg_ggzz_mass);

  THStack stack_mass_4m("stack_mass_4m", "Mass, 4m");
  if (plotZXbackground) stack_mass_4m.Add(bkg_zxss_mass_4m);
  stack_mass_4m.Add(bkg_qqzz_mass_4m);
  stack_mass_4m.Add(bkg_ggzz_mass_4m);

  THStack stack_mass_4e("stack_mass_4e", "Mass, 4e");
  if (plotZXbackground) stack_mass_4e.Add(bkg_zxss_mass_4e);
  stack_mass_4e.Add(bkg_qqzz_mass_4e);
  stack_mass_4e.Add(bkg_ggzz_mass_4e);

  THStack stack_mass_2e2m("stack_mass_2e2m", "Mass, 2e2m");
  if (plotZXbackground) stack_mass_2e2m.Add(bkg_zxss_mass_2e2m);
  stack_mass_2e2m.Add(bkg_qqzz_mass_2e2m);
  stack_mass_2e2m.Add(bkg_ggzz_mass_2e2m);

  THStack stack_allzmass("stack_allzmass", "Z candidate masses");
  if (plotZXbackground) stack_allzmass.Add(bkg_zxss_allzmass);
  stack_allzmass.Add(bkg_qqzz_allzmass);
  stack_allzmass.Add(bkg_ggzz_allzmass);

  THStack stack_allzmass_2m("stack_allzmass_2m", "Z candidate masses, 2m");
  if (plotZXbackground) stack_allzmass_2m.Add(bkg_zxss_allzmass_2m);
  stack_allzmass_2m.Add(bkg_qqzz_allzmass_2m);
  stack_allzmass_2m.Add(bkg_ggzz_allzmass_2m);

  THStack stack_allzmass_2e("stack_allzmass_2e", "Z candidate masses, 2e");
  if (plotZXbackground) stack_allzmass_2e.Add(bkg_zxss_allzmass_2e);
  stack_allzmass_2e.Add(bkg_qqzz_allzmass_2e);
  stack_allzmass_2e.Add(bkg_ggzz_allzmass_2e);

  THStack stack_mela("stack_mela", "MELA");
  if (plotZXbackground) stack_mela.Add(bkg_zxss_mela);
  stack_mela.Add(bkg_qqzz_mela);
  stack_mela.Add(bkg_ggzz_mela);

  THStack stack_mela_4m("stack_mela_4m", "MELA, 4m");
  if (plotZXbackground) stack_mela_4m.Add(bkg_zxss_mela_4m);
  stack_mela_4m.Add(bkg_qqzz_mela_4m);
  stack_mela_4m.Add(bkg_ggzz_mela_4m);

  THStack stack_mela_4e("stack_mela_4e", "MELA, 4e");
  if (plotZXbackground) stack_mela_4e.Add(bkg_zxss_mela_4e);
  stack_mela_4e.Add(bkg_qqzz_mela_4e);
  stack_mela_4e.Add(bkg_ggzz_mela_4e);

  THStack stack_mela_2e2m("stack_mela_2e2m", "MELA, 2e2m");
  if (plotZXbackground) stack_mela_2e2m.Add(bkg_zxss_mela_2e2m);
  stack_mela_2e2m.Add(bkg_qqzz_mela_2e2m);
  stack_mela_2e2m.Add(bkg_ggzz_mela_2e2m);

  // Signal
  sig_qhzz_mass->SetFillColor(46);
  sig_qhzz_mass_4m->SetFillColor(46);
  sig_qhzz_mass_4e->SetFillColor(46);
  sig_qhzz_mass_2e2m->SetFillColor(46);
  sig_qhzz_allzmass->SetFillColor(46);
  sig_qhzz_allzmass_2m->SetFillColor(46);
  sig_qhzz_allzmass_2e->SetFillColor(46);
  sig_qhzz_mela->SetFillColor(46);
  sig_qhzz_mela_4m->SetFillColor(46);
  sig_qhzz_mela_4e->SetFillColor(46);
  sig_qhzz_mela_2e2m->SetFillColor(46);

  sig_ghzz_mass->SetFillColor(44);
  sig_ghzz_mass_4m->SetFillColor(44);
  sig_ghzz_mass_4e->SetFillColor(44);
  sig_ghzz_mass_2e2m->SetFillColor(44);
  sig_ghzz_allzmass->SetFillColor(44);
  sig_ghzz_allzmass_2m->SetFillColor(44);
  sig_ghzz_allzmass_2e->SetFillColor(44);
  sig_ghzz_mela->SetFillColor(44);
  sig_ghzz_mela_4m->SetFillColor(44);
  sig_ghzz_mela_4e->SetFillColor(44);
  sig_ghzz_mela_2e2m->SetFillColor(44);

  // Adding signal to stacks
  stack_mass.Add(sig_ghzz_mass); stack_mass.Add(sig_qhzz_mass);
  stack_mass_4m.Add(sig_ghzz_mass_4m); stack_mass_4m.Add(sig_qhzz_mass_4m);
  stack_mass_4e.Add(sig_ghzz_mass_4e); stack_mass_4e.Add(sig_qhzz_mass_4e);
  stack_mass_2e2m.Add(sig_ghzz_mass_2e2m); stack_mass_2e2m.Add(sig_qhzz_mass_2e2m);
  stack_allzmass.Add(sig_ghzz_allzmass); stack_allzmass.Add(sig_qhzz_allzmass);
  stack_allzmass_2m.Add(sig_ghzz_allzmass_2m); stack_allzmass_2m.Add(sig_qhzz_allzmass_2m);
  stack_allzmass_2e.Add(sig_ghzz_allzmass_2e); stack_allzmass_2e.Add(sig_qhzz_allzmass_2e);
  stack_mela.Add(sig_ghzz_mela); stack_mela.Add(sig_qhzz_mela);
  stack_mela_4m.Add(sig_ghzz_mela_4m); stack_mela_4m.Add(sig_qhzz_mela_4m);
  stack_mela_4e.Add(sig_ghzz_mela_4e); stack_mela_4e.Add(sig_qhzz_mela_4e);
  stack_mela_2e2m.Add(sig_ghzz_mela_2e2m); stack_mela_2e2m.Add(sig_qhzz_mela_2e2m);

  // Graphic options for data
  data_mass->SetLineWidth(2);
  data_mass_4m->SetLineWidth(2);
  data_mass_4e->SetLineWidth(2);
  data_mass_2e2m->SetLineWidth(2);
  data_allzmass->SetLineWidth(2);
  data_allzmass_2m->SetLineWidth(2);
  data_allzmass_2e->SetLineWidth(2);
  data_mela->SetLineWidth(2);
  data_mela_4m->SetLineWidth(2);
  data_mela_4e->SetLineWidth(2);
  data_mela_2e2m->SetLineWidth(2);

  // Plotting and saving to files
  cv = new TCanvas("cv", "Canvas", 800, 600);
  stack_mass.Draw();
  stack_mass.SetMaximum(TMath::Max(data_mass->GetMaximum(), stack_mass.GetMaximum()) * 1.5);
  cv->Update();
  data_mass->Draw("e1same");
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectory + "doHZZHistograms_mass.gif").c_str());

  cv = new TCanvas("cv", "Canvas", 800, 600);
  stack_mass_4m.Draw();
  stack_mass_4m.SetMaximum(TMath::Max(data_mass_4m->GetMaximum(), stack_mass_4m.GetMaximum()) * 1.5);
  cv->Update();
  data_mass_4m->Draw("e1same");
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectory + "doHZZHistograms_mass_4m.gif").c_str());

  cv = new TCanvas("cv", "Canvas", 800, 600);
  stack_mass_4e.Draw();
  stack_mass_4e.SetMaximum(TMath::Max(data_mass_4e->GetMaximum(), stack_mass_4e.GetMaximum()) * 1.5);
  cv->Update();
  data_mass_4e->Draw("e1same");
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectory + "doHZZHistograms_mass_4e.gif").c_str());

  cv = new TCanvas("cv", "Canvas", 800, 600);
  stack_mass_2e2m.Draw();
  stack_mass_2e2m.SetMaximum(TMath::Max(data_mass_2e2m->GetMaximum(), stack_mass_2e2m.GetMaximum()) * 1.5);
  cv->Update();
  data_mass_2e2m->Draw("e1same");
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectory + "doHZZHistograms_mass_2e2m.gif").c_str());

  cv = new TCanvas("cv", "Canvas", 800, 600);
  stack_allzmass.Draw();
  stack_allzmass.SetMaximum(TMath::Max(data_allzmass->GetMaximum(), stack_allzmass.GetMaximum()) * 1.5);
  cv->Update();
  data_allzmass->Draw("e1same");
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectory + "doHZZHistograms_allzmass.gif").c_str());

  cv = new TCanvas("cv", "Canvas", 800, 600);
  stack_allzmass_2e.Draw();
  stack_allzmass_2e.SetMaximum(TMath::Max(data_allzmass_2e->GetMaximum(), stack_allzmass_2e.GetMaximum()) * 1.5);
  cv->Update();
  data_allzmass_2e->Draw("e1same");
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectory + "doHZZHistograms_allzmass_2e.gif").c_str());

  cv = new TCanvas("cv", "Canvas", 800, 600);
  stack_allzmass_2m.Draw();
  stack_allzmass_2m.SetMaximum(TMath::Max(data_allzmass_2m->GetMaximum(), stack_allzmass_2m.GetMaximum()) * 1.5);
  cv->Update();
  data_allzmass_2m->Draw("e1same");
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectory + "doHZZHistograms_allzmass_2m.gif").c_str());

  cv = new TCanvas("cv", "Canvas", 800, 600);
  stack_mela.Draw();
  stack_mela.SetMaximum(TMath::Max(data_mela->GetMaximum(), stack_mela.GetMaximum()) * 1.5);
  cv->Update();
  data_mela->Draw("e1same");
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectory + "doHZZHistograms_mela.gif").c_str());

  cv = new TCanvas("cv", "Canvas", 800, 600);
  stack_mela_4m.Draw();
  stack_mela_4m.SetMaximum(TMath::Max(data_mela_4m->GetMaximum(), stack_mela_4m.GetMaximum()) * 1.5);
  cv->Update();
  data_mela_4m->Draw("e1same");
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectory + "doHZZHistograms_mela_4m.gif").c_str());

  cv = new TCanvas("cv", "Canvas", 800, 600);
  stack_mela_4e.Draw();
  stack_mela_4e.SetMaximum(TMath::Max(data_mela_4e->GetMaximum(), stack_mela_4e.GetMaximum()) * 1.5);
  cv->Update();
  data_mela_4e->Draw("e1same");
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectory + "doHZZHistograms_mela_4e.gif").c_str());
  
  cv = new TCanvas("cv", "Canvas", 800, 600);
  stack_mela_2e2m.Draw();
  stack_mela_2e2m.SetMaximum(TMath::Max(data_mela_2e2m->GetMaximum(), stack_mela_2e2m.GetMaximum()) * 1.5);
  cv->Update();
  data_mela_2e2m->Draw("e1same");
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectory + "doHZZHistograms_mela_2e2m.gif").c_str());
  
  // Now plotting with different colors for the different channels
  data_mass_4m->SetLineColor(3); data_mass_4e->SetLineColor(4); data_mass_2e2m->SetLineColor(6);
  data_allzmass_2m->SetLineColor(3); data_allzmass_2e->SetLineColor(4);
  data_mela_4m->SetLineColor(3); data_mela_4e->SetLineColor(4); data_mela_2e2m->SetLineColor(6);

  // Making stacks of data 
  THStack stack_data_mass("stack_data_mass", "Mass");
  stack_data_mass.Add(data_mass_4m);
  stack_data_mass.Add(data_mass_4e);
  stack_data_mass.Add(data_mass_2e2m);

  THStack stack_data_allzmass("stack_data_allzmass", "Z candidate masses");
  stack_data_allzmass.Add(data_allzmass_2m);
  stack_data_allzmass.Add(data_allzmass_2e);

  THStack stack_data_mela("stack_data_mela", "MELA");
  stack_data_mela.Add(data_mela_4m);
  stack_data_mela.Add(data_mela_4e);
  stack_data_mela.Add(data_mela_2e2m);

  cv = new TCanvas("cv", "Canvas", 800, 600);
  stack_mass.Draw();
  stack_data_mass.Draw("e1same");
  stack_mass.SetMaximum(TMath::Max(stack_data_mass.GetMaximum(), stack_mass.GetMaximum()) * 1.5);
  cv->Update();
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectory + "doHZZHistograms_mass_color.gif").c_str());

  cv = new TCanvas("cv", "Canvas", 800, 600);
  stack_allzmass.Draw();
  stack_data_allzmass.Draw("e1same");
  stack_allzmass.SetMaximum(TMath::Max(stack_data_allzmass.GetMaximum(), stack_allzmass.GetMaximum()) * 1.5);
  cv->Update();
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectory + "doHZZHistograms_allzmass_color.gif").c_str());

  cv = new TCanvas("cv", "Canvas", 800, 600);
  stack_mela.Draw();
  stack_data_mela.Draw("e1same");
  stack_mela.SetMaximum(TMath::Max(stack_data_mela.GetMaximum(), stack_mela.GetMaximum()) * 1.5);
  cv->Update();
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectory + "doHZZHistograms_mela_color.gif").c_str());
}
