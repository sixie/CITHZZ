#include <sstream>
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TMath.h"


void doHZZHistograms_plot2D (char* filename_2011, char* filename_2012, char* outDirectory) {
  /*
     Uses the histograms stored in file_2011 and file_2012, adding them up and plotting (background/signal/data),
     saving the created canvases into files in outDirectory.
   */

  // Getting rid of stats box
  gStyle->SetOptStat(0);

  // Luminosity parameters in fb-1 (for signal intensity estimation)
  Float_t lumi_2011 = 5.05;
  Float_t lumi_2012 = 5.26;

  // Putting outDirectory into string for convenience
  std::string outDirectoryString(outDirectory);

  // Retrieving histograms from file
  TFile* file_2011 = new TFile(filename_2011);
  TFile* file_2012 = new TFile(filename_2012);

  TH2F* data_2D_2011 = (TH2F*) file_2011->Get("data_hist_2D");
  TH2F* data_2D_4m_2011 = (TH2F*) file_2011->Get("data_hist_2D_4m");
  TH2F* data_2D_4e_2011 = (TH2F*) file_2011->Get("data_hist_2D_4e");
  TH2F* data_2D_2e2m_2011 = (TH2F*) file_2011->Get("data_hist_2D_2e2m");
  TH2F* data_2D_2m2e_2011 = (TH2F*) file_2011->Get("data_hist_2D_2m2e");
  TH2F* data_2D_2012 = (TH2F*) file_2012->Get("data_hist_2D");
  TH2F* data_2D_4m_2012 = (TH2F*) file_2012->Get("data_hist_2D_4m");
  TH2F* data_2D_4e_2012 = (TH2F*) file_2012->Get("data_hist_2D_4e");
  TH2F* data_2D_2e2m_2012 = (TH2F*) file_2012->Get("data_hist_2D_2e2m");
  TH2F* data_2D_2m2e_2012 = (TH2F*) file_2012->Get("data_hist_2D_2m2e");

  TH2F* bkg_zxss_2D_2011 = (TH2F*) file_2011->Get("bkg_zxss_hist_2D");
  TH2F* bkg_zxss_2D_4m_2011 = (TH2F*) file_2011->Get("bkg_zxss_hist_2D_4m");
  TH2F* bkg_zxss_2D_4e_2011 = (TH2F*) file_2011->Get("bkg_zxss_hist_2D_4e");
  TH2F* bkg_zxss_2D_2e2m_2011 = (TH2F*) file_2011->Get("bkg_zxss_hist_2D_2e2m");
  TH2F* bkg_zxss_2D_2m2e_2011 = (TH2F*) file_2011->Get("bkg_zxss_hist_2D_2m2e");
  TH2F* bkg_zxss_2D_2012 = (TH2F*) file_2012->Get("bkg_zxss_hist_2D");
  TH2F* bkg_zxss_2D_4m_2012 = (TH2F*) file_2012->Get("bkg_zxss_hist_2D_4m");
  TH2F* bkg_zxss_2D_4e_2012 = (TH2F*) file_2012->Get("bkg_zxss_hist_2D_4e");
  TH2F* bkg_zxss_2D_2e2m_2012 = (TH2F*) file_2012->Get("bkg_zxss_hist_2D_2e2m");
  TH2F* bkg_zxss_2D_2m2e_2012 = (TH2F*) file_2012->Get("bkg_zxss_hist_2D_2m2e");

  TH2F* bkg_ggzz_2D_2011 = (TH2F*) file_2011->Get("bkg_ggzz_hist_2D");
  TH2F* bkg_ggzz_2D_4m_2011 = (TH2F*) file_2011->Get("bkg_ggzz_hist_2D_4m");
  TH2F* bkg_ggzz_2D_4e_2011 = (TH2F*) file_2011->Get("bkg_ggzz_hist_2D_4e");
  TH2F* bkg_ggzz_2D_2e2m_2011 = (TH2F*) file_2011->Get("bkg_ggzz_hist_2D_2e2m");
  TH2F* bkg_ggzz_2D_2m2e_2011 = (TH2F*) file_2011->Get("bkg_ggzz_hist_2D_2m2e");
  TH2F* bkg_ggzz_2D_2012 = (TH2F*) file_2012->Get("bkg_ggzz_hist_2D");
  TH2F* bkg_ggzz_2D_4m_2012 = (TH2F*) file_2012->Get("bkg_ggzz_hist_2D_4m");
  TH2F* bkg_ggzz_2D_4e_2012 = (TH2F*) file_2012->Get("bkg_ggzz_hist_2D_4e");
  TH2F* bkg_ggzz_2D_2e2m_2012 = (TH2F*) file_2012->Get("bkg_ggzz_hist_2D_2e2m");
  TH2F* bkg_ggzz_2D_2m2e_2012 = (TH2F*) file_2012->Get("bkg_ggzz_hist_2D_2m2e");

  TH2F* bkg_qqzz_2D_2011 = (TH2F*) file_2011->Get("bkg_qqzz_hist_2D");
  TH2F* bkg_qqzz_2D_4m_2011 = (TH2F*) file_2011->Get("bkg_qqzz_hist_2D_4m");
  TH2F* bkg_qqzz_2D_4e_2011 = (TH2F*) file_2011->Get("bkg_qqzz_hist_2D_4e");
  TH2F* bkg_qqzz_2D_2e2m_2011 = (TH2F*) file_2011->Get("bkg_qqzz_hist_2D_2e2m");
  TH2F* bkg_qqzz_2D_2m2e_2011 = (TH2F*) file_2011->Get("bkg_qqzz_hist_2D_2m2e");
  TH2F* bkg_qqzz_2D_2012 = (TH2F*) file_2012->Get("bkg_qqzz_hist_2D");
  TH2F* bkg_qqzz_2D_4m_2012 = (TH2F*) file_2012->Get("bkg_qqzz_hist_2D_4m");
  TH2F* bkg_qqzz_2D_4e_2012 = (TH2F*) file_2012->Get("bkg_qqzz_hist_2D_4e");
  TH2F* bkg_qqzz_2D_2e2m_2012 = (TH2F*) file_2012->Get("bkg_qqzz_hist_2D_2e2m");
  TH2F* bkg_qqzz_2D_2m2e_2012 = (TH2F*) file_2012->Get("bkg_qqzz_hist_2D_2m2e");

  TH2F* sig_qhzz_2D_2011 = (TH2F*) file_2011->Get("sig_qhzz_hist_2D");
  TH2F* sig_qhzz_2D_4m_2011 = (TH2F*) file_2011->Get("sig_qhzz_hist_2D_4m");
  TH2F* sig_qhzz_2D_4e_2011 = (TH2F*) file_2011->Get("sig_qhzz_hist_2D_4e");
  TH2F* sig_qhzz_2D_2e2m_2011 = (TH2F*) file_2011->Get("sig_qhzz_hist_2D_2e2m");
  TH2F* sig_qhzz_2D_2m2e_2011 = (TH2F*) file_2011->Get("sig_qhzz_hist_2D_2m2e");
  TH2F* sig_qhzz_2D_2012 = (TH2F*) file_2012->Get("sig_qhzz_hist_2D");
  TH2F* sig_qhzz_2D_4m_2012 = (TH2F*) file_2012->Get("sig_qhzz_hist_2D_4m");
  TH2F* sig_qhzz_2D_4e_2012 = (TH2F*) file_2012->Get("sig_qhzz_hist_2D_4e");
  TH2F* sig_qhzz_2D_2e2m_2012 = (TH2F*) file_2012->Get("sig_qhzz_hist_2D_2e2m");
  TH2F* sig_qhzz_2D_2m2e_2012 = (TH2F*) file_2012->Get("sig_qhzz_hist_2D_2m2e");

  TH2F* sig_ghzz_2D_2011 = (TH2F*) file_2011->Get("sig_ghzz_hist_2D");
  TH2F* sig_ghzz_2D_4m_2011 = (TH2F*) file_2011->Get("sig_ghzz_hist_2D_4m");
  TH2F* sig_ghzz_2D_4e_2011 = (TH2F*) file_2011->Get("sig_ghzz_hist_2D_4e");
  TH2F* sig_ghzz_2D_2e2m_2011 = (TH2F*) file_2011->Get("sig_ghzz_hist_2D_2e2m");
  TH2F* sig_ghzz_2D_2m2e_2011 = (TH2F*) file_2011->Get("sig_ghzz_hist_2D_2m2e");
  TH2F* sig_ghzz_2D_2012 = (TH2F*) file_2012->Get("sig_ghzz_hist_2D");
  TH2F* sig_ghzz_2D_4m_2012 = (TH2F*) file_2012->Get("sig_ghzz_hist_2D_4m");
  TH2F* sig_ghzz_2D_4e_2012 = (TH2F*) file_2012->Get("sig_ghzz_hist_2D_4e");
  TH2F* sig_ghzz_2D_2e2m_2012 = (TH2F*) file_2012->Get("sig_ghzz_hist_2D_2e2m");
  TH2F* sig_ghzz_2D_2m2e_2012 = (TH2F*) file_2012->Get("sig_ghzz_hist_2D_2m2e");

  // Now adding up all 2011 and 2012 data
  TH2F* data_2D = (TH2F*) data_2D_2011->Clone("data_2D"); data_2D->Add(data_2D_2012);
  TH2F* data_2D_4m = (TH2F*) data_2D_4m_2011->Clone("data_2D_4m"); data_2D_4m->Add(data_2D_4m_2012);
  TH2F* data_2D_4e = (TH2F*) data_2D_4e_2011->Clone("data_2D_4e"); data_2D_4e->Add(data_2D_4e_2012);
  TH2F* data_2D_2e2m = (TH2F*) data_2D_2e2m_2011->Clone("data_2D_2e2m"); data_2D_2e2m->Add(data_2D_2e2m_2012);
  TH2F* data_2D_2m2e = (TH2F*) data_2D_2m2e_2011->Clone("data_2D_2m2e"); data_2D_2m2e->Add(data_2D_2m2e_2012);

  TH2F* bkg_zxss_2D = (TH2F*) bkg_zxss_2D_2011->Clone("bkg_zxss_2D"); bkg_zxss_2D->Add(bkg_zxss_2D_2012);
  TH2F* bkg_zxss_2D_4m = (TH2F*) bkg_zxss_2D_4m_2011->Clone("bkg_zxss_2D_4m"); bkg_zxss_2D_4m->Add(bkg_zxss_2D_4m_2012);
  TH2F* bkg_zxss_2D_4e = (TH2F*) bkg_zxss_2D_4e_2011->Clone("bkg_zxss_2D_4e"); bkg_zxss_2D_4e->Add(bkg_zxss_2D_4e_2012);
  TH2F* bkg_zxss_2D_2e2m = (TH2F*) bkg_zxss_2D_2e2m_2011->Clone("bkg_zxss_2D_2e2m"); bkg_zxss_2D_2e2m->Add(bkg_zxss_2D_2e2m_2012);
  TH2F* bkg_zxss_2D_2m2e = (TH2F*) bkg_zxss_2D_2m2e_2011->Clone("bkg_zxss_2D_2m2e"); bkg_zxss_2D_2m2e->Add(bkg_zxss_2D_2m2e_2012);

  TH2F* bkg_ggzz_2D = (TH2F*) bkg_ggzz_2D_2011->Clone("bkg_ggzz_2D"); bkg_ggzz_2D->Add(bkg_ggzz_2D_2012);
  TH2F* bkg_ggzz_2D_4m = (TH2F*) bkg_ggzz_2D_4m_2011->Clone("bkg_ggzz_2D_4m"); bkg_ggzz_2D_4m->Add(bkg_ggzz_2D_4m_2012);
  TH2F* bkg_ggzz_2D_4e = (TH2F*) bkg_ggzz_2D_4e_2011->Clone("bkg_ggzz_2D_4e"); bkg_ggzz_2D_4e->Add(bkg_ggzz_2D_4e_2012);
  TH2F* bkg_ggzz_2D_2e2m = (TH2F*) bkg_ggzz_2D_2e2m_2011->Clone("bkg_ggzz_2D_2e2m"); bkg_ggzz_2D_2e2m->Add(bkg_ggzz_2D_2e2m_2012);
  TH2F* bkg_ggzz_2D_2m2e = (TH2F*) bkg_ggzz_2D_2m2e_2011->Clone("bkg_ggzz_2D_2m2e"); bkg_ggzz_2D_2m2e->Add(bkg_ggzz_2D_2m2e_2012);

  TH2F* bkg_qqzz_2D = (TH2F*) bkg_qqzz_2D_2011->Clone("bkg_qqzz_2D"); bkg_qqzz_2D->Add(bkg_qqzz_2D_2012);
  TH2F* bkg_qqzz_2D_4m = (TH2F*) bkg_qqzz_2D_4m_2011->Clone("bkg_qqzz_2D_4m"); bkg_qqzz_2D_4m->Add(bkg_qqzz_2D_4m_2012);
  TH2F* bkg_qqzz_2D_4e = (TH2F*) bkg_qqzz_2D_4e_2011->Clone("bkg_qqzz_2D_4e"); bkg_qqzz_2D_4e->Add(bkg_qqzz_2D_4e_2012);
  TH2F* bkg_qqzz_2D_2e2m = (TH2F*) bkg_qqzz_2D_2e2m_2011->Clone("bkg_qqzz_2D_2e2m"); bkg_qqzz_2D_2e2m->Add(bkg_qqzz_2D_2e2m_2012);
  TH2F* bkg_qqzz_2D_2m2e = (TH2F*) bkg_qqzz_2D_2m2e_2011->Clone("bkg_qqzz_2D_2m2e"); bkg_qqzz_2D_2m2e->Add(bkg_qqzz_2D_2m2e_2012);

  // Using the luminosity parameters to provide an estimate for the 2011 signal mc
  TH2F* sig_ghzz_2D = (TH2F*) sig_ghzz_2D_2012->Clone("sig_ghzz_2D"); sig_ghzz_2D->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH2F* sig_ghzz_2D_4m = (TH2F*) sig_ghzz_2D_4m_2012->Clone("sig_ghzz_2D_4m"); sig_ghzz_2D_4m->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH2F* sig_ghzz_2D_4e = (TH2F*) sig_ghzz_2D_4e_2012->Clone("sig_ghzz_2D_4e"); sig_ghzz_2D_4e->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH2F* sig_ghzz_2D_2e2m = (TH2F*) sig_ghzz_2D_2e2m_2012->Clone("sig_ghzz_2D_2e2m"); sig_ghzz_2D_2e2m->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH2F* sig_ghzz_2D_2m2e = (TH2F*) sig_ghzz_2D_2m2e_2012->Clone("sig_ghzz_2D_2m2e"); sig_ghzz_2D_2m2e->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  
  TH2F* sig_qhzz_2D = (TH2F*) sig_qhzz_2D_2012->Clone("sig_qhzz_2D"); sig_qhzz_2D->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH2F* sig_qhzz_2D_4m = (TH2F*) sig_qhzz_2D_4m_2012->Clone("sig_qhzz_2D_4m"); sig_qhzz_2D_4m->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH2F* sig_qhzz_2D_4e = (TH2F*) sig_qhzz_2D_4e_2012->Clone("sig_qhzz_2D_4e"); sig_qhzz_2D_4e->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH2F* sig_qhzz_2D_2e2m = (TH2F*) sig_qhzz_2D_2e2m_2012->Clone("sig_qhzz_2D_2e2m"); sig_qhzz_2D_2e2m->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  TH2F* sig_qhzz_2D_2m2e = (TH2F*) sig_qhzz_2D_2m2e_2012->Clone("sig_qhzz_2D_2m2e"); sig_qhzz_2D_2m2e->Scale((lumi_2011 + lumi_2012) / lumi_2012);
  
  // Saving all histograms individually to separate files
  cv = new TCanvas("cv", "Canvas", 800, 600);
  data_2D->SetMarkerSize(1); data_2D->SetMarkerStyle(29); data_2D->Draw("scat"); cv->SaveAs((outDirectoryString + "doHZZHistograms_data_2D.gif").c_str());
  data_2D_4m->SetMarkerSize(1); data_2D_4m->SetMarkerStyle(20); data_2D_4m->Draw("scat"); cv->SaveAs((outDirectoryString + "doHZZHistograms_data_2D_4m.gif").c_str());
  data_2D_4e->SetMarkerSize(1); data_2D_4e->SetMarkerStyle(21); data_2D_4e->Draw("scat"); cv->SaveAs((outDirectoryString + "doHZZHistograms_data_2D_4e.gif").c_str());
  data_2D_2e2m->SetMarkerSize(1); data_2D_2e2m->SetMarkerStyle(22); data_2D_2e2m->Draw("scat"); cv->SaveAs((outDirectoryString + "doHZZHistograms_data_2D_2e2m.gif").c_str());
  data_2D_2m2e->SetMarkerSize(1); data_2D_2m2e->SetMarkerStyle(23); data_2D_2m2e->Draw("scat"); cv->SaveAs((outDirectoryString + "doHZZHistograms_data_2D_2m2e.gif").c_str());

  // Setting different colors too
  data_2D_4m->SetMarkerColor(6);
  data_2D_4e->SetMarkerColor(1);
  data_2D_2e2m->SetMarkerColor(20);
  data_2D_2m2e->SetMarkerColor(28);

  bkg_zxss_2D->SetMarkerStyle(29); bkg_zxss_2D->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_zxss_2D.gif").c_str());
  bkg_zxss_2D_4m->SetMarkerStyle(29); bkg_zxss_2D_4m->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_zxss_2D_4m.gif").c_str());
  bkg_zxss_2D_4e->SetMarkerStyle(29); bkg_zxss_2D_4e->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_zxss_2D_4e.gif").c_str());
  bkg_zxss_2D_2e2m->SetMarkerStyle(29); bkg_zxss_2D_2e2m->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_zxss_2D_2e2m.gif").c_str());
  bkg_zxss_2D_2m2e->SetMarkerStyle(29); bkg_zxss_2D_2m2e->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_zxss_2D_2m2e.gif").c_str());

  bkg_ggzz_2D->SetMarkerStyle(29); bkg_ggzz_2D->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_ggzz_2D.gif").c_str());
  bkg_ggzz_2D_4m->SetMarkerStyle(29); bkg_ggzz_2D_4m->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_ggzz_2D_4m.gif").c_str());
  bkg_ggzz_2D_4e->SetMarkerStyle(29); bkg_ggzz_2D_4e->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_ggzz_2D_4e.gif").c_str());
  bkg_ggzz_2D_2e2m->SetMarkerStyle(29); bkg_ggzz_2D_2e2m->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_ggzz_2D_2e2m.gif").c_str());
  bkg_ggzz_2D_2m2e->SetMarkerStyle(29); bkg_ggzz_2D_2m2e->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_ggzz_2D_2m2e.gif").c_str());

  bkg_qqzz_2D->SetMarkerStyle(29); bkg_qqzz_2D->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_qqzz_2D.gif").c_str());
  bkg_qqzz_2D_4m->SetMarkerStyle(29); bkg_qqzz_2D_4m->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_qqzz_2D_4m.gif").c_str());
  bkg_qqzz_2D_4e->SetMarkerStyle(29); bkg_qqzz_2D_4e->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_qqzz_2D_4e.gif").c_str());
  bkg_qqzz_2D_2e2m->SetMarkerStyle(29); bkg_qqzz_2D_2e2m->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_qqzz_2D_2e2m.gif").c_str());
  bkg_qqzz_2D_2m2e->SetMarkerStyle(29); bkg_qqzz_2D_2m2e->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_bkg_qqzz_2D_2m2e.gif").c_str());

  sig_qhzz_2D->SetMarkerStyle(29); sig_qhzz_2D->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_sig_qhzz_2D.gif").c_str());
  sig_qhzz_2D_4m->SetMarkerStyle(29); sig_qhzz_2D_4m->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_sig_qhzz_2D_4m.gif").c_str());
  sig_qhzz_2D_4e->SetMarkerStyle(29); sig_qhzz_2D_4e->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_sig_qhzz_2D_4e.gif").c_str());
  sig_qhzz_2D_2e2m->SetMarkerStyle(29); sig_qhzz_2D_2e2m->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_sig_qhzz_2D_2e2m.gif").c_str());
  sig_qhzz_2D_2m2e->SetMarkerStyle(29); sig_qhzz_2D_2m2e->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_sig_qhzz_2D_2m2e.gif").c_str());

  sig_ghzz_2D->SetMarkerStyle(29); sig_ghzz_2D->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_sig_ghzz_2D.gif").c_str());
  sig_ghzz_2D_4m->SetMarkerStyle(29); sig_ghzz_2D_4m->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_sig_ghzz_2D_4m.gif").c_str());
  sig_ghzz_2D_4e->SetMarkerStyle(29); sig_ghzz_2D_4e->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_sig_ghzz_2D_4e.gif").c_str());
  sig_ghzz_2D_2e2m->SetMarkerStyle(29); sig_ghzz_2D_2e2m->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_sig_ghzz_2D_2e2m.gif").c_str());
  sig_ghzz_2D_2m2e->SetMarkerStyle(29); sig_ghzz_2D_2m2e->Draw("colz"); cv->SaveAs((outDirectoryString + "doHZZHistograms_sig_ghzz_2D_2m2e.gif").c_str());

  // Adding for total background and signal
  TH2F* bkg_2D = (TH2F*) bkg_zxss_2D->Clone("bkg_2D"); bkg_2D->Add(bkg_ggzz_2D); bkg_2D->Add(bkg_qqzz_2D);
  bkg_2D->SetNameTitle("bkg_2D", "Z1 vs. Z2 masses, all, background");
  TH2F* sig_2D = (TH2F*) sig_qhzz_2D->Clone("sig_2D"); sig_2D->Add(sig_ghzz_2D);
  sig_2D->SetNameTitle("sig_2D", "Z1 vs. Z2 masses, all, signal");
  TH2F* bkg_2D_4m = (TH2F*) bkg_zxss_2D_4m->Clone("bkg_2D_4m"); bkg_2D_4m->Add(bkg_ggzz_2D_4m); bkg_2D_4m->Add(bkg_qqzz_2D_4m);
  bkg_2D_4m->SetNameTitle("bkg_2D_4m", "Z1 vs. Z2 masses, 4m, background");
  TH2F* sig_2D_4m = (TH2F*) sig_qhzz_2D_4m->Clone("sig_2D_4m"); sig_2D_4m->Add(sig_ghzz_2D_4m);
  sig_2D_4m->SetNameTitle("sig_2D_4m", "Z1 vs. Z2 masses, 4m, signal");
  TH2F* bkg_2D_4e = (TH2F*) bkg_zxss_2D_4e->Clone("bkg_2D_4e"); bkg_2D_4e->Add(bkg_ggzz_2D_4e); bkg_2D_4e->Add(bkg_qqzz_2D_4e);
  bkg_2D_4e->SetNameTitle("bkg_2D_4e", "Z1 vs. Z2 masses, 4e, background");
  TH2F* sig_2D_4e = (TH2F*) sig_qhzz_2D_4e->Clone("sig_2D_4e"); sig_2D_4e->Add(sig_ghzz_2D_4e);
  sig_2D_4e->SetNameTitle("sig_2D_4e", "Z1 vs. Z2 masses, 4e, signal");
  TH2F* bkg_2D_2e2m = (TH2F*) bkg_zxss_2D_2e2m->Clone("bkg_2D_2e2m"); bkg_2D_2e2m->Add(bkg_ggzz_2D_2e2m); bkg_2D_2e2m->Add(bkg_qqzz_2D_2e2m);
  bkg_2D_2e2m->SetNameTitle("bkg_2D_2e2m", "Z1 vs. Z2 masses, 2e2m, background");
  TH2F* sig_2D_2e2m = (TH2F*) sig_qhzz_2D_2e2m->Clone("sig_2D_2e2m"); sig_2D_2e2m->Add(sig_ghzz_2D_2e2m);
  sig_2D_2e2m->SetNameTitle("sig_2D_2e2m", "Z1 vs. Z2 masses, 2e2m, signal");
  TH2F* bkg_2D_2m2e = (TH2F*) bkg_zxss_2D_2m2e->Clone("bkg_2D_2m2e"); bkg_2D_2m2e->Add(bkg_ggzz_2D_2m2e); bkg_2D_2m2e->Add(bkg_qqzz_2D_2m2e);
  bkg_2D_2m2e->SetNameTitle("bkg_2D_2m2e", "Z1 vs. Z2 masses, 2m2e, background");
  TH2F* sig_2D_2m2e = (TH2F*) sig_qhzz_2D_2m2e->Clone("sig_2D_2m2e"); sig_2D_2m2e->Add(sig_ghzz_2D_2m2e);
  sig_2D_2m2e->SetNameTitle("sig_2D_2m2e", "Z1 vs. Z2 masses, 2m2e, signal");

  TH2F* both_2D = (TH2F*) bkg_2D->Clone("both_2D"); both_2D->Add(sig_2D);
  TH2F* both_2D_4m = (TH2F*) bkg_2D_4m->Clone("both_2D_4m"); both_2D_4m->Add(sig_2D_4m);
  TH2F* both_2D_4e = (TH2F*) bkg_2D_4e->Clone("both_2D_4e"); both_2D_4e->Add(sig_2D_4e);
  TH2F* both_2D_2e2m = (TH2F*) bkg_2D_2e2m->Clone("both_2D_2e2m"); both_2D_2e2m->Add(sig_2D_2e2m);
  TH2F* both_2D_2m2e = (TH2F*) bkg_2D_2m2e->Clone("both_2D_2m2e"); both_2D_2m2e->Add(sig_2D_2m2e);
  
  
  cv = new TCanvas("cv", "Canvas", 800, 600);
  // Setting log scale to z axis
  //cv->SetLogz();
  bkg_2D->Draw("colz");
  data_2D_4m->Draw("scatsame");
  data_2D_4e->Draw("scatsame");
  data_2D_2e2m->Draw("scatsame");
  data_2D_2m2e->Draw("scatsame");
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectoryString + "doHZZHistograms_all_2D.gif").c_str());
  
  cv = new TCanvas("cv", "Canvas", 800, 600);
  // Setting log scale to z axis
  //cv->SetLogz();
  bkg_2D_4m->Draw("colz");
  data_2D_4m->Draw("scatsame");		// Marker style has already been set
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectoryString + "doHZZHistograms_all_2D_4m.gif").c_str());
  
  cv = new TCanvas("cv", "Canvas", 800, 600);
  // Setting log scale to z axis
  //cv->SetLogz();
  bkg_2D_4e->Draw("colz");
  data_2D_4e->Draw("scatsame");		// Marker style has already been set
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectoryString + "doHZZHistograms_all_2D_4e.gif").c_str());
  
  cv = new TCanvas("cv", "Canvas", 800, 600);
  // Setting log scale to z axis
  //cv->SetLogz();
  bkg_2D_2e2m->Draw("colz");
  data_2D_2e2m->Draw("scatsame");		// Marker style has already been set
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectoryString + "doHZZHistograms_all_2D_2e2m.gif").c_str());
  
  cv = new TCanvas("cv", "Canvas", 800, 600);
  // Setting log scale to z axis
  //cv->SetLogz();
  bkg_2D_2m2e->Draw("colz");
  data_2D_2m2e->Draw("scatsame");		// Marker style has already been set
  legend = cv->BuildLegend(0.6, 0.7, 0.9, 0.9);
  legend->Draw();
  cv->SaveAs((outDirectoryString + "doHZZHistograms_all_2D_2m2e.gif").c_str());
}
