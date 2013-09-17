/** \macro fitJPsi2SMassUnbinned.C
 *
 * $Id: fitZToMuMuGammaMassUnbinnedwtv3data.C,v 1.1 2012/10/24 15:41:26 yangyong Exp $
 *
 *
 * Macro implementing unbinned Maximum-likelihood fit of
 * the Z->mmg lineshape
 *
 * Software developed for the CMS Detector at LHC
 *
 *  \author J. Veverka - Caltech, Pasadena, USA
 *  Inspired by
 *    E. Schneider - Caltech, Pasadena, USA
 *    S. Ganzhur   - CEA/DAPNIA/SPP, Saclay
 *    Y. Yang      - Caltech, Pasadena, USA
 *
 */

#include "TTree.h"
#include "TLatex.h"
#include "TH1.h"

#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "TCanvas.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"



//#include "tdrstyle.C"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "/afs/cern.ch/user/y/yangyong/macros/effSigma.C"

using namespace RooFit;

double NormalizedIntegral(RooAbsPdf & function, RooRealVar & integrationVar, double lowerLimit, double upperLimit){

  integrationVar.setRange("integralRange",lowerLimit,upperLimit);
  RooAbsReal* integral = function.createIntegral(integrationVar,NormSet(integrationVar),Range("integralRange"));


  double normlizedIntegralValue = integral->getVal();

  //  cout<<normlizedIntegralValue<<endl;


  return normlizedIntegralValue;


}

void fitZToMuMuGammaMassUnbinnedv3data(
					float &dm_fit,
					float &dm_fitErr,
					float &a_fit,
				     float &n_fit,
					float &sigcb_fit,
					float &sigcb_fitErr,
					 //const char *filename = "ZMuMuGammaMass_Zmumu_Spring10_EB.txtwt",
					 RooDataSet *dataset,
					 RooRealVar *mass,
					 float xlowFit = 60,
					 float xhighFit = 120,
					 float a_fix = 1.0,
					 float n_fix = 5.0,
					bool usecdf = false,
					bool fix_antoMC = false,
					 const char *gifdir = "",
					 const char *gifName = ""
  )
{
  char* plotOpt = "NEU";

  //const int nbins = 60;
  const int nbins = 80;
  
  
  //gROOT->ProcessLine(".L tdrstyle.C");
  // setTDRStyle();
  gStyle->SetPadRightMargin(0.05);
    
  // Build p.d.f.
  
  ////////////////////////////////////////////////
  //             Parameters                     //
  ////////////////////////////////////////////////
  
  //  Signal p.d.f. parameters
  //  Parameters for a Gaussian and a Crystal Ball Lineshape
  RooRealVar  cbBias ("#Deltam_{CB}", "CB Bias", 0.5, -5, 5,"");
  
  //cbBias.setConstant(kTRUE);
  
  
  RooRealVar  cbSigma("#sigma_{CB}","CB Width", 1.38, 0.01, 10.0,"");
  //RooRealVar  cbCut  ("a_{CB}","CB Cut", 1.5, 0.5, 10);
  RooRealVar  cbCut  ("a_{CB}","CB Cut", 1.5, 1, 10);
  
  //RooRealVar  cbPower("n_{CB}","CB Power", 1.3, 0.5, 10.0);
  //RooRealVar  cbPower("n_{CB}","CB Power", 1.3, 1, 20);
  RooRealVar  cbPower("n_{CB}","CB Power", 1.3, 1, 10);
  
  mass->SetTitle("m_{ee} (GeV/c^{2})");
  
  mass->setRange(xlowFit,xhighFit);
  
  
  if(fix_antoMC){
    cbCut.setVal(a_fix);
    cbCut.setConstant(true);
    cbPower.setVal(n_fix);
    cbPower.setConstant(true);
  }
  

  //   cbSigma.setConstant(kTRUE);
  //   cbCut.setConstant(kTRUE);
  //   cbPower.setConstant(kTRUE);
  
  //  Parameters for Breit-Wigner
  RooRealVar bwMean("m_{Z}","BW Mean", 91.1876, "GeV/c^{2}");
  RooRealVar bwWidth("#Gamma_{Z}", "BW Width", 2.4952, "GeV/c^{2}");
  
  // Keep Breit-Wigner parameters fixed to the PDG values
  //   bwMean.setConstant(kTRUE);
  //   bwWidth.setConstant(kTRUE);
  
  
  //  Background p.d.f. parameters
  // Parameters for exponential
  RooRealVar expRate("#lambda_{exp}", "Exponential Rate", -0.119, -10, 1);
  
  
  // fraction of signal
  //  RooRealVar  frac("frac", "Signal Fraction", 0.1,0.,0.3.);
  /*  RooRealVar  nsig("N_{S}", "#signal events", 9000, 0.,10000.);
      RooRealVar  nbkg("N_{B}", "#background events", 1000,2,10000.);*/
  RooRealVar  nsig("N_{S}", "#signal events", 29300, 0.1, 10000000.);
  RooRealVar  nbkg("N_{B}", "#background events", 0, 0., 100000.);

////////////////////////////////////////////////
//               P.D.F.s                      //
////////////////////////////////////////////////
  
// Di-photon mass signal p.d.f.
  RooBreitWigner bw("bw", "bw", *mass, bwMean, bwWidth);
  //   RooGaussian    signal("signal", "A  Gaussian Lineshape", mass, m0, sigma);
  RooCBShape     cball("cball", "A  Crystal Ball Lineshape", *mass, cbBias, cbSigma, cbCut, cbPower);

  mass->setBins(100000, "fft");

  RooFFTConvPdf BWxCB("BWxCB","bw (X) crystall ball", *mass, bw, cball);
  
  
  // Di-photon mass background  p.d.f.
  //RooExponential bg("bg","bkgd exp", *mass, expRate);
  // Di-photon mass model p.d.f.
  
  ///RooAddPdf      model("model", "signal + background mass model", RooArgList(BWxCB, bg), RooArgList(nsig, nbkg));
  //RooAddPdf  *model; 
  /// p.d.f. provides expected number of events, including extended term in likelihood. the same as by doing RooExtendPdf
  //RooExtendPdf modelext("modelex","",BWxCB,nsig);
  
  //model.fitTo(*dataset,Strategy(0),FitOptions("mh"),SumW2Error(kTRUE),Optimize(0),Timer(1));

  RooExtendPdf *model = new RooExtendPdf("model", "signal + background mass model", BWxCB, nsig);
  //  if(IsData ==1){
  
  //BWxCB.fitTo(*dataset,SumW2Error(kTRUE));
  model->fitTo(*dataset,SumW2Error(kTRUE));
  
  //model->fitTo(*dataset,Strategy(2),FitOptions("mh"),SumW2Error(kTRUE),Optimize(0),Timer(1));
  
  
  TCanvas *c = new TCanvas("c","Unbinned Invariant Mass Fit", 0,0,800,600);
  // Plot the fit results
  RooPlot* plot = mass->frame(Range(xlowFit,xhighFit),Bins(nbins));
  dataset->plotOn(plot,Name("data"));
  model->plotOn(plot,Name("model"));
  
  
  
  plot->Draw();
  TLatex l;
  l.SetTextSize(0.045);
  l.SetTextColor(1);
  l.SetNDC();
  TString result = TString( Form("#sigma_{CB} = %3.2f #pm %3.2f",cbSigma.getVal(),cbSigma.getError())); 
  cout<<"result " << result <<endl; 
  l.DrawLatex(0.2,0.85,result);
  result = TString( Form("#Deltam_{CB} = %3.2f #pm %3.2f",cbBias.getVal(),cbBias.getError())); 
  cout<<"result " << result <<endl; 
  l.DrawLatex(0.2,0.8,result);
  if(!fix_antoMC){
    result = TString( Form("a_{CB} = %2.2f #pm %2.2f",cbCut.getVal(),cbCut.getError())); 
    // l.DrawLatex(0.2,0.75,result);
    result = TString( Form("n_{CB} = %2.2f #pm %2.2f",cbPower.getVal(),cbPower.getError())); 
    // l.DrawLatex(0.2,0.7,result);
  }else{
    result = TString( Form("a_{CB} = %2.2f(MC)",cbCut.getVal()));
    //l.DrawLatex(0.2,0.75,result);
    result = TString( Form("n_{CB} = %2.2f(MC)",cbPower.getVal()));
    //l.DrawLatex(0.2,0.7,result);
  }
  

  a_fit = cbCut.getVal();
  n_fit = cbPower.getVal();
  dm_fit = cbBias.getVal();
  dm_fitErr = cbBias.getError();
  sigcb_fit = cbSigma.getVal();
  sigcb_fitErr = cbSigma.getError();
  
  int npar = 6; 
  if(fix_antoMC){
    npar =-2;
  }
  
  double chi2 = plot->chiSquare("model","data",npar);
  double ndoff = plot->GetNbinsX() - npar;
  TString histname = TString("th1f_") + dataset->GetName(); 
  TH1F *hhtemp = convertRooDataSetToTH1F(dataset,mass,histname,xlowFit,xhighFit,1);
  float eff_sig = effSigma(hhtemp);


  //calculate effective sigma and corresponding interval
  if( usecdf){
   RooAbsReal *cdf = model->createCdf(*mass);
  
  float testmass = 91.1876; 
  
  float center = testmass-10.0;
  float minwidth = 999.0;
  float mlmin = 0.0;
  float mhmin = 0.0;
  
  float step=0.01;
  int Nstep = int(20/step+0.1);
  
  cout<<"cdf created " << Nstep <<endl; 

  int kkk = 0; 
  for (int i=0; i<Nstep; ++i) {

    if(i % 100 ==0) cout<<"i " << i <<endl; 
    
    float mlow = center+i*step;
      mass->setVal(mlow);
      float cdflo = cdf->getVal();
      for (int j=i+1; j<Nstep; ++j) {
        float mhigh = center+j*step;
        mass->setVal(mhigh);
        float cdfhi = cdf->getVal();
        if ( (cdfhi-cdflo)>0.683 ) {
          if ( (mhigh-mlow)<minwidth) {
            minwidth = mhigh-mlow;
            mlmin = mlow;
            mhmin = mhigh;
          }
          break;

        }
      }
  }
  
  
  float sigmaeff = minwidth/2.0;
  
  cout<<"cdf created sigmaeff done " <<endl; 

  //calculate FWHM and corresponding interval
  //fullpdf->getDependents()->Print("V");
  //return;
  ///RooRealVar *clonedmass = (RooRealVar*)model->getDependents(&RooArgSet(*mass))->find("mass");
  double mmax = 0.0;
  double pmax = 0.0;
    double mleft = 0.0;
    double mright = 0.0;
    double dpminleft = 999.;
    double dpminright = 999.;
    for (int i=0; i<Nstep; ++i) {

      double mlow = center+i*step;
      mass->setVal(mlow);
      double pval = model->getVal(&RooArgSet(*mass));
      if (pval>pmax) {
        pmax = pval;
        mmax = mlow;
      }
    }
    
    double phalf = pmax/2.0;
    for (int i=0; i<Nstep; ++i) {
      double mlow = center+i*step;
      //clonedmass->setVal(mlow);
      mass->setVal(mlow);
      double pval = model->getVal(&RooArgSet(*mass));
      double dp = fabs(pval-phalf);
      
      if (dp<dpminleft && mlow<mmax) {
        dpminleft = dp;
        mleft = mlow;
      }
      
      if (dp<dpminright && mlow>mmax) {
        dpminright = dp;
        mright = mlow;
      }      
      
    }    
    double fwhm = mright-mleft;
    
    
    cout<<"cdf created fwhm done" <<endl; 
 

 
  result = TString( Form("#sigma_{eff}= %3.2f",sigmaeff));
  // l.DrawLatex(0.2,0.65,result);
  
  result = TString( Form("FWHM = %3.2f",fwhm));
  // l.DrawLatex(0.2,0.6,result);
  }else{
    result = TString( Form("#sigma_{eff}= %3.2f",eff_sig));
    // l.DrawLatex(0.2,0.65,result);
  }
  
  
  
  result = TString( Form("#chi^{2}/d.o.f. = %3.3f",chi2/ndoff));
  //l.DrawLatex(0.2,0.5,result);
  hhtemp->Delete();
    

  TString sgifdir = TString(gifdir);
  if( sgifdir != ""){
    TString gfname = sgifdir + TString("/") + TString(gifName) + TString(".gif");
    //c->Print(gfname);
    gfname = sgifdir + TString("/") + TString(gifName) + TString(".pdf");
    c->Print(gfname);
    gfname = sgifdir + TString("/") + TString(gifName) + TString(".C");
    //c->Print(gfname);
  }
  

    
  //return BWxCB;
  
}

