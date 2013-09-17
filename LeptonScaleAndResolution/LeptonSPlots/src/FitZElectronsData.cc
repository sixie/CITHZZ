struct fitbinning {
  enum ptBins {kPtLow = 0,
               kPtHigh = 1}; 
  
  enum etaBins {kCentralEB = 0,
                kOuterEB,
                kEE};
};

//===================================
// make the RooDataSet from ROOT tree
//=================================== 
RooDataSet *getDataSet(const char *rootfile) {
  cout << "roofitting file " << rootfile << endl;

  // fit variables
  RooRealVar* mass = new RooRealVar("mass",  "Mass [GeV]", 60., 120.);

  // kinematic cut variables
  RooRealVar* pt = new RooRealVar("pt", "p_{T} [GeV]",5,1000);
  RooRealVar* eta = new RooRealVar("eta","#eta",-2.5,2.5);

  RooArgSet setkine(*mass,*pt,*eta);
  TFile *file = TFile::Open(rootfile);
  TTree *tree = (TTree*)file->Get("eleIDdir/T1");

  RooDataSet *data = new RooDataSet("T1","dataset",tree,setkine);

  return data;
}

// Set Fit Options
MLOptions GetDefaultOptions() {
  MLOptions opts;
  // Fit configuration
  opts.addBoolOption("useinvMass", "Use invMass", kTRUE);        
  return opts;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

MLFit theFit;

void myFit() {

  MLFit theFit;

  // Various fit options...
  MLOptions opts = GetDefaultOptions();
  
  // define the structure of the dataset
  RooRealVar* mass = new RooRealVar("mass",  "Mass [GeV/c^{2}]" , 60., 120.);
  
  theFit.AddFlatFileColumn(mass);
  
  // define a fit model
  theFit.addModel("myFit", "ZtoEE");
  
  // define species
  theFit.addSpecies("myFit", "sig", "Signal Component");
  theFit.addSpecies("myFit", "bkg", "Bkg   Component");
  
  // mLL PDF
  if(opts.getBoolVal("useinvMass")) {
    theFit.addPdfWName("myFit", "sig" , "mass",  "Cruijff", "sig_Mass");
    theFit.addPdfWName("myFit", "bkg" , "mass",  "Expo",   "bkg_Mass");
  }
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Fit a sample of Z events
void FitZElectrons(int ptbin=-1, int etabin=-1) {
  
  myFit();

  // Load the data
  char datasetname[200];
  sprintf(datasetname,"root://eoscms//eos/cms/store/group/phys_egamma/emanuele/elereg/eletrees/V1/DoubleElectron_Run2012AB_ZTaP.root");
  RooDataSet *data = getDataSet(datasetname);
  if(ptbin==fitbinning::kPtLow) data = (RooDataSet*)data->reduce("pt>=7 && pt<15");
  if(ptbin==fitbinning::kPtHigh) data = (RooDataSet*)data->reduce("pt>=15");
  if(etabin==fitbinning::kCentralEB) data = (RooDataSet*)data->reduce("abs(eta)<0.8");
  if(etabin==fitbinning::kOuterEB) data = (RooDataSet*)data->reduce("abs(eta)>=0.8 && abs(eta)<1.479");
  if(etabin==fitbinning::kEE) data = (RooDataSet*)data->reduce("abs(eta)>=1.479");

  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");
  
  // Initialize the fit...
  char fitconfig[200];
  //  sprintf(fitconfig,"fitconfig/fitZee-ptbin%d-etabin%d.config",ptbin,etabin);
  sprintf(fitconfig,"fitconfig/fitZee-ptbin0-etabin0.config"); 
  theFit.initialize(fitconfig);
  
  // Print Fit configuration 
  myPdf->getParameters(data)->selectByAttrib("Constant",kTRUE)->Print("V");  
  myPdf->getParameters(data)->selectByAttrib("Constant",kFALSE)->Print("V");
  
  RooFitResult *fitres =  myPdf->fitTo(*data,
                                       RooFit::ConditionalObservables(theFit.getNoNormVars("myFit")),
                                       RooFit::FitOptions("MHTER"),
                                       RooFit::NumCPU(4));
  fitres->Print("V");
  
  // write the config file corresponding to the fit minimum
  char configfilename[200];
  sprintf(configfilename, "fitres/fitResults-ptbin%d-etabin%d.config",ptbin,etabin);
  theFit.writeConfigFile(configfilename);  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PlotZElectrons(int ptbin=-1, int etabin=-1, int nbins=19) {

  myFit();

  // Load the data
  char datasetname[200];
  sprintf(datasetname,"root://eoscms//eos/cms/store/group/phys_egamma/emanuele/elereg/eletrees/V1/DoubleElectron_Run2012AB_ZTaP.root");
  RooDataSet *data = getDataSet(datasetname);
  if(ptbin==fitbinning::kPtLow) data = (RooDataSet*)data->reduce("pt>=7 && pt<15");
  if(ptbin==fitbinning::kPtHigh) data = (RooDataSet*)data->reduce("pt>=15");
  if(etabin==fitbinning::kCentralEB) data = (RooDataSet*)data->reduce("abs(eta)<0.8");
  if(etabin==fitbinning::kOuterEB) data = (RooDataSet*)data->reduce("abs(eta)>=0.8 && abs(eta)<1.479");
  if(etabin==fitbinning::kEE) data = (RooDataSet*)data->reduce("abs(eta)>=1.479");

  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");

  // Initialize the fit...
  char configfilename[200];
  sprintf(configfilename, "fitres/fitResults-ptbin%d-etabin%d.config",ptbin,etabin);
  theFit.initialize(configfilename);

  TCanvas *c = new TCanvas("c","fitResult");
  char rootfilename[200];
  sprintf(rootfilename,"fitres/mll-ptbin%d-etabin%d.root",ptbin,etabin);
  TFile *output = new TFile(rootfilename,"RECREATE");

  RooPlot* MassPlot = MakePlot("mass", &theFit, data, nbins, false);    

  MassPlot->SetYTitle("Events");
  MassPlot->Draw();
  char epsfilename[200];
  sprintf(epsfilename,"fitres/mll-ptbin%d-etabin%d.eps",ptbin,etabin);
  c->SaveAs(epsfilename);
  MassPlot->Write();
  //  output->Close();

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Make the plot for a given variable
RooPlot *MakePlot(TString VarName, MLFit* theFit, RooDataSet* theData, int nbins, bool poissonError=true)
{
  RooRealVar* Var = theFit->RealObservable(VarName);
  double min=Var->getMin();
  double max=Var->getMax();
  RooPlot *plot = Var->frame(min,max,nbins);
  
  // plot the data
  if(poissonError)
    theData->plotOn(plot);
  else 
    theData->plotOn(plot, RooFit::DataError(RooAbsData::SumW2) );

  // plot the total likelihood
  RooAbsPdf *thePdf = theFit->getPdf("myFit");
  thePdf->plotOn(plot, RooFit::LineColor(kBlack));

  double Ns = theFit->getRealPar("N_sig")->getVal();
  double Nb = theFit->getRealPar("N_bkg")->getVal();

  // plot (dashed) the bkg component
  theFit->getRealPar("N_sig")->setVal(0.);
  thePdf->plotOn(plot, RooFit::Normalization(Nb/(Ns+Nb)),RooFit::LineColor(kBlack),RooFit::LineStyle(kDashed));

  
  return plot;
}

// =========================================================
// *** Calculate the sWeights from the minimum of the LH
// =========================================================

void calcSWeight(int ptbin=-1, int etabin=-1) {

  myFit();

  // Load the data
  // reduction has to be done hard on the tree because reduce does not propagate correctly the tree
  char datasetname[1000];
  sprintf(datasetname,"root://eoscms//eos/cms/store/group/phys_egamma/emanuele/elereg/eletrees/V1/DoubleElectron_Run2012AB_ZTaP.root");
  char cut[200];
  if(ptbin==fitbinning::kPtLow) sprintf(cut,"pt>=7 && pt<15");
  if(ptbin==fitbinning::kPtHigh) sprintf(cut,"pt>=15");
  if(etabin==fitbinning::kCentralEB) sprintf(cut,"%s && abs(eta)<0.8",cut);
  if(etabin==fitbinning::kOuterEB) sprintf(cut,"%s && abs(eta)>=0.8 && abs(eta)<1.479",cut);
  if(etabin==fitbinning::kEE) sprintf(cut,"%s && abs(eta)>=1.479",cut);

  gROOT->LoadMacro("src/reduceTree.cc");
  char filesmall[200];
  sprintf(filesmall,"sPlots/electrons_ptbin%d_etabin%d.root",ptbin,etabin);
  makeSmall(datasetname,"eleIDdir/T1",cut,filesmall);
  RooDataSet *data = getDataSet(filesmall);

  // build the fit likelihood
  RooAbsPdf *myPdf = theFit.buildModel("myFit");
  
  // Initialize the fit: the first fit should be already done and the minimum config ready...
  char fitconfig[200];
  sprintf(fitconfig, "fitres/fitResults-ptbin%d-etabin%d.config",ptbin,etabin);
  theFit.initialize(fitconfig);

  // fix all parameters, float the yields and fit
  theFit._parameterSet.selectByName("*")->setAttribAll("Constant",kTRUE);
  (static_cast<RooRealVar*>(theFit._parameterSet.find("N_sig")))->setConstant(kFALSE) ;
  (static_cast<RooRealVar*>(theFit._parameterSet.find("N_bkg")))->setConstant(kFALSE) ;
  
  cout << "printing data: " << endl;
  cout << data->tree() << endl;
  data->Print("V");

  RooFitResult *fitres =  myPdf->fitTo(*data,
                                       RooFit::ConditionalObservables(theFit.getNoNormVars("myFit")),
                                       RooFit::FitOptions("MHTER"));
  fitres->Print("V");

  // add appropriate column to dataset
  RooArgList yieldsList;
  yieldsList.add(*theFit._fracList.find("N_sig"));
  yieldsList.add(*theFit._fracList.find("N_bkg"));
  cout << "number of entries in set to write: " << data->numEntries() << endl ;
  RooArgSet nonormvars;
  RooDataSet* dsnew = MLSPlot::addSWeightToData((RooSimultaneous*)(myPdf), yieldsList, *data, nonormvars) ;

  char filesweight[200];
  sprintf(filesweight,"sPlots/sweights_ptbin%d_etabin%d.root",ptbin,etabin);
  TFile sPlotsTree(filesweight,"recreate");

  TreeFillerFromRooDataSet sweights(dsnew);
  sweights.addVar("N_sig_sw");
  sweights.addVar("N_bkg_sw");
  TTree *treewithweights = sweights.getTree();
  treewithweights->Write();
  sPlotsTree.Close();
}
