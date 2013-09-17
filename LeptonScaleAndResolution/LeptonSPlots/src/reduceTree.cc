#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TEventList.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>

using namespace std;

void makeList(const char *file, const char* tree, const char* cut) {
  cout << "Reducing file " << file << endl;
  TFile *f1 = TFile::Open(file);
  TTree *ntuple = (TTree*) f1->Get(tree);
  ntuple->Draw(">>elist",cut);
  TEventList *elist = (TEventList*)gDirectory->Get("elist");
  TFile ef("elist.root","recreate");
  elist->Write();
}

void makeSmall(const char *file, const char* tree, const char* cut, const char* filesmall) {
  makeList(file,tree,cut);
  TFile *f = TFile::Open("elist.root");
  TEventList *elist = (TEventList*)f->Get("elist");
  TFile *f1 = TFile::Open(file);
  TTree *ntuple = (TTree*) f1->Get(tree);
  ntuple->SetEventList(elist);
  TFile *f2 = TFile::Open(filesmall,"recreate");
  f2->mkdir("eleIDdir");
  TTree *small = ntuple->CopyTree("");
  f2->cd("eleIDdir");
  small->Write();
  f2->Close();
}


