//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 18 16:22:12 2012 by ROOT version 5.32/00
// from TTree probe_tree/probe_tree
// found on file: hcp8.root
//////////////////////////////////////////////////////////

#ifndef X_h
#define X_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class TemplateClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         bdtScalarVsBkg_125;
   Float_t         channel;
   Float_t         elPtMin;
   Float_t         genhiggsmass;
   Float_t         genhiggsmassweight;
   Float_t         hypo;
   Float_t         intimeSimVertices;
   Float_t         jet1eta;
   Float_t         jet1phi;
   Float_t         jet1pt;
   Float_t         jet2eta;
   Float_t         jet2phi;
   Float_t         jet2pt;
   Float_t         l1bdtIso;
   Float_t         l1dz;
   Float_t         l1eta;
   Float_t         l1idNew;
   Float_t         l1idPRL;
   Float_t         l1ip2d;
   Float_t         l1pdgId;
   Float_t         l1pfIsoChHad04;
   Float_t         l1pfIsoComb04EACorr;
   Float_t         l1phi;
   Float_t         l1pt;
   Float_t         l1q;
   Float_t         l1sip3d;
   Float_t         l1trig;
   Float_t         l2bdtIso;
   Float_t         l2dz;
   Float_t         l2eta;
   Float_t         l2idNew;
   Float_t         l2idPRL;
   Float_t         l2ip2d;
   Float_t         l2pdgId;
   Float_t         l2pfIsoChHad04;
   Float_t         l2pfIsoComb04EACorr;
   Float_t         l2phi;
   Float_t         l2pt;
   Float_t         l2q;
   Float_t         l2sip3d;
   Float_t         l2trig;
   Float_t         l3bdtIso;
   Float_t         l3dz;
   Float_t         l3eta;
   Float_t         l3idNew;
   Float_t         l3idPRL;
   Float_t         l3ip2d;
   Float_t         l3pdgId;
   Float_t         l3pfIsoChHad04;
   Float_t         l3pfIsoComb04EACorr;
   Float_t         l3phi;
   Float_t         l3pt;
   Float_t         l3q;
   Float_t         l3sip3d;
   Float_t         l3trig;
   Float_t         l4bdtIso;
   Float_t         l4dz;
   Float_t         l4eta;
   Float_t         l4idNew;
   Float_t         l4idPRL;
   Float_t         l4ip2d;
   Float_t         l4pdgId;
   Float_t         l4pfIsoChHad04;
   Float_t         l4pfIsoComb04EACorr;
   Float_t         l4phi;
   Float_t         l4pt;
   Float_t         l4q;
   Float_t         l4sip3d;
   Float_t         l4trig;
   Float_t         m4l;
   Float_t         mass;
   Float_t         massErr;
   Float_t         melaCosTheta1;
   Float_t         melaCosTheta2;
   Float_t         melaCosThetaStar;
   Float_t         melaLD;
   Float_t         melaPSLD;
   Float_t         melaPhi;
   Float_t         melaPhiStar1;
   Float_t         melaPt;
   Float_t         melaPtY;
   Float_t         melaSpinOneEven;
   Float_t         melaSpinOneOdd;
   Float_t         melaSpinTwoMinimal;
   Float_t         melaY;
   Float_t         mjj;
   Float_t         mll13;
   Float_t         mll14;
   Float_t         mll23;
   Float_t         mll24;
   Float_t         muPtMin;
   Float_t         njets30;
   Float_t         numTrueInteractions;
   Float_t         pfmet;
   Float_t         pho1dr;
   Float_t         pho1eta;
   Float_t         pho1iso;
   Float_t         pho1phi;
   Float_t         pho1pt;
   Float_t         pho2dr;
   Float_t         pho2eta;
   Float_t         pho2iso;
   Float_t         pho2phi;
   Float_t         pho2pt;
   Float_t         pt;
   Float_t         qll13;
   Float_t         qll14;
   Float_t         qll23;
   Float_t         qll24;
   Float_t         rap;
   Float_t         recoVertices;
   Float_t         rho;
   Float_t         worsePairCombRelIsoBaseline;
   Float_t         z1eta;
   Float_t         z1mass;
   Float_t         z1mll;
   Float_t         z1pt;
   Float_t         z1rap;
   Float_t         z2eta;
   Float_t         z2mass;
   Float_t         z2mll;
   Float_t         z2pt;
   Float_t         z2rap;
   Int_t           fourMassCut4Any;
   Int_t           fourMassCut4AnyS;
   Int_t           fourMassCut4SF;
   Int_t           mc_4l;
   Int_t           mc_4lany;
   Int_t           mc_l1;
   Int_t           mc_l2;
   Int_t           mc_l3;
   Int_t           mc_l4;
   Int_t           mc_z1;
   Int_t           mc_z2;
   Int_t           threeMassCut12SF;
   UInt_t          run;
   UInt_t          lumi;
   UInt_t          event;
   Double_t        SuperMELA;

   // List of branches
   TBranch        *b_bdtScalarVsBkg_125;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_elPtMin;   //!
   TBranch        *b_genhiggsmass;   //!
   TBranch        *b_genhiggsmassweight;   //!
   TBranch        *b_hypo;   //!
   TBranch        *b_intimeSimVertices;   //!
   TBranch        *b_jet1eta;   //!
   TBranch        *b_jet1phi;   //!
   TBranch        *b_jet1pt;   //!
   TBranch        *b_jet2eta;   //!
   TBranch        *b_jet2phi;   //!
   TBranch        *b_jet2pt;   //!
   TBranch        *b_l1bdtIso;   //!
   TBranch        *b_l1dz;   //!
   TBranch        *b_l1eta;   //!
   TBranch        *b_l1idNew;   //!
   TBranch        *b_l1idPRL;   //!
   TBranch        *b_l1ip2d;   //!
   TBranch        *b_l1pdgId;   //!
   TBranch        *b_l1pfIsoChHad04;   //!
   TBranch        *b_l1pfIsoComb04EACorr;   //!
   TBranch        *b_l1phi;   //!
   TBranch        *b_l1pt;   //!
   TBranch        *b_l1q;   //!
   TBranch        *b_l1sip3d;   //!
   TBranch        *b_l1trig;   //!
   TBranch        *b_l2bdtIso;   //!
   TBranch        *b_l2dz;   //!
   TBranch        *b_l2eta;   //!
   TBranch        *b_l2idNew;   //!
   TBranch        *b_l2idPRL;   //!
   TBranch        *b_l2ip2d;   //!
   TBranch        *b_l2pdgId;   //!
   TBranch        *b_l2pfIsoChHad04;   //!
   TBranch        *b_l2pfIsoComb04EACorr;   //!
   TBranch        *b_l2phi;   //!
   TBranch        *b_l2pt;   //!
   TBranch        *b_l2q;   //!
   TBranch        *b_l2sip3d;   //!
   TBranch        *b_l2trig;   //!
   TBranch        *b_l3bdtIso;   //!
   TBranch        *b_l3dz;   //!
   TBranch        *b_l3eta;   //!
   TBranch        *b_l3idNew;   //!
   TBranch        *b_l3idPRL;   //!
   TBranch        *b_l3ip2d;   //!
   TBranch        *b_l3pdgId;   //!
   TBranch        *b_l3pfIsoChHad04;   //!
   TBranch        *b_l3pfIsoComb04EACorr;   //!
   TBranch        *b_l3phi;   //!
   TBranch        *b_l3pt;   //!
   TBranch        *b_l3q;   //!
   TBranch        *b_l3sip3d;   //!
   TBranch        *b_l3trig;   //!
   TBranch        *b_l4bdtIso;   //!
   TBranch        *b_l4dz;   //!
   TBranch        *b_l4eta;   //!
   TBranch        *b_l4idNew;   //!
   TBranch        *b_l4idPRL;   //!
   TBranch        *b_l4ip2d;   //!
   TBranch        *b_l4pdgId;   //!
   TBranch        *b_l4pfIsoChHad04;   //!
   TBranch        *b_l4pfIsoComb04EACorr;   //!
   TBranch        *b_l4phi;   //!
   TBranch        *b_l4pt;   //!
   TBranch        *b_l4q;   //!
   TBranch        *b_l4sip3d;   //!
   TBranch        *b_l4trig;   //!
   TBranch        *b_m4l;   //!
   TBranch        *b_mass;   //!
   TBranch        *b_massErr;   //!
   TBranch        *b_melaCosTheta1;   //!
   TBranch        *b_melaCosTheta2;   //!
   TBranch        *b_melaCosThetaStar;   //!
   TBranch        *b_melaLD;   //!
   TBranch        *b_melaPSLD;   //!
   TBranch        *b_melaPhi;   //!
   TBranch        *b_melaPhiStar1;   //!
   TBranch        *b_melaPt;   //!
   TBranch        *b_melaPtY;   //!
   TBranch        *b_melaSpinOneEven;   //!
   TBranch        *b_melaSpinOneOdd;   //!
   TBranch        *b_melaSpinTwoMinimal;   //!
   TBranch        *b_melaY;   //!
   TBranch        *b_mjj;   //!
   TBranch        *b_mll13;   //!
   TBranch        *b_mll14;   //!
   TBranch        *b_mll23;   //!
   TBranch        *b_mll24;   //!
   TBranch        *b_muPtMin;   //!
   TBranch        *b_njets30;   //!
   TBranch        *b_numTrueInteractions;   //!
   TBranch        *b_pfmet;   //!
   TBranch        *b_pho1dr;   //!
   TBranch        *b_pho1eta;   //!
   TBranch        *b_pho1iso;   //!
   TBranch        *b_pho1phi;   //!
   TBranch        *b_pho1pt;   //!
   TBranch        *b_pho2dr;   //!
   TBranch        *b_pho2eta;   //!
   TBranch        *b_pho2iso;   //!
   TBranch        *b_pho2phi;   //!
   TBranch        *b_pho2pt;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_qll13;   //!
   TBranch        *b_qll14;   //!
   TBranch        *b_qll23;   //!
   TBranch        *b_qll24;   //!
   TBranch        *b_rap;   //!
   TBranch        *b_recoVertices;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_worsePairCombRelIsoBaseline;   //!
   TBranch        *b_z1eta;   //!
   TBranch        *b_z1mass;   //!
   TBranch        *b_z1mll;   //!
   TBranch        *b_z1pt;   //!
   TBranch        *b_z1rap;   //!
   TBranch        *b_z2eta;   //!
   TBranch        *b_z2mass;   //!
   TBranch        *b_z2mll;   //!
   TBranch        *b_z2pt;   //!
   TBranch        *b_z2rap;   //!
   TBranch        *b_fourMassCut4Any;   //!
   TBranch        *b_fourMassCut4AnyS;   //!
   TBranch        *b_fourMassCut4SF;   //!
   TBranch        *b_mc_4l;   //!
   TBranch        *b_mc_4lany;   //!
   TBranch        *b_mc_l1;   //!
   TBranch        *b_mc_l2;   //!
   TBranch        *b_mc_l3;   //!
   TBranch        *b_mc_l4;   //!
   TBranch        *b_mc_z1;   //!
   TBranch        *b_mc_z2;   //!
   TBranch        *b_threeMassCut12SF;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_SuperMELA;   //!

   TemplateClass(TTree *tree=0);
   virtual ~TemplateClass();
   virtual void     Init(TTree *tree);
   virtual void     Branch(TTree *tree);
};

TemplateClass::TemplateClass(TTree *tree) : fChain(0) 
{
}

TemplateClass::~TemplateClass()
{
}

void TemplateClass::Init(TTree *tree)
{
   TTree *fChain = tree;

   fChain->SetBranchAddress("bdtScalarVsBkg_125", &bdtScalarVsBkg_125, &b_bdtScalarVsBkg_125);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("elPtMin", &elPtMin, &b_elPtMin);
   fChain->SetBranchAddress("genhiggsmass", &genhiggsmass, &b_genhiggsmass);
   fChain->SetBranchAddress("genhiggsmassweight", &genhiggsmassweight, &b_genhiggsmassweight);
   fChain->SetBranchAddress("hypo", &hypo, &b_hypo);
   fChain->SetBranchAddress("intimeSimVertices", &intimeSimVertices, &b_intimeSimVertices);
   fChain->SetBranchAddress("jet1eta", &jet1eta, &b_jet1eta);
   fChain->SetBranchAddress("jet1phi", &jet1phi, &b_jet1phi);
   fChain->SetBranchAddress("jet1pt", &jet1pt, &b_jet1pt);
   fChain->SetBranchAddress("jet2eta", &jet2eta, &b_jet2eta);
   fChain->SetBranchAddress("jet2phi", &jet2phi, &b_jet2phi);
   fChain->SetBranchAddress("jet2pt", &jet2pt, &b_jet2pt);
   fChain->SetBranchAddress("l1bdtIso", &l1bdtIso, &b_l1bdtIso);
   fChain->SetBranchAddress("l1dz", &l1dz, &b_l1dz);
   fChain->SetBranchAddress("l1eta", &l1eta, &b_l1eta);
   fChain->SetBranchAddress("l1idNew", &l1idNew, &b_l1idNew);
   fChain->SetBranchAddress("l1idPRL", &l1idPRL, &b_l1idPRL);
   fChain->SetBranchAddress("l1ip2d", &l1ip2d, &b_l1ip2d);
   fChain->SetBranchAddress("l1pdgId", &l1pdgId, &b_l1pdgId);
   fChain->SetBranchAddress("l1pfIsoChHad04", &l1pfIsoChHad04, &b_l1pfIsoChHad04);
   fChain->SetBranchAddress("l1pfIsoComb04EACorr", &l1pfIsoComb04EACorr, &b_l1pfIsoComb04EACorr);
   fChain->SetBranchAddress("l1phi", &l1phi, &b_l1phi);
   fChain->SetBranchAddress("l1pt", &l1pt, &b_l1pt);
   fChain->SetBranchAddress("l1q", &l1q, &b_l1q);
   fChain->SetBranchAddress("l1sip3d", &l1sip3d, &b_l1sip3d);
   fChain->SetBranchAddress("l1trig", &l1trig, &b_l1trig);
   fChain->SetBranchAddress("l2bdtIso", &l2bdtIso, &b_l2bdtIso);
   fChain->SetBranchAddress("l2dz", &l2dz, &b_l2dz);
   fChain->SetBranchAddress("l2eta", &l2eta, &b_l2eta);
   fChain->SetBranchAddress("l2idNew", &l2idNew, &b_l2idNew);
   fChain->SetBranchAddress("l2idPRL", &l2idPRL, &b_l2idPRL);
   fChain->SetBranchAddress("l2ip2d", &l2ip2d, &b_l2ip2d);
   fChain->SetBranchAddress("l2pdgId", &l2pdgId, &b_l2pdgId);
   fChain->SetBranchAddress("l2pfIsoChHad04", &l2pfIsoChHad04, &b_l2pfIsoChHad04);
   fChain->SetBranchAddress("l2pfIsoComb04EACorr", &l2pfIsoComb04EACorr, &b_l2pfIsoComb04EACorr);
   fChain->SetBranchAddress("l2phi", &l2phi, &b_l2phi);
   fChain->SetBranchAddress("l2pt", &l2pt, &b_l2pt);
   fChain->SetBranchAddress("l2q", &l2q, &b_l2q);
   fChain->SetBranchAddress("l2sip3d", &l2sip3d, &b_l2sip3d);
   fChain->SetBranchAddress("l2trig", &l2trig, &b_l2trig);
   fChain->SetBranchAddress("l3bdtIso", &l3bdtIso, &b_l3bdtIso);
   fChain->SetBranchAddress("l3dz", &l3dz, &b_l3dz);
   fChain->SetBranchAddress("l3eta", &l3eta, &b_l3eta);
   fChain->SetBranchAddress("l3idNew", &l3idNew, &b_l3idNew);
   fChain->SetBranchAddress("l3idPRL", &l3idPRL, &b_l3idPRL);
   fChain->SetBranchAddress("l3ip2d", &l3ip2d, &b_l3ip2d);
   fChain->SetBranchAddress("l3pdgId", &l3pdgId, &b_l3pdgId);
   fChain->SetBranchAddress("l3pfIsoChHad04", &l3pfIsoChHad04, &b_l3pfIsoChHad04);
   fChain->SetBranchAddress("l3pfIsoComb04EACorr", &l3pfIsoComb04EACorr, &b_l3pfIsoComb04EACorr);
   fChain->SetBranchAddress("l3phi", &l3phi, &b_l3phi);
   fChain->SetBranchAddress("l3pt", &l3pt, &b_l3pt);
   fChain->SetBranchAddress("l3q", &l3q, &b_l3q);
   fChain->SetBranchAddress("l3sip3d", &l3sip3d, &b_l3sip3d);
   fChain->SetBranchAddress("l3trig", &l3trig, &b_l3trig);
   fChain->SetBranchAddress("l4bdtIso", &l4bdtIso, &b_l4bdtIso);
   fChain->SetBranchAddress("l4dz", &l4dz, &b_l4dz);
   fChain->SetBranchAddress("l4eta", &l4eta, &b_l4eta);
   fChain->SetBranchAddress("l4idNew", &l4idNew, &b_l4idNew);
   fChain->SetBranchAddress("l4idPRL", &l4idPRL, &b_l4idPRL);
   fChain->SetBranchAddress("l4ip2d", &l4ip2d, &b_l4ip2d);
   fChain->SetBranchAddress("l4pdgId", &l4pdgId, &b_l4pdgId);
   fChain->SetBranchAddress("l4pfIsoChHad04", &l4pfIsoChHad04, &b_l4pfIsoChHad04);
   fChain->SetBranchAddress("l4pfIsoComb04EACorr", &l4pfIsoComb04EACorr, &b_l4pfIsoComb04EACorr);
   fChain->SetBranchAddress("l4phi", &l4phi, &b_l4phi);
   fChain->SetBranchAddress("l4pt", &l4pt, &b_l4pt);
   fChain->SetBranchAddress("l4q", &l4q, &b_l4q);
   fChain->SetBranchAddress("l4sip3d", &l4sip3d, &b_l4sip3d);
   fChain->SetBranchAddress("l4trig", &l4trig, &b_l4trig);
   fChain->SetBranchAddress("m4l", &m4l, &b_m4l);
   fChain->SetBranchAddress("mass", &mass, &b_mass);
   fChain->SetBranchAddress("massErr", &massErr, &b_massErr);
   fChain->SetBranchAddress("melaCosTheta1", &melaCosTheta1, &b_melaCosTheta1);
   fChain->SetBranchAddress("melaCosTheta2", &melaCosTheta2, &b_melaCosTheta2);
   fChain->SetBranchAddress("melaCosThetaStar", &melaCosThetaStar, &b_melaCosThetaStar);
   fChain->SetBranchAddress("melaLD", &melaLD, &b_melaLD);
   fChain->SetBranchAddress("melaPSLD", &melaPSLD, &b_melaPSLD);
   fChain->SetBranchAddress("melaPhi", &melaPhi, &b_melaPhi);
   fChain->SetBranchAddress("melaPhiStar1", &melaPhiStar1, &b_melaPhiStar1);
   fChain->SetBranchAddress("melaPt", &melaPt, &b_melaPt);
   fChain->SetBranchAddress("melaPtY", &melaPtY, &b_melaPtY);
   fChain->SetBranchAddress("melaSpinOneEven", &melaSpinOneEven, &b_melaSpinOneEven);
   fChain->SetBranchAddress("melaSpinOneOdd", &melaSpinOneOdd, &b_melaSpinOneOdd);
   fChain->SetBranchAddress("melaSpinTwoMinimal", &melaSpinTwoMinimal, &b_melaSpinTwoMinimal);
   fChain->SetBranchAddress("melaY", &melaY, &b_melaY);
   fChain->SetBranchAddress("mjj", &mjj, &b_mjj);
   fChain->SetBranchAddress("mll13", &mll13, &b_mll13);
   fChain->SetBranchAddress("mll14", &mll14, &b_mll14);
   fChain->SetBranchAddress("mll23", &mll23, &b_mll23);
   fChain->SetBranchAddress("mll24", &mll24, &b_mll24);
   fChain->SetBranchAddress("muPtMin", &muPtMin, &b_muPtMin);
   fChain->SetBranchAddress("njets30", &njets30, &b_njets30);
   fChain->SetBranchAddress("numTrueInteractions", &numTrueInteractions, &b_numTrueInteractions);
   fChain->SetBranchAddress("pfmet", &pfmet, &b_pfmet);
   fChain->SetBranchAddress("pho1dr", &pho1dr, &b_pho1dr);
   fChain->SetBranchAddress("pho1eta", &pho1eta, &b_pho1eta);
   fChain->SetBranchAddress("pho1iso", &pho1iso, &b_pho1iso);
   fChain->SetBranchAddress("pho1phi", &pho1phi, &b_pho1phi);
   fChain->SetBranchAddress("pho1pt", &pho1pt, &b_pho1pt);
   fChain->SetBranchAddress("pho2dr", &pho2dr, &b_pho2dr);
   fChain->SetBranchAddress("pho2eta", &pho2eta, &b_pho2eta);
   fChain->SetBranchAddress("pho2iso", &pho2iso, &b_pho2iso);
   fChain->SetBranchAddress("pho2phi", &pho2phi, &b_pho2phi);
   fChain->SetBranchAddress("pho2pt", &pho2pt, &b_pho2pt);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("qll13", &qll13, &b_qll13);
   fChain->SetBranchAddress("qll14", &qll14, &b_qll14);
   fChain->SetBranchAddress("qll23", &qll23, &b_qll23);
   fChain->SetBranchAddress("qll24", &qll24, &b_qll24);
   fChain->SetBranchAddress("rap", &rap, &b_rap);
   fChain->SetBranchAddress("recoVertices", &recoVertices, &b_recoVertices);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("worsePairCombRelIsoBaseline", &worsePairCombRelIsoBaseline, &b_worsePairCombRelIsoBaseline);
   fChain->SetBranchAddress("z1eta", &z1eta, &b_z1eta);
   fChain->SetBranchAddress("z1mass", &z1mass, &b_z1mass);
   fChain->SetBranchAddress("z1mll", &z1mll, &b_z1mll);
   fChain->SetBranchAddress("z1pt", &z1pt, &b_z1pt);
   fChain->SetBranchAddress("z1rap", &z1rap, &b_z1rap);
   fChain->SetBranchAddress("z2eta", &z2eta, &b_z2eta);
   fChain->SetBranchAddress("z2mass", &z2mass, &b_z2mass);
   fChain->SetBranchAddress("z2mll", &z2mll, &b_z2mll);
   fChain->SetBranchAddress("z2pt", &z2pt, &b_z2pt);
   fChain->SetBranchAddress("z2rap", &z2rap, &b_z2rap);
   fChain->SetBranchAddress("fourMassCut4Any", &fourMassCut4Any, &b_fourMassCut4Any);
   fChain->SetBranchAddress("fourMassCut4AnyS", &fourMassCut4AnyS, &b_fourMassCut4AnyS);
   fChain->SetBranchAddress("fourMassCut4SF", &fourMassCut4SF, &b_fourMassCut4SF);
   fChain->SetBranchAddress("mc_4l", &mc_4l, &b_mc_4l);
   fChain->SetBranchAddress("mc_4lany", &mc_4lany, &b_mc_4lany);
   fChain->SetBranchAddress("mc_l1", &mc_l1, &b_mc_l1);
   fChain->SetBranchAddress("mc_l2", &mc_l2, &b_mc_l2);
   fChain->SetBranchAddress("mc_l3", &mc_l3, &b_mc_l3);
   fChain->SetBranchAddress("mc_l4", &mc_l4, &b_mc_l4);
   fChain->SetBranchAddress("mc_z1", &mc_z1, &b_mc_z1);
   fChain->SetBranchAddress("mc_z2", &mc_z2, &b_mc_z2);
   fChain->SetBranchAddress("threeMassCut12SF", &threeMassCut12SF, &b_threeMassCut12SF);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("SuperMELA", &SuperMELA, &b_SuperMELA);
}

void TemplateClass::Branch(TTree *tree)
{
   TTree *fChain = tree;

   fChain->Branch("bdtScalarVsBkg_125", &bdtScalarVsBkg_125, "bdtScalarVsBkg_125/F");
   fChain->Branch("channel", &channel, "channel/F");
   fChain->Branch("elPtMin", &elPtMin, "elPtMin/F");
   fChain->Branch("genhiggsmass", &genhiggsmass, "genhiggsmass/F");
   fChain->Branch("genhiggsmassweight", &genhiggsmassweight, "genhiggsmassweight/F");
   fChain->Branch("hypo", &hypo, "hypo/F");
   fChain->Branch("intimeSimVertices", &intimeSimVertices, "intimeSimVertices/F");
   fChain->Branch("jet1eta", &jet1eta, "jet1eta/F");
   fChain->Branch("jet1phi", &jet1phi, "jet1phi/F");
   fChain->Branch("jet1pt", &jet1pt, "jet1pt/F");
   fChain->Branch("jet2eta", &jet2eta, "jet2eta/F");
   fChain->Branch("jet2phi", &jet2phi, "jet2phi/F");
   fChain->Branch("jet2pt", &jet2pt, "jet2pt/F");
   fChain->Branch("l1bdtIso", &l1bdtIso, "l1bdtIso/F");
   fChain->Branch("l1dz", &l1dz, "l1dz/F");
   fChain->Branch("l1eta", &l1eta, "l1eta/F");
   fChain->Branch("l1idNew", &l1idNew, "l1idNew/F");
   fChain->Branch("l1idPRL", &l1idPRL, "l1idPRL/F");
   fChain->Branch("l1ip2d", &l1ip2d, "l1ip2d/F");
   fChain->Branch("l1pdgId", &l1pdgId, "l1pdgId/F");
   fChain->Branch("l1pfIsoChHad04", &l1pfIsoChHad04, "l1pfIsoChHad04/F");
   fChain->Branch("l1pfIsoComb04EACorr", &l1pfIsoComb04EACorr, "l1pfIsoComb04EACorr/F");
   fChain->Branch("l1phi", &l1phi, "l1phi/F");
   fChain->Branch("l1pt", &l1pt, "l1pt/F");
   fChain->Branch("l1q", &l1q, "l1q/F");
   fChain->Branch("l1sip3d", &l1sip3d, "l1sip3d/F");
   fChain->Branch("l1trig", &l1trig, "l1trig/F");
   fChain->Branch("l2bdtIso", &l2bdtIso, "l2bdtIso/F");
   fChain->Branch("l2dz", &l2dz, "l2dz/F");
   fChain->Branch("l2eta", &l2eta, "l2eta/F");
   fChain->Branch("l2idNew", &l2idNew, "l2idNew/F");
   fChain->Branch("l2idPRL", &l2idPRL, "l2idPRL/F");
   fChain->Branch("l2ip2d", &l2ip2d, "l2ip2d/F");
   fChain->Branch("l2pdgId", &l2pdgId, "l2pdgId/F");
   fChain->Branch("l2pfIsoChHad04", &l2pfIsoChHad04, "l2pfIsoChHad04/F");
   fChain->Branch("l2pfIsoComb04EACorr", &l2pfIsoComb04EACorr, "l2pfIsoComb04EACorr/F");
   fChain->Branch("l2phi", &l2phi, "l2phi/F");
   fChain->Branch("l2pt", &l2pt, "l2pt/F");
   fChain->Branch("l2q", &l2q, "l2q/F");
   fChain->Branch("l2sip3d", &l2sip3d, "l2sip3d/F");
   fChain->Branch("l2trig", &l2trig, "l2trig/F");
   fChain->Branch("l3bdtIso", &l3bdtIso, "l3bdtIso/F");
   fChain->Branch("l3dz", &l3dz, "l3dz/F");
   fChain->Branch("l3eta", &l3eta, "l3eta/F");
   fChain->Branch("l3idNew", &l3idNew, "l3idNew/F");
   fChain->Branch("l3idPRL", &l3idPRL, "l3idPRL/F");
   fChain->Branch("l3ip2d", &l3ip2d, "l3ip2d/F");
   fChain->Branch("l3pdgId", &l3pdgId, "l3pdgId/F");
   fChain->Branch("l3pfIsoChHad04", &l3pfIsoChHad04, "l3pfIsoChHad04/F");
   fChain->Branch("l3pfIsoComb04EACorr", &l3pfIsoComb04EACorr, "l3pfIsoComb04EACorr/F");
   fChain->Branch("l3phi", &l3phi, "l3phi/F");
   fChain->Branch("l3pt", &l3pt, "l3pt/F");
   fChain->Branch("l3q", &l3q, "l3q/F");
   fChain->Branch("l3sip3d", &l3sip3d, "l3sip3d/F");
   fChain->Branch("l3trig", &l3trig, "l3trig/F");
   fChain->Branch("l4bdtIso", &l4bdtIso, "l4bdtIso/F");
   fChain->Branch("l4dz", &l4dz, "l4dz/F");
   fChain->Branch("l4eta", &l4eta, "l4eta/F");
   fChain->Branch("l4idNew", &l4idNew, "l4idNew/F");
   fChain->Branch("l4idPRL", &l4idPRL, "l4idPRL/F");
   fChain->Branch("l4ip2d", &l4ip2d, "l4ip2d/F");
   fChain->Branch("l4pdgId", &l4pdgId, "l4pdgId/F");
   fChain->Branch("l4pfIsoChHad04", &l4pfIsoChHad04, "l4pfIsoChHad04/F");
   fChain->Branch("l4pfIsoComb04EACorr", &l4pfIsoComb04EACorr, "l4pfIsoComb04EACorr/F");
   fChain->Branch("l4phi", &l4phi, "l4phi/F");
   fChain->Branch("l4pt", &l4pt, "l4pt/F");
   fChain->Branch("l4q", &l4q, "l4q/F");
   fChain->Branch("l4sip3d", &l4sip3d, "l4sip3d/F");
   fChain->Branch("l4trig", &l4trig, "l4trig/F");
   fChain->Branch("m4l", &m4l, "m4l/F");
   fChain->Branch("mass", &mass, "mass/F");
   fChain->Branch("massErr", &massErr, "massErr/F");
   fChain->Branch("melaCosTheta1", &melaCosTheta1, "melaCosTheta1/F");
   fChain->Branch("melaCosTheta2", &melaCosTheta2, "melaCosTheta2/F");
   fChain->Branch("melaCosThetaStar", &melaCosThetaStar, "melaCosThetaStar/F");
   fChain->Branch("melaLD", &melaLD, "melaLD/F");
   fChain->Branch("melaPSLD", &melaPSLD, "melaPSLD/F");
   fChain->Branch("melaPhi", &melaPhi, "melaPhi/F");
   fChain->Branch("melaPhiStar1", &melaPhiStar1, "melaPhiStar1/F");
   fChain->Branch("melaPt", &melaPt, "melaPt/F");
   fChain->Branch("melaPtY", &melaPtY, "melaPtY/F");
   fChain->Branch("melaSpinOneEven", &melaSpinOneEven, "melaSpinOneEven/F");
   fChain->Branch("melaSpinOneOdd", &melaSpinOneOdd, "melaSpinOneOdd/F");
   fChain->Branch("melaSpinTwoMinimal", &melaSpinTwoMinimal, "melaSpinTwoMinimal/F");
   fChain->Branch("melaY", &melaY, "melaY/F");
   fChain->Branch("mjj", &mjj, "mjj/F");
   fChain->Branch("mll13", &mll13, "mll13/F");
   fChain->Branch("mll14", &mll14, "mll14/F");
   fChain->Branch("mll23", &mll23, "mll23/F");
   fChain->Branch("mll24", &mll24, "mll24/F");
   fChain->Branch("muPtMin", &muPtMin, "muPtMin/F");
   fChain->Branch("njets30", &njets30, "njets30/F");
   fChain->Branch("numTrueInteractions", &numTrueInteractions, "numTrueInteractions/F");
   fChain->Branch("pfmet", &pfmet, "pfmet/F");
   fChain->Branch("pho1dr", &pho1dr, "pho1dr/F");
   fChain->Branch("pho1eta", &pho1eta, "pho1eta/F");
   fChain->Branch("pho1iso", &pho1iso, "pho1iso/F");
   fChain->Branch("pho1phi", &pho1phi, "pho1phi/F");
   fChain->Branch("pho1pt", &pho1pt, "pho1pt/F");
   fChain->Branch("pho2dr", &pho2dr, "pho2dr/F");
   fChain->Branch("pho2eta", &pho2eta, "pho2eta/F");
   fChain->Branch("pho2iso", &pho2iso, "pho2iso/F");
   fChain->Branch("pho2phi", &pho2phi, "pho2phi/F");
   fChain->Branch("pho2pt", &pho2pt, "pho2pt/F");
   fChain->Branch("pt", &pt, "pt/F");
   fChain->Branch("qll13", &qll13, "qll13/F");
   fChain->Branch("qll14", &qll14, "qll14/F");
   fChain->Branch("qll23", &qll23, "qll23/F");
   fChain->Branch("qll24", &qll24, "qll24/F");
   fChain->Branch("rap", &rap, "rap/F");
   fChain->Branch("recoVertices", &recoVertices, "recoVertices/F");
   fChain->Branch("rho", &rho, "rho/F");
   fChain->Branch("worsePairCombRelIsoBaseline", &worsePairCombRelIsoBaseline, "worsePairCombRelIsoBaseline/F");
   fChain->Branch("z1eta", &z1eta, "z1eta/F");
   fChain->Branch("z1mass", &z1mass, "z1mass/F");
   fChain->Branch("z1mll", &z1mll, "z1mll/F");
   fChain->Branch("z1pt", &z1pt, "z1pt/F");
   fChain->Branch("z1rap", &z1rap, "z1rap/F");
   fChain->Branch("z2eta", &z2eta, "z2eta/F");
   fChain->Branch("z2mass", &z2mass, "z2mass/F");
   fChain->Branch("z2mll", &z2mll, "z2mll/F");
   fChain->Branch("z2pt", &z2pt, "z2pt/F");
   fChain->Branch("z2rap", &z2rap, "z2rap/F");
   fChain->Branch("fourMassCut4Any", &fourMassCut4Any, "fourMassCut4Any/I");
   fChain->Branch("fourMassCut4AnyS", &fourMassCut4AnyS, "fourMassCut4AnyS/I");
   fChain->Branch("fourMassCut4SI", &fourMassCut4SF, "fourMassCut4SF/I");
   fChain->Branch("mc_4l", &mc_4l, "mc_4l/I");
   fChain->Branch("mc_4lany", &mc_4lany, "mc_4lany/I");
   fChain->Branch("mc_l1", &mc_l1, "mc_l1/I");
   fChain->Branch("mc_l2", &mc_l2, "mc_l2/I");
   fChain->Branch("mc_l3", &mc_l3, "mc_l3/I");
   fChain->Branch("mc_l4", &mc_l4, "mc_l4/I");
   fChain->Branch("mc_z1", &mc_z1, "mc_z1/I");
   fChain->Branch("mc_z2", &mc_z2, "mc_z2/I");
   fChain->Branch("threeMassCut12SI", &threeMassCut12SF, "threeMassCut12SF/I");
   fChain->Branch("run", &run, "run/i");
   fChain->Branch("lumi", &lumi, "lumi/i");
   fChain->Branch("event", &event, "event/i");
   fChain->Branch("SuperMELA", &SuperMELA, "SuperMELA/D");
}

#endif // #ifdef X_cxx
