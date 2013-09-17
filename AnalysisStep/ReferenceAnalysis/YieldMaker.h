#ifndef YIELDMAKER_H
#define YIELDMAKER_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <TH1D.h>
#include <TH2D.h>

#include "LeptonOrdering.h"
#include "CITHZZ/CommonCode/CommonDefs.hh"
using namespace RooFit;

class YieldMaker {

	protected :

		RooRealVar rrchannel  ;
		RooRealVar rrz1mass   ;
		RooRealVar rrz2mass   ;
		RooRealVar rrmass     ;
		RooRealVar rrmela     ;
		RooRealVar rrweight   ;
		RooRealVar rrweighterr;
		RooRealVar rrl1eta    ;
		RooRealVar rrl1phi    ;
		RooRealVar rrl1pt     ;
		RooRealVar rrl1pdgId  ;
		RooRealVar rrl2eta    ;
		RooRealVar rrl2phi    ;
		RooRealVar rrl2pt     ;
		RooRealVar rrl2pdgId  ;
		RooRealVar rrl3eta    ;
		RooRealVar rrl3phi    ;
		RooRealVar rrl3pt     ;
		RooRealVar rrl3pdgId  ;
		RooRealVar rrl4eta    ;
		RooRealVar rrl4phi    ;
		RooRealVar rrl4pt     ;
		RooRealVar rrl4pdgId  ;
		RooArgSet argset      ;
		RooArgSet argset_l1   ;			// ArgSets to hold the lepton variables
		RooArgSet argset_l2   ;
		RooArgSet argset_l3   ;
		RooArgSet argset_l4   ;
		RooDataSet dataset    ;
		RooDataSet dataset_l1 ;			// DataSets to hold the lepton variables
		RooDataSet dataset_l2 ;
		RooDataSet dataset_l3 ;
		RooDataSet dataset_l4 ;

	public :        

		YieldMaker():
			rrchannel  (RooRealVar("channel",   "channel",   0., 10.)), 
			rrz1mass   (RooRealVar("z1mass",    "z1mass",    0., 200.)), 
			rrz2mass   (RooRealVar("z2mass",    "z2mass",    0., 200.)), 
			rrmass     (RooRealVar("mass",      "mass",      0., 1000000.)),
			rrmela     (RooRealVar("mela",      "mela",      0., 1.)),
			rrweight   (RooRealVar("weight",    "weight",    -10., 10.)),
			rrweighterr(RooRealVar("weighterr", "weighterr", 0., 1000000.)),

			// Lepton parameters
			rrl1eta    (RooRealVar("l1eta",     "l1eta",     -5., 5.)),
			rrl1phi    (RooRealVar("l1phi",     "l1phi",     -5., 5.)),
			rrl1pt     (RooRealVar("l1pt",     "l1pt",      0., 1000.)),
			rrl1pdgId  (RooRealVar("l1pdgId",   "l1pdgId",   -100., 100.)),
			rrl2eta    (RooRealVar("l2eta",     "l2eta",     -5., 5.)),
			rrl2phi    (RooRealVar("l2phi",     "l2phi",     -5., 5.)),
			rrl2pt     (RooRealVar("l2pt",     "l2pt",      0., 1000.)),
			rrl2pdgId  (RooRealVar("l2pdgId",   "l2pdgId",   -100., 100.)),
			rrl3eta    (RooRealVar("l3eta",     "l3eta",     -5., 5.)),
			rrl3phi    (RooRealVar("l3phi",     "l3phi",     -5., 5.)),
			rrl3pt     (RooRealVar("l3pt",     "l3pt",      0., 1000.)),
			rrl3pdgId  (RooRealVar("l3pdgId",   "l3pdgId",   -100., 100.)),
			rrl4eta    (RooRealVar("l4eta",     "l4eta",     -5., 5.)),
			rrl4phi    (RooRealVar("l4phi",     "l4phi",     -5., 5.)),
			rrl4pt     (RooRealVar("l4pt",     "l4pt",      0., 1000.)),
			rrl4pdgId  (RooRealVar("l4pdgId",   "l4pdgId",   -100., 100.)),

			argset(RooArgSet(rrchannel, rrz1mass, rrz2mass, rrmass, rrmela, rrweight, rrweighterr, "argset")),
			argset_l1(RooArgSet(rrl1eta, rrl1phi, rrl1pt, rrl1pdgId, "argset_l1")),
			argset_l2(RooArgSet(rrl2eta, rrl2phi, rrl2pt, rrl2pdgId, "argset_l2")),
			argset_l3(RooArgSet(rrl3eta, rrl3phi, rrl3pt, rrl3pdgId, "argset_l3")),
			argset_l4(RooArgSet(rrl4eta, rrl4phi, rrl4pt, rrl4pdgId, "argset_l4")),
			dataset(RooDataSet("dataset", "dataset", argset)),
			dataset_l1(RooDataSet("dataset_l1", "dataset_l1", argset_l1)),
			dataset_l2(RooDataSet("dataset_l2", "dataset_l2", argset_l2)),
			dataset_l3(RooDataSet("dataset_l3", "dataset_l3", argset_l3)),
			dataset_l4(RooDataSet("dataset_l4", "dataset_l4", argset_l4)) {}


		void get1DHist_zmass_pairing(int channel, float z1min, float z1max, float z2min, float z2max, float m4lmin, float m4lmax, float melacut, TH1* z1hist, TH1* z2hist, int option) {

		  // Pairs up the leptons according to the given option, and then fills up z1hist and z2hist according to the invariant masses it calculates

		  for (int i = 0; i < dataset.numEntries(); i++) {
		    float mass      = dataset.get(i)->getRealValue("mass");
		    float mela      = dataset.get(i)->getRealValue("mela");
		    float weight    = dataset.get(i)->getRealValue("weight");
		    float ch        = dataset.get(i)->getRealValue("channel");

		    // Lepton parameters
		    float l1eta     = dataset_l1.get(i)->getRealValue("l1eta");
		    float l1phi     = dataset_l1.get(i)->getRealValue("l1phi");
		    float l1pt      = dataset_l1.get(i)->getRealValue("l1pt");
		    int   l1id      = (int) dataset_l1.get(i)->getRealValue("l1pdgId");
		    float l2eta     = dataset_l2.get(i)->getRealValue("l2eta");
		    float l2phi     = dataset_l2.get(i)->getRealValue("l2phi");
		    float l2pt      = dataset_l2.get(i)->getRealValue("l2pt");
		    int   l2id      = (int) dataset_l2.get(i)->getRealValue("l2pdgId");
		    float l3eta     = dataset_l3.get(i)->getRealValue("l3eta");
		    float l3phi     = dataset_l3.get(i)->getRealValue("l3phi");
		    float l3pt      = dataset_l3.get(i)->getRealValue("l3pt");
		    int   l3id      = (int) dataset_l3.get(i)->getRealValue("l3pdgId");
		    float l4eta     = dataset_l4.get(i)->getRealValue("l4eta");
		    float l4phi     = dataset_l4.get(i)->getRealValue("l4phi");
		    float l4pt      = dataset_l4.get(i)->getRealValue("l4pt");
		    int   l4id      = (int) dataset_l4.get(i)->getRealValue("l4pdgId");

		    // Creating TLorentzVectors to hold the lepton parameters
		    TLorentzVector l1, l2, l3, l4;
		    if (abs(l1id) == 11) l1.SetPtEtaPhiM(l1pt, l1eta, l1phi, ELECTRONMASS);
		    else if (abs(l1id) == 13) l1.SetPtEtaPhiM(l1pt, l1eta, l1phi, MUONMASS);
		    if (abs(l2id) == 11) l2.SetPtEtaPhiM(l2pt, l2eta, l2phi, ELECTRONMASS);
		    else if (abs(l2id) == 13) l2.SetPtEtaPhiM(l2pt, l2eta, l2phi, MUONMASS);
		    if (abs(l3id) == 11) l3.SetPtEtaPhiM(l3pt, l3eta, l3phi, ELECTRONMASS);
		    else if (abs(l3id) == 13) l3.SetPtEtaPhiM(l3pt, l3eta, l3phi, MUONMASS);
		    if (abs(l4id) == 11) l4.SetPtEtaPhiM(l4pt, l4eta, l4phi, ELECTRONMASS);
		    else if (abs(l4id) == 13) l4.SetPtEtaPhiM(l4pt, l4eta, l4phi, MUONMASS);

		    // Testing for lepton ID, and seeing whether it agrees with the correct signature
		    if ((l1id + l2id != 0) || (l3id + l4id != 0)) continue;

		    // Ordering leptons based on the given option
		    std::vector<TLorentzVector> leptons(4);
		    leptons[0] = l1; leptons[1] = l2; leptons[2] = l3; leptons[3] = l4;
		    std::vector<Int_t> leptonIDs(4);
		    leptonIDs[0] = l1id; leptonIDs[1] = l2id; leptonIDs[2] = l3id; leptonIDs[3] = l4id; 
		    order_lepton_list(leptons, leptonIDs, option);

		    // Updating channel information
		    if (abs(leptonIDs[0]) == 13 && abs(leptonIDs[2]) == 13) ch = 0;
		    if (abs(leptonIDs[0]) == 11 && abs(leptonIDs[2]) == 11) ch = 1;
		    if (abs(leptonIDs[0]) == 11 && abs(leptonIDs[2]) == 13) ch = 2;
		    if (abs(leptonIDs[0]) == 13 && abs(leptonIDs[2]) == 11) ch = 3;

		    // Filling in the masses
		    float z1mass = (leptons[0] + leptons[1]).M();
		    float z2mass = (leptons[2] + leptons[3]).M();

		    if (channel == (int)ch && z1mass>z1min && z1mass<z1max && z2mass>z2min && z2mass<z2max && mass>m4lmin && mass<m4lmax && mela>melacut) {
		      z1hist->Fill(z1mass, weight);
		      z2hist->Fill(z2mass, weight);
		    }
		  }
		}

		void get1DHist_all_zmasses(int channel, float z1min, float z1max, float z2min, float z2max, float m4lmin, float m4lmax, float melacut, TH1* zhist) {

		  // Fills zhist up with the invariant mass of all valid lepton pairs

		  for (int i = 0; i < dataset.numEntries(); i++) {
		    float mass      = dataset.get(i)->getRealValue("mass");
		    float mela      = dataset.get(i)->getRealValue("mela");
		    float weight    = dataset.get(i)->getRealValue("weight");
		    float ch        = dataset.get(i)->getRealValue("channel");

		    // Lepton parameters
		    float l1eta     = dataset_l1.get(i)->getRealValue("l1eta");
		    float l1phi     = dataset_l1.get(i)->getRealValue("l1phi");
		    float l1pt      = dataset_l1.get(i)->getRealValue("l1pt");
		    int   l1id      = (int) dataset_l1.get(i)->getRealValue("l1pdgId");
		    float l2eta     = dataset_l2.get(i)->getRealValue("l2eta");
		    float l2phi     = dataset_l2.get(i)->getRealValue("l2phi");
		    float l2pt      = dataset_l2.get(i)->getRealValue("l2pt");
		    int   l2id      = (int) dataset_l2.get(i)->getRealValue("l2pdgId");
		    float l3eta     = dataset_l3.get(i)->getRealValue("l3eta");
		    float l3phi     = dataset_l3.get(i)->getRealValue("l3phi");
		    float l3pt      = dataset_l3.get(i)->getRealValue("l3pt");
		    int   l3id      = (int) dataset_l3.get(i)->getRealValue("l3pdgId");
		    float l4eta     = dataset_l4.get(i)->getRealValue("l4eta");
		    float l4phi     = dataset_l4.get(i)->getRealValue("l4phi");
		    float l4pt      = dataset_l4.get(i)->getRealValue("l4pt");
		    int   l4id      = (int) dataset_l4.get(i)->getRealValue("l4pdgId");

		    // Creating TLorentzVectors to hold the lepton parameters
		    TLorentzVector l1, l2, l3, l4;
		    if (abs(l1id) == 11) l1.SetPtEtaPhiM(l1pt, l1eta, l1phi, ELECTRONMASS);
		    else if (abs(l1id) == 13) l1.SetPtEtaPhiM(l1pt, l1eta, l1phi, MUONMASS);
		    if (abs(l2id) == 11) l2.SetPtEtaPhiM(l2pt, l2eta, l2phi, ELECTRONMASS);
		    else if (abs(l2id) == 13) l2.SetPtEtaPhiM(l2pt, l2eta, l2phi, MUONMASS);
		    if (abs(l3id) == 11) l3.SetPtEtaPhiM(l3pt, l3eta, l3phi, ELECTRONMASS);
		    else if (abs(l3id) == 13) l3.SetPtEtaPhiM(l3pt, l3eta, l3phi, MUONMASS);
		    if (abs(l4id) == 11) l4.SetPtEtaPhiM(l4pt, l4eta, l4phi, ELECTRONMASS);
		    else if (abs(l4id) == 13) l4.SetPtEtaPhiM(l4pt, l4eta, l4phi, MUONMASS);

		    if (/*z1mass>z1min && z1mass<z1max && z2mass>z2min && z2mass<z2max &&*/ mass>m4lmin && mass<m4lmax && mela>melacut) {
		      // Going over every all possible pairs
		      if (l1id + l2id == 0) zhist->Fill( (l1 + l2).M() , weight);
		      if (l3id + l4id == 0) zhist->Fill( (l3 + l4).M() , weight);
		      if (l1id + l3id == 0) zhist->Fill( (l1 + l3).M() , weight);
		      if (l2id + l4id == 0) zhist->Fill( (l2 + l4).M() , weight);
		      if (l1id + l4id == 0) zhist->Fill( (l1 + l4).M() , weight);
		      if (l2id + l3id == 0) zhist->Fill( (l2 + l3).M() , weight);
		    }
		  }
		}

		void get1DHist_deltaR_pairing(int channel, float z1min, float z1max, float z2min, float z2max, float m4lmin, float m4lmax, float melacut, TH1* dr12hist, TH1* dr34hist, int option) {

		  // Pairs up the leptons according to the given option, and then fills dr12hist and dr34hist
		  // respectively with the delta R of the pair (lep1, lep2) and (lep3, lep4)

		  for (int i = 0; i < dataset.numEntries(); i++) {
		    float mass      = dataset.get(i)->getRealValue("mass");
		    float mela      = dataset.get(i)->getRealValue("mela");
		    float weight    = dataset.get(i)->getRealValue("weight");
		    float ch        = dataset.get(i)->getRealValue("channel");

		    // Lepton parameters
		    float l1eta     = dataset_l1.get(i)->getRealValue("l1eta");
		    float l1phi     = dataset_l1.get(i)->getRealValue("l1phi");
		    float l1pt      = dataset_l1.get(i)->getRealValue("l1pt");
		    int   l1id      = (int) dataset_l1.get(i)->getRealValue("l1pdgId");
		    float l2eta     = dataset_l2.get(i)->getRealValue("l2eta");
		    float l2phi     = dataset_l2.get(i)->getRealValue("l2phi");
		    float l2pt      = dataset_l2.get(i)->getRealValue("l2pt");
		    int   l2id      = (int) dataset_l2.get(i)->getRealValue("l2pdgId");
		    float l3eta     = dataset_l3.get(i)->getRealValue("l3eta");
		    float l3phi     = dataset_l3.get(i)->getRealValue("l3phi");
		    float l3pt      = dataset_l3.get(i)->getRealValue("l3pt");
		    int   l3id      = (int) dataset_l3.get(i)->getRealValue("l3pdgId");
		    float l4eta     = dataset_l4.get(i)->getRealValue("l4eta");
		    float l4phi     = dataset_l4.get(i)->getRealValue("l4phi");
		    float l4pt      = dataset_l4.get(i)->getRealValue("l4pt");
		    int   l4id      = (int) dataset_l4.get(i)->getRealValue("l4pdgId");

		    // Creating TLorentzVectors to hold the lepton parameters
		    TLorentzVector l1, l2, l3, l4;
		    if (abs(l1id) == 11) l1.SetPtEtaPhiM(l1pt, l1eta, l1phi, ELECTRONMASS);
		    else if (abs(l1id) == 13) l1.SetPtEtaPhiM(l1pt, l1eta, l1phi, MUONMASS);
		    if (abs(l2id) == 11) l2.SetPtEtaPhiM(l2pt, l2eta, l2phi, ELECTRONMASS);
		    else if (abs(l2id) == 13) l2.SetPtEtaPhiM(l2pt, l2eta, l2phi, MUONMASS);
		    if (abs(l3id) == 11) l3.SetPtEtaPhiM(l3pt, l3eta, l3phi, ELECTRONMASS);
		    else if (abs(l3id) == 13) l3.SetPtEtaPhiM(l3pt, l3eta, l3phi, MUONMASS);
		    if (abs(l4id) == 11) l4.SetPtEtaPhiM(l4pt, l4eta, l4phi, ELECTRONMASS);
		    else if (abs(l4id) == 13) l4.SetPtEtaPhiM(l4pt, l4eta, l4phi, MUONMASS);

		    // Testing for lepton ID, and seeing whether it agrees with the correct signature
		    if ((l1id + l2id != 0) || (l3id + l4id != 0)) continue;

		    // Ordering leptons based on the given option
		    std::vector<TLorentzVector> leptons(4);
		    leptons[0] = l1; leptons[1] = l2; leptons[2] = l3; leptons[3] = l4;
		    std::vector<Int_t> leptonIDs(4);
		    leptonIDs[0] = l1id; leptonIDs[1] = l2id; leptonIDs[2] = l3id; leptonIDs[3] = l4id; 
		    order_lepton_list(leptons, leptonIDs, option);

		    // Updating channel information
		    if (abs(leptonIDs[0]) == 13 && abs(leptonIDs[2]) == 13) ch = 0;
		    if (abs(leptonIDs[0]) == 11 && abs(leptonIDs[2]) == 11) ch = 1;
		    if (abs(leptonIDs[0]) == 11 && abs(leptonIDs[2]) == 13) ch = 2;
		    if (abs(leptonIDs[0]) == 13 && abs(leptonIDs[2]) == 11) ch = 3;

		    float z1mass = (leptons[0] + leptons[1]).M();
		    float z2mass = (leptons[2] + leptons[3]).M();

		    if (channel == (int)ch && z1mass>z1min && z1mass<z1max && z2mass>z2min && z2mass<z2max && mass>m4lmin && mass<m4lmax && mela>melacut) {
		      // Filling histograms with delta R values
		      dr12hist->Fill( (leptons[0]).DeltaR(leptons[1]) , weight);
		      dr34hist->Fill( (leptons[2]).DeltaR(leptons[3]) , weight);
		    }
		  }
		}

		void get2DHist_deltaR_pairing(int channel, float z1min, float z1max, float z2min, float z2max, float m4lmin, float m4lmax, float melacut, TH2* drhist, int option) {

		  // Pairs up the leptons according to the given option, and then fills drhist
		  // respectively with the delta R of the pair (lep1, lep2) and (lep3, lep4)

		  for (int i = 0; i < dataset.numEntries(); i++) {
		    float mass      = dataset.get(i)->getRealValue("mass");
		    float mela      = dataset.get(i)->getRealValue("mela");
		    float weight    = dataset.get(i)->getRealValue("weight");
		    float ch        = dataset.get(i)->getRealValue("channel");

		    // Lepton parameters
		    float l1eta     = dataset_l1.get(i)->getRealValue("l1eta");
		    float l1phi     = dataset_l1.get(i)->getRealValue("l1phi");
		    float l1pt      = dataset_l1.get(i)->getRealValue("l1pt");
		    int   l1id      = (int) dataset_l1.get(i)->getRealValue("l1pdgId");
		    float l2eta     = dataset_l2.get(i)->getRealValue("l2eta");
		    float l2phi     = dataset_l2.get(i)->getRealValue("l2phi");
		    float l2pt      = dataset_l2.get(i)->getRealValue("l2pt");
		    int   l2id      = (int) dataset_l2.get(i)->getRealValue("l2pdgId");
		    float l3eta     = dataset_l3.get(i)->getRealValue("l3eta");
		    float l3phi     = dataset_l3.get(i)->getRealValue("l3phi");
		    float l3pt      = dataset_l3.get(i)->getRealValue("l3pt");
		    int   l3id      = (int) dataset_l3.get(i)->getRealValue("l3pdgId");
		    float l4eta     = dataset_l4.get(i)->getRealValue("l4eta");
		    float l4phi     = dataset_l4.get(i)->getRealValue("l4phi");
		    float l4pt      = dataset_l4.get(i)->getRealValue("l4pt");
		    int   l4id      = (int) dataset_l4.get(i)->getRealValue("l4pdgId");

		    // Creating TLorentzVectors to hold the lepton parameters
		    TLorentzVector l1, l2, l3, l4;
		    if (abs(l1id) == 11) l1.SetPtEtaPhiM(l1pt, l1eta, l1phi, ELECTRONMASS);
		    else if (abs(l1id) == 13) l1.SetPtEtaPhiM(l1pt, l1eta, l1phi, MUONMASS);
		    if (abs(l2id) == 11) l2.SetPtEtaPhiM(l2pt, l2eta, l2phi, ELECTRONMASS);
		    else if (abs(l2id) == 13) l2.SetPtEtaPhiM(l2pt, l2eta, l2phi, MUONMASS);
		    if (abs(l3id) == 11) l3.SetPtEtaPhiM(l3pt, l3eta, l3phi, ELECTRONMASS);
		    else if (abs(l3id) == 13) l3.SetPtEtaPhiM(l3pt, l3eta, l3phi, MUONMASS);
		    if (abs(l4id) == 11) l4.SetPtEtaPhiM(l4pt, l4eta, l4phi, ELECTRONMASS);
		    else if (abs(l4id) == 13) l4.SetPtEtaPhiM(l4pt, l4eta, l4phi, MUONMASS);

		    // Testing for lepton ID, and seeing whether it agrees with the correct signature
		    if ((l1id + l2id != 0) || (l3id + l4id != 0)) continue;

		    // Ordering leptons based on the given option
		    std::vector<TLorentzVector> leptons(4);
		    leptons[0] = l1; leptons[1] = l2; leptons[2] = l3; leptons[3] = l4;
		    std::vector<Int_t> leptonIDs(4);
		    leptonIDs[0] = l1id; leptonIDs[1] = l2id; leptonIDs[2] = l3id; leptonIDs[3] = l4id; 
		    order_lepton_list(leptons, leptonIDs, option);


		    // Updating channel information
		    if (abs(leptonIDs[0]) == 13 && abs(leptonIDs[2]) == 13) ch = 0;
		    if (abs(leptonIDs[0]) == 11 && abs(leptonIDs[2]) == 11) ch = 1;
		    if (abs(leptonIDs[0]) == 11 && abs(leptonIDs[2]) == 13) ch = 2;
		    if (abs(leptonIDs[0]) == 13 && abs(leptonIDs[2]) == 11) ch = 3;

		    float z1mass = (leptons[0] + leptons[1]).M();
		    float z2mass = (leptons[2] + leptons[3]).M();

		    if (channel == (int)ch && z1mass>z1min && z1mass<z1max && z2mass>z2min && z2mass<z2max && mass>m4lmin && mass<m4lmax && mela>melacut) {
		      // Filling histograms with delta R values
		      drhist->Fill( (leptons[0]).DeltaR(leptons[1]) , (leptons[2]).DeltaR(leptons[3]) , weight);
		    }
		  }
		}

		void get2DHist_zmass_pairing(int channel, float z1min, float z1max, float z2min, float z2max, float m4lmin, float m4lmax, float melacut, TH2* zhist, int option) {

		  // Pairs up the leptons according to the given option, and then fills up z1hist and z2hist according to the invariant masses it calculates

		  for (int i = 0; i < dataset.numEntries(); i++) {
		    float mass      = dataset.get(i)->getRealValue("mass");
		    float mela      = dataset.get(i)->getRealValue("mela");
		    float weight    = dataset.get(i)->getRealValue("weight");
		    float ch        = dataset.get(i)->getRealValue("channel");

		    // Lepton parameters
		    float l1eta     = dataset_l1.get(i)->getRealValue("l1eta");
		    float l1phi     = dataset_l1.get(i)->getRealValue("l1phi");
		    float l1pt      = dataset_l1.get(i)->getRealValue("l1pt");
		    int   l1id      = (int) dataset_l1.get(i)->getRealValue("l1pdgId");
		    float l2eta     = dataset_l2.get(i)->getRealValue("l2eta");
		    float l2phi     = dataset_l2.get(i)->getRealValue("l2phi");
		    float l2pt      = dataset_l2.get(i)->getRealValue("l2pt");
		    int   l2id      = (int) dataset_l2.get(i)->getRealValue("l2pdgId");
		    float l3eta     = dataset_l3.get(i)->getRealValue("l3eta");
		    float l3phi     = dataset_l3.get(i)->getRealValue("l3phi");
		    float l3pt      = dataset_l3.get(i)->getRealValue("l3pt");
		    int   l3id      = (int) dataset_l3.get(i)->getRealValue("l3pdgId");
		    float l4eta     = dataset_l4.get(i)->getRealValue("l4eta");
		    float l4phi     = dataset_l4.get(i)->getRealValue("l4phi");
		    float l4pt      = dataset_l4.get(i)->getRealValue("l4pt");
		    int   l4id      = (int) dataset_l4.get(i)->getRealValue("l4pdgId");

		    // Creating TLorentzVectors to hold the lepton parameters
		    TLorentzVector l1, l2, l3, l4;
		    if (abs(l1id) == 11) l1.SetPtEtaPhiM(l1pt, l1eta, l1phi, ELECTRONMASS);
		    else if (abs(l1id) == 13) l1.SetPtEtaPhiM(l1pt, l1eta, l1phi, MUONMASS);
		    if (abs(l2id) == 11) l2.SetPtEtaPhiM(l2pt, l2eta, l2phi, ELECTRONMASS);
		    else if (abs(l2id) == 13) l2.SetPtEtaPhiM(l2pt, l2eta, l2phi, MUONMASS);
		    if (abs(l3id) == 11) l3.SetPtEtaPhiM(l3pt, l3eta, l3phi, ELECTRONMASS);
		    else if (abs(l3id) == 13) l3.SetPtEtaPhiM(l3pt, l3eta, l3phi, MUONMASS);
		    if (abs(l4id) == 11) l4.SetPtEtaPhiM(l4pt, l4eta, l4phi, ELECTRONMASS);
		    else if (abs(l4id) == 13) l4.SetPtEtaPhiM(l4pt, l4eta, l4phi, MUONMASS);

		    // Testing for lepton ID, and seeing whether it agrees with the correct signature
		    if ((l1id + l2id != 0) || (l3id + l4id != 0)) continue;

		    // Ordering leptons based on the given option
		    std::vector<TLorentzVector> leptons(4);
		    leptons[0] = l1; leptons[1] = l2; leptons[2] = l3; leptons[3] = l4;
		    std::vector<Int_t> leptonIDs(4);
		    leptonIDs[0] = l1id; leptonIDs[1] = l2id; leptonIDs[2] = l3id; leptonIDs[3] = l4id; 
		    order_lepton_list(leptons, leptonIDs, option);

		    // Updating channel information
		    if (abs(leptonIDs[0]) == 13 && abs(leptonIDs[2]) == 13) ch = 0;
		    if (abs(leptonIDs[0]) == 11 && abs(leptonIDs[2]) == 11) ch = 1;
		    if (abs(leptonIDs[0]) == 11 && abs(leptonIDs[2]) == 13) ch = 2;
		    if (abs(leptonIDs[0]) == 13 && abs(leptonIDs[2]) == 11) ch = 3;

		    // Filling in the masses
		    float z1mass = (leptons[0] + leptons[1]).M();
		    float z2mass = (leptons[2] + leptons[3]).M();
		    if (channel == (int)ch && z1mass>z1min && z1mass<z1max && z2mass>z2min && z2mass<z2max && mass>m4lmin && mass<m4lmax && mela>melacut) zhist->Fill(z1mass, z2mass, weight);
		  }

		}

		void get1DHist_z1mass(int channel, float z1min, float z1max, float z2min, float z2max, float m4lmin, float m4lmax, float melacut, TH1* hist) {
		  for (int i = 0; i < dataset.numEntries(); i++) {
		    float z1mass    = dataset.get(i)->getRealValue("z1mass");
		    float z2mass    = dataset.get(i)->getRealValue("z2mass");
		    float mass      = dataset.get(i)->getRealValue("mass");
		    float mela      = dataset.get(i)->getRealValue("mela");
		    float weight    = dataset.get(i)->getRealValue("weight");
		    float ch        = dataset.get(i)->getRealValue("channel");

		    // Lepton IDs
		    int   l1id      = (int) dataset_l1.get(i)->getRealValue("l1pdgId");
		    int   l2id      = (int) dataset_l2.get(i)->getRealValue("l2pdgId");
		    int   l3id      = (int) dataset_l3.get(i)->getRealValue("l3pdgId");
		    int   l4id      = (int) dataset_l4.get(i)->getRealValue("l4pdgId");

		    if (channel == (int)ch && z1mass>z1min && z1mass<z1max && z2mass>z2min && z2mass<z2max && mass>m4lmin && mass<m4lmax && mela>melacut) hist->Fill(z1mass, weight);
		  }
		}

		void get1DHist_z2mass(int channel, float z1min, float z1max, float z2min, float z2max, float m4lmin, float m4lmax, float melacut, TH1* hist) {
		  for (int i = 0; i < dataset.numEntries(); i++) {
		    float z1mass    = dataset.get(i)->getRealValue("z1mass");
		    float z2mass    = dataset.get(i)->getRealValue("z2mass");
		    float mass      = dataset.get(i)->getRealValue("mass");
		    float mela      = dataset.get(i)->getRealValue("mela");
		    float weight    = dataset.get(i)->getRealValue("weight");
		    float ch        = dataset.get(i)->getRealValue("channel");

		    // Lepton IDs
		    int   l1id      = (int) dataset_l1.get(i)->getRealValue("l1pdgId");
		    int   l2id      = (int) dataset_l2.get(i)->getRealValue("l2pdgId");
		    int   l3id      = (int) dataset_l3.get(i)->getRealValue("l3pdgId");
		    int   l4id      = (int) dataset_l4.get(i)->getRealValue("l4pdgId");


		    if (channel == (int)ch && z1mass>z1min && z1mass<z1max && z2mass>z2min && z2mass<z2max && mass>m4lmin && mass<m4lmax && mela>melacut) hist->Fill(z2mass, weight);
		  }
		}

		void get1DHist_mela(int channel, float z1min, float z1max, float z2min, float z2max, float m4lmin, float m4lmax, float melacut, TH1* hist) {
		  for (int i = 0; i < dataset.numEntries(); i++) {
		    float z1mass    = dataset.get(i)->getRealValue("z1mass");
		    float z2mass    = dataset.get(i)->getRealValue("z2mass");
		    float mass      = dataset.get(i)->getRealValue("mass");
		    float mela      = dataset.get(i)->getRealValue("mela");
		    float weight    = dataset.get(i)->getRealValue("weight");
		    float ch        = dataset.get(i)->getRealValue("channel");

		    // Lepton IDs
		    int   l1id      = (int) dataset_l1.get(i)->getRealValue("l1pdgId");
		    int   l2id      = (int) dataset_l2.get(i)->getRealValue("l2pdgId");
		    int   l3id      = (int) dataset_l3.get(i)->getRealValue("l3pdgId");
		    int   l4id      = (int) dataset_l4.get(i)->getRealValue("l4pdgId");

		    if (channel == (int)ch && z1mass>z1min && z1mass<z1max && z2mass>z2min && z2mass<z2max && mass>m4lmin && mass<m4lmax && mela>melacut) hist->Fill(mela, weight);
		  }
		}

		void get2DHist_zmass(int channel, float z1min, float z1max, float z2min, float z2max, float m4lmin, float m4lmax, float melacut, TH2* hist) {
		  // Puts the Z1 vs Z2 mass histogram into hist
		  // Z1 is on X axis, Z2 is on Y axis

		  for (int i = 0; i < dataset.numEntries(); i++) {
		    float z1mass    = dataset.get(i)->getRealValue("z1mass");
		    float z2mass    = dataset.get(i)->getRealValue("z2mass");
		    float mass      = dataset.get(i)->getRealValue("mass");
		    float mela      = dataset.get(i)->getRealValue("mela");
		    float weight    = dataset.get(i)->getRealValue("weight");
		    float ch        = dataset.get(i)->getRealValue("channel");

		    // Lepton IDs
		    int   l1id      = (int) dataset_l1.get(i)->getRealValue("l1pdgId");
		    int   l2id      = (int) dataset_l2.get(i)->getRealValue("l2pdgId");
		    int   l3id      = (int) dataset_l3.get(i)->getRealValue("l3pdgId");
		    int   l4id      = (int) dataset_l4.get(i)->getRealValue("l4pdgId");

		    if (channel == (int)ch && z1mass>z1min && z1mass<z1max && z2mass>z2min && z2mass<z2max && mass>m4lmin && mass<m4lmax && mela>melacut) hist->Fill(z1mass, z2mass, weight);
		  }

		}

		void get1DHist_mass(int channel, float z1min, float z1max, float z2min, float z2max, float m4lmin, float m4lmax, float melacut, TH1* hist) {

		  for (int i = 0; i < dataset.numEntries(); i++) {
		    float z1mass    = dataset.get(i)->getRealValue("z1mass");
		    float z2mass    = dataset.get(i)->getRealValue("z2mass");
		    float mass      = dataset.get(i)->getRealValue("mass");
		    float mela      = dataset.get(i)->getRealValue("mela");
		    float weight    = dataset.get(i)->getRealValue("weight");
		    float ch        = dataset.get(i)->getRealValue("channel");

		    if (channel == (int)ch && z1mass>z1min && z1mass<z1max && z2mass>z2min && z2mass<z2max && mass>m4lmin && mass<m4lmax && mela>melacut) hist->Fill(mass, weight);
		  }

		}

};


class DataYieldMaker : public YieldMaker {

	private :
		std::vector<std::pair<int, int> > runeventinfo;


	public:

		DataYieldMaker():YieldMaker(){}

		void fill(std::string filepath) {
		  TFile* file = new TFile(filepath.c_str());
		  TTree* tree = (TTree*)file->Get("zz4lTree/probe_tree");

		  TBranch *bchannel   = tree->GetBranch("channel");
		  TBranch *bz1mass    = tree->GetBranch("z1mass");
		  TBranch *bz2mass    = tree->GetBranch("z2mass");
		  TBranch *bmass      = tree->GetBranch("mass");
		  TBranch *bmela      = tree->GetBranch("melaLD");
		  TBranch *bevent     = tree->GetBranch("event");
		  TBranch *brun       = tree->GetBranch("run");
		  TBranch *bl3pt      = tree->GetBranch("l3pt");
		  TBranch *bl3eta     = tree->GetBranch("l3eta");
		  TBranch *bl3phi     = tree->GetBranch("l3phi");
		  TBranch *bl3id      = tree->GetBranch("l3idNew");
		  TBranch *bl3iso     = tree->GetBranch("l3pfIsoComb04EACorr");
		  TBranch *bl3pdgId   = tree->GetBranch("l3pdgId");
		  TBranch *bl4pt      = tree->GetBranch("l4pt");
		  TBranch *bl4eta     = tree->GetBranch("l4eta");
		  TBranch *bl4phi     = tree->GetBranch("l4phi");
		  TBranch *bl4id      = tree->GetBranch("l4idNew");
		  TBranch *bl4iso     = tree->GetBranch("l4pfIsoComb04EACorr");
		  TBranch *bl4pdgId   = tree->GetBranch("l4pdgId");
		  TBranch *bl1pt      = tree->GetBranch("l1pt");
		  TBranch *bl1eta     = tree->GetBranch("l1eta");
		  TBranch *bl1phi     = tree->GetBranch("l1phi");
		  TBranch *bl1pdgId   = tree->GetBranch("l1pdgId");
		  TBranch *bl2pt      = tree->GetBranch("l2pt");
		  TBranch *bl2eta     = tree->GetBranch("l2eta");
		  TBranch *bl2phi     = tree->GetBranch("l2phi");
		  TBranch *bl2pdgId   = tree->GetBranch("l2pdgId");

		  float channel   = 0.0;
		  float z1mass    = 0.0;
		  float z2mass    = 0.0;
		  float mass      = 0.0;
		  float mela      = 0.0;
		  int   event     = 0;
		  int   run       = 0;
		  float l3pt      = 0.0;
		  float l3eta     = 0.0;
		  float l3phi     = 0.0;
		  float l3id      = 0.0;
		  float l3iso     = 0.0;
		  float l3pdgId   = 0.0;
		  float l4pt      = 0.0;
		  float l4eta     = 0.0;
		  float l4phi     = 0.0;
		  float l4id      = 0.0;
		  float l4iso     = 0.0;
		  float l4pdgId   = 0.0;
		  float l1pt      = 0.0;
		  float l1eta     = 0.0;
		  float l1phi     = 0.0;
		  float l1pdgId   = 0.0;
		  float l2pt      = 0.0;
		  float l2eta     = 0.0;
		  float l2phi     = 0.0;
		  float l2pdgId   = 0.0;

		  bchannel   ->SetAddress(&channel);
		  bz1mass    ->SetAddress(&z1mass);
		  bz2mass    ->SetAddress(&z2mass);
		  bmass      ->SetAddress(&mass);
		  bmela      ->SetAddress(&mela);
		  bevent     ->SetAddress(&event);
		  brun       ->SetAddress(&run);
		  bl3pt      ->SetAddress(&l3pt);
		  bl3eta     ->SetAddress(&l3eta);
		  bl3phi     ->SetAddress(&l3phi);
		  bl3id      ->SetAddress(&l3id);
		  bl3iso     ->SetAddress(&l3iso);
		  bl3pdgId   ->SetAddress(&l3pdgId);
		  bl4pt      ->SetAddress(&l4pt);
		  bl4eta     ->SetAddress(&l4eta);
		  bl4phi     ->SetAddress(&l4phi);
		  bl4id      ->SetAddress(&l4id);
		  bl4iso     ->SetAddress(&l4iso);
		  bl4pdgId   ->SetAddress(&l4pdgId);
		  bl1pt      ->SetAddress(&l1pt);
		  bl1eta     ->SetAddress(&l1eta);
		  bl1phi     ->SetAddress(&l1phi);
		  bl1pdgId   ->SetAddress(&l1pdgId);
		  bl2pt      ->SetAddress(&l2pt);
		  bl2eta     ->SetAddress(&l2eta);
		  bl2phi     ->SetAddress(&l2phi);
		  bl2pdgId   ->SetAddress(&l2pdgId);

		  for (int i = 0; i < tree->GetEntries(); i++) {
		    bchannel   ->GetEvent(i);
		    bz1mass    ->GetEvent(i);
		    bz2mass    ->GetEvent(i);
		    bmass      ->GetEvent(i);
		    bevent     ->GetEvent(i);
		    brun       ->GetEvent(i);
		    bmela      ->GetEvent(i);
		    bl1pt      ->GetEvent(i);
		    bl1eta     ->GetEvent(i);
		    bl1phi     ->GetEvent(i);
		    bl1pdgId   ->GetEvent(i);
		    bl2pt      ->GetEvent(i);
		    bl2eta     ->GetEvent(i);
		    bl2phi     ->GetEvent(i);
		    bl2pdgId   ->GetEvent(i);
		    bl3pt      ->GetEvent(i);
		    bl3eta     ->GetEvent(i);
		    bl3phi     ->GetEvent(i);
		    bl3id      ->GetEvent(i);
		    bl3iso     ->GetEvent(i);
		    bl3pdgId   ->GetEvent(i);
		    bl4pt      ->GetEvent(i);
		    bl4eta     ->GetEvent(i);
		    bl4phi     ->GetEvent(i);
		    bl4id      ->GetEvent(i);
		    bl4iso     ->GetEvent(i);
		    bl4pdgId   ->GetEvent(i);

		    bool existsAlready = false;
		    for (std::size_t k = 0; k < runeventinfo.size(); k++) {
		      if (run == runeventinfo[k].first && event == runeventinfo[k].second) existsAlready = true;
		    }

		    if (!existsAlready) {
		      argset.setRealValue("z1mass",    z1mass);
		      argset.setRealValue("z2mass",    z2mass);
		      argset.setRealValue("mass",      mass);
		      argset.setRealValue("mela",      mela);
		      argset.setRealValue("channel",   channel);
		      argset.setRealValue("weight",    1.0);
		      argset.setRealValue("weighterr", 0.0);
		      argset_l1.setRealValue("l1eta", l1eta);
		      argset_l1.setRealValue("l1phi", l1phi);
		      argset_l1.setRealValue("l1pt", l1pt);
		      argset_l1.setRealValue("l1pdgId", l1pdgId);
		      argset_l2.setRealValue("l2eta", l2eta);
		      argset_l2.setRealValue("l2phi", l2phi);
		      argset_l2.setRealValue("l2pt", l2pt);
		      argset_l2.setRealValue("l2pdgId", l2pdgId);
		      argset_l3.setRealValue("l3eta", l3eta);
		      argset_l3.setRealValue("l3phi", l3phi);
		      argset_l3.setRealValue("l3pt", l3pt);
		      argset_l3.setRealValue("l3pdgId", l3pdgId);
		      argset_l4.setRealValue("l4eta", l4eta);
		      argset_l4.setRealValue("l4phi", l4phi);
		      argset_l4.setRealValue("l4pt", l4pt);
		      argset_l4.setRealValue("l4pdgId", l4pdgId);
		      runeventinfo.push_back(std::pair<int, int>(run, event));
		      dataset.add(argset);
		      dataset_l1.add(argset_l1);
		      dataset_l2.add(argset_l2);
		      dataset_l3.add(argset_l3);
		      dataset_l4.add(argset_l4);
		    }
		  }
		}
};


class ZXYieldMaker : public YieldMaker {

	private :
		std::vector<std::pair<int, int> > runeventinfo;    

		RooDataSet data2p2f;
		RooDataSet data3p2f;



	public:

		ZXYieldMaker():YieldMaker(){}

		void fill(std::string filepath, float wgt, FakeRateCalculator FR, bool doSS) {

		  TFile* file = new TFile(filepath.c_str());
		  TTree* tree = (TTree*)file->Get("zxTree/probe_tree");

		  TBranch *bchannel   = tree->GetBranch("channel");
		  TBranch *bz1mass    = tree->GetBranch("z1mass");
		  TBranch *bz2mass    = tree->GetBranch("z2mass");
		  TBranch *bmass      = tree->GetBranch("mass");
		  TBranch *bl3pt      = tree->GetBranch("l3pt");
		  TBranch *bl3eta     = tree->GetBranch("l3eta");
		  TBranch *bl3phi     = tree->GetBranch("l3phi");
		  TBranch *bl3id      = tree->GetBranch("l3idNew");
		  TBranch *bl3iso     = tree->GetBranch("l3pfIsoComb04EACorr");
		  TBranch *bl3pdgId   = tree->GetBranch("l3pdgId");
		  TBranch *bl4pt      = tree->GetBranch("l4pt");
		  TBranch *bl4eta     = tree->GetBranch("l4eta");
		  TBranch *bl4phi     = tree->GetBranch("l4phi");
		  TBranch *bl4id      = tree->GetBranch("l4idNew");
		  TBranch *bl4iso     = tree->GetBranch("l4pfIsoComb04EACorr");
		  TBranch *bl4pdgId   = tree->GetBranch("l4pdgId");
		  TBranch *bl1pt      = tree->GetBranch("l1pt");
		  TBranch *bl1eta     = tree->GetBranch("l1eta");
		  TBranch *bl1phi     = tree->GetBranch("l1phi");
		  TBranch *bl1pdgId   = tree->GetBranch("l1pdgId");
		  TBranch *bl2pt      = tree->GetBranch("l2pt");
		  TBranch *bl2eta     = tree->GetBranch("l2eta");
		  TBranch *bl2phi     = tree->GetBranch("l2phi");
		  TBranch *bl2pdgId   = tree->GetBranch("l2pdgId");
		  TBranch *bmela      = tree->GetBranch("melaLD");
		  TBranch *bevent     = tree->GetBranch("event");
		  TBranch *brun       = tree->GetBranch("run");

		  float channel   = 0.0;
		  float z1mass    = 0.0;
		  float z2mass    = 0.0;
		  float mass      = 0.0;
		  float l3pt      = 0.0;
		  float l3eta     = 0.0;
		  float l3phi     = 0.0;
		  float l3id      = 0.0;
		  float l3iso     = 0.0;
		  float l3pdgId   = 0.0;
		  float l4pt      = 0.0;
		  float l4eta     = 0.0;
		  float l4phi     = 0.0;
		  float l4id      = 0.0;
		  float l4iso     = 0.0;
		  float l4pdgId   = 0.0;
		  float l1pt      = 0.0;
		  float l1eta     = 0.0;
		  float l1phi     = 0.0;
		  float l1pdgId   = 0.0;
		  float l2pt      = 0.0;
		  float l2eta     = 0.0;
		  float l2phi     = 0.0;
		  float l2pdgId   = 0.0;
		  float mela      = 0.0;
		  int   event     = 0;
		  int   run       = 0;

		  bchannel   ->SetAddress(&channel); 
		  bz1mass    ->SetAddress(&z1mass);
		  bz2mass    ->SetAddress(&z2mass);
		  bmass      ->SetAddress(&mass);
		  bl3pt      ->SetAddress(&l3pt);
		  bl3eta     ->SetAddress(&l3eta);
		  bl3phi     ->SetAddress(&l3phi);
		  bl3id      ->SetAddress(&l3id);
		  bl3iso     ->SetAddress(&l3iso);
		  bl3pdgId   ->SetAddress(&l3pdgId);
		  bl4pt      ->SetAddress(&l4pt);
		  bl4eta     ->SetAddress(&l4eta);
		  bl4phi     ->SetAddress(&l4phi);
		  bl4id      ->SetAddress(&l4id);
		  bl4iso     ->SetAddress(&l4iso);
		  bl4pdgId   ->SetAddress(&l4pdgId);
		  bl1pt      ->SetAddress(&l1pt);
		  bl1eta     ->SetAddress(&l1eta);
		  bl1phi     ->SetAddress(&l1phi);
		  bl1pdgId   ->SetAddress(&l1pdgId);
		  bl2pt      ->SetAddress(&l2pt);
		  bl2eta     ->SetAddress(&l2eta);
		  bl2phi     ->SetAddress(&l2phi);
		  bl2pdgId   ->SetAddress(&l2pdgId);
		  bmela      ->SetAddress(&mela);
		  bevent     ->SetAddress(&event);
		  brun       ->SetAddress(&run);

		  for (int i = 0; i < tree->GetEntries(); i++) {
		    bchannel   ->GetEvent(i);
		    bz1mass    ->GetEvent(i);
		    bz2mass    ->GetEvent(i);
		    bmass      ->GetEvent(i);
		    bl3pt      ->GetEvent(i);
		    bl3eta     ->GetEvent(i);
		    bl3phi     ->GetEvent(i);
		    bl3id      ->GetEvent(i);
		    bl3iso     ->GetEvent(i);
		    bl3pdgId   ->GetEvent(i);
		    bl4pt      ->GetEvent(i);
		    bl4eta     ->GetEvent(i);
		    bl4phi     ->GetEvent(i);
		    bl4id      ->GetEvent(i);
		    bl4iso     ->GetEvent(i);
		    bl4pdgId   ->GetEvent(i);
		    bevent     ->GetEvent(i);
		    brun       ->GetEvent(i);
		    bl1pt      ->GetEvent(i);
		    bl1eta     ->GetEvent(i);
		    bl1phi     ->GetEvent(i);
		    bl1pdgId   ->GetEvent(i);
		    bl2pt      ->GetEvent(i);
		    bl2eta     ->GetEvent(i);
		    bl2phi     ->GetEvent(i);
		    bl2pdgId   ->GetEvent(i);
		    bmela      ->GetEvent(i);


		    bool existsAlready = false;
		    for (std::size_t k = 0; k < runeventinfo.size(); k++) {
		      if (run == runeventinfo[k].first && event == runeventinfo[k].second) existsAlready = true;
		    }

		    if (!existsAlready) {
		      argset.setRealValue("z1mass",    z1mass);
		      argset.setRealValue("z2mass",    z2mass);
		      argset.setRealValue("mass",      mass);
		      argset.setRealValue("mela",      mela);
		      argset.setRealValue("channel",   channel);
		      argset.setRealValue("weight",    1.0);
		      argset.setRealValue("weighterr", 0.0);
		      argset_l1.setRealValue("l1eta", l1eta);
		      argset_l1.setRealValue("l1phi", l1phi);
		      argset_l1.setRealValue("l1pt", l1pt);
		      argset_l1.setRealValue("l1pdgId", l1pdgId);
		      argset_l2.setRealValue("l2eta", l2eta);
		      argset_l2.setRealValue("l2phi", l2phi);
		      argset_l2.setRealValue("l2pt", l2pt);
		      argset_l2.setRealValue("l2pdgId", l2pdgId);
		      argset_l3.setRealValue("l3eta", l3eta);
		      argset_l3.setRealValue("l3phi", l3phi);
		      argset_l3.setRealValue("l3pt", l3pt);
		      argset_l3.setRealValue("l3pdgId", l3pdgId);
		      argset_l4.setRealValue("l4eta", l4eta);
		      argset_l4.setRealValue("l4phi", l4phi);
		      argset_l4.setRealValue("l4pt", l4pt);
		      argset_l4.setRealValue("l4pdgId", l4pdgId);
		      runeventinfo.push_back(std::pair<int, int>(run, event));
		      dataset_l1.add(argset_l1);
		      dataset_l2.add(argset_l2);
		      dataset_l3.add(argset_l3);
		      dataset_l4.add(argset_l4);

		      if (!doSS && l3pdgId == -l4pdgId) { 
			/*
			   float f1    = FR.getFakeRate(l3pt, l3eta, l3pdgId);
			   float f2    = FR.getFakeRate(l4pt, l4eta, l4pdgId);
			   float p1    = FR.getPromptRate(l3pt, l3eta, l3pdgId);
			   float p2    = FR.getPromptRate(l4pt, l4eta, l4pdgId);

			   float eps1  = f1/(1.0-f1);
			   float eps2  = f2/(1.0-f2);
			   float eta1  = (1.0-p1)/p1;
			   float eta2  = (1.0-p2)/p2;
			   float deno = (1-(eps1*eta1))*(1-(eps2*eta2));
			   float weight = 0.0;

			   if      ((l3id==0 || l3iso/l3pt>0.4) && (l4id==0 || l4iso/l4pt>0.4)) weight = -eps1*eps2*wgt;
			   else if ((l3id==0 || l3iso/l3pt>0.4) && (l4id==1 && l4iso/l4pt<0.4)) weight = eps1*wgt;
			   else if ((l3id==1 && l3iso/l3pt<0.4) && (l4id==0 || l4iso/l4pt>0.4)) weight = eps2*wgt;
			   else if ((l3id==1 && l3iso/l3pt<0.4) && (l4id==1 && l4iso/l4pt<0.4)) weight = (-eps1*eta1 - eps2*eta2 + eps1*eps2*eta1*eta2)*wgt;

			   weight /= deno;
			   argset.setRealValue("weight", weight);
			   argset.setRealValue("weighterr", 0.0);
			 */

			argset.setRealValue("weight", 1.0);
			argset.setRealValue("weighterr", 1.0);
			dataset.add(argset);

		      }


		      else if (doSS && l3pdgId == l4pdgId) {
			float f1    = FR.getFakeRate(l3pt, l3eta, l3pdgId);
			float f2    = FR.getFakeRate(l4pt, l4eta, l4pdgId);
			float f1_up = FR.getFakeRate(l3pt, l3eta, l3pdgId) + FR.getFakeRateErr(l3pt, l3eta, l3pdgId);
			float f2_up = FR.getFakeRate(l4pt, l4eta, l4pdgId) + FR.getFakeRateErr(l4pt, l4eta, l4pdgId);
			float f1_dn = FR.getFakeRate(l3pt, l3eta, l3pdgId) - FR.getFakeRateErr(l3pt, l3eta, l3pdgId);
			float f2_dn = FR.getFakeRate(l4pt, l4eta, l4pdgId) - FR.getFakeRateErr(l4pt, l4eta, l4pdgId);
			float sf = f1*f2;
			if (channel == 0) sf*= 1.28;
			if (channel == 1) sf*= 0.93;
			if (channel == 2 || channel == 3) sf*= 0.94;

			float sf_up = f1_up*f2_up;
			if (channel == 0) sf_up*= 1.28;
			if (channel == 1) sf_up*= 0.93;
			if (channel == 2 || channel == 3) sf_up*= 0.94;

			float sf_dn = f1_dn*f2_dn;
			if (channel == 0) sf_dn*= 1.28;
			if (channel == 1) sf_dn*= 0.93;
			if (channel == 2 || channel == 3) sf_dn*= 0.94;

			float weight    = sf*wgt;
			float weight_up = sf_up*wgt;
			float weight_dn = sf_dn*wgt;

			argset.setRealValue("weight", weight);
			argset.setRealValue("weighterr", std::max<float>(fabs(weight_up - weight), fabs(weight_dn - weight)));
			dataset.add(argset);

		      }
		    }
		  }
		}
};


class ZZYieldMaker : public YieldMaker {

	private :
		std::vector<std::pair<int, int> > runeventinfo;    

	public:

		ZZYieldMaker():YieldMaker(){}

		void fill(std::string filepath, float wgt, float wgterr, bool isSignal) {
		  if (runeventinfo.size()>0) runeventinfo.clear();       

		  TFile* file = new TFile(filepath.c_str());
		  TTree* tree = (TTree*)file->Get("zz4lTree/probe_tree");

		  TBranch *bchannel   = tree->GetBranch("channel");
		  TBranch *bz1mass    = tree->GetBranch("z1mass");
		  TBranch *bz2mass    = tree->GetBranch("z2mass");
		  TBranch *bmass      = tree->GetBranch("mass");
		  TBranch *bl3pt      = tree->GetBranch("l3pt");
		  TBranch *bl3eta     = tree->GetBranch("l3eta");
		  TBranch *bl3phi     = tree->GetBranch("l3phi");
		  TBranch *bl3id      = tree->GetBranch("l3idNew");
		  TBranch *bl3iso     = tree->GetBranch("l3pfIsoComb04EACorr");
		  TBranch *bl3pdgId   = tree->GetBranch("l3pdgId");
		  TBranch *bl4pt      = tree->GetBranch("l4pt");
		  TBranch *bl4eta     = tree->GetBranch("l4eta");
		  TBranch *bl4phi     = tree->GetBranch("l4phi");
		  TBranch *bl4id      = tree->GetBranch("l4idNew");
		  TBranch *bl4iso     = tree->GetBranch("l4pfIsoComb04EACorr");
		  TBranch *bl4pdgId   = tree->GetBranch("l4pdgId");
		  TBranch *bl1pt      = tree->GetBranch("l1pt");
		  TBranch *bl1eta     = tree->GetBranch("l1eta");
		  TBranch *bl1phi     = tree->GetBranch("l1phi");
		  TBranch *bl1pdgId   = tree->GetBranch("l1pdgId");
		  TBranch *bl2pt      = tree->GetBranch("l2pt");
		  TBranch *bl2eta     = tree->GetBranch("l2eta");
		  TBranch *bl2phi     = tree->GetBranch("l2phi");
		  TBranch *bl2pdgId   = tree->GetBranch("l2pdgId");
		  TBranch *bmela      = tree->GetBranch("melaLD");
		  TBranch *bevent     = tree->GetBranch("event");
		  TBranch *brun       = tree->GetBranch("run");
		  TBranch *bnumsim    = tree->GetBranch("numTrueInteractions");

		  float channel   = 0.0;
		  float z1mass    = 0.0;
		  float z2mass    = 0.0;
		  float mass      = 0.0;
		  float l3pt      = 0.0;
		  float l3eta     = 0.0;
		  float l3phi     = 0.0;
		  float l3id      = 0.0;
		  float l3iso     = 0.0;
		  float l3pdgId   = 0.0;
		  float l4pt      = 0.0;
		  float l4eta     = 0.0;
		  float l4phi     = 0.0;
		  float l4id      = 0.0;
		  float l4iso     = 0.0;
		  float l4pdgId   = 0.0;
		  float l1pt      = 0.0;
		  float l1eta     = 0.0;
		  float l1phi     = 0.0;
		  float l1pdgId   = 0.0;
		  float l2pt      = 0.0;
		  float l2eta     = 0.0;
		  float l2phi     = 0.0;
		  float l2pdgId   = 0.0;
		  float mela      = 0.0;
		  float numsim    = 0.0;
		  int   event     = 0;
		  int   run       = 0;

		  bchannel   ->SetAddress(&channel); 
		  bz1mass    ->SetAddress(&z1mass);
		  bz2mass    ->SetAddress(&z2mass);
		  bmass      ->SetAddress(&mass);
		  bl3pt      ->SetAddress(&l3pt);
		  bl3eta     ->SetAddress(&l3eta);
		  bl3phi     ->SetAddress(&l3phi);
		  bl3id      ->SetAddress(&l3id);
		  bl3iso     ->SetAddress(&l3iso);
		  bl3pdgId   ->SetAddress(&l3pdgId);
		  bl4pt      ->SetAddress(&l4pt);
		  bl4eta     ->SetAddress(&l4eta);
		  bl4phi     ->SetAddress(&l4phi);
		  bl4id      ->SetAddress(&l4id);
		  bl4iso     ->SetAddress(&l4iso);
		  bl4pdgId   ->SetAddress(&l4pdgId);
		  bl1pt      ->SetAddress(&l1pt);
		  bl1eta     ->SetAddress(&l1eta);
		  bl1phi     ->SetAddress(&l1phi);
		  bl1pdgId   ->SetAddress(&l1pdgId);
		  bl2pt      ->SetAddress(&l2pt);
		  bl2eta     ->SetAddress(&l2eta);
		  bl2phi     ->SetAddress(&l2phi);
		  bl2pdgId   ->SetAddress(&l2pdgId);
		  bmela      ->SetAddress(&mela);
		  bevent     ->SetAddress(&event);
		  brun       ->SetAddress(&run);
		  bnumsim    ->SetAddress(&numsim);

		  for (int i = 0; i < tree->GetEntries(); i++) {
		    bchannel   ->GetEvent(i);
		    bz1mass    ->GetEvent(i);
		    bz2mass    ->GetEvent(i);
		    bmass      ->GetEvent(i);
		    bnumsim    ->GetEvent(i);
		    bl3pt      ->GetEvent(i);
		    bl3phi     ->GetEvent(i);
		    bl3eta     ->GetEvent(i);
		    bl3id      ->GetEvent(i);
		    bl3iso     ->GetEvent(i);
		    bl3pdgId   ->GetEvent(i);
		    bl4pt      ->GetEvent(i);
		    bl4phi     ->GetEvent(i);
		    bl4eta     ->GetEvent(i);
		    bl4id      ->GetEvent(i);
		    bl4iso     ->GetEvent(i);
		    bl4pdgId   ->GetEvent(i);
		    bevent     ->GetEvent(i);
		    brun       ->GetEvent(i);
		    bl1pt      ->GetEvent(i);
		    bl1phi     ->GetEvent(i);
		    bl1eta     ->GetEvent(i);
		    bl1pdgId   ->GetEvent(i);
		    bl2pt      ->GetEvent(i);
		    bl2phi     ->GetEvent(i);
		    bl2eta     ->GetEvent(i);
		    bl2pdgId   ->GetEvent(i);
		    bmela      ->GetEvent(i);

		    bool existsAlready = false;
		    for (std::size_t k = 0; k < runeventinfo.size(); k++) {
		      if (run == runeventinfo[k].first && event == runeventinfo[k].second) existsAlready = true;
		    }

		    if (!existsAlready) {
		      argset.setRealValue("z1mass", z1mass);
		      argset.setRealValue("z2mass", z2mass);
		      argset.setRealValue("mass",   mass);
		      argset.setRealValue("mela",   mela);
		      argset.setRealValue("channel",channel);

		      argset_l1.setRealValue("l1eta", l1eta);
		      argset_l1.setRealValue("l1phi", l1phi);
		      argset_l1.setRealValue("l1pt", l1pt);
		      argset_l1.setRealValue("l1pdgId", l1pdgId);
		      argset_l2.setRealValue("l2eta", l2eta);
		      argset_l2.setRealValue("l2phi", l2phi);
		      argset_l2.setRealValue("l2pt", l2pt);
		      argset_l2.setRealValue("l2pdgId", l2pdgId);
		      argset_l3.setRealValue("l3eta", l3eta);
		      argset_l3.setRealValue("l3phi", l3phi);
		      argset_l3.setRealValue("l3pt", l3pt);
		      argset_l3.setRealValue("l3pdgId", l3pdgId);
		      argset_l4.setRealValue("l4eta", l4eta);
		      argset_l4.setRealValue("l4phi", l4phi);
		      argset_l4.setRealValue("l4pt", l4pt);
		      argset_l4.setRealValue("l4pdgId", l4pdgId);


		      dataset_l1.add(argset_l1);
		      dataset_l2.add(argset_l2);
		      dataset_l3.add(argset_l3);
		      dataset_l4.add(argset_l4);
		      runeventinfo.push_back(std::pair<int, int>(run, event));

		      float weight  = wgt * getPUWeight((int)numsim);
		      float weighterr  = wgterr * getPUWeight((int)numsim);
		      weight *= getSF(l1pt, l1eta, l1pdgId);
		      weight *= getSF(l2pt, l2eta, l2pdgId);
		      weight *= getSF(l3pt, l3eta, l3pdgId);
		      weight *= getSF(l4pt, l4eta, l4pdgId);

		      weighterr *= getSF(l1pt, l1eta, l1pdgId);
		      weighterr *= getSF(l2pt, l2eta, l2pdgId);
		      weighterr *= getSF(l3pt, l3eta, l3pdgId);
		      weighterr *= getSF(l4pt, l4eta, l4pdgId);

		      //if (isSignal) weight    *= 0.5 + 0.5*erf((mass-80.85)/50.42);
		      //if (isSignal) weighterr *= 0.5 + 0.5*erf((mass-80.85)/50.42);

		      argset.setRealValue("weight", weight);
		      argset.setRealValue("weighterr", weighterr);
		      dataset.add(argset);
		    }
		  }
		}

};

#endif
