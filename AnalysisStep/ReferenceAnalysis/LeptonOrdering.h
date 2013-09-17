#ifndef LEPTONORDERING_H
#define LEPTONORDERING_H

#include "TVector3.h"
#include "TLorentzVector.h"
#include <cmath>

// Function to calculate "proximity" of two leptons based on various criteria
Float_t proximity (TLorentzVector l1, TLorentzVector l2, Int_t option = 1) {
  /*
     Options
     1  use delta R to measure proximity
     2  use dot product of unit 3-momentum vectors
     3  use closeness of reconstructed invariant mass to on-shell mass of Z boson
     4  use sum of absolute values of transverse momentum
   */

  if (option == 1)
    return -l1.DeltaR(l2);
  if (option == 2) {
    TVector3 l1_3p(l1.Px(), l1.Py(), l1.Pz());
    TVector3 l2_3p(l2.Px(), l2.Py(), l2.Pz());
    return l1_3p.Dot(l2_3p);
  }
  if (option == 3) {
    Float_t Z_onshell_mass = 91.2;
    return -fabs((l1 + l2).M() - Z_onshell_mass);
  }
  if (option == 4) {
    return fabs(l1.Pt()) + fabs(l2.Pt());
  }
}

// Function to order four leptons based on different proximity criteria
void order_lepton_list(std::vector<TLorentzVector>& leptons, std::vector<Int_t>& leptonIDs, Int_t prox_option = 0) {

  /*
     Orders a set of four leptons (leptons and leptonIDs) based on the desired proximity option (prox_option)
     Ordered set has the following structure: 
     (l1, l2) and (l3, l4) form charge- and lepton number- neutral pairs (l1id + l2id == l3id + l4id == 0)

     Options
     0		keep ordering of leptons
     1		(l1, l2) has smaller delta R than (l3, l4)
     2		(l1, l2) has larger 3-momenta dot product than (l3, l4)
     3		(l1, l2) has reconstructed mass closer to the Z pole than (l3, l4)
     4		(l1, l2) has larger sum of absolute values of transverse momentum than (l3, l4)
     5		pairing is such that Z1 and Z2 are most back-to-back, i.e. have the biggest delta phi
     6		pairing is such that the mass difference between the two reconstructed bosons is reduced
     7		pairing is such that, in the center of mass frame, the momenta of the Z's is minimized

   */


  TLorentzVector l1(leptons[0]); Int_t l1id = leptonIDs[0];
  TLorentzVector l2(leptons[1]); Int_t l2id = leptonIDs[1];
  TLorentzVector l3(leptons[2]); Int_t l3id = leptonIDs[2];
  TLorentzVector l4(leptons[3]); Int_t l4id = leptonIDs[3];
  Float_t prox = 0.; Float_t biggest_prox = -1000000.;

  // If option is zero then use the ordering given in tree
  if (prox_option == 0) return;

  // If option is five, use different pairing scheme, not based on two-lepton proximity functions
  if (prox_option == 5) {
    if (l1id + l2id == 0 && l3id + l4id == 0) {
      prox = fabs((l1 + l2).DeltaPhi((l3 + l4)));
      if (prox > biggest_prox){
	biggest_prox = prox;
	leptons[0] = l1; leptonIDs[0] = l1id;
	leptons[1] = l2; leptonIDs[1] = l2id;
	leptons[2] = l3; leptonIDs[2] = l3id;
	leptons[3] = l4; leptonIDs[3] = l4id;
      }
    }

    if (l1id + l3id == 0 && l2id + l4id == 0) {
      prox = fabs((l1 + l3).DeltaPhi((l2 + l4)));
      if (prox > biggest_prox){
	biggest_prox = prox;
	leptons[0] = l1; leptonIDs[0] = l1id;
	leptons[1] = l3; leptonIDs[1] = l3id;
	leptons[2] = l2; leptonIDs[2] = l2id;
	leptons[3] = l4; leptonIDs[3] = l4id;
      }
    }

    if (l1id + l4id == 0 && l2id + l3id == 0) {
      prox = fabs((l1 + l4).DeltaPhi((l2 + l3)));
      if (prox > biggest_prox){
	biggest_prox = prox;
	leptons[0] = l1; leptonIDs[0] = l1id;
	leptons[1] = l4; leptonIDs[1] = l4id;
	leptons[2] = l2; leptonIDs[2] = l2id;
	leptons[3] = l3; leptonIDs[3] = l3id;
      }
    }
  }

  // If option is six, also do not use two-electron procimity function
  if (prox_option == 6) {
    if (l1id + l2id == 0 && l3id + l4id == 0) {
      prox = - fabs((l1 + l2).M() - (l3 + l4).M());
      if (prox > biggest_prox){
	biggest_prox = prox;
	leptons[0] = l1; leptonIDs[0] = l1id;
	leptons[1] = l2; leptonIDs[1] = l2id;
	leptons[2] = l3; leptonIDs[2] = l3id;
	leptons[3] = l4; leptonIDs[3] = l4id;
      }
    }

    if (l1id + l3id == 0 && l2id + l4id == 0) {
      prox = - fabs((l1 + l3).M() - (l2 + l4).M());
      if (prox > biggest_prox){
	biggest_prox = prox;
	leptons[0] = l1; leptonIDs[0] = l1id;
	leptons[1] = l3; leptonIDs[1] = l3id;
	leptons[2] = l2; leptonIDs[2] = l2id;
	leptons[3] = l4; leptonIDs[3] = l4id;
      }
    }

    if (l1id + l4id == 0 && l2id + l3id == 0) {
      prox = - fabs((l1 + l4).M() - (l2 + l3).M());
      if (prox > biggest_prox){
	biggest_prox = prox;
	leptons[0] = l1; leptonIDs[0] = l1id;
	leptons[1] = l4; leptonIDs[1] = l4id;
	leptons[2] = l2; leptonIDs[2] = l2id;
	leptons[3] = l3; leptonIDs[3] = l3id;
      }
    }
  }

  // If option is seven, boost to the center of mass frame and choose the pairing where the Zs have the smallest momentum
  if (prox_option == 7) {
    // Boost everything to the center of mass frame
    TVector3 totalBoost = - (l1 + l2 + l3 + l4).BoostVector();
    l1.Boost(totalBoost); l2.Boost(totalBoost); l3.Boost(totalBoost); l4.Boost(totalBoost);

    if (l1id + l2id == 0 && l3id + l4id == 0) {
      prox = - ((l1 + l2).Rho() + (l3 + l4).Rho());
      if (prox > biggest_prox){
	biggest_prox = prox;
	leptons[0] = l1; leptonIDs[0] = l1id;
	leptons[1] = l2; leptonIDs[1] = l2id;
	leptons[2] = l3; leptonIDs[2] = l3id;
	leptons[3] = l4; leptonIDs[3] = l4id;
      }
    }

    if (l1id + l3id == 0 && l2id + l4id == 0) {
      prox = - ((l1 + l3).Rho() + (l2 + l4).Rho());
      if (prox > biggest_prox){
	biggest_prox = prox;
	leptons[0] = l1; leptonIDs[0] = l1id;
	leptons[1] = l3; leptonIDs[1] = l3id;
	leptons[2] = l2; leptonIDs[2] = l2id;
	leptons[3] = l4; leptonIDs[3] = l4id;
      }
    }

    if (l1id + l4id == 0 && l2id + l3id == 0) {
      prox = - ((l1 + l4).Rho() + (l2 + l3).Rho());
      if (prox > biggest_prox){
	biggest_prox = prox;
	leptons[0] = l1; leptonIDs[0] = l1id;
	leptons[1] = l4; leptonIDs[1] = l4id;
	leptons[2] = l2; leptonIDs[2] = l2id;
	leptons[3] = l3; leptonIDs[3] = l3id;
      }
    }
  }

  // For options that use two-lepton proximity criteria
  if (prox_option > 0 && prox_option < 5) {	
    if (l1id + l2id == 0 && l3id + l4id == 0) {
      prox = proximity(l1, l2, prox_option);
      if (prox > biggest_prox) {
	biggest_prox = prox;
	leptons[0] = l1; leptonIDs[0] = l1id;
	leptons[1] = l2; leptonIDs[1] = l2id;
	leptons[2] = l3; leptonIDs[2] = l3id;
	leptons[3] = l4; leptonIDs[3] = l4id;
      }

      prox = proximity(l3, l4, prox_option);
      if (prox > biggest_prox) {
	biggest_prox = prox;
	leptons[0] = l3; leptonIDs[0] = l3id;
	leptons[1] = l4; leptonIDs[1] = l4id;
	leptons[2] = l1; leptonIDs[2] = l1id;
	leptons[3] = l2; leptonIDs[3] = l2id;
      }
    }

    if (l1id + l3id == 0 && l2id + l4id == 0) {
      prox = proximity(l1, l3, prox_option);
      if (prox > biggest_prox) {
	biggest_prox = prox;
	leptons[0] = l1; leptonIDs[0] = l1id;
	leptons[1] = l3; leptonIDs[1] = l3id;
	leptons[2] = l2; leptonIDs[2] = l2id;
	leptons[3] = l4; leptonIDs[3] = l4id;
      }

      prox = proximity(l2, l4, prox_option);
      if (prox > biggest_prox) {
	biggest_prox = prox;
	leptons[0] = l2; leptonIDs[0] = l2id;
	leptons[1] = l4; leptonIDs[1] = l4id;
	leptons[2] = l1; leptonIDs[2] = l1id;
	leptons[3] = l3; leptonIDs[3] = l3id;
      }
    }

    if (l1id + l4id == 0 && l2id + l3id == 0) {
      prox = proximity(l1, l4, prox_option);
      if (prox > biggest_prox) {
	biggest_prox = prox;
	leptons[0] = l1; leptonIDs[0] = l1id;
	leptons[1] = l4; leptonIDs[1] = l4id;
	leptons[2] = l2; leptonIDs[2] = l2id;
	leptons[3] = l3; leptonIDs[3] = l3id;
      }

      prox = proximity(l2, l3, prox_option);
      if (prox > biggest_prox) {
	biggest_prox = prox;
	leptons[0] = l2; leptonIDs[0] = l2id;
	leptons[1] = l3; leptonIDs[1] = l3id;
	leptons[2] = l1; leptonIDs[2] = l1id;
	leptons[3] = l4; leptonIDs[3] = l4id;
      }
    }
  }
}

#endif
