//---------------------------------------------------------------------------
// AttachLikelihood.cpp
// Author: Yi Chen (11231)
//
// This program reads in a file, get tree named "probe_tree",
// evaluate likelihood using a 6D map and make a new branch with
// the log likelihood for each entry
//
//
// To run it, compile as usual
//
// > g++ AttachLikelihood.cpp -o AttachBranch -I/Path/To/Header/Files \
//    `root-config --cflags` `root-config --glibs`
//
//
// And then execute the program
//
// > time ./AttachBranch RootFile.root MapFile.map logLikelihood
//
// where the arguments are the root file containing probe_tree,
// the 6D map we want to use, and last one is the new name of the
// branch we want to create.
//---------------------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
using namespace std;
//---------------------------------------------------------------------------
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
//---------------------------------------------------------------------------
#include "TauHelperFunctions2.h"
#include "ProgressBar.h"
//---------------------------------------------------------------------------
#include "AngleConversion.h"
//---------------------------------------------------------------------------
class Index;
int main(int argc, char *argv[]);
//---------------------------------------------------------------------------
class Index
{
public:
   int Indices[6];
public:
   Index();
   ~Index() {}
   bool operator <(const Index &other) const;
   bool operator ==(const Index &other) const;
   void Print();
};
//---------------------------------------------------------------------------
Index::Index()
{
   for(int i = 0; i < 6; i++)
      Indices[i] = 0;
}
//---------------------------------------------------------------------------
bool Index::operator <(const Index &other) const
{
   for(int i = 0; i < 6; i++)
   {
      if(Indices[i] < other.Indices[i])   return true;
      if(Indices[i] > other.Indices[i])   return false;
   }

   return false;
}
//---------------------------------------------------------------------------
bool Index::operator ==(const Index &other) const
{
   for(int i = 0; i < 6; i++)
      if(Indices[i] != other.Indices[i])
         return false;
   return true;
}
//---------------------------------------------------------------------------
void Index::Print()
{
   cout << "(";
   for(int i = 0; i < 5; i++)
      cout << Indices[i] << " ";
   cout << Indices[5] << ")";
}
//---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
   if(argc < 4)
   {
      cerr << "Usage: " << argv[0] << " FileWithTree MapFile BranchName" << endl;
      return -1;
   }

   // Input parameters
   string TreeFile = "HZZEventNtuple_HZZ125.root";
   string MapFile = "All_PseudoScalar.map";
   string BranchName = "logLikelihood_pseudoscalar";

   TreeFile = argv[1];
   MapFile = argv[2];
   BranchName = argv[3];

   // Open file in update mode, and get tree
   TFile File(TreeFile.c_str(), "UPDATE");
   TTree *Tree = (TTree *)File.Get("probe_tree");

   // Set branch addresses
   float l1pt, l1eta, l1phi;
   float l2pt, l2eta, l2phi;
   float l3pt, l3eta, l3phi;
   float l4pt, l4eta, l4phi;
   bool passFullSelection;
   Tree->SetBranchAddress("l1pt", &l1pt);
   Tree->SetBranchAddress("l1eta", &l1eta);
   Tree->SetBranchAddress("l1phi", &l1phi);
   Tree->SetBranchAddress("l2pt", &l2pt);
   Tree->SetBranchAddress("l2eta", &l2eta);
   Tree->SetBranchAddress("l2phi", &l2phi);
   Tree->SetBranchAddress("l3pt", &l3pt);
   Tree->SetBranchAddress("l3eta", &l3eta);
   Tree->SetBranchAddress("l3phi", &l3phi);
   Tree->SetBranchAddress("l4pt", &l4pt);
   Tree->SetBranchAddress("l4eta", &l4eta);
   Tree->SetBranchAddress("l4phi", &l4phi);
   Tree->SetBranchAddress("passFullSelection", &passFullSelection);

   // Get ranges and bin counts from map
   double HMass = -1, ZMass = -1;
   int NBinsMass = -1, NBinsPhi0 = -1, NBinsTheta0 = -1, NBinsPhi = -1, NBinsTheta1 = -1, NBinsTheta2 = -1;
   double Phi0Min = 0, Phi0Max = 0;
   double Theta0Min = 0, Theta0Max = 0;
   double PhiMin = 0, PhiMax = 0;
   double Theta1Min = 0, Theta1Max = 0;
   double Theta2Min = 0, Theta2Max = 0;
   double MassMin = 0, MassMax = 0;
  
   ifstream in(MapFile.c_str());
   in >> HMass >> ZMass;
   in >> NBinsPhi0 >> NBinsTheta0 >> NBinsPhi >> NBinsTheta1 >> NBinsTheta2 >> NBinsMass;
   in >> Phi0Min >> Phi0Max;
   in >> Theta0Min >> Theta0Max;
   in >> PhiMin >> PhiMax;
   in >> Theta1Min >> Theta1Max;
   in >> Theta2Min >> Theta2Max;
   in >> MassMin >> MassMax;

   // First loop over the entries and create index
   map<int, Index> EventIndices;
   vector<Index> AvailableIndices;
   map<Index, double> Likelihoods;

   int EntryCount = Tree->GetEntries();
   
   cout << "Loop over tree and collect indices" << endl;
   ProgressBar Bar1(cout, EntryCount);
   Bar1.SetStyle(1);

   for(int iEvent = 0; iEvent < EntryCount; iEvent++)
   {
      Tree->GetEntry(iEvent);

      if(iEvent % 1000 == 0)
      {
         Bar1.Update(iEvent + 1);
         Bar1.Print();
      }

      if(passFullSelection == false)
         continue;

      LeptonVectors Leptons;
      Leptons.Lepton11.SetPtEtaPhi(l1pt, l1eta, l1phi);
      Leptons.Lepton12.SetPtEtaPhi(l2pt, l2eta, l2phi);
      Leptons.Lepton21.SetPtEtaPhi(l3pt, l3eta, l3phi);
      Leptons.Lepton22.SetPtEtaPhi(l4pt, l4eta, l4phi);

      EventParameters Angles = ConvertVectorsToAngles(Leptons);

      Index index;
      index.Indices[0] = (Angles.Phi0 - Phi0Min) / (Phi0Max - Phi0Min) * NBinsPhi0;
      index.Indices[1] = (cos(Angles.Theta0) - Theta0Min) / (Theta0Max - Theta0Min) * NBinsTheta0;
      index.Indices[2] = (Angles.Phi - PhiMin) / (PhiMax - PhiMin) * NBinsPhi;
      index.Indices[3] = (cos(Angles.Theta1) - Theta1Min) / (Theta1Max - Theta1Min) * NBinsTheta1;
      index.Indices[4] = (cos(Angles.Theta2) - Theta2Min) / (Theta2Max - Theta2Min) * NBinsTheta2;
      index.Indices[5] = (Angles.Z2Mass - MassMin) / (MassMax - MassMin) * NBinsMass;

      if(index.Indices[0] >= NBinsPhi0)
         index.Indices[0] = NBinsPhi0 - 1;
      if(index.Indices[1] >= NBinsTheta0)
         index.Indices[1] = NBinsTheta0 - 1;
      if(index.Indices[2] >= NBinsPhi)
         index.Indices[2] = NBinsPhi - 1;
      if(index.Indices[3] >= NBinsTheta1)
         index.Indices[3] = NBinsTheta1 - 1;
      if(index.Indices[4] >= NBinsTheta2)
         index.Indices[4] = NBinsTheta2 - 1;
      if(index.Indices[5] >= NBinsMass)
         index.Indices[5] = NBinsMass - 1;

      EventIndices.insert(pair<int, Index>(iEvent, index));

      AvailableIndices.push_back(index);
      if(iEvent % 10000 == 0)
      {
         sort(AvailableIndices.begin(), AvailableIndices.end());
         AvailableIndices.erase(unique(AvailableIndices.begin(), AvailableIndices.end()), AvailableIndices.end());
      }
   }
         
   Bar1.Update(EntryCount);
   Bar1.Print();
   Bar1.PrintLine();
   cout << endl;

   sort(AvailableIndices.begin(), AvailableIndices.end());
   AvailableIndices.erase(unique(AvailableIndices.begin(), AvailableIndices.end()), AvailableIndices.end());

   // Evaluate likelihoods - calculate integral along the way
   bool ErrorMessagePrinted = false;
   double Integral = 0;
   int VectorCount = 0;

   cout << "Starting evaluating likelihoods..." << endl;
   ProgressBar Bar2(cout, NBinsPhi0 * NBinsTheta0);
   Bar2.SetStyle(1);

   for(int iPhi0 = 0; iPhi0 < NBinsPhi0; iPhi0++)
   {
      for(int iTheta0 = 0; iTheta0 < NBinsTheta0; iTheta0++)
      {
         Bar2.Update(iPhi0 * NBinsPhi0 + iTheta0 + 1);
         Bar2.Print();

         for(int iPhi = 0; iPhi < NBinsPhi; iPhi++)
         {
            for(int iTheta1 = 0; iTheta1 < NBinsTheta1; iTheta1++)
            {
               for(int iTheta2 = 0; iTheta2 < NBinsTheta2; iTheta2++)
               {
                  for(int iZ2Mass = 0; iZ2Mass < NBinsMass; iZ2Mass++)
                  {
                     double Temp = -1;
                     in >> Temp;
                     
                     Integral = Integral + Temp;

                     Index ThisEntry;
                     ThisEntry.Indices[0] = iPhi0;
                     ThisEntry.Indices[1] = iTheta0;
                     ThisEntry.Indices[2] = iPhi;
                     ThisEntry.Indices[3] = iTheta1;
                     ThisEntry.Indices[4] = iTheta2;
                     ThisEntry.Indices[5] = iZ2Mass;

                     if(VectorCount >= AvailableIndices.size())
                        continue;

                     if(ThisEntry == AvailableIndices[VectorCount])
                     {
                        Likelihoods.insert(pair<Index, double>(ThisEntry, Temp));
                        VectorCount = VectorCount + 1;
                     }
                     else if(AvailableIndices[VectorCount] < ThisEntry && ErrorMessagePrinted == false)
                     {
                        cerr << "This should not happen!  Contact Yi for more detail and to debug!" << endl;
                        ThisEntry.Print();
                        cout << " ";
                        AvailableIndices[VectorCount].Print();
                        cout << endl;
                        ErrorMessagePrinted = true;
                     }
                  }
               }
            }
         }
      }
   }
   AvailableIndices.erase(AvailableIndices.begin(), AvailableIndices.end());

   Bar2.Update(NBinsPhi0 * NBinsTheta0);
   Bar2.Print();
   Bar2.PrintLine();
   cout << endl;

   // Add a new branch to tree
   float logLikelihood = 0;
   TBranch *bLikelihood = Tree->Branch(BranchName.c_str(), &logLikelihood, (BranchName + "/F").c_str());

   // Loop over entries and include the new branch
   cout << "Attaching a new branch" << endl;

   for(int iEntry = 0; iEntry < EntryCount; iEntry++)
   {
      if(iEntry % 1000 == 0)
      {
         Bar1.Update(iEntry + 1);
         Bar1.Print();
      }
      
      if(EventIndices.find(iEntry) == EventIndices.end())
         logLikelihood = 0;
      else
         logLikelihood = log(Likelihoods[EventIndices[iEntry]]) - log(Integral);
         
      bLikelihood->Fill();
   }
   
   Bar1.Update(EntryCount);
   Bar1.Print();
   Bar1.PrintLine();
   cout << endl;

   // Overwrite the tree
   Tree->Write("", TObject::kOverwrite);

   // Clean-up
   in.close();
   File.Close();

   return 0;
}
//---------------------------------------------------------------------------





