// Testing code for checking whether the lepton <=> angle conversion makes sense
// Author: Yi Chen 11223

#include <iostream>
#include <fstream>
using namespace std;

#include "ProgressBar.h"

#include "AngleConversion.h"

int main()
{
   EventParameters TestParameters;
   TestParameters.Phi0 = 0.524;
   TestParameters.Theta0 = 0.653;
   TestParameters.Phi = 2.6623;
   TestParameters.Theta1 = 2.1123;
   TestParameters.Theta2 = 1.95352;
   TestParameters.HMass = 126.5235;
   TestParameters.ZMass = 91.186;
   TestParameters.Z2Mass = 25.3253;
   TestParameters.PhiH = 0.61253;

   cout << "===== Input parameters =====" << endl;
   cout << TestParameters << endl;

   LeptonVectors Leptons = ConvertAnglesToVectors(TestParameters, 30, 0.4);

   cout << "===== Lepton vectors =====" << endl;
   cout << Leptons << endl;

   cout << "Z1 Mass = " << (Leptons.Lepton11 + Leptons.Lepton12).GetMass() << endl;
   cout << "Z2 Mass = " << (Leptons.Lepton21 + Leptons.Lepton22).GetMass() << endl;
   cout << "4L Mass = " << (Leptons.Lepton11 + Leptons.Lepton12 + Leptons.Lepton21 + Leptons.Lepton22).GetMass() << endl;

   EventParameters OutputParameter = ConvertVectorsToAngles(Leptons);

   cout << endl;
   cout << "===== Output parameters =====" << endl;
   cout << OutputParameter << endl;
   cout << endl;

   cout << endl;
   cout << "==> Now the program will loop over all angles (10 bins each), and verify" << endl;
   cout << "   that for ALL cases the output parameter is the same as input parameter" << endl;
   cout << "   In case of mismatch, they will be summarized in the file \"MisMatch.txt\"" << endl;
   cout << "   and also at the end of this program" << endl;
   cout << "   This part should run for 1-2 minutes" << endl;
   cout << endl;

   // now the real deal

   ProgressBar Bar(cout, 100);
   Bar.SetStyle(1);

   ofstream out("MisMatch.txt");

   int FailCount = 0;

   for(double Phi0 = 0; Phi0 < 2 * PI; Phi0 = Phi0 + PI / 5)
   {
      for(double Theta0 = PI / 40; Theta0 < PI; Theta0 = Theta0 + PI / 10)
      {
         Bar.Update(Phi0 * 5 / PI * 10 + Theta0 * 10 / PI);
         Bar.Print();

         for(double Phi = 0; Phi < 2 * PI; Phi = Phi + PI / 5)
         {
            for(double Theta1 = PI / 40; Theta1 < PI; Theta1 = Theta1 + PI / 10)
            {
               for(double Theta2 = PI / 40; Theta2 < PI; Theta2 = Theta2 + PI / 10)
               {
                  for(double Z2Mass = 0; Z2Mass < 30; Z2Mass = Z2Mass + 1)
                  {
                     TestParameters.Phi0 = Phi0;
                     TestParameters.Theta0 = Theta0;
                     TestParameters.Phi = Phi;
                     TestParameters.Theta1 = Theta1;
                     TestParameters.Theta2 = Theta2;
                     TestParameters.Z2Mass = Z2Mass;

                     Leptons = ConvertAnglesToVectors(TestParameters, 40, 1.3);

                     OutputParameter = ConvertVectorsToAngles(Leptons);

                     if(!(OutputParameter == TestParameters))
                     {
                        out << "Mismatch:" << endl;
                        out << endl;
                        out << TestParameters << endl;
                        out << OutputParameter << endl;
                        out << endl;
                        out << endl;

                        FailCount = FailCount + 1;
                     }
                  }
               }
            }
         }
      }
   }

   out.close();

   Bar.Update(100);
   Bar.Print();
   Bar.PrintLine();

   if(FailCount > 0)
      cout << "Failed count = " << FailCount << ".  Contact Yi to debug." << endl;
   else
      cout << "All tests cleared.  Everything matches." << endl;
   cout << endl;

   return 0;
}






