#include <cmath>
#include <iostream>
using namespace std;

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "PlotHelper2.h"

int main();
void AddPlot(TFile &FData7TeV, TFile &FData8TeV, TFile &FMC7TeV, TFile &FMC8TeV,
   TFile &FSignal7TeV, TFile &FSignal8TeV, TFile &FZX7TeV, TFile &FZX8TeV,
   double Lumi7TeV, double Lumi8TeV, string HistogramName, PsFileHelper &PsFile,
   double LX1 = 0.15, double LY1 = 0.65, double LX2 = 0.35, double LY2 = 0.85);

void AddPlotWithPSSignal(
  TFile &FData7TeV, TFile &FData8TeV, TFile &FMC7TeV, TFile &FMC8TeV,
  TFile &FSignal7TeV, TFile &FSignal8TeV, TFile &FZX7TeV, TFile &FZX8TeV,
  TFile &FPSSignal,
  double Lumi7TeV, double Lumi8TeV, string HistogramName, PsFileHelper &PsFile,
  double LX1 = 0.15, double LY1 = 0.65, double LX2 = 0.35, double LY2 = 0.85);


int main()
{
   double Lumi7TeV = 5.05;
   double Lumi8TeV = 12.1;

   TFile FData7TeV("Plots7TeV_DATA.root", "READ");
   TFile FData8TeV("Plots8TeV_DATA.root", "READ");
   TFile FMC7TeV("All_Plots7TeV.root", "READ");
   TFile FMC8TeV("All_Plots8TeV.root", "READ");
   TFile FSignal7TeV("All_Signal7.root", "READ");
   TFile FSignal8TeV("All_Signal8.root", "READ");
   TFile FZX7TeV("Plots7TeV_ZX.root", "READ");
   TFile FZX8TeV("Plots8TeV_ZX.root", "READ");
   TFile FPSSignal("All_PSSignal8.root", "READ");

   PsFileHelper PsFile("FinalPlots.ps");
   PsFile.AddTextPage("Plots!");

   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HMass1", PsFile, 0.65, 0.65, 0.85, 0.85);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HMass2", PsFile, 0.65, 0.65, 0.85, 0.85);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HMass3", PsFile, 0.65, 0.65, 0.85, 0.85);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HZ1Mass", PsFile);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HZ1MassRestricted1", PsFile);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HZ1MassRestricted2", PsFile);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HZ1MassRestricted3", PsFile);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HZ1MassRestricted4", PsFile);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HZ2Mass", PsFile, 0.40, 0.65, 0.60, 0.85);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HZ2MassRestricted1", PsFile, 0.65, 0.65, 0.85, 0.85);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HZ2MassRestricted2", PsFile, 0.65, 0.65, 0.85, 0.85);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HPT1", PsFile, 0.65, 0.65, 0.85, 0.85);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HPT1Restricted1", PsFile, 0.65, 0.65, 0.85, 0.85);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HPT1Restricted2", PsFile, 0.65, 0.65, 0.85, 0.85);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HPT2", PsFile, 0.65, 0.65, 0.85, 0.85);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HPT2Restricted1", PsFile, 0.65, 0.65, 0.85, 0.85);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HPT2Restricted2", PsFile, 0.65, 0.65, 0.85, 0.85);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HRapidity", PsFile);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HRapidityRestricted1", PsFile);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HRapidityRestricted2", PsFile);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HMELA", PsFile);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HMELARestricted1", PsFile);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HMELARestricted2", PsFile);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HPSMELA", PsFile);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HPSMELARestricted1", PsFile);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HPSMELA_AfterMELACut_Mass110To140", PsFile);
   AddPlot(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV,
      Lumi7TeV, Lumi8TeV, "HPSMELA_AfterMELACut_Mass120To130", PsFile);

   
   AddPlotWithPSSignal( FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV, 
                        FPSSignal, Lumi7TeV, Lumi8TeV, "HZ1MassRestricted1", PsFile);
   AddPlotWithPSSignal(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV, FPSSignal,
                       Lumi7TeV, Lumi8TeV, "HZ1MassRestricted2", PsFile);
   AddPlotWithPSSignal(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV, FPSSignal,
                       Lumi7TeV, Lumi8TeV, "HZ1MassRestricted3", PsFile);
   AddPlotWithPSSignal(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV, FPSSignal,
                       Lumi7TeV, Lumi8TeV, "HZ1MassRestricted4", PsFile);
   AddPlotWithPSSignal(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV, FPSSignal,
                       Lumi7TeV, Lumi8TeV, "HPSMELA", PsFile);
   AddPlotWithPSSignal(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV, FPSSignal,
                       Lumi7TeV, Lumi8TeV, "HPSMELARestricted1", PsFile);
   AddPlotWithPSSignal(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV, FPSSignal,
                       Lumi7TeV, Lumi8TeV, "HPSMELA_AfterMELACut_Mass110To140", PsFile);
   AddPlotWithPSSignal(FData7TeV, FData8TeV, FMC7TeV, FMC8TeV, FSignal7TeV, FSignal8TeV, FZX7TeV, FZX8TeV, FPSSignal,
                       Lumi7TeV, Lumi8TeV, "HPSMELA_AfterMELACut_Mass120To130", PsFile);
   

   PsFile.AddTimeStampPage();
   PsFile.Close();

   FZX8TeV.Close();
   FZX7TeV.Close();
   FSignal8TeV.Close();
   FSignal7TeV.Close();
   FMC8TeV.Close();
   FMC7TeV.Close();
   FData8TeV.Close();
//   FData7TeV.Close();
   FPSSignal.Close();

   return 0;
}


void AddPlot(TFile &FData7TeV, TFile &FData8TeV, TFile &FMC7TeV, TFile &FMC8TeV,
             TFile &FSignal7TeV, TFile &FSignal8TeV, TFile &FZX7TeV, TFile &FZX8TeV,
             double Lumi7TeV, double Lumi8TeV, string HistogramName, PsFileHelper &PsFile,
             double LX1, double LY1, double LX2, double LY2)
{
   TH1D *HData7TeV = (TH1D *)((TH1D *)FData7TeV.Get(HistogramName.c_str()))->Clone("Data7TeV");
   TH1D *HData8TeV = (TH1D *)((TH1D *)FData8TeV.Get(HistogramName.c_str()))->Clone("Data8TeV");
   TH1D *HMC7TeV = (TH1D *)((TH1D *)FMC7TeV.Get(HistogramName.c_str()))->Clone("MC7TeV");
   TH1D *HMC8TeV = (TH1D *)((TH1D *)FMC8TeV.Get(HistogramName.c_str()))->Clone("MC8TeV");
   TH1D *HSignal7TeV = (TH1D *)((TH1D *)FSignal7TeV.Get(HistogramName.c_str()))->Clone("Signal7TeV");
   TH1D *HSignal8TeV = (TH1D *)((TH1D *)FSignal8TeV.Get(HistogramName.c_str()))->Clone("Signal8TeV");
   TH1D *HZX7TeV = (TH1D *)((TH1D *)FZX7TeV.Get(HistogramName.c_str()))->Clone("ZX7TeV");
   TH1D *HZX8TeV = (TH1D *)((TH1D *)FZX8TeV.Get(HistogramName.c_str()))->Clone("ZX8TeV");

   HMC7TeV->Scale(Lumi7TeV);
   HMC8TeV->Scale(Lumi8TeV);
   HSignal7TeV->Scale(Lumi7TeV);
   HSignal8TeV->Scale(Lumi8TeV);

   HData8TeV->Add(HData7TeV);
   HMC8TeV->Add(HMC7TeV);
   HSignal8TeV->Add(HSignal7TeV);
   HZX8TeV->Add(HZX7TeV);

   HMC8TeV->Add(HZX8TeV);
   HSignal8TeV->Add(HMC8TeV);

   TCanvas C;

   HSignal8TeV->SetStats(0);
   HMC8TeV->SetStats(0);
   HZX8TeV->SetStats(0);
   HData8TeV->SetStats(0);

   HSignal8TeV->SetLineColor(kOrange + 10);
   HSignal8TeV->SetLineWidth(2);
   HMC8TeV->SetFillColor(kAzure - 9);
   HZX8TeV->SetFillColor(kGreen - 5);
   HData8TeV->SetLineWidth(2);
   HData8TeV->SetLineColor(kBlack);
   HData8TeV->SetMarkerStyle(11);

   double Maximum = max(HSignal8TeV->GetMaximum(), HData8TeV->GetMaximum());
   Maximum = Maximum + sqrt(Maximum) * 1.5;
   HSignal8TeV->SetMaximum(Maximum);
   
   HSignal8TeV->Draw();
   HMC8TeV->Draw("same");
   HZX8TeV->Draw("same");
   HData8TeV->Draw("same error");
   HSignal8TeV->Draw("axis same");

   TLegend Legend(LX1, LY1, LX2, LY2);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   Legend.SetTextFont(42);
   Legend.AddEntry(HData8TeV, "Data");
   Legend.AddEntry(HSignal8TeV, "Signal");
   Legend.AddEntry(HMC8TeV, "ZZ");
   Legend.AddEntry(HZX8TeV, "ZX");
   Legend.Draw();

   HSignal8TeV->SetMinimum(0);
   C.SetLogy(false);
   PsFile.AddCanvas(C);
   C.SaveAs(("Plots/" + HistogramName + "_Normal.png").c_str());
   C.SaveAs(("Plots/" + HistogramName + "_Normal.C").c_str());
   C.SaveAs(("Plots/" + HistogramName + "_Normal.eps").c_str());
   C.SaveAs(("Plots/" + HistogramName + "_Normal.pdf").c_str());

   HSignal8TeV->SetMinimum(0.31);
   C.SetLogy(true);
   PsFile.AddCanvas(C);
   C.SaveAs(("Plots/" + HistogramName + "_Log.png").c_str());
   C.SaveAs(("Plots/" + HistogramName + "_Log.C").c_str());
   C.SaveAs(("Plots/" + HistogramName + "_Log.eps").c_str());
   C.SaveAs(("Plots/" + HistogramName + "_Log.pdf").c_str());
}


// void AddPlotWithPSSignal( //TFile &FPSSignal, 
//   TFile &FData7TeV, TFile &FData8TeV, TFile &FMC7TeV, TFile &FMC8TeV,
//   TFile &FSignal7TeV, TFile &FSignal8TeV, TFile &FZX7TeV, TFile &FZX8TeV,
//   TFile &FPSSignal,
//   double Lumi7TeV, double Lumi8TeV, string HistogramName, PsFileHelper &PsFile,
//   double LX1, double LY1, double LX2, double LY2)
// {
//   TH1D *HData7TeV = (TH1D *)((TH1D *)FData7TeV.Get(HistogramName.c_str()))->Clone("Data7TeV");
//   TH1D *HData8TeV = (TH1D *)((TH1D *)FData8TeV.Get(HistogramName.c_str()))->Clone("Data8TeV");
//   TH1D *HMC7TeV = (TH1D *)((TH1D *)FMC7TeV.Get(HistogramName.c_str()))->Clone("MC7TeV");
//   TH1D *HMC8TeV = (TH1D *)((TH1D *)FMC8TeV.Get(HistogramName.c_str()))->Clone("MC8TeV");
//   TH1D *HSignal7TeV = (TH1D *)((TH1D *)FSignal7TeV.Get(HistogramName.c_str()))->Clone("Signal7TeV");
//   TH1D *HSignal8TeV = (TH1D *)((TH1D *)FSignal8TeV.Get(HistogramName.c_str()))->Clone("Signal8TeV");
//   TH1D *HZX7TeV = (TH1D *)((TH1D *)FZX7TeV.Get(HistogramName.c_str()))->Clone("ZX7TeV");
//   TH1D *HZX8TeV = (TH1D *)((TH1D *)FZX8TeV.Get(HistogramName.c_str()))->Clone("ZX8TeV");
//   TH1D *HPSSignal = (TH1D *)((TH1D *)FPSSignal.Get(HistogramName.c_str()))->Clone("PSSignal");

//    cout << HistogramName << " : " << FPSSignal.GetName() << endl;
//    cout << HPSSignal->Integral() << endl;
//    cout << HSignal8TeV->Integral() << endl;
//    cout << HData7TeV->Integral() << endl;


//    HMC7TeV->Scale(Lumi7TeV);
//    HMC8TeV->Scale(Lumi8TeV);
//    HSignal7TeV->Scale(Lumi7TeV);
//    HSignal8TeV->Scale(Lumi8TeV);

//    HData8TeV->Add(HData7TeV);
//    HMC8TeV->Add(HMC7TeV);
//    HSignal8TeV->Add(HSignal7TeV);
//    HZX8TeV->Add(HZX7TeV);

//    HMC8TeV->Add(HZX8TeV);
//    HSignal8TeV->Add(HMC8TeV);

//    TCanvas C;

//    HSignal8TeV->SetStats(0);
//    HMC8TeV->SetStats(0);
//    HZX8TeV->SetStats(0);
//    HData8TeV->SetStats(0);

//    HSignal8TeV->SetLineColor(kOrange + 10);
//    HSignal8TeV->SetLineWidth(2);
//    HMC8TeV->SetFillColor(kAzure - 9);
//    HZX8TeV->SetFillColor(kGreen - 5);
//    HData8TeV->SetLineWidth(2);
//    HData8TeV->SetLineColor(kBlack);
//    HData8TeV->SetMarkerStyle(11);

//    double Maximum = max(HSignal8TeV->GetMaximum(), HData8TeV->GetMaximum());
//    Maximum = Maximum + sqrt(Maximum) * 1.5;
//    HSignal8TeV->SetMaximum(Maximum);
   
//    HSignal8TeV->Draw();
//    HMC8TeV->Draw("same");
//    HZX8TeV->Draw("same");
//    HData8TeV->Draw("same error");
//    HSignal8TeV->Draw("axis same");

//    TLegend Legend(LX1, LY1, LX2, LY2);
//    Legend.SetFillStyle(0);
//    Legend.SetBorderSize(0);
//    Legend.SetTextFont(42);
//    Legend.AddEntry(HData8TeV, "Data");
//    Legend.AddEntry(HSignal8TeV, "Signal");
//    Legend.AddEntry(HMC8TeV, "ZZ");
//    Legend.AddEntry(HZX8TeV, "ZX");
//    Legend.Draw();

//    HSignal8TeV->SetMinimum(0);
//    C.SetLogy(false);
//    PsFile.AddCanvas(C);
//    C.SaveAs(("Plots/" + HistogramName + "_Normal.png").c_str());
//    C.SaveAs(("Plots/" + HistogramName + "_Normal.C").c_str());
//    C.SaveAs(("Plots/" + HistogramName + "_Normal.eps").c_str());
//    C.SaveAs(("Plots/" + HistogramName + "_Normal.pdf").c_str());

//    HSignal8TeV->SetMinimum(0.31);
//    C.SetLogy(true);
//    PsFile.AddCanvas(C);
//    C.SaveAs(("Plots/" + HistogramName + "_Log.png").c_str());
//    C.SaveAs(("Plots/" + HistogramName + "_Log.C").c_str());
//    C.SaveAs(("Plots/" + HistogramName + "_Log.eps").c_str());
//    C.SaveAs(("Plots/" + HistogramName + "_Log.pdf").c_str());
// }




void AddPlotWithPSSignal(TFile &FData7TeV, TFile &FData8TeV, TFile &FMC7TeV, TFile &FMC8TeV,
                         TFile &FSignal7TeV, TFile &FSignal8TeV, TFile &FZX7TeV, TFile &FZX8TeV,
                         TFile &FPSSignal,
                         double Lumi7TeV, double Lumi8TeV, string HistogramName, PsFileHelper &PsFile,
                         double LX1, double LY1, double LX2, double LY2)
{
  TH1D *HData7TeV = (TH1D *)((TH1D *)FData7TeV.Get(HistogramName.c_str()))->Clone("Data7TeV");
  TH1D *HData8TeV = (TH1D *)((TH1D *)FData8TeV.Get(HistogramName.c_str()))->Clone("Data8TeV");
  TH1D *HMC7TeV = (TH1D *)((TH1D *)FMC7TeV.Get(HistogramName.c_str()))->Clone("MC7TeV");
  TH1D *HMC8TeV = (TH1D *)((TH1D *)FMC8TeV.Get(HistogramName.c_str()))->Clone("MC8TeV");
  TH1D *HSignal7TeV = (TH1D *)((TH1D *)FSignal7TeV.Get(HistogramName.c_str()))->Clone("Signal7TeV");
  TH1D *HSignal8TeV = (TH1D *)((TH1D *)FSignal8TeV.Get(HistogramName.c_str()))->Clone("Signal8TeV");
  TH1D *HZX7TeV = (TH1D *)((TH1D *)FZX7TeV.Get(HistogramName.c_str()))->Clone("ZX7TeV");
  TH1D *HZX8TeV = (TH1D *)((TH1D *)FZX8TeV.Get(HistogramName.c_str()))->Clone("ZX8TeV");
  TH1D *HPSSignal = (TH1D *)((TH1D *)FPSSignal.Get(HistogramName.c_str()))->Clone("PSSignal");

  HMC7TeV->Scale(Lumi7TeV);
  HMC8TeV->Scale(Lumi8TeV);
  HSignal7TeV->Scale(Lumi7TeV);
  HSignal8TeV->Scale(Lumi8TeV);

  HData8TeV->Add(HData7TeV);
  HMC8TeV->Add(HMC7TeV);
  HSignal8TeV->Add(HSignal7TeV);
  HZX8TeV->Add(HZX7TeV);


  //Scale PS Signal to Sum of 7TeV + 8TeV SM Higgs signal
  HPSSignal->Scale( HSignal8TeV->Integral() / HPSSignal->Integral() );

  //make stacks ourselves
  HMC8TeV->Add(HZX8TeV);
  HSignal8TeV->Add(HMC8TeV);
  HPSSignal->Add(HMC8TeV);

  TCanvas C;

  HSignal8TeV->SetStats(0);
  HMC8TeV->SetStats(0);
  HZX8TeV->SetStats(0);
  HData8TeV->SetStats(0);
  HPSSignal->SetStats(0);

  HSignal8TeV->SetLineColor(kOrange + 10);
  HSignal8TeV->SetLineWidth(2);
  HPSSignal->SetLineColor(kMagenta);
  HPSSignal->SetLineWidth(2);
  HMC8TeV->SetFillColor(kAzure - 9);
  HZX8TeV->SetFillColor(kGreen - 5);
  HData8TeV->SetLineWidth(2);
  HData8TeV->SetLineColor(kBlack);
  HData8TeV->SetMarkerStyle(11);

  double Maximum = max(HSignal8TeV->GetMaximum(), HData8TeV->GetMaximum());
  Maximum = Maximum + sqrt(Maximum) * 1.5;
  HSignal8TeV->SetMaximum(Maximum);
   
  HSignal8TeV->Draw();
  HMC8TeV->Draw("same");
  HZX8TeV->Draw("same");
  HData8TeV->Draw("same error");
  HSignal8TeV->Draw("axis same");
  HPSSignal->Draw("same"); 

  TLegend Legend(LX1, LY1, LX2, LY2);
  Legend.SetFillStyle(0);
  Legend.SetBorderSize(0);
  Legend.SetTextFont(42);
  Legend.AddEntry(HData8TeV, "Data");
  Legend.AddEntry(HSignal8TeV, "SM Higgs");
  Legend.AddEntry(HPSSignal, "PS Higgs");
  Legend.AddEntry(HMC8TeV, "ZZ");
  Legend.AddEntry(HZX8TeV, "ZX");
  Legend.Draw();

  HSignal8TeV->SetMinimum(0);
  C.SetLogy(false);
  PsFile.AddCanvas(C);
  C.SaveAs(("Plots/" + HistogramName + "_WithPSSignal_Normal.png").c_str());
  C.SaveAs(("Plots/" + HistogramName + "_WithPSSignal_Normal.C").c_str());
  C.SaveAs(("Plots/" + HistogramName + "_WithPSSignal_Normal.eps").c_str());
  C.SaveAs(("Plots/" + HistogramName + "_WithPSSignal_Normal.pdf").c_str());

  HSignal8TeV->SetMinimum(0.31);
  C.SetLogy(true);
  PsFile.AddCanvas(C);
  C.SaveAs(("Plots/" + HistogramName + "_WithPSSignal_Log.png").c_str());
  C.SaveAs(("Plots/" + HistogramName + "_WithPSSignal_Log.C").c_str());
  C.SaveAs(("Plots/" + HistogramName + "_WithPSSignal_Log.eps").c_str());
  C.SaveAs(("Plots/" + HistogramName + "_WithPSSignal_Log.pdf").c_str());
}



