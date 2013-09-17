#########################
## Script to make TeX table with fit values
#########################

root -l fitWork/fitResultsToTex.C+'("/afs/cern.ch/work/a/atakeda/public/release/CMSSW_5_2_4_patch1/src/ZeeFits_2012Data/ZeeFits_2012Data", "/afs/cern.ch/work/a/atakeda/public/release/CMSSW_5_2_4_patch1/src/ZeeFits_MC/ZeeFits_MC", "/afs/cern.ch/work/a/atakeda/public/release/CMSSW_5_2_4_patch1/src/testTable.tex")'


root -l CITHZZ/LeptonScaleAndResolution/Electrons/fitResultsToTex.C+'("Notes/ElectronEnergyRegression/ZeeFits/FromAlex/ZeeFitsNew_2012DataRestricted/ZeeFits_2012DataRestricted", "Notes/ElectronEnergyRegression/ZeeFits/FromAlex/ZeeFitsNew_MC/ZeeFits_MC", "ElectronEnergyScaleAndResolutionTable.tex")'



