
###########################################################
# Do Fits For Regression Energy
###########################################################


#53X
foreach etype( 0 1 2 3 )
   foreach bin ( 0 1 2 )
     root -l -b -q  CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/FitZMassScaleAndResolution.C+\(\"/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_Summer12DY53X.root\",\"Summer12_DY_53X_EnergyType${etype}_CategoryBin${bin}\",$etype,$bin,false,\"\"\)
     root -l -b -q  CITHZZ/scripts/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/FitZMassScaleAndResolution.C+\(\"/afs/cern.ch/work/s/sixie/public/LeptonScaleAndResolution/ZeeEvents/ZeeNtuple_HCP2012.root\",\"ZMassScaleAndResolutionFits_HCP2012_DY.root\",$etype,$bin,false,\"\"\)
  end
end


###########################################################
# Make Plots
###########################################################
foreach etype( 0 1 2 3 )
   foreach bin ( 0 1 2 )

    root -l -b -q CITHZZ/LeptonScaleAndResolution/Electrons/PlotZMassScaleAndResolutionFit.C+\(\"Notes/ElectronEnergyRegression/ZeeFits/FromAlex/ZeeFitsNew_2012DataRestricted/ZeeFits_2012DataRestricted_EnergyType${etype}_CategoryBin${bin}.root\",\"/afs/cern.ch/work/a/atakeda/public/LeptonScaleAndResolution/ZeeEventsJoined/ZeeNtuple.2012Data.root\",\"ZeeFits_2012Data_EnergyType${etype}_CategoryBin${bin}\",${etype},${bin}\)

    root -l -b -q CITHZZ/LeptonScaleAndResolution/Electrons/PlotZMassScaleAndResolutionFit.C+\(\"Notes/ElectronEnergyRegression/ZeeFits/FromAlex/ZeeFitsNew_MC/ZeeFits_MC_EnergyType${etype}_CategoryBin${bin}.root\",\"/afs/cern.ch/work/a/atakeda/public/LeptonScaleAndResolution/ZeeEventsJoined/ZeeNtuple.s12-zllm50-2-v9.root\",\"ZeeFits_MC_EnergyType${etype}_CategoryBin${bin}\",${etype},${bin}\)

  end
end

