##################################################################
## Applies the regression weights for the four versions		##
## replace the filenames with the ones you got from		##
## trainElectronEnergyRegression.C				##
##################################################################

root -l -b -q CITHZZ/LeptonScaleAndResolution/Electrons/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/applyElectronEnergyRegression.C+'("/afs/cern.ch/work/a/atakeda/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12HZZ4l.Testing.root", "weights_Summer12HZZ4l/weightFile_V00.root", "targets_Summer12HZZ4l/", "V00")'

root -l -b -q CITHZZ/LeptonScaleAndResolution/Electrons/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/applyElectronEnergyRegression.C+'("/afs/cern.ch/work/a/atakeda/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12HZZ4l.Testing.root", "weights_Summer12HZZ4l/weightFile_V01.root", "targets_Summer12HZZ4l/", "V01")'

root -l -b -q CITHZZ/LeptonScaleAndResolution/Electrons/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/applyElectronEnergyRegression.C+'("/afs/cern.ch/work/a/atakeda/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12HZZ4l.Testing.root", "weights_Summer12HZZ4l/weightFile_V10.root", "targets_Summer12HZZ4l/", "V10")'

root -l -b -q CITHZZ/LeptonScaleAndResolution/Electrons/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/applyElectronEnergyRegression.C+'("/afs/cern.ch/work/a/atakeda/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12HZZ4l.Testing.root", "weights_Summer12HZZ4l/weightFile_V11.root", "targets_Summer12HZZ4l/", "V11")'

