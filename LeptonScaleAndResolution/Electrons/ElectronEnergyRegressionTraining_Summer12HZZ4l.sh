##################################################
## Doing the training for different options     ##
## replace the directory names here with        ##
## you own exisiting directories		##
##################################################
root -l -b -q CITHZZ/LeptonScaleAndResolution/Electrons/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/trainElectronEnergyRegression.C+'("/afs/cern.ch/work/a/atakeda/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12HZZ4l.Training.root", "weights_Summer12HZZ4l/weightFile_V00.root","V00")'

root -l -b -q CITHZZ/LeptonScaleAndResolution/Electrons/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/trainElectronEnergyRegression.C+'("/afs/cern.ch/work/a/atakeda/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12HZZ4l.Training.root", "weights_Summer12HZZ4l/weightFile_V01.root","V01")'

root -l -b -q CITHZZ/LeptonScaleAndResolution/Electrons/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/trainElectronEnergyRegression.C+'("/afs/cern.ch/work/a/atakeda/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12HZZ4l.Training.root", "weights_Summer12HZZ4l/weightFile_V10.root","V10")'

root -l -b -q CITHZZ/LeptonScaleAndResolution/Electrons/rootlogon.C CITHZZ/LeptonScaleAndResolution/Electrons/trainElectronEnergyRegression.C+'("/afs/cern.ch/work/a/atakeda/public/LeptonScaleAndResolution/Electrons/ElectronNtuples.Summer12HZZ4l.Training.root", "weights_Summer12HZZ4l/weightFile_V11.root","V11")'
