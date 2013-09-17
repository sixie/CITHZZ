#!/bin/csh

which root

set barrelOrEndcap = AAA
set energyRegrVer = BBB
set eleID = CCC
set testcat = DDD
set fixscale = EEE
set fitlow = FFF
set fithigh = GGG
set binwidth = HHH
set Ntrial = III


set workdir = /afs/cern.ch/work/y/yangyong/CMSSW/CMSSW_5_2_4/src/CITHZZ/LeptonScaleAndResolution/Electrons


root -b <<EOF
gSystem->Load("$workdir/testEnergySmearZee_C.so")
testEnergySmearZee($barrelOrEndcap,$energyRegrVer,$eleID,$testcat,$fixscale,$fitlow,$fithigh,$binwidth,$Ntrial)
.qqq
EOF


cp testEnergySmearZee*root $workdir/ressmear
