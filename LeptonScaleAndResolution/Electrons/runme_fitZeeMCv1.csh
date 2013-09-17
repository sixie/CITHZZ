#!/bin/csh

which root

set energyRegrVer = AAA
set eleID = BBB
set fitmin = CCC
set fitmax = DDD

set workdir = /afs/cern.ch/work/y/yangyong/CMSSW/CMSSW_5_2_4/src/CITHZZ/LeptonScaleAndResolution/Electrons

root -b <<EOF
gSystem->Load("$workdir/fitZeeMCv1_C.so")
fitZeeMCv1($energyRegrVer,$eleID,$fitmin,$fitmax)
.qqq
EOF

cp *fitres $workdir/fitres/
cp *pdf $workdir/fitres/
