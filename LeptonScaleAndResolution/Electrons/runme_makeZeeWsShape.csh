#!/bin/csh

which root

set dataOrMC = AAA
set energyRegrVer = BBB
set eleID = CCC
set applyEScale = DDD

set workdir = /afs/cern.ch/work/y/yangyong/CMSSW/CMSSW_5_2_4/src/CITHZZ/LeptonScaleAndResolution/Electrons

root -b <<EOF
gSystem->Load("$workdir/makeZeeWsShape_C.so")
makeZeeWsShape($dataOrMC,$energyRegrVer,$eleID,$applyEScale)
.qqq
EOF

cp makeZeeWsShape*root $workdir/zeeWS
