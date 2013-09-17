#! /bin/bash

WORK=$PWD
export SCRAM_ARCH=slc5_amd64_gcc462
SRC=/afs/cern.ch/work/e/emanuele/hzz4l/CMSSW_5_3_3_patch3/src
if [[ "${LSB_JOBID}" != "" ]]; then
    export TMPDIR=$WORK;
fi;

cd $SRC;
eval $(scramv1 runtime -sh);
cd $WORK

mkdir -p $SRC/out/${LSB_JOBID}
echo combine $* 2>&1
combine $* 
mv *.root $SRC/out/${LSB_JOBID}/ -v;
