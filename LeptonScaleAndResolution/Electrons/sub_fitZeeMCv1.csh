#!/bin/csh

if ( ${#argv} != 5 ) then
echo "Usage: $0 energyRegrVer eleID fitmin fitmax runlocally";
exit
endif

set energyRegrVer = $1
set eleID = $2
set fitmin = $3
set fitmax = $4
set runLocal = $5

set jobfile = jobs/runme_fitZeeMCv1.$energyRegrVer.$eleID.$fitmin.$fitmax.csh

cp runme_fitZeeMCv1.csh $jobfile

perl ReplaceString.pl AAA $energyRegrVer $jobfile
perl ReplaceString.pl BBB $eleID $jobfile
perl ReplaceString.pl CCC $fitmin $jobfile
perl ReplaceString.pl DDD $fitmax $jobfile

if( $runLocal == 1) then

else
bsub -q cmscaf1nh -J test <  $jobfile
endif

