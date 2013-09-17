#!/bin/csh

if ( ${#argv} != 5 ) then
echo "Usage: $0 dataOrMC energyRegrVer eleID applyEScale runLocally";
exit
endif

set dataOrMC = $1
set energyRegrVer = $2
set eleID = $3
set applyEScale = $4
set runLocal = $5


if($dataOrMC == 2 && $applyEScale > 0)  then
echo "MC applyEScale=0 only!"
exit
endif

set output = zeeWS/makeZeeWsShape_dm${dataOrMC}_regver${energyRegrVer}_eleID${eleID}_escale${applyEScale}.root
if(-e $output) then
echo "done already "$output
exit
endif


set jobfile = jobs/runme_makeZeeWsShape.$dataOrMC.$energyRegrVer.$eleID.$applyEScale.csh

cp runme_makeZeeWsShape.csh $jobfile
perl ReplaceString.pl AAA $dataOrMC $jobfile
perl ReplaceString.pl BBB $energyRegrVer $jobfile
perl ReplaceString.pl CCC $eleID $jobfile
perl ReplaceString.pl DDD $applyEScale $jobfile


if( $runLocal == 1) then

else
bsub -q cmscaf1nh -J test <  $jobfile
endif

