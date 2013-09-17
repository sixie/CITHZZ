#!/bin/csh

if ( ${#argv} != 10 ) then
echo "Usage: $0 barrelOrEndcap energyRegrVer eleID testcat fixscale fitlow fithigh binwidth Ntrial runlocally";
exit
endif

set barrelOrEndcap = $1
set energyRegrVer = $2
set eleID = $3
set testcat = $4
set fixscale = $5
set fitlow = $6
set fithigh = $7
set binwidth = $8
set Ntrial = $9
set runLocal = $10

set jobfile = jobs/runme_testEnergySmearZee.$barrelOrEndcap.$energyRegrVer.$eleID.$testcat.$fixscale.$fitlow.$fithigh.$binwidth.$Ntrial.csh

cp runme_testEnergySmearZee.csh $jobfile

perl ReplaceString.pl AAA $barrelOrEndcap $jobfile
perl ReplaceString.pl BBB $energyRegrVer $jobfile
perl ReplaceString.pl CCC $eleID $jobfile
perl ReplaceString.pl DDD $testcat $jobfile
perl ReplaceString.pl EEE $fixscale $jobfile
perl ReplaceString.pl FFF $fitlow $jobfile
perl ReplaceString.pl GGG $fithigh $jobfile
perl ReplaceString.pl HHH $binwidth $jobfile
perl ReplaceString.pl III $Ntrial $jobfile


if( $runLocal == 1) then

else
bsub -q cmscaf1nd -J test <  $jobfile
endif
