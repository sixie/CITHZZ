#!/bin/csh


set Ntrial = 50

set barrelOrEndcap = 2

set energyRegrVer = 1
set eleID = 2
set testcat = -1
set fixscale = 1

set fitlow = 75
set fithigh = 105
set binwidth = 1.0

set runLocal = 0

csh sub_testEnergySmearZee.csh $barrelOrEndcap $energyRegrVer $eleID $testcat $fixscale $fitlow $fithigh $binwidth $Ntrial $runLocal

