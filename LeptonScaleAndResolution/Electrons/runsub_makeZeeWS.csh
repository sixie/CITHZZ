#!/bin/csh

set dataOrMC = 2
set eleID = 2

set applyEScale = 0
set runLocally = 0

set energyRegrVer = 0
while ($energyRegrVer <= 2)
csh sub_makeZeeWsShape.csh $dataOrMC $energyRegrVer $eleID $applyEScale $runLocally
@ energyRegrVer ++
end
