#!/usr/bin/python

# execute with: python runfits.py
import sys
import os

curdir=os.getcwd()

from ROOT import gROOT
gROOT.Macro(curdir+'/RooLogon.C')
gROOT.LoadMacro('src/FitZElectronsData.cc')

from ROOT import FitZElectrons,PlotZElectrons,calcSWeight

for pt in range(0,1):
    print 'fitting pt bin ', pt
    for eta in range(0,2):
        print '     and eta bin ', eta
        FitZElectrons(pt,eta)
        PlotZElectrons(pt,eta)
        calcSWeight(pt,eta)
print 'done.'

        


