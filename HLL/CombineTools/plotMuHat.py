#! /usr/bin/env python

from ROOT import gROOT, TLegend, TCanvas, TFile, TGraph, TMultiGraph
import sys, os, math, argparse
import rootlogonTDR
from array import array

from drawFunction import *

parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),description="analysis")
parser.add_argument('-t',dest='tag',action='store',required=True,help='tag')
parser.add_argument('-d',dest='dir',action='store',required=True,help='dir')
args=parser.parse_args()
if args.tag not in ['8TeV_4mu','8TeV_4e','8TeV_2e2mu','7TeV','8TeV','7p8TeV']: args.tag='8TeV'

masses = range(110,160,1)
masses1 = range(162,290,2)
masses2 = range(295,350,5)
masses3 = range(360,400,10)
masses4 = range(420,1000,20)

masses.extend(masses1)
masses.extend(masses2)
masses.extend(masses3)
masses.extend(masses4)

print "plotting mu for masses "+str(masses)

limits = {}
for mass in masses:
    limits[mass] = []
    print             args.dir+'/'+args.tag+'/higgsCombineTest.MaxLikelihoodFit.mH'+str(mass)+'.root'
    inputFile = TFile(args.dir+'/'+args.tag+'/higgsCombineTest.MaxLikelihoodFit.mH'+str(mass)+'.root')
    tree = inputFile.Get('limit')
    # expected
    tree.GetEntry(0)
    limits[mass].append(tree.GetLeaf('limit').GetValue())
    #-s
    tree.GetEntry(1)
    limits[mass].append(tree.GetLeaf('limit').GetValue())
    #+s
    tree.GetEntry(2)
    limits[mass].append(tree.GetLeaf('limit').GetValue())
    # observed
    tree.GetEntry(3)
    limits[mass].append(tree.GetLeaf('limit').GetValue())

#print limits

drawBestFit(limits, args.tag)

