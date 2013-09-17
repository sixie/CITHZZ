#! /usr/bin/env python

from ROOT import TFile, TH1D, TCanvas, TChain, THStack, TLegend, TGraph, TLine, TMultiGraph, TLatex
from array import array
import rootlogonTDR

def drawLimits(limits,tag):

    masses = []
    reverseMasses = []
    massesForSigma = []
    medians = []
    sigma1l = []
    sigma2l = []
    sigma1 = []
    sigma2 = []
    observed = []
    
    #-2sigma -1sigma  mean  +1sigma  +2sigma,  median
    for mass in sorted(limits):
        masses.append(float(mass))
        medians.append(float(limits[mass][2]))

        reverseMasses.append(float(mass))
        massesForSigma.append(float(mass))
        sigma1l.append(float(limits[mass][1]))
        sigma1.append(float(limits[mass][3]))
        sigma2l.append(float(limits[mass][0]))
        sigma2.append(float(limits[mass][4]))

        observed.append(float(limits[mass][5]))

    for m in reverseMasses:
        print masses[masses.index(m)],observed[masses.index(m)],medians[masses.index(m)],medians[masses.index(m)],sigma2l[masses.index(m)],sigma1l[masses.index(m)],sigma1[masses.index(m)],sigma2[masses.index(m)]

    reverseMasses.reverse()

    for m in reverseMasses:
        massesForSigma.append(m)
        sigma1.append(sigma1l[masses.index(m)])
        sigma2.append(sigma2l[masses.index(m)])

    multiGraph = TMultiGraph()
    expectedLimits = TGraph(len(masses),array('f',masses),array('f', medians))
    expectedLimits.SetLineWidth(1)
    expectedLimits.SetLineStyle(2)

    band1Sigma = TGraph(len(massesForSigma),array('f',massesForSigma),array('f',sigma1))
    band1Sigma.SetFillColor(416)
    band1Sigma.SetLineColor(416)

    band2Sigma = TGraph(len(massesForSigma),array('f',massesForSigma),array('f',sigma2))
    band2Sigma.SetFillColor(400)
    band2Sigma.SetLineColor(400)

    observedLimitsGraph = TGraph(len(masses),array('f',masses),array('f', observed))
    observedLimitsGraph.SetLineWidth(1)
    observedLimitsGraph.SetLineStyle(2)


    legend = TLegend(0.3,0.7,0.6,0.9)
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    legend.SetTextAlign(12);
    legend.SetTextFont(42);
    legend.SetTextSize(0.03);
    legend.AddEntry(observedLimitsGraph,"95% CL exclusion: observed", "l")
    legend.AddEntry(expectedLimits,"95% CL exclusion: median", "l")
    legend.AddEntry(band1Sigma,"95% CL exclusion: 1 #sigma band", "f")
    legend.AddEntry(band2Sigma,"95% CL exclusion: 2 #sigma band", "f")

    lineOne = TLine(114,1, masses[-1]+1, 1);
    lineOne.SetLineWidth(2);
    lineOne.SetLineStyle(1);


    cLimits = TCanvas('limits')
    frame = TH1D('frame','',1000,114,masses[-1]+1)
    frame.SetMinimum(0)
    frame.SetMaximum(7.5)
    frame.SetDirectory(0)
    frame.SetStats(0)
    frame.SetFillColor(63)
    frame.SetLineStyle(0)
    frame.SetMarkerStyle(20)
    frame.GetYaxis().SetLabelSize(0.05);
    frame.GetYaxis().SetTitle('95% CL Limit on #sigma/#sigma_{SM}')
    frame.GetXaxis().SetTitle('m_{H} (GeV)')
    frame.Draw('  ')
    
    #cLimits.DrawFrame(114,0,masses[-1]+1,7.5,'')
    band2Sigma.GetXaxis().SetTitle('M_{H} [GeV/c^{2}]')
    band2Sigma.GetYaxis().SetTitle('95% CL Limit on #sigma/#sigma_{SM}')
    band2Sigma.Draw('f,same')
    band1Sigma.GetXaxis().SetTitle('M_{H} [GeV/c^{2}]')
    band1Sigma.GetYaxis().SetTitle('95% CL Limit on #sigma/#sigma_{SM}')
    band1Sigma.Draw('f,same')
    expectedLimits.GetXaxis().SetTitle('M_{H} [GeV/c^{2}]')
    expectedLimits.GetYaxis().SetTitle('95% CL Limit on #sigma/#sigma_{SM}')
    expectedLimits.Draw('L,same')
    observedLimitsGraph.SetLineColor(4)
    observedLimitsGraph.SetMarkerColor(4)
    observedLimitsGraph.SetLineStyle(1)
    observedLimitsGraph.GetXaxis().SetTitle('M_{H} [GeV/c^{2}]')
    observedLimitsGraph.GetYaxis().SetTitle('95% CL Limit on #sigma/#sigma_{SM}')
    observedLimitsGraph.Draw('LP,same')

    texCMS = TLatex(120, 0.9*frame.GetMaximum(), 'CMS Preliminary')
    texCMS.SetTextSize(14)
    #texCMS.Draw()
    if tag.find('4mu')>-1: modeLabel = 'H#rightarrow ZZ #rigtharrow 4#mu, (2D)'
    elif tag.find('4e')>-1: modeLabel = 'H#rightarrow ZZ #rigtharrow 4e, (2D)'
    elif tag.find('2e2mu')>-1: modeLabel = 'H#rightarrow ZZ #rigtharrow 2e2#mu, (2D)'
    else: modeLabel = 'H#rightarrow ZZ #rigtharrow 4l, (2D)'
    texMode = TLatex(120, 0.8*frame.GetMaximum(), modeLabel)    
    texMode.SetTextSize(14)
    #texMode.Draw()
    if tag.find('7TeV')>-1: lumi = 'L=5.1 fb^{-1}'
    elif tag.find('8TeV')>-1: lumi = 'L=12.1 fb^{-1}'
    elif tag.find('7p8TeV')>-1: lumi='L=5.1 fb^{-1} (7 TeV) + 12.1 fb^{-1} (8 TeV)'
    else: lumi='L=5.1 fb^{-1} (7 TeV) + 12.1 fb^{-1} (8 TeV)'
    texLumi = TLatex(120, 0.7*frame.GetMaximum(), lumi)
    texLumi.SetTextSize(14)
    #texLumi.Draw()

    legend.Draw()
    lineOne.Draw('same')
    cLimits.Update()

    #cLimits.SaveAs('limits'+tag+'.png')

    output = TFile('limits.root','RECREATE')
    cLimits.Write()

    raw_input('type whatever to quit')


def drawBestFit(limits,tag):

    masses = []
    reverseMasses = []
    massesForSigma = []
    medians = []
    sigma1l = []
    sigma1 = []
    observed = []

    #mean -1sigma  +1sigma mean
    for mass in sorted(limits):
        masses.append(float(mass))
        medians.append(float(limits[mass][3]))

        reverseMasses.append(float(mass))
        massesForSigma.append(float(mass))
        sigma1l.append(float(limits[mass][1]))
        sigma1.append(float(limits[mass][2]))


    reverseMasses.reverse()

    print medians

    for m in reverseMasses:
        massesForSigma.append(m)
        sigma1.append(sigma1l[masses.index(m)])

    multiGraph = TMultiGraph()
    expectedLimits = TGraph(len(masses),array('f',masses),array('f', medians))
    expectedLimits.SetLineWidth(1)
    expectedLimits.SetLineStyle(2)

    band1Sigma = TGraph(len(massesForSigma),array('f',massesForSigma),array('f',sigma1))
    band1Sigma.SetFillColor(416)
    band1Sigma.SetLineColor(416)

    legend = TLegend(0.5,0.4 ,0.92,0.7)
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    legend.AddEntry(expectedLimits,"median", "l")
    legend.AddEntry(band1Sigma,"1 #sigma band", "f")

    lineOne = TLine(115,1, masses[-1]+1, 1);
    lineOne.SetLineWidth(2);
    lineOne.SetLineStyle(1);
    lineOne.SetLineColor(2);
    
    cLimits = TCanvas('limits')
    cLimits.SetGridy()
    frame = TH1D('frame','',1000,115,masses[-1]+1)
    frame.SetMinimum(-2)
    frame.SetMaximum(3)
    frame.SetDirectory(0)
    frame.SetStats(0)
    frame.SetFillColor(63)
    frame.SetLineStyle(0)
    frame.SetMarkerStyle(20)
    frame.GetYaxis().SetLabelSize(0.05);
    frame.GetYaxis().SetTitle('best fit for #sigma/#sigma_{SM}')
    frame.GetXaxis().SetTitle('m_{H} (GeV)')
    frame.Draw('  ')
    
    #cLimits.DrawFrame(114,0,masses[-1]+1,7.5,'')
    band1Sigma.GetXaxis().SetTitle('M_{H} [GeV/c^{2}]')
    band1Sigma.GetYaxis().SetTitle('95% CL Limit on #sigma/#sigma_{SM}')
    band1Sigma.Draw('f,same')
    expectedLimits.SetLineWidth(2)
    expectedLimits.GetXaxis().SetTitle('M_{H} [GeV/c^{2}]')
    expectedLimits.GetYaxis().SetTitle('95% CL Limit on #sigma/#sigma_{SM}')
    expectedLimits.Draw('L,same')

    xCorner=140
    yCorner=2.75
    text = [TLatex(xCorner,yCorner,'CMS Preliminary')]
    if tag == '7TeV':
        text.append(TLatex(xCorner,yCorner-0.25,'#sqrt{s} = 7 TeV, L = 5.1 fb^{-1}'))
    elif tag == '8TeV':
        text.append(TLatex(xCorner,yCorner-0.25,'#sqrt{s} = 8 TeV, L = 12.1 fb^{-1}'))
    else:
        text.append(TLatex(xCorner,yCorner-0.25,'#sqrt{s} = 7 TeV, L = 5.1 fb^{-1}'))
        text.append(TLatex(xCorner,yCorner-0.5,'#sqrt{s} = 8 TeV, L = 12.1 fb^{-1}'))

                    
    #legend.Draw()
    for t in text:
        t.SetTextFont(42)
        t.SetTextSize(0.05)
        t.Draw('same')
    lineOne.Draw('same')
    cLimits.Update()

    cLimits.SaveAs('mlf'+tag+'.png')
    cLimits.SaveAs('mlf'+tag+'.pdf')

    output = TFile('mlf.root','RECREATE')
    cLimits.Write()

    raw_input('type whatever to quit')

