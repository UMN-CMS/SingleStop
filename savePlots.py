#!/usr/bin/env python
import os
import ROOT
import argparse

parser = argparse.ArgumentParser(description='Save plots from ROOT file')
parser.add_argument('--input',type=str,required=True,help='Path to input files')
parser.add_argument('--output',type=str,required=True,help='Path to output plots')
parser.add_argument('--sample',type=str,required=True,choices=['overlaid','signal','TT2018','QCD2018','QCDInclusive2018','ZQQ2018','ST2018','WQQ2018','ZNuNu2018','Diboson2018'],help='Sample to scale')
parser.add_argument('--TH1',action='store_true',default=False,help='Only save 1D histograms')
parser.add_argument('--TH2',action='store_true',default=False,help='Only save 2D histograms')
parser.add_argument('--norm',action='store_true',default=False,help='Normalize histograms')
args = parser.parse_args()

if (args.TH1 and args.TH2) or (not args.TH1 and not args.TH2): print('Saving all plots (1D and 2D)...')

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

if args.sample == 'overlaid':
  outputDir = '{}/overlaid{}'.format(args.output,'_norm' if args.norm else '')
  if not os.path.exists(outputDir): os.makedirs(outputDir)
  ROOT.gStyle.SetPalette(60) # kBlueRedYellow = 60
  ROOT.gStyle.SetOptTitle(0)
  samples = ['QCD2018','QCDInclusive2018','TT2018','signal_1200_400','signal_1500_900','signal_2000_1900']
  plotNames = ['cutflow','pT4ETFracUncomp','nQLHE','pTb3','dEtabb','dR35','dPhi24','dPhi23','dRRecoUncomp','HT','dR45','pT1Frac','dEta14','pT3ETFracComp','pT2FracChiComp','dR25','dR24','nWTight','pT4Frac','dR23','dPhiRecoComp','ntMedium','dPhi34','pTSDComp','pT5FracChiComp','pT5Frac','dEtaRecoComp','pT5ETFracComp','dR12','dR14','dR15','ntDeepWP4','ntDeepWP1','ntDeepWP3','ntDeepWP2','nWMedium','dPhi25','dR34','pTb4','nbMedium','pTb2','m3NoLead','pT1ETFracUncomp','dEta12','dEta13','dEta15','pTMeanUncomp','dPhi15','dPhi14','dPhi13','pT4FracChiUncomp','pT3ETFracUncomp','pT1FracChiComp','eta3','eta2','m3NoLeadOrSub','eta4','m4','pT2Frac','m3','pT4ETFracComp','dPhi12','nbTight','dR13','nJLHE','pT2ETFracUncomp','pT4FracChiComp','pT3Frac','dEta34','dEta35','dRRecoComp','dPhiRecoUncomp','nWDeepWP4','nWDeepWP1','nWDeepWP2','nWDeepWP3','ntTight','mAll','pT2FracChiUncomp','ntLoose','dRbb','dEta23','dEta25','dEta24','pT5ETFracUncomp','HTLHE','eta1','dPhibb','pT1FracChiUncomp','nWLoose','pTSDUncomp','dPhi35','pT4','dPhi45','nGLHE','pT1','nbLoose','pT3','dEtaRecoUncomp','pT5FracChiUncomp','pT3FracChiComp','pTMeanComp','nJets','MET','pT1ETFracComp','pT2','dEta45','pT2ETFracComp','pT3FracChiUncomp','pTb1','pTSDMeanFracComp','pTSDMeanFracUncomp']
  for plotName in plotNames:
    c1 = ROOT.TCanvas()
    stack = ROOT.THStack()
    legend = ROOT.TLegend(0.68,0.89,0.89,0.65)
    for sample in samples:
      if not os.path.isfile('{}/{}.root'.format(args.input,sample)): continue
      f = ROOT.TFile.Open('{}/{}.root'.format(args.input,sample))
      plot = f.Get(plotName)
      if plot.GetEntries() == 0: continue
      plot.SetTitle(sample)
      plot.SetLineWidth(3)
      plot.SetMarkerSize(0)
      if args.norm: 
        if plotName == 'cutflow': plot.Scale(1/plot.GetBinContent(1))
        else: plot.Scale(1/plot.Integral())
      ROOT.gROOT.cd()
      plot = plot.Clone()
      legend.AddEntry(plot,plot.GetTitle(),'l')
      stack.Add(plot)
    stack.Draw('HIST PLC NOSTACK')
    stack.GetXaxis().SetTitle(plot.GetXaxis().GetTitle()) 
    if args.norm: stack.GetYaxis().SetTitle('Fraction of Events')
    else:         stack.GetYaxis().SetTitle('Events')
    legend.Draw()
    ROOT.gPad.SetLogy()
    c1.SaveAs('{}/{}.pdf'.format(outputDir,plotName))

elif args.sample != 'signal':

  outputDir = '{}/{}{}'.format(args.output,args.sample,'_norm' if args.norm else '')
  if not os.path.exists(outputDir): os.makedirs(outputDir)
  f = ROOT.TFile.Open('{}/{}.root'.format(args.input,args.sample))
  keys = [key.GetName() for key in f.GetListOfKeys()]
  for key in keys:
    h = f.Get(key)
    if h.GetEntries() == 0: continue
    if args.TH1 and h.InheritsFrom('TH2'): continue
    if args.TH2 and not h.InheritsFrom('TH2'): continue
    c1 = ROOT.TCanvas()
    c1.SetRightMargin(0.15)
    if args.norm: h.Scale(1.0 / h.Integral())
    h.Draw() if not h.InheritsFrom('TH2') else h.Draw('COLZ')
    c1.SaveAs('{}/{}.pdf'.format(outputDir,h.GetName()))

else:

  masses = ['1000_400','1000_900','1500_600','1500_1400','2000_900','2000_1900','1000_600','1500_400','2000_400','2000_1400','1500_900','2000_600']
  for m in masses:
    outputDir = '{}/{}_{}{}'.format(args.output,args.sample,m,'_norm' if args.norm else '')
    if not os.path.exists(outputDir): os.makedirs(outputDir)
    f = ROOT.TFile.Open('{}/signal_{}.root'.format(args.input,m))
    keys = [key.GetName() for key in f.GetListOfKeys()]
    for key in keys:
      h = f.Get(key)
      if h.GetEntries() == 0: continue
      if args.TH1 and h.InheritsFrom('TH2'): continue
      if args.TH2 and not h.InheritsFrom('TH2'): continue
      c1 = ROOT.TCanvas()
      c1.SetRightMargin(0.15)
      if args.norm: h.Scale(1.0 / h.Integral())
      h.Draw() if not h.InheritsFrom('TH2') else h.Draw('COLZ')
      c1.SaveAs('{}/{}.pdf'.format(outputDir,h.GetName()))
