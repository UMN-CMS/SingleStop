#!/usr/bin/env python
import os
import ROOT
import argparse

parser = argparse.ArgumentParser(description='Save plots from ROOT file')
parser.add_argument('--input',type=str,required=True,help='Path to input files')
parser.add_argument('--output',type=str,required=True,help='Path to output plots')
parser.add_argument('--sample',type=str,required=True,choices=['overlaid','signal','TT2018','QCD2018','ZQQ2018','ST2018','WQQ2018','ZNuNu2018','Diboson2018'],help='Sample to scale')
parser.add_argument('--TH1',action='store_true',default=False,help='Only save 1D histograms')
parser.add_argument('--TH2',action='store_true',default=False,help='Only save 2D histograms')
parser.add_argument('--norm',action='store_true',default=False,help='Normalize histograms')
args = parser.parse_args()

if (args.TH1 and args.TH2) or (not args.TH1 and not args.TH2): print('Saving all plots (1D and 2D)...')

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

if args.sample == 'overlaid':
  if not os.path.exists('{}/overlaid'.format(args.output)): os.makedirs('{}/overlaid'.format(args.output))
  ROOT.gStyle.SetPalette(60) # kBlueRedYellow = 60
  ROOT.gStyle.SetOptTitle(0)
  samples = ['QCD2018','signal_1000_400','signal_1200_400','signal_1400_400']
  plotNames = ['m4']
  for plotName in plotNames:
    c1 = ROOT.TCanvas()
    stack = ROOT.THStack()
    legend = ROOT.TLegend(0.68,0.89,0.89,0.65)
    for sample in samples:
      f = ROOT.TFile.Open('{}/{}.root'.format(args.input,sample))
      plot = f.Get(plotName)
      plot.SetTitle(sample)
      plot.SetLineWidth(3)
      plot.SetMarkerSize(0)
      if args.norm: plot.Scale(1/plot.Integral())
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
    c1.SaveAs('{}/overlaid/{}.pdf'.format(args.output,plotName))

elif args.sample != 'signal':

  if not os.path.exists('{}/{}'.format(args.output,args.sample)): os.makedirs('{}/{}'.format(args.output,args.sample))
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
    c1.SaveAs('{}/{}/{}.pdf'.format(args.output,args.sample,h.GetName()))

else:

  masses = ['1000_400','1000_900','1500_600','1500_1400','2000_900','2000_1900','1000_600','1500_400','2000_400','2000_1400','1500_900','2000_600']
  for m in masses:
    if not os.path.exists('{}/{}_{}'.format(args.output,args.sample,m)): os.makedirs('{}/{}_{}'.format(args.output,args.sample,m))
    f = ROOT.TFile.Open('{}/signal_{}.root'.format(args.input,m))
    keys = [key.GetName() for key in f.GetListOfKeys()]
    for key in keys:
      h = f.Get(key)
      if args.TH1 and h.InheritsFrom('TH2'): continue
      if args.TH2 and not h.InheritsFrom('TH2'): continue
      c1 = ROOT.TCanvas()
      c1.SetRightMargin(0.15)
      if args.norm: h.Scale(1.0 / h.Integral())
      h.Draw() if not h.InheritsFrom('TH2') else h.Draw('COLZ')
      c1.SaveAs('{}/{}_{}/{}.png'.format(args.output,args.sample,m,h.GetName()))
