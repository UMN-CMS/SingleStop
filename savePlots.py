#!/usr/bin/env python
import os
import ROOT

doPlots = 1

tag = 'Run3_22-12-06'
masses = ['200_100','300_100','300_200','500_100','500_200','500_400','700_100','700_400','700_600','1000_100','1000_400','1000_900','1500_100','1500_600','1500_1400','2000_100','2000_900','2000_1900','700_200','1000_200','1500_200','1500_400','2000_200','2000_400']
#masses = ['500_100']

if doPlots:
  if not os.path.exists('plots/{}'.format(tag)): os.makedirs('plots/{}'.format(tag))
  for m in masses:
    if not os.path.exists('plots/{}/{}'.format(tag,m)): os.makedirs('plots/{}/{}'.format(tag,m))
    f = ROOT.TFile.Open('output/{}/hist_{}.root'.format(tag,m))
    plots = f.GetDirectory('plots')
    keys = [key.GetName() for key in plots.GetListOfKeys()]
    for key in keys:
      h = plots.Get(key)
      c1 = ROOT.TCanvas()
      h.Draw() if not h.InheritsFrom('TH2') else h.Draw('COLZ')
      c1.SaveAs('plots/{}/{}/{}.png'.format(tag,m,h.GetName()))
