import os, argparse, math, ROOT

path = 'root-test-files'
files = ['QCD.root', 'TT.root', 'signal_1500_600.root']
colors = [2, 3, 4]
h_stack = ROOT.THStack()

ROOT.gStyle.SetPalette(ROOT.kOcean)

for i, f in enumerate(files):
	filepath = '{}/{}'.format(path, f)
	newfile = ROOT.TFile.Open(filepath, "READ")
	h_pT1Gen = newfile.Get('pT1Gen')
	#ROOT.gStyle.SetPalette(3)	
	h_stack.Add(h_pT1Gen)	

c1 = ROOT.TCanvas()
h_stack.Draw('HIST')
c1.SaveAs('plots/test-stack.png')
