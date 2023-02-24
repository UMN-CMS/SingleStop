import os, argparse, math, ROOT

path = 'root-test-files'
files = ['QCD.root', 'TT.root', 'signal_1500_600.root']
colors = [2, 3, 4]
h_stack = ROOT.THStack()
for i, f in enumerate(files):
	filepath = '{}/{}'.format(path, f)
	newfile = ROOT.TFile.Open(filepath, "READ")
	h_pT1Gen = newfile.Get('pT1Gen')
	h_pT1Gen.SetFillColor(colors[i])
	h_pT1Gen.SetMarkerColor(colors[i])
	h_pT1Gen.SetLineColor(colors[i])
	h_stack.Add(h_pT1Gen)		
	
c1 = ROOT.TCanvas()
h_stack.Draw()
c1.SaveAs('plots/test-stack.png')




 

