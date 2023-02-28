import os, argparse, math, ROOT

path = 'root-test-files'
files = ['QCD.root', 'TT.root', 'signal_1500_600.root']
colors = [2, 3, 4]


ROOT.gStyle.SetPalette(ROOT.kOcean)

file1 = ROOT.TFile.Open('root-test-files/signal_1500_600.root', "READ")
param_list = [key.GetName() for key in file1.GetListOfKeys()]
ROOT.TFile.Close(file1)

for param in param_list:
	h_stack = ROOT.THStack()
	for i, f in enumerate(files):
		filepath = '{}/{}'.format(path, f)
		newfile = ROOT.TFile.Open(filepath, "READ")	
		h_new = newfile.Get(param)
		ROOT.TFile.Close(newfile)
		#ROOT.gStyle.SetPalette(3)	
		h_stack.Add(h_new)	

	c1 = ROOT.TCanvas()
	h_stack.Draw('HIST')
	c1.SaveAs('plots/{}.png'.format(param))
