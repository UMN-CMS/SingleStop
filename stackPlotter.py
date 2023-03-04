import os, argparse, math, ROOT
ROOT.gROOT.SetBatch(True)
path = 'root-test-files'
files = ['QCD.root', 'TT.root', 'signal_1500_600.root']
colors = [2, 3, 4]


ROOT.gStyle.SetPalette(ROOT.kOcean)

file1 = ROOT.TFile.Open('root-test-files/signal_1500_600.root', "READ")
param_list = [key.GetName() for key in file1.GetListOfKeys()]
ROOT.TFile.Close(file1)

for param in param_list:
	h_stack = ROOT.THStack(param, "{};{};Events".format(param, param))
	print(param)
	for i, f in enumerate(files):
		filepath = '{}/{}'.format(path, f)
		newfile = ROOT.TFile.Open(filepath, "READ")	
		h_new = newfile.Get(param)
		if not h_new:
			break
		#ROOT.TFile.Close(newfile)
		#ROOT.gStyle.SetPalette(3)	
		h_stack.Add(h_new)	
	if not h_new or not h_new.InheritsFrom('TH1'):
		continue
	c1 = ROOT.TCanvas()
	ROOT.gPad.SetLogy()
	h_stack.Draw('HIST')
	c1.SaveAs('plots/{}.png'.format(param))
