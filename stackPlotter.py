import os, argparse, math, ROOT
ROOT.gROOT.SetBatch(True)
path = 'output/signal_2000_1900_1500_600_1000_400'
files = ['signal_2000_1900.root', 'signal_1500_600.root', 'signal_1000_400.root']
colors = [2, 3, 4]


ROOT.gStyle.SetPalette(ROOT.kOcean)

<<<<<<< HEAD
file1 = ROOT.TFile.Open('{}/{}'.format(path, files[0]), "READ")
=======
file1 = ROOT.TFile.Open('root-test-files/signal_1500_600.root', "READ")
>>>>>>> Stack Plotter as of Feb. 28 Meeting
param_list = [key.GetName() for key in file1.GetListOfKeys()]
ROOT.TFile.Close(file1)
files = {fname: ROOT.TFile.Open('{}/{}'.format(path, fname), "READ") for fname in files}

for param in param_list:
	h_stack = ROOT.THStack(param, "{};{};Events".format(param, param))
	print(param)
	for i, f in enumerate(files):
		newfile = files[f]	
		h_new = newfile.Get(param)
		if h_new == None or not h_new.InheritsFrom('TH1F'):
			continue	
		h_stack.Add(h_new)
	c1 = ROOT.TCanvas()
	h_stack.Draw('HIST')
	ROOT.gPad.SetLogy()
	h_stack.SetMinimum(1)	
	h_stack.SetMaximum(10**8)
	#h_stack.Draw('HIST')
	c1.SaveAs('plots/313/{}.png'.format(param))
