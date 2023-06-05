import ROOT
ROOT.gROOT.SetBatch(True)
dirs = ['0B', '1B', 'RestrictedpT', '300']
backgroundFile = 'QCD2018.root'
ROOT.gStyle.SetPalette(1)

files = {dirName: ROOT.TFile.Open('QCD2018_{}/output/{}'.format(dirName, backgroundFile), "READ") for dirName in dirs}

#param_list = [key.GetName() for key in files['0B'].GetListOfKeys()]
#print(param_list)
param_list = ['m3', 'm4']

for param in param_list:
	h_stack = ROOT.THStack(param, "{};{};Events".format(param, param))
	legend = ROOT.TLegend(0.68, 0.89, 0.89, 0.65)
	#print(param)
	for i, f in enumerate(files):
		newfile = files[f]
		h_new = newfile.Get(param)
		if h_new == None or not (h_new.InheritsFrom('TH1F') or h_new.InheritsFrom('TH1D')) or h_new.Integral() == 0.0: continue
		if param == 'cutflow': 
			h_new.Scale(1 / h_new.GetMaximum())
			print(h_new.GetBinContent(6), h_new.GetBinContent(7))
			print(1 - (h_new.GetBinContent(7) / h_new.GetBinContent(6)))
		else:
			h_new.Scale(1 / h_new.Integral())			
		h_stack.Add(h_new)
		#print(h_new.GetMaximum())
		legend.AddEntry(h_new, f)
	c1 = ROOT.TCanvas()
	h_stack.Draw('HIST NOSTACK PLC')
	legend.Draw()
	#h_stack.Draw('HIST')
	print(h_stack)
	c1.SaveAs('plots/QCDOverlaid/ControlRegions/{}.png'.format(param))
