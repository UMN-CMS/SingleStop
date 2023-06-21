import ROOT, os
ROOT.gROOT.SetBatch(True)
dirs = ['200']#, '1B', '200']
backgroundFiles = ['QCD2018.root', 'QCDInclusive2018.root']
ROOT.gStyle.SetPalette(ROOT.kSouthWest)

files = {dirName: ROOT.TFile.Open('QCD2018_{}/output/{}'.format(dirName, backgroundFiles[0]), "READ") for dirName in dirs}
files.update( {'{}Inc'.format(dirName): ROOT.TFile.Open('QCD2018Inc_{}/output/{}'.format(dirName, backgroundFiles[1]), "READ") for dirName in dirs} )
param_list = [key.GetName() for key in files[dirs[0]].GetListOfKeys()]
param_list.remove('pTGenBQCD')
param_list.remove('pTRecoB')

for param in param_list:
	h_stack = ROOT.THStack(param, "{};{};Events".format(param, files[dirs[0]].Get(param).GetXaxis().GetTitle()))
	legend = ROOT.TLegend(0.68, 0.89, 0.89, 0.65)
	#print(param)
	for i, f in enumerate(files):
		newfile = files[f]
		h_new = newfile.Get(param)
		if h_new == None or not (h_new.InheritsFrom('TH1F') or h_new.InheritsFrom('TH1D')) or h_new.Integral() == 0.0: continue
		if param == 'pTRecoBMatchSuccess':
			h_new.Divide(files[f].Get('pTRecoB'))
		elif param == 'pTGenBMatchSuccess':
			h_new.Divide(files[f].Get('pTGenBQCD'))
		if param == 'cutflow': 
			h_new.Scale(1 / h_new.GetMaximum())
			print(h_new.GetBinContent(6), h_new.GetBinContent(7))
			print(1 - (h_new.GetBinContent(7) / h_new.GetBinContent(6)))
		elif param != 'pTRecoBMatchSuccess' and param != 'pTGenBMatchSuccess':
			h_new.Scale(1 / h_new.Integral())			
		h_new.SetLineWidth(3)
		h_stack.Add(h_new)
		#print(h_new.GetMaximum())
		label = 'QCD2018 [total]' if ( 'Inc' in f ) else 'QCD2018 [b-enriched]'
		legend.AddEntry(h_new, label)
	if h_new == None or h_new.InheritsFrom('TH2D') or h_new.GetEntries() == 0: continue	
	c1 = ROOT.TCanvas()
	h_stack.Draw('HIST NOSTACK PLC E')
	legend.Draw()
	#h_stack.Draw('HIST')
	print(h_stack)
	outputDir = 'plots/QCDOverlaid/ControlRegions/{}'.format(dirs[0], param)
	if not os.path.exists(outputDir): os.makedirs(outputDir)
	c1.SaveAs('{}.png'.format(outputDir))



