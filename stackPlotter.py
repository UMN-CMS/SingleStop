import ROOT
ROOT.gROOT.SetBatch(True)
bTags = '1B'
signalPath = 'output/313/{}'.format(bTags)
signalFiles = ['signal_2000_1900.root', 'signal_1500_900.root', 'signal_1000_400.root']
colors = [2, 3, 4]
backgroundPath = 'QCD2018_{}/output'.format(bTags)
backgroundFile = 'QCD2018.root'


ROOT.gStyle.SetPalette(ROOT.kSouthWest)

file1 = ROOT.TFile.Open('{}/{}'.format(signalPath, signalFiles[0]), "READ")
#file1 = ROOT.TFile.Open('output/313/signal_1500_600.root', "READ")
files = {fname: ROOT.TFile.Open('{}/{}'.format(signalPath, fname), "READ") for fname in signalFiles}
files.update( {backgroundFile: ROOT.TFile.Open('{}/{}'.format(backgroundPath, backgroundFile), "READ") } )

file1.cd("plots")
param_list = [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
print(param_list)
ROOT.TFile.Close(file1)

for param in param_list:
	h_stack = ROOT.THStack(param, "{};{};Events".format(param, param))
	legend = ROOT.TLegend(0.68, 0.89, 0.89, 0.65)
	#print(param)
	for i, f in enumerate(files):
		newfile = files[f]
		if f != "QCD2018.root": 
			newfile.cd("plots")	
			h_new = ROOT.gDirectory.Get(param)
		else:
			h_new = newfile.Get(param)
		if h_new == None or not (h_new.InheritsFrom('TH1F') or h_new.InheritsFrom('TH1D')): continue
		if param == 'cutflow': 
			h_new.Scale(1 / h_new.GetMaximum())
			print(h_new.GetBinContent(6), h_new.GetBinContent(7))
			print(1 - (h_new.GetBinContent(7) / h_new.GetBinContent(6)))			
		h_stack.Add(h_new)
		#print(h_new.GetMaximum())
		legend.AddEntry(h_new, f)
	c1 = ROOT.TCanvas()
	h_stack.Draw('HIST NOSTACK PLC')
	if param == 'cutflow':
		h_stack.SetMinimum(0)
		h_stack.SetMaximum(1)
	else:
		ROOT.gPad.SetLogy()
		h_stack.SetMinimum(10**-2)	
		h_stack.SetMaximum(10**8)
	legend.Draw()
	#h_stack.Draw('HIST')
	c1.SaveAs('plots/QCDOverlaid/{}/{}.png'.format(bTags, param))
