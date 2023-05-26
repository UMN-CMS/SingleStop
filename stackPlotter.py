import ROOT, argparse
ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description = 'Stack Plotter')
parser.add_argument('--bTags', type = str, default = '300')
args = parser.parse_args()
bTags = args.bTags

signalPath = 'output/313/{}/scaled'.format(bTags)
signalFiles = ['signal_2000_1900.root', 'signal_1500_900.root', 'signal_1000_400.root']
colors = [2, 3, 4]
backgroundPath = 'QCD2018_{}/output'.format(bTags)
backgroundFile = 'QCD2018.root'


ROOT.gStyle.SetPalette(ROOT.kRainBow)

file1 = ROOT.TFile.Open('{}/{}'.format(signalPath, signalFiles[0]), "READ")
files = {fname: ROOT.TFile.Open('{}/{}'.format(signalPath, fname), "READ") for fname in signalFiles}
files.update( {backgroundFile: ROOT.TFile.Open('{}/{}'.format(backgroundPath, backgroundFile), "READ") } )

param_list = [key.GetName() for key in file1.GetListOfKeys()]
print(param_list)
ROOT.TFile.Close(file1)

bb2DHists = []
for param in param_list:
	h_stack = ROOT.THStack(param, "{};{};Events".format(param, param))
	legend = ROOT.TLegend(0.68, 0.89, 0.89, 0.65)
	for i, f in enumerate(files):
		newfile = files[f]
		h_new = newfile.Get(param)
		if 'bb' in param and h_new.InheritsFrom('TH2D'): bb2DHists.append(param)
		if h_new == None or not (h_new.InheritsFrom('TH1F') or h_new.InheritsFrom('TH1D')): continue
		if param == 'cutflow': 
			h_new.Scale(1 / h_new.GetMaximum())
		#elif h_new.Integral() != 0: h_new.Scale(1 / h_new.Integral())
		h_new.SetLineWidth(3)
		h_stack.Add(h_new)
		legend.AddEntry(h_new, f)
	if h_new == None or h_new.InheritsFrom('TH2D') or h_new.GetEntries() == 0: continue
	c1 = ROOT.TCanvas()
	h_stack.Draw('HIST NOSTACK PLC')
	if param == 'cutflow':
		h_stack.SetMinimum(0)
		h_stack.SetMaximum(1)
	else: ROOT.gPad.SetLogy()
	legend.Draw()
	c1.SaveAs('plots/QCDOverlaid/{}/{}.png'.format(bTags, param))
'''
bb2DHists = list(set(bb2DHists))
for param in bb2DHists:
	for f in files:
		newfile = files[f]
		if f != 'QCD2018.root':
			newfile.cd('plots')
			h_new = ROOT.gDirectory.Get(param)
		else:
			h_new = newfile.Get(param)
		h_new.SetTitle('{}'.format(f[:len(f) - 5]))
		if h_new.Integral() != 0: h_new.Scale(1 / h_new.Integral())
		h_new.SetStats(0)
		c1 = ROOT.TCanvas()
		c1.SetRightMargin(0.15)
		h_new.Draw('COLZ1')
		c1.SaveAs('plots/QCDOverlaid/{}/{}_{}.png'.format(bTags, param, f[:len(f) - 5]))
'''		
