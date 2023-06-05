import ROOT, argparse
ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description = 'Stack Plotter')
parser.add_argument('--bTags', type = str, default = '300')
args = parser.parse_args()
bTags = args.bTags

signalPath = 'output/313/{}/scaled'.format(bTags)
signalFiles = ['signal_2000_1900.root', 'signal_1500_900.root', 'signal_1000_400.root']
colors = [2, 3, 4]
QCDPath = 'QCD2018_{}/output'.format(bTags)
QCDFile = 'QCD2018.root'
TTPath = 'TT2018_{}/output'.format(bTags)
TTFile = 'TT2018.root'


ROOT.gStyle.SetPalette(ROOT.kRainBow)

file1 = ROOT.TFile.Open('{}/{}'.format(signalPath, signalFiles[0]), "READ")
files = {fname: ROOT.TFile.Open('{}/{}'.format(signalPath, fname), "READ") for fname in signalFiles}
files.update( {QCDFile: ROOT.TFile.Open('{}/{}'.format(QCDPath, QCDFile), "READ") } )
#files.update( {TTFile: ROOT.TFile.Open('{}/{}'.format(TTPath, TTFile), "READ") } )

param_list = [key.GetName() for key in file1.GetListOfKeys()]
#param_list = ['cutflow', 'dEtabb12', 'dEtabb13', 'dEtabb23', 'dPhibb12', 'dPhibb13', 'dPhibb23', 'dRbb12', 'dRbb13', 'dRbb23']
#print(param_list)
ROOT.TFile.Close(file1)
#labels = {'QCD2018.root': 'QCD 2018', 'TT2018.root': r'$t \overline{t}$ 2018', 'signal_1000_400.root': r'$\tilde{t}$ = 1000, $\tilde{chi}^{\pm}$ = 400', 'signal_1500_900.root': r'$\tilde{t}$ = 1500, $\tilde{chi}^{\pm}$ = 900', 'signal_2000_1900.root': r'$\tilde{t}$ = 2000, $\tilde{chi}^{\pm}$ = 1900'}

bb2DHists = []
for param in param_list:
	h_stack = ROOT.THStack(param, "{};{};Events".format(param, files['signal_2000_1900.root'].Get(param).GetXaxis().GetTitle()))
	legend = ROOT.TLegend(0.68, 0.89, 0.89, 0.65)
	for i, f in enumerate(files):
		newfile = files[f]
		h_new = newfile.Get(param)
		#if 'bb' in param and h_new.InheritsFrom('TH2D'): bb2DHists.append(param)
		if 'bb' not in param and h_new.InheritsFrom('TH2D'): bb2DHists.append(param)
		if h_new == None or not (h_new.InheritsFrom('TH1F') or h_new.InheritsFrom('TH1D')): continue
		if param == 'cutflow': 
			print(f, h_new.GetBinContent(6), h_new.GetBinContent(7), 1 - h_new.GetBinContent(7) / h_new.GetBinContent(6))				
			h_new.Scale(1 / h_new.GetMaximum())
		elif h_new.Integral() != 0: h_new.Scale(1 / h_new.Integral())
		h_new.SetLineWidth(3)
		h_stack.Add(h_new)
		legend.AddEntry(h_new, f)
	if h_new == None or h_new.InheritsFrom('TH2D') or h_new.GetEntries() == 0: continue
#'''	
	c1 = ROOT.TCanvas()
	h_stack.Draw('HIST NOSTACK PLC')
	if param == 'cutflow':
		h_stack.SetMinimum(0)
		h_stack.SetMaximum(1)
	legend.Draw()
	c1.SaveAs('plots/QCDOverlaid/{}/{}.png'.format(bTags, param))
'''
bb2DHists = list(set(bb2DHists))
print(bb2DHists)
for param in bb2DHists:
	for f in files:
		newfile = files[f]
		h_new = newfile.Get(param)
		h_new.SetTitle('{}'.format(f[:len(f) - 5]))
		if h_new.Integral() != 0: h_new.Scale(1 / h_new.Integral())
		if param == 'dPhi13Vs34': print(f, h_new.Integral(h_new.GetXaxis().FindBin(0), h_new.GetXaxis().FindBin(0.7999), h_new.GetYaxis().FindBin(2.6), h_new.GetYaxis().FindBin(3.1999)))
		h_new.SetStats(0)
		c1 = ROOT.TCanvas()
		c1.SetRightMargin(0.15)
		h_new.Draw('COLZ1')
		c1.SaveAs('plots/QCDOverlaid/{}/{}_{}.png'.format(bTags, param, f[:len(f) - 5]))
'''		
