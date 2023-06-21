import ROOT, argparse, os, math
ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description = 'Stack Plotter')
parser.add_argument('--bTags', type = str, default = '300')
parser.add_argument('--TH1', default = True)
parser.add_argument('--TH2', default = False)

args = parser.parse_args()
bTags = args.bTags
TH1 = args.TH1
TH2 = args.TH2

signalPath = 'output/313/{}/scaled'.format(bTags)
signalFiles = ['signal_2000_1900.root', 'signal_1500_900.root', 'signal_1000_400.root']
colors = [2, 3, 4]
QCDPath = 'QCD2018_{}/output'.format(bTags)
QCDFile = 'QCD2018.root'
TTPath = 'TT2018_{}/output'.format(bTags)
TTFile = 'TT2018.root'
DataPath = 'Data2018_{}/output'.format(bTags)
DataFile = 'Data2018.root'

if not os.path.exists('plots/QCDOverlaid/{}/313'.format(bTags)): os.makedirs('plots/QCDOverlaid/{}/313'.format(bTags))

ROOT.gStyle.SetPalette(ROOT.kRainBow)

file1 = ROOT.TFile.Open('{}/{}'.format(signalPath, signalFiles[0]), "READ")
files = {fname: ROOT.TFile.Open('{}/{}'.format(signalPath, fname), "READ") for fname in signalFiles}
#files.update( { QCDFile: ROOT.TFile.Open('{}/{}'.format(QCDPath, QCDFile), "READ") } )
#files.update( { TTFile: ROOT.TFile.Open('{}/{}'.format(TTPath, TTFile), "READ") } )
files.update( { 'QCDInclusive2018.root': ROOT.TFile.Open('QCD2018Inc_{}/output/QCDInclusive2018.root'.format(bTags), "READ") } )
files.update( { 'Data2018.root': ROOT.TFile.Open('Data2018_{}/output/Data2018.root'.format(bTags), "READ") } )

#param_list = [key.GetName() for key in file1.GetListOfKeys()]
#ROOT.TFile.Close(file1)
param_list = ['m3', 'm4', 'm3NoLead']
#param_list = ['HT', 'pT1']
#labels = {'QCD2018.root': 'QCD 2018', 'TT2018.root': r'$t \overline{t}$ 2018', 'signal_1000_400.root': r'$\tilde{t}$ = 1000, $\tilde{chi}^{\pm}$ = 400', 'signal_1500_900.root': r'$\tilde{t}$ = 1500, $\tilde{chi}^{\pm}$ = 900', 'signal_2000_1900.root': r'$\tilde{t}$ = 2000, $\tilde{chi}^{\pm}$ = 1900'}
#param_list = ['bJetMatchingEff', 'etaGenBQCD', 'genBMatchingRate', 'matchedPartonFlavour', 'nBsPerQCDEvent', 'phiGenBQCD', 'pTRecoBMatchSuccess', 'QCDGenpdgId']

bb2DHists = []
for param in param_list:
	h_stack = ROOT.THStack(param, "{};{};Events".format(param, files['signal_2000_1900.root'].Get(param).GetXaxis().GetTitle()))
	legend = ROOT.TLegend(0.68, 0.89, 0.89, 0.65)
	for i, f in enumerate(files):
		newfile = files[f]
		h_new = newfile.Get(param)
		if 'pTRecoB' in param or 'pTGenB' in param:
			h_new.SetLineWidth(3)
			continue
		#if 'bb' in param and h_new.InheritsFrom('TH2D'): bb2DHists.append(param)
		if 'bb' not in param and h_new != None and h_new.InheritsFrom('TH2D'): bb2DHists.append(param)
		if h_new == None or not (h_new.InheritsFrom('TH1F') or h_new.InheritsFrom('TH1D')): continue
		if param == 'cutflow': 
			print(f, h_new.GetBinContent(6), h_new.GetBinContent(7), 1 - h_new.GetBinContent(7) / h_new.GetBinContent(6))				
			h_new.Scale(1 / h_new.GetMaximum())
		#elif h_new.Integral() != 0: h_new.Scale(1 / h_new.Integral())
		h_new.SetLineWidth(3)
		h_stack.Add(h_new)
		legend.AddEntry(h_new, f[:-5])
	if h_new == None or h_new.InheritsFrom('TH2D') or h_new.GetEntries() == 0: continue
	if args.TH1:		
		c1 = ROOT.TCanvas()
		h_stack.Draw('HIST E NOSTACK PLC')
		if param == 'cutflow':
			h_stack.SetMinimum(0)
			h_stack.SetMaximum(1)
		legend.Draw()
		if param != 'cutflow':
			ROOT.gPad.SetLogy()
			h_stack.SetMaximum(10**7)
			h_stack.SetMinimum(10**0)
		c1.SaveAs('plots/QCDOverlaid/{}/313/{}.png'.format(bTags, param))

for param in ['m3', 'm4', 'm3NoLead']:
	for f in signalFiles:
		legend = ROOT.TLegend(0.68, 0.89, 0.89, 0.75)
		h_signal = files[f].Get(param)
		h_background = files['QCDInclusive2018.root'].Get(param)
		for i in range(h_background.GetNcells()):
			if h_background.GetBinContent(i) != 0: h_signal.SetBinContent(i, h_signal.GetBinContent(i) / h_background.GetBinContent(i))
			else: h_signal.SetBinContent(i, 0)
		c1 = ROOT.TCanvas()
		h_signal.SetYTitle('s / b')
		h_signal.SetStats(0)
		legend.AddEntry(h_new, f)
		h_signal.Draw('HIST')
		h_signal.Rebin(4)
		legend.Draw()
		c1.SaveAs('plots/QCDOverlaid/{}/313/{}Significance_{}.png'.format(bTags, param, f[7:-5]))	
for param in ['HT', 'pT1']:
	for f in signalFiles:
		legend = ROOT.TLegend(0.68, 0.89, 0.89, 0.65)

		h_signal = files[f].Get(param)
		h_background = files['QCDInclusive2018.root'].Get(param)
		h_background.Add(h_signal)
		h_data = files['Data2018.root'].Get(param)

		h_stack = ROOT.THStack(param, "{} (Data Vs. MC {});{};Events".format(param, f[7:-5], files[f].Get(param).GetXaxis().GetTitle()))
		h_stack.Add(h_background)
		h_stack.Add(h_data)
		legend.AddEntry(h_background, 'MC ({})'.format(f[7:-5]))
		legend.AddEntry(h_data, 'Data2018')


		c1 = ROOT.TCanvas()
		h_stack.Draw('HIST E NOSTACK PLC')
		ROOT.gPad.SetLogy()
		h_stack.SetMaximum(10**7)
		h_stack.SetMinimum(10**0)
		legend.Draw()
		c1.SaveAs('plots/QCDOverlaid/{}/313/{}DataVsMC_{}.png'.format(bTags, param, f[7:-5]))

'''
c1 = ROOT.TCanvas()
files['QCD2018.root'].Get('pTRecoBMatchSuccess').Divide(files['QCD2018.root'].Get('pTRecoB'))
files['QCD2018.root'].Get('pTRecoBMatchSuccess').Draw('HIST E PLC')
files['QCD2018.root'].Get('pTRecoBMatchSuccess').SetYTitle("Matched Events / Total Events")
c1.SaveAs('plots/QCDOverlaid/{}/pTRecoBMatchSuccess.png'.format(bTags))	

files['QCD2018.root'].Get('pTGenBMatchSuccess').Divide(files['QCD2018.root'].Get('pTGenBQCD'))
files['QCD2018.root'].Get('pTGenBMatchSuccess').Draw('HIST E PLC')
files['QCD2018.root'].Get('pTGenBMatchSuccess').SetYTitle("Matched Events / Total Events")
c1.SaveAs('plots/QCDOverlaid/{}/pTGenBMatchSuccess.png'.format(bTags))	
'''

if args.TH2:
	bb2DHists = list(set(bb2DHists))
	print(bb2DHists)
	for param in bb2DHists:
		for f in files:
			newfile = files[f]
			h_new = newfile.Get(param)
			h_new.SetTitle('{}'.format(f[:len(f) - 5]))
			if h_new.Integral() != 0: h_new.Scale(1 / h_new.Integral())
			h_new.SetStats(0)
			c1 = ROOT.TCanvas()
			c1.SetRightMargin(0.15)
			h_new.Draw('COLZ1')
			c1.SaveAs('plots/QCDOverlaid/{}/{}_{}.png'.format(bTags, param, f[:len(f) - 5]))
