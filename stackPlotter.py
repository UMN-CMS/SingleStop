from collections import OrderedDict
import ROOT, argparse, os, math
ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description = 'Stack Plotter')
parser.add_argument('--bTags', type = str, default = '300')
parser.add_argument('--coupling', type = str, default = '313')
parser.add_argument('--TH1', default = True)
parser.add_argument('--TH2', default = False)

args = parser.parse_args()
bTags = args.bTags
coupling = args.coupling
TH1 = args.TH1
TH2 = args.TH2

signalPath = 'output/{}/{}/scaled'.format(coupling, bTags)
signalFiles = ['signal_2000_1900.root', 'signal_1500_900.root', 'signal_1200_400.root'] if coupling == '312' else ['signal_2000_1900.root', 'signal_1500_900.root', 'signal_1000_400.root']
colors = [2, 3, 4]
QCDPath = 'QCD2018_{}/output'.format(bTags)
QCDFile = 'QCD2018.root'
TTPath = 'TT2018_{}/output'.format(bTags)
TTFile = 'TT2018.root'
DataPath = 'Data2018_{}/output'.format(bTags)
DataFile = 'Data2018.root'

if not os.path.exists('plots/QCDOverlaid/{}/{}'.format(bTags, coupling)): os.makedirs('plots/QCDOverlaid/{}/{}'.format(bTags, coupling))

ROOT.gStyle.SetPalette(ROOT.kCMYK)

files = OrderedDict()
files = {fname: ROOT.TFile.Open('{}/{}'.format(signalPath, fname), "READ") for fname in signalFiles}
#files.update( { QCDFile: ROOT.TFile.Open('{}/{}'.format(QCDPath, QCDFile), "READ") } )
#files.update( { TTFile: ROOT.TFile.Open('{}/{}'.format(TTPath, TTFile), "READ") } )
files.update( { 'QCDInclusive2018.root': ROOT.TFile.Open('QCDInc2018_{}/output/QCDInclusive2018.root'.format(bTags), "READ") } )
#files.update( { 'Data2018.root': ROOT.TFile.Open('Data2018_{}/output/Data2018.root'.format(bTags), "READ") } )

'''
HTBins = [50, 100, 200, 300, 500, 700, 1000, 1500, 2000]
for HT in HTBins:
	files.update( { 'QCDInclusive2018-HT{}.root'.format(HT): ROOT.TFile.Open('QCD2018Inc_{}/output/QCDInclusive2018-HT{}.root'.format(bTags, HT), "READ") } )
'''
param_list = [key.GetName() for key in files['QCDInclusive2018.root'].GetListOfKeys() if ('dEtabb' in key.GetName() or 'dPhibb' in key.GetName() or 'dRbb' in key.GetName())]
#param_list = ['m3', 'm3NoLead', 'm4']
#bTagTitles = {'0B': '0 bs', '1B': '1 b', '020': '2 medium bs', '003': '3 tight bs'}

bb2DHists = []
for param in param_list:
	h_stack = ROOT.THStack(param, "{};{};Events".format(param, files['QCDInclusive2018.root'].Get(param).GetXaxis().GetTitle()))
	legend = ROOT.TLegend(0.68, 0.89, 0.89, 0.65)
	c1 = ROOT.TCanvas()
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
		elif param != 'cutflow' and h_new.Integral() != 0: h_new.Scale(1 / h_new.Integral())
		h_new.SetLineWidth(3)
		h_stack.Add(h_new)
		legend.AddEntry(h_new, f[:-5])
	#if h_stack.GetHists() == 0 or h_new.InheritsFrom('TH2D') or h_new.GetEntries() == 0: continue
	if args.TH1 and 'Vs' not in param:		
		h_stack.Draw('HIST PLC NOSTACK E')
		if param == 'dRbb12':
			h_stack.SetTitle(r'#Delta R between the leading two bs')
			if '312' in signalPath:	h_stack.GetXaxis().SetTitle(r'#Delta R_{b_{1}, b_{2}} (medium)')	
			else:	h_stack.GetXaxis().SetTitle(r'#Delta R_{b_{1}, b_{2}} (tight)')
		#h_stack.Draw('SAME E')
		if param == 'cutflow':
			h_stack.SetMinimum(0)
			h_stack.SetMaximum(1)
		legend.Draw()
		'''
		if param != 'cutflow':
			ROOT.gPad.SetLogy()
			h_stack.SetMaximum(10**7)
			h_stack.SetMinimum(10**0)
		'''
		c1.SaveAs('plots/QCDOverlaid/{}/{}/{}.png'.format(bTags, coupling, param))

'''
for param in ['m3', 'm4', 'm3NoLead']:
	for f in signalFiles:
		legend = ROOT.TLegend(0.68, 0.89, 0.89, 0.75)
		h_signal = files[f].Get(param)
		h_background = files['QCDInclusive2018.root'].Get(param)
		for i in range(h_background.GetNcells()):
			if h_background.GetBinContent(i) != 0: 
				h_signal.SetBinContent(i, h_signal.GetBinContent(i) / math.sqrt(h_background.GetBinContent(i)))
				#if '1200_400' in f: print(h_signal.GetBinContent(i), h_background.GetBinContent(i))
			else: h_signal.SetBinContent(i, 0)
		c1 = ROOT.TCanvas()
		h_signal.SetYTitle('s / root(b)')
		h_signal.SetStats(0)
		h_signal.Draw('HIST')
		h_signal.Rebin(4)
		legend.AddEntry(h_signal, f)
		legend.Draw()
		c1.SaveAs('plots/QCDOverlaid/{}/{}/{}Significance_{}.png'.format(bTags, coupling, param, f[7:-5]))	
'''

'''
for param in ['HT', 'pTFat1', 'mSoftFat1']:
	h_data = files['Data2018.root'].Get(param)
	h_background = files['QCDInclusive2018.root'].Get(param)
	for i in range(h_background.GetNcells()):
		if h_background.GetBinContent(i) != 0: h_data.SetBinContent(i, h_data.GetBinContent(i) / h_background.GetBinContent(i))
		else: h_data.SetBinContent(i, 0)
	c1 = ROOT.TCanvas()
	h_data.SetYTitle('Data / MC')
	h_data.SetTitle('{} (Only HT Trigger)'.format(param))
	h_data.SetStats(0)
	h_data.SetLineWidth(2)
	h_data.Draw('HIST E')
	c1.SaveAs('plots/QCDOverlaid/{}/{}/{}_DataOverMCRatio.png'.format(bTags, coupling, param))	
'''

'''
for param in param_list:
	legend = ROOT.TLegend(0.68, 0.89, 0.89, 0.65)

	h_background = files['QCDInclusive2018.root'].Get(param)
	h_data = files['Data2018.root'].Get(param)
	if h_background != None: h_background.SetLineWidth(3)
	if h_data != None: h_data.SetLineWidth(3)	
	if h_data == None or h_data.InheritsFrom('TH2F') or h_data.InheritsFrom('TH2D') or h_data.GetEntries() == 0: continue	
	if param == 'cutflow': 
		h_data.Scale(1 / h_data.GetMaximum())
		h_background.Scale(1 / h_background.GetMaximum())
		h_data.Fill(7, h_data.GetBinContent(7))
		h_data.SetBinError(8, h_data.GetBinError(7))		

	h_stack = ROOT.THStack(param, "{} (Data vs. MC);{};Events".format(param, files['Data2018.root'].Get(param).GetXaxis().GetTitle()))
	h_stack.Add(h_background)
	h_stack.Add(h_data)	
	legend.AddEntry(h_background, 'MC (QCD Inclusive 2018)')
	legend.AddEntry(h_data, 'Data2018')
	

	c1 = ROOT.TCanvas()
	h_stack.Draw('HIST E NOSTACK PLC')
	if param != 'cutflow':
		ROOT.gPad.SetLogy()
		h_stack.SetMaximum(1.05 * h_background.GetMaximum())
		h_stack.SetMinimum(10**0)
	legend.Draw()
	if param == 'cutflow':
		h_stack.GetXaxis().SetBinLabel(1, "Triggers")
		h_stack.GetXaxis().SetBinLabel(2, "pT1 > 300")
		h_stack.GetXaxis().SetBinLabel(3, "4 < nJets < 6")
		h_stack.GetXaxis().SetBinLabel(4, "No electrons or muons")
		h_stack.GetXaxis().SetBinLabel(5, "2 < dR12 < 4")
		h_stack.GetXaxis().SetBinLabel(6, "0 loose bs".format(bTags))
		h_stack.GetXaxis().SetBinLabel(7, "dRbb12 < 1")
	c1.SaveAs('plots/QCDOverlaid/{}/{}/{}_DataVsMC.png'.format(bTags, coupling, param))
'''

'''
c1 = ROOT.TCanvas()
files['QCD2018.root'].Get('pTRecoBMatchSuccess').Divide(files['QCD2018.root'].Get('pTRecoB'))
files['QCD2018.root'].Get('pTRecoBMatchSuccess').Draw('HIST E PLC')
files['QCD2018.root'].Get('pTRecoBMatchSuccess').SetYTitle("Matched Events / Total Events")
c1.SaveAs('plots/QCDOverlaid/{}/pTRecoBMatchSuccess.png'.format(bTags))	

files['QCD2018.root'].Get('pTGenBMatchSuccess').Divide(files['QCD2018.root'].Get('pTGenB'))
files['QCD2018.root'].Get('pTGenBMatchSuccess').Draw('HIST E PLC')
files['QCD2018.root'].Get('pTGenBMatchSuccess').SetYTitle("Matched Events / Total Events")
c1.SaveAs('plots/QCDOverlaid/{}/pTGenBMatchSuccess.png'.format(bTags))	
'''

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
'''
