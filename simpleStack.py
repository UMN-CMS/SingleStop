import ROOT
ROOT.gROOT.SetBatch(True)
path = 'root-test-files'
files = ['QCD.root', 'TT.root', 'signal_1500_600.root']
open_files = {'QCD.root': ROOT.TFile.Open('{}/QCD.root'.format(path), "READ"), 'TT.root': ROOT.TFile.Open('{}/TT.root'.format(path), "READ"), 'signal_1500_600.root': ROOT.TFile.Open('{}/signal_1500_600.root'.format(path), "READ")}
colors = [2, 3, 4]
h_stack = ROOT.THStack()
c1 = ROOT.TCanvas()
ROOT.gStyle.SetPalette(ROOT.kOcean)
for i, f in enumerate(files):
	file_path = '{}/{}'.format(path, f)
	new_file = open_files[f]
	h_pT1Gen = new_file.Get('pT1')
	#h_pT1Gen.Draw('HIST')
	#h_pT1Gen.SetLineColor(colors[i])
	h_stack.Add(h_pT1Gen)
h_stack.Draw('HIST')
ROOT.gPad.SetLogy()
c1.SaveAs('plots/pT1Gen-stack.png')
