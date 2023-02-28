import ROOT

path = 'root-test-files'
files = ['QCD.root', 'TT.root', 'signal_1500_600.root']
colors = [2, 3, 4]
h_stack = ROOT.THStack()

for i, f in enumerate(files):
	file_path = '{}/{}'.format(path, f)
	new_file = ROOT.TFile.Open(file_path, "READ")
	h_pT1Gen = new_file.Get('pT1Gen')
	h_stack.Add(h_pT1Gen)

c1 = ROOT.TCanvas()
h_stack.Draw('HIST PFC')
c1.SaveAs('plots/simple-test-stack.png')
