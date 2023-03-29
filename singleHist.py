import ROOT
ROOT.gROOT.SetBatch(True)

file_path = 'root-test-files/QCD.root'
new_file = ROOT.TFile.Open(file_path, "READ")
h_new = new_file.Get('HT')
h_stack = ROOT.THStack()
h_stack.Add(h_new)
c1 = ROOT.TCanvas()

file_path = 'root-test-files/signal_1500_600.root'
new_file_2 = ROOT.TFile.Open(file_path, "READ")
h_new_2 = new_file_2.Get('HT')
#h_stack = ROOT.THStack()
h_stack.Add(h_new_2)
c1 = ROOT.TCanvas()
ROOT.gPad.SetLogy()
h_stack.SetMinimum(1)
h_stack.SetMaximum(10**7)
h_stack.Draw('HIST')

c1.SaveAs('plots/QCD_test.png')
