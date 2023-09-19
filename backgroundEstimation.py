import ROOT, array
import matplotlib.pyplot as plt
ROOT.gROOT.SetBatch(True)

samples = ['TT2018', 'QCDInc2018', 'ZQQ2018', 'ST2018', 'WQQ2018', 'ZNuNu2018', 'Diboson2018']
h_values = ROOT.THStack()
legend = ROOT.TLegend(0.812, 0.57, 0.995, 0.85)
ROOT.gStyle.SetPalette(ROOT.kRainBow)
ROOT.gStyle.SetPaintTextFormat("5.2e")

class Sample:
	def __init__(self, name, count):
		self.name = name
		self.count = count

values312 = {}
values313 = {}
totalCount312 = 0
totalCount313 = 0

for sample in samples:
	if sample != 'QCDInc2018':
		f = ROOT.TFile.Open('{0}_Cutflows/output/{0}.root'.format(sample))
		h_312 = f.Get('cutflow312')
		values312[sample] = Sample(sample, h_312.GetBinContent(8))
		h_313 = f.Get('cutflow313')
		values313[sample] = Sample(sample, h_313.GetBinContent(8))
		totalCount312 += values312[sample].count
		totalCount313 += values313[sample].count
	else:
		f = ROOT.TFile.Open('{}_Cutflows/output/QCDInclusive2018.root'.format(sample))
		h_312 = f.Get('cutflow312')
		values312[sample] = Sample(sample, h_312.GetBinContent(8))
		h_313 = f.Get('cutflow313')
		values313[sample] = Sample(sample, h_313.GetBinContent(8))
		totalCount312 += values312[sample].count
		totalCount313 += values313[sample].count

sampleList = list(values312.values())
print([sample.name for sample in sampleList])
sampleList = sorted(sampleList, key = lambda sample: sample.count)
print([sample.name for sample in sampleList])

labels = {'TT2018': 'tt', 'QCDInc2018': 'Inclusive QCD', 'ZQQ2018': 'Z to qq', 'ST2018': 'Single top', 'WQQ2018': 'W to qq', 'ZNuNu2018': 'Z to nu nu', 'Diboson2018': 'Diboson'}

for s in sampleList:
	h_new = ROOT.TH1F('', '', 2, 0, 2)
	h_new.Fill(0, s.count / totalCount312)
	print(s.name, 100 * s.count / totalCount312)
	h_new.Fill(1, values313[s.name].count / totalCount313)
	print(s.name, 100 * values313[s.name].count / totalCount313)
	h_new.SetLineWidth(3)
	h_new.SetMarkerColor(ROOT.kWhite)
	h_values.Add(h_new)
	legend.AddEntry(h_new, labels[s.name])

c1 = ROOT.TCanvas()
c1.SetRightMargin(0.2)

h_values.Draw('HIST PLC PFC')
h_values.GetXaxis().SetBinLabel(1, '')
h_values.GetXaxis().SetBinLabel(2, '')
h_values.GetYaxis().SetTitle('Proportion of Total Background')

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.DrawLatex(0.14, 0.06, "\lambda_{312}^{\prime \prime} \ \mathrm{Signal \ Region}")
latex.DrawLatex(0.5, 0.06, "\lambda_{313}^{\prime \prime} \ \mathrm{Signal \ Region}")

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.035)
latex.SetTextColor(ROOT.kWhite)
latex.DrawLatex(0.22, 0.75, "0.9236")
latex.DrawLatex(0.22, 0.63, "0.0487")
latex.DrawLatex(0.22, 0.57, "0.0124")
latex.DrawLatex(0.22, 0.51, "0.00933")
latex.DrawLatex(0.22, 0.39, "0.00643")
latex.DrawLatex(0.22, 0.27, "0.00047")
latex.DrawLatex(0.22, 0.15, "0.00022")
latex.DrawLatex(0.58, 0.74, "0.9335")
latex.DrawLatex(0.58, 0.62, "0.0368")
latex.DrawLatex(0.58, 0.58, "0.00785")
latex.DrawLatex(0.58, 0.50, "0.0180")
latex.DrawLatex(0.58, 0.37, "0.00304")
latex.DrawLatex(0.58, 0.28, "0.00042")
latex.DrawLatex(0.58, 0.16, "0.00033")

legend.Draw()

ROOT.gPad.SetLogy()
h_values.SetMinimum(10 ** -4)

c1.Draw()
c1.SaveAs('plots/CutflowPieCharts/backgroundStack.pdf')	

