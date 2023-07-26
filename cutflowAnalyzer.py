#!/usr/bin/env python
from framework.eventloop import Module
from framework.datamodel import Collection
from framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
import glob
import argparse
import array
from math import sqrt,fabs,copysign
from numpy import mean,std
from itertools import combinations
sys.path.append('~/nobackup/SingleStop/CMSSW_10_6_19_patch2/src/PhysicsTools/NanoAODTools/python/postprocessing/SingleStop/')
from jetMatching import jetMatcher, checker, dR, pT
ROOT.PyConfig.IgnoreCommandLineOptions = True

class ExampleAnalysis(Module):
    def __init__(self,isSignal,MCCampaign,isSkimmed,coupling,isQCD,bAlgo,isData,isCR0b):
			self.writeHistFile = True
			self.isSignal = isSignal
			self.MCCampaign = MCCampaign
			self.isSkimmed = isSkimmed
			self.coupling = coupling
			self.isQCD = isQCD
			self.bAlgo = bAlgo
			self.isData = isData
			self.isCR0b = isCR0b

    def beginJob(self, histFile=None, histDirName=None):
			Module.beginJob(self, histFile, histDirName)
			self.h_cutflow312                  = ROOT.TH1D('cutflow312',    	';Cut',      				10,     0,      10      )
			self.h_cutflow313                  = ROOT.TH1D('cutflow313',    	';Cut',      				10,     0,      10      )

			for h in list(vars(self)):
				if h[0:2] == 'h_': self.addObject(getattr(self, h))
			
    def analyze(self, event):

			dRMatch = 0.1

			# Get MC weight
			if not self.isData: genWeight = 1 if event.genWeight > 0 else -1
			else: genWeight = 1

			# Fill histograms before selections (but after any pre-selections)

			# Set b tagging WPs
			if   self.MCCampaign == 'UL2016preVFP': 	bTagWPs = [0.0508,0.2598,0.6502]
			elif self.MCCampaign == 'UL2016postVFP': 	bTagWPs = [0.0480,0.2489,0.6377]
			elif self.MCCampaign == 'UL2017':     		bTagWPs = [0.0532,0.3040,0.7476]
			elif self.MCCampaign == 'UL2018':     		bTagWPs = [0.0490,0.2783,0.7100]

			# Get event collections
			jets           = filter(lambda x: x.pt > 30 and abs(x.eta) < 2.4,list(Collection(event,"Jet")))
			fatJets        = list(Collection(event,"FatJet"))
			looseTs        = filter(lambda x: x.particleNet_TvsQCD > 0.58,fatJets)
			mediumTs       = filter(lambda x: x.particleNet_TvsQCD > 0.80,fatJets)
			tightTs        = filter(lambda x: x.particleNet_TvsQCD > 0.97,fatJets)
			looseWs        = filter(lambda x: x.particleNet_WvsQCD > 0.70,fatJets)
			mediumWs       = filter(lambda x: x.particleNet_WvsQCD > 0.94,fatJets)
			tightWs        = filter(lambda x: x.particleNet_WvsQCD > 0.98,fatJets)
			deepWP1Ts      = filter(lambda x: x.deepTag_TvsQCD > 0.436,fatJets)
			deepWP2Ts      = filter(lambda x: x.deepTag_TvsQCD > 0.802,fatJets)
			deepWP3Ts      = filter(lambda x: x.deepTag_TvsQCD > 0.922,fatJets)
			deepWP4Ts      = filter(lambda x: x.deepTag_TvsQCD > 0.989,fatJets)
			deepWP1Ws      = filter(lambda x: x.deepTag_WvsQCD > 0.458,fatJets)
			deepWP2Ws      = filter(lambda x: x.deepTag_WvsQCD > 0.762,fatJets)
			deepWP3Ws      = filter(lambda x: x.deepTag_WvsQCD > 0.918,fatJets)
			deepWP4Ws      = filter(lambda x: x.deepTag_WvsQCD > 0.961,fatJets)
			looseBs        = filter(lambda x: x.btagDeepFlavB > bTagWPs[0],jets)
			mediumBs       = filter(lambda x: x.btagDeepFlavB > bTagWPs[1],jets)
			tightBs        = filter(lambda x: x.btagDeepFlavB > bTagWPs[2],jets)
			MET            = ROOT.TVector2()
			MET.SetMagPhi(event.MET_pt,event.MET_phi)

			if not self.isData:
				genParts       = list(Collection(event,"GenPart"))
				genAK4Jets     = filter(lambda x: x.pt > 30 and abs(x.eta) < 2.4, list(Collection(event,"GenJet")))
			if self.isSignal:
				# Get only outgoing particles of the hardest subprocess
				gens = filter(lambda x: (((x.statusFlags >> 13) & 1) and ((x.statusFlags >> 8) & 1)) and not (((abs(x.pdgId) == 1) or (abs(x.pdgId) == 5)) and ((x.statusFlags >> 11) & 1)), genParts)
			elif self.isQCD:
				# 7: FromHardProcess, 13: IsLastCopy (gen Bs)
				gens = filter(lambda x: (x.statusFlags >> 13) & 1, genParts)
				genBJets = filter(lambda x: abs(x.partonFlavour) == 5, genAK4Jets)
			if not self.isSkimmed:
				goodElectrons  = filter(lambda x: x.cutBased == 4 and x.miniPFRelIso_all < 0.1 and x.pt > 30 and abs(x.eta) < 2.4,list(Collection(event,"Electron")))
				goodMuons      = filter(lambda x: x.mediumId and x.miniPFRelIso_all < 0.2 and x.pt > 30 and abs(x.eta) < 2.4,list(Collection(event,"Muon")))

			# Cuts
			if not self.isSkimmed:
				self.h_cutflow312.Fill(0,genWeight)
				self.h_cutflow313.Fill(0,genWeight)

				if not (event.HLT_PFHT1050 or event.HLT_AK8PFJet400_TrimMass30): return False
				self.h_cutflow312.Fill(1,genWeight)
				self.h_cutflow313.Fill(1,genWeight)

				if len(jets) > 0 and jets[0].pt < 300: return False
				self.h_cutflow312.Fill(2,genWeight)
				self.h_cutflow313.Fill(2,genWeight)

				if len(jets) < 4 or len(jets) > 6: return False
				self.h_cutflow312.Fill(3,genWeight)
				self.h_cutflow313.Fill(3,genWeight)

				if len(goodElectrons) != 0 or len(goodMuons) != 0: return False
				self.h_cutflow312.Fill(4,genWeight)
				self.h_cutflow313.Fill(4,genWeight)

				if len(jets) > 1 and not 2 < abs(jets[0].p4().DeltaR(jets[1].p4())) < 4: return False
				self.h_cutflow312.Fill(5,genWeight)
				self.h_cutflow313.Fill(5,genWeight)

				if len(tightBs) < 3:
					if (len(mediumBs) < 2 or len(tightBs) < 1): return False
					else: 
						self.h_cutflow312.Fill(6, genWeight)
						if abs(mediumBs[0].p4().DeltaR(mediumBs[1].p4())) >= 1: self.h_cutflow312.Fill(7, genWeight)
				else:
					self.h_cutflow312.Fill(6, genWeight)
					self.h_cutflow313.Fill(6, genWeight)
					if abs(mediumBs[0].p4().DeltaR(mediumBs[1].p4())) >= 1: self.h_cutflow312.Fill(7, genWeight)
					if abs(tightBs[0].p4().DeltaR(tightBs[1].p4())) >= 1: self.h_cutflow313.Fill(7, genWeight)
			return True

parser = argparse.ArgumentParser(description='Single Stop Analyzer')
parser.add_argument('--sample',type=str,default='signal',choices=['Data2018','signal','TT','TT2018','QCD','QCD2018','QCDInclusive2018','ZQQ2018','ST2018','WQQ2018','ZNuNu2018','Diboson2018'],help='Sample to run over')
parser.add_argument('--tag',type=str,default='test',help='Tag for output label')
parser.add_argument('-n',type=int,default=1,help='Sample index to run over for backgrounds')
parser.add_argument('--points',type=str,default='all',help='Signal point(s) to run over, comma separated in MSTOP_MCHI format; "all" to run over all available points')
parser.add_argument('--useskim',action='store_true',default=False,help='Flag to use NANOAODs skimmed with the nominal selections')
parser.add_argument('--coupling', type = str, default = '313', choices = ['312', '313'])
parser.add_argument('--bAlgo', type = str, default = 'medium', choices = ['tight', 'medium'])
parser.add_argument('--CR0b',action='store_true',default=False,help='Flag to use the 0b CR selection')
args = parser.parse_args()
outputPath = 'output/{}'.format(args.tag)
if not os.path.exists(outputPath):
  os.makedirs(outputPath)

if   args.sample == 'TT':               sampleFile = 'TTToHadronic.txt'
elif args.sample == 'TT2018':           sampleFile = 'TTToHadronic2018.txt'
elif args.sample == 'QCD':              sampleFile = 'QCDBEnriched.txt'
elif args.sample == 'QCD2018':          sampleFile = 'QCDBEnriched2018.txt'
elif args.sample == 'QCDInclusive2018': sampleFile = 'QCDInclusive2018.txt'
elif args.sample == 'ZQQ2018':          sampleFile = 'ZJetsToQQ2018.txt'
elif args.sample == 'ST2018':           sampleFile = 'STHadronic2018.txt'
elif args.sample == 'WQQ2018':          sampleFile = 'WJetsToQQ2018.txt'
elif args.sample == 'ZNuNu2018':        sampleFile = 'ZJetsToNuNu2018.txt'
elif args.sample == 'Diboson2018':      sampleFile = 'Diboson2018.txt'
elif args.sample == 'Data2018':         sampleFile = 'Data2018.txt'
elif args.sample != 'signal': print('ERROR: Unexpected sample argument')

preselection = (
		'(Jet_pt[3] > 30) &&'
		'(Jet_pt[0] > 300) &&'
		'HLT_PFHT1050 || HLT_AK8PFJet400_TrimMass30'
               )

if args.sample == 'Data2018':

  print('Running over {} files'.format(args.sample))
  print('Blinding data, only using CR cuts (0 loose b\'s)')
  print('Using UL 2018 campaign working points')

  files = open('samples/{}'.format(sampleFile)).read().split('\n')
  files = [['root://cmsxrootd.fnal.gov/' + x.replace('/eos/uscms','') for x in files][:-1][args.n - 1]]
  if len(files) != 1: print('WARNING: Multiple files selected. All must be from the same MC campaign.')
  if 'UL2016' in files[0]:
    if 'preVFP' in files[0]:    MCCampaign = 'UL2016preVFP'
    else:                       MCCampaign = 'UL2016postVFP'
  elif 'UL2017' in files[0]:      MCCampaign = 'UL2017'
  elif 'UL2018' in files[0]:      MCCampaign = 'UL2018'
  else:
    print('ERROR: Unable to determine campaign of {}'.format(files[0]))
    sys.exit()
  p = PostProcessor(".", files, cut=preselection, branchsel=None,
                    modules=[ExampleAnalysis(isData=1,isSignal=0,MCCampaign=MCCampaign,isSkimmed=False, coupling = args.coupling, bAlgo = args.bAlgo, isQCD = 0, isCR0b = args.CR0b)],
                    noOut=True, histFileName='{}/{}-{}.root'.format(outputPath,args.sample,args.n), histDirName="plots",
                    maxEntries=None)
  p.run()

elif args.sample == 'signal':

  allPoints = ['1200_400', '1200_1100', '1500_600', '1500_900', '1500_1400', '2000_1900'] if args.coupling == '312' else ['1000_400', '1000_900', '1500_600', '1500_900', '1500_1400', '2000_1900']
							 #,'1000_600','1000_900',
               #'1200_400','1200_600','1200_1100',
               #'1300_400','1300_600','1300_1200',
               #'1400_400','1400_600','1400_1300',
               #'1500_400','1500_600','1500_900','1500_1400',
               #'2000_400','2000_600','2000_900','2000_1400','2000_1900']
  if args.points == 'all': 
    print('Running over all available signal points...')
    points = allPoints
  elif '_' in args.points: 
    points = args.points.split(',')
    if False in [x in allPoints for x in points]: 
      print('ERROR: Invalid --points argument provided. Available signal points: {}'.format(allPoints))
      sys.exit()
  else: 
    print('ERROR: Invalid --points arguement format provided')
    sys.exit()
  for masses in points:
    if args.coupling == '313':
    	files = glob.glob('/eos/uscms/store/user/dmahon/condor/RPVSingleStopMC313/NANOAOD-ALL/NANOAOD-{}.root'.format(masses))
    else:
			files = glob.glob('/eos/uscms/store/user/dmahon/condor/RPVSingleStopMC/NANOAOD-ALL/NANOAOD-{}.root'.format(masses))
			#files = glob.glob('/eos/uscms/store/user/dmahon/condor/RPVSingleStopMC/NANOAOD/NANOAOD-{}-*.root'.format(masses))
			files = ['root://cmsxrootd.fnal.gov/' + x.replace('/eos/uscms','') for x in files]
			#files = ['file:/uscms_data/d3/dmahon/RPVSingleStopRun3Patched/NANOAOD/CMSSW_12_4_5/test_2000_100-1.root']
			#files = ['/uscms_data/d3/dmahon/RPVSingleStopRun3Patched/NANOAOD/files/NANOAOD-{}.root'.format(masses)]
    p = PostProcessor(".", files, cut=preselection, branchsel=None,
											modules=[ExampleAnalysis(isCR0b = args.CR0b,isSignal=1,MCCampaign='UL2018',isSkimmed=False, coupling = args.coupling, isQCD = 0, bAlgo = args.bAlgo, isData = 0)],
											noOut=True, histFileName='{}/{}_{}.root'.format(outputPath,args.sample,masses), histDirName="plots",
											maxEntries=None)
    p.run()

elif args.useskim: 

  print('Running over skimmed {} files'.format(args.sample))
  print('Using UL 2018 MC campaign working points')

  files = ['root://cmsxrootd.fnal.gov//store/user/ckapsiak/SingleStop/Skims/Skim_2023_23_03/{}.root'.format(args.sample)]
  if len(files) != 1: print('WARNING: Multiple files selected. All must be from the same MC campaign.')
  p = PostProcessor(".", files, cut='', branchsel=None,
                    modules=[ExampleAnalysis(isCR0b = args.CR0b, isSignal=0,MCCampaign='UL2018',isSkimmed=True, coupling = args.coupling, bAlgo = args.bAlgo, isData = 0)],
                    noOut=True, histFileName='{}/{}-ALL.root'.format(outputPath,args.sample), histDirName="plots",
                    maxEntries=None)
  p.run()

else:

  print('Running over file index {} for sample {}'.format(args.n,args.sample))

  files = open('samples/{}'.format(sampleFile)).read().split('\n')
  files = [['root://cmsxrootd.fnal.gov/' + x.replace('/eos/uscms','') for x in files][:-1][args.n - 1]]
  if len(files) != 1: print('WARNING: Multiple files selected. All must be from the same MC campaign.')
  if 'UL16' in files[0]: 
    if 'preVFP' in files[0]:	MCCampaign = 'UL2016preVFP'
    else: 			MCCampaign = 'UL2016postVFP'
  elif 'UL17' in files[0]: 	MCCampaign = 'UL2017'
  elif 'UL18' in files[0]:	MCCampaign = 'UL2018'
  else: 
    print('ERROR: Unable to determine MC campaign of {}'.format(files[0]))
    sys.exit()
  p = PostProcessor(".", files, cut=preselection, branchsel=None, 
                    modules=[ExampleAnalysis(isSignal=0,MCCampaign=MCCampaign,isSkimmed=False, coupling = args.coupling, 
										isQCD = (args.sample == 'QCD2018' or args.sample == 'QCDInclusive2018'), bAlgo = args.bAlgo, isData = 0, isCR0b = args.CR0b)],
										noOut=True, histFileName= '{}/{}-{}.root'.format(outputPath,args.sample,args.n), histDirName = "plots",
										maxEntries = None)
  p.run() 
