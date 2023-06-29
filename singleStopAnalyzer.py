#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
import glob
import argparse
from math import sqrt,fabs,copysign
from numpy import mean,std
from itertools import combinations
import pandas as pd
ROOT.PyConfig.IgnoreCommandLineOptions = True

def jetMatcher(genJets, recoJets, dR_match = 0.4, pT_match = 3):
  def getp4(obj):
  	return obj.p4()
  def dR(obj1, obj2):
  	return obj1.p4().DeltaR(obj2.p4())
  def pT(obj):
  	return obj.p4().Pt()
  matched = []
  remaining_idxs = set(range(len(genJets)))
  used_idxs = set()
  reverse_matched = {}
  for i, j in enumerate(recoJets):
  	diff = lambda x: dR(j, genJets[x])
  	if remaining_idxs: 
  		smallest = min(remaining_idxs, key = diff)
  		if diff(smallest) < dR_match and ( abs(pT(j) - pT(genJets[smallest])) / pT(genJets[smallest]) ) <  pT_match:
  			if smallest in used_idxs:
  				if dR(genJets[smallest], j) < dR(genJets[smallest], recoJets[reverse_matched[smallest]]):
  					matched.remove((reverse_matched[smallest], smallest))
  					matched.append((i, smallest))
  					reverse_matched[smallest] = i
  			used_idxs.add(smallest)
  			matched.append((i, smallest))
  			reverse_matched[smallest] = i
  	else: matched.append((None, None))
  return matched

class ExampleAnalysis(Module):

    def __init__(self,isData,isSignal,MCCampaign,isSkimmed,isCR0b):
        self.writeHistFile = True
        self.isSignal = isSignal
        self.MCCampaign = MCCampaign
        self.isSkimmed = isSkimmed
        self.isData = isData
        self.isCR0b = isCR0b

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

    def analyze(self, event):

        dRMatch = 0.1

        # Get MC weight
        if not self.isData: genWeight = 1 if event.genWeight > 0 else -1
        else: genWeight = 1

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

        if self.isSignal:
          genParts       = list(Collection(event,"GenPart"))
          genAK4Jets     = list(Collection(event,"GenJet"))
          # Get only outgoing particles of the hardest subprocess
          gens = filter(lambda x: (((x.statusFlags >> 13) & 1) and ((x.statusFlags >> 8) & 1)) and not (((abs(x.pdgId) == 1) or (abs(x.pdgId) == 3)) and ((x.statusFlags >> 11) & 1)), genParts)

        if not self.isSkimmed:
          goodElectrons  = filter(lambda x: x.cutBased == 4 and x.miniPFRelIso_all < 0.1 and x.pt > 30 and abs(x.eta) < 2.4,list(Collection(event,"Electron")))
          goodMuons      = filter(lambda x: x.mediumId and x.miniPFRelIso_all < 0.2 and x.pt > 30 and abs(x.eta) < 2.4,list(Collection(event,"Muon")))

        # Cuts
        if not self.isSkimmed:
          if not (event.HLT_PFHT1050 or event.HLT_AK8PFJet400_TrimMass30): return False
          if len(jets) > 0 and not jets[0].pt > 300: return False
          if len(jets) < 4 or len(jets) > 6: return False
          if len(goodElectrons) != 0 or len(goodMuons) != 0: return False
          if not 2 < abs(jets[0].p4().DeltaR(jets[1].p4())) < 4: return False
          if not self.isCR0b: 
            if not self.isData and len(looseBs) < 2: return False
            elif self.isData and len(looseBs) != 0: return False
            if not self.isData and abs(looseBs[0].p4().DeltaR(looseBs[1].p4())) < 1: return False
          elif len(looseBs) != 0: return False

        if self.isSignal:

          # Match gens to particle
          stopPlus = True
          if True in (g.pdgId == 1000006 for g in gens):
            for g in gens:
              if g.pdgId == 1000006: genStop = g
              elif g.pdgId == 1000024: genChi = g
              elif g.pdgId == 5: genBStop = g
              elif g.pdgId == -5: genBChi = g
              elif g.pdgId == -1: genD = g
              elif g.pdgId == -3: genS = g
              else: print('WARNING: Unexpected particle with pdgId {}'.format(g.pdgId))
            genStopPlus = genStop
            genBStopPlus = genBStop
            genChiPlus = genChi
            genBChiPlus = genBChi
          elif True in (g.pdgId == -1000006 for g in gens):
            stopPlus = False
            for g in gens:
              if g.pdgId == -1000006: genStop = g
              elif g.pdgId == -1000024: genChi = g
              elif g.pdgId == -5: genBStop = g
              elif g.pdgId == 5: genBChi = g
              elif g.pdgId == 1: genD = g
              elif g.pdgId == 3: genS = g
              else: print('WARNING: Unexpected particle with pdgId {}'.format(g.pdgId))
            genStopMinus = genStop
            genBStopMinus = genBStop
            genChiMinus = genChi
            genBChiMinus = genBChi
          else: print('WARNING: No stop found in event')
          genQuarks = [genBStop,genBChi,genD,genS]

        #-----------------------------------------------------------------------
        # RECO
        #-----------------------------------------------------------------------

        sumJet             = ROOT.TLorentzVector()
        sumJet4            = ROOT.TLorentzVector()
        sumJet3            = ROOT.TLorentzVector()
        sumJet3NoLead      = ROOT.TLorentzVector()
        sumJet3NoLeadOrSub = ROOT.TLorentzVector()

        # Loop over jets
        HT = 0; HT3 = 0; HT3NoLead = 0; nbLoose = 0; nbMedium = 0; nbTight=0
        for i,j in enumerate(jets):
          sumJet += j.p4()
          HT += j.pt
          if i < 4:           
            sumJet4 += j.p4()
          if i < 3:           
            sumJet3 += j.p4()
            HT3 += j.pt
          if i > 0 and i < 4: 
            sumJet3NoLead += j.p4()
            HT3NoLead += j.pt
          if i > 1 and i < 5: sumJet3NoLeadOrSub += j.p4()

        matches = jetMatcher(genQuarks,jets)
        if not len(matches) == 4: return False
        #for match in matches:
        #  print('Jet {} matches to pdgId {}, so it comes from {}'.format(match[0],genQuarks[match[1]].pdgId,'stop' if ((genQuarks[match[1]].pdgId == 5 and genD.pdgId < 0) or (genQuarks[match[1]].pdgId == -5 and genD.pdgId > 0)) else 'chargino'))
 
        with open('output/exportJetInfo/jets_1500_900.csv','a') as f:
          for i,j in enumerate(jets):
            isMatched = 1 if i in [x[0] for x in matches] else 0
            matchedQuark = genQuarks[[x[1] for x in matches if x[0] == i][0]] if isMatched else None
            isStopMatched, isChiMatched, isOther = -1,-1,-1
            if isMatched and ((matchedQuark.pdgId == 5 and genD.pdgId < 0) or (matchedQuark.pdgId == -5 and genD.pdgId > 0)): isStopMatched, isChiMatched, isOther = 1,0,0
            elif isMatched and ((genD.pdgId < 0 and matchedQuark.pdgId in [-1,-3,-5]) or (genD.pdgId > 0 and matchedQuark.pdgId in [1,3,5])): isStopMatched, isChiMatched, isOther = 0,1,0
            else: isStopMatched, isChiMatched, isOther = 0,0,1
            f.write('{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(
                     i,j.pt,j.eta,j.phi,j.btagDeepFlavB,
                     sumJet3NoLead.M(),sumJet3NoLead.Pt(),sumJet3NoLead.Eta(),sumJet3NoLead.Phi(),
                     sumJet4.M(),sumJet4.Pt(),sumJet4.Eta(),sumJet4.Phi(),
                     isStopMatched,isChiMatched,isOther
                    ))

          with open('output/exportJetInfo/baseline.txt','a') as f:
            for i,j in enumerate(jets):
              isMatched = 1 if i in [x[0] for x in matches] else 0
              matchedQuark = genQuarks[[x[1] for x in matches if x[0] == i][0]] if isMatched else None
              isStopMatched, isChiMatched, isOther = -1,-1,-1
              if isMatched and ((matchedQuark.pdgId == 5 and genD.pdgId < 0) or (matchedQuark.pdgId == -5 and genD.pdgId > 0)): isStopMatched, isChiMatched, isOther = 1,0,0
              elif isMatched and ((genD.pdgId < 0 and matchedQuark.pdgId in [-1,-3,-5]) or (genD.pdgId > 0 and matchedQuark.pdgId in [1,3,5])): isStopMatched, isChiMatched, isOther = 0,1,0
              else: isStopMatched, isChiMatched, isOther = 0,0,1
              matchIdx = isChiMatched + 2*isOther
              if i == 0:         f.write('{},0\n'.format(matchIdx))
              elif i in [1,2,3]: f.write('{},1\n'.format(matchIdx))
              else:              f.write('{},2\n'.format(matchIdx))

        return True

parser = argparse.ArgumentParser(description='Single Stop Analyzer')
parser.add_argument('--sample',type=str,default='signal',choices=['Data2018','signal','TT','TT2018','QCD','QCD2018','QCDInclusive2018','ZQQ2018','ST2018','WQQ2018','ZNuNu2018','Diboson2018'],help='Sample to run over')
parser.add_argument('--tag',type=str,default='test',help='Tag for output label')
parser.add_argument('-n',type=int,default=1,help='Sample index to run over for backgrounds')
parser.add_argument('--points',type=str,default='all',help='Signal point(s) to run over, comma separated in MSTOP_MCHI format; "all" to run over all available points')
parser.add_argument('--useskim',action='store_true',default=False,help='Flag to use NANOAODs skimmed with the nominal selections')
parser.add_argument('--CR0b',action='store_true',default=False,help='Flag to use the 0b CR selection')
args = parser.parse_args()

with open('output/exportJetInfo/jets_1500_900.csv','w') as f:
  f.write('jetOrdinality,jetPT,jetEta,jetPhi,jetBScore,m3M,m3PT,m3Eta,m3Phi,m4M,m4PT,m4Eta,m4Phi,isStopMatched,isChiMatched,isOther\n')

os.remove('output/exportJetInfo/baseline.txt')

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
		'(HLT_PFHT1050 || HLT_AK8PFJet400_TrimMass30)'
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
                    modules=[ExampleAnalysis(isData=1,
                             isSignal=0,
                             MCCampaign=MCCampaign,
                             isSkimmed=False,
                             isCR0b=args.CR0b)],
                    noOut=True, histFileName='{}/{}-{}.root'.format(outputPath,args.sample,args.n), histDirName="plots",
                    maxEntries=None)
  p.run()

elif args.sample == 'signal':

  allPoints = ['1000_400','1000_600','1000_900',
               '1200_400','1200_600','1200_1100',
               '1300_400','1300_600','1300_1200',
               '1400_400','1400_600','1400_1300',
               '1500_400','1500_600','1500_900','1500_1400',
               '2000_400','2000_600','2000_900','2000_1400','2000_1900']
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
    files = glob.glob('/eos/uscms/store/user/dmahon/condor/RPVSingleStopMC/NANOAOD-ALL/NANOAOD-{}.root'.format(masses))
    #files = glob.glob('/eos/uscms/store/user/dmahon/condor/RPVSingleStopMC/NANOAOD/NANOAOD-{}-*.root'.format(masses))
    files = ['root://cmsxrootd.fnal.gov/' + x.replace('/eos/uscms','') for x in files]
    #files = ['file:/uscms_data/d3/dmahon/RPVSingleStopRun3Patched/NANOAOD/CMSSW_12_4_5/test_2000_100-1.root']
    #files = ['/uscms_data/d3/dmahon/RPVSingleStopRun3Patched/NANOAOD/files/NANOAOD-{}.root'.format(masses)]
    p = PostProcessor(".", files, cut=preselection, branchsel=None,
                      modules=[ExampleAnalysis(isData=0,
                                               isSignal=1,
                                               MCCampaign='UL2018',
                                               isSkimmed=False,
                                               isCR0b=args.CR0b)],
                      noOut=True, histFileName='{}/{}_{}.root'.format(outputPath,args.sample,masses), histDirName="plots",
                      maxEntries=None)
    p.run()

elif args.useskim: 

  print('Running over skimmed {} files'.format(args.sample))
  print('Using UL 2018 MC campaign working points')

  files = ['root://cmsxrootd.fnal.gov//store/user/ckapsiak/SingleStop/Skims/Skim_2023_23_03/{}.root'.format(args.sample)]
  if len(files) != 1: print('WARNING: Multiple files selected. All must be from the same MC campaign.')
  p = PostProcessor(".", files, cut='', branchsel=None,
                    modules=[ExampleAnalysis(isData=0,
                                             isSignal=0,
                                             MCCampaign='UL2018',
                                             isSkimmed=True,
                                             isCR0b=args.CR0b)],
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
                    modules=[ExampleAnalysis(isData=0,
                                             isSignal=0,
                                             MCCampaign=MCCampaign,
                                             isSkimmed=False,
                                             isCR0b=args.CR0b)], 
                    noOut=True, histFileName='{}/{}-{}.root'.format(outputPath,args.sample,args.n), histDirName="plots",
                    maxEntries=None)
  p.run() 
