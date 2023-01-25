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
ROOT.PyConfig.IgnoreCommandLineOptions = True

class ExampleAnalysis(Module):

    def __init__(self,isSignal):
        self.writeHistFile = True
        self.isSignal = isSignal

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        #-----------------------------------------------------------------------
        # GEN
        #-----------------------------------------------------------------------

        # Jet kinematics
        self.h_pT1Gen          		= ROOT.TH1F('pT1Gen',		';p_{T,1}^{gen} [GeV]',			150,	0, 	1500 	)
        self.h_pT2Gen          		= ROOT.TH1F('pT2Gen',		';p_{T,2}^{gen} [GeV]',			150,	0,     	1500 	)
        self.h_pT3Gen            	= ROOT.TH1F('pT3GEn',		';p_{T,3}^{gen} [GeV]',			150,	0,     	1500 	)
        self.h_pT4Gen            	= ROOT.TH1F('pT4Gen',		';p_{T,4}^{gen} [GeV]',			150,	0,     	1500 	)
        self.h_eta1Gen     	       	= ROOT.TH1F('eta1Gen',		';#eta_{1}^{gen}',			80,	-8,	8	)
        self.h_eta2Gen  	       	= ROOT.TH1F('eta2Gen',		';#eta_{2}^{gen}',			80,	-8,    	8	)
        self.h_eta3Gen		       	= ROOT.TH1F('eta3Gen',		';#eta_{3}^{gen}',			80,	-8,   	8	)
        self.h_eta4Gen       		= ROOT.TH1F('eta4Gen',		';#eta_{4}^{gen}',			80,	-8,   	8	)

        # SUSY particle kinematics
        self.h_pTStop              	= ROOT.TH1F('pTStop', 		';p_{T,{#tilde{t}}  [GeV]',		150,	0,    	1500	)
        self.h_pTStopPlus          	= ROOT.TH1F('pTStopPlus', 	';p_{T,{#tilde{t}^{+2/3}}  [GeV]',	150,	0,    	1500	)
        self.h_pTStopMinus         	= ROOT.TH1F('pTStopMinus', 	';p_{T,#tilde{t}^{-2/3}} [GeV]',	150,	0,    	1500	)
        self.h_pTChi               	= ROOT.TH1F('pTChi',		';p_{T,#chi^{#pm}}  [GeV]',		150,	0,    	1500	)
        self.h_pTChiPlus           	= ROOT.TH1F('pTChiPlus', 	';p_{T,#chi^{+}} [GeV]',		150,	0,    	1500	)
        self.h_pTChiMinus          	= ROOT.TH1F('pTChiMinus', 	';p_{T,#chi^{-}} [GeV]',		150,	0,    	1500	)
        self.h_pTBStop             	= ROOT.TH1F('pTBStop', 		';p_{T,b from #tilde{t}}  [GeV]',	150,	0,    	1500	)
        self.h_pTBStopPlus         	= ROOT.TH1F('pTBStopPlus', 	';p_{T,b from #tilde{t}^{+2/3}} [GeV]',	150,	0, 	1500	)
        self.h_pTBStopMinus    	   	= ROOT.TH1F('pTBStopMinus', 	';p_{T,b from #tilde{t}^{-2/3}} [GeV]',	150,	0, 	1500	)
        self.h_pTBChi              	= ROOT.TH1F('pTBChi', 		';p_{T,b from #chi^{#pm}}  [GeV]',	150,	0, 	1500	)
        self.h_pTBChiPlus          	= ROOT.TH1F('pTBChiPlus', 	';p_{T,b from #chi^{+}} [GeV]',		150,	0, 	1500	)
        self.h_pTBChiMinus         	= ROOT.TH1F('pTBChiMinus', 	';p_{T,b from #chi^{-}} [GeV]',		150,	0, 	1500	) 

        self.h_etaStop             	= ROOT.TH1F('etaStop', 		';#eta_{#tilde{t}',			80,	-8,	8	)
        self.h_etaStopPlus         	= ROOT.TH1F('etaStopPlus', 	';#eta_{#tilde{t}^{+2/3}',		80,	-8,	8	)
        self.h_etaStopMinus        	= ROOT.TH1F('etaStopMinus', 	';eta_{#tilde{t}^{-2/3}',		80,	-8,	8	) 
        self.h_etaChi              	= ROOT.TH1F('etaChi', 		';#eta_{#tilde{#chi^{#pm}}',		80,	-8,	8	)
        self.h_etaChiPlus          	= ROOT.TH1F('etaChiPlus',	';#eta_{#chi^{+}}',			80,	-8,	8	)
        self.h_etaChiMinus         	= ROOT.TH1F('etaChiMinus',	';#eta_{#chi^{-}}',			80,	-8,	8	)
        self.h_etaBStop            	= ROOT.TH1F('etaBStop',		';eta_{b from #tilde{t}}',		80,	-8,	8	)       
        self.h_etaBStopPlus        	= ROOT.TH1F('etaBStopPlus',	';eta_{b from #tilde{t}^{+2/3}}',	80,	-8,	8	)
        self.h_etaBStopMinus       	= ROOT.TH1F('etaBStopMinus',	';eta_{b from #tilde{t}^{-2/3}}',	80,	-8,	8	)
        self.h_etaBChi             	= ROOT.TH1F('etaBChi' , 	';eta_{b from #chi^{#pm}}',		80,	-8,	8	)
        self.h_etaBChiPlus         	= ROOT.TH1F('etaBChiPlus', 	';eta_{b from #chi^{+}}',		80,	-8,	8	)
        self.h_etaBChiMinus        	= ROOT.TH1F('etaBChiMinus', 	';eta_{b from #chi^{-}}',		80,	-8,	8	)

        self.h_dEtaBChi         	= ROOT.TH1F('dEtaBChi',		';|#Delta#eta_{b,#chi^{#pm}}|',		50,	0,	5.0   	)
        self.h_dPhiBChi         	= ROOT.TH1F('dPhiBChi',		';|#Delta#phi_{b,#chi^{#pm}}|',		50,	0,	5.0   	)
        self.h_nJetsChiMerged   	= ROOT.TH1F('nJetsChiMerged',	';N_{j} matched with #chi^{#pm}',	4,	0,	4   	)
        self.h_dRBB             	= ROOT.TH1F('dRBB',		';#DeltaR_{b,b}',			50,	0,	7    	)
        self.h_dEtaBB           	= ROOT.TH1F('dEtaBB',		';|#Delta#eta_{b,b}|',			50,	0,	5.0   	)
        self.h_dPhiBB           	= ROOT.TH1F('dPhiBB',		';|#Delta#phi_{b,b}|',			50,	0,	5.0   	)
        self.h_passDijet        	= ROOT.TH1F('passDijet',	';Pass Dijet Search Requirements',	2,	0,	2	)
        self.h_dEtaWJs			= ROOT.TH1F('dEtaWJs',		';|#Delta#eta_{j,j}| (wide jets)',	50,	0,	5.0	)
        self.h_mWJs			= ROOT.TH1F('mWJs',		';m_jj (wide jets) [GeV]',		150,    0,      3000    )	
        self.h_dRChiMax         	= ROOT.TH1F('dRChiMax',		';#Delta R_{#chi^{#pm},max}',		50,	0,	8   	)
        self.h_dRBChi         		= ROOT.TH1F('dRBChi',		';#Delta R_{b,#chi^{#pm}}',		50,	0,	7  	)

        self.h_pTBVsChi		 	= ROOT.TH2F('pTBVsChi',		';p_{T,#tilde{#chi^{#pm}}} [GeV];p_{T,b} [GeV]',		75,	0,	1500,	75,	0,	1500	)
        self.h_dEtaVsPTStop      	= ROOT.TH2F('dEtaVsPTStop',	';p_{T,#tilde{t}} [GeV];|#Delta#eta_{b,#tilde{#chi^{#pm}}}|',	75,	0,	1500,	80,	0,	8	)
        self.h_dEtaVsPTStopRatio 	= ROOT.TH2F('dEtaVsPTStopRatio',';p_{T,#tilde{t}} / #Delta m_{#tilde{t},#tilde{#chi^{#pm}}};|#Delta#eta_{b,#tilde{#chi^{#pm}}}|',40,0,2,80,0,8)

        # Matching
        self.h_jetMatch          	= ROOT.TH1F('jetMatch',		';Jet matched',				2,	0,	2   	)
        self.h_leadJetMatch      	= ROOT.TH1F('leadJetMatch',	';Correct lead jet matched',		2,	0,	2   	)

        # Trigger
        self.h_L1_HTT450er	 	= ROOT.TH1F('L1_HTT450er',	';L1_HTT450er',				2,	0,	2	)

        #-----------------------------------------------------------------------
        # RECO
        #-----------------------------------------------------------------------
        
        self.h_HT 			= ROOT.TH1F('HT',	 	';H_{T} [GeV]',			  	150,	0,	3000	)
        self.h_nJets	  		= ROOT.TH1F('nJets',  	  	';N_{j}',  				20, 	0,	20  	)
        self.h_nb  			= ROOT.TH1F('nb',  	  	';n_{b} (tight)',  			5,	0,	5  	)
        self.h_mAll             	= ROOT.TH1F('mAll',       	';m_{#sum j} [GeV]',                   	150,	0,	3000	)
        self.h_m4   			= ROOT.TH1F('m4',   	  	';m_{4j} [GeV]',   			150,	0,	3000	)
        self.h_m3   			= ROOT.TH1F('m3',   	  	';m_{3j} [GeV]',   			150,	0,	3000	)
        self.h_pT1          		= ROOT.TH1F('pT1',           	';p_{T,1} [GeV]',                    	150,	0,	1500 	)
        self.h_pT2          		= ROOT.TH1F('pT2',           	';p_{T,2} [GeV]',                    	150,	0,	1500 	)
        self.h_pT3          		= ROOT.TH1F('pT3',           	';p_{T,3} [GeV]',                    	150,	0,	1500 	)
        self.h_pT4          		= ROOT.TH1F('pT4',           	';p_{T,4} [GeV]',                    	150,	0,	1500 	)
        self.h_eta1          		= ROOT.TH1F('eta1',       	';#eta_{1}',     			80,	-8,	8	)
        self.h_eta2          		= ROOT.TH1F('eta2',       	';#eta_{2}',     			80,	-8,	8	)
        self.h_eta3          		= ROOT.TH1F('eta3',       	';#eta_{3}',     			80,	-8,	8	)
        self.h_eta4          		= ROOT.TH1F('eta4',       	';#eta_{4}',     			80,	-8,	8	)
        self.h_m3NoLead      		= ROOT.TH1F('m3NoLead',   	';m_{3j} (excl. leading) [GeV]',	150,	0,	3000	)
        self.h_dEta12			= ROOT.TH1F('dEta12',           ';|#Delta#eta_{1,2}|',                  50,     0,     	5     	)
        self.h_dPhi12			= ROOT.TH1F('dPhi12',           ';|#Delta#phi_{1,2}|',                  50,     0,      5     	)
        self.h_dR12			= ROOT.TH1F('dR12',           	';#Delta R_{1,2}',                    	50,     0,      7     	)

        # Add histograms to analysis object
        for h in list(vars(self)):
          if h[0:2] == 'h_':  self.addObject(getattr(self,h))

    def analyze(self, event):

        dRMatch = 0.1

        jets       = list(Collection(event,"Jet"))
        genParts   = list(Collection(event,"GenPart"))
        genAK4Jets = list(Collection(event,"GenJet"))

        # Reco jet cuts
        jets = filter(lambda x: x.pt > 30 and abs(x.eta) < 2.4,jets)

        # Get only outgoing particles of the hardest subprocess
        #gens = list(filter(lambda x: ((x.statusFlags >> 13) & 1) and (not((x.statusFlags >> 7)) & 1) and ((x.statusFlags >> 8) & 1), genParts)) 
        gens = filter(lambda x: (((x.statusFlags >> 13) & 1) and ((x.statusFlags >> 8) & 1)) and not (((abs(x.pdgId) == 1) or (abs(x.pdgId) == 3)) and ((x.statusFlags >> 11) & 1)), genParts)
        #gens = filter(lambda x: ((x.statusFlags >> 13) & 1),genParts)
        #genAK4Jets = filter(lambda x: abs(x.pdgId) <= 5,gens)
        #genAK4Jets = sorted(genAK4Jets, key=lambda x: x.pt, reverse=True)

        #print(len(genAK4Jets))
        #tot = ROOT.TLorentzVector()
        #for g in genAK4Jets:
        #  #if g.pdgId == -5 and ((g.statusFlags >> 13) &1):
        #  #  print('b from stop')
        #  #  tot += g.p4()
        #  #elif abs(g.pdgId) == 1000024 and not ((g.statusFlags >> 13) &1):
        #  #  print('chi')
        #  #  tot += g.p4()
        #  #if g.genPartIdxMother >= 0 and genParts[g.genPartIdxMother].pdgId == -1000006 and g.pdgId != -1000006:
        #  #  print(g.pdgId)
        #  #  tot += g.p4()
        #  tot += g.p4()
        #print(tot.M())
        #return True

        if isSignal:

          # Triggers
          self.h_L1_HTT450er.Fill(event.L1_HTT450er)

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
            genQuarks = [genBStop,genBChi,genD,genS]
          else: print('WARNING: No stop found in event')

          #-----------------------------------------------------------------------
          # GEN
          #-----------------------------------------------------------------------

          # Gen kinematics
          self.h_pTStop.Fill(genStop.pt) 
          self.h_pTChi.Fill(genChi.pt)
          self.h_pTBStop.Fill(genBStop.pt)
          self.h_pTBChi.Fill(genBChi.pt)
          if stopPlus:
            self.h_pTBStopPlus.Fill(genBStopPlus.pt)
            self.h_etaBStopPlus.Fill(genBStopPlus.eta)
            self.h_pTBChiPlus.Fill(genBChiPlus.pt)
            self.h_etaBChiPlus.Fill(genBChiPlus.eta)
            self.h_pTStopPlus.Fill(genStopPlus.pt)
            self.h_etaStopPlus.Fill(genStopPlus.eta)
            self.h_pTChiPlus.Fill(genChiPlus.pt)
            self.h_etaChiPlus.Fill(genChiPlus.eta)
          else:
            self.h_pTBStopMinus.Fill(genBStopMinus.pt)
            self.h_etaBStopMinus.Fill(genBStopMinus.eta)
            self.h_pTBChiMinus.Fill(genBChiMinus.pt)
            self.h_etaBChiMinus.Fill(genBChiMinus.eta)
            self.h_pTStopMinus.Fill(genStopMinus.pt)
            self.h_etaStopMinus.Fill(genStopMinus.eta)
            self.h_pTChiMinus.Fill(genChiMinus.pt)
            self.h_etaChiMinus.Fill(genChiMinus.eta)
          dEtaBChi = abs(genBStop.eta - genChi.eta)
          self.h_dEtaBChi.Fill(dEtaBChi)
          self.h_dPhiBChi.Fill(abs(genBStop.p4().DeltaPhi(genChi.p4())))
          self.h_dRBChi.Fill(abs(genBStop.p4().DeltaR(genChi.p4())))
          dRChiMax = max(genBChi.p4().DeltaR(genD.p4()),
                                   genBChi.p4().DeltaR(genS.p4()),
                                   genD.p4().DeltaR(genS.p4()))
          self.h_dRChiMax.Fill(dRChiMax)
          self.h_dRBB.Fill(genBChi.p4().DeltaR(genBStop.p4()))
          self.h_dEtaBB.Fill(abs(genBChi.eta - genBStop.eta))
          self.h_dPhiBB.Fill(abs(genBChi.p4().DeltaPhi(genBStop.p4())))
          self.h_pTBVsChi.Fill(genChi.pt,genBStop.pt) 
          self.h_dEtaVsPTStop.Fill(genChi.pt,abs(genChi.eta - genBStop.eta))
          self.h_dEtaVsPTStopRatio.Fill(genChi.pt / (genStop.p4().M() - genChi.p4().M()),abs(genChi.eta - genBStop.eta))
          self.h_passDijet.Fill(1 if (dRChiMax < 1.1 and dEtaBChi < 1.1) else 0)

          #pT and eta of the gen AK4 jets
          for i,j in enumerate(genAK4Jets):
           if i == 0:
            self.h_pT1Gen.Fill(j.pt)
            self.h_eta1Gen.Fill(j.eta)
           if i == 1:
            self.h_pT2Gen.Fill(j.pt)
            self.h_eta2Gen.Fill(j.eta)
           if i == 2:
            self.h_pT3Gen.Fill(j.pt)
            self.h_eta3Gen.Fill(j.eta)
           if i == 3:
            self.h_pT4Gen.Fill(j.pt)
            self.h_eta4Gen.Fill(j.eta)

          # Wide jets
          if len(genAK4Jets) >= 2:
            genJ1 = genAK4Jets[0].p4()
            genJ2 = genAK4Jets[1].p4()
            genWJ1 = genJ1
            genWJ2 = genJ2
            #passDijet = 1
            for j in genAK4Jets[2:]:
              j = j.p4()
              dR1 = j.DeltaR(genJ1)
              dR2 = j.DeltaR(genJ2)
              if dR1 < 1.1 and dR1 < dR2: genWJ1 += j
              elif dR2 < 1.1 and dR2 < dR1: genWJ2 += j
            self.h_mWJs.Fill((genWJ1 + genWJ2).M())
            dEtaWJ = abs(genWJ1.Eta() - genWJ2.Eta())
            self.h_dEtaWJs.Fill(dEtaWJ)
            #if dEtaWJ > 1.1: passDijet = 0
            #self.h_passDijet.Fill(passDijet)

        #-----------------------------------------------------------------------
        # RECO
        #-----------------------------------------------------------------------

        sumJet          = ROOT.TLorentzVector()
        sumJet4         = ROOT.TLorentzVector()
        sumJet3         = ROOT.TLorentzVector()
        sumJet3NoLead   = ROOT.TLorentzVector()

        # n jets
        self.h_nJets.Fill(len(jets))

        # 4 leading jet pTs
        HT = 0; nb = 0
        for i,j in enumerate(jets):
          sumJet += j.p4()
          HT += j.pt
          if j.btagDeepB > 0.7665: nb += 1
          if i < 4:           sumJet4 += j.p4()
          if i < 3:           sumJet3 += j.p4()
          if i > 0 and i < 4: sumJet3NoLead += j.p4()
          if i == 0:
           self.h_pT1.Fill(j.pt)
           self.h_eta1.Fill(j.eta)
          if i == 1:
           self.h_pT2.Fill(j.pt)
           self.h_eta2.Fill(j.eta)
          if i == 2:
           self.h_pT3.Fill(j.pt)
           self.h_eta3.Fill(j.eta)
          if i == 3:
           self.h_pT4.Fill(j.pt)
           self.h_eta4.Fill(j.eta)
        if len(jets) >= 4:
          self.h_m4.Fill(sumJet4.M())
          self.h_m3NoLead.Fill(sumJet3NoLead.M())
        if len(jets) >= 3:
          self.h_m3.Fill(sumJet3.M())
        if len(jets) >= 2:
          self.h_dEta12.Fill(abs(jets[0].eta - jets[1].eta))
          self.h_dPhi12.Fill(abs(jets[0].p4().DeltaPhi(jets[1].p4())))
          self.h_dR12.Fill(abs(jets[0].p4().DeltaR(jets[1].p4())))
        if len(jets) >= 1:
          self.h_HT.Fill(HT)
          self.h_nb.Fill(nb)
          self.h_mAll.Fill(sumJet.M())

        return True

parser = argparse.ArgumentParser(description='Single Stop Analyzer')
parser.add_argument('--sample',type=str,default='signal',choices=['signal','TT','QCD'],help='Sample to run over')
parser.add_argument('--tag',type=str,default='test',help='Tag for output label')
parser.add_argument('-n',type=int,default=1,help='Sample index to run over for backgrounds')
args = parser.parse_args()

outputPath = 'output/Run3/{}'.format(args.tag)
if not os.path.exists(outputPath):
  os.makedirs(outputPath)

if args.sample == 'TT': sampleFile = 'TTToHadronic.txt'
elif args.sample == 'QCD': sampleFile = 'QCDBEnriched.txt'
elif args.sample != 'signal': print('ERROR: Unexpected sample argument')

if args.sample == 'signal':

  isSignal = 1
  preselection = '' #"Jet_pt[0] > 250"
  #files = [
  #"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAOD/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/NANOAODSIM/PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/40000/2CE738F9-C212-E811-BD0E-EC0D9A8222CE.root"
  #"file:/eos/uscms/store/user/dmahon/condor/RPVSingleStopMC/NANOAOD/NANOAOD-300_200-?.root"
  #]
  points = ['200_100','300_100','300_200','500_100','500_200','500_400','700_100','700_400','700_600','1000_100','1000_400','1000_900','1500_100','1500_600','1500_1400','2000_100','2000_900','2000_1900','700_200','1000_200','1500_200','1500_400','2000_200','2000_400']
  #points = ['2000_100']
  for masses in points:
    files = glob.glob('/eos/uscms/store/user/dmahon/condor/RPVSingleStopRun3MC/NANOAOD-ALL/NANOAOD-{}.root'.format(masses))
    #files = glob.glob('/eos/uscms/store/user/dmahon/condor/RPVSingleStopMC/NANOAOD/NANOAOD-{}-*.root'.format(masses))
    files = ['root://cmsxrootd.fnal.gov/' + x.replace('/eos/uscms','') for x in files]
    #files = ['file:/uscms_data/d3/dmahon/RPVSingleStopRun3Patched/NANOAOD/CMSSW_12_4_5/test_2000_100-1.root']
    #files = ['/uscms_data/d3/dmahon/RPVSingleStopRun3Patched/NANOAOD/files/NANOAOD-{}.root'.format(masses)]
    p = PostProcessor(".", files, cut=preselection, branchsel=None, modules=[
                    ExampleAnalysis(isSignal)], noOut=True, histFileName='{}/{}_{}.root'.format(outputPath,args.sample,masses), histDirName="plots",
                    maxEntries=None)
    p.run()

else:

  isSignal = 0
  preselection = ''
  files = open('samples/{}'.format(sampleFile)).read().split('\n')
  #files = glob.glob('/eos/uscms/store/user/dmahon/condor/RPVSingleStopMC/NANOAOD/NANOAOD-{}-*.root'.format(masses))
  files = [['root://cmsxrootd.fnal.gov/' + x.replace('/eos/uscms','') for x in files][:-1][args.n - 1]]
  p = PostProcessor(".", files, cut=preselection, branchsel=None, modules=[
                  ExampleAnalysis(isSignal)], noOut=True, histFileName='{}/{}-{}.root'.format(outputPath,args.sample,args.n), histDirName="plots",
                  maxEntries=None)
  p.run() 
