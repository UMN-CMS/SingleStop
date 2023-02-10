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
from match_algos import matcher

ROOT.PyConfig.IgnoreCommandLineOptions = True
DR_MATCH=0.2
JET_COMBINATORICS_COUNT=6





def leadJetFromStop(histo, jet, stop):
    return stop.p4().DeltaR(jet.p4()) < DR_MATCH

def doJetMatching(jets, particles):
    best_permutation = matcher(jets, particles)
    for i in range(len(particles)):
        print("Particle {} matched with jet {} with distance {}".format(particles[i].pdgId , 
            best_permitation[i], particles[i].p4().DeltaR(jets[best_permutation[i]].p4())
            ))









class ExampleAnalysis(Module):

    def __init__(self,isSignal,MCCampaign):
        self.writeHistFile = True
        self.isSignal = isSignal
        self.MCCampaign = MCCampaign

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
        self.h_pTStop              	= ROOT.TH1F('pTStop', 		';p_{T,#tilde{t}}}  [GeV]',		150,	0,    	1500	)
        self.h_pTStopPlus          	= ROOT.TH1F('pTStopPlus', 	';p_{T,#tilde{t}^{+2/3}}  [GeV]',	150,	0,    	1500	)
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

        self.h_etaStop             	= ROOT.TH1F('etaStop', 		';#eta_{#tilde{t}}',			80,	-8,	8	)
        self.h_etaStopPlus         	= ROOT.TH1F('etaStopPlus', 	';#eta_{#tilde{t}^{+2/3}}',		80,	-8,	8	)
        self.h_etaStopMinus        	= ROOT.TH1F('etaStopMinus', 	';eta_{#tilde{t}^{-2/3}}',		80,	-8,	8	) 
        self.h_etaChi              	= ROOT.TH1F('etaChi', 		';#eta_{#tilde{#chi^{#pm}}}',		80,	-8,	8	)
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
        self.h_dRBB             	= ROOT.TH1F('dRBB',		';#DeltaR_{b,b}',			35,	0,	7    	)
        self.h_dEtaBB           	= ROOT.TH1F('dEtaBB',		';|#Delta#eta_{b,b}|',			50,	0,	5.0   	)
        self.h_dPhiBB           	= ROOT.TH1F('dPhiBB',		';|#Delta#phi_{b,b}|',			50,	0,	5.0   	)
        self.h_passDijet        	= ROOT.TH1F('passDijet',	';Pass Dijet Search Requirements',	2,	0,	2	)
        self.h_dEtaWJs			= ROOT.TH1F('dEtaWJs',		';|#Delta#eta_{j,j}| (wide jets)',	50,	0,	5.0	)
        self.h_mWJs			= ROOT.TH1F('mWJs',		';m_jj (wide jets) [GeV]',		150,    0,      3000    )	
        self.h_dRChiMax         	= ROOT.TH1F('dRChiMax',		';#Delta R_{#chi^{#pm},max}',		40,	0,	8   	)
        self.h_dRBChi         		= ROOT.TH1F('dRBChi',		';#Delta R_{b,#chi^{#pm}}',		35,	0,	7  	)

        self.h_pTBVsChi		 	= ROOT.TH2F('pTBVsChi',		';p_{T,#tilde{#chi^{#pm}}} [GeV];p_{T,b} [GeV]',		75,	0,	1500,	75,	0,	1500	)
        self.h_dEtaVsPTStop      	= ROOT.TH2F('dEtaVsPTStop',	';p_{T,#tilde{t}} [GeV];|#Delta#eta_{b,#tilde{#chi^{#pm}}}|',	75,	0,	1500,	80,	0,	8	)
        self.h_dEtaVsPTStopRatio 	= ROOT.TH2F('dEtaVsPTStopRatio',';p_{T,#tilde{t}} / #Delta m_{#tilde{t},#tilde{#chi^{#pm}}};|#Delta#eta_{b,#tilde{#chi^{#pm}}}|',40,0,2,80,0,8)

        # Matching
        self.h_jetMatch          	= ROOT.TH1F('jetMatch',		';Jet matched',				2,	0,	2   	)
        self.h_leadJetMatch      	= ROOT.TH1F('leadJetMatch',	';Correct lead jet matched',		2,	0,	2   	)

        self.h_leadBFromStop              = ROOT.TH1F("leadBFromStop" , ";leadBFromStop", 2, 0, 2)
        self.h_leadBFromChi              = ROOT.TH1F("leadBFromChi" , ";leadBFromChi", 2, 0, 2)
        self.h_leadJetFromStop             = ROOT.TH1F("leadJetFromStop" , ";leadJetFromStop", 2, 0, 2)
        self.h_leadGenJetFromStop             = ROOT.TH1F("leadGenJetFromStop" , ";leadGenJetFromStop", 2, 0, 2)

        self.h_stopBGenJetMatchOrdinal     = ROOT.TH1F("stopBGenJetMatchOrdinal", ";stopBJetMatchOrdinal",            10,     0,      10      )
        self.h_chiBGenJetMatchOrdinal      = ROOT.TH1F("chiBGenJetMatchOrdinal", ";chiBJetMatchOrdinal",            10,     0,      10      )
        self.h_chiqOneGenJetMatchOrdinal      = ROOT.TH1F("chiqOneGenJetMatchOrdinal", ";chiqOneGenJetMatchOrdinal",            10,     0,      10      )
        self.h_chiqAnyGenJetMatchOrdinal      = ROOT.TH1F("chiqAnyGenJetMatchOrdinal", ";chiqAnyGenJetMatchOrdinal",            10,     0,      10      )

        self.h_stopBJetMatchOrdinal     = ROOT.TH1F("stopBJetMatchOrdinal", ";stopBJetMatchOrdinal",            10,     0,      10      )
        self.h_chiBJetMatchOrdinal      = ROOT.TH1F("chiBJetMatchOrdinal", ";chiBJetMatchOrdinal",            10,     0,      10      )
        self.h_chiqOneJetMatchOrdinal      = ROOT.TH1F("chiqOneJetMatchOrdinal", ";chiqOneJetMatchOrdinal",            10,     0,      10      )
        self.h_chiqTwoJetMatchOrdinal      = ROOT.TH1F("chiqTwoJetMatchOrdinal", ";chiqTwoJetMatchOrdinal",            10,     0,      10      )
        self.h_chiqAnyJetMatchOrdinal      = ROOT.TH1F("chiqAnyJetMatchOrdinal", ";chiqAnyJetMatchOrdinal",            10,     0,      10      )


        self.h_numTopFourGenMatched      = ROOT.TH1F("numTopFourGenMatched", "numTopFourGenMatched",            4,     0,      4     )
        self.h_numTopThreeGenMatched      = ROOT.TH1F("numTopThreeGenMatched", "numTopThreeGenMatched",            4,     0,      4     )

        self.h_numTopFourRecoMatched      = ROOT.TH1F("numTopFourRecoMatched", "numTopFourRecoMatched",            4,     0,      4     )
        self.h_numTopThreeRecoMatched      = ROOT.TH1F("numTopThreeRecoMatched", "numTopThreeRecoMatched",            4,     0,      4     )


        # Trigger
        self.h_L1_HTT450er	 	= ROOT.TH1F('L1_HTT450er',	';L1_HTT450er',				2,	0,	2	)

        #-----------------------------------------------------------------------
        # RECO
        #-----------------------------------------------------------------------
        
        self.h_HT 			= ROOT.TH1F('HT',	 	';H_{T} [GeV]',			  	150,	0,	3000	)
        self.h_MET                      = ROOT.TH1F('MET',              ';p^{miss}_{T} [GeV]',                  50,     0,      1000    )
        self.h_dPhiMET1                 = ROOT.TH1F('dPhiMET1',         ';|#Delta#phi_{p^{miss}_{T},1}|',       50,     0,      5       )
        self.h_dPhiMET2                 = ROOT.TH1F('dPhiMET2',         ';|#Delta#phi_{p^{miss}_{T},1}|',       50,     0,      5       )
        self.h_dPhiMET3                 = ROOT.TH1F('dPhiMET3',         ';|#Delta#phi_{p^{miss}_{T},1}|',       50,     0,      5       )
        self.h_dPhiMET4                 = ROOT.TH1F('dPhiMET4',         ';|#Delta#phi_{p^{miss}_{T},1}|',       50,     0,      5       )
        self.h_nJets	  		= ROOT.TH1F('nJets',  	  	';N_{j}',  				20, 	0,	20  	)
        self.h_nbLoose 			= ROOT.TH1F('nbLoose', 	  	';n_{b} (loose)',  			7,	0,	7  	)
        self.h_nbMedium                 = ROOT.TH1F('nbMedium',         ';n_{b} (medium)',                      7,      0,      7       )
        self.h_nbTight                  = ROOT.TH1F('nbTight',          ';n_{b} (tight)',                       7,      0,      7       )
        self.h_ntLoose                  = ROOT.TH1F('ntLoose',          ';n_{t} (loose)',                       7,      0,      7       )
        self.h_ntMedium                 = ROOT.TH1F('ntMedium',         ';n_{t} (medium)',                      7,      0,      7       )
        self.h_ntTight                  = ROOT.TH1F('ntTight',          ';n_{t} (tight)',                       7,      0,      7       )
        self.h_mAll             	= ROOT.TH1F('mAll',       	';m_{#sum j} [GeV]',                   	150,	0,	3000	)
        self.h_m4   			= ROOT.TH1F('m4',   	  	';m_{4j} [GeV]',   			150,	0,	3000	)
        self.h_m3   			= ROOT.TH1F('m3',   	  	';m_{3j} [GeV]',   			150,	0,	3000	)
        self.h_m3NoLead                 = ROOT.TH1F('m3NoLead',         ';m_{3j} (excl. leading) [GeV]',        150,    0,      3000    )
        self.h_m3NoLeadOrSub            = ROOT.TH1F('m3NoLeadOrSub',    ';m_{3j} (excl. lead and sub) [GeV]',   150,    0,      3000    )
        self.h_pT1          		= ROOT.TH1F('pT1',           	';p_{T,1} [GeV]',                    	150,	0,	1500 	)
        self.h_pT2          		= ROOT.TH1F('pT2',           	';p_{T,2} [GeV]',                    	150,	0,	1500 	)
        self.h_pT3          		= ROOT.TH1F('pT3',           	';p_{T,3} [GeV]',                    	150,	0,	1500 	)
        self.h_pT4          		= ROOT.TH1F('pT4',           	';p_{T,4} [GeV]',                    	150,	0,	1500 	)
        self.h_eta1          		= ROOT.TH1F('eta1',       	';#eta_{1}',     			80,	-8,	8	)
        self.h_eta2          		= ROOT.TH1F('eta2',       	';#eta_{2}',     			80,	-8,	8	)
        self.h_eta3          		= ROOT.TH1F('eta3',       	';#eta_{3}',     			80,	-8,	8	)
        self.h_eta4          		= ROOT.TH1F('eta4',       	';#eta_{4}',     			80,	-8,	8	)
        self.h_dEta12			= ROOT.TH1F('dEta12',           ';|#Delta#eta_{1,2}|',                  50,     0,     	5     	)
        self.h_dPhi12			= ROOT.TH1F('dPhi12',           ';|#Delta#phi_{1,2}|',                  50,     0,      5     	)
        self.h_dR12			= ROOT.TH1F('dR12',           	';#Delta R_{1,2}',                    	35,     0,      7     	)
        self.h_m3Vsm4			= ROOT.TH2F('m3Vsm4',		';m_{4j} [GeV];m_{3j} [GeV]',		150,0,3000,150,0,3000   )
        self.h_m3NoLeadVsm4		= ROOT.TH2F('m3NoLeadVsm4',     ';m_{4j} [GeV];m_{3j} (excl. leading) [GeV]',150,0,3000,150,0,3000)

        # Add histograms to analysis object
        for h in list(vars(self)):
          if h[0:2] == 'h_':  self.addObject(getattr(self,h))

    def analyze(self, event):

        dRMatch = 0.1

        # Set b tagging WPs
        if self.MCCampaign == 'UL2016preVFP': 	bTagWPs = [0.0508,0.2598,0.6502]
        elif self.MCCampaign == 'UL2016postVFP': 	bTagWPs = [0.0480,0.2489,0.6377]
        elif self.MCCampaign == 'UL2017':     	bTagWPs = [0.0532,0.3040,0.7476]
        elif self.MCCampaign == 'UL2018':     	bTagWPs = [0.0490,0.2783,0.7100]

        # Get event collections
        jets           = filter(lambda x: x.pt > 30 and abs(x.eta) < 2.4,list(Collection(event,"Jet")))
        fatJets        = list(Collection(event,"FatJet"))
        looseTs        = filter(lambda x: x.particleNet_TvsQCD > 0.58,fatJets)
        mediumTs       = filter(lambda x: x.particleNet_TvsQCD > 0.80,fatJets)
        tightTs        = filter(lambda x: x.particleNet_TvsQCD > 0.97,fatJets)
        looseBs        = filter(lambda x: x.btagDeepFlavB > bTagWPs[0],jets)
        mediumBs       = filter(lambda x: x.btagDeepFlavB > bTagWPs[1],jets)
        tightBs        = filter(lambda x: x.btagDeepFlavB > bTagWPs[2],jets)
        genParts       = list(Collection(event,"GenPart"))
        genAK4Jets     = list(Collection(event,"GenJet"))
        goodElectrons  = filter(lambda x: x.cutBased == 4 and x.miniPFRelIso_all < 0.1 and x.pt > 30 and abs(x.eta) < 2.4,list(Collection(event,"Electron")))
        goodMuons      = filter(lambda x: x.mediumId and x.miniPFRelIso_all < 0.2 and x.pt > 30 and abs(x.eta) < 2.4,list(Collection(event,"Muon")))
        MET            = ROOT.TVector2()
        MET.SetMagPhi(event.MET_pt,event.MET_phi)

        # Get only outgoing particles of the hardest subprocess
        gens = filter(lambda x: (((x.statusFlags >> 13) & 1) and ((x.statusFlags >> 8) & 1)) and not (((abs(x.pdgId) == 1) or (abs(x.pdgId) == 3)) and ((x.statusFlags >> 11) & 1)), genParts)

        # Cuts
        if len(jets) < 4 or len(jets) > 5: return False
        if not jets[0].pt > 300: return False
        if not (event.HLT_PFHT1050 or event.HLT_AK8PFJet360_TrimMass30): return False
        if len(goodElectrons) != 0: return False
        if len(goodMuons) != 0 : return False
        if len(mediumBs) < 1: return False
        if not 2 < abs(jets[0].p4().DeltaR(jets[1].p4())) < 4: return False
        if len(tightTs) != 0: return False

        if self.isSignal:

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






          doJetMatching(genAK4Jets, genQuarks)



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

        sumJet             = ROOT.TLorentzVector()
        sumJet4            = ROOT.TLorentzVector()
        sumJet3            = ROOT.TLorentzVector()
        sumJet3NoLead      = ROOT.TLorentzVector()
        sumJet3NoLeadOrSub = ROOT.TLorentzVector()

        # n jets
        self.h_nJets.Fill(len(jets))
        self.h_ntLoose.Fill(len(looseTs))
        self.h_ntMedium.Fill(len(mediumTs))
        self.h_ntTight.Fill(len(tightTs))
        self.h_nbLoose.Fill(len(looseBs))
        self.h_nbMedium.Fill(len(mediumBs))
        self.h_nbTight.Fill(len(tightBs))

        # Event characteristics
        self.h_MET.Fill(MET.Mod())

        # 4 leading jet pTs
        HT = 0; nbLoose = 0; nbMedium = 0; nbTight=0
        for i,j in enumerate(jets):
          sumJet += j.p4()
          HT += j.pt
          if i < 4:           sumJet4 += j.p4()
          if i < 3:           sumJet3 += j.p4()
          if i > 0 and i < 4: sumJet3NoLead += j.p4()
          if i > 1 and i < 5: sumJet3NoLeadOrSub += j.p4()
          if i == 0:
           self.h_pT1.Fill(j.pt)
           self.h_eta1.Fill(j.eta)
           self.h_dPhiMET1.Fill(abs(j.p4().EtaPhiVector().DeltaPhi(MET)))
          if i == 1:
           self.h_pT2.Fill(j.pt)
           self.h_eta2.Fill(j.eta)
           self.h_dPhiMET2.Fill(abs(j.p4().EtaPhiVector().DeltaPhi(MET)))
          if i == 2:
           self.h_pT3.Fill(j.pt)
           self.h_eta3.Fill(j.eta)
           self.h_dPhiMET3.Fill(abs(j.p4().EtaPhiVector().DeltaPhi(MET)))
          if i == 3:
           self.h_pT4.Fill(j.pt)
           self.h_eta4.Fill(j.eta)
           self.h_dPhiMET4.Fill(abs(j.p4().EtaPhiVector().DeltaPhi(MET)))
        if len(jets) >= 5:
          self.h_m3NoLeadOrSub.Fill(sumJet3NoLeadOrSub.M())
        if len(jets) >= 4:
          self.h_m4.Fill(sumJet4.M())
          self.h_m3NoLead.Fill(sumJet3NoLead.M())
          self.h_m3Vsm4.Fill(sumJet4.M(),sumJet3.M())
          self.h_m3NoLeadVsm4.Fill(sumJet4.M(),sumJet3NoLead.M())
        if len(jets) >= 3:
          self.h_m3.Fill(sumJet3.M())
        if len(jets) >= 2:
          self.h_dEta12.Fill(abs(jets[0].eta - jets[1].eta))
          self.h_dPhi12.Fill(abs(jets[0].p4().DeltaPhi(jets[1].p4())))
          self.h_dR12.Fill(abs(jets[0].p4().DeltaR(jets[1].p4())))
        if len(jets) >= 1:
          self.h_HT.Fill(HT)
          self.h_mAll.Fill(sumJet.M())

        return True

parser = argparse.ArgumentParser(description='Single Stop Analyzer')
parser.add_argument('--sample',type=str,default='signal',choices=['signal','TT','TT2018','QCD','QCD2018','ZQQ2018','ST2018','WQQ2018','ZNuNu2018','Diboson2018'],help='Sample to run over')
parser.add_argument('--tag',type=str,default='test',help='Tag for output label')
parser.add_argument('-n',type=int,default=1,help='Sample index to run over for backgrounds')
parser.add_argument('--points',type=str,default='all',help='Signal point(s) to run over, comma separated in MSTOP_MCHI format; `"all`" to run over all available points')
args = parser.parse_args()

outputPath = 'output/{}'.format(args.tag)
if not os.path.exists(outputPath):
  os.makedirs(outputPath)

if   args.sample == 'TT':          sampleFile = 'TTToHadronic.txt'
elif args.sample == 'TT2018':      sampleFile = 'TTToHadronic2018.txt'
elif args.sample == 'QCD':         sampleFile = 'QCDBEnriched.txt'
elif args.sample == 'QCD2018':     sampleFile = 'QCDBEnriched2018.txt'
elif args.sample == 'ZQQ2018':     sampleFile = 'ZJetsToQQ2018.txt'
elif args.sample == 'ST2018':      sampleFile = 'STHadronic2018.txt'
elif args.sample == 'WQQ2018':     sampleFile = 'WJetsToQQ2018.txt'
elif args.sample == 'ZNuNu2018':   sampleFile = 'ZJetsToNuNu2018.txt'
elif args.sample == 'Diboson2018': sampleFile = 'Diboson2018.txt'
elif args.sample != 'signal': print('ERROR: Unexpected sample argument')

preselection = (
		'(Jet_pt[3] > 30) &&'
		'(Jet_pt[0] > 300) &&'
		'(HLT_PFHT1050 || HLT_AK8PFJet360_TrimMass30)'
               )

if args.sample == 'signal':

  allPoints = ['1000_400','1000_900','1500_600','1500_1400','2000_900','2000_1900','1000_600','1500_400','2000_400','2000_1400','1500_900','2000_600']
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
    p = PostProcessor(".", files, cut=preselection, branchsel=None, modules=[
                    ExampleAnalysis(isSignal=1,MCCampaign='UL2018')], noOut=True, histFileName='{}/{}_{}.root'.format(outputPath,args.sample,masses), histDirName="plots",
                    maxEntries=None)
    p.run()

else:

  files = open('samples/{}'.format(sampleFile)).read().split('\n')
  #files = glob.glob('/eos/uscms/store/user/dmahon/condor/RPVSingleStopMC/NANOAOD/NANOAOD-{}-*.root'.format(masses))
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
 
  p = PostProcessor(".", files, cut=preselection, branchsel=None, modules=[ ExampleAnalysis(isSignal=0,MCCampaign=MCCampaign)], noOut=True, histFileName='{}/{}-{}.root'.format(outputPath,args.sample,args.n), histDirName="plots", maxEntries=None)
  p.run() 
