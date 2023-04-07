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
ROOT.PyConfig.IgnoreCommandLineOptions = True

class ExampleAnalysis(Module):

    def __init__(self,isSignal,MCCampaign,isSkimmed):
        self.writeHistFile = True
        self.isSignal = isSignal
        self.MCCampaign = MCCampaign
        self.isSkimmed = isSkimmed

    def create2DHists(self,titles,labels,*args):
        for i,title in enumerate(titles):
          setattr(self,'h_{}'.format(title),ROOT.TH2D(title,'{}'.format(labels[i]),*args))

    def fill2DHists(self,titles,xVals,yVals,genWeight):
        for i,title in enumerate(titles):
          getattr(self,'h_{}'.format(title)).Fill(xVals[i],yVals[i],genWeight)

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        self.h_nEventsPostPre           = ROOT.TH1D('nEventsPostPre',   ';N_{events} (post-preselection)',      1,      0,      1       )
        self.h_cutflow                  = ROOT.TH1D('cutflow',    	';Cut',      				10,     0,      10      )

        #-----------------------------------------------------------------------
        # GEN
        #-----------------------------------------------------------------------

        # Gen-level (before parton showering)
        self.h_HTLHE                    = ROOT.TH1F('HTLHE',            ';H_{T,LHE} [GeV]',                     150,    0,      3000    )
        self.h_nQLHE                    = ROOT.TH1D('nQLHE',            ';N_{q} (LHE)',                         7,      0,      7       )
        self.h_nGLHE                    = ROOT.TH1D('nGLHE',            ';N_{g} (LHE)',                         7,      0,      7       )
        self.h_nJLHE                    = ROOT.TH1D('nJLHE',            ';N_{j} (LHE)',                         7,      0,      7       )

        # Jet kinematics
        self.h_pT1Gen          		= ROOT.TH1F('pT1Gen',		';p_{T,1}^{gen} [GeV]',			75,	0, 	1500 	)
        self.h_pT2Gen          		= ROOT.TH1F('pT2Gen',		';p_{T,2}^{gen} [GeV]',			75,	0,     	1500 	)
        self.h_pT3Gen            	= ROOT.TH1F('pT3GEn',		';p_{T,3}^{gen} [GeV]',			75,	0,     	1500 	)
        self.h_pT4Gen            	= ROOT.TH1F('pT4Gen',		';p_{T,4}^{gen} [GeV]',			75,	0,     	1500 	)
        self.h_eta1Gen     	       	= ROOT.TH1F('eta1Gen',		';#eta_{1}^{gen}',			80,	-8,	8	)
        self.h_eta2Gen  	       	= ROOT.TH1F('eta2Gen',		';#eta_{2}^{gen}',			80,	-8,    	8	)
        self.h_eta3Gen		       	= ROOT.TH1F('eta3Gen',		';#eta_{3}^{gen}',			80,	-8,   	8	)
        self.h_eta4Gen       		= ROOT.TH1F('eta4Gen',		';#eta_{4}^{gen}',			80,	-8,   	8	)

        # SUSY particle kinematics
        self.h_pTStop              	= ROOT.TH1F('pTStop', 		';p_{T,#tilde{t}}  [GeV]',		75,	0,    	1500	)
        self.h_pTStopPlus          	= ROOT.TH1F('pTStopPlus', 	';p_{T,#tilde{t}^{+2/3}}  [GeV]',	75,	0,    	1500	)
        self.h_pTStopMinus         	= ROOT.TH1F('pTStopMinus', 	';p_{T,#tilde{t}^{-2/3}} [GeV]',	75,	0,    	1500	)
        self.h_pTChi               	= ROOT.TH1F('pTChi',		';p_{T,#tilde{#chi}^{#pm}}  [GeV]',     75,	0,    	1500	)
        self.h_pTChiPlus           	= ROOT.TH1F('pTChiPlus', 	';p_{T,#chi^{+}} [GeV]',		75,	0,    	1500	)
        self.h_pTChiMinus          	= ROOT.TH1F('pTChiMinus', 	';p_{T,#chi^{-}} [GeV]',		75,	0,    	1500	)
        self.h_pTBStop             	= ROOT.TH1F('pTBStop', 		';p_{T,b from #tilde{t}}  [GeV]',	75,	0,    	1500	)
        self.h_pTBStopPlus         	= ROOT.TH1F('pTBStopPlus', 	';p_{T,b from #tilde{t}^{+2/3}} [GeV]',	75,	0, 	1500	)
        self.h_pTBStopMinus    	   	= ROOT.TH1F('pTBStopMinus', 	';p_{T,b from #tilde{t}^{-2/3}} [GeV]',	75,	0, 	1500	)
        self.h_pTBChi              	= ROOT.TH1F('pTBChi', 		';p_{T,b from #tilde{#chi}^{#pm}}  [GeV]',75,	0, 	1500	)
        self.h_pTBChiPlus          	= ROOT.TH2F('pTBChiPlus', 	';p_{T,b from #chi^{+}_{1}} [GeV];p_{T, b from #chi^{+}_{2}} [GeV]', 75,	0, 	1500, 	75,	0, 	1500	)
        self.h_pTBChiMinus         	= ROOT.TH2F('pTBChiMinus', 	';p_{T,b from #chi^{-}_{1}} [GeV];p_{T, b from #chi^{-}_{2}} [GeV]', 75,	0, 	1500,	75, 	0, 	1500	) 

        self.h_etaStop             	= ROOT.TH1F('etaStop', 		';#eta_{#tilde{t}}',			80,	-8,	8	)
        self.h_etaStopPlus         	= ROOT.TH1F('etaStopPlus', 	';#eta_{#tilde{t}^{+2/3}}',		80,	-8,	8	)
        self.h_etaStopMinus        	= ROOT.TH1F('etaStopMinus', 	';eta_{#tilde{t}^{-2/3}}',		80,	-8,	8	) 
        self.h_etaChi              	= ROOT.TH1F('etaChi', 		';#eta_{#tilde{#chi}^{#pm}}',		80,	-8,	8	)
        self.h_etaChiPlus          	= ROOT.TH1F('etaChiPlus',	';#eta_{#chi^{+}}',			80,	-8,	8	)
        self.h_etaChiMinus         	= ROOT.TH1F('etaChiMinus',	';#eta_{#chi^{-}}',			80,	-8,	8	)
        self.h_etaBStop            	= ROOT.TH1F('etaBStop',		';eta_{b from #tilde{t}}',		80,	-8,	8	)       
        self.h_etaBStopPlus        	= ROOT.TH1F('etaBStopPlus',	';eta_{b from #tilde{t}^{+2/3}}',	80,	-8,	8	)
        self.h_etaBStopMinus       	= ROOT.TH1F('etaBStopMinus',	';eta_{b from #tilde{t}^{-2/3}}',	80,	-8,	8	)
        self.h_etaBChi             	= ROOT.TH1F('etaBChi' , 	';eta_{b from #tilde{#chi}^{#pm}}',	80,	-8,	8	)
        self.h_etaBChiPlus         	= ROOT.TH2F('etaBChiPlus', 	';eta_{b from #chi^{+}_{1}};eta_{b from #chi&{+}_{2}}',	80,	-8,	8, 	80,	-8, 	8	)
        self.h_etaBChiMinus        	= ROOT.TH2F('etaBChiMinus', 	';eta_{b from #chi^{-}_{1}};eta_{b from #chi^{-}_{2}}',	80,	-8,	8, 	80, 	-8, 	8	)

        self.h_dEtaBChi         	= ROOT.TH1F('dEtaBChi',		';|#Delta#eta_{b,#tilde{#chi}^{#pm}}|',	50,	0,	5.0   	)
        self.h_dPhiBChi         	= ROOT.TH1F('dPhiBChi',		';|#Delta#phi_{b,#tilde{#chi}^{#pm}}|',	50,	0,	5.0   	)
        self.h_nJetsChiMerged   	= ROOT.TH1D('nJetsChiMerged',	';N_{j} matched with #tilde{#chi}^{#pm}',4,	0,	4   	)
        self.h_dRBB             	= ROOT.TH1F('dRBB',		';#DeltaR_{b,b}',			35,	0,	7    	)
        self.h_dEtaBB           	= ROOT.TH1F('dEtaBB',		';|#Delta#eta_{b,b}|',			50,	0,	5.0   	)
        self.h_dPhiBB           	= ROOT.TH1F('dPhiBB',		';|#Delta#phi_{b,b}|',			50,	0,	5.0   	)
        self.h_passDijet        	= ROOT.TH1D('passDijet',	';Pass Dijet Search Requirements',	2,	0,	2	)
        self.h_dEtaWJs			= ROOT.TH1F('dEtaWJs',		';|#Delta#eta_{j,j}| (wide jets)',	50,	0,	5.0	)
        self.h_mWJs			= ROOT.TH1F('mWJs',		';m_jj (wide jets) [GeV]',		150,    0,      3000    )	
        self.h_dRChiMax         	= ROOT.TH1F('dRChiMax',		';#Delta R_{#tilde{#chi}^{#pm},max}',	40,	0,	8   	)
        self.h_dRBChi         		= ROOT.TH1F('dRBChi',		';#Delta R_{b,#tilde{#chi}^{#pm}}',	35,	0,	7  	)

        self.h_pTBVsChi		 	= ROOT.TH2F('pTBVsChi',		';p_{T,#tilde{#tilde{#chi}^{#pm}}} [GeV];p_{T,b} [GeV]',	75,	0,	1500,	75,	0,	1500	)
        self.h_dEtaVsPTStop      	= ROOT.TH2F('dEtaVsPTStop',	';p_{T,#tilde{t}} [GeV];|#Delta#eta_{b,#tilde{#chi}^{#pm}}|',	75,	0,	1500,	80,	0,	8	)
        self.h_dEtaVsPTStopRatio 	= ROOT.TH2F('dEtaVsPTStopRatio',';p_{T,#tilde{t}} / #Delta m_{#tilde{t},#tilde{#chi}^{#pm}};|#Delta#eta_{b,#tilde{#chi}^{#pm}}|',40,0,2,80,0,8)

        # Trigger
        self.h_L1_HTT450er	 	= ROOT.TH1D('L1_HTT450er',	';L1_HTT450er',				2,	0,	2	)

        #-----------------------------------------------------------------------
        # RECO
        #-----------------------------------------------------------------------
        
        self.h_HT 			= ROOT.TH1F('HT',	 	';H_{T} [GeV]',			  	60,	0,	3000	)
        self.h_MET                      = ROOT.TH1F('MET',              ';p^{miss}_{T} [GeV]',                  50,     0,      1000    )
        self.h_nJets	  		= ROOT.TH1D('nJets',  	  	';N_{j}',  				20, 	0,	20  	)
        self.h_nbLoose 			= ROOT.TH1D('nbLoose', 	  	';n_{b} (loose)',  			7,	0,	7  	)
        self.h_nbMedium                 = ROOT.TH1D('nbMedium',         ';n_{b} (medium)',                      7,      0,      7       )
        self.h_nbTight                  = ROOT.TH1D('nbTight',          ';n_{b} (tight)',                       7,      0,      7       )
        self.h_ntLoose                  = ROOT.TH1D('ntLoose',          ';n_{t} (loose)',                       7,      0,      7       )
        self.h_ntMedium                 = ROOT.TH1D('ntMedium',         ';n_{t} (medium)',                      7,      0,      7       )
        self.h_ntTight                  = ROOT.TH1D('ntTight',          ';n_{t} (tight)',                       7,      0,      7       )
        self.h_nWLoose                  = ROOT.TH1D('nWLoose',          ';n_{W} (loose)',                       7,      0,      7       )
        self.h_nWMedium                 = ROOT.TH1D('nWMedium',         ';n_{W} (medium)',                      7,      0,      7       )
        self.h_nWTight                  = ROOT.TH1D('nWTight',          ';n_{W} (tight)',                       7,      0,      7       )
        self.h_ntDeepWP1                = ROOT.TH1D('ntDeepWP1',        ';n_{t} (DeepJet WP 1)',                7,      0,      7       )
        self.h_ntDeepWP2                = ROOT.TH1D('ntDeepWP2',        ';n_{t} (DeepJet WP 2)',                7,      0,      7       )
        self.h_ntDeepWP3                = ROOT.TH1D('ntDeepWP3',        ';n_{t} (DeepJet WP 3)',                7,      0,      7       )
        self.h_ntDeepWP4                = ROOT.TH1D('ntDeepWP4',        ';n_{t} (DeepJet WP 4)',                7,      0,      7       )
        self.h_nWDeepWP1                = ROOT.TH1D('nWDeepWP1',        ';n_{W} (DeepJet WP 1)',                7,      0,      7       )
        self.h_nWDeepWP2                = ROOT.TH1D('nWDeepWP2',        ';n_{W} (DeepJet WP 2)',                7,      0,      7       )
        self.h_nWDeepWP3                = ROOT.TH1D('nWDeepWP3',        ';n_{W} (DeepJet WP 3)',                7,      0,      7       )
        self.h_nWDeepWP4                = ROOT.TH1D('nWDeepWP4',        ';n_{W} (DeepJet WP 4)',                7,      0,      7       )
        self.h_mAll             	= ROOT.TH1F('mAll',       	';m_{#sum j} [GeV]',                   	60,	0,	3000	)
        self.h_m4   			= ROOT.TH1F('m4',   	  	';m_{4j} [GeV]',   			60,	0,	3000	)
        self.h_m3   			= ROOT.TH1F('m3',   	  	';m_{3j} [GeV]',   			60,	0,	3000	)
        self.h_m3NoLead                 = ROOT.TH1F('m3NoLead',         ';m_{3j} (excl. leading) [GeV]',        60,     0,      3000    )
        self.h_m3NoLeadOrSub            = ROOT.TH1F('m3NoLeadOrSub',    ';m_{3j} (excl. lead and sub) [GeV]',   60,     0,      3000    )
        self.h_pT1          		= ROOT.TH1F('pT1',           	';p_{T,1} [GeV]',                    	75,	0,	1500 	)
        self.h_pT2          		= ROOT.TH1F('pT2',           	';p_{T,2} [GeV]',                    	75,	0,	1500 	)
        self.h_pT3          		= ROOT.TH1F('pT3',           	';p_{T,3} [GeV]',                    	75,	0,	1500 	)
        self.h_pT4          		= ROOT.TH1F('pT4',           	';p_{T,4} [GeV]',                    	75,	0,	1500 	)
        self.h_pTb1                     = ROOT.TH1F('pTb1',             ';p_{T,b_{1}} (loose) [GeV]',           75,     0,      1500    )
        self.h_pTb2                     = ROOT.TH1F('pTb2',             ';p_{T,b_{2}} (loose) [GeV]',           75,     0,      1500    )
        self.h_pTb3                     = ROOT.TH1F('pTb3',             ';p_{T,b_{3}} (loose) [GeV]',           75,     0,      1500    )
        self.h_pTb4                     = ROOT.TH1F('pTb4',             ';p_{T,b_{4}} (loose) [GeV]',           75,     0,      1500    )
        self.h_eta1          		= ROOT.TH1F('eta1',       	';#eta_{1}',     			80,	-8,	8	)
        self.h_eta2          		= ROOT.TH1F('eta2',       	';#eta_{2}',     			80,	-8,	8	)
        self.h_eta3          		= ROOT.TH1F('eta3',       	';#eta_{3}',     			80,	-8,	8	)
        self.h_eta4          		= ROOT.TH1F('eta4',       	';#eta_{4}',     			80,	-8,	8	)

        self.h_pT1Frac                  = ROOT.TH1F('pT1Frac',          ';p_{T,1} / H_{T}',                     20,     0,      1       )
        self.h_pT2Frac                  = ROOT.TH1F('pT2Frac',          ';p_{T,2} / H_{T}',                     20,     0,      1       )
        self.h_pT3Frac                  = ROOT.TH1F('pT3Frac',          ';p_{T,3} / H_{T}',                     20,     0,      1       )
        self.h_pT4Frac                  = ROOT.TH1F('pT4Frac',          ';p_{T,4} / H_{T}',                     20,     0,      1       )
        self.h_pT5Frac                  = ROOT.TH1F('pT5Frac',          ';p_{T,5} / H_{T}',                     20,     0,      1       )

        self.h_pT1ETFracComp           = ROOT.TH1F('pT1ETFracComp',   ';p_{T,1} / E_{T,#tilde{#chi}^{#pm}} (compressed)',    60,     0,      3       )
        self.h_pT2ETFracComp           = ROOT.TH1F('pT2ETFracComp',   ';p_{T,2} / E_{T,#tilde{#chi}^{#pm}} (compressed)',    60,     0,      3       )
        self.h_pT3ETFracComp           = ROOT.TH1F('pT3ETFracComp',   ';p_{T,3} / E_{T,#tilde{#chi}^{#pm}} (compressed)',    60,     0,      3       )
        self.h_pT4ETFracComp           = ROOT.TH1F('pT4ETFracComp',   ';p_{T,4} / E_{T,#tilde{#chi}^{#pm}} (compressed)',    60,     0,      3       )
        self.h_pT5ETFracComp           = ROOT.TH1F('pT5ETFracComp',   ';p_{T,5} / E_{T,#tilde{#chi}^{#pm}} (compressed)',    60,     0,      3       )

        self.h_pT1ETFracUncomp         = ROOT.TH1F('pT1ETFracUncomp', ';p_{T,1} / E_{T,#tilde{#chi}^{#pm}} (uncompressed)',  60,     0,      3       )
        self.h_pT2ETFracUncomp         = ROOT.TH1F('pT2ETFracUncomp', ';p_{T,2} / E_{T,#tilde{#chi}^{#pm}} (uncompressed)',  60,     0,      3       )
        self.h_pT3ETFracUncomp         = ROOT.TH1F('pT3ETFracUncomp', ';p_{T,3} / E_{T,#tilde{#chi}^{#pm}} (uncompressed)',  60,     0,      3       )
        self.h_pT4ETFracUncomp         = ROOT.TH1F('pT4ETFracUncomp', ';p_{T,4} / E_{T,#tilde{#chi}^{#pm}} (uncompressed)',  60,     0,      3       )
        self.h_pT5ETFracUncomp         = ROOT.TH1F('pT5ETFracUncomp', ';p_{T,5} / E_{T,#tilde{#chi}^{#pm}} (uncompressed)',  60,     0,      3       )

        self.h_pT1FracChiComp           = ROOT.TH1F('pT1FracChiComp',   ';p_{T,1} / H_{T,#tilde{#chi}^{#pm}} (compressed)',  40,     0,      2       )
        self.h_pT2FracChiComp           = ROOT.TH1F('pT2FracChiComp',   ';p_{T,2} / H_{T,#tilde{#chi}^{#pm}} (compressed)',  40,     0,      2       )
        self.h_pT3FracChiComp           = ROOT.TH1F('pT3FracChiComp',   ';p_{T,3} / H_{T,#tilde{#chi}^{#pm}} (compressed)',  40,     0,      2       )
        self.h_pT4FracChiComp           = ROOT.TH1F('pT4FracChiComp',   ';p_{T,4} / H_{T,#tilde{#chi}^{#pm}} (compressed)',  40,     0,      2       )
        self.h_pT5FracChiComp           = ROOT.TH1F('pT5FracChiComp',   ';p_{T,5} / H_{T,#tilde{#chi}^{#pm}} (compressed)',  40,     0,      2       )

        self.h_pT1FracChiUncomp         = ROOT.TH1F('pT1FracChiUncomp', ';p_{T,1} / H_{T,#tilde{#chi}^{#pm}} (uncompressed)',40,     0,      2       )
        self.h_pT2FracChiUncomp         = ROOT.TH1F('pT2FracChiUncomp', ';p_{T,2} / H_{T,#tilde{#chi}^{#pm}} (uncompressed)',40,     0,      2       )
        self.h_pT3FracChiUncomp         = ROOT.TH1F('pT3FracChiUncomp', ';p_{T,3} / H_{T,#tilde{#chi}^{#pm}} (uncompressed)',40,     0,      2       )
        self.h_pT4FracChiUncomp         = ROOT.TH1F('pT4FracChiUncomp', ';p_{T,4} / H_{T,#tilde{#chi}^{#pm}} (uncompressed)',40,     0,      2       )
        self.h_pT5FracChiUncomp         = ROOT.TH1F('pT5FracChiUncomp', ';p_{T,5} / H_{T,#tilde{#chi}^{#pm}} (uncompressed)',40,     0,      2       )


        self.h_pTMeanComp               = ROOT.TH1F('pTMeanComp',       	';#mu(p_{T,1},p_{T,2},p_{T,3}) [GeV]',               50,    0,      1000    )
        self.h_pTSDComp                 = ROOT.TH1F('pTSDComp',     		';#sigma(p_{T,1},p_{T,2},p_{T,3}) [GeV]',            50,    0,      1000    )
        self.h_pTSDMeanFracComp         = ROOT.TH1F('pTSDMeanFracComp', 	';#frac{#sigma}{#mu}(p_{T,1},p_{T,2},p_{T,3})',      30,    0,      1.5     )
        self.h_pTMeanUncomp             = ROOT.TH1F('pTMeanUncomp',     	';#mu(p_{T,2},p_{T,3},p_{T,4}) [GeV]',               50,    0,      1000    )
        self.h_pTSDUncomp               = ROOT.TH1F('pTSDUncomp',   		';#sigma(p_{T,2},p_{T,3},p_{T,4}) [GeV]',            50,    0,      1000    )
        self.h_pTSDMeanFracUncomp       = ROOT.TH1F('pTSDMeanFracUncomp',	';#frac{#sigma}{#mu}(p_{T,2},p_{T,3},p_{T,4})',      30,    0,      1.5     )

        self.h_dEta12			= ROOT.TH1F('dEta12',           ';|#Delta#eta_{1,2}|',                  50,     0,     	5     	)
        self.h_dEta13                   = ROOT.TH1F('dEta13',           ';|#Delta#eta_{1,3}|',                  50,     0,      5       )
        self.h_dEta14                   = ROOT.TH1F('dEta14',           ';|#Delta#eta_{1,4}|',                  50,     0,      5       )
        self.h_dEta15                   = ROOT.TH1F('dEta15',           ';|#Delta#eta_{1,5}|',                  50,     0,      5       )
        self.h_dEta23                   = ROOT.TH1F('dEta23',           ';|#Delta#eta_{2,3}|',                  50,     0,      5       )
        self.h_dEta24                   = ROOT.TH1F('dEta24',           ';|#Delta#eta_{2,4}|',                  50,     0,      5       )
        self.h_dEta25                   = ROOT.TH1F('dEta25',           ';|#Delta#eta_{2,5}|',                  50,     0,      5       )
        self.h_dEta34                   = ROOT.TH1F('dEta34',           ';|#Delta#eta_{3,4}|',                  50,     0,      5       )
        self.h_dEta35                   = ROOT.TH1F('dEta35',           ';|#Delta#eta_{3,5}|',                  50,     0,      5       )
        self.h_dEta45                   = ROOT.TH1F('dEta45',           ';|#Delta#eta_{4,5}|',                  50,     0,      5       )

        self.h_dPhi12			= ROOT.TH1F('dPhi12',           ';|#Delta#phi_{1,2}|',                  50,     0,      5     	)
        self.h_dPhi13                   = ROOT.TH1F('dPhi13',           ';|#Delta#phi_{1,3}|',                  50,     0,      5       )
        self.h_dPhi14                   = ROOT.TH1F('dPhi14',           ';|#Delta#phi_{1,4}|',                  50,     0,      5       )
        self.h_dPhi15                   = ROOT.TH1F('dPhi15',           ';|#Delta#phi_{1,5}|',                  50,     0,      5       )
        self.h_dPhi23                   = ROOT.TH1F('dPhi23',           ';|#Delta#phi_{2,3}|',                  50,     0,      5       )
        self.h_dPhi24                   = ROOT.TH1F('dPhi24',           ';|#Delta#phi_{2,4}|',                  50,     0,      5       )
        self.h_dPhi25                   = ROOT.TH1F('dPhi25',           ';|#Delta#phi_{2,5}|',                  50,     0,      5       )
        self.h_dPhi34                   = ROOT.TH1F('dPhi34',           ';|#Delta#phi_{3,4}|',                  50,     0,      5       )
        self.h_dPhi35                   = ROOT.TH1F('dPhi35',           ';|#Delta#phi_{3,5}|',                  50,     0,      5       )
        self.h_dPhi45                   = ROOT.TH1F('dPhi45',           ';|#Delta#phi_{4,5}|',                  50,     0,      5       )

        jetCombos2D = [[''.join(map(str,i)),''.join(map(str,j))] for i,j in combinations(combinations([1,2,3,4],2),2)]
        self.create2DHists(['dPhi{}Vs{}'.format(i,j) for [i,j] in jetCombos2D],
                           [';|#Delta#phi_{{{}}}|;|#Delta#phi_{{{}}}|'.format(j,i) for [i,j] in jetCombos2D],
                           16,0,3.2,16,0,3.2)

        self.h_dR12			= ROOT.TH1F('dR12',           	';#Delta R_{1,2}',                    	35,     0,      7     	)
        self.h_dR13                     = ROOT.TH1F('dR13',             ';#Delta R_{1,3}',                      35,     0,      7       )
        self.h_dR14                     = ROOT.TH1F('dR14',             ';#Delta R_{1,4}',                      35,     0,      7       )
        self.h_dR15                     = ROOT.TH1F('dR15',             ';#Delta R_{1,5}',                      35,     0,      7       )
        self.h_dR23                     = ROOT.TH1F('dR23',             ';#Delta R_{2,3}',                      35,     0,      7       )
        self.h_dR24                     = ROOT.TH1F('dR24',             ';#Delta R_{2,4}',                      35,     0,      7       )
        self.h_dR25                     = ROOT.TH1F('dR25',             ';#Delta R_{2,5}',                      35,     0,      7       )
        self.h_dR34                     = ROOT.TH1F('dR34',             ';#Delta R_{3,4}',                      35,     0,      7       )
        self.h_dR35                     = ROOT.TH1F('dR35',             ';#Delta R_{3,5}',                      35,     0,      7       )
        self.h_dR45                     = ROOT.TH1F('dR45',             ';#Delta R_{4,5}',                      35,     0,      7       )

        self.h_dEtabb                   = ROOT.TH1F('dEtabb',           ';|#Delta#eta_{b_{1},b_{2}}| (loose)',  50,     0,      5       )
        self.h_dPhibb                   = ROOT.TH1F('dPhibb',           ';|#Delta#phi_{b_{1},b_{2}}| (loose)',  50,     0,      5       )
        self.h_dRbb                     = ROOT.TH1F('dRbb',             ';#Delta R_{b_{1},b_{2}} (loose)',      35,     0,      7       )

        self.h_dEtaRecoComp             = ROOT.TH1F('dEtaRecoComp',     ';|#Delta#eta_{#tilde{t},#tilde{#chi}^{#pm}}| (compressed)',  50,     0,      5       )
        self.h_dPhiRecoComp             = ROOT.TH1F('dPhiRecoComp',     ';|#Delta#phi_{#tilde{t},#tilde{#chi}^{#pm}}| (compressed)',  50,     0,      5       )
        self.h_dRRecoComp               = ROOT.TH1F('dRRecoComp',       ';#Delta R_{#tilde{t},#tilde{#chi}^{#pm}} (compressed)',      35,     0,      7       )
        self.h_dEtaRecoUncomp           = ROOT.TH1F('dEtaRecoUncomp',   ';|#Delta#eta_{#tilde{t},#tilde{#chi}^{#pm}}| (uncompressed)',50,     0,      5       )
        self.h_dPhiRecoUncomp           = ROOT.TH1F('dPhiRecoUncomp',   ';|#Delta#phi_{#tilde{t},#tilde{#chi}^{#pm}}| (uncompressed)',50,     0,      5       )
        self.h_dRRecoUncomp             = ROOT.TH1F('dRRecoUncomp',     ';#Delta R_{#tilde{t},#tilde{#chi}^{#pm}} (uncompressed)',    35,     0,      7       )

        self.h_m3Vsm4			= ROOT.TH2F('m3Vsm4',		';m_{4j} [GeV];m_{3j} [GeV]',		     150,0,3000,150,0,3000   )
        self.h_m3NoLeadVsm4		= ROOT.TH2F('m3NoLeadVsm4',     ';m_{4j} [GeV];m_{3j} (excl. leading) [GeV]',150,0,3000,150,0,3000)

        # Add histograms to analysis object
        for h in list(vars(self)):
          if h[0:2] == 'h_':  self.addObject(getattr(self,h))

    def analyze(self, event):

        dRMatch = 0.1

        # Get MC weight
        genWeight = 1 if event.genWeight > 0 else -1

        # Fill histograms before selections (but after any pre-selections)
        self.h_nEventsPostPre.Fill(0,genWeight)
        try: self.h_HTLHE.Fill(event.LHE_HT,genWeight)
        except RuntimeError: pass

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
          gens = filter(lambda x: (((x.statusFlags >> 13) & 1) and ((x.statusFlags >> 8) & 1)) and not (((abs(x.pdgId) == 1) or (abs(x.pdgId) == 5)) and ((x.statusFlags >> 11) & 1)), genParts)

        if not self.isSkimmed:
          goodElectrons  = filter(lambda x: x.cutBased == 4 and x.miniPFRelIso_all < 0.1 and x.pt > 30 and abs(x.eta) < 2.4,list(Collection(event,"Electron")))
          goodMuons      = filter(lambda x: x.mediumId and x.miniPFRelIso_all < 0.2 and x.pt > 30 and abs(x.eta) < 2.4,list(Collection(event,"Muon")))

        # Cuts
        if not self.isSkimmed:
          self.h_cutflow.Fill(0,genWeight)
          if not (event.HLT_PFHT1050 or event.HLT_AK8PFJet360_TrimMass30): return False
          self.h_cutflow.Fill(1,genWeight)
          if len(jets) > 0 and not jets[0].pt > 300: return False
          self.h_cutflow.Fill(2,genWeight)
          if len(jets) < 4 or len(jets) > 6: return False
          self.h_cutflow.Fill(3,genWeight)
          if len(goodElectrons) != 0 or len(goodMuons) != 0: return False
          self.h_cutflow.Fill(4,genWeight)
          if not 2 < abs(jets[0].p4().DeltaR(jets[1].p4())) < 4: return False
          self.h_cutflow.Fill(5,genWeight)
          if len(looseBs) < 2: return False
          self.h_cutflow.Fill(6,genWeight)
          if abs(looseBs[0].p4().DeltaR(looseBs[1].p4())) < 1: return False
          self.h_cutflow.Fill(7,genWeight)

        try: 
          self.h_nQLHE.Fill(event.LHE_Nuds + event.LHE_Nc + event.LHE_Nb,genWeight)
          self.h_nGLHE.Fill(event.LHE_Nglu,genWeight)
          self.h_nJLHE.Fill(event.LHE_Njets,genWeight)
        except RuntimeError: pass

        if self.isSignal:

          # Triggers
          self.h_L1_HTT450er.Fill(event.L1_HTT450er,genWeight)

          # Match gens to particle
          stopPlus = True
          if True in (g.pdgId == 1000006 for g in gens):
		  oneBChi = False
		  for g in gens:
		      if g.pdgId == 1000006: genStop = g
		      elif g.pdgId == 1000024: genChi = g
		      elif g.pdgId == 5: genBStop = g
		      elif g.pdgId == -5: 
			if not oneBChi:
				genBChi1 = g
				oneBChi = True
			else:
				genBChi2 = g
		      elif g.pdgId == -1: genD = g
		      else: print('WARNING: Unexpected particle with pdgId {}'.format(g.pdgId))
		  if genBChi2.pt > genBChi1.pt: genBChi1, genBChi2 = genBChi2, genBChi1
		  genStopPlus = genStop
	          genBStopPlus = genBStop
		  genChiPlus = genChi
		  genBChiPlus1 = genBChi1
		  genBChiPlus2 = genBChi2
          elif True in (g.pdgId == -1000006 for g in gens):
	      stopPlus = False
	      oneBChi = False
	      for g in gens:
		      if g.pdgId == -1000006: genStop = g
		      elif g.pdgId == -1000024: genChi = g
		      elif g.pdgId == -5: genBStop = g
		      elif g.pdgId == 5:
			      if not oneBChi:
				      genBChi1 = g
				      oneBChi = True
			      else:
				      genBChi2 = g
		      elif g.pdgId == 1: genD = g
		      else: print('WARNING: Unexpected particle with pdgId {}'.format(g.pdgId))
	      if genBChi2.pt > genBChi1.pt: genBChi1, genBChi2 = genBChi2, genBChi1  
	      genStopMinus = genStop
	      genBStopMinus = genBStop
	      genChiMinus = genChi
	      genBChiMinus1 = genBChi1
	      genBChiMinus2 = genBChi2
	      genQuarks = [genBStop,genBChi1,genD,genBChi2]

          else: print('WARNING: No stop found in event')

          #-----------------------------------------------------------------------
          # GEN
          #-----------------------------------------------------------------------

          # Gen kinematics
          self.h_pTStop.Fill(genStop.pt,genWeight)
          self.h_pTChi.Fill(genChi.pt,genWeight)
          self.h_pTBStop.Fill(genBStop.pt,genWeight)
          self.h_pTBChi.Fill(genBChi1.pt,genWeight)
          if stopPlus:
            self.h_pTBStopPlus.Fill(genBStopPlus.pt,genWeight)
            self.h_etaBStopPlus.Fill(genBStopPlus.eta,genWeight)
            self.h_pTBChiPlus.Fill(genBChiPlus1.pt, genBChiPlus2.pt, genWeight)
            self.h_etaBChiPlus.Fill(genBChiPlus1.eta, genBChiPlus2.eta, genWeight)
            self.h_pTStopPlus.Fill(genStopPlus.pt,genWeight)
            self.h_etaStopPlus.Fill(genStopPlus.eta,genWeight)
            self.h_pTChiPlus.Fill(genChiPlus.pt,genWeight)
            self.h_etaChiPlus.Fill(genChiPlus.eta,genWeight)
          else:
            self.h_pTBStopMinus.Fill(genBStopMinus.pt,genWeight)
            self.h_etaBStopMinus.Fill(genBStopMinus.eta,genWeight)
            self.h_pTBChiMinus.Fill(genBChiMinus1.pt, genBChiMinus2.pt, genWeight)
            self.h_etaBChiMinus.Fill(genBChiMinus1.eta, genBChiMinus2.eta, genWeight)
            self.h_pTStopMinus.Fill(genStopMinus.pt,genWeight)
            self.h_etaStopMinus.Fill(genStopMinus.eta,genWeight)
            self.h_pTChiMinus.Fill(genChiMinus.pt,genWeight)
            self.h_etaChiMinus.Fill(genChiMinus.eta,genWeight)
          dEtaBChi = abs(genBStop.eta - genChi.eta)
          self.h_dEtaBChi.Fill(dEtaBChi,genWeight)
          self.h_dPhiBChi.Fill(abs(genBStop.p4().DeltaPhi(genChi.p4())),genWeight)
          self.h_dRBChi.Fill(abs(genBStop.p4().DeltaR(genChi.p4())),genWeight)
          dRChiMax = max(genBChi1.p4().DeltaR(genD.p4()),
                                   genBChi1.p4().DeltaR(genBChi2.p4()),
                                   genD.p4().DeltaR(genBChi2.p4()))
          self.h_dRChiMax.Fill(dRChiMax,genWeight)
          self.h_dRBB.Fill(genBChi1.p4().DeltaR(genBStop.p4()),genWeight)
          self.h_dEtaBB.Fill(abs(genBChi1.eta - genBStop.eta),genWeight)
          self.h_dPhiBB.Fill(abs(genBChi1.p4().DeltaPhi(genBStop.p4())),genWeight)
          self.h_pTBVsChi.Fill(genChi.pt,genBStop.pt,genWeight)
          self.h_dEtaVsPTStop.Fill(genChi.pt,abs(genChi.eta - genBStop.eta),genWeight)
          self.h_dEtaVsPTStopRatio.Fill(genChi.pt / (genStop.p4().M() - genChi.p4().M()),abs(genChi.eta - genBStop.eta),genWeight)
          self.h_passDijet.Fill(1 if (dRChiMax < 1.1 and dEtaBChi < 1.1) else 0,genWeight)

          #pT and eta of the gen AK4 jets
          for i,j in enumerate(genAK4Jets):
            if i == 0:
              self.h_pT1Gen.Fill(j.pt,genWeight)
              self.h_eta1Gen.Fill(j.eta,genWeight)
            if i == 1:
              self.h_pT2Gen.Fill(j.pt,genWeight)
              self.h_eta2Gen.Fill(j.eta,genWeight)
            if i == 2:
              self.h_pT3Gen.Fill(j.pt,genWeight)
              self.h_eta3Gen.Fill(j.eta,genWeight)
            if i == 3:
              self.h_pT4Gen.Fill(j.pt,genWeight)
              self.h_eta4Gen.Fill(j.eta,genWeight)

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
            self.h_mWJs.Fill((genWJ1 + genWJ2).M(),genWeight)
            dEtaWJ = abs(genWJ1.Eta() - genWJ2.Eta())
            self.h_dEtaWJs.Fill(dEtaWJ,genWeight)

        #-----------------------------------------------------------------------
        # RECO
        #-----------------------------------------------------------------------

        sumJet             = ROOT.TLorentzVector()
        sumJet4            = ROOT.TLorentzVector()
        sumJet3            = ROOT.TLorentzVector()
        sumJet3NoLead      = ROOT.TLorentzVector()
        sumJet3NoLeadOrSub = ROOT.TLorentzVector()

        # n jets
        self.h_nJets.Fill(len(jets),genWeight)
        self.h_ntLoose.Fill(len(looseTs),genWeight)
        self.h_ntMedium.Fill(len(mediumTs),genWeight)
        self.h_ntTight.Fill(len(tightTs),genWeight)
        self.h_nWLoose.Fill(len(looseWs),genWeight)
        self.h_nWMedium.Fill(len(mediumWs),genWeight)
        self.h_nWTight.Fill(len(tightWs),genWeight)
        self.h_nbLoose.Fill(len(looseBs),genWeight)
        self.h_nbMedium.Fill(len(mediumBs),genWeight)
        self.h_nbTight.Fill(len(tightBs),genWeight)
        self.h_ntDeepWP1.Fill(len(deepWP1Ts),genWeight)
        self.h_ntDeepWP2.Fill(len(deepWP2Ts),genWeight)
        self.h_ntDeepWP3.Fill(len(deepWP3Ts),genWeight)
        self.h_ntDeepWP4.Fill(len(deepWP4Ts),genWeight)
        self.h_nWDeepWP1.Fill(len(deepWP1Ws),genWeight)
        self.h_nWDeepWP2.Fill(len(deepWP2Ws),genWeight)
        self.h_nWDeepWP3.Fill(len(deepWP3Ws),genWeight)
        self.h_nWDeepWP4.Fill(len(deepWP4Ws),genWeight)

        # Event characteristics
        self.h_MET.Fill(MET.Mod(),genWeight)

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
          if i == 0:
           self.h_pT1.Fill(j.pt,genWeight)
           self.h_eta1.Fill(j.eta,genWeight)
          if i == 1:
           self.h_pT2.Fill(j.pt,genWeight)
           self.h_eta2.Fill(j.eta,genWeight)
          if i == 2:
           self.h_pT3.Fill(j.pt,genWeight)
           self.h_eta3.Fill(j.eta,genWeight)
          if i == 3:
           self.h_pT4.Fill(j.pt,genWeight)
           self.h_eta4.Fill(j.eta,genWeight)

        if len(jets) >= 5:
          self.h_m3NoLeadOrSub.Fill(sumJet3NoLeadOrSub.M(),genWeight)
          self.h_dEta15.Fill(abs(jets[0].eta - jets[4].eta),genWeight)
          self.h_dPhi15.Fill(abs(jets[0].p4().DeltaPhi(jets[4].p4())),genWeight)
          self.h_dR15.Fill(abs(jets[0].p4().DeltaR(jets[4].p4())),genWeight)
          self.h_dEta25.Fill(abs(jets[1].eta - jets[4].eta),genWeight)
          self.h_dPhi25.Fill(abs(jets[1].p4().DeltaPhi(jets[4].p4())),genWeight)
          self.h_dR25.Fill(abs(jets[1].p4().DeltaR(jets[4].p4())),genWeight)
          self.h_dEta35.Fill(abs(jets[2].eta - jets[4].eta),genWeight)
          self.h_dPhi35.Fill(abs(jets[2].p4().DeltaPhi(jets[4].p4())),genWeight)
          self.h_dR35.Fill(abs(jets[2].p4().DeltaR(jets[4].p4())),genWeight)
          self.h_dEta45.Fill(abs(jets[3].eta - jets[4].eta),genWeight)
          self.h_dPhi45.Fill(abs(jets[3].p4().DeltaPhi(jets[4].p4())),genWeight)
          self.h_dR45.Fill(abs(jets[3].p4().DeltaR(jets[4].p4())),genWeight)
          self.h_pT5Frac.Fill(jets[4].pt / HT,genWeight)
          self.h_pT5ETFracComp.Fill(jets[4].pt / sumJet3.Et(),genWeight)
          self.h_pT5ETFracUncomp.Fill(jets[4].pt / sumJet3NoLead.Et(),genWeight)
          self.h_pT5FracChiComp.Fill(jets[4].pt / HT3,genWeight)
          self.h_pT5FracChiUncomp.Fill(jets[4].pt / HT3NoLead,genWeight)
        if len(jets) >= 4:
          self.h_m4.Fill(sumJet4.M(),genWeight)
          self.h_m3NoLead.Fill(sumJet3NoLead.M(),genWeight)
          self.h_m3Vsm4.Fill(sumJet4.M(),sumJet3.M(),genWeight)
          self.h_m3NoLeadVsm4.Fill(sumJet4.M(),sumJet3NoLead.M(),genWeight)
          self.h_dEta14.Fill(abs(jets[0].eta - jets[3].eta),genWeight)
          self.h_dPhi14.Fill(abs(jets[0].p4().DeltaPhi(jets[3].p4())),genWeight)
          self.h_dR14.Fill(abs(jets[0].p4().DeltaR(jets[3].p4())),genWeight)
          self.h_dEta24.Fill(abs(jets[1].eta - jets[3].eta),genWeight)
          self.h_dPhi24.Fill(abs(jets[1].p4().DeltaPhi(jets[3].p4())),genWeight)
          self.h_dR24.Fill(abs(jets[1].p4().DeltaR(jets[3].p4())),genWeight)
          self.h_dEta34.Fill(abs(jets[2].eta - jets[3].eta),genWeight)
          self.h_dPhi34.Fill(abs(jets[2].p4().DeltaPhi(jets[3].p4())),genWeight)
          self.h_dR34.Fill(abs(jets[2].p4().DeltaR(jets[3].p4())),genWeight)
          self.h_pT4Frac.Fill(jets[3].pt / HT,genWeight)
          self.h_pT4ETFracComp.Fill(jets[3].pt / sumJet3.Et(),genWeight)
          self.h_pT4ETFracUncomp.Fill(jets[3].pt / sumJet3NoLead.Et(),genWeight)
          self.h_pT4FracChiComp.Fill(jets[3].pt / HT3,genWeight)
          self.h_pT4FracChiUncomp.Fill(jets[3].pt / HT3NoLead,genWeight)
          pTMeanComp   = mean([jets[0].pt,jets[1].pt,jets[2].pt])
          pTSDComp     = std([jets[0].pt,jets[1].pt,jets[2].pt])
          pTMeanUncomp = mean([jets[1].pt,jets[2].pt,jets[3].pt])
          pTSDUncomp   = std([jets[1].pt,jets[2].pt,jets[3].pt])
          self.h_pTMeanComp.Fill(pTMeanComp,genWeight)
          self.h_pTSDComp.Fill(pTSDComp,genWeight)
          self.h_pTSDMeanFracComp.Fill(pTSDComp / pTMeanComp,genWeight)
          self.h_pTMeanUncomp.Fill(pTMeanUncomp,genWeight)
          self.h_pTSDUncomp.Fill(pTSDUncomp,genWeight)
          self.h_pTSDMeanFracUncomp.Fill(pTSDUncomp / pTMeanUncomp,genWeight)
        if len(jets) >= 3:
          self.h_m3.Fill(sumJet3.M(),genWeight)
          self.h_dEta13.Fill(abs(jets[0].eta - jets[2].eta),genWeight)
          self.h_dPhi13.Fill(abs(jets[0].p4().DeltaPhi(jets[2].p4())),genWeight)
          self.h_dR13.Fill(abs(jets[0].p4().DeltaR(jets[2].p4())),genWeight)
          self.h_dEta23.Fill(abs(jets[1].eta - jets[2].eta),genWeight)
          self.h_dPhi23.Fill(abs(jets[1].p4().DeltaPhi(jets[2].p4())),genWeight)
          self.h_dR23.Fill(abs(jets[1].p4().DeltaR(jets[2].p4())),genWeight)
          self.h_pT3Frac.Fill(jets[2].pt / HT,genWeight)
          self.h_pT3ETFracComp.Fill(jets[2].pt / sumJet3.Et(),genWeight)
          self.h_pT3ETFracUncomp.Fill(jets[2].pt / sumJet3NoLead.Et(),genWeight)
          self.h_pT3FracChiComp.Fill(jets[2].pt / HT3,genWeight)
          self.h_pT3FracChiUncomp.Fill(jets[2].pt / HT3NoLead,genWeight)
        if len(jets) >= 2:
          self.h_dEta12.Fill(abs(jets[0].eta - jets[1].eta),genWeight)
          self.h_dPhi12.Fill(abs(jets[0].p4().DeltaPhi(jets[1].p4())),genWeight)
          self.h_dR12.Fill(abs(jets[0].p4().DeltaR(jets[1].p4())),genWeight)
          self.h_pT2Frac.Fill(jets[1].pt / HT,genWeight)
          self.h_pT2ETFracComp.Fill(jets[1].pt / sumJet3.Et(),genWeight)
          self.h_pT2ETFracUncomp.Fill(jets[1].pt / sumJet3NoLead.Et(),genWeight)
          self.h_pT2FracChiComp.Fill(jets[1].pt / HT3,genWeight)
          self.h_pT2FracChiUncomp.Fill(jets[1].pt / HT3NoLead,genWeight)
        if len(jets) >= 1:
          self.h_HT.Fill(HT,genWeight)
          self.h_mAll.Fill(sumJet.M(),genWeight)
          self.h_pT1Frac.Fill(jets[0].pt / HT,genWeight)
          self.h_pT1ETFracComp.Fill(jets[0].pt / sumJet3.Et(),genWeight)
          self.h_pT1ETFracUncomp.Fill(jets[0].pt / sumJet3NoLead.Et(),genWeight)
          self.h_pT1FracChiComp.Fill(jets[0].pt / HT3,genWeight)
          self.h_pT1FracChiUncomp.Fill(jets[0].pt / HT3NoLead,genWeight)

        jetCombos2D = [[''.join(map(str,i)),''.join(map(str,j))] for i,j in combinations(combinations([1,2,3,4],2),2)]
        self.fill2DHists(['dPhi{}Vs{}'.format(i,j) for [i,j] in jetCombos2D],
                          [abs(jets[int(i)-1].p4().DeltaPhi(jets[int(j)-1].p4())) for i,j in [x[1] for x in jetCombos2D]],
                          [abs(jets[int(i)-1].p4().DeltaPhi(jets[int(j)-1].p4())) for i,j in [x[0] for x in jetCombos2D]],
                          genWeight)

        # b jets
        for i,b in enumerate(looseBs):
          if i == 0: self.h_pTb1.Fill(b.pt,genWeight)
          elif i == 1: self.h_pTb2.Fill(b.pt,genWeight)
          elif i == 2: self.h_pTb3.Fill(b.pt,genWeight)
          elif i == 3: self.h_pTb4.Fill(b.pt,genWeight)

        if len(looseBs) >= 2: 
          self.h_dEtabb.Fill(abs(looseBs[0].eta - looseBs[1].eta),genWeight)
          self.h_dPhibb.Fill(abs(looseBs[0].p4().DeltaPhi(looseBs[1].p4())),genWeight)
          self.h_dRbb.Fill(abs(looseBs[0].p4().DeltaR(looseBs[1].p4())),genWeight)

        self.h_dEtaRecoComp.Fill(abs(sumJet3.Eta() - sumJet4.Eta()),genWeight)
        self.h_dPhiRecoComp.Fill(abs(sumJet3.DeltaPhi(sumJet4)),genWeight)
        self.h_dRRecoComp.Fill(abs(sumJet3.DeltaR(sumJet4)),genWeight)
        self.h_dEtaRecoUncomp.Fill(abs(sumJet3NoLead.Eta() - sumJet4.Eta()),genWeight)
        self.h_dPhiRecoUncomp.Fill(abs(sumJet3NoLead.DeltaPhi(sumJet4)),genWeight)
        self.h_dRRecoUncomp.Fill(abs(sumJet3NoLead.DeltaR(sumJet4)),genWeight)

        return True

parser = argparse.ArgumentParser(description='Single Stop Analyzer')
parser.add_argument('--sample',type=str,default='signal',choices=['signal','TT','TT2018','QCD','QCD2018','ZQQ2018','ST2018','WQQ2018','ZNuNu2018','Diboson2018'],help='Sample to run over')
parser.add_argument('--tag',type=str,default='test',help='Tag for output label')
parser.add_argument('-n',type=int,default=1,help='Sample index to run over for backgrounds')
parser.add_argument('--points',type=str,default='all',help='Signal point(s) to run over, comma separated in MSTOP_MCHI format; "all" to run over all available points')
parser.add_argument('--useskim',action='store_true',default=False,help='Flag to use NANOAODs skimmed with the nominal selections')
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
    files = glob.glob('/eos/uscms/store/user/dmahon/condor/RPVSingleStopMC313/NANOAOD-ALL/NANOAOD-{}.root'.format(masses))
    #files = glob.glob('/eos/uscms/store/user/dmahon/condor/RPVSingleStopMC/NANOAOD/NANOAOD-{}-*.root'.format(masses))
    files = ['root://cmsxrootd.fnal.gov/' + x.replace('/eos/uscms','') for x in files]
    #files = ['file:/uscms_data/d3/dmahon/RPVSingleStopRun3Patched/NANOAOD/CMSSW_12_4_5/test_2000_100-1.root']
    #files = ['/uscms_data/d3/dmahon/RPVSingleStopRun3Patched/NANOAOD/files/NANOAOD-{}.root'.format(masses)]
    p = PostProcessor(".", files, cut=preselection, branchsel=None,
                      modules=[ExampleAnalysis(isSignal=1,MCCampaign='UL2018',isSkimmed=False)],
                      noOut=True, histFileName='{}/{}_{}.root'.format(outputPath,args.sample,masses), histDirName="plots",
                      maxEntries=None)
    p.run()

elif args.useskim: 

  print('Running over skimmed {} files'.format(args.sample))
  print('Using UL 2018 MC campaign working points')

  files = ['root://cmsxrootd.fnal.gov//store/user/ckapsiak/SingleStop/Skims/Skim_2023_23_03/{}.root'.format(args.sample)]
  if len(files) != 1: print('WARNING: Multiple files selected. All must be from the same MC campaign.')
  p = PostProcessor(".", files, cut='', branchsel=None,
                    modules=[ExampleAnalysis(isSignal=0,MCCampaign='UL2018',isSkimmed=True)],
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
                    modules=[ExampleAnalysis(isSignal=0,MCCampaign=MCCampaign,isSkimmed=False)], 
                    noOut=True, histFileName='{}/{}-{}.root'.format(outputPath,args.sample,args.n), histDirName="plots",
                    maxEntries=None)
  p.run() 
