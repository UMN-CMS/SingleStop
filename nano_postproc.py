#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from importlib import import_module
import os
import sys
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

class ExampleAnalysis(Module):

    def __init__(self,isSignal,MCCampaign):
        self.writeHistFile = True
        self.isSignal = isSignal
        self.MCCampaign = MCCampaign

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

    def analyze(self, event):

        # Set b tagging WPs
        if   self.MCCampaign == 'UL2016preVFP':         bTagWPs = [0.0508,0.2598,0.6502]
        elif self.MCCampaign == 'UL2016postVFP':        bTagWPs = [0.0480,0.2489,0.6377]
        elif self.MCCampaign == 'UL2017':               bTagWPs = [0.0532,0.3040,0.7476]
        elif self.MCCampaign == 'UL2018':               bTagWPs = [0.0490,0.2783,0.7100]

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
        genParts       = list(Collection(event,"GenPart"))
        genAK4Jets     = list(Collection(event,"GenJet"))
        goodElectrons  = filter(lambda x: x.cutBased == 4 and x.miniPFRelIso_all < 0.1 and x.pt > 30 and abs(x.eta) < 2.4,list(Collection(event,"Electron")))
        goodMuons      = filter(lambda x: x.mediumId and x.miniPFRelIso_all < 0.2 and x.pt > 30 and abs(x.eta) < 2.4,list(Collection(event,"Muon")))
        MET            = ROOT.TVector2()
        MET.SetMagPhi(event.MET_pt,event.MET_phi)

        # Cuts
        if len(jets) < 4 or len(jets) > 5: return False
        if not jets[0].pt > 300: return False
        if not (event.HLT_PFHT1050 or event.HLT_AK8PFJet360_TrimMass30): return False
        if len(goodElectrons) != 0: return False
        if len(goodMuons) != 0 : return False
        if len(looseBs) < 2: return False
        if not 2 < abs(jets[0].p4().DeltaR(jets[1].p4())) < 4: return False
        #if len(tightTs) != 0: return False

        return True

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] outputDir inputFiles")
    parser.add_option('-n',type=int,default=1,help='Sample index to run over for backgrounds')
    parser.add_option("-s", "--postfix", dest="postfix", type="string", default=None,
                      help="Postfix which will be appended to the file name (default: _Friend for friends, _Skim for skims)")
    parser.add_option("-J", "--json", dest="json", type="string",
                      default=None, help="Select events using this JSON file")
    parser.add_option("-c", "--cut", dest="cut", type="string",
                      default=None, help="Cut string")
    parser.add_option("-b", "--branch-selection", dest="branchsel",
                      type="string", default=None, help="Branch selection")
    parser.add_option("--bi", "--branch-selection-input", dest="branchsel_in",
                      type="string", default=None, help="Branch selection input")
    parser.add_option("--bo", "--branch-selection-output", dest="branchsel_out",
                      type="string", default=None, help="Branch selection output")
    parser.add_option("--friend", dest="friend", action="store_true", default=False,
                      help="Produce friend trees in output (current default is to produce full trees)")
    parser.add_option("--full", dest="friend", action="store_false", default=False,
                      help="Produce full trees in output (this is the current default)")
    parser.add_option("--noout", dest="noOut", action="store_true",
                      default=False, help="Do not produce output, just run modules")
    parser.add_option("-P", "--prefetch", dest="prefetch", action="store_true", default=False,
                      help="Prefetch input files locally instead of accessing them via xrootd")
    parser.add_option("--long-term-cache", dest="longTermCache", action="store_true", default=False,
                      help="Keep prefetched files across runs instead of deleting them at the end")
    parser.add_option("-N", "--max-entries", dest="maxEntries", type="long", default=None,
                      help="Maximum number of entries to process from any single given input tree")
    parser.add_option("--first-entry", dest="firstEntry", type="long", default=0,
                      help="First entry to process in the three (to be used together with --max-entries)")
    parser.add_option("--justcount", dest="justcount", default=False,
                      action="store_true", help="Just report the number of selected events")
    parser.add_option("-I", "--import", dest="imports", type="string", default=[], action="append",
                      nargs=2, help="Import modules (python package, comma-separated list of ")
    parser.add_option("-z", "--compression", dest="compression", type="string",
                      default=("LZMA:9"), help="Compression: none, or (algo):(level) ")
    parser.add_option("--output",dest="output",type="string",default="output",help="Output directory name")
    parser.add_option('--sample',type=str,default='QCD2018',help='Sample to run over')

    (options,args) = parser.parse_args()
    args = options

    if args.friend:
        if args.cut or args.json:
            raise RuntimeError(
                "Can't apply JSON or cut selection when producing friends")

    outdir = args.output
   
    if   args.sample == 'TT':          sampleFile = 'TTToHadronic.txt'
    elif args.sample == 'TT2018':      sampleFile = 'TTToHadronic2018.txt'
    elif args.sample == 'QCD':         sampleFile = 'QCDBEnriched.txt'
    elif args.sample == 'QCD2018':     sampleFile = 'QCDBEnriched2018.txt'
    elif args.sample == 'ZQQ2018':     sampleFile = 'ZJetsToQQ2018.txt'
    elif args.sample == 'ST2018':      sampleFile = 'STHadronic2018.txt'
    elif args.sample == 'WQQ2018':     sampleFile = 'WJetsToQQ2018.txt'
    elif args.sample == 'ZNuNu2018':   sampleFile = 'ZJetsToNuNu2018.txt'
    elif args.sample == 'Diboson2018': sampleFile = 'Diboson2018.txt'
    elif args.sample == 'signal': sampleFile = 'RPV.txt'
    elif args.sample != 'signal': print('ERROR: Unexpected sample argument')

    files = open('samples/{}'.format(sampleFile)).read().split('\n') 
    #files = ['root://cmsxrootd.fnal.gov/' + x.replace('/eos/uscms','') for x in files]
    files = [['root://cmsxrootd.fnal.gov/' + x.replace('/eos/uscms','') for x in files][:-1][args.n - 1]]

    modules = []
    for mod, names in args.imports:
        import_module(mod)
        obj = sys.modules[mod]
        selnames = names.split(",")
        mods = dir(obj)
        for name in selnames:
            if name in mods:
                print("Loading %s from %s " % (name, mod))
                if type(getattr(obj, name)) == list:
                    for mod in getattr(obj, name):
                        modules.append(mod())
                else:
                    modules.append(getattr(obj, name)())
    if args.noOut:
        if len(modules) == 0:
            raise RuntimeError(
                "Running with --noout and no modules does nothing!")
    if args.branchsel != None:
        args.branchsel_in = args.branchsel
        args.branchsel_out = args.branchsel
    if 'UL16' in files[0]:
      if 'preVFP' in files[0]:    MCCampaign = 'UL2016preVFP'
      else:                      MCCampaign = 'UL2016postVFP'
    elif 'UL17' in files[0]:      MCCampaign = 'UL2017'
    elif 'UL18' in files[0]:      MCCampaign = 'UL2018'
    else:
      print('ERROR: Unable to determine MC campaign of {}'.format(files[0]))
      sys.exit()
    modules.append(ExampleAnalysis(isSignal=0,MCCampaign=MCCampaign))
    outputPath = 'output/{}/'.format(args.sample)
    print(outputPath)
    preselection = ( '(Jet_pt[3] > 30) &&'
                    '(Jet_pt[0] > 300) &&'
                    '(HLT_PFHT1050 || HLT_AK8PFJet360_TrimMass30)'
                   )
    p = PostProcessor(outputPath, files,
                      cut=preselection,
                      branchsel=args.branchsel_in,
                      modules=modules,
                      compression=args.compression,
                      friend=args.friend,
                      postfix=args.postfix,
                      jsonInput=args.json,
                      noOut=args.noOut,
                      justcount=args.justcount,
                      prefetch=args.prefetch,
                      longTermCache=args.longTermCache,
                      maxEntries=args.maxEntries,
                      firstEntry=args.firstEntry,
                      outputbranchsel=args.branchsel_out)
    p.run()
