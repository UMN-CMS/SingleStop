#!/usr/bin/env python

import os, sys, ROOT, argparse, glob
import numpy
from array import array
from functools import wraps
import itertools as it
import re


def zeroFill(tree, brName, brObj, allowNonBool=False):
    # typename: (numpy type code, root type code)
    branch_type_dict = {
        "Bool_t": ("?", "O"),
        "Float_t": ("f4", "F"),
        "UInt_t": ("u4", "i"),
        "Long64_t": ("i8", "L"),
        "Double_t": ("f8", "D"),
    }
    brType = brObj.GetLeaf(brName).GetTypeName()
    if (not allowNonBool) and (brType != "Bool_t"):
        print(
            (
                "Did not expect to back fill non-boolean branches",
                tree,
                brName,
                brObj.GetLeaf(br).GetTypeName(),
            )
        )
    else:
        if brType not in branch_type_dict:
            raise RuntimeError("Impossible to backfill branch of type %s" % brType)
        buff = numpy.zeros(1, dtype=numpy.dtype(branch_type_dict[brType][0]))
        b = tree.Branch(brName, buff, brName + "/" + branch_type_dict[brType][1])
        # be sure we do not trigger flushing
        b.SetBasketSize(tree.GetEntries() * 2)
        for x in range(0, tree.GetEntries()):
            b.Fill()
        b.ResetAddress()


def addBranchToTree(tree, name, val):
    buff = array("f", [val])
    b = tree.Branch(name, buff, name + "/f")
    b.SetBasketSize(tree.GetEntries() * 100)
    for x in range(0, tree.GetEntries()):
        b.Fill()
    b.ResetAddress()


def haddNano(ofname, files):
    fileHandles = []
    goFast = True
    for fn in files:
        print("Adding file" + str(fn))
        fileHandles.append((fn, ROOT.TFile.Open(fn)))
        if (
            fileHandles[-1][1].GetCompressionSettings()
            != fileHandles[0][1].GetCompressionSettings()
        ):
            goFast = False
            print("Disabling fast merging as inputs have different compressions")
    of = ROOT.TFile(ofname, "recreate")
    if goFast:
        of.SetCompressionSettings(fileHandles[0][1].GetCompressionSettings())
    of.cd()

    for e in fileHandles[0][1].GetListOfKeys():
        print(e)
        scalevar = None
        name = e.GetName()
        obj = e.ReadObj()
        cl = ROOT.TClass.GetClass(e.GetClassName())
        inputs = ROOT.TList()
        isTree = obj.IsA().InheritsFrom(ROOT.TTree.Class())
        scalevar = None
        if isTree:
            obj = obj.CloneTree(-1, "fast" if goFast else "")
            branchNames = set([x.GetName() for x in obj.GetListOfBranches()])

        for fht in fileHandles[1:]:
            fn = fht[0]
            fh = fht[1]
            print(fh)
            otherObj = fh.GetListOfKeys().FindObject(name).ReadObj()
            inputs.Add(otherObj)
            if isTree and obj.GetName() == "Events":
                otherObj.SetAutoFlush(0)
                otherBranches = set([x.GetName() for x in otherObj.GetListOfBranches()])
                missingBranches = list(branchNames - otherBranches)
                additionalBranches = list(otherBranches - branchNames)
                for br in missingBranches:
                    # fill "Other"
                    zeroFill(otherObj, br, obj.GetListOfBranches().FindObject(br))
                for br in additionalBranches:
                    # fill main
                    branchNames.add(br)
                    zeroFill(obj, br, otherObj.GetListOfBranches().FindObject(br))
                # merge immediately for trees
            if isTree and obj.GetName() == "Runs":
                otherObj.SetAutoFlush(0)
                otherBranches = set([x.GetName() for x in otherObj.GetListOfBranches()])
                missingBranches = list(branchNames - otherBranches)
                additionalBranches = list(otherBranches - branchNames)
                for br in missingBranches:
                    # fill "Other"
                    zeroFill(
                        otherObj,
                        br,
                        obj.GetListOfBranches().FindObject(br),
                        allowNonBool=True,
                    )
                for br in additionalBranches:
                    # fill main
                    branchNames.add(br)
                    zeroFill(
                        obj,
                        br,
                        otherObj.GetListOfBranches().FindObject(br),
                        allowNonBool=True,
                    )
                # merge immediately for trees
            if isTree:
                obj.Merge(inputs, "fast" if goFast else "")
                inputs.Clear()

        if isTree:
            obj.Write()
        elif obj.IsA().InheritsFrom(ROOT.TH1.Class()):
            obj.Merge(inputs)
            obj.Write()
        elif obj.IsA().InheritsFrom(ROOT.TObjString.Class()):
            for st in inputs:
                if st.GetString() != obj.GetString():
                    print("Strings are not matching")
            obj.Write()
        else:
            print("Cannot handle " + str(obj.IsA().GetName()))


def associateFiles(sample_file, input_dir):
    fnames = []
    with open(sample_file, "r") as f:
        fnames = [line.strip() for line in f]
    g = glob.glob("{}/*.root".format(input_dir))
    fmap = {
        rootfile: (
            next(f for f in fnames if os.path.split(rootfile)[1].split("_")[0] in f)
        )
        for rootfile in g
    }
    return fmap


def qcdMatcher(x):
    m = re.search(r'(RunIISummer.+NanoAODv9).+(QCD_HT\d+to(?:\d+|Inf))', x)
    return "{}_{}".format(m.group(1), m.group(2))

def makeGroups(fmap, matcher, xsec_getter):
    pairs = list(fmap.items())
    pairs = sorted(pairs, key=lambda x: x[1])
    groups = it.groupby(pairs, lambda x: matcher(x[1]))
    return {x:(xsec_getter(x), list(y[0] for y in z)) for x,z in groups}

def getQCDXSec(fname):
    xsecs = { (1,50)   : 186100000.0 * 1E3, 
        (1,100)  : 23630000.0  * 1E3, 
        (1,200)  : 1554000.0   * 1E3, 
        (1,300)  : 323800.0    * 1E3, 
        (1,500)  : 30280.0     * 1E3, 
        (1,700)  : 6392.0      * 1E3, 
        (1,1000) : 1118.0      * 1E3, 
        (1,1500) : 108.9       * 1E3, 
        (1,2000) : 21.93       * 1E3, 
        (3,50)   : 187300000.0 * 1E3, 
        (3,100)  : 23590000.0  * 1E3, 
        (3,200)  : 1555000.0   * 1E3, 
        (3,300)  : 324500.0    * 1E3, 
        (3,500)  : 30310.0     * 1E3, 
        (3,700)  : 6444.0      * 1E3, 
        (3,1000) : 1127.0      * 1E3, 
        (3,1500) : 109.8       * 1E3, 
        (3,2000) : 21.98       * 1E3 
    }
    period = -1 
    HT = -1
    if   'preVFP' in fname: period = 0
    elif 'UL16' in fname:   period = 1
    elif 'UL17' in fname:   period = 2
    elif 'UL18' in fname:   period = 3
    else: print('ERROR: Data period could not be determined for {}'.format(fname))

    if   '50to100' in fname:    HT = 50	
    elif '100to200' in fname:   HT = 100
    elif '200to300' in fname:   HT = 200    
    elif '300to500' in fname:   HT = 300    
    elif '500to700' in fname:   HT = 500    
    elif '700to1000' in fname:  HT = 700   
    elif '1000to1500' in fname: HT = 1000  
    elif '1500to2000' in fname: HT = 1500  
    elif '2000toInf' in fname:  HT = 2000
    else: print('ERROR: HT range could not be determined for {}'.format(fname))
    
    if period == -1 or HT == -1: sys.exit()
    else: x = xsecs[(period,HT)]
    return x
#
#
#

def ttMatcher(x):
    return "TT"

def getTTXSec(fname):
    #lumiTarget = 59.8
    #lumiSample = 331506194 / (831.8 * 0.457) * 1e-3
    #SF = lumiTarget / lumiSample
    return (831.8 * 0.457) * 1e-3

def getSignalXSec(fname):
    points = ['1000_400','1000_600','1000_900',
              '1200_400','1200_600','1200_1100',
              '1300_400','1300_600','1300_1200',
              '1400_400','1400_600','1400_1300',
              '1500_400','1500_600','1500_900','1500_1400',
              '2000_400','2000_600','2000_900','2000_1400','2000_1900'] 
    lambdapp = 0.1
    SFs = {
      1000:48000 * lambdapp**2,
      1200:21000 * lambdapp**2,
      1300:15000 * lambdapp**2,
      1400:10000 * lambdapp**2,
      1500:7300  * lambdapp**2,
      2000:1600  * lambdapp**2,
    }
    part = fname.split("-")[1]
    mStop = int(part.split("_")[0])
    return SFs[mStop]



def zqqMatcher(x):
    m = re.search(r'(RunIISummer.+NanoAODv9)/(.+\d+to(?:\d+|Inf))', x)
    return "{}_{}".format(m.group(1), m.group(2))

def getZQQ(fname):
    SFs = {
        # 2018
        (3, 200):  1012.0 * 1e3 ,
        (3, 400):  114.2 * 1e3 ,
        (3, 600):  25.34 * 1e3 ,
        (3, 800):  12.99 * 1e3 
    }
    if "preVFP" in fname:
        period = 0
    elif "UL16" in fname:
        period = 1
    elif "UL17" in fname:
        period = 2
    elif "UL18" in fname:
        period = 3
    else:
        print("ERROR: Data period could not be determined for {}".format(fname))

    if "200to400" in fname:
        HT = 200
    elif "400to600" in fname:
        HT = 400
    elif "600to800" in fname:
        HT = 600
    elif "800toInf" in fname:
        HT = 800
    else:
        print("ERROR: HT range could not be determined for {}".format(fname))

    return SFs[(period, HT)]

def zNuNuMatcher(x):
    m = re.search(r'(RunIISummer.+NanoAODv9)/(.+\d+to(?:\d+|Inf))', x)
    return "{}_{}".format(m.group(1), m.group(2))

def getZNuNu(fname):
    SFs = {
        # 2018
        (3, 100): 267.0 * 1.1347 * 1e3 ,
        (3, 200): 73.08 * 1.1347 * 1e3 ,
        (3, 400): 9.921 * 1.1347 * 1e3 ,
        (3, 600): 2.409 * 1.1347 * 1e3 ,
        (3, 800): 1.078 * 1.1347 * 1e3 ,
        (3, 1200): 0.2514 * 1.1347 * 1e3 ,
        (3, 2500): 0.005614 * 1.1347 * 1e3 ,
    }

    if "preVFP" in fname:
        period = 0
    elif "UL16" in fname:
        period = 1
    elif "UL17" in fname:
        period = 2
    elif "UL18" in fname:
        period = 3
    else:
        print("ERROR: Data period could not be determined for {}".format(fname))

    if "100To200" in fname:
        HT = 100
    elif "200To400" in fname:
        HT = 200
    elif "400To600" in fname:
        HT = 400
    elif "600To800" in fname:
        HT = 600
    elif "800To1200" in fname:
        HT = 800
    elif "1200To2500" in fname:
        HT = 1200
    elif "2500ToInf" in fname:
        HT = 2500
    else:
        print("ERROR: HT range could not be determined for {}".format(fname))

    SF = SFs[(period, HT)]
    return SF
#
#

def wqqMatcher(x):
    m = re.search(r'(RunIISummer.+NanoAODv9)/(.+\d+to(?:\d+|Inf))', x)

    return "{}_{}".format(m.group(1), m.group(2))
def getWQQ(fname):
    SFs = {
        # 2018
        (3, 200): 59.8 * 2549.0 * 1e3 ,
        (3, 400): 59.8 * 276.5 * 1e3 ,
        (3, 600): 59.8 * 59.25 * 1e3 ,
        (3, 800): 59.8 * 28.75 * 1e3 ,
    }

    if "preVFP" in fname:
        period = 0
    elif "UL16" in fname:
        period = 1
    elif "UL17" in fname:
        period = 2
    elif "UL18" in fname:
        period = 3
    else:
        print("ERROR: Data period could not be determined for {}".format(fname))

    if "200to400" in fname:
        HT = 200
    elif "400to600" in fname:
        HT = 400
    elif "600to800" in fname:
        HT = 600
    elif "800toInf" in fname:
        HT = 800
    else:
        print("ERROR: HT range could not be determined for {}".format(fname))

    SF = SFs[(period, HT)]
    return SF

def dibosonMatcher(x):
    m = re.search(r'(RunIISummer.+NanoAODv9)/(.+13TeV)', x)
    return "{}_{}".format(m.group(1), m.group(2))

def getDiboson(fname):
    SFs = {
        # 2018
        (3, "WW"): 59.8 * 118.7 * 1e3 ,
        (3, "WZ"): 59.8 * 47.13 * 1e3 ,
        (3, "ZZ"): 59.8 * 16.523 * 1e3 ,
    }

    if "preVFP" in fname:
        period = 0
    elif "UL16" in fname:
        period = 1
    elif "UL17" in fname:
        period = 2
    elif "UL18" in fname:
        period = 3
    else:
        print("ERROR: Data period could not be determined for {}".format(fname))

    if "WW" in fname:
        mode = "WW"
    elif "WZ" in fname:
        mode = "WZ"
    elif "ZZ" in fname:
        mode = "ZZ"
    else:
        print("ERROR: HT range could not be determined for {}".format(fname))

    SF = SFs[(period, mode)]
    return SF
#
#
def stMatcher(x):
    m = re.search(r'(RunIISummer.+NanoAODv9)/(.+13TeV)', x)
    return "{}_{}".format(m.group(1), m.group(2))

def getST(fname):
    SFs = {
        # 2018
        (3, "s-channel"): 59.8 * 11.03 * 0.457 * 1e3 ,
        (3, "t-channel_antitop"): 59.8 * 80.95 * 1e3 ,
        (3, "t-channel_top"): 59.8 * 136.02 * 1e3 ,
        (3, "tW_antitop"): 59.8 * 35.85 * 1e3 ,
        (3, "tW_top"): 59.8 * 35.85 * 1e3 ,
    }
    if "preVFP" in fname:
        period = 0
    elif "UL16" in fname:
        period = 1
    elif "UL17" in fname:
        period = 2
    elif "UL18" in fname:
        period = 3
    else:
        print("ERROR: Data period could not be determined for {}".format(fname))

    if "s-channel" in fname:
        mode = "s-channel"
    elif "t-channel_antitop" in fname:
        mode = "t-channel_antitop"
    elif "t-channel_top" in fname:
        mode = "t-channel_top"
    elif "tW_antitop" in fname:
        mode = "tW_antitop"
    elif "tW_top" in fname:
        mode = "tW_top"
    else:
        print("ERROR: HT range could not be determined for {}".format(fname))

    if period == -1 or mode == -1:
        sys.exit()
    SF = SFs[(period, mode)]
    return SF
funcs = {
    #"QCD2018": scaleQCDBEnchriched,
    #"QCDBEnriched2018": scaleQCDBEnchriched,
    "QCDInclusive2018": (qcdMatcher, getQCDXSec),
    "TT2018": (ttMatcher, getTTXSec),
    "ZQQ2018": (zqqMatcher, getZQQ),
    "WQQ2018": (wqqMatcher, getWQQ),
    "WQQ2018": (wqqMatcher, getWQQ),
    #"ZJetsToQQ2018": scaleZQQ,
    #"WQQ2018": scaleWQQ,
    #"WJetsToQQ": scaleWQQ,
    #"WJetsToQQ2018": scaleWQQ,
    "Diboson2018": (dibosonMatcher, getDiboson),
    "ST2018": (stMatcher, getST),
    "ZNuNu2018": (zNuNuMatcher, getZNuNu),
    #"signal": scaleSignal,
}


def haddGroups(outdir, groups):
    for group,files in groups.items():
        out = os.path.join(outdir, group) + ".root"
        print("Creating group {} at {}".format(group, out))
        haddNano(out , files)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MC Scaler")
    parser.add_argument("--sample", type=str, default=None)
    parser.add_argument(
        "--input-dir", type=str, required=True, help="Path to input files"
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Path to output scaled file"
    )
    parser.add_argument("--sample-dir", type=str, default="samples")

    args = parser.parse_args()
    sample = args.sample
    input_dir = args.input_dir
    output = args.output

    is_signal = False

    if args.sample == "TT":
        sample_file = "TTToHadronic.txt"
    elif args.sample == "TT2018":
        sample_file = "TTToHadronic2018.txt"
    elif args.sample == "QCD2018":
        sample_file = "QCDBEnriched2018.txt"
    elif args.sample == "QCDInclusive2018":
        sample_file = "QCDInclusive2018.txt"
    elif args.sample == "ZQQ2018":
        sample_file = "ZJetsToQQ2018.txt"
    elif args.sample == "ST2018":
        sample_file = "STHadronic2018.txt"
    elif args.sample == "WQQ2018":
        sample_file = "WJetsToQQ2018.txt"
    elif args.sample == "ZNuNu2018":
        sample_file = "ZJetsToNuNu2018.txt"
    elif args.sample == "Diboson2018":
        sample_file = "Diboson2018.txt"
    elif args.sample == "signal":
        is_signal = True
        sample_file = "RPV.txt"
    else:
        raise ValueError("Cannot scale {}".format(args.sample))

    sample_file = "{}/{}".format(args.sample_dir, sample_file)

    try:
        os.remove(output)
    except OSError:
        pass

    if is_signal:
        raw_files = glob.glob("{}/*.root".format(input_dir))
        for f in raw_files:
            p = os.path.split(f)[-1].split('-')[1]
            #haddNano(os.path.join(output, "signal_{}".format(p)), [f], sf_func)
    else:
        afiles = associateFiles(sample_file, input_dir)
        #print(afiles)
        raw_files = [f for f in afiles]
        matcher, xsec = funcs[sample]
        groups = makeGroups(afiles, matcher, xsec)
        print(list(groups))
        haddGroups(output, {x:y[1] for x,y in groups.items()})

        #haddNano(output, raw_files, sf_func)
