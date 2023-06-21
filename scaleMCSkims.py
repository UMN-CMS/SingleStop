#!/usr/bin/env python

import os, sys, ROOT, argparse, glob
import numpy
from array import array
from functools import wraps


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


def haddNano(ofname, files, scale_func):
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
        scalevar = None
        name = e.GetName()
        obj = e.ReadObj()
        cl = ROOT.TClass.GetClass(e.GetClassName())
        inputs = ROOT.TList()
        isTree = obj.IsA().InheritsFrom(ROOT.TTree.Class())
        scalevar = None
        if isTree:
            obj = obj.CloneTree(-1, "fast" if goFast else "")
            if scale_func and obj.GetName() == "Events":
                addBranchToTree(obj, "MCScaleWeight", scale_func(fileHandles[0][0]))
            branchNames = set([x.GetName() for x in obj.GetListOfBranches()])

        for fht in fileHandles[1:]:
            fn = fht[0]
            fh = fht[1]
            print("Adding file {}".format(fh))
            otherObj = fh.GetListOfKeys().FindObject(name).ReadObj()
            inputs.Add(otherObj)
            if isTree and obj.GetName() == "Events":
                if scale_func:
                    addBranchToTree(otherObj, "MCScaleWeight", scale_func(fn))
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


def scaled_lumi(lumi_target, lumi_sample):
    def decorator(func):
        @wraps(func)
        def retFunc(fname):
            x = func(fname) * lumi_target / lumi_sample
            return x

        return retFunc

    return decorator


def translate_files(mapping):
    def decorator(func):
        @wraps(func)
        def retFunc(fname, *args, **kwargs):
            x = func(mapping[fname], *args, **kwargs)
            return x

        return retFunc

    return decorator


def scaleQCDBEnchriched(fname):
    SFs = {
        # 2016preVFP
        (0, 100): 19.5 * 1.122e06 * 1e3 / 17657456,
        (0, 200): 19.5 * 8.006e04 * 1e3 / 8886507,
        (0, 300): 19.5 * 1.672e04 * 1e3 / 4978755,
        (0, 500): 19.5 * 1.496e03 * 1e3 / 4433560,
        (0, 700): 19.5 * 3.001e02 * 1e3 / 979344,
        (0, 1000): 19.5 * 4.768e01 * 1e3 / 591966,
        (0, 1500): 19.5 * 4.037e00 * 1e3 / 675657,
        (0, 2000): 19.5 * 6.951e-01 * 1e3 / 668223,
        # 2016postVFP
        (1, 100): 16.8 * 1.124e06 * 1e3 / 19202473,
        (1, 200): 16.8 * 8.040e04 * 1e3 / 9328147,
        (1, 300): 16.8 * 1.668e04 * 1e3 / 5612374,
        (1, 500): 16.8 * 1.502e03 * 1e3 / 4616176,
        (1, 700): 16.8 * 2.995e02 * 1e3 / 903293,
        (1, 1000): 16.8 * 4.756e01 * 1e3 / 663922,
        (1, 1500): 16.8 * 4.024e00 * 1e3 / 698469,
        (1, 2000): 16.8 * 6.963e-01 * 1e3 / 684942,
        # 2017
        (2, 100): 41.5 * 1.125e06 * 1e3 / 37427427,
        (2, 200): 41.5 * 8.013e04 * 1e3 / 19844424,
        (2, 300): 41.5 * 1.669e04 * 1e3 / 11312350,
        (2, 500): 41.5 * 1.506e03 * 1e3 / 10203561,
        (2, 700): 41.5 * 2.998e02 * 1e3 / 1881618,
        (2, 1000): 41.5 * 4.771e01 * 1e3 / 1385631,
        (2, 1500): 41.5 * 4.016e00 * 1e3 / 1458069,
        (2, 2000): 41.5 * 6.979e-01 * 1e3 / 1408971,
        # 2018
        (3, 100): 59.8 * 1.121e06 * 1e3 / 36118282,
        (3, 200): 59.8 * 8.015e04 * 1e3 / 18462183,
        (3, 300): 59.8 * 1.674e04 * 1e3 / 11197722,
        (3, 500): 59.8 * 1.496e03 * 1e3 / 9246898,
        (3, 700): 59.8 * 3.000e02 * 1e3 / 1844165,
        (3, 1000): 59.8 * 4.755e01 * 1e3 / 1330829,
        (3, 1500): 59.8 * 4.030e00 * 1e3 / 1431254,
        (3, 2000): 59.8 * 6.984e-01 * 1e3 / 1357334,
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

    if "100to200" in fname:
        HT = 100
    elif "200to300" in fname:
        HT = 200
    elif "300to500" in fname:
        HT = 300
    elif "500to700" in fname:
        HT = 500
    elif "700to1000" in fname:
        HT = 700
    elif "1000to1500" in fname:
        HT = 1000
    elif "1500to2000" in fname:
        HT = 1500
    elif "2000toInf" in fname:
        HT = 2000
    else:
        print("ERROR: HT range could not be determined for {}".format(fname))

    SF = SFs[(period, HT)]
    return SF

def scaleQCDInclusive(fname):
    SFs = {
        # QCD 2016
        (1,50)   : 59.8 * 186100000.0 * 1E3 / 35474117,
        (1,100)  : 59.8 * 23630000.0  * 1E3 / 73506112,
        (1,200)  : 59.8 * 1554000.0   * 1E3 / 43280518,
        (1,300)  : 59.8 * 323800.0    * 1E3 / 46335846,
        (1,500)  : 59.8 * 30280.0     * 1E3 / 52661606,
        (1,700)  : 59.8 * 6392.0      * 1E3 / 41664730,
        (1,1000) : 59.8 * 1118.0      * 1E3 / 12254238,
        (1,1500) : 59.8 * 108.9       * 1E3 / 9376965,
        (1,2000) : 59.8 * 21.93       * 1E3 / 4867995,
        # QCD 2018
        (3,50)   : 59.8 * 187300000.0 * 1E3 / 38599389,
        (3,100)  : 59.8 * 23590000.0  * 1E3 / 84434559,
        (3,200)  : 59.8 * 1555000.0   * 1E3 / 57336623,
        (3,300)  : 59.8 * 324500.0    * 1E3 / 61609663,
        (3,500)  : 59.8 * 30310.0     * 1E3 / 49184771,
        (3,700)  : 59.8 * 6444.0      * 1E3 / 48506751,
        (3,1000) : 59.8 * 1127.0      * 1E3 / 14394786,
        (3,1500) : 59.8 * 109.8       * 1E3 / 10411831,
        (3,2000) : 59.8 * 21.98       * 1E3 / 5374711,
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
    else: SF = SFs[(period,HT)]
    return SF



def scaleTT(fname):
    lumiTarget = 59.8
    lumiSample = 331506194 / (831.8 * 0.457) * 1e-3
    SF = lumiTarget / lumiSample
    return SF


# -------------------------------------------------
# Signal
# -------------------------------------------------


def scaleSignal(fname):
    points = ['1000_400','1000_600','1000_900',
              '1200_400','1200_600','1200_1100',
              '1300_400','1300_600','1300_1200',
              '1400_400','1400_600','1400_1300',
              '1500_400','1500_600','1500_900','1500_1400',
              '2000_400','2000_600','2000_900','2000_1400','2000_1900'] 

    lumiTarget = 59.8
    NEventsGen = 10000
    lambdapp = 0.1
    
    SFs = {
      1000:lumiTarget * 48000 * lambdapp**2 / NEventsGen,
      1200:lumiTarget * 21000 * lambdapp**2 / NEventsGen,
      1300:lumiTarget * 15000 * lambdapp**2 / NEventsGen,
      1400:lumiTarget * 10000 * lambdapp**2 / NEventsGen,
      1500:lumiTarget * 7300  * lambdapp**2 / NEventsGen,
      2000:lumiTarget * 1600  * lambdapp**2 / NEventsGen,
    }


    part = fname.split("-")[1]
    mStop = int(part.split("_")[0])
    return SFs[mStop]

def scaleSignal313(fname):
    lumiTarget = 59.8
    NEventsGen = 10000
    lambdapp = 0.1
    points = ['1000_400','1000_600','1000_900',
              '1500_400','1500_600','1500_900','1500_1400',
              '2000_400','2000_600','2000_900','2000_1400','2000_1900']

    SFs = {
      1000:lumiTarget * 27000 * lambdapp**2 / NEventsGen,
      1500:lumiTarget * 3800  * lambdapp**2 / NEventsGen,
      2000:lumiTarget * 760   * lambdapp**2 / NEventsGen,
    }

    part = fname.split("-")[1]
    mStop = int(part.split("_")[0])
    return SFs[mStop]


def scaleZQQ(fname):
    SFs = {
        # 2018
        (3, 200): 59.8 * 1012.0 * 1e3 / 15002757,
        (3, 400): 59.8 * 114.2 * 1e3 / 13930474,
        (3, 600): 59.8 * 25.34 * 1e3 / 12029507,
        (3, 800): 59.8 * 12.99 * 1e3 / 9681521,
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


def scaleZNuNu(fname):
    SFs = {
        # 2018
        (3, 100): 59.8 * 267.0 * 1.1347 * 1e3 / 28876062,
        (3, 200): 59.8 * 73.08 * 1.1347 * 1e3 / 22749608,
        (3, 400): 59.8 * 9.921 * 1.1347 * 1e3 / 19676607,
        (3, 600): 59.8 * 2.409 * 1.1347 * 1e3 / 5968910,
        (3, 800): 59.8 * 1.078 * 1.1347 * 1e3 / 2129122,
        (3, 1200): 59.8 * 0.2514 * 1.1347 * 1e3 / 381695,
        (3, 2500): 59.8 * 0.005614 * 1.1347 * 1e3 / 268224,
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


def scaleWQQ(fname):
    SFs = {
        # 2018
        (3, 200): 59.8 * 2549.0 * 1e3 / 14494966,
        (3, 400): 59.8 * 276.5 * 1e3 / 9335298,
        (3, 600): 59.8 * 59.25 * 1e3 / 13633226,
        (3, 800): 59.8 * 28.75 * 1e3 / 13581343,
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


def scaleDiboson(fname):
    SFs = {
        # 2018
        (3, "WW"): 59.8 * 118.7 * 1e3 / 15679000,
        (3, "WZ"): 59.8 * 47.13 * 1e3 / 7940000,
        (3, "ZZ"): 59.8 * 16.523 * 1e3 / 3526000,
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


def scaleST(fname):
    SFs = {
        # 2018
        (3, "s-channel"): 59.8 * 11.03 * 0.457 * 1e3 / 10592646,
        (3, "t-channel_antitop"): 59.8 * 80.95 * 1e3 / 90022642,
        (3, "t-channel_top"): 59.8 * 136.02 * 1e3 / 167111718,
        (3, "tW_antitop"): 59.8 * 35.85 * 1e3 / 7748690,
        (3, "tW_top"): 59.8 * 35.85 * 1e3 / 7955614,
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


sf_funcs = {
    "QCD2018": scaleQCDBEnchriched,
    "QCDBEnriched2018": scaleQCDBEnchriched,
    "QCDInclusive2018": scaleQCDInclusive,
    "TT2018": scaleTT,
    "ZQQ2018": scaleZQQ,
    "ZJetsToQQ2018": scaleZQQ,
    "WQQ2018": scaleWQQ,
    "WJetsToQQ": scaleWQQ,
    "WJetsToQQ2018": scaleWQQ,
    "Diboson2018": scaleDiboson,
    "ST2018": scaleST,
    "ZNuNu2018": scaleZNuNu,
    "signal": scaleSignal,
}

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

    sf_func = sf_funcs[sample]
    sf_func = scaled_lumi(137.62, 59.8)(sf_func)
    if is_signal:
        raw_files = glob.glob("{}/*.root".format(input_dir))
        for f in raw_files:
            p = os.path.split(f)[-1].split('-')[1]
            haddNano(os.path.join(output, "signal_{}".format(p)), [f], sf_func)
    else:
        afiles = associateFiles(sample_file, input_dir)
        sf_func = translate_files(afiles)(sf_func)
        raw_files = [f for f in afiles]
        haddNano(output, raw_files, sf_func)
