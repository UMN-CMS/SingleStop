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


def haddNano(ofname, files, max_events=None):
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

    n_events = 0
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
            branchNames = set([x.GetName() for x in obj.GetListOfBranches()])

        for fht in fileHandles[1:]:
            fn = fht[0]
            fh = fht[1]
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


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def chunkMax(lst, f, target):
    target *= 1024**3
    ret = []
    curr = []
    accum = 0
    for x in lst:
        curr.append(x)
        accum += f(x)
        print(f(x))
        if accum >= target:
            ret.append(curr)
            curr = []
            accum = 0
    ret.append(curr)
    return ret


def qcdMatcher(x):
    m = re.search(r"(RunIISummer.+NanoAODv9).+(QCD_HT\d+to(?:\d+|Inf))", x)
    return "{}_{}".format(m.group(1), m.group(2))


def ttMatcher(x):
    return "TT"


def zqqMatcher(x):
    m = re.search(r"(RunIISummer.+NanoAODv9)/(.+\d+to(?:\d+|Inf))", x)
    return "{}_{}".format(m.group(1), m.group(2))


def zNuNuMatcher(x):
    m = re.search(r"(RunIISummer.+NanoAODv9)/(.+HT-\d+To(?:\d+|Inf))", x)
    return "{}_{}".format(m.group(1), m.group(2))


def wqqMatcher(x):
    m = re.search(r"(RunIISummer.+NanoAODv9)/(.+\d+to(?:\d+|Inf))", x)
    return "{}_{}".format(m.group(1), m.group(2))


def dibosonMatcher(x):
    m = re.search(r"(RunIISummer.+NanoAODv9)/(.+13TeV)", x)
    return "{}_{}".format(m.group(1), m.group(2))


def stMatcher(x):
    m = re.search(r"(RunIISummer.+NanoAODv9)/(.+13TeV)", x)
    return "{}_{}".format(m.group(1), m.group(2))


matchers = {
    "QCDInclusive2018": qcdMatcher,
    "TT2018": ttMatcher,
    "ZQQ2018": zqqMatcher,
    "WQQ2018": wqqMatcher,
    "Diboson2018": dibosonMatcher,
    "ST2018": stMatcher,
    "ZNuNu2018": zNuNuMatcher,
}


def associateFiles(sample_file, input_dir, split_num=None):
    fnames = []
    with open(sample_file, "r") as f:
        fnames = [line.strip() for line in f]
    files = glob.glob("{}/*.root".format(input_dir))

    def g(rootfile):
        return next(f for f in fnames if os.path.split(rootfile)[1].split("_")[0] in f)

    fmap = {rootfile: (g(rootfile)) for rootfile in files}
    return fmap


def makeGroups(fmap, matcher):
    pairs = list(fmap.items())
    pairs = sorted(pairs, key=lambda x: x[1])
    groups = it.groupby(pairs, lambda x: matcher(x[1]))
    return {x: list(y[0] for y in z) for x, z in groups}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="hadd")
    parser.add_argument(
        "--input-dir", type=str, required=True, help="Path to input files"
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Path to output scaled file"
    )
    parser.add_argument("--max-size-gb", type=float, required=True, help="Max size gb")
    args = parser.parse_args()

    if "TT2018" in args.input_dir:
        sample_file = "TTToHadronic2018.txt"
    elif "QCD2018" in args.input_dir:
        sample_file = "QCDBEnriched2018.txt"
    elif "QCDInclusive2018" in args.input_dir:
        sample_file = "QCDInclusive2018.txt"
    elif "ZQQ2018" in args.input_dir:
        sample_file = "ZJetsToQQ2018.txt"
    elif "ST2018" in args.input_dir:
        sample_file = "STHadronic2018.txt"
    elif "WQQ2018" in args.input_dir:
        sample_file = "WJetsToQQ2018.txt"
    elif "ZNuNu2018" in args.input_dir:
        sample_file = "ZJetsToNuNu2018.txt"
    elif "Diboson2018" in args.input_dir:
        sample_file = "Diboson2018.txt"

    if "Data" not in args.input_dir:
        sample_dir = "samples"
        sample_file = "{}/{}".format(sample_dir, sample_file)
        afiles = associateFiles(sample_file, args.input_dir, 20)
        matcher = next(x for y, x in matchers.items() if y in args.input_dir)
        sample_groups =  makeGroups(afiles, matcher)
        sample_groups = {x : chunkMax(g, os.path.getsize, args.max_size_gb) for x,g in sample_groups.items()}
    else:
        sample_groups = {"Data2018" : chunkMax(glob.glob("{}/*.root".format(args.input_dir)), os.path.getsize, args.max_size_gb)}

        
    # files = glob.glob("{}/*.root".format(args.input_dir))
    #print(sample_groups)
    for sample_name, groups in sample_groups.items():
        for i, group in enumerate(groups):
            #print("{}/{}_{}.root".format(args.output, sample_name, i))
            haddNano("{}/{}_{}.root".format(args.output, sample_name, i), group)

    #    #haddNano(output, raw_files, sf_func)
