import itertools
import ROOT

def massRecoTrivial(jets, bwp):
    return sum(x.p4() for x in jets[:4]).M()

def massRecoPtMin(jets, bwp, to_consider=6):
    vec = None
    for perm in itertools.combinations(jets, min(to_consider, len(jets))):
        new_jet = sum((x.p4() for x in perm[:4]),ROOT.TLorentzVector())
        if vec is None or new_jet.Pt() < vec.Pt():
            vec = new_jet
    return vec.M()

def massRecoBSeed(jets, bwp, nbs = 1):
    return sum(x.p4() for x in jets).M()

