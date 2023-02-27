import itertools


def orderedMatcher(jets, particles, dr_match=0.3):
    def get4(obj):
        return obj.p4()
    def DR(self,other):
        return self.DeltaR(other)
    matched = []
    remaining_idxs = set(range(len(jets)))
    used_idxs = set()
    for i,p in enumerate(particles):
        diff=lambda x : DR(get4(p),get4(jets[x]))
        if remaining_idxs: smallest = min(remaining_idxs, key=diff)
        if remaining_idxs and diff(smallest) < dr_match:
            remaining_idxs.remove(smallest)
            used_idxs.add(smallest)
            matched.append((smallest,None))
        elif used_idxs:
            smallest = min(used_idxs, key=diff)
            if diff(smallest) < dr_match: matched.append((None,smallest))
            else: matched.append((None,None))
        else: matched.append((None,None))
    return matched

