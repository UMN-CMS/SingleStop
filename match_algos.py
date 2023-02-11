import itertools

def globalMatcher(jets, particles, dr_match=0.3, combo_count=8, require_all_matches=False):
    def get4(obj):
        return obj.p4()

    def DR(self,other):
        return self.DeltaR(other)
        
    jets_to_use = jets[:combo_count]
    num_to_test = min(len(particles), len(jets_to_use))
    njets = len(jets_to_use)
    loss = None 
    best_match = None
    for perm in itertools.permutations(range(njets), mum_to_test):
        #this_perm = [jets[x] for x in perm]
        closeness = [DR(get4(particles[i]),get4(jets[perm[i]]))**2 for i in range(num_to_test)]
        if require_all_matches and any(x > dr_match for x in closeness):
            continue
        new_loss = sum(closeness)
#       print("Permutation {} has loss {}".format(perm, new_loss))
        if loss is None  or new_loss < loss:
            loss = new_loss
            best_match = perm
    return best_match

def orderedMatcher(jets, particles, dr_match=0.3):
    def get4(obj):
        return obj.p4()
    def DR(self,other):
        return self.DeltaR(other)
    matched = []
    remaining_idxs = set(range(len(jets)))

    for i,p in enumerate(particles):
        if not remaining_idxs:
            matched.append(None)
            continue
        smallest = min(remaining_idxs, key=lambda x : DR(get4(p),get4(jets.__getitem__(x))))
        if DR(get4(p),get4(jets.__getitem__(smallest))) > dr_match:
            matched.append(None)
            continue
        remaining_idxs.remove(smallest)
        matched.append(smallest)

    return matched

def priorityOrderedMatcher(jets, particles, dr_match=0.3):
    def get4(obj):
        return obj.p4()
    def DR(self,other):
        return self.DeltaR(other)
    matched = []
    remaining_idxs = list(range(len(jets)))

    for i,p in enumerate(particles):
        smallest = next((x for x in remaining_idxs if DR(get4(p),get4(jets.__getitem__(x))) < dr_match), None)
        if not smallest:
            smallest = None
        else:
            remaining_idxs.remove(smallest)
        matched.append(smallest)
        if not remaining_idxs:
            break

    return matched
        
