import itertools

def matcher(jets, particles, dr_match=0.2, combo_count=6):
    def get4(obj):
        return obj.p4()

    def DR(self,other):
        return self.DeltaR(other)
        
    jets_to_use = jets[:combo_count]
    loss = None 
    best_match = None
    for perm in itertools.permutations(range(len(jets_to_use))):
        this_perm = [jets[x] for x in perm]
        closeness = [DR(get4(particles[i]),get4(this_perm[i])) for i in range(len(jets_to_use))]
        if any(x > dr_match for x in closeness):
            continue
        new_loss = sum(closeness)
#        print("Permutation {} has loss {}".format(perm, new_loss))
        if loss is None  or new_loss < loss:
            loss = new_loss
            best_match = perm
        
    return best_match

