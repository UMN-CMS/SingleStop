def bJetMatcher(genJets, recoJets, dR_match = 0.4, pT_match = 3):
	def getp4(obj):
		return obj.p4()
	def dR(obj1, obj2):
		return obj1.p4().DeltaR(obj2.p4())
	def pT(obj):
		return obj.p4().Pt()
	matched = []
	remaining_idxs = set(range(len(genJets)))
	used_idxs = set()
	reverse_matched = {}
	for i, j in enumerate(recoJets):
		diff = lambda x: dR(j, genJets[x])
		if remaining_idxs: 
			smallest = min(remaining_idxs, key = diff)
			if diff(smallest) < dR_match and ( (1 / pT_match) * pT(genJets[smallest]) < pT(j) < (pT_match) * pT(genJets[smallest]) ):
				if smallest in used_idxs:
					if dR(genJets[smallest], j) < dR(genJets[smallest], recoJets[reverse_matched[smallest]]):
						matched.remove((reverse_matched[smallest], smallest))
						matched.append((i, smallest))
						reverse_matched[smallest] = i
				used_idxs.add(smallest)
				matched.append((i, smallest))
				reverse_matched[smallest] = i
		else: matched.append((None, None))
	return matched
