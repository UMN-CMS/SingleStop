from itertools import product
def dR(obj1, obj2):
	return obj1.p4().DeltaR(obj2.p4())

def pT(obj):
	return obj.p4().Pt()

def jetMatcher(jetsA, jetsB, dRMatch = 0.4, pTMatch = 3):
	matched = {}
	combos = list(product(range(len(jetsA)), range(len(jetsB))))
	for (a, b) in combos:
		jetA = jetsA[a]
		jetB = jetsB[b]
		if dR(jetA, jetB) <= dRMatch and ( abs(pT(jetA) - pT(jetB)) / pT(jetB) ) <= pTMatch:
			matched[a] = matched[a] + [b] if a in matched else [b]
	multipleMatches = 0
	for a in matched:
		possibleMatches = matched[a]
		sorted(possibleMatches, key = lambda b: dR(jetsA[a], jetsB[b]))
		matched[a] = possibleMatches
		if len(matched[a]) > 1: multipleMatches += 1
	return checker(matched, jetsA, jetsB), multipleMatches
	
def checker(matched, jetsA, jetsB):
	primaryMatches = {a: 0 for a in matched}
	noMatches = []
	for a in matched:
		if matched[a] == []: 
			noMatches.append(a)
			continue
		primaryMatches[a] = matched[a][0]
	for m in noMatches:
		matched.pop(m)
	if len(primaryMatches.values()) == len(set(primaryMatches.values())): 
		for a in matched:
			matched[a] = matched[a][0]
		return matched	

	sameMatches = {}
	for a, b in primaryMatches.items():
		sameMatches[b] = sameMatches[b] + [a] if b in sameMatches else [a]
	for b in sameMatches:
		if len(sameMatches[b]) <= 1: continue
		closestMatch = min(sameMatches[b], key = lambda a: dR(jetsA[a], jetsB[b]))
		for a in sameMatches[b]:
			if a in matched and a != closestMatch: matched[a].remove(b)
	return checker(matched, jetsA, jetsB)	
