import collections
targetFile = 'samples/QCDInclusive2018.txt'
counts = collections.OrderedDict()
with open(targetFile, 'r') as f:
	lines = f.read().split('\n')
	for line in lines:
		if line == '': continue
		splitLine = list(line.split('_'))[1]
		HTBinStr = list(splitLine[2:].split('to'))
		bins = (HTBinStr[0], HTBinStr[1])
		counts[bins] = counts[bins] + 1 if bins in counts else 1
for count in counts:
	print(count)
