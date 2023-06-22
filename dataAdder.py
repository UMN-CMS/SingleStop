import ROOT, argparse, os, subprocess
parser = argparse.ArgumentParser(description = 'Adding data files')
parser.add_argument('--inputDir', type = str, required = True, help = 'Input Path')
parser.add_argument('--outputDir', type = str, required = True, help = 'Output Path')
args = parser.parse_args()
inputDir = args.inputDir
outputDir = args.outputDir

numFiles = subprocess.check_output('ls {} | wc -l'.format(args.inputDir), shell = True)
fileList = []
for i in range(int(numFiles)):
  fileList.append('{}/Data2018-{}-temp.root'.format(outputDir, i + 1))
  fTemp = ROOT.TFile.Open('{}/Data2018-{}-temp.root'.format(outputDir, i + 1), 'RECREATE')
  f = ROOT.TFile.Open('{}/Data2018-{}.root'.format(inputDir, i + 1), 'READ')
  plots = f.GetDirectory('plots')
  keys = [key.GetName() for key in plots.GetListOfKeys()]
  for key in keys:
    h = plots.Get(key)
    fTemp.WriteObject(h, h.GetName())
  fTemp.Write()
  fTemp.Close()

os.system('hadd -f {}/Data2018.root {}'.format(outputDir, ' '.join(fileList)))
os.system('rm {}/*temp.root'.format(outputDir))
