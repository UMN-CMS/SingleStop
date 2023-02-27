#!/usr/bin/env python

import os, sys, ROOT, argparse

parser = argparse.ArgumentParser(description='MC Scaler')
parser.add_argument('--sample',type=str,required=True,choices=['signal','TT2018','QCD','QCD2018','ZQQ2018','ST2018','WQQ2018','ZNuNu2018','Diboson2018'],help='Sample to scale')
parser.add_argument('--input',type=str,required=True,help='Path to input files')
parser.add_argument('--output',type=str,required=True,help='Path to output scaled file')
args = parser.parse_args()

scale2018toRun2 = True if '2018' in args.sample else False

sample = args.sample
inputDir = args.input
outputDir = args.output

if not os.path.exists(outputDir):
  os.makedirs(outputDir)

try: os.remove('{}/{}.root'.format(outputDir,sample))
except OSError: pass

if   args.sample == 'TT':          sampleFile = 'TTToHadronic.txt'
elif args.sample == 'TT2018':      sampleFile = 'TTToHadronic2018.txt'
elif args.sample == 'QCD':         sampleFile = 'QCDBEnriched.txt'
elif args.sample == 'QCD2018':     sampleFile = 'QCDBEnriched2018.txt'
elif args.sample == 'ZQQ2018':     sampleFile = 'ZJetsToQQ2018.txt'
elif args.sample == 'ST2018':      sampleFile = 'STHadronic2018.txt'
elif args.sample == 'WQQ2018':     sampleFile = 'WJetsToQQ2018.txt'
elif args.sample == 'ZNuNu2018':   sampleFile = 'ZJetsToNuNu2018.txt'
elif args.sample == 'Diboson2018': sampleFile = 'Diboson2018.txt'
elif args.sample != 'signal': print('ERROR: Unexpected sample argument')

#-------------------------------------------------
# QCD
#-------------------------------------------------

if 'QCD' in sample:

  SFs = {
  	# 2016preVFP
  	(0,100) : 19.5 * 1.122e+06 * 1E3 / 17657456,
  	(0,200) : 19.5 * 8.006e+04 * 1E3 / 8886507,
  	(0,300) : 19.5 * 1.672e+04 * 1E3 / 4978755,
  	(0,500) : 19.5 * 1.496e+03 * 1E3 / 4433560,
  	(0,700) : 19.5 * 3.001e+02 * 1E3 / 979344,
  	(0,1000): 19.5 * 4.768e+01 * 1E3 / 591966,
  	(0,1500): 19.5 * 4.037e+00 * 1E3 / 675657,
  	(0,2000): 19.5 * 6.951e-01 * 1E3 / 668223,
  	# 2016postVFP
        (1,100) : 19.5 * 1.124e+06 * 1E3 / 19202473, 
        (1,200) : 19.5 * 8.040e+04 * 1E3 / 9328147,
        (1,300) : 19.5 * 1.668e+04 * 1E3 / 5612374,
        (1,500) : 19.5 * 1.502e+03 * 1E3 / 4616176,
        (1,700) : 19.5 * 2.995e+02 * 1E3 / 903293,
        (1,1000): 19.5 * 4.756e+01 * 1E3 / 663922,
        (1,1500): 19.5 * 4.024e+00 * 1E3 / 698469,
        (1,2000): 19.5 * 6.963e-01 * 1E3 / 684942,
  	# 2017
        (2,100) : 19.5 * 1.125e+06 * 1E3 / 37427427,
        (2,200) : 19.5 * 8.013e+04 * 1E3 / 19844424,
        (2,300) : 19.5 * 1.669e+04 * 1E3 / 11312350,
        (2,500) : 19.5 * 1.506e+03 * 1E3 / 10203561,
        (2,700) : 19.5 * 2.998e+02 * 1E3 / 1881618,
        (2,1000): 19.5 * 4.771e+01 * 1E3 / 1385631,
        (2,1500): 19.5 * 4.016e+00 * 1E3 / 1458069,
        (2,2000): 19.5 * 6.979e-01 * 1E3 / 1408971,
  	# 2018
        (3,100) : 19.5 * 1.121e+06 * 1E3 / 36118282,
        (3,200) : 19.5 * 8.015e+04 * 1E3 / 18462183,
        (3,300) : 19.5 * 1.674e+04 * 1E3 / 11197722,
        (3,500) : 19.5 * 1.496e+03 * 1E3 / 9246898,
        (3,700) : 19.5 * 3.000e+02 * 1E3 / 1844165,
        (3,1000): 19.5 * 4.755e+01 * 1E3 / 1330829,
        (3,1500): 19.5 * 4.030e+00 * 1E3 / 1431254,
        (3,2000): 19.5 * 6.984e-01 * 1E3 / 1470453,
  }
  
  if scale2018toRun2:
  
    lumiTarget = 137.62
    lumiSample = 59.83
  
    for key,value in SFs.items():
      if key[0] == 3: SFs[key] *= lumiTarget / lumiSample
  
  with open('../samples/{}'.format(sampleFile)) as f:
    for i,line in enumerate(f):
     
      period = -1
      HT = -1
  
      if   'preVFP' in line: period = 0
      elif 'UL16' in line: period = 1
      elif 'UL17' in line: period = 2
      elif 'UL18' in line: period = 3
      else: print('ERROR: Data period could not be determined for {}'.format(line))
  
      if   '100to200' in line:   HT = 100	
      elif '200to300' in line:   HT = 200    
      elif '300to500' in line:   HT = 300    
      elif '500to700' in line:   HT = 500    
      elif '700to1000' in line:  HT = 700   
      elif '1000to1500' in line: HT = 1000  
      elif '1500to2000' in line: HT = 1500  
      elif '2000toInf' in line:  HT = 2000
      else: print('ERROR: HT range could not be determined for {}'.format(line))
  
      if period == -1 or HT == -1: sys.exit()
      else: SF = SFs[(period,HT)]
  
      fTemp = ROOT.TFile.Open('{}/{}-{}-temp.root'.format(outputDir,sample,i + 1),'RECREATE')
      f = ROOT.TFile.Open('{}/{}-{}.root'.format(inputDir,sample,i + 1))
      plots = f.GetDirectory('plots')
      keys = [key.GetName() for key in plots.GetListOfKeys()]
      for key in keys:
        h = plots.Get(key)
        h.Scale(SF)
        fTemp.WriteObject(h,h.GetName())
  
      fTemp.Write()
      fTemp.Close()
  
  fileList = []
  for j in range(i + 1):
    fileList.append('{}/{}-{}-temp.root'.format(outputDir,sample,j + 1))
  os.system('hadd -f {}/{}.root {}'.format(outputDir,sample,' '.join(fileList)))
  os.system('rm {}/*temp.root'.format(outputDir))

#-------------------------------------------------
# TT2018
#-------------------------------------------------

if sample == 'TT2018':

  lumiTarget = 137.62
  lumiSample = 334206000 / (831.8 * 0.457) * 1E-3
  
  SF = lumiTarget / lumiSample
  
  os.system('hadd -f {}/{}-temp.root {}/{}*.root'.format(outputDir,sample,inputDir,sample))

  fScaled = ROOT.TFile.Open('{}/{}.root'.format(outputDir,sample),'RECREATE')
  fTemp = ROOT.TFile.Open('{}/{}-temp.root'.format(outputDir,sample),'READ')
  plots = fTemp.GetDirectory('plots')
  keys = [key.GetName() for key in plots.GetListOfKeys()]
  for key in keys:
    h = plots.Get(key)
    h.Scale(SF)
    fScaled.WriteObject(h,h.GetName())
  
  fScaled.Write()
  fScaled.Close()

  os.system('rm {}/*temp.root'.format(outputDir))

#-------------------------------------------------
# Signal
#-------------------------------------------------

if sample == 'signal':

  points = ['1000_400','1000_900','1500_600','1500_1400','2000_900','2000_1900','1000_600','1500_400','2000_400','2000_1400','1500_900','2000_600']
  
  lumiTarget = 137.62
  NEventsGen = 10000
  lambdapp312 = 0.1
  
  SFs = {
    1000:lumiTarget * 47000 * lambdapp312**2 / NEventsGen,
    1500:lumiTarget * 7200  * lambdapp312**2 / NEventsGen,
    2000:lumiTarget * 1600  * lambdapp312**2 / NEventsGen,
  }
  
  for point in points:
  
    mStop = int(point.split('_')[0])
  
    fScaled = ROOT.TFile.Open('{}/{}_{}.root'.format(outputDir,sample,point),'RECREATE')
    f = ROOT.TFile.Open('{}/{}_{}.root'.format(inputDir,sample,point),'READ')
    plots = f.GetDirectory('plots')
    keys = [key.GetName() for key in plots.GetListOfKeys()]
    for key in keys:
      h = plots.Get(key)
      h.Scale(SFs[mStop])
      fScaled.WriteObject(h,h.GetName())
  
    fScaled.Write()
    fScaled.Close()

#-------------------------------------------------
# ZQQ2018
#-------------------------------------------------

if sample == 'ZQQ2018':

  SFs = {
  	# 2018
        (3,200): 59.83 * 1012.0 * 1E3 / 15002757,
        (3,400): 59.83 * 114.2  * 1E3 / 13930474,
        (3,600): 59.83 * 25.34  * 1E3 / 12029507,
        (3,800): 59.83 * 12.99  * 1E3 / 9681521,
  }
  
  if scale2018toRun2:
  
    lumiTarget = 137.62
    lumiSample = 59.83
  
    for key,value in SFs.items():
      if key[0] == 3: SFs[key] *= lumiTarget / lumiSample
  
  with open('../samples/{}'.format(sampleFile)) as f:
    for i,line in enumerate(f):
     
      period = -1
      HT = -1
  
      if   'preVFP' in line: period = 0
      elif 'UL16' in line: period = 1
      elif 'UL17' in line: period = 2
      elif 'UL18' in line: period = 3
      else: print('ERROR: Data period could not be determined for {}'.format(line))
  
      if   '200to400' in line:   HT = 200	
      elif '400to600' in line:   HT = 400    
      elif '600to800' in line:   HT = 600    
      elif '800toInf' in line:   HT = 800    
      else: print('ERROR: HT range could not be determined for {}'.format(line))
  
      if period == -1 or HT == -1: sys.exit()
      else: SF = SFs[(period,HT)]
  
      fTemp = ROOT.TFile.Open('{}/{}-{}-temp.root'.format(outputDir,sample,i + 1),'RECREATE')
      f = ROOT.TFile.Open('{}/{}-{}.root'.format(inputDir,sample,i + 1))
      plots = f.GetDirectory('plots')
      keys = [key.GetName() for key in plots.GetListOfKeys()]
      for key in keys:
        h = plots.Get(key)
        h.Scale(SF)
        fTemp.WriteObject(h,h.GetName())
  
      fTemp.Write()
      fTemp.Close()
  
  fileList = []
  for j in range(i + 1):
    fileList.append('{}/{}-{}-temp.root'.format(outputDir,sample,j + 1))
  os.system('hadd -f {}/{}.root {}'.format(outputDir,sample,' '.join(fileList)))
  os.system('rm {}/*temp.root'.format(outputDir))

#-------------------------------------------------
# ZNuNu2018
#-------------------------------------------------

if sample == 'ZNuNu2018':

  SFs = {
  	# 2018
        (3,100):  59.83 * 267.0    * 1.1347 * 1E3 / 28876062,
        (3,200):  59.83 * 73.08    * 1.1347 * 1E3 / 22749608,
        (3,400):  59.83 * 9.921    * 1.1347 * 1E3 / 19676607,
        (3,600):  59.83 * 2.409    * 1.1347 * 1E3 / 5968910,
        (3,800):  59.83 * 1.078    * 1.1347 * 1E3 / 2129122,
        (3,1200): 59.83 * 0.2514   * 1.1347 * 1E3 / 381695,
        (3,2500): 59.83 * 0.005614 * 1.1347 * 1E3 / 268224,
  }
  
  if scale2018toRun2:
  
    lumiTarget = 137.62
    lumiSample = 59.83
  
    for key,value in SFs.items():
      if key[0] == 3: SFs[key] *= lumiTarget / lumiSample
  
  with open('../samples/{}'.format(sampleFile)) as f:
    for i,line in enumerate(f):
     
      period = -1
      HT = -1
  
      if   'preVFP' in line: period = 0
      elif 'UL16' in line: period = 1
      elif 'UL17' in line: period = 2
      elif 'UL18' in line: period = 3
      else: print('ERROR: Data period could not be determined for {}'.format(line))
  
      if   '100To200' in line:   HT = 100	
      elif '200To400' in line:   HT = 200    
      elif '400To600' in line:   HT = 400
      elif '600To800' in line:   HT = 400 
      elif '800To1200' in line:  HT = 800    
      elif '1200To2500' in line: HT = 1200
      elif '2500ToInf' in line:  HT = 2500
      else: print('ERROR: HT range could not be determined for {}'.format(line))
  
      if period == -1 or HT == -1: sys.exit()
      else: SF = SFs[(period,HT)]
  
      fTemp = ROOT.TFile.Open('{}/{}-{}-temp.root'.format(outputDir,sample,i + 1),'RECREATE')
      f = ROOT.TFile.Open('{}/{}-{}.root'.format(inputDir,sample,i + 1))
      plots = f.GetDirectory('plots')
      keys = [key.GetName() for key in plots.GetListOfKeys()]
      for key in keys:
        h = plots.Get(key)
        h.Scale(SF)
        fTemp.WriteObject(h,h.GetName())
  
      fTemp.Write()
      fTemp.Close()
  
  fileList = []
  for j in range(i + 1):
    fileList.append('{}/{}-{}-temp.root'.format(outputDir,sample,j + 1))
  os.system('hadd -f {}/{}.root {}'.format(outputDir,sample,' '.join(fileList)))
  os.system('rm {}/*temp.root'.format(outputDir))

#-------------------------------------------------
# WQQ2018
#-------------------------------------------------

if sample == 'WQQ2018':

  SFs = {
  	# 2018
        (3,200): 59.83 * 2549.0 * 1E3 / 14494966,
        (3,400): 59.83 * 276.5  * 1E3 / 9335298,
        (3,600): 59.83 * 59.25  * 1E3 / 13633226,
        (3,800): 59.83 * 28.75  * 1E3 / 13581343,
  }
  
  if scale2018toRun2:
  
    lumiTarget = 137.62
    lumiSample = 59.83
  
    for key,value in SFs.items():
      if key[0] == 3: SFs[key] *= lumiTarget / lumiSample
  
  with open('../samples/{}'.format(sampleFile)) as f:
    for i,line in enumerate(f):
     
      period = -1
      HT = -1
  
      if   'preVFP' in line: period = 0
      elif 'UL16' in line: period = 1
      elif 'UL17' in line: period = 2
      elif 'UL18' in line: period = 3
      else: print('ERROR: Data period could not be determined for {}'.format(line))
  
      if   '200to400' in line:   HT = 200	
      elif '400to600' in line:   HT = 400    
      elif '600to800' in line:   HT = 600    
      elif '800toInf' in line:   HT = 800    
      else: print('ERROR: HT range could not be determined for {}'.format(line))
  
      if period == -1 or HT == -1: sys.exit()
      else: SF = SFs[(period,HT)]
  
      fTemp = ROOT.TFile.Open('{}/{}-{}-temp.root'.format(outputDir,sample,i + 1),'RECREATE')
      f = ROOT.TFile.Open('{}/{}-{}.root'.format(inputDir,sample,i + 1))
      plots = f.GetDirectory('plots')
      keys = [key.GetName() for key in plots.GetListOfKeys()]
      for key in keys:
        h = plots.Get(key)
        h.Scale(SF)
        fTemp.WriteObject(h,h.GetName())
  
      fTemp.Write()
      fTemp.Close()
  
  fileList = []
  for j in range(i + 1):
    fileList.append('{}/{}-{}-temp.root'.format(outputDir,sample,j + 1))
  os.system('hadd -f {}/{}.root {}'.format(outputDir,sample,' '.join(fileList)))
  os.system('rm {}/*temp.root'.format(outputDir))

#-------------------------------------------------
# Diboson2018
#-------------------------------------------------

if sample == 'Diboson2018':

  SFs = {
  	# 2018
        (3,'WW'): 59.83 * 118.7  * 1E3 / 15679000,
        (3,'WZ'): 59.83 * 47.13  * 1E3 / 7940000,
        (3,'ZZ'): 59.83 * 16.523 * 1E3 / 3526000,
  }
  
  if scale2018toRun2:
  
    lumiTarget = 137.62
    lumiSample = 59.83
  
    for key,value in SFs.items():
      if key[0] == 3: SFs[key] *= lumiTarget / lumiSample
  
  with open('../samples/{}'.format(sampleFile)) as f:
    for i,line in enumerate(f):
     
      period = -1
      mode = -1
  
      if   'preVFP' in line: period = 0
      elif 'UL16' in line: period = 1
      elif 'UL17' in line: period = 2
      elif 'UL18' in line: period = 3
      else: print('ERROR: Data period could not be determined for {}'.format(line))
  
      if   'WW' in line:   mode = 'WW'	
      elif 'WZ' in line:   mode = 'WZ'    
      elif 'ZZ' in line:   mode = 'ZZ' 
      else: print('ERROR: HT range could not be determined for {}'.format(line))
  
      if period == -1 or mode == -1: sys.exit()
      else: SF = SFs[(period,mode)]
  
      fTemp = ROOT.TFile.Open('{}/{}-{}-temp.root'.format(outputDir,sample,i + 1),'RECREATE')
      f = ROOT.TFile.Open('{}/{}-{}.root'.format(inputDir,sample,i + 1))
      plots = f.GetDirectory('plots')
      keys = [key.GetName() for key in plots.GetListOfKeys()]
      for key in keys:
        h = plots.Get(key)
        h.Scale(SF)
        fTemp.WriteObject(h,h.GetName())
  
      fTemp.Write()
      fTemp.Close()
  
  fileList = []
  for j in range(i + 1):
    fileList.append('{}/{}-{}-temp.root'.format(outputDir,sample,j + 1))
  os.system('hadd -f {}/{}.root {}'.format(outputDir,sample,' '.join(fileList)))
  os.system('rm {}/*temp.root'.format(outputDir))

#-------------------------------------------------
# ST2018
#-------------------------------------------------

if sample == 'ST2018':

  SFs = {
  	# 2018
        (3,'s-channel'):         59.83 * 11.03 * 0.457 * 1E3 / 16259000,
        (3,'t-channel_antitop'): 59.83 * 80.95         * 1E3 / 95627000,
        (3,'t-channel_top'):     59.83 * 136.02        * 1E3 / 178336000,
        (3,'tW_antitop'):        59.83 * 35.85         * 1E3 / 7749000,
        (3,'tW_top'):            59.83 * 35.85         * 1E3 / 7956000,
  }
  
  if scale2018toRun2:
  
    lumiTarget = 137.62
    lumiSample = 59.83
  
    for key,value in SFs.items():
      if key[0] == 3: SFs[key] *= lumiTarget / lumiSample
  
  with open('../samples/{}'.format(sampleFile)) as f:
    for i,line in enumerate(f):
     
      period = -1
      mode = -1
  
      if   'preVFP' in line: period = 0
      elif 'UL16' in line: period = 1
      elif 'UL17' in line: period = 2
      elif 'UL18' in line: period = 3
      else: print('ERROR: Data period could not be determined for {}'.format(line))
  
      if   's-channel' in line:   mode = 's-channel'
      elif 't-channel_antitop' in line:   mode = 't-channel_antitop'    
      elif 't-channel_top' in line:   mode = 't-channel_top'
      elif 'tW_antitop' in line:   mode = 'tW_antitop'
      elif 'tW_top' in line:   mode = 'tW_top'
      else: print('ERROR: HT range could not be determined for {}'.format(line))
  
      if period == -1 or mode == -1: sys.exit()
      else: SF = SFs[(period,mode)]
  
      fTemp = ROOT.TFile.Open('{}/{}-{}-temp.root'.format(outputDir,sample,i + 1),'RECREATE')
      f = ROOT.TFile.Open('{}/{}-{}.root'.format(inputDir,sample,i + 1))
      plots = f.GetDirectory('plots')
      keys = [key.GetName() for key in plots.GetListOfKeys()]
      for key in keys:
        h = plots.Get(key)
        h.Scale(SF)
        fTemp.WriteObject(h,h.GetName())
  
      fTemp.Write()
      fTemp.Close()
  
  fileList = []
  for j in range(i + 1):
    fileList.append('{}/{}-{}-temp.root'.format(outputDir,sample,j + 1))
  os.system('hadd -f {}/{}.root {}'.format(outputDir,sample,' '.join(fileList)))
  os.system('rm {}/*temp.root'.format(outputDir))

