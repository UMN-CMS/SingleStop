## Setup nanoAOD-tools
```
cmsrel CMSSW_10_6_19_patch2
cd CMSSW_10_6_19_patch2/src
cmsenv
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
scram b
```
## Setup grid
```
voms-proxy-init -voms cms
```
## Setup single stop repository
```
cd PhysicsTools/NanoAODTools/python/postprocessing/
git clone git@github.com:UMN-CMS/SingleStop.git
```
## Set architecture
This should be done every time you log in.
```
source setup.sh
``` 
## Running the analyzer
```
python singleStopAnalyzer.py --sample [signal/TT/QCD] --tag [OUTPUT TAG] -n [BKG FILE INDEX]
```
### Input
For running over signal, the desired signal points can be specified in the `points` list. 
For QCD and TT backgrounds, the corresponding files are specified in the respective files in the `samples` directory.
The `-n` argument must be used for background in order to specify the index of the files to run over. 
### Output
Output histograms are added in labelled directories within `output`. 
### Condor
#### QCD
First, copy over the `singleStopAnalyzer.py` script and in the import section at the top, remove all `PhysicsTools.NanoAODTools.postprocessing.` prefixes.
Then run `doSubmitQCD.sh` to subnmit jobs for all files in `samples/QCDBEnriched.txt`.
