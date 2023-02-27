## Setup nanoAOD-tools
Make sure you've set up the default CMS software environment (add this to your `~/.bash_profile` if it's not already). 
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
```
Then, setup CMSSW, clone the `NanoAODTools` repository, and compile.
```
cmsrel CMSSW_10_6_19_patch2
cd CMSSW_10_6_19_patch2/src
cmsenv
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
scram b
```
## Setup single stop repository
Checkout the `SingleStop` repository.
```
cd PhysicsTools/NanoAODTools/python/postprocessing/
git clone git@github.com:UMN-CMS/SingleStop.git
```
## Setup grid
You will need a valid grid certificate in order to run over any signal or bkg files.
The command below open a grid proxy for 8 days (the maximum length).
```
voms-proxy-init --rfc --voms cms -valid 192:00
```
## Set architecture
CMSSW will complain if you don't specify the correct architecture.
The `setup.sh` script simply sets the correct one.
This should be done every time you log in.
```
source setup.sh
``` 
## Running the analyzer
```
python singleStopAnalyzer.py --sample [SAMPLE] --tag [TAG] -n [BKG FILE INDEX] --points [MSTOP_MCHI]
```
### Input
For running over signal, the desired signal points can be specified with the `--points` argument (`MSTOP_MCHI`), comma-separated. 
For signal, all points will be used unless specified.
For all backgrounds, the corresponding files are specified in the respective files in the `samples` directory.
Use `--help` to list all bkg options.
The `-n` argument must be used for background in order to specify the index of the files to run over, i.e. only one can be run at a time.
### Output
Output histograms are sent to `output`, inside of which will be a directory labelled with whatever was provided to `--tag` (`test` by default), which contains all output files, labelled by the sample name. 
For signal, the masses are also appended to the output files. 
For bkg, the sample file index is appended, which corresponds to the line number in the relevant sample file within `samples`. 
### Condor
```
cp -r condor [CONDOR LABEL]
cd [CONDOR LABEL]
./doSubmit.sh [SAMPLE] # TT, TT2018, QCD, QCD2018, ZQQ2018, ST2018, WQQ2018, ZNuNu2018, Diboson2018 
```
The outputs from condor will be send to `output/test/` within the new directory you made, labelled with the sample and indexed by line number in the corresponsing sample file within `../samples`. 
Note: `doSubmit.sh` will simply copy over the `singleStopAnalyzer.py` script from the parent directory, as well as the relevant samples lies from `samples`. 
