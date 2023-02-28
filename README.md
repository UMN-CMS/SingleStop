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
### Locally
```
python singleStopAnalyzer.py --sample [SAMPLE] --tag [TAG] -n [BKG FILE INDEX] --points [MSTOP_MCHI]
```
#### Input
For running over signal, the desired signal points can be specified with the `--points` argument (`MSTOP_MCHI`), comma-separated. 
All points will be used unless specified.
For backgrounds, the corresponding files are specified in their respective files in the `samples` directory.
Use `--help` to list all bkg options.
The `-n` argument must be used for background in order to specify the index of the files to run over, i.e. only one file can be run at a time locally. 
See the "Condor" section for details on how  to submit jobs to run over multiple files.
#### Output
Output histograms are sent to `output`, inside of which will be a directory labelled with whatever was provided to `--tag` (`test` by default), which contains all output files, labelled by the sample name. 
For signal, the masses are also appended to the output files. 
For bkg, the sample file index is appended, which corresponds to the line number in the relevant sample file within `samples`. 
#### Signal example
```
python singleStopAnalyzer.py --sample signal --tag signal_23-02-27 --points 2000_900,1000_400
```
#### Background example
```
python singleStopAnalyzer.py --sample TT2018 --tag TT2018_23-02-27 -n 206
```
### Condor
To run over multiple files (those specified in `samples`), follow the instricutions below. 
```
cp -r condor [CONDOR LABEL]
cd [CONDOR LABEL]
./doSubmit.sh [SAMPLE] # TT, TT2018, QCD, QCD2018, ZQQ2018, ST2018, WQQ2018, ZNuNu2018, Diboson2018 
```
The outputs from condor will be send to `output/test/` within the new directory you made, labelled with the sample and indexed by line number in the corresponsing sample file within `../samples`. 
Note: `doSubmit.sh` will simply copy over the `singleStopAnalyzer.py` script from the parent directory, as well as the relevant file list from `samples`. 
### Scaling MC
The analyzer does not apply any MC weights, i.e. all histograms are filled with a weight of 1.
The `scaleMC.py` script can be run over the output of the analyzer to scale all histograms to the full Run 2 luminosity (any 2018-only samples are automatically scaled to the full Run 2 luminosity).
This script also calls `hadd` to combine all outputs, e.g. all those from condor jobs in the `--input` directory. 
The output will be a single ROOT file labelled with the sample name.
```
python scaleMC.py --input [PATH TO ANALYZER ROOT OUTPUTS] --output [OUTPUT DIRECTORY] --sample [SAMPLE]
```
`[SAMPLE]` is the same as that used by the analyzer. 
Use `--help` to see all options.
Note: `scaleMC.py` extracts all histograms from the `TDirectory` named `plots` and write them directly to the output ROOT file without the `plots` directory.
