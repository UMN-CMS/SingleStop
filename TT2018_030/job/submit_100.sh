#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc820
export XRD_STREAMTIMEOUT=600
scramv1 project CMSSW CMSSW_10_6_19_patch2
echo "Running ls -alrth:"
ls -alrth
cd CMSSW_10_6_19_patch2/src/
eval `scramv1 runtime -sh`
echo $CMSSW_BASE "is the CMSSW we created on the local worker node"
cd ${_CONDOR_SCRATCH_DIR}
python singleStopAnalyzer.py --sample TT2018 -n 100 --coupling 313
echo "Running pwd:"
pwd
echo "Running ls -alrth:"
ls -alrth
