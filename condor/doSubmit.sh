#!/bin/bash

if [ -z "$1" ]
  then
    echo "ERROR: No sample argument supplied. Usage: ./doSubmit.sh [SAMPLE]"
    exit 1
fi

sample=$1

case $sample in
  "QCD") sampleFile="QCDBEnriched.txt";;
  "QCD2018") sampleFile="QCDBEnriched2018.txt";;
  "TT") sampleFile="TTToHadronic.txt";;
  "TT2018") sampleFile="TTToHadronic2018.txt";;
  "ZQQ2018") sampleFile="ZJetsToQQ2018.txt";;
  "ST2018") sampleFile="STHadronic2018.txt";;
  "WQQ2018") sampleFile="WJetsToQQ2018.txt";;
  "ZNuNu2018") sampleFile="ZJetsToNuNu2018.txt";;
  "Diboson2018") sampleFile="Diboson2018.txt";;
  *) echo "ERROR: Invalid sample argument.";exit 1;;
esac

rm -rf "job/$sample" "out/$sample" "err/$sample" "log/$sample" samples/$sample/samples
mkdir -p "job/$sample" "out/$sample" "err/$sample" "log/$sample" samples/$sample/samples

cp ../samples/"$sampleFile" samples/$sample/samples
cp ../singleStopAnalyzer.py .
cp ../match_algos.py .
cp ../mass_reco.py .

sed -i "2,4s/PhysicsTools.NanoAODTools.postprocessing.//g" singleStopAnalyzer.py

i=1
while read file; do
#if (($i > 2)); then continue; fi

jobfilename=job/"$sample"/submit_${i}.sh

cat << EOT >> ${jobfilename} 
#!/bin/bash
echo "Starting job on " \`date\` #Date/time of start of job
echo "Running on: \`uname -a\`" #Condor job is running on this node
echo "System software: \`cat /etc/redhat-release\`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc820
export XRD_STREAMTIMEOUT=600
scramv1 project CMSSW CMSSW_10_6_19_patch2
echo "Running ls -alrth:"
ls -alrth
cd CMSSW_10_6_19_patch2/src/
eval \`scramv1 runtime -sh\`
echo \$CMSSW_BASE "is the CMSSW we created on the local worker node"
cd \${_CONDOR_SCRATCH_DIR}
python singleStopAnalyzer.py --sample $sample -n $i
echo "Running pwd:"
pwd
echo "Running ls -alrth:"
ls -alrth
EOT

i=$((i+1))
done <samples/$sample/samples/"$sampleFile"

submitfilename=job/"$sample"/submit.sub
cat << EOT >> $submitfilename
Executable	= job/$sample/submit_\$(ijobname).sh
Output		= out/$sample/submit_\$(ijobname).out
Error		= err/$sample/submit_\$(ijobname).err
Log		= log/$sample/submit_\$(ijobname).log
transfer_input_files = samples/$sample/samples,framework,singleStopAnalyzer.py,match_algos.py,mass_reco.py
transfer_output_files = output
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
universe	= vanilla
queue ijobname from (
EOT

i=$((i-1))
for j in $(seq 1 ${i})
do
echo "${j}" >> $submitfilename
done
echo ")" >> $submitfilename

condor_submit $submitfilename
