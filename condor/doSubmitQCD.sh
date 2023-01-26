#!/bin/bash
rm -rf job out err log samples
mkdir -p job out err log samples

cp ../samples/TTToHadronic.txt ../samples/QCDBEnriched.txt samples

i=1
while read file; do
#if (($i > 2)); then continue; fi

jobfilename=job/submit_${i}.sh

cat << EOT >> ${jobfilename} 
#!/bin/bash
echo "Starting job on " \`date\` #Date/time of start of job
echo "Running on: \`uname -a\`" #Condor job is running on this node
echo "System software: \`cat /etc/redhat-release\`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc820
scramv1 project CMSSW CMSSW_10_6_19_patch2
echo "Running ls -alrth:"
ls -alrth
cd CMSSW_10_6_19_patch2/src/
eval \`scramv1 runtime -sh\`
echo \$CMSSW_BASE "is the CMSSW we created on the local worker node"
cd \${_CONDOR_SCRATCH_DIR}
python singleStopAnalyzer.py --sample QCD -n $i
echo "Running pwd:"
pwd
echo "Running ls -alrth:"
ls -alrth
EOT

i=$((i+1))
done <samples/QCDBEnriched.txt

submitfilename=job/submit.sub
cat << EOT >> $submitfilename
Executable	= job/submit_\$(ijobname).sh
Output		= out/submit_\$(ijobname).out
Error		= err/submit_\$(ijobname).err
Log		= log/submit_\$(ijobname).log
transfer_input_files = samples,framework,singleStopAnalyzer.py
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
