#!/bin/bash                                                                                                                                                 
echo "Starting job on " `date`
echo "Running on: `uname -a`"
echo "System software: `cat /etc/redhat-release`"
source /cvmfs/cms.cern.ch/cmsset_default.sh

OUTDIR=/eos/cms/store/group/phys_egamma/ec/rsalvatico/MiniHoE/test/

[ ! -d "${OUTDIR}" ] && mkdir -p "${OUTDIR}"

echo "======="
ls
echo "======"

echo $PWD
eval `scramv1 project CMSSW CMSSW_14_0_14`
cd CMSSW_14_0_14/src/
# set cmssw environment                                                                                                                                     
eval `scram runtime -sh`
cd -
echo "========================="
echo "==> List all files..."
ls
echo "+=============================="
echo "==> Running the desired process"

cmsRun ${3}
echo "List all root files = "
ls *.root
cp stepMINI.root stepMINI_${1}_${2}.root
echo "========================="
echo "==> List all files..."
ls *.root
echo "+=============================="
echo "xrdcp output for condor"
cp -f stepMINI_${1}_${2}.root ${OUTDIR}/stepMINI_${1}_${2}.root
rm stepMINI_${1}_${2}.root
echo "========================="
echo "========================="
echo "==> List all files..."
ls *.root
echo "+=============================="
date
echo "+=============================="
