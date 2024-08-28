# EGEfficiency

To set up the repository:
```
cmsrel CMSSW_14_0_14
cd CMSSW_14_0_14/src
cmsenv
voms-proxy-init --voms cms
git cms-init
git clone git@github.com:RSalvatico/EGEfficiency.git EGEfficiency
scram b -j 8
```

To run the EfficiencyCalculator, do:
```
cd EGEfficiency/TrigTools/test
cmsRun run_EfficiencyCalculator.py
```

This produces a set of numerator, denominator, and occupancy histograms. The content of run_EfficiencyCalculator.py can easily be changed to pick a different list of input files.

If you want to produce MINIAOD out of RAW-RECO input and run the EfficiencyCalculator on those, you can use the following `cmsDriver` command:
```
cmsDriver.py stepMINI -s PAT --conditions auto:phase1_2024_realistic --datatier MINIAOD -n -1 --eventcontent MINIAOD --python_filename makeMini_cfg.py --geometry DB:Extended --era Run3_2024 --filein file:/afs/cern.ch/work/r/rverma/public/cms-egamma-hlt/CMSHLT-3313/test_RAW2DIGI_L1REPACK_HLT.root --fileout file:stepMINI.root --hltProcess MYHLT --no_exec
```
If you remove the option `--hltProcess MYHLT` from the configuration, the standard `HLT` collections will be used to make PAT objects, such as `pat::TriggerObjectStandAlone`.

There is a script to run `makeMini_cfg.py` on HTCondor. From `test/`, you can do:
```
python3 prepareCondorJobs.py
```
This will create `Jobs/Job_*` directories, one for every input file. To send the jobs to the scheduler, do:
```
./submit_all.jobb
```

There is also a customizable script to make efficiencies out of the histos produced by EfficiencyCalculator, and that is in:
```
TrigTools/test/compareEfficiencies.py
```