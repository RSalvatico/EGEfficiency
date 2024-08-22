# EGEfficiency

To set up the repository:
```
cmsrel CMSSW_13_0_10
cd CMSSW_13_0_10/src
cmsenv
voms-proxy-init --voms cms
git cms-init
git clone git@github.com:RSalvatico/EGEfficiency.git EGTools
scram b -j 8
```

To run the EfficiencyCalculator, do:
```
cd EGTools/TrigTools/test
cmsRun run_EfficiencyCalculator.py
```