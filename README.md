#CompactSkim-Examples

This package is mean to be run after the BPH CompactSkim.

* **Onia2MuMuRootupler** - Simple rootupler of the Onia2MuMu branch 
* **PsiTrakRootupler**   - Reconstrcut B+ -> J/psi K+, with a kinematical fit and rootuple

* Setup: (has being tested with 74X and 75X, but should run with latest releases, 
  where CompactSkim has run)

```
export SCRAM_ARCH=slc6_amd64_gcc491
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_7_4_5
cd CMSSW_7_4_5/src/
cmsenv
git clone https://github.com/alberto-sanchez/CompactSkim-Examples.git CompactSkim/Examples
scram b
```

* Run: (use your favorite input sample)

```
vi CompactSkim/Examples/test/runOnia2MuMuRootupler.py
cmsRun CompactSkim/Examples/test/runOnia2MuMuRootupler.py
```

Look for other configurations in *test* directory
