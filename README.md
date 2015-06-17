This package is mean to be run after the BPH CompactSkim.

Onia2MuMuRootupler - Simple rootupler of the Onia2MuMu branch 


* Setup:

export SCRAM_ARCH=slc6_amd64_gcc491

source /cvmfs/cms.cern.ch/cmsset_default.sh

cmsrel CMSSW_7_4_5

cd CMSSW_7_4_5/src/

cmsenv

git clone https://github.com/alberto-sanchez/CompactSkim-Examples.git CompactSkim/Examples

scram b

* Run:

vi CompactSkim/Examples/test/runOnia2MuMuRootupler.py

cmsRun CompactSkim/Examples/test/runOnia2MuMuRootupler.py
