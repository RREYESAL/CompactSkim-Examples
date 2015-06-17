
import FWCore.ParameterSet.Config as cms
process = cms.Process('PSIK')

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc') 

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:inputfile')
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load("CompactSkim.Examples.PsiTrakRootupler_cfi")

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('rootuple-PsiTrakRootupler.root'),
)

process.p = cms.Path(process.PsiTrakSequence)
