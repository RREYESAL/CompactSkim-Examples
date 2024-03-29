import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:infilename')
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('file:outfilename'),
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load('CompactSkim.Examples.Onia2MuMuRootupler_cfi')
process.p = cms.Path(process.rootuple)
