import FWCore.ParameterSet.Config as cms

PsiTrakProducer = cms.EDProducer('OniaTrakProducer',
    Onia = cms.InputTag('onia2MuMuPAT'),
    Trak = cms.InputTag('cleanPatTrackCands'),
    OniaMassCuts = cms.vdouble(2.947,3.247),
    OniaTrakMassCuts = cms.vdouble(5.0,5.7),
    OnlyBest = cms.bool(True)    
)

PsiTrakFitter = cms.EDProducer('PsiTrakKinematicFit',
    PsiTrak         = cms.InputTag('PsiTrakProducer','OniaTrakCandidates'),
    mass_constraint = cms.double(3.0969),              # J/psi mass in GeV
    product_name    = cms.string('PsiTrakCandidates')
)

rootuple = cms.EDAnalyzer('PsiTrakRootupler',
    oniat_cand = cms.InputTag("PsiTrakProducer","OniaTrakCandidates"),
    oniat_rf_cand = cms.InputTag("PsiTrakFitter","PsiTrakCandidates"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    oniat_pdgid = cms.uint32(521),
    onia_pdgid = cms.uint32(443),
    trak_pdgid = cms.uint32(321),
    isMC = cms.bool(True),
    OnlyBest = cms.bool(True)                          
)

PsiTrakSequence = cms.Sequence(PsiTrakProducer*PsiTrakFitter*rootuple)
