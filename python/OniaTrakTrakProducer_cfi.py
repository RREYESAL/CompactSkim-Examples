import FWCore.ParameterSet.Config as cms

PsiPhiProducer = cms.EDProducer('OniaTrakTrakProducer',
    Onia = cms.InputTag('onia2MuMuPAT'),
    Trak = cms.InputTag('cleanPatTrackCands'),
    OniaMassCuts = cms.vdouble(2.946916,3.246916),      # J/psi mass window 3.096916 +/- 0.150
    TrakTrakMassCuts = cms.vdouble(1.004461,1.034461),  # phi mass window 1.019461 +/- .015
    OniaTrakTrakMassCuts = cms.vdouble(5.0,5.7),            # b-hadron mass window
    MassTraks = cms.vdouble(0.493677,0.493677),         # traks masses
    OnlyBest  = cms.bool(True)    
)

PsiPhiFitter = cms.EDProducer('PsiTrakTrakKinematicFit',
    PsiTrakTrak     = cms.InputTag('PsiPhiProducer','OniaTrakTrakCandidates'),
    mass_constraint = cms.double(3.096916),              # J/psi mass in GeV
    product_name    = cms.string('PsiPhiCandidates')
)

PsiPhiSequence = cms.Sequence(PsiPhiProducer*PsiPhiFitter)
