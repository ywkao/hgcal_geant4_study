import FWCore.ParameterSet.Config as cms

generator = cms.EDProducer("CloseByParticleGunProducer",
    PGunParameters = cms.PSet(PartID = cms.vint32(-13),
        EnMin = cms.double(99.99),
        EnMax = cms.double(100.01),
        RMin = cms.double(35),
        RMax = cms.double(140),
        ZMin = cms.double(320),
        ZMax = cms.double(321),
        Delta = cms.double(10),
        Pointing = cms.bool(True),
        Overlapping = cms.bool(False),
        RandomShoot = cms.bool(False),
        NParticles = cms.int32(1),
        MaxEta = cms.double(3.1),
        MinEta = cms.double(1.3),
        MaxPhi = cms.double(3.14159265359),
        MinPhi = cms.double(-3.14159265359),
    ),
    Verbosity = cms.untracked.int32(0),

    psethack = cms.string('random particles in phi and r windows'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)
