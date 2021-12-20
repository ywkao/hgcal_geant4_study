import FWCore.ParameterSet.Config as cms


# from Configuration.Eras.Era_Phase2C11_cff import Phase2C11
# process = cms.Process('PROD',Phase2C11)

#process.load("SimGeneral.HepPDTESSource.pdt_cfi")
#process.load("Configuration.Geometry.GeometryExtended2026D76Reco_cff")

#from Configuration.Eras.Modifier_phase2_hgcalV12_cff import phase2_hgcalV12
#process = cms.Process('PROD',phase2_hgcalV12)

# from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
# process = cms.Process('PROD',Phase2C9)

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
process = cms.Process('PROD',Phase2C11I13M9)

#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('Configuration.Geometry.GeometryExtended2026D83Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')

#process.load("IOMC.RandomEngine.IOMC_cff")
#process.RandomNumberGeneratorService.generator.initialSeed = 456789

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/home/mikumar/t3store3/workarea/HGCAL_Validation/CMSSW_12_2_X_2021-11-29-2300/src/Muon_gun/Mip_peak/step2.root')

)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))


process.prodEE = cms.EDAnalyzer('simHit',
                             simtrack = cms.untracked.InputTag("g4SimHits"),
                             simhits = cms.untracked.InputTag("g4SimHits","HGCHitsEE", "SIM"),
                             #Detector   = cms.string("HGCalEESensitive"),
                         )


process.prodHEF = process.prodEE.clone(
    simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEfront", "SIM"),
    #Detector   = cms.string("HGCalHESiliconSensitive"),
)

process.prodHEB = process.prodHEF.clone(
    simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEback", "SIM"),
    #Detector   = cms.string("HGCalHEScintillatorSensitive"),
)

#process.Tracer = cms.Service("Tracer")

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('geantoutput.root')
 )

process.p = cms.Path(process.prodEE*process.prodHEF*process.prodHEB)
