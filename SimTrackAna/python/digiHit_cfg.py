import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
process = cms.Process('PROD',Phase2C11I13M9)

process.load('FWCore.MessageService.MessageLogger_cfi')


process.source = cms.Source("PoolSource",
        #fileNames = cms.untracked.vstring('file:/home/mikumar/t3store3/workarea/CMSSW_9_4_9/src/tmp/step2_1.root')
        fileNames = cms.untracked.vstring('file:/eos/user/y/ykao/www/HGCAL_Geant4_project/testbeam_positron_D86_R80To100_E100/step2.root')
        )

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10))


process.prodEE_DigiSim = cms.EDAnalyzer('DigiSim',
        #simtrack = cms.untracked.InputTag("g4SimHits"),
        #digihits = cms.untracked.InputTag("hgcalDigis","EE"),# "HLT"),
        simhits = cms.untracked.InputTag("g4SimHits","HGCHitsEE", "SIM"),
        digihits = cms.untracked.InputTag("simHGCalUnsuppressedDigis","EE"),# "HLT"),
        Detector   = cms.string("HGCalEESensitive"),
        ifNose = cms.untracked.bool(False),
        Verbosity = cms.untracked.int32(0),
        SampleIndx = cms.untracked.int32(2),
        mightGet = cms.optional.untracked.vstring,

        #HGCEEdigiCollection = cms.InputTag('hgcalDigis:EE'),
        #HGCHEFdigiCollection = cms.InputTag('hgcalDigis:HEfront'),
        #HGCHEBdigiCollection = cms.InputTag('hgcalDigis:HEback'),
        #HGCHFNosedigiCollection = cms.InputTag('hfnoseDigis:HFNose'),

        #HGCEEdigiCollection = cms.InputTag('simHGCalUnsuppressedDigis:EE'),
        #HGCHEFdigiCollection = cms.InputTag('simHGCalUnsuppressedDigis:HEfront'),
        #HGCHEBdigiCollection = cms.InputTag('simHGCalUnsuppressedDigis:HEback'),
        algo = cms.string("HGCalUncalibRecHitWorkerWeights")
        )

process.prodHEF_DigiSim = process.prodEE_DigiSim.clone(
        simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEfront", "SIM"),
        Detector   = cms.string("HGCalHESiliconSensitive"),
        )

process.prodHEB_DigiSim = process.prodEE_DigiSim.clone(
        simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEback", "SIM"),
        Detector   = cms.string("HGCalHEScintillatorSensitive"),
        )


process.TFileService = cms.Service("TFileService",
        fileName = cms.string('geantoutput.root')
        )

process.p = cms.Path(process.prodEE_DigiSim*process.prodHEF_DigiSim*process.prodHEB_DigiSim)
