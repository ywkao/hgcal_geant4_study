import FWCore.ParameterSet.Config as cms
process = cms.Process("myEvent")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100) 
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),SkipEvent = cms.untracked.vstring('ProductNotFound'))
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
            'file:/eos/user/y/ykao/www/HGCAL_Geant4_project/testbeam_positron_D86_R80To100_E100/step3.root'
            #'file:/eos/user/y/ykao/www/HGCAL_Geant4_project/testbeam_positron_D86_R120To140_E100/step3.root'
            #'file:/eos/user/y/ykao/www/HGCAL_Geant4_project/testbeam_positron_D86_R35To60_E100/step3.root'
            #'file:/eos/user/y/ykao/www/HGCAL_Geant4_project/testbeam_positron_D83_R120To140_E100/step3.root'
            #'file:/eos/user/y/ykao/www/HGCAL_Geant4_project/testbeam_positron_D83_R35To60_E100/step3.root'
            )
        )
process.analyzer = cms.EDAnalyzer('clue_performance',
        layerclusters = cms.untracked.InputTag('hgcalLayerClusters'),
        tracks = cms.untracked.InputTag('generalTracks'),
        photons = cms.untracked.InputTag('photonsFromMultiCl'),
        multis = cms.untracked.InputTag('hgcalMultiClusters'),
        pfcandsticl = cms.untracked.InputTag('pfTICL'),        
        pfcands = cms.untracked.InputTag('particleFlow') ,
        tracksters = cms.untracked.InputTag('ticlTrackstersMerge')
        )

process.TFileService = cms.Service("TFileService", fileName = cms.string("mix_density_D86_R80To100_E100.root"))
#process.TFileService = cms.Service("TFileService", fileName = cms.string("mix_density_D86_R120To140_E100.root"))
#process.TFileService = cms.Service("TFileService", fileName = cms.string("mix_density_D86_R35To60_E100.root"))
#process.TFileService = cms.Service("TFileService", fileName = cms.string("mix_density_D83_R120To140_E100.root"))
#process.TFileService = cms.Service("TFileService", fileName = cms.string("mix_density_D83_R35To60_E100.root"))
process.p = cms.Path(process.analyzer)
