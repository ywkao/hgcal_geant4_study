import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
process = cms.Process('PROD',Phase2C11I13M9)

process.load('FWCore.MessageService.MessageLogger_cfi')


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/m/mikumar/work/private/HGCAL_Validation/CMSSW_12_2_X_2021-12-10-2300/src/step2.root')#/home/mikumar/t3store3/workarea/HGCAL_Validation/CMSSW_12_2_X_2021-11-29-2300/src/Muon_gun/Events_10k/step2.root')#Mip_peak/

)

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
  				mightGet = cms.optional.untracked.vstring
                         )

"""process.prodHEF_DigiSim = process.prodEE_DigiSim.clone(
    simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEfront", "SIM"),
    Detector   = cms.string("HGCalHESiliconSensitive"),
)

process.prodHEB_DigiSim = process.prodEE_DigiSim.clone(
    simhits = cms.untracked.InputTag("g4SimHits","HGCHitsHEback", "SIM"),
    Detector   = cms.string("HGCalHEScintillatorSensitive"),
)"""


process.TFileService = cms.Service("TFileService",
     fileName = cms.string('geantoutput.root')
 )

process.p = cms.Path(process.prodEE_DigiSim)#*process.prodHEF_DigiSim*process.prodHEB_DigiSim)
