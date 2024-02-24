import FWCore.ParameterSet.Config as cms
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import hgceeDigitizer, hgchefrontDigitizer, hgchebackDigitizer, hfnoseDigitizer

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
process = cms.Process('PROD',Phase2C11I13M9)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )

path = "file:/eos/user/y/ykao/www/HGCAL_Geant4_project"
tag = "D86_R90To130_E100"

process.source = cms.Source("PoolSource",
        #fileNames = cms.untracked.vstring('file:/home/mikumar/t3store3/workarea/CMSSW_9_4_9/src/tmp/step2_1.root')
        #fileNames = cms.untracked.vstring('file:/eos/user/y/ykao/www/HGCAL_Geant4_project/testbeam_positron_D86_R80To100_E100/step2.root')
        fileNames = cms.untracked.vstring( path + '/testbeam_positron_' + tag + '/step2.root')
        #fileNames = cms.untracked.vstring('file:/eos/cms/store/group/dpg_hgcal/tb_hgcal/2022/sps_oct2022/pion_beam_150_320fC/beam_run/run_20221011_130925/beam_run0.root')
        )

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1))

process.prodEE_DigiSim = cms.EDAnalyzer('Calibration',
        simhits = cms.untracked.InputTag("g4SimHits","HGCHitsEE", "SIM"),
        digihits = cms.untracked.InputTag("simHGCalUnsuppressedDigis","EE"),
        generatorSmeared = cms.untracked.InputTag("generatorSmeared"),
        genParticles = cms.untracked.InputTag("genParticles"),
        Detector   = cms.string("HGCalEESensitive"),
        ifNose = cms.untracked.bool(False),
        Verbosity = cms.untracked.int32(0),
        SampleIndx = cms.untracked.int32(2),
        mightGet = cms.optional.untracked.vstring,

        HGCEEConfig = cms.PSet(
            isSiFE = cms.bool(True),
            # adc information
            adcNbits      = hgceeDigitizer.digiCfg.feCfg.adcNbits,
            adcSaturation = hgceeDigitizer.digiCfg.feCfg.adcSaturation_fC,
            #tdc information
            tdcNbits      = hgceeDigitizer.digiCfg.feCfg.tdcNbits,
            tdcSaturation = hgceeDigitizer.digiCfg.feCfg.tdcSaturation_fC,
            tdcOnset      = hgceeDigitizer.digiCfg.feCfg.tdcOnset_fC,
            toaLSB_ns     = hgceeDigitizer.digiCfg.feCfg.toaLSB_ns,
            fCPerMIP      = cms.vdouble(1.25,2.57,3.88) #100um, 200um, 300um
            #fCPerMIP      = cms.vdouble(24.,36.,50.) #120um, 200um, 300um; from Mintu
            ),

        HGCHEFConfig = cms.PSet(
            isSiFE = cms.bool(True),
            # adc information
            adcNbits      = hgchefrontDigitizer.digiCfg.feCfg.adcNbits,
            adcSaturation = hgchefrontDigitizer.digiCfg.feCfg.adcSaturation_fC,
            #tdc information
            tdcNbits      = hgchefrontDigitizer.digiCfg.feCfg.tdcNbits,
            tdcSaturation = hgchefrontDigitizer.digiCfg.feCfg.tdcSaturation_fC,
            tdcOnset      = hgchefrontDigitizer.digiCfg.feCfg.tdcOnset_fC,
            toaLSB_ns     = hgchefrontDigitizer.digiCfg.feCfg.toaLSB_ns,
            fCPerMIP      = cms.vdouble(1.25,2.57,3.88) #100um, 200um, 300um
            ),

        HGCHEBConfig = cms.PSet(
            isSiFE  = cms.bool(True),
            # adc information
            adcNbits      = hgchebackDigitizer.digiCfg.feCfg.adcNbits,
            adcSaturation = hgchebackDigitizer.digiCfg.feCfg.adcSaturation_fC,
            # tdc information
            tdcNbits      = hgchebackDigitizer.digiCfg.feCfg.tdcNbits,
            tdcSaturation = hgchebackDigitizer.digiCfg.feCfg.tdcSaturation_fC,
            tdcOnset      = hgchebackDigitizer.digiCfg.feCfg.tdcOnset_fC,
            toaLSB_ns     = hgchebackDigitizer.digiCfg.feCfg.toaLSB_ns,
            fCPerMIP      = cms.vdouble(1.0,1.0,1.0) #dummy values, it's scintillator
            ),

        HGCHFNoseConfig = cms.PSet(
                isSiFE = cms.bool(False),
                # adc information
                adcNbits      = hfnoseDigitizer.digiCfg.feCfg.adcNbits,
                adcSaturation = hfnoseDigitizer.digiCfg.feCfg.adcSaturation_fC,
                #tdc information
                tdcNbits      = hfnoseDigitizer.digiCfg.feCfg.tdcNbits,
                tdcSaturation = hfnoseDigitizer.digiCfg.feCfg.tdcSaturation_fC,
                tdcOnset      = hfnoseDigitizer.digiCfg.feCfg.tdcOnset_fC,
                toaLSB_ns     = hfnoseDigitizer.digiCfg.feCfg.toaLSB_ns,
                fCPerMIP      = cms.vdouble(1.25,2.57,3.88) #100um, 200um, 300um
                ),

        algo = cms.string("HGCalUncalibRecHitWorkerWeights")
        )

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('rootfiles/geantoutput_' + tag + '.root')
        )

process.p = cms.Path(process.prodEE_DigiSim)
