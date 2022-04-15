#include <memory>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <iostream>
// user include files {{{
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "CoralBase/Exception.h"

#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TMath.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"

//#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
//#include "DQMServices/Core/interface/DQMStore.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToModule.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToROC.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalUncalibRecHitWorkerFactory.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalUncalibRecHitRecWeightsAlgo.h"
//}}}
#include "../interface/toolbox.h"
#include <TNtuple.h>
#include <TVector3.h>
namespace tb = toolbox;

void configureIt(const edm::ParameterSet& conf, HGCalUncalibRecHitRecWeightsAlgo<HGCalDataFrame>& maker) //{{{
{
    //https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/HGCalRecProducers/plugins/HGCalUncalibRecHitWorkerWeights.cc#L10-L58
    constexpr char isSiFE[] = "isSiFE";
    constexpr char adcNbits[] = "adcNbits";
    constexpr char adcSaturation[] = "adcSaturation";
    constexpr char tdcNbits[] = "tdcNbits";
    constexpr char tdcSaturation[] = "tdcSaturation";
    constexpr char tdcOnset[] = "tdcOnset";
    constexpr char toaLSB_ns[] = "toaLSB_ns";
    constexpr char fCPerMIP[] = "fCPerMIP";

    if (conf.exists(isSiFE)) {
        maker.set_isSiFESim(conf.getParameter<bool>(isSiFE));
    } else {
        maker.set_isSiFESim(false);
    }

    if (conf.exists(adcNbits)) {
        uint32_t nBits = conf.getParameter<uint32_t>(adcNbits);
        double saturation = conf.getParameter<double>(adcSaturation);
        float adcLSB = saturation / pow(2., nBits);
        maker.set_ADCLSB(adcLSB);
    } else {
        maker.set_ADCLSB(-1.);
    }

    if (conf.exists(tdcNbits)) {
        uint32_t nBits = conf.getParameter<uint32_t>(tdcNbits);
        double saturation = conf.getParameter<double>(tdcSaturation);
        double onset = conf.getParameter<double>(tdcOnset);  // in fC
        float tdcLSB = saturation / pow(2., nBits);
        maker.set_TDCLSB(tdcLSB);
        maker.set_tdcOnsetfC(onset);
    } else {
        maker.set_TDCLSB(-1.);
        maker.set_tdcOnsetfC(-1.);
    }

    if (conf.exists(toaLSB_ns)) {
        maker.set_toaLSBToNS(conf.getParameter<double>(toaLSB_ns));
    } else {
        maker.set_toaLSBToNS(-1.);
    }

    if (conf.exists(fCPerMIP)) {
        maker.set_fCPerMIP(conf.getParameter<std::vector<double> >(fCPerMIP));
    } else {
        maker.set_fCPerMIP(std::vector<double>({1.0}));
    }
} //}}}
class DigiSim : public edm::one::EDAnalyzer<edm::one::SharedResources> { //{{{
    public:
        //Implemented following Validation/HGCalValidation/plugins/HGCalSimHitValidation.cc

        explicit DigiSim(const edm::ParameterSet&);
        ~DigiSim();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
        struct energysum {
            energysum() {
                etotal = 0;
                for (int i = 0; i < 6; ++i)
                    eTime[i] = 0.;
            }
            double eTime[6], etotal;
        };
        struct adcinfo {
            adcinfo() {  adc = 0;}
            uint32_t adc;
        };  
        /*struct waferinfo {
          waferinfo() {      
          layer = u = v = type = -999;
          }
          int layer, u, v, type;
          };*/

        struct hitsinfo {
            hitsinfo() {
                u_cor = v_cor = type = layer = 0;
                hitid = nhits = 0;
            }
            int u_cor, v_cor, type, layer;
            unsigned int hitid, nhits;
        };

        struct digisinfo {
            digisinfo() {
                u_cor = v_cor = type = layer = 0;
                x_pos = y_pos = z_pos = eta = phi = 0.;
                hitid = ndigis = 0;
            }
            int u_cor, v_cor, type, layer;
            float x_pos, y_pos, z_pos, eta, phi;
            unsigned int hitid, ndigis;
        };

        struct myDigis {
            digisinfo dinfo;
            adcinfo ainfo;
            double amplitude;
        };

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        //const std::string name;

        const std::string nameDetector_; 
        const bool ifNose_;
        const int verbosity_, SampleIndx_;
        // ----------member data ---------------------------
        const edm::ESGetToken<HGCalGeometry, IdealGeometryRecord> tok_hgcalg_;
        int firstLayer_; 
        edm::EDGetTokenT<edm::PCaloHitContainer> tSimCaloHitContainer; 
        //edm::EDGetTokenT<HGCalDigiCollection> eeDigiCollection_;     // collection of HGCEE digis
        //edm::EDGetTokenT<HGCalDigiCollection> hefDigiCollection_;    // collection of HGCHEF digis
        //edm::EDGetTokenT<HGCalDigiCollection> hebDigiCollection_;    // collection of HGCHEB digis
        //edm::EDGetTokenT<HGCalDigiCollection> hfnoseDigiCollection_; // collection of HGCHFNose digis
        //std::unique_ptr<HGCalUncalibRecHitWorkerBaseClass> worker_;

        edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
        hgcal::RecHitTools rhtools_;
        edm::EDGetToken digiSource_;
        //edm::ConsumesCollector iC;
        edm::ESGetToken<HGCalGeometry, IdealGeometryRecord> ee_geometry_token_;
        HGCalUncalibRecHitRecWeightsAlgo<HGCalDataFrame> uncalibMaker_ee_;

        TH1D *hELossEE;TH1D *hELossEEF;TH1D *hELossEECN;TH1D *hELossEECK;
        TH1D *hELossHEF;TH1D *hELossHEFF;TH1D *hELossHEFCN;TH1D *hELossHEFCK;
        
        TH1D *hEta;
        TH1D *hPhi;

        // hit
        std::vector<TH1D*> vechist;   
        TNtuple *nt_hit_position;
        TNtuple *nt_total_[26];
        TNtuple *nt_120mum_[26];
        TNtuple *nt_200mum_[26];
        TNtuple *nt_300mum_[26];
        TH1D *ADC_total_[26];
        TH1D *ADC_120mum_[26];
        TH1D *ADC_200mum_[26];
        TH1D *ADC_300mum_[26];
        TH1D *MIP_total_[26];
        TH1D *MIP_120mum_[26];
        TH1D *MIP_200mum_[26];
        TH1D *MIP_300mum_[26];
        TH1D *SIM_total_[26];
        TH1D *SIM_120mum_[26];
        TH1D *SIM_200mum_[26];
        TH1D *SIM_300mum_[26];
        TProfile *adc_sim_total_[26];
        TProfile *adc_sim_120mum_[26];
        TProfile *adc_sim_200mum_[26];
        TProfile *adc_sim_300mum_[26];
        TProfile *adc_mip_total_[26];
        TProfile *adc_mip_120mum_[26];
        TProfile *adc_mip_200mum_[26];
        TProfile *adc_mip_300mum_[26];
        TProfile *mip_sim_total_[26];
        TProfile *mip_sim_120mum_[26];
        TProfile *mip_sim_200mum_[26];
        TProfile *mip_sim_300mum_[26];

        // total energy
        TH1D *total_ADC_total_[26];
        TH1D *total_ADC_120mum_[26];
        TH1D *total_ADC_200mum_[26];
        TH1D *total_ADC_300mum_[26];
        TH1D *total_MIP_total_[26];
        TH1D *total_MIP_coarse_[26];
        TH1D *total_MIP_fine_[26];
        TH1D *total_MIP_120mum_[26];
        TH1D *total_MIP_200mum_[26];
        TH1D *total_MIP_300mum_[26];
        TH1D *total_SIM_total_[26];
        TH1D *total_SIM_120mum_[26];
        TH1D *total_SIM_200mum_[26];
        TH1D *total_SIM_300mum_[26];

        // multiplicity
        TH1D *multiplicity_digis_total_[26];
        TH1D *multiplicity_digis_120mum_[26];
        TH1D *multiplicity_digis_200mum_[26];
        TH1D *multiplicity_digis_300mum_[26];
        TH1D *multiplicity_simhits_total_[26];
        TH1D *multiplicity_simhits_coarse_[26];
        TH1D *multiplicity_simhits_fine_[26];
        TH1D *multiplicity_simhits_120mum_[26];
        TH1D *multiplicity_simhits_200mum_[26];
        TH1D *multiplicity_simhits_300mum_[26];

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
        edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif

};
// }}}
DigiSim::DigiSim(const edm::ParameterSet& iconfig) : //{{{
        //auto temp = iConfig.getUntrackedParameter<edm::InputTag>("digihits");
        nameDetector_(iconfig.getParameter<std::string>("Detector")),
        ifNose_(iconfig.getUntrackedParameter<bool>("ifNose")),
        verbosity_(iconfig.getUntrackedParameter<int>("Verbosity", 0)),
        SampleIndx_(iconfig.getUntrackedParameter<int>("SampleIndx", 0)),
        tok_hgcalg_(esConsumes<HGCalGeometry, IdealGeometryRecord>(edm::ESInputTag{"", nameDetector_})), 
        firstLayer_(1), 
        tSimCaloHitContainer(consumes<edm::PCaloHitContainer>(iconfig.getUntrackedParameter<edm::InputTag>("simhits"))),
        //eeDigiCollection_(consumes<HGCalDigiCollection>(iconfig.getParameter<edm::InputTag>("HGCEEdigiCollection"))),
        //hefDigiCollection_(consumes<HGCalDigiCollection>(iconfig.getParameter<edm::InputTag>("HGCHEFdigiCollection"))),
        //hebDigiCollection_(consumes<HGCalDigiCollection>(iconfig.getParameter<edm::InputTag>("HGCHEBdigiCollection"))),
        //hfnoseDigiCollection_(consumes<HGCalDigiCollection>(iconfig.getParameter<edm::InputTag>("HGCHFNosedigiCollection"))),
        //worker_{HGCalUncalibRecHitWorkerFactory::get()->create(
        //        iconfig.getParameter<std::string>("algo"), iconfig, consumesCollector())},
        ee_geometry_token_(consumesCollector().esConsumes(edm::ESInputTag("", "HGCalEESensitive")))
{
    caloGeomToken_ = esConsumes<CaloGeometry, CaloGeometryRecord>();

    auto temp = iconfig.getUntrackedParameter<edm::InputTag>("digihits");
    if ((nameDetector_ == "HGCalEESensitive") || (nameDetector_ == "HGCalHESiliconSensitive") ||
            (nameDetector_ == "HGCalHEScintillatorSensitive") || (nameDetector_ == "HGCalHFNoseSensitive")) {
        digiSource_ = consumes<HGCalDigiCollection>(temp);
        //digiSource_=consumes<HGCalDigiCollection>(iconfig.getUntrackedParameter<edm::InputTag>("digihits"));
        //tSimCaloHitContainer=consumes<edm::PCaloHitContainer>(iconfig.getUntrackedParameter<edm::InputTag>("simhits"));
    } else {
        throw cms::Exception("BadHGCDigiSource") << "HGCal DetectorName given as " << nameDetector_ << " must be: "
            << "\"HGCalEESensitive\", \"HGCalHESiliconSensitive\", or "
            << "\"HGCalHEScintillatorSensitive\", \"HGCalHFNoseSensitive\"!";
    }

    const edm::ParameterSet& ee_cfg = iconfig.getParameterSet("HGCEEConfig");
    configureIt(ee_cfg, uncalibMaker_ee_);

    // tSimCaloHitContainer(consumes<edm::PCaloHitContainer>(iconfig.getUntrackedParameter<edm::InputTag>("simhits")))
    //now do what ever initialization is needed
    //name = iconfig.getParameter<std::string>("Detector");

    usesResource("TFileService");
    edm::Service<TFileService> fs; 
    hELossEE    = fs->make<TH1D>("hELossEE"    , "hELossEE"    , 1000 , 0. , 1000.);
    hELossEEF   = fs->make<TH1D>("hELossEEF"   , "hELossEEF"   , 1000 , 0. , 1000.);
    hELossEECN  = fs->make<TH1D>("hELossEECN"  , "hELossEECN"  , 1000 , 0. , 1000.);
    hELossEECK  = fs->make<TH1D>("hELossEECK"  , "hELossEECK"  , 1000 , 0. , 1000.);
    hELossHEF   = fs->make<TH1D>("hELossHEF"   , "hELossHEF"   , 1000 , 0. , 1000.);
    hELossHEFF  = fs->make<TH1D>("hELossHEFF"  , "hELossHEFF"  , 1000 , 0. , 1000.);
    hELossHEFCN = fs->make<TH1D>("hELossHEFCN" , "hELossHEFCN" , 1000 , 0. , 1000.);
    hELossHEFCK = fs->make<TH1D>("hELossHEFCK" , "hELossHEFCK" , 1000 , 0. , 1000.);
    nt_hit_position = fs->make<TNtuple>("nt_hit_position","nt_hit_position", "r:z:is_Silicon_w120:is_Silicon_w200:is_Silicon_w300:is_Scintillator");
    //hEta = fs->make<TH1D>("hEta" , "hEta" , 20 , -5. , 5.);
    hEta = fs->make<TH1D>("hEta" , "hEta" , 20 ,  1. , 3.);
    hPhi = fs->make<TH1D>("hPhi" , "hPhi" , 20 , -3. , 3.);
    std::ostringstream hnamestr (std::ostringstream::ate);
    for(int i=0;i<26;i++) {
        // total energy & multiplicity
        tb::set_string(hnamestr, "total_ADC_total_", i+1);
        total_ADC_total_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 5000.);
        tb::set_string(hnamestr, "total_MIP_total_", i+1);
        total_MIP_total_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 5000.);
        tb::set_string(hnamestr, "total_MIP_coarse_", i+1);
        total_MIP_coarse_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 5000.);
        tb::set_string(hnamestr, "total_MIP_fine_", i+1);
        total_MIP_fine_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 5000.);
        tb::set_string(hnamestr, "total_SIM_total_", i+1);
        total_SIM_total_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 5000.);
        tb::set_string(hnamestr, "multiplicity_digis_total_", i+1);
        multiplicity_digis_total_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 2000.);
        tb::set_string(hnamestr, "multiplicity_simhits_total_", i+1);
        multiplicity_simhits_total_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 200.);
        tb::set_string(hnamestr, "multiplicity_simhits_coarse_", i+1);
        multiplicity_simhits_coarse_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 200.);
        tb::set_string(hnamestr, "multiplicity_simhits_fine_", i+1);
        multiplicity_simhits_fine_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 200.);

        tb::set_string(hnamestr, "total_ADC_120mum_", i+1);
        total_ADC_120mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 5000.);
        tb::set_string(hnamestr, "total_MIP_120mum_", i+1);
        total_MIP_120mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 5000.);
        tb::set_string(hnamestr, "total_SIM_120mum_", i+1);
        total_SIM_120mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 5000.);
        tb::set_string(hnamestr, "multiplicity_digis_120mum_", i+1);
        multiplicity_digis_120mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 200.);
        tb::set_string(hnamestr, "multiplicity_simhits_120mum_", i+1);
        multiplicity_simhits_120mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 200.);

        tb::set_string(hnamestr, "total_ADC_200mum_", i+1);
        total_ADC_200mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 5000.);
        tb::set_string(hnamestr, "total_MIP_200mum_", i+1);
        total_MIP_200mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 5000.);
        tb::set_string(hnamestr, "total_SIM_200mum_", i+1);
        total_SIM_200mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 5000.);
        tb::set_string(hnamestr, "multiplicity_digis_200mum_", i+1);
        multiplicity_digis_200mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 200.);
        tb::set_string(hnamestr, "multiplicity_simhits_200mum_", i+1);
        multiplicity_simhits_200mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 200.);

        tb::set_string(hnamestr, "total_ADC_300mum_", i+1);
        total_ADC_300mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 5000.);
        tb::set_string(hnamestr, "total_MIP_300mum_", i+1);
        total_MIP_300mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 5000.);
        tb::set_string(hnamestr, "total_SIM_300mum_", i+1);
        total_SIM_300mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 5000.);
        tb::set_string(hnamestr, "multiplicity_digis_300mum_", i+1);
        multiplicity_digis_300mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 200.);
        tb::set_string(hnamestr, "multiplicity_simhits_300mum_", i+1);
        multiplicity_simhits_300mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 200.);

        // individual hits
        tb::set_string(hnamestr, "ADC_total_layer_", i+1);
        ADC_total_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 80, 0, 800.);
        tb::set_string(hnamestr, "MIP_total_layer_", i+1);
        MIP_total_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 150, 0, 150.);
        tb::set_string(hnamestr, "SIM_total_layer_", i+1);
        SIM_total_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 80, 0, 800.);
        tb::set_string(hnamestr, "ADC_SimhitE_total_layer_", i+1);
        adc_sim_total_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 800., 0, 800., 0., 800.);
        tb::set_string(hnamestr, "ADC_MIP_total_layer_", i+1);
        adc_mip_total_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 800., 0, 800., 0., 150.);
        tb::set_string(hnamestr, "MIP_SimhitE_total_layer_", i+1);
        mip_sim_total_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 150., 0, 150., 0., 800.);

        tb::set_string(hnamestr, "ADC_120mum_layer_", i+1);
        ADC_120mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 80, 0, 800.);
        tb::set_string(hnamestr, "MIP_120mum_layer_", i+1);
        MIP_120mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 150, 0, 150.);
        tb::set_string(hnamestr, "SIM_120mum_layer_", i+1);
        SIM_120mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 80, 0, 800.);
        tb::set_string(hnamestr, "ADC_SimhitE_120mum_layer_", i+1);
        adc_sim_120mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 800., 0, 800., 0., 800.);
        tb::set_string(hnamestr, "ADC_MIP_120mum_layer_", i+1);
        adc_mip_120mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 800., 0, 800., 0., 150.);
        tb::set_string(hnamestr, "MIP_SimhitE_120mum_layer_", i+1);
        mip_sim_120mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 150., 0, 150., 0., 800.);

        tb::set_string(hnamestr, "ADC_200mum_layer_", i+1);
        ADC_200mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 80, 0, 800.);
        tb::set_string(hnamestr, "MIP_200mum_layer_", i+1);
        MIP_200mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 150, 0, 150.);
        tb::set_string(hnamestr, "SIM_200mum_layer_", i+1);
        SIM_200mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 80, 0, 800.);
        tb::set_string(hnamestr, "ADC_SimhitE_200mum_layer_", i+1);
        adc_sim_200mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 800., 0, 800., 0., 800.);
        tb::set_string(hnamestr, "ADC_MIP_200mum_layer_", i+1);
        adc_mip_200mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 800., 0, 800., 0., 150.);
        tb::set_string(hnamestr, "MIP_SimhitE_200mum_layer_", i+1);
        mip_sim_200mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 150., 0, 150., 0., 800.);

        tb::set_string(hnamestr, "ADC_300mum_layer_", i+1);
        ADC_300mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 80, 0, 800.);
        tb::set_string(hnamestr, "MIP_300mum_layer_", i+1);
        MIP_300mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 150, 0, 150.);
        tb::set_string(hnamestr, "SIM_300mum_layer_", i+1);
        SIM_300mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 80, 0, 800.);
        tb::set_string(hnamestr, "ADC_SimhitE_300mum_layer_", i+1);
        adc_sim_300mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 800., 0, 800., 0., 800.);
        tb::set_string(hnamestr, "ADC_MIP_300mum_layer_", i+1);
        adc_mip_300mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 800., 0, 800., 0., 150.);
        tb::set_string(hnamestr, "MIP_SimhitE_300mum_layer_", i+1);
        mip_sim_300mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 150., 0, 150., 0., 800.);

        tb::set_string(hnamestr, "nt_total_layer_", i+1);
        nt_total_[i] = fs->make<TNtuple>(hnamestr.str().c_str(),hnamestr.str().c_str(), "adc:mip:simhitE");
        tb::set_string(hnamestr, "nt_120mum_layer_", i+1);
        nt_120mum_[i] = fs->make<TNtuple>(hnamestr.str().c_str(),hnamestr.str().c_str(), "adc:mip:simhitE");
        tb::set_string(hnamestr, "nt_200mum_layer_", i+1);
        nt_200mum_[i] = fs->make<TNtuple>(hnamestr.str().c_str(),hnamestr.str().c_str(), "adc:mip:simhitE");
        tb::set_string(hnamestr, "nt_300mum_layer_", i+1);
        nt_300mum_[i] = fs->make<TNtuple>(hnamestr.str().c_str(),hnamestr.str().c_str(), "adc:mip:simhitE");
    }

#ifdef this_is_an_eventsetup_example
    setupdatatoken_ = esConsumes<setupdata, setuprecord>();
#endif
} //}}}
DigiSim::~DigiSim(){}

// ------------ method called for each event  ------------
void DigiSim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // init {{{
    using namespace edm;

    int counter = 0;
    std::vector<double> total_energy_adc_total ;
    std::vector<double> total_energy_mip_total ;
    std::vector<double> total_energy_mip_coarse ;
    std::vector<double> total_energy_mip_fine ;
    std::vector<double> total_energy_sim_total ;
    std::vector<int>    num_digis_total        ;
    std::vector<int>    num_simhits_total      ;
    std::vector<int>    num_simhits_coarse      ;
    std::vector<int>    num_simhits_fine      ;

    std::vector<double> total_energy_adc_120mum ;
    std::vector<double> total_energy_mip_120mum ;
    std::vector<double> total_energy_sim_120mum ;
    std::vector<int>    num_digis_120mum        ;
    std::vector<int>    num_simhits_120mum      ;

    std::vector<double> total_energy_adc_200mum ;
    std::vector<double> total_energy_mip_200mum ;
    std::vector<double> total_energy_sim_200mum ;
    std::vector<int>    num_digis_200mum        ;
    std::vector<int>    num_simhits_200mum      ;

    std::vector<double> total_energy_adc_300mum ;
    std::vector<double> total_energy_mip_300mum ;
    std::vector<double> total_energy_sim_300mum ;
    std::vector<int>    num_digis_300mum        ;
    std::vector<int>    num_simhits_300mum      ;

    for(int idx=0; idx<26; ++idx) {
        total_energy_adc_total .push_back(0);
        total_energy_mip_total .push_back(0);
        total_energy_mip_coarse.push_back(0);
        total_energy_mip_fine  .push_back(0);
        total_energy_sim_total .push_back(0);
        num_digis_total        .push_back(0);
        num_simhits_total      .push_back(0);
        num_simhits_coarse     .push_back(0);
        num_simhits_fine       .push_back(0);

        total_energy_adc_120mum .push_back(0);
        total_energy_mip_120mum .push_back(0);
        total_energy_sim_120mum .push_back(0);
        num_digis_120mum        .push_back(0);
        num_simhits_120mum      .push_back(0);

        total_energy_adc_200mum .push_back(0);
        total_energy_mip_200mum .push_back(0);
        total_energy_sim_200mum .push_back(0);
        num_digis_200mum        .push_back(0);
        num_simhits_200mum      .push_back(0);

        total_energy_adc_300mum .push_back(0);
        total_energy_mip_300mum .push_back(0);
        total_energy_sim_300mum .push_back(0);
        num_digis_300mum        .push_back(0);
        num_simhits_300mum      .push_back(0);
    }
    
    //double Z_[47] = {322.155,323.149,325.212,326.206,328.269,329.263,331.326,332.32,334.383,335.377,337.44,338.434,340.497,341.491,343.554,344.548,346.611,347.605,349.993,350.987,353.375,354.369,356.757,357.751,360.139,361.133,367.976,374.281,380.586,386.891,393.196,399.501,405.806,412.111,418.416,424.721,431.026,439.251,447.476,455.701,463.926,472.151,480.376,488.601,496.826,505.051,513.276};

    //int geomType(0);
    //const HGCalGeometry* geom0 = &iSetup.getData(tok_hgcalg_);
    //std::cout<<geom0->topology().waferHexagon8()<<std::endl;
    
    /*********************************************************************************
     * Tool to convert Digis to Uncalibrated RecHits (amplitude in unit of MIPs)
     * $CMSSW_RELEASE_BASE/src/DataFormats/HGCRecHit/interface/HGCUncalibratedRecHit.h
     * $CMSSW_RELEASE_BASE/src/DataFormats/HGCRecHit/interface/HGCRecHitCollections.h
     *********************************************************************************/
    if (uncalibMaker_ee_.isSiFESim()) uncalibMaker_ee_.setGeometry(&iSetup.getData(ee_geometry_token_));

    const CaloGeometry &geomCalo = iSetup.getData(caloGeomToken_);
    rhtools_.setGeometry(geomCalo);
    
    //----------------------------------------------------------------------------------------------------}}}
    // SimHit Handle {{{
    //----------------------------------------------------------------------------------------------------
    std::map<uint32_t, std::pair<hitsinfo, energysum> > map_Simhits;
    map_Simhits.clear();
    unsigned int nofSiHits = 0;
    Handle<PCaloHitContainer> simhit;
    iEvent.getByToken(tSimCaloHitContainer, simhit);
    for(PCaloHitContainer::const_iterator itHit= simhit->begin(); itHit!= simhit->end(); ++itHit) {
        DetId detId = static_cast<DetId>(itHit->id());

        bool store_hit_position_info = true;
        if(store_hit_position_info) {
            // Warn: Might need to set proper sub detector configuration.
            // Currently the code is for CEE analayis
            GlobalPoint gp = rhtools_.getPosition(detId);
            float x = gp.x();
            float y = gp.y();
            float z = gp.z();
            float Rxy = sqrt(pow(x, 2) + pow(y, 2));
            bool isSilicon = rhtools_.isSilicon(detId);
            bool isScintillator = rhtools_.isScintillator(detId);
            int wafer_type = HGCSiliconDetId(detId).type();

            bool is_Silicon_w120 = false;
            bool is_Silicon_w200 = false;
            bool is_Silicon_w300 = false;
            bool is_Scintillator = isScintillator;

            if(isSilicon) {
                if(wafer_type==0) is_Silicon_w120 = true;
                if(wafer_type==1) is_Silicon_w200 = true;
                if(wafer_type==2) is_Silicon_w300 = true;
            }

            nt_hit_position->Fill(Rxy, z, is_Silicon_w120, is_Silicon_w200, is_Silicon_w300, is_Scintillator);

            bool debug = false;
            if(debug) {
                tb::print_debug_info("x"               , x                      );
                tb::print_debug_info("y"               , y                      );
                tb::print_debug_info("z"               , z                      );
                tb::print_debug_info("Rxy"             , Rxy                    );
                tb::print_debug_info("isSilicon"       , isSilicon              );
                tb::print_debug_info("isScintillator"  , isScintillator         );
                tb::print_debug_info("wafer_type"      , wafer_type             );
                tb::print_debug_info("is_Silicon_w120" , is_Silicon_w120        );
                tb::print_debug_info("is_Silicon_w200" , is_Silicon_w200        );
                tb::print_debug_info("is_Silicon_w300" , is_Silicon_w300        );
                tb::print_debug_info("is_Scintillator" , is_Scintillator , true );
            }
        }

        if(rhtools_.isSilicon(detId)){
            HGCSiliconDetId id(itHit->id());
            double energy = itHit->energy()*1.e6; // in kev

            if(nameDetector_ == "HGCalEESensitive"){
                hELossEE->Fill(energy);
                if(id.type()==HGCSiliconDetId::HGCalFine)        hELossEEF ->Fill(energy);
                if(id.type()==HGCSiliconDetId::HGCalCoarseThin)  hELossEECN->Fill(energy);
                if(id.type()==HGCSiliconDetId::HGCalCoarseThick) hELossEECK->Fill(energy);
            }

            if(nameDetector_ == "HGCalHESiliconSensitive"){
                hELossHEF->Fill(energy);
                if(id.type()==HGCSiliconDetId::HGCalFine)        hELossHEFF ->Fill(energy);
                if(id.type()==HGCSiliconDetId::HGCalCoarseThin)  hELossHEFCN->Fill(energy);
                if(id.type()==HGCSiliconDetId::HGCalCoarseThick) hELossHEFCK->Fill(energy);
            }  

            //GlobalPoint gp = rhtools_.getPosition(detId);
            //std::cout << "check: " << gp.x() << ", " << gp.y() << std::endl;

            //std::cout<<" hit energy = "<<itHit->energy()<<std::endl;
            uint32_t id_ = itHit->id();
            energysum esum;
            hitsinfo hinfo;

            if (map_Simhits.count(id_) != 0) {
                hinfo = map_Simhits[id_].first;
                esum = map_Simhits[id_].second;
            } else {
                hinfo.hitid = nofSiHits;
                hinfo.u_cor = rhtools_.getCell(detId).first ;
                hinfo.v_cor = rhtools_.getCell(detId).second ;
                hinfo.type  = id.type();
                hinfo.layer = rhtools_.getLayerWithOffset(detId);
            }	
            esum.etotal += energy;
            esum.eTime[0] = energy;
            map_Simhits[id_] = std::pair<hitsinfo, energysum>(hinfo, esum);

            bool debug = false;
            if(debug) {
                tb::print_debug_info("hinfo.hitid"   , hinfo.hitid         );
                tb::print_debug_info("hinfo.layer"   , hinfo.layer         );
                tb::print_debug_info("hinfo.type"    , hinfo.type          );
                tb::print_debug_info("hinfo.u_cor"   , hinfo.u_cor         );
                tb::print_debug_info("hinfo.v_cor"   , hinfo.v_cor         );
                tb::print_debug_info("esum.etotal"   , esum.etotal         );
                tb::print_debug_info("esum.eTime[0]" , esum.eTime[0], true );
            }
        } // end of isSilicon bool
    } // end of simhit loop
    //std::cout<<">>> simhit map size = "<< map_Simhits.size()<<std::endl;

    //----------------------------------------------------------------------------------------------------}}}
    // Digi Handle {{{
    //----------------------------------------------------------------------------------------------------
    Handle<HGCalDigiCollection> digicollection;
    iEvent.getByToken(digiSource_, digicollection);
    
    auto UncalibRechits = std::make_unique<HGCUncalibratedRecHitCollection>();
    UncalibRechits->reserve(digicollection->size());

    std::map<uint32_t, myDigis > my_map_digihits;
    std::map<uint32_t, std::pair<digisinfo,adcinfo > > map_digihits;
    if (digicollection.isValid()) {
        for (const auto& it : *(digicollection.product())) {
            DetId detId = it.id();
            if(rhtools_.isSilicon(detId)){
                uint32_t id_digi = uint32_t(it.id());

                const HGCSample& hgcSample = it.sample(SampleIndx_);
                uint16_t adc_ = hgcSample.data();
                //uint16_t gain = hgcSample.toa();

                GlobalPoint gp = rhtools_.getPosition(detId);
                float eta = rhtools_.getEta(detId, 0);
                float phi = rhtools_.getPhi(detId);
                
                digisinfo dinfo;
                adcinfo ainfo;
                if (map_digihits.count(id_digi) == 0) {
                    dinfo.u_cor = HGCSiliconDetId(detId).cellU();
                    dinfo.v_cor = HGCSiliconDetId(detId).cellV();
                    dinfo.type  = HGCSiliconDetId(detId).type();
                    dinfo.layer = HGCSiliconDetId(detId).layer();

                    dinfo.x_pos = gp.x();
                    dinfo.y_pos = gp.y();
                    dinfo.z_pos = gp.z();
                    dinfo.eta = eta;
                    dinfo.phi = phi;

                    ainfo.adc   = adc_;
                } else {
                    dinfo = map_digihits[id_digi].first;
                    ainfo = map_digihits[id_digi].second;
                }

                HGCUncalibratedRecHit rechit = uncalibMaker_ee_.makeRecHit(it);
                UncalibRechits->push_back( rechit );
                map_digihits[id_digi] = std::pair<digisinfo, adcinfo>(dinfo, ainfo);
                my_map_digihits[id_digi].dinfo = dinfo;
                my_map_digihits[id_digi].ainfo = ainfo;
                my_map_digihits[id_digi].amplitude = rechit.amplitude();

                int idx = dinfo.layer-1;
                if(dinfo.layer <= 26) num_digis_total[idx] += 1;
                if(dinfo.layer <= 26 && dinfo.type==0) num_digis_120mum[idx] += 1;
                if(dinfo.layer <= 26 && dinfo.type==1) num_digis_200mum[idx] += 1;
                if(dinfo.layer <= 26 && dinfo.type==2) num_digis_300mum[idx] += 1;

                bool debug = false;
                if(debug) {
                    tb::print_debug_info("Id_digi"    , id_digi                                  );
                    tb::print_debug_info("x_pos"      , my_map_digihits[id_digi].dinfo.x_pos     );
                    tb::print_debug_info("y_pos"      , my_map_digihits[id_digi].dinfo.y_pos     );
                    tb::print_debug_info("z_pos"      , my_map_digihits[id_digi].dinfo.z_pos     );
                    tb::print_debug_info("eta"        , my_map_digihits[id_digi].dinfo.eta       );
                    tb::print_debug_info("phi"        , my_map_digihits[id_digi].dinfo.phi       );
                    tb::print_debug_info("layer"      , my_map_digihits[id_digi].dinfo.layer     );
                    tb::print_debug_info("wafer type" , my_map_digihits[id_digi].dinfo.type      );
                    tb::print_debug_info("adc_"       , my_map_digihits[id_digi].ainfo.adc       );
                    tb::print_debug_info("amplitude"  , my_map_digihits[id_digi].amplitude, true );
                }

            } // end of isSilicon bool
        } // end of digicollection loop
    } // end of digicollection valid

    // loop over HGCEE digis
    if(false) {
        printf(">>> check UncalibRechits->size(): %ld\n", UncalibRechits->size());
        counter = 0;
        for (auto itdg = UncalibRechits->begin(); itdg != UncalibRechits->end(); ++itdg) {
          if(counter>9) continue;
          std::cout<<"id_digi = "<<(*itdg).id().rawId();
          std::cout<<", amplitude = "<< (*itdg).amplitude()<<std::endl;
          counter++;
        }
    }

    //----------------------------------------------------------------------------------------------------}}}
    // Matching digiHits to simHits {{{
    //----------------------------------------------------------------------------------------------------
    std::map<uint32_t, myDigis>::iterator itr_mydigi;
    std::map<uint32_t, std::pair<digisinfo, adcinfo>>::iterator itr_digi;
    std::map<uint32_t, std::pair<hitsinfo, energysum> >::iterator itr_sim;
    Double_t max_energy=0; 
    for (itr_sim = map_Simhits.begin(); itr_sim != map_Simhits.end(); ++itr_sim) {
        energysum esum = (*itr_sim).second.second;
        if(max_energy<esum.eTime[0]) max_energy=esum.eTime[0];
    }

    for (itr_sim = map_Simhits.begin(); itr_sim != map_Simhits.end(); ++itr_sim) {
        uint32_t id_simhit = (*itr_sim).first;
        energysum esum = (*itr_sim).second.second;
        if(esum.etotal>0) // && esum.eTime[0]==max_energy)
        {
            for (itr_mydigi = my_map_digihits.begin(); itr_mydigi != my_map_digihits.end(); ++itr_mydigi) {
                uint32_t  id_digihit = (*itr_mydigi).first;
                digisinfo dinfo      = (*itr_mydigi).second.dinfo;
                uint32_t  adc        = (*itr_mydigi).second.ainfo.adc;
                double    eta        = dinfo.eta;
                double    phi        = dinfo.phi;
                double    amplitude  = (*itr_mydigi).second.amplitude;
                double    energy     = esum.eTime[0];
                int       idx        = dinfo.layer-1;
                int       cellType   = (*itr_sim).second.first.type;
                bool      is_coarse  = cellType==HGCSiliconDetId::HGCalCoarseThin || cellType==HGCSiliconDetId::HGCalCoarseThick;
                bool      is_fine    = cellType==HGCSiliconDetId::HGCalFine;

                bool selection_on_mips = amplitude > 0.5;
                if(!selection_on_mips) continue;

                if(id_simhit==id_digihit){
                    bool debug = false;
                    if(debug) {
                        tb::print_debug_info("Id_digi"   , id_digihit   );
                        tb::print_debug_info("layer"     , dinfo.layer  );
                        tb::print_debug_info("wafer type", dinfo.type   );
                        tb::print_debug_info("cell type" , cellType   );
                        tb::print_debug_info("adc_"      , adc          );
                        tb::print_debug_info("amplitude" , amplitude    );
                        tb::print_debug_info("energy"    , energy, true );
                    }
                    hEta->Fill(eta);
                    hPhi->Fill(phi);
                    // individual hit info
                    if(dinfo.layer <= 26) {
                        ADC_total_     [idx] -> Fill(adc);
                        MIP_total_     [idx] -> Fill(amplitude);
                        SIM_total_     [idx] -> Fill(energy);
                        adc_sim_total_ [idx] -> Fill(adc,energy);
                        adc_mip_total_ [idx] -> Fill(adc,amplitude);
                        mip_sim_total_ [idx] -> Fill(amplitude,energy);
                        nt_total_      [idx] -> Fill(adc,amplitude,energy);

                        total_energy_adc_total [idx] += adc;
                        total_energy_mip_total [idx] += amplitude;
                        if(is_coarse) total_energy_mip_coarse[idx] += amplitude;
                        if(is_fine)   total_energy_mip_fine  [idx] += amplitude;
                        total_energy_sim_total [idx] += energy;
                        num_simhits_total      [idx] += 1;
                        if(is_coarse) num_simhits_coarse     [idx] += 1;
                        if(is_fine)   num_simhits_fine       [idx] += 1;
                    }
                    if(dinfo.layer <= 26 && dinfo.type==0) {
                        ADC_120mum_     [idx] -> Fill(adc);
                        MIP_120mum_     [idx] -> Fill(amplitude);
                        SIM_120mum_     [idx] -> Fill(energy);
                        adc_sim_120mum_ [idx] -> Fill(adc,energy);
                        adc_mip_120mum_ [idx] -> Fill(adc,amplitude);
                        mip_sim_120mum_ [idx] -> Fill(amplitude,energy);
                        nt_120mum_      [idx] -> Fill(adc,amplitude,energy);

                        total_energy_adc_120mum [idx] += adc;
                        total_energy_mip_120mum [idx] += amplitude;
                        total_energy_sim_120mum [idx] += energy;
                        num_simhits_120mum      [idx] += 1;
                    }
                    if(dinfo.layer <= 26 && dinfo.type==1) {
                        ADC_200mum_     [idx] -> Fill(adc);
                        MIP_200mum_     [idx] -> Fill(amplitude);
                        SIM_200mum_     [idx] -> Fill(energy);
                        adc_sim_200mum_ [idx] -> Fill(adc,energy);
                        adc_mip_200mum_ [idx] -> Fill(adc,amplitude);
                        mip_sim_200mum_ [idx] -> Fill(amplitude,energy);
                        nt_200mum_      [idx] -> Fill(adc,amplitude,energy);

                        total_energy_adc_200mum [idx] += adc;
                        total_energy_mip_200mum [idx] += amplitude;
                        total_energy_sim_200mum [idx] += energy;
                        num_simhits_200mum      [idx] += 1;
                    }
                    if(dinfo.layer <= 26 && dinfo.type==2) {
                        ADC_300mum_     [idx] -> Fill(adc);
                        MIP_300mum_     [idx] -> Fill(amplitude);
                        SIM_300mum_     [idx] -> Fill(energy);
                        adc_sim_300mum_ [idx] -> Fill(adc,energy);
                        adc_mip_300mum_ [idx] -> Fill(adc,amplitude);
                        mip_sim_300mum_ [idx] -> Fill(amplitude,energy);
                        nt_300mum_      [idx] -> Fill(adc,amplitude,energy);

                        total_energy_adc_300mum [idx] += adc;
                        total_energy_mip_300mum [idx] += amplitude;
                        total_energy_sim_300mum [idx] += energy;
                        num_simhits_300mum      [idx] += 1;
                    }
                }
            } // end of digihits for loop
        }
    } // end of simhits for loop

    // Fill information of an event
    for(int idx=0; idx<26; ++idx) {
        //tb::print_debug_info("num_digis_total"  , num_digis_total  [idx]        );
        //tb::print_debug_info("num_digis_120mum" , num_digis_120mum [idx]        );
        //tb::print_debug_info("num_digis_200mum" , num_digis_200mum [idx]        );
        //tb::print_debug_info("num_digis_300mum" , num_digis_300mum [idx] , true );

        total_ADC_total_             [idx] -> Fill( total_energy_adc_total [idx] );
        total_MIP_total_             [idx] -> Fill( total_energy_mip_total [idx] );
        total_MIP_coarse_            [idx] -> Fill( total_energy_mip_coarse [idx] );
        total_MIP_fine_              [idx] -> Fill( total_energy_mip_fine [idx] );
        total_SIM_total_             [idx] -> Fill( total_energy_sim_total [idx] );
        multiplicity_digis_total_    [idx] -> Fill( num_digis_total        [idx] );
        multiplicity_simhits_total_  [idx] -> Fill( num_simhits_total      [idx] );
        multiplicity_simhits_coarse_ [idx] -> Fill( num_simhits_coarse     [idx] );
        multiplicity_simhits_fine_   [idx] -> Fill( num_simhits_fine       [idx] );

        total_ADC_120mum_            [idx] -> Fill( total_energy_adc_120mum [idx] );
        total_MIP_120mum_            [idx] -> Fill( total_energy_mip_120mum [idx] );
        total_SIM_120mum_            [idx] -> Fill( total_energy_sim_120mum [idx] );
        multiplicity_digis_120mum_   [idx] -> Fill( num_digis_120mum        [idx] );
        multiplicity_simhits_120mum_ [idx] -> Fill( num_simhits_120mum      [idx] );

        total_ADC_200mum_            [idx] -> Fill( total_energy_adc_200mum [idx] );
        total_MIP_200mum_            [idx] -> Fill( total_energy_mip_200mum [idx] );
        total_SIM_200mum_            [idx] -> Fill( total_energy_sim_200mum [idx] );
        multiplicity_digis_200mum_   [idx] -> Fill( num_digis_200mum        [idx] );
        multiplicity_simhits_200mum_ [idx] -> Fill( num_simhits_200mum      [idx] );

        total_ADC_300mum_            [idx] -> Fill( total_energy_adc_300mum [idx] );
        total_MIP_300mum_            [idx] -> Fill( total_energy_mip_300mum [idx] );
        total_SIM_300mum_            [idx] -> Fill( total_energy_sim_300mum [idx] );
        multiplicity_digis_300mum_   [idx] -> Fill( num_digis_300mum        [idx] );
        multiplicity_simhits_300mum_ [idx] -> Fill( num_simhits_300mum      [idx] );
    }
    //}}}
} // end of analyze

// others {{{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
// if the SetupData is always needed
auto setup = iSetup.getData(setupToken_);
// if need the ESHandle to check if the SetupData was there or not
auto pSetup = iSetup.getHandle(setupToken_);
#endif


// ------------ method called once each job just before starting event loop  ------------
void DigiSim::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void DigiSim::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DigiSim::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);

    //Specify that only 'tracks' is allowed
    //To use, remove the default given above and uncomment below
    //ParameterSetDescription desc;
    //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
    //descriptions.addDefault(desc);
}
//}}}
//define this as a plug-in
DEFINE_FWK_MODULE(DigiSim);
