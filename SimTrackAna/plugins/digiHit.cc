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
                hitid = ndigis = 0;
            }
            int u_cor, v_cor, type, layer;
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

        hgcal::RecHitTools rhtools_;
        edm::EDGetToken digiSource_;
        //edm::ConsumesCollector iC;
        edm::ESGetToken<HGCalGeometry, IdealGeometryRecord> ee_geometry_token_;
        HGCalUncalibRecHitRecWeightsAlgo<HGCalDataFrame> uncalibMaker_ee_;

        TH1D *hELossEE;TH1D *hELossEEF;TH1D *hELossEECN;TH1D *hELossEECK;
        TH1D *hELossHEF;TH1D *hELossHEFF;TH1D *hELossHEFCN;TH1D *hELossHEFCK;
        std::vector<TH1D*> vechist;   
        TH1D *ADC_120mum_[26];
        TH1D *ADC_200mum_[26];
        TH1D *ADC_300mum_[26];
        TH1D *MIP_120mum_[26];
        TH1D *MIP_200mum_[26];
        TH1D *MIP_300mum_[26];
        TH1D *SIM_120mum_[26];
        TH1D *SIM_200mum_[26];
        TH1D *SIM_300mum_[26];
        TProfile *adc_sim_120mum_[26];
        TProfile *adc_sim_200mum_[26];
        TProfile *adc_sim_300mum_[26];
        TProfile *adc_mip_120mum_[26];
        TProfile *adc_mip_200mum_[26];
        TProfile *adc_mip_300mum_[26];
        TProfile *mip_sim_120mum_[26];
        TProfile *mip_sim_200mum_[26];
        TProfile *mip_sim_300mum_[26];

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
    std::ostringstream hnamestr (std::ostringstream::ate);
    for(int i=0;i<26;i++){
        tb::set_string(hnamestr, "ADC_120mum_layer_", i+1);
        ADC_120mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 1024, 0, 1024.);
        tb::set_string(hnamestr, "MIP_120mum_layer_", i+1);
        MIP_120mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 500, 0, 500.);
        tb::set_string(hnamestr, "SIM_120mum_layer_", i+1);
        SIM_120mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 1000, 0, 1000.);

        tb::set_string(hnamestr, "ADC_SimhitE_120mum_layer_", i+1);
        adc_sim_120mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),1024., 0, 1024.,0., 1000.);
        tb::set_string(hnamestr, "ADC_MIP_120mum_layer_", i+1);
        adc_mip_120mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),1024., 0, 1024.,0., 500.);
        tb::set_string(hnamestr, "MIP_SimhitE_120mum_layer_", i+1);
        mip_sim_120mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),500., 0, 500.,0., 1000.);

        tb::set_string(hnamestr, "ADC_200mum_layer_", i+1);
        ADC_200mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 1024, 0, 1024.);
        tb::set_string(hnamestr, "MIP_200mum_layer_", i+1);
        MIP_200mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 500, 0, 500.);
        tb::set_string(hnamestr, "SIM_200mum_layer_", i+1);
        SIM_200mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 1000, 0, 1000.);

        tb::set_string(hnamestr, "ADC_SimhitE_200mum_layer_", i+1);
        adc_sim_200mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),1024., 0, 1024.,0., 1000.);
        tb::set_string(hnamestr, "ADC_MIP_200mum_layer_", i+1);
        adc_mip_200mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),1024., 0, 1024.,0., 500.);
        tb::set_string(hnamestr, "MIP_SimhitE_200mum_layer_", i+1);
        mip_sim_200mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),500., 0, 500.,0., 1000.);

        tb::set_string(hnamestr, "ADC_300mum_layer_", i+1);
        ADC_300mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 1024, 0, 1024.);
        tb::set_string(hnamestr, "MIP_300mum_layer_", i+1);
        MIP_300mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 500, 0, 500.);
        tb::set_string(hnamestr, "SIM_300mum_layer_", i+1);
        SIM_300mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 1000, 0, 1000.);

        tb::set_string(hnamestr, "ADC_SimhitE_300mum_layer_", i+1);
        adc_sim_300mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),1024., 0, 1024.,0., 1000.);
        tb::set_string(hnamestr, "ADC_MIP_300mum_layer_", i+1);
        adc_mip_300mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),1024., 0, 1024.,0., 500.);
        tb::set_string(hnamestr, "MIP_SimhitE_300mum_layer_", i+1);
        mip_sim_300mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),500., 0, 500.,0., 1000.);

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
    int counter = 0;
    using namespace edm;
    //int geomType(0);
    //const HGCalGeometry* geom0 = &iSetup.getData(tok_hgcalg_);
    //std::cout<<geom0->topology().waferHexagon8()<<std::endl;
    
    /*********************************************************************************
     * Tool to convert Digis to Uncalibrated RecHits (amplitude in unit of MIPs)
     * $CMSSW_RELEASE_BASE/src/DataFormats/HGCRecHit/interface/HGCUncalibratedRecHit.h
     * $CMSSW_RELEASE_BASE/src/DataFormats/HGCRecHit/interface/HGCRecHitCollections.h
     *********************************************************************************/
    if (uncalibMaker_ee_.isSiFESim()) uncalibMaker_ee_.setGeometry(&iSetup.getData(ee_geometry_token_));
    
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
        std::cout<<"valid"<<std::endl;
        for (const auto& it : *(digicollection.product())) {
            DetId detId = it.id();
            if(rhtools_.isSilicon(detId)){
                uint32_t id_digi = uint32_t(it.id());

                const HGCSample& hgcSample = it.sample(SampleIndx_);
                uint16_t adc_ = hgcSample.data();
                //uint16_t gain = hgcSample.toa();

                digisinfo dinfo;
                adcinfo ainfo;
                if (map_digihits.count(id_digi) != 0) {
                    dinfo = map_digihits[id_digi].first;
                    ainfo = map_digihits[id_digi].second;
                } else {
                    dinfo.u_cor = HGCSiliconDetId(detId).cellU() ;
                    dinfo.v_cor = HGCSiliconDetId(detId).cellV() ;
                    dinfo.type  = HGCSiliconDetId(detId).type();
                    dinfo.layer = HGCSiliconDetId(detId).layer();
                    ainfo.adc   = adc_;
                }

                HGCUncalibratedRecHit rechit = uncalibMaker_ee_.makeRecHit(it);
                UncalibRechits->push_back( rechit );
                map_digihits[id_digi] = std::pair<digisinfo, adcinfo>(dinfo, ainfo);
                my_map_digihits[id_digi].dinfo = dinfo;
                my_map_digihits[id_digi].ainfo = ainfo;
                my_map_digihits[id_digi].amplitude = rechit.amplitude();

                bool debug = false;
                if(debug) {
                    tb::print_debug_info("Id_digi"   , id_digi                                  );
                    tb::print_debug_info("layer"     , my_map_digihits[id_digi].dinfo.layer     );
                    tb::print_debug_info("wafer type", my_map_digihits[id_digi].dinfo.type      );
                    tb::print_debug_info("adc_"      , my_map_digihits[id_digi].ainfo.adc       );
                    tb::print_debug_info("amplitude" , my_map_digihits[id_digi].amplitude, true );
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

    printf(">>> start matching!\n");
    for (itr_sim = map_Simhits.begin(); itr_sim != map_Simhits.end(); ++itr_sim) {
        uint32_t id_simhit = (*itr_sim).first;
        energysum esum = (*itr_sim).second.second;
        if(esum.etotal>0) // && esum.eTime[0]==max_energy)
        {
            for (itr_mydigi = my_map_digihits.begin(); itr_mydigi != my_map_digihits.end(); ++itr_mydigi) {
                uint32_t  id_digihit = (*itr_mydigi).first;
                digisinfo dinfo      = (*itr_mydigi).second.dinfo;
                adcinfo   ainfo      = (*itr_mydigi).second.ainfo;
                double    amplitude  = (*itr_mydigi).second.amplitude;
                double    energy     = esum.eTime[0];
                int       idx        = dinfo.layer-1;

                if(id_simhit==id_digihit){
                    bool debug = false;
                    if(debug) {
                        tb::print_debug_info("Id_digi"   , id_digihit );
                        tb::print_debug_info("layer"     , dinfo.layer         );
                        tb::print_debug_info("wafer type", dinfo.type          );
                        tb::print_debug_info("adc_"      , ainfo.adc           );
                        tb::print_debug_info("amplitude" , amplitude           );
                        tb::print_debug_info("energy"    , energy, true        );
                    }
                    if(dinfo.layer <= 26 && dinfo.type==0) {
                        ADC_120mum_[idx]->Fill(ainfo.adc);
                        MIP_120mum_[idx]->Fill(amplitude);
                        SIM_120mum_[idx]->Fill(energy);
                        adc_sim_120mum_[idx]->Fill(ainfo.adc,energy);
                        adc_mip_120mum_[idx]->Fill(ainfo.adc,amplitude);
                        mip_sim_120mum_[idx]->Fill(amplitude,energy);
                    }
                    if(dinfo.layer <= 26 && dinfo.type==1) {
                        ADC_200mum_[idx]->Fill(ainfo.adc);
                        MIP_200mum_[idx]->Fill(amplitude);
                        SIM_200mum_[idx]->Fill(energy);
                        adc_sim_200mum_[idx]->Fill(ainfo.adc,energy);
                        adc_mip_200mum_[idx]->Fill(ainfo.adc,amplitude);
                        mip_sim_200mum_[idx]->Fill(amplitude,energy);
                    }
                    if(dinfo.layer <= 26 && dinfo.type==2) {
                        ADC_300mum_[idx]->Fill(ainfo.adc);
                        MIP_300mum_[idx]->Fill(amplitude);
                        SIM_300mum_[idx]->Fill(energy);
                        adc_sim_300mum_[idx]->Fill(ainfo.adc,energy);
                        adc_mip_300mum_[idx]->Fill(ainfo.adc,amplitude);
                        mip_sim_300mum_[idx]->Fill(amplitude,energy);
                    }
                }
            } // end of digihits for loop
        }
    } // end of simhits for loop
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
