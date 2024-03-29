#include <memory>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <iostream>
// user include files {{{

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToModule.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToROC.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"

#include "CoralBase/Exception.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalUncalibRecHitWorkerFactory.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalUncalibRecHitRecWeightsAlgo.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TVector3.h"
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

class DigiSim : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        //Implemented following Validation/HGCalValidation/plugins/HGCalSimHitValidation.cc
        // structures {{{
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
                x_pos = y_pos = z_pos = 0.;
                hitid = nhits = 0;
                is_Silicon_w120 = is_Silicon_w200 = is_Silicon_w300 = is_Scintillator = false;
            }
            int u_cor, v_cor, type, layer;
            float x_pos, y_pos, z_pos;
            unsigned int hitid, nhits;
            bool is_Silicon_w120, is_Silicon_w200, is_Silicon_w300, is_Scintillator;
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

        struct myVector {
            // expected hits for 26 layers (energy, id, vector of x, y, z positions)
            std::vector<float> tr_vx;
            std::vector<float> tr_vy;
            std::vector<float> tr_vz;
            std::vector<float> tr_ve;
        };
        //}}}
        explicit DigiSim(const edm::ParameterSet&);
        ~DigiSim();
        double get_corrected_energy_from_dEdx_method(int layer, double amplitude, TString tag);
        double get_additional_correction(int layer);
        bool is_this_in_set1(int layer);
        bool is_this_in_set2(int layer);
        double convert_amplitude_to_total_energy_pedro(int type, double amplitude);
        void reset_tree_variables();
        void reset_per_event_counters();
        void reset_expected_hit_containers(myVector &mv);
        void calculate_efficiency();
        void fill_event_info();
        GlobalPoint projectHitPositionAt(float z,float eta,float phi);
        float get_distance_from_expected_hit(double x, double y, double z, double eta, double phi);
        float get_distance_from_expected_hit(double x, double y, double x0, double y0);
        int get_signal_region(float d);
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        //const std::string name;
        // details {{{
        const std::string nameDetector_; 
        const bool ifNose_;
        const int verbosity_, SampleIndx_;
        // ----------member data ---------------------------
        const edm::ESGetToken<HGCalGeometry, IdealGeometryRecord> tok_hgcalg_;
        int firstLayer_; 
        edm::EDGetTokenT<edm::PCaloHitContainer> tSimCaloHitContainer; 
        edm::EDGetTokenT<edm::HepMCProduct> mc_;
        edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;
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
        TTree   *tr_hits;
        TTree   *tr_positron;
        TTree   *tr_max_cell;
        TTree   *tr_energy_weighted;
        TTree   *tr_linear_trajectory;
        int tr_evtNo;
        int tr_layerNo;
        int tr_signal_region_max_cell;
        int tr_signal_region_energy_weighted;
        int tr_signal_region_linear_track;
        float tr_x;
        float tr_y;
        float tr_z;
        float tr_e;
        float tr_r;
        float tr_d_max_cell;
        float tr_d_energy_weighted;
        float tr_d_linear_track;
        float tr_eta;
        float tr_phi;
        bool tr_is_Silicon_w120;
        bool tr_is_Silicon_w200;
        bool tr_is_Silicon_w300;
        bool tr_is_Scintillator;

        // expected hits for 26 layers (energy, vector of x, y, z positions)
        myVector mv_max_cell;
        myVector mv_energy_weighted;
        myVector mv_linear_track;

        // hiistograms
        std::vector<TH1D*> vechist;   
        TNtuple *nt_total_[26];
        TH1D *ADC_total_[26];
        TH1D *MIP_total_[26];
        TH1D *SIM_total_[26];
        TProfile *adc_sim_total_[26];
        TProfile *adc_mip_total_[26];
        TProfile *mip_sim_total_[26];

        TH1D *efficiency_linear_track_[26];

        // total energy
        TH1D *total_ADC_total_[26];
        TH1D *total_MIP_total_[26];
        TH1D *total_MIP_coarse_[26];
        TH1D *total_MIP_fine_[26];
        TH1D *total_SIM_total_[26];

        TH1D *total_MIP_odd;
        TH1D *total_MIP_even;
        TH1D *total_MIP_set0;
        TH1D *total_MIP_set1;
        TH1D *total_MIP_set2;

        TH1D *total_SIM_odd;
        TH1D *total_SIM_even;
        TH1D *total_SIM_set0;
        TH1D *total_SIM_set1;
        TH1D *total_SIM_set2;

        //TH1D *total_ENE_odd;
        //TH1D *total_ENE_even;
        TH1D *total_ENE_set0;
        TH1D *total_ENE_set1;
        TH1D *total_ENE_set2;

        // multiplicity
        TH1D *multiplicity_digis_total_[26];
        TH1D *multiplicity_simhits_total_[26];
        TH1D *multiplicity_simhits_coarse_[26];
        TH1D *multiplicity_simhits_fine_[26];

        // dE/dx weights from https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/HGCalRecProducers/python/HGCalRecHit_cfi.py#L12-L60
        std::vector<double> weightsPerLayer_V16 = { 0., 5.55, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 13.54, 12.86, 13.54, 12.86, 13.54, 12.86, 13.54, 12.86,
                                                    58.63, 60.7, 60.7, 60.7, 60.7, 60.7, 60.7, 60.7, 60.7, 60.7, 60.7, 83.08, 83.08, 83.43, 83.61, 83.61, 83.61, 83.61, 83.61, 83.61, 83.61 };
        std::vector<double> calibration_weights_set0 = {0.00,  9.21, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 11.13, 13.20, 13.20, 13.20, 13.20, 13.20, 13.20, 13.20, 13.20, 35.75};
        std::vector<double> calibration_weights_set1 = {0.00,  13.91,  0.00, 22.26,  0.00, 22.26,  0.00, 22.26,  0.00, 22.26,  0.00, 22.26,  0.00, 22.26,  0.00, 22.26,  0.00, 24.33,  0.00, 26.40,  0.00, 26.40,  0.00, 26.40,  0.00, 48.95,  0.00};
        std::vector<double> calibration_weights_set2 = {0.00,  13.91,  0.00, 22.26,  0.00, 28.69,  0.00,  0.00, 28.69,  0.00, 22.26,  0.00, 22.26,  0.00, 15.83, 15.83,  0.00, 24.33,  0.00, 26.40,  0.00, 26.40,  0.00, 26.40,  0.00, 48.95,  0.00};

        std::vector<double> x_D86 = { 0., 0.564,1.567,2.547,3.549,4.528,5.531,6.509,7.512,8.49,9.493,10.472,11.474,12.453,13.455,14.434,15.437,16.415,17.418,18.975,19.978,21.536,22.538,24.096,25.099,26.656,27.659 };
        //double Z_[47] = {322.155,323.149,325.212,326.206,328.269,329.263,331.326,332.32,334.383,335.377,337.44,338.434,340.497,341.491,343.554,344.548,346.611,347.605,349.993,350.987,353.375,354.369,356.757,357.751,360.139,361.133,367.976,374.281,380.586,386.891,393.196,399.501,405.806,412.111,418.416,424.721,431.026,439.251,447.476,455.701,463.926,472.151,480.376,488.601,496.826,505.051,513.276};
        //}}}
        // per event counters {{{
        int counter = 0;
        double total_energy_mip_odd  = 0.;
        double total_energy_mip_even = 0.;
        double total_energy_mip_set0 = 0.;
        double total_energy_mip_set1 = 0.;
        double total_energy_mip_set2 = 0.;
        double total_energy_sim_odd  = 0.;
        double total_energy_sim_even = 0.;
        double total_energy_sim_set0 = 0.;
        double total_energy_sim_set1 = 0.;
        double total_energy_sim_set2 = 0.;
        double total_corrected_energy_odd  = 0.;
        double total_corrected_energy_even = 0.;
        double total_corrected_energy_set0 = 0.;
        double total_corrected_energy_set1 = 0.;
        double total_corrected_energy_set2 = 0.;

        std::vector<double> efficiency_numerators   ;
        std::vector<double> efficiency_denominators ;
        std::vector<double> efficiency_signal_region_linear_track ;

        std::vector<double> total_energy_adc_total  ;
        std::vector<double> total_energy_mip_total  ;
        std::vector<double> total_energy_mip_coarse ;
        std::vector<double> total_energy_mip_fine   ;
        std::vector<double> total_energy_sim_total  ;
        std::vector<int>    num_digis_total         ;
        std::vector<int>    num_simhits_total       ;
        std::vector<int>    num_simhits_coarse      ;
        std::vector<int>    num_simhits_fine        ;

        //}}}

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
        edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

DigiSim::DigiSim(const edm::ParameterSet& iconfig) : //{{{
        //auto temp = iConfig.getUntrackedParameter<edm::InputTag>("digihits");
        nameDetector_(iconfig.getParameter<std::string>("Detector")),
        ifNose_(iconfig.getUntrackedParameter<bool>("ifNose")),
        verbosity_(iconfig.getUntrackedParameter<int>("Verbosity", 0)),
        SampleIndx_(iconfig.getUntrackedParameter<int>("SampleIndx", 0)),
        tok_hgcalg_(esConsumes<HGCalGeometry, IdealGeometryRecord>(edm::ESInputTag{"", nameDetector_})), 
        firstLayer_(1), 
        tSimCaloHitContainer(consumes<edm::PCaloHitContainer>(iconfig.getUntrackedParameter<edm::InputTag>("simhits"))),
        mc_( consumes<edm::HepMCProduct>(edm::InputTag("generatorSmeared")) ),
        genParticles_( consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles")) ),
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

    // tree of hits
    tr_hits = fs->make<TTree>("tr_hits","");
    tr_hits -> Branch("evtNo"   , &tr_evtNo   );
    tr_hits -> Branch("layerNo" , &tr_layerNo );
    tr_hits -> Branch("x"       , &tr_x       );
    tr_hits -> Branch("y"       , &tr_y       );
    tr_hits -> Branch("z"       , &tr_z       );
    tr_hits -> Branch("e"       , &tr_e       );
    tr_hits -> Branch("r"       , &tr_r       );
    tr_hits -> Branch("eta"     , &tr_eta     );
    tr_hits -> Branch("phi"     , &tr_phi     );

    tr_hits -> Branch("d_max_cell"                    , &tr_d_max_cell                    );
    tr_hits -> Branch("d_energy_weighted"             , &tr_d_energy_weighted             );
    tr_hits -> Branch("d_linear_track"                , &tr_d_linear_track                );

    tr_hits -> Branch("signal_region_linear_track"    , &tr_signal_region_linear_track    );
    tr_hits -> Branch("signal_region_max_cell"        , &tr_signal_region_max_cell        );
    tr_hits -> Branch("signal_region_energy_weighted" , &tr_signal_region_energy_weighted );

    tr_hits -> Branch("is_Silicon_w120", &tr_is_Silicon_w120);
    tr_hits -> Branch("is_Silicon_w200", &tr_is_Silicon_w200);
    tr_hits -> Branch("is_Silicon_w300", &tr_is_Silicon_w300);
    tr_hits -> Branch("is_Scintillator", &tr_is_Scintillator);

    tr_evtNo = 0;

    // tree of truth positron
    tr_positron = fs->make<TTree>("tr_positron","");
    tr_positron -> Branch("evtNo"   , &tr_evtNo   );
    tr_positron -> Branch("e"       , &tr_e       );
    tr_positron -> Branch("eta"     , &tr_eta     );
    tr_positron -> Branch("phi"     , &tr_phi     );

    tr_max_cell = fs->make<TTree>("tr_max_cell","");
    tr_max_cell -> Branch("evtNo"   , &tr_evtNo  );
    tr_max_cell -> Branch("vx"      , &mv_max_cell.tr_vx     );
    tr_max_cell -> Branch("vy"      , &mv_max_cell.tr_vy     );
    tr_max_cell -> Branch("vz"      , &mv_max_cell.tr_vz     );
    tr_max_cell -> Branch("ve"      , &mv_max_cell.tr_ve     );

    tr_energy_weighted = fs->make<TTree>("tr_energy_weighted","");
    tr_energy_weighted -> Branch("evtNo"   , &tr_evtNo  );
    tr_energy_weighted -> Branch("vx"      , &mv_energy_weighted.tr_vx     );
    tr_energy_weighted -> Branch("vy"      , &mv_energy_weighted.tr_vy     );
    tr_energy_weighted -> Branch("vz"      , &mv_energy_weighted.tr_vz     );
    tr_energy_weighted -> Branch("ve"      , &mv_energy_weighted.tr_ve     );

    tr_linear_trajectory = fs->make<TTree>("tr_linear_trajectory","");
    tr_linear_trajectory -> Branch("evtNo"   , &tr_evtNo  );
    tr_linear_trajectory -> Branch("vx"      , &mv_linear_track.tr_vx     );
    tr_linear_trajectory -> Branch("vy"      , &mv_linear_track.tr_vy     );
    tr_linear_trajectory -> Branch("vz"      , &mv_linear_track.tr_vz     );
    tr_linear_trajectory -> Branch("ve"      , &mv_linear_track.tr_ve     );

    // histograms
    hELossEE    = fs->make<TH1D>("hELossEE"    , "hELossEE"    , 1000 , 0. , 1000.);
    hELossEEF   = fs->make<TH1D>("hELossEEF"   , "hELossEEF"   , 1000 , 0. , 1000.);
    hELossEECN  = fs->make<TH1D>("hELossEECN"  , "hELossEECN"  , 1000 , 0. , 1000.);
    hELossEECK  = fs->make<TH1D>("hELossEECK"  , "hELossEECK"  , 1000 , 0. , 1000.);
    hELossHEF   = fs->make<TH1D>("hELossHEF"   , "hELossHEF"   , 1000 , 0. , 1000.);
    hELossHEFF  = fs->make<TH1D>("hELossHEFF"  , "hELossHEFF"  , 1000 , 0. , 1000.);
    hELossHEFCN = fs->make<TH1D>("hELossHEFCN" , "hELossHEFCN" , 1000 , 0. , 1000.);
    hELossHEFCK = fs->make<TH1D>("hELossHEFCK" , "hELossHEFCK" , 1000 , 0. , 1000.);
    hEta = fs->make<TH1D>("hEta" , "hEta" , 100 , -5. , 5.); // [1., 3.]
    hPhi = fs->make<TH1D>("hPhi" , "hPhi" , 20 , -3. , 3.);

    //E_set0 = all CEE layers
    //E_set1 = E1+E3+E5+E7+E9+E11+E13+E15+E17+...E25
    //E_set2 = E1+E3+E5+E8+E10+E12+E14+E15+E17+...E25
    
    total_MIP_odd  = fs->make<TH1D>("total_MIP_odd"  , "total_MIP_odd"  , 600  , 0. ,  30000.);
    total_MIP_even = fs->make<TH1D>("total_MIP_even" , "total_MIP_even" , 600  , 0. ,  30000.);
    total_MIP_set0 = fs->make<TH1D>("total_MIP_set0" , "total_MIP_set0" , 1200 , 0. ,  60000.);
    total_MIP_set1 = fs->make<TH1D>("total_MIP_set1" , "total_MIP_set1" , 600  , 0. ,  30000.);
    total_MIP_set2 = fs->make<TH1D>("total_MIP_set2" , "total_MIP_set2" , 600  , 0. ,  30000.);

    total_SIM_odd  = fs->make<TH1D>("total_SIM_odd"  , "total_SIM_odd"  , 200 , 0. , 200.);
    total_SIM_even = fs->make<TH1D>("total_SIM_even" , "total_SIM_even" , 200 , 0. , 200.);
    total_SIM_set0 = fs->make<TH1D>("total_SIM_set0" , "total_SIM_set0" , 400 , 0. , 400.);
    total_SIM_set1 = fs->make<TH1D>("total_SIM_set1" , "total_SIM_set1" , 200 , 0. , 200.);
    total_SIM_set2 = fs->make<TH1D>("total_SIM_set2" , "total_SIM_set2" , 200 , 0. , 200.);

    // energy projected from MIP to SIM_set0
    //total_ENE_odd  = fs->make<TH1D>("total_ENE_odd"  , "total_ENE_odd" , 5000 , 0. , 500.);
    //total_ENE_even = fs->make<TH1D>("total_ENE_even" , "total_ENE_even" , 5000 , 0. , 500.);
    total_ENE_set0 = fs->make<TH1D>("total_ENE_set0" , "total_ENE_set0" , 5000 , 0. , 500.);
    total_ENE_set1 = fs->make<TH1D>("total_ENE_set1" , "total_ENE_set1" , 5000 , 0. , 500.);
    total_ENE_set2 = fs->make<TH1D>("total_ENE_set2" , "total_ENE_set2" , 5000 , 0. , 500.);

    std::ostringstream hnamestr (std::ostringstream::ate);
    for(int i=0;i<26;i++) {
        tb::set_string(hnamestr, "efficiency_linear_track_", i+1);
        efficiency_linear_track_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 100, 0, 1.);

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
        total_SIM_total_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 10.);
        tb::set_string(hnamestr, "multiplicity_digis_total_", i+1);
        multiplicity_digis_total_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 2000.);
        tb::set_string(hnamestr, "multiplicity_simhits_total_", i+1);
        multiplicity_simhits_total_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 200.);
        tb::set_string(hnamestr, "multiplicity_simhits_coarse_", i+1);
        multiplicity_simhits_coarse_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 200.);
        tb::set_string(hnamestr, "multiplicity_simhits_fine_", i+1);
        multiplicity_simhits_fine_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 50, 0, 200.);

        // individual hits
        tb::set_string(hnamestr, "ADC_total_layer_", i+1);
        ADC_total_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 80, 0, 800.);
        tb::set_string(hnamestr, "MIP_total_layer_", i+1);
        MIP_total_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 150, 0, 150.);
        tb::set_string(hnamestr, "SIM_total_layer_", i+1);
        SIM_total_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 100, 0, 1.);
        tb::set_string(hnamestr, "ADC_SimhitE_total_layer_", i+1);
        adc_sim_total_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 800., 0, 800., 0., 800.);
        tb::set_string(hnamestr, "ADC_MIP_total_layer_", i+1);
        adc_mip_total_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 800., 0, 800., 0., 150.);
        tb::set_string(hnamestr, "MIP_SimhitE_total_layer_", i+1);
        mip_sim_total_[i] = fs->make<TProfile>(hnamestr.str().c_str(), hnamestr.str().c_str(), 150., 0, 150., 0., 800.);

        tb::set_string(hnamestr, "nt_total_layer_", i+1);
        nt_total_[i] = fs->make<TNtuple>(hnamestr.str().c_str(),hnamestr.str().c_str(), "adc:mip:simhitE");
    }

#ifdef this_is_an_eventsetup_example
    setupdatatoken_ = esConsumes<setupdata, setuprecord>();
#endif
} 
//}}}
DigiSim::~DigiSim(){}

void DigiSim::reset_tree_variables()
{
    tr_evtNo = 0;
    tr_layerNo = 0;
    tr_signal_region_max_cell = -1;
    tr_signal_region_energy_weighted = -1;
    tr_signal_region_linear_track = -1;
    tr_x = 0.;
    tr_y = 0.;
    tr_z = 0.;
    tr_e = 0.;
    tr_r = 0.;
    tr_d_max_cell = 0.;
    tr_d_energy_weighted = 0.;
    tr_d_linear_track = 0.;
    tr_eta = 0.;
    tr_phi = 0.;
    tr_is_Silicon_w120 = false;
    tr_is_Silicon_w200 = false;
    tr_is_Silicon_w300 = false;
    tr_is_Scintillator = false;
}

void DigiSim::reset_expected_hit_containers(myVector &mv)
{
    mv.tr_vx.clear();
    mv.tr_vy.clear();
    mv.tr_vz.clear();
    mv.tr_ve.clear();

    for(int idx=0; idx<26; ++idx) {
        mv.tr_vx.push_back(0.);
        mv.tr_vy.push_back(0.);
        mv.tr_vz.push_back(0.);
        mv.tr_ve.push_back(0.);
    }
}

void DigiSim::reset_per_event_counters()
{
    counter = 0;
    total_energy_mip_odd  = 0.;
    total_energy_mip_even = 0.;
    total_energy_mip_set0 = 0.;
    total_energy_mip_set1 = 0.;
    total_energy_mip_set2 = 0.;
    total_energy_sim_odd  = 0.;
    total_energy_sim_even = 0.;
    total_energy_sim_set0 = 0.;
    total_energy_sim_set1 = 0.;
    total_energy_sim_set2 = 0.;
    total_corrected_energy_odd  = 0.;
    total_corrected_energy_even = 0.;
    total_corrected_energy_set0 = 0.;
    total_corrected_energy_set1 = 0.;
    total_corrected_energy_set2 = 0.;

    efficiency_numerators   .clear();
    efficiency_denominators .clear();
    efficiency_signal_region_linear_track .clear();

    total_energy_adc_total  .clear();
    total_energy_mip_total  .clear();
    total_energy_mip_coarse .clear();
    total_energy_mip_fine   .clear();
    total_energy_sim_total  .clear();
    num_digis_total         .clear();
    num_simhits_total       .clear();
    num_simhits_coarse      .clear();
    num_simhits_fine        .clear();

    for(int idx=0; idx<26; ++idx) {
        efficiency_numerators  .push_back(0);
        efficiency_denominators.push_back(0);
        efficiency_signal_region_linear_track  .push_back(0);

        total_energy_adc_total .push_back(0);
        total_energy_mip_total .push_back(0);
        total_energy_mip_coarse.push_back(0);
        total_energy_mip_fine  .push_back(0);
        total_energy_sim_total .push_back(0);
        num_digis_total        .push_back(0);
        num_simhits_total      .push_back(0);
        num_simhits_coarse     .push_back(0);
        num_simhits_fine       .push_back(0);
    }
}

double DigiSim::get_corrected_energy_from_dEdx_method(int layer, double amplitude, TString tag)
{
    if(layer>26) printf("[WARNING] get_corrected_energy_from_dEdx_method::layer = %d is outside CEE\n", layer);

    double output = 0.;
    if(tag=="set0")
        output = calibration_weights_set0[layer] * amplitude;
        
    else if(tag=="set1")
        output = calibration_weights_set1[layer] * amplitude;

    else if(tag=="set2")
        output = calibration_weights_set2[layer] * amplitude;

    else
        printf("[WARNING] get_corrected_energy_from_dEdx_method::tag = %s\n is not defined\n", tag.Data());

    return output;
}

double DigiSim::get_additional_correction(int layer)
{
    // no correction
    return 1.;

    //double correction = 1.;
    // start the correction from 3rd layer & consider only odd layers
    //if( layer%2==1 && layer>2 ) {
    //    correction = 1. - 0.5 * (1. - weightsPerLayer_V16[layer] / weightsPerLayer_V16[layer-1] );
    //    return correction;
    //} else {
    //    return correction;
    //}

    if(layer<=26) {
        double width = x_D86[layer] - x_D86[layer-1];
        double dEdx_weight = weightsPerLayer_V16[layer]; 
        //double correction = 1. / (width*dEdx_weight);
        //double correction = dEdx_weight / width;
        //double correction = 1. / width;
        double correction = dEdx_weight;
        return correction;
    } else {
        return 1.;
    }
}

bool DigiSim::is_this_in_set1(int layer)
{
    //E_set1 = E1+E3+E5+E7+E9+E11+E13+E15+E17+...E25
    if(layer==1)       return true;
    else if(layer==3)  return true;
    else if(layer==5)  return true;
    else if(layer==7)  return true;
    else if(layer==9)  return true;
    else if(layer==11) return true;
    else if(layer==13) return true;
    else if(layer==15) return true;
    else if(layer==17) return true;
    else if(layer==19) return true;
    else if(layer==21) return true;
    else if(layer==23) return true;
    else if(layer==25) return true;
    else               return false;
}

bool DigiSim::is_this_in_set2(int layer)
{
    //E_set2 = E1+E3+E5+E8+E10+E12+E14+E15+E17+...E25
    if(layer==1)       return true;
    else if(layer==3)  return true;
    else if(layer==5)  return true;
    else if(layer==8)  return true;
    else if(layer==10) return true;
    else if(layer==12) return true;
    else if(layer==14) return true;
    else if(layer==15) return true;
    else if(layer==17) return true;
    else if(layer==19) return true;
    else if(layer==21) return true;
    else if(layer==23) return true;
    else if(layer==25) return true;
    else               return false;
}

double DigiSim::convert_amplitude_to_total_energy_pedro(int type, double amplitude)
{
    double corrected_energy = 0.;

    // convert E_set0 / E_set1 / E_set2 (corresponding to type 0 / 1 / 2 respectively)
    
    // R80to150
    //if(type==0) corrected_energy = 2.31546e+04 + 4.74438*amplitude;
    //if(type==1) corrected_energy = 2.28611e+04 + 9.03174*amplitude;
    //if(type==2) corrected_energy = 2.36739e+04 + 9.89*amplitude;

    // R80to130
    //if(type==0) corrected_energy = 2.18562e+04 + 4.86902e+00*amplitude;
    //if(type==1) corrected_energy = 2.15429e+04 + 9.28686e+00*amplitude;
    //if(type==2) corrected_energy = 2.22060e+04 + 1.02199e+01*amplitude;

    // R90to130 (MIPs to MeV)
    if(type==0) corrected_energy = 2.17293e+01 + 4.96546e-03*amplitude;
    if(type==1) corrected_energy = 2.15038e+01 + 9.45215e-03*amplitude;
    if(type==2) corrected_energy = 2.22265e+01 + 1.03710e-02*amplitude;

    // R90to130 (Anne-Marie algo + linear track)
    //if(type==0) corrected_energy = 1.70404e+01 + 1.01105e-03*amplitude;
    //if(type==1) corrected_energy = 1.70133e+01 + 1.94038e-03*amplitude;
    //if(type==2) corrected_energy = 1.72140e+01 + 2.11144e-03*amplitude;

    return corrected_energy;
}

GlobalPoint DigiSim::projectHitPositionAt(float z,float eta,float phi)
{
  float theta=2*TMath::ATan(exp(-eta));
  float rho=z*TMath::Tan(theta);
  GlobalPoint xyz(rho*TMath::Cos(phi),rho*TMath::Sin(phi),z);
  return xyz;
}

float DigiSim::get_distance_from_expected_hit(double x, double y, double z, double eta, double phi)
{
    TVector2 xy(x,y);
    GlobalPoint xyzExp = projectHitPositionAt(z, eta, phi);
    TVector2 xyExp(xyzExp.x(),xyzExp.y());
    float d = (xyExp-xy).Mod();
    return d;
}

float DigiSim::get_distance_from_expected_hit(double x, double y, double x0, double y0)
{
    TVector2 xy(x,y);
    TVector2 xyExp(x0,y0);
    float d = (xyExp-xy).Mod();
    return d;
}

int DigiSim::get_signal_region(float d)
{
    int signal_region = -1;
    if(d<=1.3)      signal_region = 1;
    else if(d<=2.6) signal_region = 2;
    else if(d<=5.3) signal_region = 3;
    else            signal_region = -1;
    return signal_region;
}

void DigiSim::calculate_efficiency()
{
    for(int idx=0; idx<26; ++idx) {
        if(efficiency_denominators[idx] > 0.)
            efficiency_signal_region_linear_track[idx] = efficiency_numerators[idx] / efficiency_denominators[idx];
        else
            efficiency_signal_region_linear_track[idx] = -1.;
    }
}

void DigiSim::fill_event_info()
{
    // Fill information of an event
    for(int idx=0; idx<26; ++idx) {
        //tb::print_debug_info("num_digis_total"  , num_digis_total  [idx]        );
        efficiency_linear_track_ [idx] -> Fill( efficiency_signal_region_linear_track[idx] );

        total_ADC_total_             [idx] -> Fill( total_energy_adc_total  [idx] );
        total_MIP_total_             [idx] -> Fill( total_energy_mip_total  [idx] );
        total_MIP_coarse_            [idx] -> Fill( total_energy_mip_coarse [idx] );
        total_MIP_fine_              [idx] -> Fill( total_energy_mip_fine   [idx] );
        total_SIM_total_             [idx] -> Fill( total_energy_sim_total  [idx] );
        multiplicity_digis_total_    [idx] -> Fill( num_digis_total         [idx] );
        multiplicity_simhits_total_  [idx] -> Fill( num_simhits_total       [idx] );
        multiplicity_simhits_coarse_ [idx] -> Fill( num_simhits_coarse      [idx] );
        multiplicity_simhits_fine_   [idx] -> Fill( num_simhits_fine        [idx] );

    }

    total_MIP_odd  -> Fill(total_energy_mip_odd);
    total_MIP_even -> Fill(total_energy_mip_even);
    total_MIP_set0 -> Fill(total_energy_mip_set0);
    total_MIP_set1 -> Fill(total_energy_mip_set1);
    total_MIP_set2 -> Fill(total_energy_mip_set2);

    total_SIM_odd  -> Fill(total_energy_sim_odd);
    total_SIM_even -> Fill(total_energy_sim_even);
    total_SIM_set0 -> Fill(total_energy_sim_set0);
    total_SIM_set1 -> Fill(total_energy_sim_set1);
    total_SIM_set2 -> Fill(total_energy_sim_set2);

    // linear fit (not preferred because deposited energy in passive layers is not considered)
    //total_corrected_energy_set0 = convert_amplitude_to_total_energy_pedro(0, total_energy_mip_set0);
    //total_corrected_energy_set1 = convert_amplitude_to_total_energy_pedro(1, total_energy_mip_set1);
    //total_corrected_energy_set2 = convert_amplitude_to_total_energy_pedro(2, total_energy_mip_set2);

    // store energy in unit of GeV instead of MeV
    //total_ENE_odd  -> Fill(total_corrected_energy_odd);
    //total_ENE_even -> Fill(total_corrected_energy_even);
    total_ENE_set0 -> Fill(total_corrected_energy_set0 / 1000.);
    total_ENE_set1 -> Fill(total_corrected_energy_set1 / 1000.);
    total_ENE_set2 -> Fill(total_corrected_energy_set2 / 1000.);

    //tb::print_debug_info("total_energy_sim_set0", total_energy_sim_set0);
    //tb::print_debug_info("total_energy_mip_set1", total_energy_mip_set1);
    //tb::print_debug_info("total_corrected_energy_set1", total_corrected_energy_set1, true);
}
