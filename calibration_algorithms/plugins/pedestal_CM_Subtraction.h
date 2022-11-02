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

class Calibration : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit Calibration(const edm::ParameterSet&);
        ~Calibration();

        virtual void InitTree(TTree *tree=0);
        virtual Long64_t LoadTree(Long64_t entry);
        void reset_tree_variables();
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

        TTree *fChain;   //!pointer to the analyzed TTree or TChain
        Int_t  fCurrent; //!current Tree number in a TChain

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;


        // load test beam data
        int event;
        int corruption;
        int chip;
        int half;
        int channel;
        int adc;
        int toa;
        int tot;
        int trigtime;


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

        edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
        hgcal::RecHitTools rhtools_;
        edm::EDGetToken digiSource_;
        edm::ESGetToken<HGCalGeometry, IdealGeometryRecord> ee_geometry_token_;
        HGCalUncalibRecHitRecWeightsAlgo<HGCalDataFrame> uncalibMaker_ee_;

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

Calibration::Calibration(const edm::ParameterSet& iconfig) : //{{{
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

    InitTree(0);

#ifdef this_is_an_eventsetup_example
    setupdatatoken_ = esConsumes<setupdata, setuprecord>();
#endif
} 
//}}}
Calibration::~Calibration()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void Calibration::InitTree(TTree *tree)
{
    //*** the code in this function is the same as the following two lines ***//
    // TFile *file = TFile::Open(rootfile.Data(), "R");
    // TTree *raw_tree = (TTree*)file->Get("unpacker_data/hgcroc");

    // load tree
    if (tree == 0) {
       //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/group/dpg_hgcal/tb_hgcal/2022/sps_oct2022/pion_beam_150_320fC/beam_run/run_20221011_130925/beam_run0.root");
       //if (!f || !f->IsOpen()) {
       //   f = new TFile("/eos/cms/store/group/dpg_hgcal/tb_hgcal/2022/sps_oct2022/pion_beam_150_320fC/beam_run/run_20221011_130925/beam_run0.root");
       //}
       TString rootfile = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2022/sps_oct2022/pion_beam_150_320fC/beam_run/run_20221011_130925/beam_run0.root";
       TFile *f = TFile::Open(rootfile.Data(), "R");
       TDirectory * dir = (TDirectory*)f->Get("unpacker_data");
       dir->GetObject("hgcroc",tree);
    }

    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("event", &event);
    fChain->SetBranchAddress("corruption", &corruption);
    fChain->SetBranchAddress("chip", &chip);
    fChain->SetBranchAddress("half", &half);
    fChain->SetBranchAddress("channel", &channel);
    fChain->SetBranchAddress("adc", &adc);
    fChain->SetBranchAddress("toa", &toa);
    fChain->SetBranchAddress("tot", &tot);
    fChain->SetBranchAddress("trigtime", &trigtime);
}

Long64_t Calibration::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}

void Calibration::reset_tree_variables()
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

