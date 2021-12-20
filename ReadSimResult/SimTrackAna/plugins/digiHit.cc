#include <memory>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <iostream>

// user include files
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
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"

#include "CoralBase/Exception.h"

#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToModule.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToROC.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"

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
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"

//#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
//#include "DQMServices/Core/interface/DQMStore.h"

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

class DigiSim : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
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
      u_cor = v_cor = type = layer = eta = 0;
      hitid = nhits = 0;
    }
    int u_cor, v_cor, type, layer,eta;
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
  hgcal::RecHitTools rhtools_;
  edm::EDGetToken digiSource_;
  //edm::ConsumesCollector iC;

  TH1D *hELossEE;TH1D *hELossEEF;TH1D *hELossEECN;TH1D *hELossEECK;
  TH1D *hELossHEF;TH1D *hELossHEFF;TH1D *hELossHEFCN;TH1D *hELossHEFCK;
  std::vector<TH1D*> vechist;   
  TH1D *ADC_120mum_[26];
  TH1D *ADC_200mum_[26];
  TH1D *ADC_300mum_[26];
  TProfile *adc_Simhit_120mum_[26];
  TProfile *adc_Simhit_200mum_[26];
  TProfile *adc_Simhit_300mum_[26];
  TProfile *Simhit_adc_120mum_[26];
  TProfile *Simhit_adc_200mum_[26];
  TProfile *Simhit_adc_300mum_[26];
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif

};


DigiSim::DigiSim(const edm::ParameterSet& iconfig)
  :
  //auto temp = iConfig.getUntrackedParameter<edm::InputTag>("digihits");
  nameDetector_(iconfig.getParameter<std::string>("Detector")),
  ifNose_(iconfig.getUntrackedParameter<bool>("ifNose")),
  verbosity_(iconfig.getUntrackedParameter<int>("Verbosity", 0)),
  SampleIndx_(iconfig.getUntrackedParameter<int>("SampleIndx", 0)),
  tok_hgcalg_(esConsumes<HGCalGeometry, IdealGeometryRecord>(edm::ESInputTag{"", nameDetector_})), 
  firstLayer_(1), 
  tSimCaloHitContainer(consumes<edm::PCaloHitContainer>(iconfig.getUntrackedParameter<edm::InputTag>("simhits")))
{
  
  auto temp = iconfig.getUntrackedParameter<edm::InputTag>("digihits");
  if ((nameDetector_ == "HGCalEESensitive") || (nameDetector_ == "HGCalHESiliconSensitive") ||
      (nameDetector_ == "HGCalHEScintillatorSensitive") || (nameDetector_ == "HGCalHFNoseSensitive")) {
      digiSource_ = consumes<HGCalDigiCollection>(temp);
      //digiSource_=consumes<HGCalDigiCollection>(iconfig.getUntrackedParameter<edm::InputTag>("digihits"));
      //tSimCaloHitContainer=consumes<edm::PCaloHitContainer>(iconfig.getUntrackedParameter<edm::InputTag>("simhits"));
  } 
  else {
    throw cms::Exception("BadHGCDigiSource") << "HGCal DetectorName given as " << nameDetector_ << " must be: "
                                             << "\"HGCalEESensitive\", \"HGCalHESiliconSensitive\", or "
                                             << "\"HGCalHEScintillatorSensitive\", \"HGCalHFNoseSensitive\"!";
  }
  //}

 // tSimCaloHitContainer(consumes<edm::PCaloHitContainer>(iconfig.getUntrackedParameter<edm::InputTag>("simhits")))
//{
  //now do what ever initialization is needed
  //name = iconfig.getParameter<std::string>("Detector");
  usesResource("TFileService");
  edm::Service<TFileService> fs; 
  hELossEE = fs->make<TH1D>("hELossEE","hELossEE", 1000, 0., 1000.); 
  hELossEEF = fs->make<TH1D>("hELossEEF","hELossEEF", 1000, 0., 1000.); 
  hELossEECN = fs->make<TH1D>("hELossEECN","hELossEECN", 1000, 0., 1000.); 
  hELossEECK = fs->make<TH1D>("hELossEECK","hELossEECK", 1000, 0., 1000.);
  hELossHEF = fs->make<TH1D>("hELossHEF","hELossHEF", 1000, 0., 1000.);
  hELossHEFF = fs->make<TH1D>("hELossHEFF","hELossHEFF", 1000, 0., 1000.);
  hELossHEFCN = fs->make<TH1D>("hELossHEFCN","hELossHEFCN", 1000, 0., 1000.);
  hELossHEFCK = fs->make<TH1D>("hELossHEFCK","hELossHEFCK", 1000, 0., 1000.);
  std::ostringstream hnamestr (std::ostringstream::ate);
  for(int i=0;i<26;i++){
    hnamestr.str("ADC_120mum_layer_");
    hnamestr<<i+1;
    ADC_120mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 1024., 0, 1024.);
    hnamestr.clear();

    hnamestr.str("ADC_SimhitE_120mum_layer_");
    hnamestr<<i+1;
    adc_Simhit_120mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),1024., 0, 1024.,0., 1000.);
    hnamestr.clear();

    hnamestr.str("SimhitE_ADC_120mum_layer_");
    hnamestr<<i+1;
    Simhit_adc_120mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),1000., 0., 1000.,0, 1024.);
    hnamestr.clear();


    hnamestr.str("ADC_200mum_layer_");
    hnamestr<<i+1;
    ADC_200mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 1024., 0, 1024.);
    hnamestr.clear();

    hnamestr.str("ADC_SimhitE_200mum_layer_");
    hnamestr<<i+1;
    adc_Simhit_200mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),1024., 0, 1024.,0., 1000.);
    hnamestr.clear();

    hnamestr.str("SimhitE_ADC_200mum_layer_");
    hnamestr<<i+1;
    Simhit_adc_200mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),1000., 0., 1000.,0, 1024.);
    hnamestr.clear();


    hnamestr.str("ADC_300mum_layer_");
    hnamestr<<i+1;
    ADC_300mum_[i] = fs->make<TH1D>(hnamestr.str().c_str(),hnamestr.str().c_str(), 1024., 0, 1024.);
    hnamestr.clear();		    

    hnamestr.str("ADC_SimhitE_300mum_layer_");
    hnamestr<<i+1;
    adc_Simhit_300mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),1024., 0, 1024.,0., 1000.);
    hnamestr.clear();

    hnamestr.str("SimhitE_ADC_300mum_layer_");
    hnamestr<<i+1;
    Simhit_adc_300mum_[i] = fs->make<TProfile>(hnamestr.str().c_str(),hnamestr.str().c_str(),1000., 0., 1000.,0, 1024.);
    hnamestr.clear();


  }
  //ADC_ = fs->make<TH1D>("ADC_","ADC", 1024., 0, 1024.);
  //for(int i=0; i<5;i++){
  //vechist.push_back(iB.book1D(histoname.str().c_str(), "ADCDigiOccupancy_300mum", 1024, 0, 1024));

#ifdef this_is_an_eventsetup_example
  setupdatatoken_ = esConsumes<setupdata, setuprecord>();
#endif
}

DigiSim::~DigiSim()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}















//
// member functions
//

// ------------ method called for each event  ------------
void
DigiSim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout<<"----------------------------------  new track  -----------------"<<std::endl;
  using namespace edm;
  //int geomType(0);
  //const HGCalGeometry* geom0 = &iSetup.getData(tok_hgcalg_);
  //std::cout<<geom0->topology().waferHexagon8()<<std::endl;
////////////////////////////////////////////////////----------- Simhit Hnadle -----------------///////////////
  std::map<uint32_t, std::pair<hitsinfo, energysum> > map_Simhits;
  map_Simhits.clear();
  unsigned int nofSiHits = 0;
  Handle<PCaloHitContainer> simhit;
  iEvent.getByToken(tSimCaloHitContainer, simhit);
  for(PCaloHitContainer::const_iterator itHit= simhit->begin(); itHit!= simhit->end(); ++itHit) {
    DetId detId = static_cast<DetId>(itHit->id());
    if(rhtools_.isSilicon(detId)){
      HGCSiliconDetId id(itHit->id());
      if(nameDetector_ == "HGCalEESensitive"){
        hELossEE->Fill(itHit->energy()*1.e6);
        if(id.type()==HGCSiliconDetId::HGCalFine)
	  hELossEEF->Fill(itHit->energy()*1.e6); //in keV
        if(id.type()==HGCSiliconDetId::HGCalCoarseThin)
	  hELossEECN->Fill(itHit->energy()*1.e6); //in keV
        if(id.type()==HGCSiliconDetId::HGCalCoarseThick)
          hELossEECK->Fill(itHit->energy()*1.e6); //in keV
      }

      if(nameDetector_ == "HGCalHESiliconSensitive"){
        hELossHEF->Fill(itHit->energy()*1.e6);
        if(id.type()==HGCSiliconDetId::HGCalFine)
	  hELossHEFF->Fill(itHit->energy()*1.e6); //in keV
        if(id.type()==HGCSiliconDetId::HGCalCoarseThin)
	  hELossHEFCN->Fill(itHit->energy()*1.e6); //in keV
        if(id.type()==HGCSiliconDetId::HGCalCoarseThick)
	  hELossHEFCK->Fill(itHit->energy()*1.e6); //in keV
      }  
    
      //std::cout<<" hit energy = "<<itHit->energy()<<std::endl;
      uint32_t id_ = itHit->id();
      double energy = itHit->energy();
      energysum esum;
      hitsinfo hinfo;
    	
      if (map_Simhits.count(id_) != 0) {
        hinfo = map_Simhits[id_].first;
        esum = map_Simhits[id_].second;
      } 
      else {
        hinfo.hitid =  nofSiHits;
        hinfo.u_cor = rhtools_.getCell(detId).first ;
        hinfo.v_cor = rhtools_.getCell(detId).second ;
        hinfo.type = id.type();
        hinfo.layer = rhtools_.getLayerWithOffset(detId);
	//hinfo.eta = itHit->momentum().eta()
      }	
      esum.etotal += energy;
      esum.eTime[0] = energy;
      map_Simhits[id_] = std::pair<hitsinfo, energysum>(hinfo, esum);
    }
  }
     //std::cout<<"map size = "<< map_Simhits.size()<<std::endl;
 







//////////////////////////////////////////////////////------------ Digi Handle ------------------////////////
  std::map<uint32_t, std::pair<digisinfo,adcinfo > > map_digihits;
  map_digihits.clear();
  Handle<HGCalDigiCollection> digicollection;
  iEvent.getByToken(digiSource_, digicollection);
  if (digicollection.isValid()) {
	//std::cout<<"valid"<<std::endl;
    for (const auto& it : *(digicollection.product())) {
      DetId detId = it.id();
      //std::cout<<"isSilicon = "<<rhtools_.isSilicon(detId)<<std::endl;
      if(rhtools_.isSilicon(detId)){
        //std::cout<<"entered in the condition"<<std::endl;
        uint32_t id_digi = uint32_t(it.id());
        //  int layer = ((geomType == 0)   ? HGCalDetId(detId).layer()
        //               : (geomType == 1) ? HGCSiliconDetId(detId).layer()
        //                                 : HFNoseDetId(detId).layer());
        //int waferType = ((geomType == 1) ? HGCSiliconDetId(detId).type()
        //                                   : HFNoseDetId(detId).type());
        const HGCSample& hgcSample = it.sample(SampleIndx_);
        //uint16_t gain = hgcSample.toa();
        uint16_t adc_ = hgcSample.data();
	
        digisinfo dinfo;
        adcinfo ainfo;
        if (map_digihits.count(id_digi) != 0) {
          dinfo = map_digihits[id_digi].first;
          ainfo = map_digihits[id_digi].second;
        }
        else {
          dinfo.u_cor = HGCSiliconDetId(detId).cellU() ;
          dinfo.v_cor = HGCSiliconDetId(detId).cellV() ;
          dinfo.type =  HGCSiliconDetId(detId).type();
          dinfo.layer = HGCSiliconDetId(detId).layer();
          ainfo.adc = adc_;
        }


        map_digihits[id_digi] = std::pair<digisinfo, adcinfo>(dinfo, ainfo);
	//std::cout<<" Id_digi : "<<id_digi<<"; adc_ : "<<ainfo.adc<<std::endl;
        //double charge = gain;
        //bool totmode = hgcSample.mode();
        //bool zerothreshold = hgcSample.threshold();
	//std::cout<<" layer = "<<HGCalDetId(detId).layer()<< " gain = "<<gain<<" adc "<<adc<<std::endl;
        //digiValidation(detId, geom0, layer, waferType, adc, charge, totmode, zerothreshold);
      }
    }
  }
////////////////////////////////////////////////////////////////////////////////////////////




  std::map<uint32_t, std::pair<hitsinfo, energysum> > map_Simhits_matched;
  std::map<uint32_t, std::pair<digisinfo, adcinfo>>::iterator itr_digi;
  std::map<uint32_t, std::pair<hitsinfo, energysum> >::iterator itr_sim;
  Double_t max_energy=0; 
  uint32_t matched_id=0;
  for (itr_sim = map_Simhits.begin(); itr_sim != map_Simhits.end(); ++itr_sim) {
    energysum esum = (*itr_sim).second.second;
    hitsinfo simhitinfo = (*itr_sim).second.first;
    if(simhitinfo.layer==1){
        if(max_energy<esum.etotal*1.e6){
	  max_energy=esum.etotal*1.e6;
	  if(matched_id!=0){
		map_Simhits_matched.erase(matched_id);
	  }
	  matched_id = (*itr_sim).first;
	  map_Simhits_matched[matched_id] = std::pair<hitsinfo, energysum>(simhitinfo,esum);
        }
    }
  }

  /*std::map<uint32_t, std::pair<hitsinfo, energysum> >::iterator itr_sim_matched;
  for (itr_sim_matched = map_Simhits_matched.begin(); itr_sim_matched != map_Simhits_matched.end(); ++itr_sim_matched) {
        //energysum esum = (*itr_sim).second.second;
        //hitsinfo simhitinfo = (*itr_sim).second.first;
	std::cout<<"id_digi from new map = "<<(*itr_sim_matched).first<<std::endl;
  }*/
  //std::cout<<"for loop for id removel end "<<std::endl;

  for (itr_sim = map_Simhits.begin(); itr_sim != map_Simhits.end(); ++itr_sim) {
    energysum esum = (*itr_sim).second.second;
    //hitsinfo simhitinfo = (*itr_sim).second.first;
    if(esum.etotal>0 && (*itr_sim).first==matched_id){//simhitinfo.layer==1){// && esum.etotal*1.e6==max_energy){
       //std::cout<<max_energy<<std::endl; 
       //std::cout<<"maxenergy = "<<esum.etotal*1.e6<<std::endl;
       for (itr_digi = map_digihits.begin(); itr_digi != map_digihits.end(); ++itr_digi) {
	  digisinfo dinfo = (*itr_digi).second.first;
          adcinfo ainfo = (*itr_digi).second.second;
	  if((*itr_sim).first==(*itr_digi).first){
                std::cout<<"id_digi = "<<(*itr_digi).first<<" adc = "<<ainfo.adc<<" wafer type = "<<dinfo.type<<" layer = "<<dinfo.layer<<std::endl;
	     if(dinfo.layer < 26 && dinfo.type==0){
                //std::cout<<"id_digi = "<<(*itr_digi).first<<" adc = "<<ainfo.adc<<" wafer type = "<<dinfo.type<<" layer = "<<dinfo.layer<<std::endl;
	        ADC_120mum_[dinfo.layer-1]->Fill(ainfo.adc);
		adc_Simhit_120mum_[dinfo.layer-1]->Fill(ainfo.adc,esum.eTime[0]*1.e6);
		Simhit_adc_120mum_[dinfo.layer-1]->Fill(esum.eTime[0]*1.e6,ainfo.adc);
	     }
	     if(dinfo.layer < 26 && dinfo.type==1){
                //std::cout<<"id_digi = "<<(*itr_digi).first<<" adc = "<<ainfo.adc<<" wafer type = "<<dinfo.type<<" layer = "<<dinfo.layer<<std::endl;
	        ADC_200mum_[dinfo.layer-1]->Fill(ainfo.adc);
		adc_Simhit_200mum_[dinfo.layer-1]->Fill(ainfo.adc,esum.eTime[0]*1.e6);
                Simhit_adc_200mum_[dinfo.layer-1]->Fill(esum.eTime[0]*1.e6,ainfo.adc);
	     }
	     if(dinfo.layer < 26 && dinfo.type==2){
                //std::cout<<"id_digi = "<<(*itr_digi).first<<" adc = "<<ainfo.adc<<" wafer type = "<<dinfo.type<<" layer = "<<dinfo.layer<<std::endl;
	        ADC_300mum_[dinfo.layer-1]->Fill(ainfo.adc);
		adc_Simhit_300mum_[dinfo.layer-1]->Fill(ainfo.adc,esum.eTime[0]*1.e6);
                Simhit_adc_300mum_[dinfo.layer-1]->Fill(esum.eTime[0]*1.e6,ainfo.adc);
	     }
          }
       }
    }
  }

  //std::cout<<"for loop for sim digi maching end "<<std::endl;





  #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
  #endif
}  




// ------------ method called once each job just before starting event loop  ------------
void
DigiSim::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DigiSim::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DigiSim::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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


//define this as a plug-in
DEFINE_FWK_MODULE(DigiSim);
