// -*- C++ -*-
//
// Package:    Demo/simHit
// Class:      simHit
//
/**\class simHit TrackAnalyzer.cc Track/TrackAnalyzer/plugins/TrackAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Indranil Das
//         Created:  Wed, 25 Aug 2021 06:18:11 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <fstream>
#include <string>

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
#include <TGraph.h>
#include <TMath.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class simHit : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  //Implemented following Validation/HGCalValidation/plugins/HGCalSimHitValidation.cc
  
  explicit simHit(const edm::ParameterSet&);
  ~simHit();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  //edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<edm::SimTrackContainer> tSimTrackContainer; 
  edm::EDGetTokenT<edm::PCaloHitContainer> tSimCaloHitContainer; 
  std::string name;
  //edm::ESGetToken<HGCalGeometry, IdealGeometryRecord> geomToken_;
  

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_; 
  hgcal::RecHitTools rhtools_;
  //edm::ConsumesCollector iC;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif

};

// constructors and destructor
//
simHit::simHit(const edm::ParameterSet& iconfig)
  :
  tSimTrackContainer(consumes<edm::SimTrackContainer>(iconfig.getUntrackedParameter<edm::InputTag>("simtrack"))),
  tSimCaloHitContainer(consumes<edm::PCaloHitContainer>(iconfig.getUntrackedParameter<edm::InputTag>("simhits")))
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs; 

  caloGeomToken_ = esConsumes<CaloGeometry, CaloGeometryRecord>();

#ifdef this_is_an_eventsetup_example
  setupdatatoken_ = esConsumes<setupdata, setuprecord>();
#endif
}


simHit::~simHit()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
simHit::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  Handle<SimTrackContainer> simtrack;
  iEvent.getByToken(tSimTrackContainer, simtrack);
  int itrk = 0;
  //double muonpt = 0.0;
  SimTrackContainer::const_iterator itTrack;
  for(itTrack = simtrack->begin(); itTrack != simtrack->end(); ++itTrack) {
    int charge = itTrack->charge();  
    std::cout << "charge " << charge <<std::endl;
    if(!itTrack->noGenpart()){
      std::cout << "pt " << itTrack->momentum().pt() << ", eta "<< itTrack->momentum().eta() << ", phi" << itTrack->momentum().phi() <<std::endl;
    }
      std::cout << "PDG " << itTrack->type() <<std::endl;
    if(itTrack->noGenpart())
      std::cout << "pt NoGen " << itTrack->momentum().pt() <<std::endl;
  itrk++;
  }
  #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
     auto setup = iSetup.getData(setupToken_);
  //     // if need the ESHandle to check if the SetupData was there or not
     auto pSetup = iSetup.getHandle(setupToken_);
         #endif
}  




// ------------ method called once each job just before starting event loop  ------------
void
simHit::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
simHit::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
simHit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(simHit);
