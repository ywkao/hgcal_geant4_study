// -*- C++ -*-
//
// Package:    clue_performance/clue_performance
// Class:      clue_performance
//
/**\class clue_performance clue_performance.cc clue_performance/clue_performance/plugins/clue_performance.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zheng-Gang Chen
//         Created:  Tue, 29 Jun 2021 07:55:39 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"


#include "TTree.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
const int kmax =10000;

using reco::TrackCollection;
using reco::PhotonCollection;
using reco::PFCandidateCollection;
class clue_performance : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit clue_performance(const edm::ParameterSet&);
  ~clue_performance();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  int evt =0;
  hgcal::RecHitTools tool;
  TTree *mytree = new TTree("ntuple","ntuple");
  // ----------member data ---------------------------
  int layercluster_number;
  int thickness_type[kmax];
  int layercluster_layer[kmax];
  double layercluster_energy[kmax];
  double layercluster_phi[kmax];
  double layercluster_eta[kmax];
  double layercluster_z[kmax];
  // edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  
  edm::EDGetTokenT<reco::CaloClusterCollection> layerclustersToken_;
  edm::EDGetTokenT<reco::PhotonCollection> photonsToken_;
  edm::EDGetTokenT<std::vector<reco::HGCalMultiCluster>> multisToken_;
  edm::EDGetTokenT< reco::PFCandidateCollection> pfcandsticlToken_;
  //edm::EDGetTokenT< std::vector<TICLCandidate>> pfcandsticlToken_;
  edm::EDGetTokenT< reco::PFCandidateCollection> pfcandsToken_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> trackstersToken_;
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
clue_performance::clue_performance(const edm::ParameterSet& iConfig)
      :layerclustersToken_(consumes<reco::CaloClusterCollection>(iConfig.getUntrackedParameter<edm::InputTag>("layerclusters"))),
      photonsToken_(consumes<reco::PhotonCollection>(iConfig.getUntrackedParameter<edm::InputTag>("photons"))),
      multisToken_(consumes<std::vector<reco::HGCalMultiCluster>>(iConfig.getUntrackedParameter<edm::InputTag>("multis"))),
      pfcandsticlToken_(consumes<reco::PFCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfcandsticl"))),
      pfcandsToken_(consumes<reco::PFCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfcands"))),
      trackstersToken_(consumes<std::vector<ticl::Trackster>>(iConfig.getUntrackedParameter<edm::InputTag>("tracksters")))
      
{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

clue_performance::~clue_performance() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void clue_performance::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  ++evt;

  // RecHit Collection
  Handle<reco::PCaloHitCollection> hgcalSimHitEE;
  iEvent.getByToken(hgcalSimHitEEToken_, hgcalSimHitEE);

  // LayerCluster Collection
  Handle<reco::CaloClusterCollection> layerclusters;
  iEvent.getByToken(layerclustersToken_,layerclusters);
  layercluster_number = 0;
  for(reco::CaloClusterCollection::const_iterator c = layerclusters->begin(); c != layerclusters->end(); c++) {
      if (TMath::Abs(c->z())>367.699 && tool.getLayer(c->seed())<28) {
          layercluster_layer[layercluster_number] = tool.getLayer(c->seed())+28;
      } else {
          layercluster_layer[layercluster_number] = tool.getLayer(c->seed());
      }
      thickness_type[layercluster_number] = tool.getSiThickIndex(c->seed());
      layercluster_energy[layercluster_number] = c->energy();
      layercluster_phi[layercluster_number] = c->phi();
      layercluster_eta[layercluster_number] = c->eta();
      layercluster_z[layercluster_number] = c->z();
      layercluster_number++;
  }
  mytree->Fill();

//vector<reco::CaloCluster>             "hgcalLayerClusters"        ""                "RECO"    
//vector<PCaloHit>                      "g4SimHits"                 "HGCHitsEE"       "SIM"     
//edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >    "HGCalRecHit"               "HGCEERecHits"    "RECO"    
//edm::SortedCollection<HGCDataFrame<DetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<DetId,HGCSample> > >    "simHGCalUnsuppressedDigis"   "HEfront"         "HLT"     

  /*
  Handle<reco::PFCandidateCollection> ticlcands;
  iEvent.getByToken(pfcandsticlToken_,ticlcands);
  std::cout<<"pfTICL Collection:"<<std::endl;
  for(reco::PFCandidateCollection::const_iterator cand = ticlcands->begin(); cand != ticlcands->end(); cand++){
     std::cout<<"PDGID : "<<cand->pdgId()<<std::endl;
     std::cout<<"Energy: "<<cand->energy()<<std::endl;
     std::cout<<"phi   : "<<cand->phi()<<std::endl;
     std::cout<<"eta   : "<<cand->eta()<<std::endl;
     std::cout<<std::endl;
  }
  Handle<std::vector<ticl::Trackster>> tracksters;
  iEvent.getByToken(trackstersToken_,tracksters);
  std::cout<<"tracksters Collection:"<<std::endl;
  for(std::vector<ticl::Trackster>::const_iterator cand = tracksters->begin(); cand != tracksters->end(); cand++){
     std::cout<<"PDGID : "<<cand->pdgId()<<std::endl;
     std::cout<<"Energy: "<<cand->energy()<<std::endl;
     std::cout<<"phi   : "<<cand->phi()<<std::endl;
     std::cout<<"eta   : "<<cand->eta()<<std::endl;
     std::cout<<std::endl;
  }

  Handle<reco::PFCandidateCollection> cands;
  iEvent.getByToken(pfcandsToken_,cands);
  std::cout<<"ParticleFlow Collection:"<<std::endl;
  for(reco::PFCandidateCollection::const_iterator cand = cands->begin(); cand != cands->end(); cand++){
     if(!cand->elementsInBlocks().empty()){
        reco::PFCandidate::ElementsInBlocks eleInBlocks = cand->elementsInBlocks();
        for(unsigned i = 0; i < eleInBlocks.size();i++)
        {
            reco::PFBlockRef& blockRef = eleInBlocks[i].first;
            
            std::cout<<"PFBlockRef: "<<i<<std::endl;
            for(const auto& elem : blockRef->elements())
            {
                std::cout<<"PFBlockElements: "<<elem.index()<<std::endl;
                std::cout<<"Type           :"<<elem.type()<<std::endl;
                
                reco::PFClusterRef clusterRef = elem.clusterRef();
                std::cout<<"Layer: "<<clusterRef.get()->layer()<<std::endl;
                std::cout<<"Eta  : "<<clusterRef.get()->eta()<<std::endl;
            }
        }
     }

     else
     {
         std::cout<<"No PFBlock!"<<std::endl;
         std::cout<<"PDGID : "<<cand->pdgId()<<std::endl;
         std::cout<<"Energy: "<<cand->energy()<<std::endl;
         std::cout<<"phi   : "<<cand->phi()<<std::endl;
         std::cout<<"eta   : "<<cand->eta()<<std::endl;
         std::cout<<std::endl;
     }
  }
  Handle<reco::PhotonCollection> photons;
  iEvent.getByToken(photonsToken_,photons);
  int ph =1;
  std::cout<<std::endl;
  std::cout<<"Reco::PhotonCollection: "<<std::endl;
  for(reco::PhotonCollection::const_iterator p = photons->begin(); p != photons->end(); p++)
  {
      std::cout<<"Photon: "<<ph<<std::endl;
      std::cout<<"Energy: "<<p->energy()<<std::endl;
      std::cout<<"Phi   : "<<p->phi()<<std::endl;
      std::cout<<"Eta   : "<<p->eta()<<std::endl;
      std::cout<<"LayerCluster in SuperCluster"<<std::endl;
      std::cout<<"Number of LayerCluster: "<<p->superCluster()->clustersSize()<<std::endl;
      for(auto cluster : p->superCluster()->clusters())
      {
        std::cout<<"Layer : "<<tool.getLayer(cluster->seed())<<std::endl;
        std::cout<<"Energy: "<<cluster->energy()<<std::endl;
        std::cout<<"Phi   : "<<cluster->phi()<<std::endl;
        std::cout<<"Eta   : "<<cluster->eta()<<std::endl;
        std::cout<<"AlgoID: "<<cluster->algoID()<<std::endl;
        std::cout<<"Z     : "<<cluster->z()<<std::endl;
      }
      std::cout<<std::endl;
      ph++;
  }
  Handle<std::vector<reco::HGCalMultiCluster>> multis;
  iEvent.getByToken(multisToken_,multis);
  int mc = 0 ;
  double energy = 0;
  for(auto m= multis->begin();m!=multis->end();m++)
  {
      std::cout<<"Multilcluster: "<<mc<<std::endl;
      std::cout<<"Eta    :"<<m->eta()<<std::endl;
      std::cout<<"Phi    :"<<m->phi()<<std::endl;
      std::cout<<"Energy :"<<m->energy()<<std::endl;
      energy += m->energy();
      for(auto cluster: m->clusters())
      {
            std::cout<<"AlgoID: "<<cluster->algoID()<<std::endl;
            std::cout<<"Energy: "<<cluster->energy()<<std::endl;
            std::cout<<"Phi   : "<<cluster->phi()<<std::endl;
            std::cout<<"Eta   : "<<cluster->eta()<<std::endl;
      }
      mc++;
  }
  std::cout<<"Total Energy: "<<energy<<std::endl;
  std::cout<<std::endl;
  */

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void clue_performance::beginJob() {
    mytree->Branch("layercluster_number",&layercluster_number,"layercluster_number/I");
    mytree->Branch("thickness_type",thickness_type,"thickness_type[layercluster_number]/I");
    mytree->Branch("layercluster_energy",layercluster_energy,"layercluster_energy[layercluster_number]/D");
    mytree->Branch("layercluster_phi",layercluster_phi,"layercluster_phi[layercluster_number]/D");
    mytree->Branch("layercluster_eta",layercluster_eta,"layercluster_eta[layercluster_number]/D");
    mytree->Branch("layercluster_layer",layercluster_layer,"layercluster_layer[layercluster_number]/I");
}

// ------------ method called once each job just after ending the event loop  ------------
void clue_performance::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void clue_performance::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(clue_performance);
