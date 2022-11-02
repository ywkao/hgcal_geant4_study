#include "pedestal_CM_Subtraction.h"

// ------------ method called for each event  ------------ //
void Calibration::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // init {{{
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

    const CaloGeometry &geomCalo = iSetup.getData(caloGeomToken_);
    rhtools_.setGeometry(geomCalo);
    
    //----------------------------------------------------------------------------------------------------}}}
    
    //*** loading -> to be placed in constructor ***//
    //TString rootfile = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2022/sps_oct2022/pion_beam_150_320fC/beam_run/run_20221011_130925/beam_run0.root";
    //TString treename = "unpacker_data/hgcroc";
    //printf("Input file: %s\n", rootfile.Data());
    //TFile *file = TFile::Open(rootfile.Data(), "R");
    //TTree *raw_tree = (TTree*)file->Get("unpacker_data/hgcroc");
    
    printf("Hello World!\n");

    // loop
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    printf(">>> check nentries = %lld\n", nentries);

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       // if (Cut(ientry) < 0) continue;
    }

    // masked channels (TBI)
    
    // pedestal subtraction
    
    tr_evtNo += 1;
} // end of analyze

// others {{{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
// if the SetupData is always needed
auto setup = iSetup.getData(setupToken_);
// if need the ESHandle to check if the SetupData was there or not
auto pSetup = iSetup.getHandle(setupToken_);
#endif


// ------------ method called once each job just before starting event loop  ------------
void Calibration::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void Calibration::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Calibration::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(Calibration);
