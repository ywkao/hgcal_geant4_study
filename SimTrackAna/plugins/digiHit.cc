#include "digiHit.h"

// ------------ method called for each event  ------------ //
void DigiSim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // init {{{
    using namespace edm;

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
    double total_corrected_energy_set0 = 0.;
    double total_corrected_energy_set1 = 0.;
    double total_corrected_energy_set2 = 0.;
    std::vector<double> total_energy_adc_total  ;
    std::vector<double> total_energy_mip_total  ;
    std::vector<double> total_energy_mip_coarse ;
    std::vector<double> total_energy_mip_fine   ;
    std::vector<double> total_energy_sim_total  ;
    std::vector<int>    num_digis_total         ;
    std::vector<int>    num_simhits_total       ;
    std::vector<int>    num_simhits_coarse      ;
    std::vector<int>    num_simhits_fine        ;

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
    // Primary vertex {{{
    //----------------------------------------------------------------------------------------------------
    // Not used becasue it is checked that (px, py, pz) of gen-particle can be used for linear trajectory
    edm::Handle<edm::HepMCProduct> mcHandle;
    iEvent.getByToken(mc_, mcHandle);
    HepMC::GenVertex *primaryVertex = *(mcHandle)->GetEvent()->vertices_begin();

    if(false) {
        // PV (x, y, z) in unit of mm
        tb::print_debug_info("pv.x()" , primaryVertex->position().x()       );
        tb::print_debug_info("pv.y()" , primaryVertex->position().y()       );
        tb::print_debug_info("pv.z()" , primaryVertex->position().z()       );
        tb::print_debug_info("pv.t()" , primaryVertex->position().t(), true );
    }
    //----------------------------------------------------------------------------------------------------}}}
    // Gen particles {{{
    //----------------------------------------------------------------------------------------------------
    std::vector<reco::GenParticle> truth_positron;
    Handle<std::vector<reco::GenParticle> > genParticleHandle;
    iEvent.getByToken(genParticles_, genParticleHandle);
    for(size_t i = 0; i < genParticleHandle->size(); ++i )  {    
        const reco::GenParticle &p = (*genParticleHandle)[i];
        if(fabs(p.pdgId())!=11) continue;
        if(!p.isPromptFinalState()) continue;    
        if(fabs(p.eta())<1.5 || fabs(p.eta())>2.9) continue;
        
        truth_positron.push_back(p);

        if(false) {
            tb::print_debug_info("p.pdgId()" , p.pdgId()     );
            tb::print_debug_info("p.pt()"    , p.pt()        );
            tb::print_debug_info("p.energy()", p.energy()    );
            tb::print_debug_info("p.px()"    , p.px()        );
            tb::print_debug_info("p.py()"    , p.py()        );
            tb::print_debug_info("p.px()"    , p.pz()        );
            tb::print_debug_info("p.eta()"   , p.eta(), true );
        }
    }

    // Efficiency is 100% for the cases of R90To130 positron beam
    if(truth_positron.size()!=1) return;

    // Record gen-info
    double gen_eta = truth_positron[0].eta();
    double gen_phi = truth_positron[0].phi();
    tr_e = truth_positron[0].energy(); // GeV
    tr_eta = gen_eta;
    tr_phi = gen_phi;
    tr_positron->Fill();

    // CloseByParticle gun: PV (x, y, z) is parallel with gen-particle momentum (px, py, pz)
    if(false) {
        tb::print_debug_info("p.px()"    , truth_positron[0].px()        );
        tb::print_debug_info("p.py()"    , truth_positron[0].py()        );
        tb::print_debug_info("p.pz()"    , truth_positron[0].pz(), true  );

        tb::print_debug_info("r.px()"    , primaryVertex->position().x()/truth_positron[0].px()        );
        tb::print_debug_info("r.py()"    , primaryVertex->position().y()/truth_positron[0].py()        );
        tb::print_debug_info("r.pz()"    , primaryVertex->position().z()/truth_positron[0].pz(), true  );
    }

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
        bool isSilicon = rhtools_.isSilicon(detId);
        bool isScintillator = rhtools_.isScintillator(detId);
        int wafer_type = HGCSiliconDetId(detId).type();

        bool store_hit_position_info = false; //{{{
        if(store_hit_position_info) {
            // Warn: Might need to set proper sub detector configuration.
            // Currently the code is for CEE analayis
            GlobalPoint gp = rhtools_.getPosition(detId);
            float x = gp.x();
            float y = gp.y();
            float z = gp.z();
            float Rxy = sqrt(pow(x, 2) + pow(y, 2));

            bool is_Silicon_w120 = false;
            bool is_Silicon_w200 = false;
            bool is_Silicon_w300 = false;
            bool is_Scintillator = isScintillator;

            if(isSilicon) {
                if(wafer_type==0) is_Silicon_w120 = true;
                if(wafer_type==1) is_Silicon_w200 = true;
                if(wafer_type==2) is_Silicon_w300 = true;
            }

            //tr_hits->Fill(Rxy, z, is_Silicon_w120, is_Silicon_w200, is_Silicon_w300, is_Scintillator);

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
        } //}}}

        if(isSilicon){
            HGCSiliconDetId id(itHit->id());
            double energy = itHit->energy()*1.e3; // from GeV to MeV

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

            uint32_t id_ = itHit->id();
            energysum esum;
            hitsinfo hinfo;

            if (map_Simhits.count(id_) == 0) {
                hinfo.hitid = nofSiHits;
                hinfo.u_cor = rhtools_.getCell(detId).first ;
                hinfo.v_cor = rhtools_.getCell(detId).second ;
                hinfo.type  = id.type(); // wafer_type
                hinfo.layer = rhtools_.getLayerWithOffset(detId);
                hinfo.is_Silicon_w120 = (wafer_type==0);
                hinfo.is_Silicon_w200 = (wafer_type==1);
                hinfo.is_Silicon_w300 = (wafer_type==2);
                hinfo.is_Scintillator = isScintillator;
            } else {
                hinfo = map_Simhits[id_].first;
                esum = map_Simhits[id_].second;
            }	
            esum.etotal += energy;
            esum.eTime[0] = energy;
            map_Simhits[id_] = std::pair<hitsinfo, energysum>(hinfo, esum);

            bool debug = false;
            if(debug) {
                tb::print_debug_info("hinfo.hitid"           , hinfo.hitid           );
                tb::print_debug_info("hinfo.layer"           , hinfo.layer           );
                tb::print_debug_info("hinfo.type (wafer)"    , hinfo.type            );
                tb::print_debug_info("hinfo.u_cor"           , hinfo.u_cor           );
                tb::print_debug_info("hinfo.v_cor"           , hinfo.v_cor           );
                tb::print_debug_info("hinfo.is_Silicon_w120" , hinfo.is_Silicon_w120 );
                tb::print_debug_info("hinfo.is_Silicon_w200" , hinfo.is_Silicon_w200 );
                tb::print_debug_info("hinfo.is_Silicon_w300" , hinfo.is_Silicon_w300 );
                tb::print_debug_info("hinfo.is_Scintillator" , hinfo.is_Scintillator );
                tb::print_debug_info("esum.etotal"           , esum.etotal           );
                tb::print_debug_info("esum.eTime[0]"         , esum.eTime[0], true   );
                counter += 1;
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
    //Double_t max_energy=0; 
    //for (itr_sim = map_Simhits.begin(); itr_sim != map_Simhits.end(); ++itr_sim) {
    //    energysum esum = (*itr_sim).second.second;
    //    if(max_energy<esum.eTime[0]) max_energy=esum.eTime[0];
    //}

    for (itr_sim = map_Simhits.begin(); itr_sim != map_Simhits.end(); ++itr_sim) {
        //if(counter>9) break;
        uint32_t id_simhit = (*itr_sim).first;
        energysum esum = (*itr_sim).second.second;
        hitsinfo hinfo = (*itr_sim).second.first;
        if(esum.etotal>0) // && esum.eTime[0]==max_energy)
        {
            for (itr_mydigi = my_map_digihits.begin(); itr_mydigi != my_map_digihits.end(); ++itr_mydigi) {
                uint32_t  id_digihit = (*itr_mydigi).first;
                digisinfo dinfo      = (*itr_mydigi).second.dinfo;
                uint32_t  adc        = (*itr_mydigi).second.ainfo.adc;
                double    eta        = dinfo.eta;
                double    phi        = dinfo.phi;
                double    amplitude  = (*itr_mydigi).second.amplitude;
                double    energy     = esum.eTime[0]; // MeV
                int       idx        = dinfo.layer-1;
                int       cellType   = (*itr_sim).second.first.type;
                bool      is_coarse  = cellType==HGCSiliconDetId::HGCalCoarseThin || cellType==HGCSiliconDetId::HGCalCoarseThick;
                bool      is_fine    = cellType==HGCSiliconDetId::HGCalFine;

                bool selection_on_mips = amplitude > 0.5;
                if(!selection_on_mips) continue;

                //bool selectoin_on_eta = eta > 1.8;
                //bool selectoin_on_eta = eta > 1.6 && eta < 2.0;
                //if(!selectoin_on_eta) continue;

                if(id_simhit==id_digihit){
                    double dEdx_weights = get_additional_correction(idx+1); // layer = idx+1
                    amplitude = amplitude * dEdx_weights;

                    //reset_tree_variables();
                    //tr_evtNo = 0;
                    tr_layerNo = idx+1;
                    tr_x = dinfo.x_pos;
                    tr_y = dinfo.y_pos;
                    tr_z = dinfo.z_pos;
                    tr_e = energy; // MeV
                    tr_r = sqrt(pow(tr_x, 2) + pow(tr_y, 2));
                    tr_eta = eta;
                    tr_phi = phi;
                    tr_is_Silicon_w120 = hinfo.is_Silicon_w120;
                    tr_is_Silicon_w200 = hinfo.is_Silicon_w200;
                    tr_is_Silicon_w300 = hinfo.is_Silicon_w300;
                    tr_is_Scintillator = hinfo.is_Scintillator;

                    tr_d = get_distance_from_expected_hit(tr_x, tr_y, tr_z, gen_eta, gen_phi);

                    if(tr_d<=1.3)      tr_signal_region = 1;
                    else if(tr_d<=2.6) tr_signal_region = 2;
                    else if(tr_d<=5.3) tr_signal_region = 3;
                    else               tr_signal_region = -1;

                    tr_hits->Fill();

                    if(false) {
                        counter += 1;
                        tb::print_debug_info("tr_d" , tr_d    );
                        tb::print_debug_info("tr_signal_region" , tr_signal_region, true );
                        continue;
                    }

                    if(false) {
                        counter += 1;
                        tb::print_debug_info("Id_digi"   , id_digihit   );
                        tb::print_debug_info("layer"     , dinfo.layer  );
                        tb::print_debug_info("wafer type", dinfo.type   );
                        tb::print_debug_info("cell type" , cellType     );
                        tb::print_debug_info("adc_"      , adc          );
                        tb::print_debug_info("amplitude" , amplitude    );
                        tb::print_debug_info("energy"    , energy, true );
                        continue;
                    }
                    hEta->Fill(eta);
                    hPhi->Fill(phi);

                    //--------------------------------------------------
                    // info for a set of layers
                    //--------------------------------------------------
                    total_energy_mip_set0 += amplitude;
                    total_energy_sim_set0 += energy;

                    double corrected_energy = convert_amplitude_to_total_energy_pedro(0, amplitude);
                    total_corrected_energy_set0 += corrected_energy;

                    //E_set1 = E1+E3+E5+E7+E9+E11+E13+E15+E17+...E25
                    //E_set2 = E1+E3+E5+E8+E10+E12+E14+E15+E17+...E25
                    bool is_in_set1 = is_this_in_set1(idx+1); // layer = idx+1
                    bool is_in_set2 = is_this_in_set2(idx+1); // layer = idx+1
    
                    if(is_in_set1) {
                        total_energy_mip_set1 += amplitude;
                        total_energy_sim_set1 += energy;
                    }

                    if(is_in_set2) {
                        total_energy_mip_set2 += amplitude;
                        total_energy_sim_set2 += energy;
                    }

                    if(idx%2==0) {
                        total_energy_mip_odd += amplitude;
                        total_energy_sim_odd += energy;
                    } else {
                        total_energy_mip_even += amplitude;
                        total_energy_sim_even += energy;
                    }

                    //--------------------------------------------------
                    // info for each layer
                    //--------------------------------------------------
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

    total_corrected_energy_set0 = convert_amplitude_to_total_energy_pedro(0, total_energy_mip_set0);
    total_corrected_energy_set1 = convert_amplitude_to_total_energy_pedro(1, total_energy_mip_set1);
    total_corrected_energy_set2 = convert_amplitude_to_total_energy_pedro(2, total_energy_mip_set2);

    total_ENE_set0 -> Fill(total_corrected_energy_set0);
    total_ENE_set1 -> Fill(total_corrected_energy_set1);
    total_ENE_set2 -> Fill(total_corrected_energy_set2);

    //tb::print_debug_info("total_energy_sim_set0", total_energy_sim_set0);
    //tb::print_debug_info("total_energy_mip_set1", total_energy_mip_set1);
    //tb::print_debug_info("total_corrected_energy_set1", total_corrected_energy_set1, true);
    //}}}
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
