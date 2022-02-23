#define hits_cxx
#include "hits.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void hits::Loop()
{
    fChain->SetBranchStatus("*",0);  // disable all branches
    fChain->SetBranchStatus("rechit_layer",1);  // activate branchname
    fChain->SetBranchStatus("rechit_energy",1);  // activate branchname
    fChain->SetBranchStatus("rechit_energy_noHG",1);  // activate branchname
    fChain->SetBranchStatus("rechit_amplitudeHigh",1);  // activate branchname
    fChain->SetBranchStatus("rechit_amplitudeLow",1);  // activate branchname

    fChain->SetBranchStatus("pdgID",1);  // activate branchname
    fChain->SetBranchStatus("beamEnergy",1);  // activate branchname

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    myEntries = (double) nentries;

    bool debug = false;
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;

        nb = fChain->GetEntry(jentry);   nbytes += nb;

        // check pdgId and beam energy for the root file
        mean_pdgID += pdgID;
        mean_energy += beamEnergy;

        // init counters for 28 layers
        std::vector<double> total_amplitude_high;
        std::vector<double> total_amplitude_low;
        for (int i=0; i<28; ++i) {
            total_amplitude_high.push_back(0.);
            total_amplitude_low .push_back(0.);
        }

        // loop over all hits
        for (int i=0; i<rechit_layer->size(); ++i) {
            unsigned int layer = rechit_layer->at(i);
            h_rechit_layer         -> Fill(rechit_layer         -> at(i));
            h_rechit_energy        -> Fill(rechit_energy        -> at(i));
            h_rechit_energy_noHG   -> Fill(rechit_energy_noHG   -> at(i));
            h_rechit_amplitudeHigh -> Fill(rechit_amplitudeHigh -> at(i));
            h_rechit_amplitudeLow  -> Fill(rechit_amplitudeLow  -> at(i));

            if(layer <= 28) {
                total_amplitude_high[layer-1] += rechit_amplitudeHigh -> at(i);
                total_amplitude_low [layer-1] += rechit_amplitudeLow  -> at(i);
            }

            if(debug) printf("%2d, value = %u\n", i, rechit_layer->at(i));
        }

        // fill total amplitude (MIPs, I think) per layer for an event
        for (int i=0; i<28; ++i) {
            //printf("%.2f, ", total_amplitude_high[i]);
            //printf("%.2f, ", total_amplitude_low[i]);

            h_rechit_amplitudeHigh_layers[i] -> Fill(total_amplitude_high[i]);
            h_rechit_amplitudeLow_layers[i]  -> Fill(total_amplitude_low[i]);
        }

    } // end of event loop


}
