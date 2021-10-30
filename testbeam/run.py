#!/usr/bin/env python
import subprocess
import mylist
path = "/Users/ywkao/Desktop/Geant4/test_beam"
path = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v17"

def run():
    files = [
        "ntuple_496.root",
        #"ntuple_1240.root",
        #"ntuple_sim_config22_pdgID11_beamMomentum120_listFTFP_BERT_EMN.root",
    ]

    files = mylist.test_beam_data
    files = mylist.test_beam_data_e300
    
    for f in files:
        subprocess.call("./examine %s/%s %s" % (path, f, f.split('_')[1].split('.')[0]), shell=True)
    
    print(">>> finished!")

if __name__ == "__main__":
    subprocess.call("make", shell=True)
    run()
