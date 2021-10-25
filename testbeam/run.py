#!/usr/bin/env python
import subprocess
path = "/Users/ywkao/Desktop/Geant4/test_beam"

def run():
    files = [
        "ntuple_496.root",
        #"ntuple_1240.root",
        #"ntuple_sim_config22_pdgID11_beamMomentum120_listFTFP_BERT_EMN.root",
    ]
    
    for f in files:
        subprocess.call("./examine %s/%s" % (path, f), shell=True)
    
    print(">>> finished!")

if __name__ == "__main__":
    subprocess.call("make", shell=True)
    run()
