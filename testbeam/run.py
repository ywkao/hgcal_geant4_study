#!/usr/bin/env python
import subprocess
import mylist
import parallel_utils
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-e", help = "Enable to execute", action="store_true")
args = parser.parse_args()

path = "/Users/ywkao/Desktop/Geant4/test_beam"
path = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v17"
#"/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v17_patch2/ntuple_654.root"

command_list = []
def register(command):
    global command_list
    command_list.append(command)

def run():
    global command_list
    if args.e:
        parallel_utils.submit_jobs(command_list, 10)
    else:
        for command in command_list: print command

if __name__ == "__main__":
    subprocess.call("make", shell=True)

    files = [ "ntuple_496.root", "ntuple_sim_config22_pdgID11_beamMomentum120_listFTFP_BERT_EMN.root" ]
    files = mylist.test_beam_data_e300
    files = mylist.test_beam_data
    
    for f in files:
        #subprocess.call("./examine %s/%s %s" % (path, f, f.split('_')[1].split('.')[0]), shell=True)
        register( "./examine %s/%s %s" % (path, f, f.split('_')[1].split('.')[0]) )
        #break
    
    run()

    print(">>> finished!")
