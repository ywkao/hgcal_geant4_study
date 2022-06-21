#!/usr/bin/env python
import subprocess
import glob
import argparse
import parallel_utils
parser = argparse.ArgumentParser()
parser.add_argument("-e", help = "Enable to execute", action="store_true")
args = parser.parse_args()

command_list = []
def register(command):
    global command_list
    command_list.append(command)

def run():
    global command_list
    if args.e:
        parallel_utils.submit_jobs(command_list, 10)
    else:
        for command in command_list:
            print command
            #subprocess.call(command, shell=True)

if __name__ == "__main__":

    subprocess.call("mkdir -p tmp", shell=True)

    configs = [
        #"python/run_digiHit_cfg_PCB.py",
        #"python/run_digiHit_cfg_D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_50mm.py",
        #"python/run_digiHit_cfg_D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_100mm.py",
        #"python/run_digiHit_cfg_turnOffCompton.py",

        "python/run_digiHit_cfg_nominal.py",
        "python/run_digiHit_cfg_PCB.py",
        "python/run_digiHit_cfg_turnOffCompton.py",
        "python/run_digiHit_cfg_D86_R80To100_E100_airPCB_turnOffComptionScattering.py",
        #"python/run_digiHit_cfg_D86_R80To100_E100_airPCB_turnOffComptionScattering_ProdCut_electron_100mm.py",
        #"python/run_digiHit_cfg_D86_R80To100_E100_airPCB_turnOffComptionScattering_ProdCut_electron_10mm.py",
        #"python/run_digiHit_cfg_D86_R80To100_E100_airPCB_turnOffComptionScattering_ProdCut_electron_5mm.py",
        #"python/run_digiHit_cfg_D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_1000mm.py",
        #"python/run_digiHit_cfg_D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_100mm.py",
        #"python/run_digiHit_cfg_D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_10mm.py",
        #"python/run_digiHit_cfg_D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_1mm.py",
        #"python/run_digiHit_cfg_D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_50mm.py",
        #"python/run_digiHit_cfg_D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_5mm.py",
        #"python/run_digiHit_cfg_ProdCut_egamma_1000mm.py",
        #"python/run_digiHit_cfg_ProdCut_egamma_100mm.py",
        #"python/run_digiHit_cfg_ProdCut_egamma_1mm.py",
        #"python/run_digiHit_cfg_ProdCut_electron_1000mm.py",
        #"python/run_digiHit_cfg_ProdCut_electron_100mm.py",
        #"python/run_digiHit_cfg_ProdCut_electron_1mm.py",
        #"python/run_digiHit_cfg_ProdCut_photon_1000mm.py",
        #"python/run_digiHit_cfg_ProdCut_photon_100mm.py",
        #"python/run_digiHit_cfg_ProdCut_photon_1mm.py",
    ]

    for f in configs:
    #for f in glob.glob("python/*air*.py"):
    #for f in glob.glob("python/*turnOffComptionScattering_ProdCut_electron*.py"):
    #for f in glob.glob("python/run_digiHit_cfg_turnOffCompton.py"):
    #for f in glob.glob("python/run_digiHit_cfg_*1000*"):
    #for f in glob.glob("python/test_digiHit_cfg*"):
    #for f in glob.glob("python/test_digiHit_cfg*muon*"):
        tag = f.split('_cfg_')[1].split('.')[0]
        command = "time cmsRun %s 2>&1 | tee tmp/log_%s.txt" % (f, tag)
        #command = "sed -n '18p' %s" % f
        register(command)

    run()


#python/test_digiHit_cfg_D86_R35To60_E100.py
#python/test_digiHit_cfg_D86_R35To60_E20.py
#python/test_digiHit_cfg_D86_R35To60_E300.py
#python/test_digiHit_cfg_D86_R80To100_E100.py
#python/test_digiHit_cfg_D86_R80To100_E20.py
#python/test_digiHit_cfg_D86_R80To100_E300.py

