#!/usr/bin/env python
import subprocess
import glob
import argparse
import toolbox.parallel_utils as pu
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
        pu.submit_jobs(command_list, 10)
    else:
        for command in command_list:
            print command
            #subprocess.call(command, shell=True)

if __name__ == "__main__":

    subprocess.call("mkdir -p tmp", shell=True)

    configs = [
        "python/run_digiHit_cfg_D86_R80To110_E100.py",
        "python/run_digiHit_cfg_D86_R80To120_E100.py",
        "python/run_digiHit_cfg_D86_R80To130_E100.py",
        "python/run_digiHit_cfg_D86_R80To140_E100.py",
    ]

    #for f in configs:
    for f in glob.glob("python/*130*"):
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

