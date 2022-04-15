#!/usr/bin/env python
import subprocess
import glob
import argparse
import parallel_utils
parser = argparse.ArgumentParser()
parser.add_argument("-e", help = "Enable to execute", action="store_true")
args = parser.parse_args()

command_list = []
def regester(command):
    global command_list
    command_list.append(command)

def run():
    global command_list
    if args.e:
        parallel_utils.submit_jobs(command_list, 10)
    else:
        for command in command_list: print command

if __name__ == "__main__":

    subprocess.call("mkdir -p tmp", shell=True)

    for f in glob.glob("python/run_digiHit_cfg_*"):
    #for f in glob.glob("python/test_digiHit_cfg*"):
    #for f in glob.glob("python/test_digiHit_cfg*muon*"):
        tag = f.split('_cfg_')[1].split('.')[0]
        command = "time cmsRun %s 2>&1 | tee tmp/log_%s.txt" % (f, tag)
        regester(command)

    run()


#python/test_digiHit_cfg_D86_R35To60_E100.py
#python/test_digiHit_cfg_D86_R35To60_E20.py
#python/test_digiHit_cfg_D86_R35To60_E300.py
#python/test_digiHit_cfg_D86_R80To100_E100.py
#python/test_digiHit_cfg_D86_R80To100_E20.py
#python/test_digiHit_cfg_D86_R80To100_E300.py

