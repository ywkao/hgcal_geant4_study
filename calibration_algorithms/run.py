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
    command_list.append(command)

def run():
    if args.e:
        pu.submit_jobs(command_list, 10)
    else:
        for command in command_list: print command

if __name__ == "__main__":

    subprocess.call("mkdir -p tmp", shell=True)

    for f in glob.glob("python/*90*"):
        tag = f.split('_cfg_')[1].split('.')[0]
        command = "time cmsRun %s 2>&1 | tee tmp/log_%s.txt" % (f, tag)
        register(command)

    run()

