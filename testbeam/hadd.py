#!/usr/bin/env python
import subprocess
import mylist as ml
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-e", help = "exe", action="store_true")
args = parser.parse_args()

def create_command(output, rootfiles):
    mylist = ""
    for root in rootfiles:
        produced = "./output/" + root.replace("ntuple", "output")
        mylist += produced + " "
        #subprocess.call("ls %s" % produced, shell=True)
    
    command = "hadd -f %s %s" % (output, mylist)
    if args.e:
        subprocess.call(command, shell=True)
    else:
        print command

if __name__ == "__main__":
    create_command( "output/my_test_beam_data_e300.root", ml.test_beam_data_e300 )
    create_command( "output/my_test_beam_data_e100.root", ml.test_beam_data_e100 )
    create_command( "output/my_test_beam_data_e20.root" , ml.test_beam_data_e20  )
