#!/usr/bin/env python
import os
import subprocess
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-a', help = 'add packages', action='store_true')
parser.add_argument('-r', help = 'fetch cms release', action='store_true')
parser.add_argument('-t', help = 'print ED tutorial', action='store_true')
args = parser.parse_args()

cmssw = "CMSSW_12_1_X_2021-10-16-1100"
cmssw = "CMSSW_12_2_X_2021-10-29-2300"

arch = "slc7_amd64_gcc700"
arch = "slc7_amd64_gcc900"

def fetch_cms_release(): #{{{
    print "scram list -a"
    print "scram list CMSSW"
    print "export SCRAM_ARCH=%s" % arch
    print "cmsrel %s" % cmssw
    print "cd %s/src" % cmssw
    print "cmsenv"
    print "git cms-init"
    print "cd $CMSSW_BASE/src"

#}}}

def get_cmssw_path(): #{{{
    log = "tmp_cmssw_path.txt"
    subprocess.call("echo $CMSSW_BASE > %s" % log, shell=True)

    with open(log, 'r') as f:
        for line in f.readlines():
            return line.strip()
#}}}
def add_pkg(): #{{{
    path = get_cmssw_path()
    os.chdir(path + '/src')
    print ">>>", os.getcwd()

    subprocess.call('pwd', shell=True)
    subprocess.call("git cms-addpkg Geometry/HGCalCommonData", shell=True)
    subprocess.call("git cms-addpkg Geometry/HGCalGeometry", shell=True)
    subprocess.call("git cms-addpkg Geometry/HGCalSimData", shell=True)
    subprocess.call("git cms-addpkg DataFormats/ForwardDetId", shell=True)
    subprocess.call("git cms-addpkg Geometry/CMSCommonData", shell=True)
    subprocess.call("git cms-addpkg Validation/HGCalValidation", shell=True) # newly added, to be checked
#}}}
def print_common_xml(): #{{{
    path = get_cmssw_path()
    os.chdir(path + '/src')
    print ">>>", os.getcwd()
    print ">>> print out xml files specific to HGCal"

    # 9 xml files in Geometry/HGCalCommonData/data
    subprocess.call("ls Geometry/HGCalCommonData/data/hgcal/v15/hgcal.xml", shell=True)
    subprocess.call("ls Geometry/HGCalCommonData/data/hgcalCons/v14/hgcalCons.xml", shell=True)
    subprocess.call("ls Geometry/HGCalCommonData/data/hgcalConsData/v15/hgcalConsData.xml", shell=True)
    subprocess.call("ls Geometry/HGCalCommonData/data/hgcalEE/v15/hgcalEE.xml", shell=True)
    subprocess.call("ls Geometry/HGCalCommonData/data/hgcalHEmix/v15/hgcalHEmix.xml", shell=True)
    subprocess.call("ls Geometry/HGCalCommonData/data/hgcalHEsil/v15/hgcalHEsil.xml", shell=True)
    subprocess.call("ls Geometry/HGCalCommonData/data/hgcalcell/v15/hgcalcell.xml", shell=True)
    subprocess.call("ls Geometry/HGCalCommonData/data/hgcalwafer/v15/hgcalwafer.xml", shell=True)
    subprocess.call("ls Geometry/HGCalCommonData/data/hgcalMaterial/v2/hgcalMaterial.xml", shell=True)

    # 3 xml files in Geometry/HGCalCommonData/data/TB*
    # TB160, TB161, TB170, TB180, TB181
    subprocess.call("ls Geometry/HGCalCommonData/data/TB160/hgcalsense.xml", shell=True)
    subprocess.call("ls Geometry/HGCalCommonData/data/TB160/hgcProdCuts.xml", shell=True)
    subprocess.call("ls Geometry/HGCalCommonData/data/TB160/cms.xml", shell=True)

    # 2 xml files in Geometry/CMSCommonData/data
    subprocess.call("ls Geometry/CMSCommonData/data/caloBase/2026/v5/caloBase.xml", shell=True)
    subprocess.call("ls Geometry/CMSCommonData/data/materials/2021/v2/materials.xml", shell=True)
#}}}

def print_tutorial():
    print ""
    print "--- ED Producer ---"
    print "cd $CMSSW_BASE/src"
    print "cmsenv"
    print "mkdir ProdTutorial"
    print "cd ProdTutorial"
    print "mkedprod ProducerTest"
    print ""
    print "--- ED Analyzer ---"
    print "git clone git@github.com:ZhengGang85129/CLUE_FLATTUPLE.git"
    print ""

if __name__ == "__main__":

    if args.r:
        fetch_cms_release()

    if args.t:
        print_tutorial()

    if args.a:
        add_pkg()
        print_common_xml()
        #grep DDHGCalSiliconModule Geometry/HGCalCommonData/data/*/*/*xml
        #grep DDHGCalMixedLayer Geometry/HGCalCommonData/data/*/*/*xml

    print ">>> finished!"
