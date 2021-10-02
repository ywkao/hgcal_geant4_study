#!/usr/bin/env python
import os
import subprocess

cmssw = "CMSSW_12_1_X_2021-09-19-2300"

def add_pkg(): #{{{
    subprocess.call('pwd', shell=True)
    subprocess.call("git cms-addpkg Geometry/HGCalGeometry", shell=True)
    subprocess.call("git cms-addpkg Geometry/HGCalCommonData", shell=True)
    subprocess.call("git cms-addpkg Geometry/HGCalSimData", shell=True)
    subprocess.call("git cms-addpkg Geometry/CMSCommonData", shell=True)
    subprocess.call("git cms-addpkg DataFormats/ForwardDetId", shell=True)
#}}}
def get_cmssw_path(): #{{{
    log = "tmp_cmssw_path.txt"
    subprocess.call("echo $CMSSW_BASE > %s" % log, shell=True)

    with open(log, 'r') as f:
        for line in f.readlines():
            return line.strip()
#}}}
def fetch_cms_release(): #{{{
    subprocess.call("/cvmfs/cms.cern.ch/common/scramv1 %s" % cmssw, shell=True) # cmsrel

    #os.chdir("%s/src" % cmssw)
    #subprocess.call("eval `scramv1 runtime -sh`", shell=True) # cmsenv
    #subprocess.call("git cms-init", shell=True)

    #path = get_cmssw_path()
    #os.chdir(path + '/src')

    #add_pkg()
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

if __name__ == "__main__":
    #subprocess.call("scram list -a | grep CMSSW_12_1", shell=True)

    #fetch_cms_release()
    add_pkg()
    #print_common_xml()

    #grep DDHGCalSiliconModule Geometry/HGCalCommonData/data/*/*/*xml
    #grep DDHGCalMixedLayer Geometry/HGCalCommonData/data/*/*/*xml
