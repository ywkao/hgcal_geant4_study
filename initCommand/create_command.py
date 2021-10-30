#!/usr/bin/env python
import os
import subprocess
import argparse
import my_metadata_sets
cwd = os.getcwd()
parser = argparse.ArgumentParser()
parser.add_argument("-e", help = "enable execution", action="store_true")
args = parser.parse_args()

# exmples {{{
## quick start
example_step1 = " cmsDriver.py CloseByParticle_Photon_ERZRanges_cfi -s GEN,SIM -n 10 --conditions auto:phase2_realistic_T15 --beamspot HGCALCloseBy --datatier GEN-SIM --eventcontent FEVTDEBUG --geometry Extended2026D49 --era Phase2C9 --relval 9000,100 --fileout file:step1.root  > step1_CloseByParticleGun+2026D49+CloseByParticle_Photon_ERZRanges_GenSimHLBeamSpotHGCALCloseBy+DigiTrigger+RecoGlobal+HARVESTGlobal.log  2>&1"
example_step2 = " cmsDriver.py step2  -s DIGI:pdigi_valid,L1TrackTrigger,L1,DIGI2RAW,HLT:@fake2 --conditions auto:phase2_realistic_T15 --datatier GEN-SIM-DIGI-RAW -n 10 --eventcontent FEVTDEBUGHLT --geometry Extended2026D49 --era Phase2C9 --filein  file:step1.root  --fileout file:step2.root  > step2_CloseByParticleGun+2026D49+CloseByParticle_Photon_ERZRanges_GenSimHLBeamSpotHGCALCloseBy+DigiTrigger+RecoGlobal+HARVESTGlobal.log  2>&1"
example_step3 = " cmsDriver.py step3  -s RAW2DIGI,L1Reco,RECO,RECOSIM,PAT,VALIDATION:@phase2Validation+@miniAODValidation,DQM:@phase2+@miniAODDQM --conditions auto:phase2_realistic_T15 --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO -n 10 --eventcontent FEVTDEBUGHLT,MINIAODSIM,DQM --geometry Extended2026D49 --era Phase2C9 --filein  file:step2.root  --fileout file:step3.root  > step3_CloseByParticleGun+2026D49+CloseByParticle_Photon_ERZRanges_GenSimHLBeamSpotHGCALCloseBy+DigiTrigger+RecoGlobal+HARVESTGlobal.log  2>&1"
example_step4 = " cmsDriver.py step4  -s HARVESTING:@phase2Validation+@phase2+@miniAODValidation+@miniAODDQM --conditions auto:phase2_realistic_T15 --mc  --geometry Extended2026D49 --scenario pp --filetype DQM --era Phase2C9 -n 100  --filein file:step3_inDQM.root --fileout file:step4.root  > step4_CloseByParticleGun+2026D49+CloseByParticle_Photon_ERZRanges_GenSimHLBeamSpotHGCALCloseBy+DigiTrigger+RecoGlobal+HARVESTGlobal.log  2>&1"

## indra
example_init = " cmsDriver.py SingleMuFlatPt2To100_cfi.py -s GEN,SIM --conditions auto:mc --datatier GEN-SIM-RAW --eventcontent RAWSIM -n 100000 --no_exec --python_filename SingleMuFlatPt2To100_cfi_py_GEN_geo_default.py --fileout file:SingleMuFlatPt2To100_cfi_py_GEN_geo_default_Phase2C11_Extended2026D83_higheta.root --era Phase2C11 --geometry Extended2026D83 --nThreads 8"
example_D83 = " cmsDriver.py SingleMuFlatPt2To100_cfi -s GEN,SIM -n 100000 --conditions auto:phase2_realistic_T21 --beamspot HLLHC14TeV --datatier GEN-SIM --eventcontent FEVTDEBUG --geometry Extended2026D83 --era Phase2C11I13M9 --fileout file:SingleMuFlatPt2To100_D83_step1.root --nThreads 8 "
example_D86 = " cmsDriver.py SingleMuFlatPt2To100_cfi -s GEN,SIM -n 100000 --conditions auto:phase2_realistic_T21 --beamspot HLLHC14TeV --datatier GEN-SIM --eventcontent FEVTDEBUG --geometry Extended2026D86 --era Phase2C11I13M9 --fileout file:SingleMuFlatPt2To100_D86_step1.root --nThreads 8 "
#}}}
# templates {{{
template_step1 = ''' cmsDriver.py {PARTICLE_GUN} -s GEN,SIM -n {NEVENTS} \
--conditions {CONDITIONS} \
--beamspot {BEAMSPOT} \
--datatier {DATATIER} \
--eventcontent {EVENTCONTENT} \
--geometry {GEOMETRY} \
--era {ERA} \
--fileout {FILEOUT} \
--nThreads {NTHREADS} \
--no_exec'''

template_step2 = ''' cmsDriver.py step2  -s DIGI:pdigi_valid,L1TrackTrigger,L1,DIGI2RAW,HLT:@fake2 \
--conditions {CONDITIONS} \
--datatier {DATATIER} \
-n {NEVENTS} \
--eventcontent {EVENTCONTENT} \
--geometry {GEOMETRY} \
--era {ERA} \
--filein  {FILEIN}  \
--fileout {FILEOUT} \
--no_exec'''

template_step3 = ''' cmsDriver.py step3  -s RAW2DIGI,L1Reco,RECO,RECOSIM,PAT,VALIDATION:@phase2Validation+@miniAODValidation,DQM:@phase2+@miniAODDQM \
--conditions {CONDITIONS} \
--datatier {DATATIER} \
-n {NEVENTS} \
--eventcontent {EVENTCONTENT} \
--geometry {GEOMETRY} \
--era {ERA} \
--filein  {FILEIN}  \
--fileout {FILEOUT} \
--nThreads {NTHREADS} \
--no_exec'''

template_step4 = ''' cmsDriver.py step4  -s HARVESTING:@phase2Validation+@phase2+@miniAODValidation+@miniAODDQM \
--conditions {CONDITIONS} \
--mc  \
--geometry {GEOMETRY} \
--scenario {SCENARIO} \
--filetype {FILETYPE} \
--era {ERA} \
-n {NEVENTS}  \
--filein {FILEIN} \
--fileout {FILEOUT} \
--nThreads {NTHREADS} \
--no_exec'''
#}}}
python_file_fragments = [ #{{{
    "SingleMuFlatPt2To100_cfi_GEN_SIM.py", # the 1st string will be updated in create_command(par) 
    "step2_DIGI_L1TrackTrigger_L1_DIGI2RAW_HLT.py",
    "step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PAT_VALIDATION_DQM.py",
    "step4_HARVESTING.py",
] #}}}

directory = ""
command_list = []
def register(command): #{{{
    global command_list
    command_list.append(command)
#}}}
def prepare_scripts(): #{{{
    global directory
    subprocess.call("mkdir -p %s" % directory, shell=True)
    os.chdir(directory)

    script = "cmd.sh"
    with open(script, 'w') as fout:
        fout.write('#!/bin/bash\n\n')
        for command in command_list:
            fout.write(command + '\n\n')

    if args.e:
        # create python files automatically
        subprocess.call("chmod +x %s" % script, shell=True)
        subprocess.call('time ./%s' % script, shell=True)

        # prepare scripts for cmsRun
        subprocess.call('sed -i "s/ cmsDriver/# cmsDriver/g" %s' % script, shell=True)
        with open(script, 'a') as fout:
            for py in python_file_fragments:
                log = 'log_' + py.replace('.py', '.txt')
                command = 'time cmsRun ' + py + ' > ' + log + ' 2>&1'
                fout.write( command + '\n' )
                subprocess.call('echo "double check: `ls %s`"' % py, shell=True)

        print ''
        subprocess.call('cat ./%s' % script, shell=True)
        print ''
        print '>>> cd %s' % directory
        print '>>> time ./%s' % script

    else:
        subprocess.call('cat ./%s' % script, shell=True)
        print ''
        print '>>> ./create_command.py -e'

    os.chdir(cwd)
#}}}

def create_command(par):
    ''' Input parameters: my_metadata_sets.py
        Command for each step are created by register()
        Configuration files for cmsRun are prepared by prepare_scripts() '''
    # init {{{
    global command_list, directory, python_file_fragments
    command_list = []
    directory = par["directory"]
    python_file_fragments[0] = par["particle_gun"] + '_GEN_SIM.py'

    particle_gun = par["particle_gun"]
    nevents      = par["nevents"]
    conditions   = par["conditions"]
    beamspot     = par["beamspot"]
    geometry     = par["geometry"]
    era          = par["era"]
    #}}}
    register( template_step1.format( #{{{
                PARTICLE_GUN = particle_gun,
                NEVENTS      = nevents,
                CONDITIONS   = conditions,
                BEAMSPOT     = beamspot,
                DATATIER     = "GEN-SIM",
                EVENTCONTENT = "FEVTDEBUG",
                GEOMETRY     = geometry,
                ERA          = era,
                FILEOUT      = "file:step1.root",
                NTHREADS     = "8"
            )) #}}}
    register( template_step2.format( #{{{
                CONDITIONS   = conditions,
                DATATIER     = "GEN-SIM-DIGI-RAW",
                NEVENTS      = "-1",
                EVENTCONTENT = "FEVTDEBUGHLT",
                GEOMETRY     = geometry,
                ERA          = era,
                FILEIN       = "file:step1.root",
                FILEOUT      = "file:step2.root",
                NTHREADS     = "8"
            )) #}}}
    register( template_step3.format( #{{{
                CONDITIONS   = conditions,
                DATATIER     = "GEN-SIM-RECO,MINIAODSIM,DQMIO",
                NEVENTS      = "-1",
                EVENTCONTENT = "FEVTDEBUGHLT,MINIAODSIM,DQM",
                GEOMETRY     = geometry,
                ERA          = era,
                FILEIN       = "file:step2.root",
                FILEOUT      = "file:step3.root",
                NTHREADS     = "8"
            )) #}}}
    register( template_step4.format( #{{{
                CONDITIONS = conditions,
                GEOMETRY   = geometry,
                SCENARIO   = "pp",
                FILETYPE   = "DQM",
                ERA        = era,
                NEVENTS    = "-1",
                FILEIN     = "file:step3_inDQM.root",
                FILEOUT    = "file:step4.root",
                NTHREADS     = "8"
            )) #}}}
    prepare_scripts()

if __name__ == "__main__":
    create_command( my_metadata_sets.parameters['set1'] )
    #create_command( my_metadata_sets.parameters['set2'] )
    #create_command( my_metadata_sets.parameters['set3'] )
    #create_command( my_metadata_sets.parameters['set4'] )

# relval 9000,100

