#!/bin/bash
#To be run on remote machine
#Take input arguments as an array
myArray=( "$@" )
#Array: Size=$#, an element=$1, all element = $@

printf "Start Running Histogramming at ";/bin/date
printf "Worker node hostname ";/bin/hostname
CMSVER=CMSSW_12_2_X_2021-12-10-2300

if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then 
    echo "Running Interactively" ; 
else
    echo "Running In Batch"
    echo ${_CONDOR_SCRATCH_DIR}
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    SCRAM_ARCH=slc7_amd64_gcc900
    scramv1 project CMSSW $CMSVER
    cd $CMSVER/src
    eval `scramv1 runtime -sh`
    #git cms-addpkg Configuration/Generator
    #cd ../..
    cp ../../generator.tar.gz .
fi

pwd
ls -la
tar --strip-components=0 -zxvf generator.tar.gz
cp ReadSimResult/GenConfig/CloseByParticle_Photon_ERZRanges_cfi.py Configuration/Generator/python/CloseByParticle_Photon_ERZRanges_cfi.py

if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then 
    echo "Running Interactively" ; 
else
    pwd
    ls -la
    scram b -j 4
fi
#Run for Base, Signal region

echo "All arguements: "$@
echo "Number of arguements: "$#
geom=$1
index=$2

cmsDriver.py CloseByParticle_Photon_ERZRanges_cfi -s GEN,SIM -n 10000 --conditions auto:phase2_realistic_T21 --beamspot HGCALCloseBy --datatier GEN-SIM --eventcontent FEVTDEBUG --geometry $geom --era Phase2C11I13M9 --relval 9000,100 --fileout file:step1_${index}.root  --customise_commands process.RandomNumberGeneratorService.generator.initialSeed="cms.untracked.uint32($RANDOM)" --no_exec  --nThreads 4 > step1_${index}_CloseByParticleGun.log  2>&1

grep -n "process.RandomNumberGeneratorService.generator.initialSeed" CloseByParticle_Photon_ERZRanges_cfi.py

cmsRun CloseByParticle_Photon_ERZRanges_cfi_GEN_SIM.py

cmsDriver.py step2  -s DIGI:pdigi_valid,L1TrackTrigger,L1,DIGI2RAW,HLT:@fake2 --conditions auto:phase2_realistic_T21 --datatier GEN-SIM-DIGI-RAW -n 10000 --eventcontent FEVTDEBUGHLT --geometry Extended2026D86 --era Phase2C11I13M9  --filein  file:step1_${index}.root  --fileout file:step2_${index}.root  --nThreads 4 > step2_${index}_CloseByParticleGun.log  2>&1


ls -ltr

pwd

printf "Simulation completed at ";/bin/date
#---------------------------------------------
#Copy the ouput root files
#---------------------------------------------
condorOutDir1=/eos/user/m/mikumar/HGCAl_Validation/$geom
#condorOutDir=/cms/store/user/idas/SimOut/DeltaPt/$geom
if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then
    echo "Running Interactively" ;
else
    #xrdcp -f ${sample}_tree_*.root root://se01.indiacms.res.in:1094/${condorOutDir}/${year} 
    xrdcp -f CloseByParticle_Photon_ERZRanges_cfi_GEN_SIM.py root://eosuser.cern.ch/${condorOutDir1}
    xrdcp -f step1_${index}.root root://eosuser.cern.ch/${condorOutDir1}
    xrdcp -f step1_${index}_CloseByParticleGun.log root://eosuser.cern.ch/${condorOutDir1}
    xrdcp -f step2_${index}.root root://eosuser.cern.ch/${condorOutDir1}
    xrdcp -f step2_${index}_CloseByParticleGun.log root://eosuser.cern.ch/${condorOutDir1}

    echo "Cleanup"
    cd ../../
    rm -rf $CMSVER
    rm *.root
fi
printf "Done ";/bin/date
