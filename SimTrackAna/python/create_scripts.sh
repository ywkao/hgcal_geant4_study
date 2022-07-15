#!/bin/bash

function create()
{
    echo ">>> tag: $1"
    target="run_digiHit_cfg_${1}.py"
    cp run_digiHit_cfg_ProdCut_electron_1mm.py $target
    echo $1 >> $target
    vim $target

    return 

    echo ">>> tag: $1"
    target="test_digiHit_cfg_${1}.py"
    cp digiHit_cfg.py $target
    vim $target
}

function copy()
{
    cp $1 $2
    vim $2
}
#++++++++++++++++++++++++++++++
# main
#++++++++++++++++++++++++++++++

#create "D86_R80To150_E60"
#create "D86_R80To150_E175"
#create "D86_R80To150_E225"
#create "D86_R80To150_E20"
#create "D86_R80To150_E100"
#create "D86_R80To150_E300"

#create D86_R80To100_E100_airPCB_turnOffComptionScattering_ProdCut_electron_10mm
#create D86_R80To100_E100_airPCB_turnOffComptionScattering_ProdCut_electron_5mm
#create D86_R80To100_E100_airPCB_turnOffComptionScattering_ProdCut_electron_100mm
#create D86_R80To100_E100_airPCB_turnOffComptionScattering

#create D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_1000mm
#create D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_10mm
#create D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_1mm
#create D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_5mm

#create D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_50mm
#create D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_100mm

#copy run_digiHit_cfg_ProdCut_egamma_1mm.py run_digiHit_cfg_PCB.py

exit

copy run_digiHit_cfg_ProdCut_egamma_1mm.py run_digiHit_cfg_ProdCut_egamma_1000mm.py
copy run_digiHit_cfg_ProdCut_electron_1mm.py run_digiHit_cfg_ProdCut_electron_1000mm.py
copy run_digiHit_cfg_ProdCut_photon_1mm.py run_digiHit_cfg_ProdCut_photon_1000mm.py

exit

create "D86_R35To60_E20"
create "D86_R35To60_E100"
create "D86_R35To60_E300"
create "D86_R80To100_E20"
create "D86_R80To100_E100"
create "D86_R80To100_E300"
