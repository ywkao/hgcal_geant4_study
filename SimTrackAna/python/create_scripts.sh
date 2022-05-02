#!/bin/bash

function create()
{
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

