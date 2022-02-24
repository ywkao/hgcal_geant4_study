#!/bin/bash

function create()
{
    echo ">>> tag: $1"
    target="test_digiHit_cfg_${1}.py"
    cp digiHit_cfg.py $target
    vim $target
}

exit

create "D86_R35To60_E20"
create "D86_R35To60_E100"
create "D86_R35To60_E300"
create "D86_R80To100_E20"
create "D86_R80To100_E100"
create "D86_R80To100_E300"
