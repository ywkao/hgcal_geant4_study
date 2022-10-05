#!/usr/bin/env python
import os, subprocess
import ROOT

import toolbox.plot_utils as pu
import toolbox.MetaData as m
import toolbox.runners as ru

eos = "./eos"

def perform_unclustered_study():
    m.type_resolution = "resolution_unclustered"
    m.json_fit_parameters = "./toolbox/unclustered_fit_parameters.json"
    m.json_fit_parameters_goodness = "./toolbox/unclustered_fit_parameters_goodness.json"

    target_directory = "R90To130_v1p6"
    target_directory = "R90To130_unclustered"
    output_directory = eos + "/" + target_directory
    m.specified_directory = output_directory
    pu.create_directory( m.specified_directory )

    if run_full_commands:
        tags = ["E300", "E100", "E20"]
        fit_constraints = m.fit_constraints_v1p1
        ru.run_manager( m.input_files["R90To130_v1p1"], tags, fit_constraints )

        tags = ["E225", "E175", "E60"]
        fit_constraints = m.fit_constraints_v1p2
        ru.run_manager( m.input_files["R90To130_v1p2"], tags, fit_constraints )
        ru.run_register_fit_parameters()

    ru.run_fitters_and_summary()

def perform_clustered_study():
    m.type_resolution = "resolution_clustered"
    m.json_fit_parameters = "./toolbox/clustered_fit_parameters.json"
    m.json_fit_parameters_goodness = "./toolbox/clustered_fit_parameters_goodness.json"

    target_directory = "R90To130_v2p6"
    target_directory = "R90To130_clustered"
    target_directory = "20221005"
    output_directory = eos + "/" + target_directory
    m.specified_directory = output_directory
    pu.create_directory( m.specified_directory )

    if run_full_commands:
        tags = ["E300", "E100", "E20"]
        fit_constraints = m.fit_constraints_v2p1
        ru.run_manager( m.input_files["R90To130_v2p1"], tags, fit_constraints )

        tags = ["E225", "E175", "E60"]
        fit_constraints = m.fit_constraints_v2p2
        ru.run_manager( m.input_files["R90To130_v2p2"], tags, fit_constraints )
        ru.run_register_fit_parameters()

    ru.run_fitters_and_summary()

#----------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    m.enable_check_odd_even = False
    run_full_commands = False 
    run_full_commands = True
    #perform_unclustered_study()
    perform_clustered_study()

