#!/usr/bin/env python
import os, subprocess
import ROOT

import toolbox.plot_utils as pu
import toolbox.MetaData as m
import toolbox.runners as ru

eos = "./eos"

def perform_lognitudinal_only():
    m.type_resolution = "resolution_unclustered"
    m.bias = "bias_unclustered"
    m.json_fit_parameters = "./toolbox/unclustered_fit_parameters.json"
    m.json_fit_parameters_goodness = "./toolbox/unclustered_fit_parameters_goodness.json"

    target_directory = "20240224"
    output_directory = eos + "/" + target_directory
    m.specified_directory = output_directory
    pu.create_directory( m.specified_directory )

    m.tags = ["E300", "E100", "E20"]
    fit_constraints = m.fit_constraints_v1p1
    ru.run_manager( m.input_files["R90To130_v1p1"], m.tags, fit_constraints )

    #----------------------------------------------------------------------------------------------------

    target_directory = "20240224_pcb"
    output_directory = eos + "/" + target_directory
    m.specified_directory = output_directory
    pu.create_directory( output_directory )
    m.tags = ["nominal", "Turn_off_Compton", "TOC+airPCB", "airPCB"]
    ru.run_manager( m.input_files["pcbStudy"], m.tags, fit_constraints )

    # tags = ["E300", "E100", "E20"]
    # for tag in tags: label[tag] = tag.split("E")[1] + " GeV"
    # run( m.input_files["R80To100"], eos + "/" + "R80To100_v4p3" )
    # exit()

    # tags = ["nominal", "with_dEdx_weight"]
    # for tag in tags: label[tag] = tag
    # run( m.input_files["dEdxStudy"]    , eos + "/" + "R80To100_study_with_dEdx_weights"    )
    # exit()

    # #tags = ["nominal", "Turn_off_Compton", "Replace_PCB_by_Air", "TurnOffCompton_CutEle5mm", "TurnOffCompton_CutEle10mm"]
    # #tags = ["nominal", "Turn_off_Compton", "TOC+airPCB", "TOC+airPCB+CutEle5mm", "TOC+airPCB+CutEle10mm"]
    # tags = ["nominal", "Turn_off_Compton", "TOC+airPCB", "airPCB"]
    # for tag in tags: label[tag] = tag
    # run( m.input_files["pcbStudy"]    , eos + "/" + "R80To100_PCB_study_with_dEdx_weights_v2"    )
    # exit()

    # tags = ["nominal", "Turn_off_Compton", "TurnOffCompton_CutEle1mm", "TurnOffCompton_CutEle10mm", "TurnOffCompton_CutEle100mm", "TurnOffCompton_CutEle100mm"]
    # for tag in tags: label[tag] = tag
    # run( m.input_files["turnOffCompton_ProdCutElectron"]    , eos + "/" + "R80To100_turnOffCompton_ProdCutElectron_v3"    )
    # exit()

    # tags = ["nominal", "Turn_off_Compton", "TurnOffCompton_CutEle1mm", "TurnOffCompton_CutEle5mm", "TurnOffCompton_CutEle10mm", "TurnOffCompton_CutEle50mm", "TurnOffCompton_CutEle100mm"]
    # for tag in tags: label[tag] = tag
    # run( m.input_files["turnOffCompton_ProdCutElectron"]    , eos + "/" + "R80To100_turnOffCompton_ProdCutElectron_v2"    )
    # exit()

    # # extra
    # tags = ["nominal", "Turn_off_Compton", "Turn_off_Conversion", "Turn_off_both"]
    # for tag in tags: label[tag] = tag
    # run( m.input_files["extraStudy"]    , eos + "/" + "R80To100_extra_study"    )
    # exit()

    # # compton
    # tags = ["nominal", "Turn_off_Compton", "Turn_off_both", "ProdCut_electron"]
    # for tag in tags: label[tag] = tag
    # run( m.input_files["turnOffCompton"]    , eos + "/" + "R80To100_TurnOffCompton_comparison"    )
    # exit()

    # # prodCut study
    # tags = ["nominal", "ProdCut_electron", "ProdCut_photon", "ProdCut_egamma"]
    # for tag in tags: label[tag] = tag
    # run( m.input_files["ProdCut1mm"]    , eos + "/" + "R80To100_ProdCut1mm"    )
    # run( m.input_files["ProdCut100mm"]  , eos + "/" + "R80To100_ProdCut100mm"  )
    # run( m.input_files["ProdCut1000mm"] , eos + "/" + "R80To100_ProdCut1000mm" )

def perform_unclustered_study():
    m.type_resolution = "resolution_unclustered"
    m.bias = "bias_unclustered"
    m.json_fit_parameters = "./toolbox/unclustered_fit_parameters.json"
    m.json_fit_parameters_goodness = "./toolbox/unclustered_fit_parameters_goodness.json"

    target_directory = "R90To130_v1p6"
    target_directory = "R90To130_unclustered"
    target_directory = "20240224"
    output_directory = eos + "/" + target_directory
    m.specified_directory = output_directory
    pu.create_directory( m.specified_directory )

    #m.tags = ["E300", "E100", "E20"]
    #fit_constraints = m.fit_constraints_v1p1
    #ru.run_manager( m.input_files["R90To130_v1p1"], m.tags, fit_constraints )
    #return

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
    m.bias = "bias_clustered"
    m.json_fit_parameters = "./toolbox/clustered_fit_parameters.json"
    m.json_fit_parameters_goodness = "./toolbox/clustered_fit_parameters_goodness.json"

    target_directory = "R90To130_v2p6"
    target_directory = "20221005"
    target_directory = "R90To130_clustered"
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
    perform_lognitudinal_only()
    #perform_unclustered_study()
    #perform_clustered_study()

