#!/usr/bin/env python
import os, subprocess
import ROOT

import toolbox.plot_utils as pu
import toolbox.MetaData as m
import toolbox.analyzer as an
import toolbox.plotter as pl

eos = "./eos"

# runners {{{

#----------------------------------------------------------------------------------------------------

def run_logitudinal_profile():
    output_directory = pl.specified_directory + "/longitudinal_profiles"
    pu.create_directory(output_directory)

    thickness = ["120mum", "200mum", "300mum", "total"]
    thickness = ["total", "coarse", "fine"] # consider 120, 200, 300 altogether
    thickness = ["total"]
    for t in thickness:
        pl.make_plot( "multiplicity_simhits_%s" % t , output_directory, True  )
        pl.make_plot( "total_MIP_%s"            % t , output_directory, True  )
        pl.make_plot( "total_SIM_%s"            % t , output_directory, True  )

        continue

        pl.make_plot( "total_ADC_%s"            % t , output_directory, True  )
        pl.make_plot( "total_MIP_%s"            % t , output_directory, True  )
        pl.make_plot( "total_SIM_%s"            % t , output_directory, True  )
        pl.make_plot( "multiplicity_digis_%s"   % t , output_directory, True  )
        pl.make_plot( "multiplicity_simhits_%s" % t , output_directory, True  )

        pl.make_plot( "ADC_SimhitE_%s_layer"    % t , output_directory, False )
        pl.make_plot( "ADC_MIP_%s_layer"        % t , output_directory, False )
        pl.make_plot( "MIP_SimhitE_%s_layer"    % t , output_directory, False )

        pl.make_plot( "ADC_%s"                  % t , output_directory, True  )
        pl.make_plot( "MIP_%s"                  % t , output_directory, True  )
        pl.make_plot( "SIM_%s"                  % t , output_directory, True  )

def run_energy_resolution():
    output_directory = pl.specified_directory + "/energy_resolution"
    pu.create_directory(output_directory)

    #pl.make_simple_plot("MIP", output_directory, "odd_even")
    #pl.make_simple_plot("SIM", output_directory, "odd_even")

    pl.make_simple_plot("MIP", output_directory, "set1_set2")
    pl.make_simple_plot("MIP", output_directory, "set0")

    pl.make_simple_plot("SIM", output_directory, "set1_set2")
    pl.make_simple_plot("SIM", output_directory, "set0")

    pl.make_simple_plot("ENE", output_directory, "set1_set2")
    pl.make_simple_plot("ENE", output_directory, "set0")

def run_hit_distribution():
    output_directory = pl.specified_directory + "/hit_distributions"
    pu.create_directory(output_directory)

    pl.make_plot( "hEta", output_directory, False )
    pl.make_plot( "hPhi", output_directory, False )

    analyzer = an.HitAnalyzer(pl.myRootfiles, output_directory, tags)
    analyzer.loop()

def run_fitters_and_summary():
    output_directory = pl.specified_directory + "/linear_fit"
    pu.create_directory(output_directory)

    fit_result = pu.fit_result

    # report
    for energyType in ["MIP", "SIM"]:
        for tag in ["E20", "E60", "E100", "E175", "E225", "E300"]:
            print energyType, tag
            #print fit_result[energyType]["odd" ][tag]
            #print fit_result[energyType]["even"][tag]
            print fit_result[energyType]["set0"][tag]
            print fit_result[energyType]["set1"][tag]
            print fit_result[energyType]["set2"][tag]
            print ""

    # linear fit
    #pl.run_linear_fit(output_directory, "even", fit_result["MIP"]["odd" ], fit_result["SIM"]["odd" ])
    #pl.run_linear_fit(output_directory, "odd" , fit_result["MIP"]["even"], fit_result["SIM"]["even"])
    pl.run_linear_fit(output_directory, "set0", fit_result["MIP"]["set0"], fit_result["SIM"]["set0"])
    pl.run_linear_fit(output_directory, "set0_set1", fit_result["MIP"]["set1"], fit_result["SIM"]["set0"])
    pl.run_linear_fit(output_directory, "set0_set2", fit_result["MIP"]["set2"], fit_result["SIM"]["set0"])

    # summary
    labels = ["#sigma#left(E_{set1}#right) / #bar{E}_{set1}", "#sigma#left(E_{set2}#right) / #bar{E}_{set2}"]
    labels = ["Resolution of E_set1", "Resolution of E_set2"]
    pl.run_resolution_summary(pl.specified_directory, labels, fit_result["ENE"]["set1"], fit_result["ENE"]["set2"])

#----------------------------------------------------------------------------------------------------

#}}}

def run_manager(myfin, mydin):
    # set parameters
    pl.tags = tags
    pl.myRootfiles = myfin
    pl.specified_directory = mydin
    pl.fit_constraints = fit_constraints
    for tag in tags: pl.label[tag] = tag.split("E")[1] + " GeV"

    # prepare the main output directory
    pu.create_directory( pl.specified_directory )

    # runners
    #run_hit_distribution()
    #run_logitudinal_profile()
    run_energy_resolution()

#----------------------------------------------------------------------------------------------------

if __name__ == "__main__":

    tags = ["E300", "E100", "E20"]
    fit_constraints = m.fit_constraints_v1
    run_manager( m.input_files["R90To130"], eos + "/" + "R90To130_v2p1" )

    tags = ["E225", "E175", "E60"]
    fit_constraints = m.fit_constraints_v2
    run_manager( m.input_files["R90To130_v2"], eos + "/" + "R90To130_v2p1" )

    run_fitters_and_summary()

