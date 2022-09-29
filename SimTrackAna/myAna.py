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
        #pl.make_plot( "total_SIM_%s"            % t , output_directory, True  )

        pl.make_plot( "efficiency_linear_track"     , output_directory, True  )

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

    if enable_check_odd_even:
        pl.make_simple_plot("MIP", output_directory, "odd_even")
        pl.make_simple_plot("SIM", output_directory, "odd_even")

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

def run_hit_analyzer():
    output_directory = pl.specified_directory + "/hit_distributions"
    pu.create_directory(output_directory)

    analyzer = an.HitAnalyzer(pl.myRootfiles, output_directory, tags)
    analyzer.loop()

def run_fitters_and_summary():
    output_directory = pl.specified_directory + "/linear_fit"
    pu.create_directory(output_directory)

    fit_result = pu.fit_result
    fit_result_goodness = pu.fit_result_goodness

    # report
    if False:
        for energyType in ["MIP", "SIM"]:
            for tag in ["E20", "E60", "E100", "E175", "E225", "E300"]:
                print energyType, tag
                print fit_result[energyType]["set0"][tag]
                print fit_result[energyType]["set1"][tag]
                print fit_result[energyType]["set2"][tag]
                print ""

    # linear fit
    pl.run_linear_fit(output_directory, "set0", fit_result["MIP"]["set0"], fit_result["SIM"]["set0"])
    pl.run_linear_fit(output_directory, "set0_set1", fit_result["MIP"]["set1"], fit_result["SIM"]["set0"])
    pl.run_linear_fit(output_directory, "set0_set2", fit_result["MIP"]["set2"], fit_result["SIM"]["set0"])

    # summary
    if enable_check_odd_even:
        m.energy_type = "odd_even_MIPs"
        m.labels = ["E_odd_MIP", "E_even_MIP"]
        pl.run_summary("pvalue", fit_result_goodness["MIP"]["odd"], fit_result_goodness["MIP"]["even"])
        pl.run_summary("chi2ndf", fit_result_goodness["MIP"]["odd"], fit_result_goodness["MIP"]["even"])
        pl.run_summary(m.type_resolution, fit_result["MIP"]["odd"], fit_result["MIP"]["even"])

    m.energy_type = "set1_set2_MeV"
    m.labels = ["E_set1_MeV", "E_set2_MeV"]
    pl.run_summary("pvalue", fit_result_goodness["ENE"]["set1"], fit_result_goodness["ENE"]["set2"])
    pl.run_summary("chi2ndf", fit_result_goodness["ENE"]["set1"], fit_result_goodness["ENE"]["set2"])
    pl.run_summary(m.type_resolution, fit_result["ENE"]["set1"], fit_result["ENE"]["set2"])

    m.energy_type = "set1_set2_MIPs"
    m.labels = ["E_set1_MIP", "E_set2_MIP"]
    pl.run_summary("pvalue", fit_result_goodness["MIP"]["set1"], fit_result_goodness["MIP"]["set2"])
    pl.run_summary("chi2ndf", fit_result_goodness["MIP"]["set1"], fit_result_goodness["MIP"]["set2"])
    pl.run_summary(m.type_resolution, fit_result["MIP"]["set1"], fit_result["MIP"]["set2"])

#----------------------------------------------------------------------------------------------------

#}}}

def run_manager(myfin, mydin, tags, fit_constraints):
    # set parameters
    pl.tags = tags
    pl.myRootfiles = myfin
    pl.specified_directory = mydin
    pl.fit_constraints = fit_constraints
    for tag in tags: pl.label[tag] = tag.split("E")[1] + " GeV"

    # prepare the main output directory
    pu.create_directory( pl.specified_directory )

    # runners
    #run_hit_analyzer()
    #run_hit_distribution()
    #run_logitudinal_profile()
    run_energy_resolution()

#----------------------------------------------------------------------------------------------------

def perform_unclustered_study():
    m.type_resolution = "resolution_unclustered"
    target_directory = "R90To130_v1p6"
    target_directory = "R90To130_unclustered"
    output_directory = eos + "/" + target_directory

    tags = ["E300", "E100", "E20"]
    fit_constraints = m.fit_constraints_v1p1
    run_manager( m.input_files["R90To130_v1p1"], output_directory, tags, fit_constraints )

    tags = ["E225", "E175", "E60"]
    fit_constraints = m.fit_constraints_v1p2
    run_manager( m.input_files["R90To130_v1p2"], output_directory, tags, fit_constraints )

    run_fitters_and_summary()

def perform_clustered_study():
    m.type_resolution = "resolution_clustered"
    target_directory = "R90To130_v2p6"
    target_directory = "R90To130_clustered"
    output_directory = eos + "/" + target_directory

    tags = ["E300", "E100", "E20"]
    fit_constraints = m.fit_constraints_v2p1
    run_manager( m.input_files["R90To130_v2p1"], output_directory, tags, fit_constraints )

    tags = ["E225", "E175", "E60"]
    fit_constraints = m.fit_constraints_v2p2
    run_manager( m.input_files["R90To130_v2p2"], output_directory, tags, fit_constraints )

    run_fitters_and_summary()

#----------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    enable_check_odd_even = False
    perform_unclustered_study()
    #perform_clustered_study()

