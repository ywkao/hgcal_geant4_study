#!/usr/bin/env python2
import json
import toolbox.plot_utils as pu
import toolbox.MetaData as m
import toolbox.analyzer as an
import toolbox.plotter as pl

eos = "./eos"

#----------------------------------------------------------------------------------------------------

def run_logitudinal_profile():
    output_directory = m.specified_directory + "/longitudinal_profiles"
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
    """ extract histograms, perform gaussian fit, and store fit result """
    output_directory = m.specified_directory + "/energy_resolution"
    pu.create_directory(output_directory)

    if m.enable_check_odd_even:
        pl.make_simple_plot("MIP", output_directory, "odd_even")
        pl.make_simple_plot("SIM", output_directory, "odd_even")

    pl.make_simple_plot("MIP", output_directory, "set1_set2")
    pl.make_simple_plot("MIP", output_directory, "set0")

    pl.make_simple_plot("SIM", output_directory, "set1_set2")
    pl.make_simple_plot("SIM", output_directory, "set0")

    pl.make_simple_plot("ENE", output_directory, "set1_set2")
    pl.make_simple_plot("ENE", output_directory, "set0")

def run_hit_distribution():
    output_directory = m.specified_directory + "/hit_distributions"
    pu.create_directory(output_directory)

    pl.make_plot( "hEta", output_directory, False )
    pl.make_plot( "hPhi", output_directory, False )

def run_hit_analyzer():
    output_directory = m.specified_directory + "/hit_distributions"
    pu.create_directory(output_directory)

    analyzer = an.HitAnalyzer(pl.myRootfiles, output_directory, pl.tags)
    analyzer.loop()

def run_fitters_and_summary():
    output_directory = m.specified_directory + "/linear_fit"
    pu.create_directory(output_directory)

    #fit_result = pu.fit_result
    #fit_result_goodness = pu.fit_result_goodness

    fit_result, fit_result_goodness = {}, {}
    with open(m.json_fit_parameters, 'r') as f: fit_result = json.load(f)
    with open(m.json_fit_parameters_goodness, 'r') as f: fit_result_goodness = json.load(f)

    #--------------------------------------------------
    # report
    #--------------------------------------------------
    if False:
        print "[INFO] report fit result:"
        for energyType in ["MIP", "SIM", "ENE"]:
            for tag in ["E20", "E60", "E100", "E175", "E225", "E300"]:
                print energyType, tag
                print fit_result[energyType]["set0"][tag]
                print fit_result[energyType]["set1"][tag]
                print fit_result[energyType]["set2"][tag]
                print ""

    #--------------------------------------------------
    # linear fit
    #--------------------------------------------------
    # x = E[MIP], y = E[GeV]
    #pl.run_linear_fit(output_directory, fit_result, "set0"     , fit_result["MIP"]["set0"], fit_result["ENE"]["set0"])
    #pl.run_linear_fit(output_directory, fit_result, "set0_set1", fit_result["MIP"]["set1"], fit_result["ENE"]["set1"])
    #pl.run_linear_fit(output_directory, fit_result, "set0_set2", fit_result["MIP"]["set2"], fit_result["ENE"]["set2"])

    # x = E_beam, y = E[MIP]
    pl.run_linear_fit(output_directory, fit_result, "set0_beam", 0, fit_result["MIP"]["set0"])
    pl.run_linear_fit(output_directory, fit_result, "set1_beam", 0, fit_result["MIP"]["set1"])
    pl.run_linear_fit(output_directory, fit_result, "set2_beam", 0, fit_result["MIP"]["set2"])

    m.json_fit_parameters = "./toolbox/test.json"
    run_register_fit_parameters(fit_result) # record fit result

    #--------------------------------------------------
    # summary
    #--------------------------------------------------
    if m.enable_check_odd_even:
        m.energy_type = "odd_even_MIPs"
        m.labels = ["E_odd_MIP", "E_even_MIP"]
        pl.run_summary("pvalue", fit_result_goodness["MIP"]["odd"], fit_result_goodness["MIP"]["even"])
        pl.run_summary("chi2ndf", fit_result_goodness["MIP"]["odd"], fit_result_goodness["MIP"]["even"])
        pl.run_summary(m.type_resolution, fit_result["MIP"]["odd"], fit_result["MIP"]["even"])

    # Remark: less meaningful to use resolution of E[MIP] directly
    # Resolution of E reco from E[MIP] / slope
    m.energy_type = "set1_set2_MIPs_slope"
    m.labels = ["E_default_slope_method", "E_alternative_slope_method"]
    m.resolution["set1"][m.energy_type] = {}
    m.resolution["set2"][m.energy_type] = {}
    pl.run_summary("pvalue", fit_result_goodness["MIP"]["set1"], fit_result_goodness["MIP"]["set2"])
    pl.run_summary("chi2ndf", fit_result_goodness["MIP"]["set1"], fit_result_goodness["MIP"]["set2"])
    pl.run_summary(m.type_resolution, fit_result["E_slope_method"]["set1"], fit_result["E_slope_method"]["set2"])
    pl.run_summary(m.bias+"_MIP", fit_result["E_slope_method"]["set1"], fit_result["E_slope_method"]["set2"])

    # Resolution of E reco from the dE/dx method 
    m.energy_type = "set1_set2_MeV"
    m.labels = ["E_default_dEdx_method", "E_alternative_dEdx_method"]
    m.resolution["set1"][m.energy_type] = {}
    m.resolution["set2"][m.energy_type] = {}
    pl.run_summary("pvalue", fit_result_goodness["ENE"]["set1"], fit_result_goodness["ENE"]["set2"])
    pl.run_summary("chi2ndf", fit_result_goodness["ENE"]["set1"], fit_result_goodness["ENE"]["set2"])
    pl.run_summary(m.type_resolution, fit_result["ENE"]["set1"], fit_result["ENE"]["set2"])
    pl.run_summary(m.bias+"_MeV", fit_result["ENE"]["set1"], fit_result["ENE"]["set2"])

    m.energy_type = "changes_in_resolution"
    m.labels = ["E_default", "E_alternative"]
    pl.run_summary("changes_in_resolution", m.resolution["set1"], m.resolution["set2"])

    #pl.run_summary("linearity", m.resolution["set1"], m.resolution["set2"])

#----------------------------------------------------------------------------------------------------

def run_register_fit_parameters(fit_result, fit_result_goodness=""):

    with open(m.json_fit_parameters, 'w') as f:
        json.dump(fit_result, f, sort_keys=True, indent=4)
        f.write("\n")

    if isinstance(fit_result_goodness, dict):
        with open(m.json_fit_parameters_goodness, 'w') as f:
            json.dump(fit_result_goodness, f, sort_keys=True, indent=4)
            f.write("\n")

#----------------------------------------------------------------------------------------------------

def run_manager(myfin, tags, fit_constraints):
    # set parameters
    pl.tags = tags
    pl.myRootfiles = myfin
    pl.fit_constraints = fit_constraints
    for tag in tags: pl.label[tag] = tag.split("E")[1] + " GeV" if "E" in tag else tag

    # runners
    # run_hit_analyzer()
    # run_hit_distribution()
    # run_logitudinal_profile()
    run_energy_resolution()

