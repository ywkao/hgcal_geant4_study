#!/usr/bin/env python2

#====================================================================================================
# containers for parameters to be set
#====================================================================================================
type_resolution = ""
energy_type = "" # for run_summary
labels = [] # for run_summary

resolution = {}
resolution = {"set1":{}, "set2":{}}

#====================================================================================================
# toolbox/plotter.py: run_summary
#====================================================================================================
# title and range
draw_options_for_run_summary = {
    "pvalue": {
        "ytitle"         : "p-value of gaussian fit",
        "yrange"         : [1e-3, 1.0],
        "leg_pos"        : [0.49, 0.25, 0.89, 0.45],
        "leg_option"     : "p",
        "draw_goodness"  : True,
        "draw_lower_pad" : False,
        "reference_line" : 0.05,
        "useLog"         : 1,
    },

    "chi2ndf": {
        "ytitle"         : "chi2/ndf of gaussian fit",
        "yrange"         : [0., 4.0],
        "leg_pos"        : [0.49, 0.65, 0.89, 0.85],
        "leg_option"     : "p",
        "draw_goodness"  : True,
        "draw_lower_pad" : False,
        "reference_line" : 1.0,
        "useLog"         : 0,
    },

    "resolution_unclustered": {
        "ytitle"            : "#sigma#left(E#right) / #LTE#GT",
        #"yrange"            : [0.018, 0.043],
        #"yrange"            : [0.010, 0.050],
        "yrange"            : [0.000, 0.160],
        "leg_pos"           : [0.49, 0.65, 0.89, 0.85],
        "leg_option"        : "ep",
        "draw_goodness"     : True,
        "draw_lower_pad"    : True,
        "reference_line"    : 1.0,
        "useLog"            : 0,

        "linear_fit_yrange" : [0, 400],
        #"use_bin_width"     : True,
        "use_bin_width"     : False,
    },

    "resolution_clustered": {
        "ytitle"            : "#sigma#left(E#right) / #LTE#GT",
        #"yrange"            : [0.000, 0.040],
        "yrange"            : [0.000, 0.160],
        "leg_pos"           : [0.49, 0.65, 0.89, 0.85],
        "leg_option"        : "ep",
        "draw_goodness"     : True,
        "draw_lower_pad"    : True,
        "reference_line"    : 1.0,
        "useLog"            : 0,

        "linear_fit_yrange" : [0, 400],
        "use_bin_width"     : False,
    },

    "resolution_unclustered_FIT": {
        "ytitle"            : "#sigma#left(E#right) / #LTE#GT",
        "yrange"            : [0.00, 0.15],
        "leg_pos"           : [0.37, 0.67, 0.87, 0.87],
        #"leg_pos"           : [0.15, 0.47, 0.40, 0.87],
        "leg_option"        : "ep",
        "draw_goodness"     : True,
        "draw_lower_pad"    : True,
        "reference_line"    : 1.0,
        "useLog"            : 0,

        "linear_fit_yrange" : [0, 400],
        "use_bin_width"     : False,
    },

    "resolution_clustered_FIT": {
        "ytitle"            : "#sigma#left(E#right) / #LTE#GT",
        "yrange"            : [0.00, 0.15],
        "leg_pos"           : [0.37, 0.67, 0.87, 0.87],
        "leg_option"        : "ep",
        "draw_goodness"     : True,
        "draw_lower_pad"    : True,
        "reference_line"    : 1.0,
        "useLog"            : 0,

        "linear_fit_yrange" : [0, 400],
        "use_bin_width"     : False,
    },

    "bias_unclustered_MIP": {
        #"ytitle"            : "Relative difference of E[MIP]",
        #"ytitle"            : "#frac{#LTE/m#GT#minusE_{beam}}{E_{beam}}",
        "ytitle"            : "(#LTE/m#GT#minusE_{beam})/E_{beam}",
        "yrange"            : [-0.10, 0.15],
        "leg_pos"           : [0.49, 0.65, 0.89, 0.85],
        "leg_option"        : "ep",
        "draw_goodness"     : True,
        "draw_lower_pad"    : False,
        "reference_line"    : 0.0,
        "useLog"            : 0,
    },

    "bias_unclustered_MeV": {
        #"ytitle"            : "#frac{E_{reco.}#minusE_{beam}}{E_{beam}}",
        "ytitle"            : "(E_{reco.}#minusE_{beam})/E_{beam}",
        "yrange"            : [-0.06, 0.14],
        "leg_pos"           : [0.49, 0.65, 0.89, 0.85],
        "leg_option"        : "ep",
        "draw_goodness"     : True,
        "draw_lower_pad"    : False,
        "reference_line"    : 0.0,
        "useLog"            : 0,
    },

    "bias_clustered_MIP": {
        #"ytitle"            : "#frac{#LTE/m#GT#minusE_{beam}}{E_{beam}}",
        "ytitle"            : "(#LTE/m#GT#minusE_{beam})/E_{beam}",
        "yrange"            : [-0.3, 0.3],
        "leg_pos"           : [0.49, 0.65, 0.89, 0.85],
        "leg_option"        : "ep",
        "draw_goodness"     : True,
        "draw_lower_pad"    : False,
        "reference_line"    : 0.0,
        "useLog"            : 0,
    },

    "bias_clustered_MeV": {
        #"ytitle"            : "#frac{E_{reco.}#minusE_{beam}}{E_{beam}}",
        "ytitle"            : "(E_{reco.}#minusE_{beam})/E_{beam}",
        "yrange"            : [-0.3, 0.3],
        "leg_pos"           : [0.49, 0.65, 0.89, 0.85],
        "leg_option"        : "ep",
        "draw_goodness"     : True,
        "draw_lower_pad"    : False,
        "reference_line"    : 0.0,
        "useLog"            : 0,
    },

    "changes_in_resolution": {
        #"ytitle"            : "#frac{Res.(after) #minus Res.(before)}{Res.(before)}",
        "ytitle"            : "(Res.(after) #minus Res.(before))/Res.(before)",
        "yrange"            : [-0.5, 0.5],
        "leg_pos"           : [0.49, 0.65, 0.89, 0.85],
        "leg_option"        : "ep",
        "draw_goodness"     : True,
        "draw_lower_pad"    : False,
        "reference_line"    : 0.0,
        "useLog"            : 0,

        "linear_fit_yrange" : [0, 400],
        "use_bin_width"     : False,
    },
}


#====================================================================================================
# fit parameters
#====================================================================================================

linear_fit_parameter = {
    "set0" : {
        "linear_fit_xrange" : [0, 30000],
        "linear_fit_yrange" : [0, 400],
    },

    "set0_set1" : {
        "linear_fit_xrange" : [0, 20000],
        "linear_fit_yrange" : [0, 400],
    },

    "set0_set2" : {
        "linear_fit_xrange" : [0, 20000],
        "linear_fit_yrange" : [0, 400],
    },

    "set0_beam" : {
        "linear_fit_xrange" : [0, 350], # GeV
        "linear_fit_yrange" : [0, 30000], # MIPs
    },

    "set1_beam" : {
        "linear_fit_xrange" : [0, 350], # GeV
        "linear_fit_yrange" : [0, 15000], # MIPs
    },

    "set2_beam" : {
        "linear_fit_xrange" : [0, 350], # GeV
        "linear_fit_yrange" : [0, 15000], # MIPs
    },
}

rebin_factor = {
    "resolution_unclustered":{
        "set0":{
            "E300" : 5,
            "E225" : 5,
            "E175" : 5,
            "E100" : 5,
            "E60"  : 2,
            "E20"  : 2,
        },
        "set1_set2":{
            "E300" : 10,
            "E225" : 10,
            "E175" : 10,
            "E100" : 10,
            "E60"  : 4,
            "E20"  : 4,
        },
    },
    "resolution_clustered":{
        "set0":{
            "E300" : 5,
            "E225" : 5,
            "E175" : 5,
            "E100" : 5,
            "E60"  : 2,
            "E20"  : 2,
        },
        "set1_set2":{
            "E300" : 10,
            "E225" : 10,
            "E175" : 10,
            "E100" : 10,
            "E60"  : 4,
            "E20"  : 4,
        },
    },
}

# E = 300, 100, 20
fit_constraints_v0 = {
    "set0" : {
        "MIP" : { "xRanges" : [[0, 30000], [0, 30000], [0, 30000], [0, 30000]] },
        "SIM" : { "xRanges" : [[0, 300], [0, 300], [0, 300]] },
        "ENE" : { "xRanges" : [[0, 300], [0, 300], [0, 300]] },
    },

    "set1set2" : {
        "MIP" : { "xRanges" : [[0, 30000], [0, 30000], [0, 30000], [0, 30000]] },
        "SIM" : { "xRanges" : [[0, 300], [0, 300], [0, 300]] },
        "ENE" : { "xRanges" : [[0, 300], [0, 300], [0, 300]] },
    },
}

fit_constraints_set0 = {
    "MIP" : { "xRanges" : [0, 30000] },
    "SIM" : { "xRanges" : [0, 300] },
    "ENE" : { "xRanges" : [0, 300] },
}

# E = 300, 100, 20
fit_constraints_v1p1 = {
    "set0" : {
        "MIP" : { "xRanges" : [[24000, 40000], [7000, 15000], [1000, 5000]] },
        "SIM" : { "xRanges" : [[100, 300], [40, 240], [0, 200]] },
        "ENE" : { "xRanges" : [[270, 400], [80, 180], [10, 60]] },
        #"ENE" : { "xRanges" : [[0, 300], [0, 300], [0, 300]] },
    },

    "set1set2" : {
        "MIP" : { "xRanges" : [[11500, 20500], [3000, 12000], [0, 9000]], },
        "SIM" : { "xRanges" : [[0, 300], [0, 300], [0, 300]] },
        "ENE" : { "xRanges" : [[270, 440], [80, 180], [10, 60]] },
        #"ENE" : { "xRanges" : [[0, 400000], [0, 400000], [0, 400000]] },
        #"ENE" : { "xRanges" : [[140, 200], [58, 98], [25, 65]] },
    },
}

# E = 225, 175, 60
fit_constraints_v1p2 = {
    "set0" : {
        "MIP" : { "xRanges" : [[17000, 30000], [12000, 25000], [3500, 10000]] },
        "SIM" : { "xRanges" : [[80, 280], [60, 260], [0, 300]] },
        "ENE" : { "xRanges" : [[200, 350], [150, 300], [45, 120]] },
        #"ENE" : { "xRanges" : [[0, 300], [0, 300], [0, 300]] },
    },

    "set1set2" : {
        "MIP" : { "xRanges" : [[8000, 17000], [6000, 15000], [0, 9000]], },
        "SIM" : { "xRanges" : [[0, 300], [0, 300], [0, 300]] },
        "ENE" : { "xRanges" : [[200, 350], [150, 300], [45, 120]] },
        #"ENE" : { "xRanges" : [[0, 400000], [0, 400000], [0, 400000]] },
        #"ENE" : { "xRanges" : [[110, 160], [90, 130], [40, 80]] },
    },
}

# E = 300, 100, 20
fit_constraints_v2p1 = {
    "set0" : {
        "MIP" : { "xRanges" : [[20000, 30000], [6500, 11500], [1000, 5000]] },
        "SIM" : { "xRanges" : [[0, 300], [0, 300], [0, 300]] },
        "ENE" : { "xRanges" : [[220, 400], [70, 170], [10, 60]] },
        #"ENE" : { "xRanges" : [[0, 400000], [0, 400000], [0, 400000]] },
        #"ENE" : { "xRanges" : [[0, 300], [0, 300], [0, 300]] },
    },

    "set1set2" : {
        "MIP" : { "xRanges" : [[9000, 20000], [2500, 13000], [400, 2400]] },
        #"MIP" : { "xRanges" : [[9000, 15000], [3000, 6000], [0, 2000]] },
        "SIM" : { "xRanges" : [[0, 300], [0, 300], [0, 300]] },
        "ENE" : { "xRanges" : [[220, 400], [70, 170], [10, 60]] },
        #"ENE" : { "xRanges" : [[0, 400000], [0, 400000], [0, 400000]] },
        #"ENE" : { "xRanges" : [[37, 52], [22, 37], [16, 31]] },
    },
}

# E = 225, 175, 60
fit_constraints_v2p2 = {
    "set0" : {
        "MIP" : { "xRanges" : [[15000, 30000], [11000, 21000], [3000, 8000]] },
        "SIM" : { "xRanges" : [[0, 300], [0, 300], [0, 300]] },
        "ENE" : { "xRanges" : [[170, 320], [130, 230], [30, 130]] },
        #"ENE" : { "xRanges" : [[0, 400000], [0, 400000], [0, 400000]] },
        #"ENE" : { "xRanges" : [[0, 300], [0, 300], [0, 300]] },
    },

    "set1set2" : {
        "MIP" : { "xRanges" : [[6500, 14500], [5000, 13000], [1500, 6000]] },
        "SIM" : { "xRanges" : [[0, 300], [0, 300], [0, 300]] },
        "ENE" : { "xRanges" : [[170, 320], [130, 230], [30, 130]] },
        #"ENE" : { "xRanges" : [[0, 400000], [0, 400000], [0, 400000]] },
        #"ENE" : { "xRanges" : [[31, 46], [27, 42], [18, 33]] },
    },
}


#====================================================================================================
# input files
#====================================================================================================
input_files = {
    "R35To60" : [
        "rootfiles/geantoutput_D86_R35To60_E300.root",
        "rootfiles/geantoutput_D86_R35To60_E100.root",
        "rootfiles/geantoutput_D86_R35To60_E20.root",
    ],

    "R80To100" : [
        "rootfiles_original/geantoutput_D86_R80To100_E300.root",
        "rootfiles_original/geantoutput_D86_R80To100_E100.root",
        "rootfiles_original/geantoutput_D86_R80To100_E20.root",
    ],

    "R80To100_v2" : [
        "rootfiles_eta_cutFree_v4/geantoutput_D86_R80To100_E300.root",
        "rootfiles_eta_cutFree_v4/geantoutput_D86_R80To100_E100.root",
        "rootfiles_eta_cutFree_v4/geantoutput_D86_R80To100_E20.root",
    ],

    "R90To130_v1p1" : [
        "rootfiles_hit_position/20240522_R90To130/geantoutput_D86_R90To130_E300.root",
        "rootfiles_hit_position/20240522_R90To130/geantoutput_D86_R90To130_E100.root",
        "rootfiles_hit_position/20240522_R90To130/geantoutput_D86_R90To130_E20.root",
        #"rootfiles_hit_position/20221011_R90To130/geantoutput_D86_R90To130_E300.root",
        #"rootfiles_hit_position/20221011_R90To130/geantoutput_D86_R90To130_E100.root",
        #"rootfiles_hit_position/20221011_R90To130/geantoutput_D86_R90To130_E20.root",
    ],

    "R90To130_v1p2" : [
        "rootfiles_hit_position/20240522_R90To130/geantoutput_D86_R90To130_E225.root",
        "rootfiles_hit_position/20240522_R90To130/geantoutput_D86_R90To130_E175.root",
        "rootfiles_hit_position/20240522_R90To130/geantoutput_D86_R90To130_E60.root",
        #"rootfiles_hit_position/20221011_R90To130/geantoutput_D86_R90To130_E225.root",
        #"rootfiles_hit_position/20221011_R90To130/geantoutput_D86_R90To130_E175.root",
        #"rootfiles_hit_position/20221011_R90To130/geantoutput_D86_R90To130_E60.root",
    ],

    # v2: SR = linear track
    "R90To130_v2p1" : [
        "rootfiles_hit_position/20240522_R90To130_AM_algorithm/geantoutput_D86_R90To130_E300.root",
        "rootfiles_hit_position/20240522_R90To130_AM_algorithm/geantoutput_D86_R90To130_E100.root",
        "rootfiles_hit_position/20240522_R90To130_AM_algorithm/geantoutput_D86_R90To130_E20.root",
        #"rootfiles_hit_position/20220929_R90To130_AM_algorithm/geantoutput_D86_R90To130_E300.root",
        #"rootfiles_hit_position/20220929_R90To130_AM_algorithm/geantoutput_D86_R90To130_E100.root",
        #"rootfiles_hit_position/20220929_R90To130_AM_algorithm/geantoutput_D86_R90To130_E20.root",
    ],

    "R90To130_v2p2" : [
        "rootfiles_hit_position/20240522_R90To130_AM_algorithm/geantoutput_D86_R90To130_E225.root",
        "rootfiles_hit_position/20240522_R90To130_AM_algorithm/geantoutput_D86_R90To130_E175.root",
        "rootfiles_hit_position/20240522_R90To130_AM_algorithm/geantoutput_D86_R90To130_E60.root",
        #"rootfiles_hit_position/20220929_R90To130_AM_algorithm/geantoutput_D86_R90To130_E225.root",
        #"rootfiles_hit_position/20220929_R90To130_AM_algorithm/geantoutput_D86_R90To130_E175.root",
        #"rootfiles_hit_position/20220929_R90To130_AM_algorithm/geantoutput_D86_R90To130_E60.root",
    ],

    "R80To130" : [
        "rootfiles_R80To130/geantoutput_D86_R80To130_E300.root",
        "rootfiles_R80To130/geantoutput_D86_R80To130_E100.root",
        "rootfiles_R80To130/geantoutput_D86_R80To130_E20.root",
    ],

    "R80To130_v2" : [
        "rootfiles_R80To130/geantoutput_D86_R80To130_E225.root",
        "rootfiles_R80To130/geantoutput_D86_R80To130_E175.root",
        "rootfiles_R80To130/geantoutput_D86_R80To130_E60.root",
    ],

    "R80To130_v3" : [
        "rootfiles_hit_position/geantoutput_D86_R80To130_E300.root",
        "rootfiles_hit_position/geantoutput_D86_R80To130_E100.root",
        "rootfiles_hit_position/geantoutput_D86_R80To130_E20.root",
    ],

    "R80To150" : [
        "rootfiles_eta_cutFree_v4/geantoutput_D86_R80To150_E300.root",
        "rootfiles_eta_cutFree_v4/geantoutput_D86_R80To150_E100.root",
        "rootfiles_eta_cutFree_v4/geantoutput_D86_R80To150_E20.root",
        #"rootfiles_eta_cutFree_v4/geantoutput_D86_R80To150_E300.root",
        #"rootfiles_eta_cutFree_v4/geantoutput_D86_R80To150_E100.root",
        #"rootfiles_eta_cutFree_v4/geantoutput_D86_R80To150_E20.root",
    ],

    "R80To150_v2" : [
        "rootfiles_eta_cutFree_v4/geantoutput_D86_R80To150_E225.root",
        "rootfiles_eta_cutFree_v4/geantoutput_D86_R80To150_E175.root",
        "rootfiles_eta_cutFree_v4/geantoutput_D86_R80To150_E60.root",
    ],

    "eta_scanning" : [
        "rootfiles_eta_scanning/geantoutput_D86_R80To110_E100.root",
        "rootfiles_eta_scanning/geantoutput_D86_R80To120_E100.root",
        "rootfiles_eta_scanning/geantoutput_D86_R80To130_E100.root",
        "rootfiles_eta_scanning/geantoutput_D86_R80To140_E100.root",
    ],

    "muon" : [
        "rootfiles/geantoutput_D86_muon_E100.root",
    ],

    "ProdCut1mm" : [
        "rootfiles/geantoutput_D86_R80To100_E100_nominal.root",
        "rootfiles/geantoutput_D86_R80To100_E100_ProdCut_electron_1mm.root",
        "rootfiles/geantoutput_D86_R80To100_E100_ProdCut_photon_1mm.root",
        "rootfiles/geantoutput_D86_R80To100_E100_ProdCut_egamma_1mm.root",
    ],

    "ProdCut100mm" : [
        "rootfiles/geantoutput_D86_R80To100_E100_nominal.root",
        "rootfiles/geantoutput_D86_R80To100_E100_ProdCut_electron_100mm.root",
        "rootfiles/geantoutput_D86_R80To100_E100_ProdCut_photon_100mm.root",
        "rootfiles/geantoutput_D86_R80To100_E100_ProdCut_egamma_100mm.root",
    ],

    "ProdCut1000mm" : [
        "rootfiles/geantoutput_D86_R80To100_E100_nominal.root",
        "rootfiles/geantoutput_D86_R80To100_E100_ProdCut_electron_1000mm.root",
        "rootfiles/geantoutput_D86_R80To100_E100_ProdCut_photon_1000mm.root",
        "rootfiles/geantoutput_D86_R80To100_E100_ProdCut_egamma_1000mm.root",
    ],

    "turnOffCompton" : [
        "rootfiles/geantoutput_D86_R80To100_E100_nominal.root",
        "rootfiles/geantoutput_D86_R80To100_E100_turnOffComptionScattering_v2.root",
        "rootfiles/geantoutput_D86_R80To100_E100_turnOffComptionScattering_turnOffGammaConversion.root",
        "rootfiles/geantoutput_D86_R80To100_E100_ProdCut_electron_100mm.root",
    ],

    "extraStudy" : [
        "rootfiles/geantoutput_D86_R80To100_E100_nominal.root",
        "rootfiles/geantoutput_D86_R80To100_E100_turnOffComptionScattering_v2.root",
        "rootfiles/geantoutput_D86_R80To100_E100_turnOffGammaConversion.root",
        "rootfiles/geantoutput_D86_R80To100_E100_turnOffComptionScattering_turnOffGammaConversion.root",
    ],

    "turnOffCompton_ProdCutElectron" : [
        "rootfiles/geantoutput_D86_R80To100_E100_nominal.root",
        "rootfiles/geantoutput_D86_R80To100_E100_turnOffComptionScattering_v3.root",
        "rootfiles/geantoutput_D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_1mm.root",
        #"rootfiles/geantoutput_D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_5mm.root",
        "rootfiles/geantoutput_D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_10mm.root",
        #"rootfiles/geantoutput_D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_50mm.root",
        "rootfiles/geantoutput_D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_100mm.root",
        "rootfiles/geantoutput_D86_R80To100_E100_turnOffComptionScattering_ProdCut_electron_1000mm.root",
    ],

    "pcbStudy" : [
        "rootfiles/geantoutput_D86_R80To100_E100_nominal.root",
        "rootfiles/geantoutput_D86_R80To100_E100_turnOffComptionScattering_v3.root",
        "rootfiles/geantoutput_D86_R80To100_E100_airPCB_turnOffComptionScattering.root",
        "rootfiles/geantoutput_D86_R80To100_E100_PCB.root",
        #"rootfiles/geantoutput_D86_R80To100_E100_airPCB_turnOffComptionScattering_ProdCut_electron_5mm.root",
        #"rootfiles/geantoutput_D86_R80To100_E100_airPCB_turnOffComptionScattering_ProdCut_electron_10mm.root",
        #"rootfiles/geantoutput_D86_R80To100_E100_airPCB_turnOffComptionScattering_ProdCut_electron_100mm.root",
    ],

    "dEdxStudy" : [
        "rootfiles/geantoutput_D86_R80To100_E100_nominal.root",
        "rootfiles_dEdx_weights/geantoutput_D86_R80To100_E100_nominal.root",
    ],

    "X0_corrections" : [
        "rootfiles/geantoutput_D86_R80To100_E100_nominal.root",
        "rootfiles_dEdx_weights/geantoutput_D86_R80To100_E100_nominal_weighted_dEdx_X0.root",
        "rootfiles_dEdx_weights/geantoutput_D86_R80To100_E100_nominal_weighted_dEdx_divided_X0.root",
        "rootfiles_dEdx_weights/geantoutput_D86_R80To100_E100_nominal_weighted_X0.root",
        "rootfiles_dEdx_weights/geantoutput_D86_R80To100_E100_nominal_weighted_dEdx.root",
    ],
}

#====================================================================================================
# test beam result
#====================================================================================================
test_beam_result = {
    # layer depth is approximated from eyes
    "layer_depth" : [0.98  , 1.97  , 2.80  , 3.80  , 4.80  , 5.75  , 6.75    ,
                     7.75  , 8.60  , 9.60  , 10.50 , 11.50 , 12.40 , 13.40   ,
                     14.40 , 15.25 , 16.25 , 17.25 , 18.20 , 19.20 , 20.10   ,
                     21.25 , 22.20 , 23.40 , 24.25 , 25.60 , 26.50 , 27.80 ] ,

    "mips" : {
        "data" : {
                20: [ 27.97529793   , 56.94511032   , 109.28077698  , 136.3278656   , 185.52850342  , 183.97288513  , 190.09320068  ,
                    169.14602661    , 172.43197632  , 125.23709869  , 123.59194946  , 86.94829559   , 80.3841095    , 55.09900284   ,
                    49.03165436     , 30.44026566   , 28.33568954   , 18.1529274    , 21.78116608   , 9.80546951    , 7.89422131    ,
                    4.74877453      , 3.71464896    , 2.26393962    , 2.68312931    , 1.84798074    , 2.96802163    , 2.07398152]   ,
                100: [ 46.63975525  , 114.66628265  , 265.69363403  , 386.59579468  , 621.76507568  , 755.52404785  , 828.36853027  ,
                    850.06390381    , 961.36834717  , 762.60217285  , 810.68963623  , 615.45013428  , 609.97467041  , 436.13134766  ,
                    409.57852173    , 271.23117065  , 254.13935852  , 169.89068604  , 159.09002686  , 98.19945526   , 78.01111603   ,
                    48.89661789     , 38.81033707   , 23.93457222   , 23.25452232   , 14.55828953   , 13.5385313    , 8.12277794]   ,
                300: [  61.28678131 , 170.39982605  , 445.28921509  , 696.44915771  , 1232.91577148 , 1706.63549805 , 2004.63244629 ,
                    2117.43017578   , 2636.14453125 , 2165.76635742 , 2491.61181641 , 2010.32653809 , 2153.66015625 , 1509.47131348 ,
                    1532.20666504   , 1083.24182129 , 985.58740234  , 675.88018799  , 626.57421875  , 408.73651123  , 326.54840088  ,
                    208.37161255    , 165.11035156  , 102.31463623  , 98.6309967    , 62.07058716   , 52.98619843   , 30.69275284]
                },
        
        "mc" : {
                20: [ 21.50681686  , 52.96774292   , 105.72846985  , 139.59152222  , 190.9624176   , 188.5632019   , 212.36528015  ,
                    178.33622742   , 179.30996704  , 136.88346863  , 128.64570618  , 92.97825623   , 83.23855591   , 57.88684082   ,
                    50.43962479    , 34.07664871   , 29.20438576   , 19.41445923   , 16.39099503   , 10.66606712   , 7.91118002    ,
                    5.2383852      , 3.97883224    , 2.6746645     , 2.32645464    , 1.55437088    , 1.33136487    , 0.88430882]   ,
                100: [ 32.89741516 , 100.46727753  , 243.90615845  , 387.94775391  , 623.9732666   , 719.44714355  , 921.58850098  ,
                    876.69610596   , 977.50720215  , 823.13946533  , 837.11102295  , 652.08312988  , 622.23828125  , 457.28369141  ,
                    417.88778687   , 294.85919189  , 260.6965332   , 178.53686523  , 153.69523621  , 102.70782471  , 76.91766357   ,
                    50.60368347    , 38.71577835   , 25.14162064   , 21.87405968   , 13.90151215   , 11.96078873   , 7.24694681]   ,
                300: [  41.1002121 , 141.27975464  , 386.300354    , 688.69451904  , 1221.99047852 , 1552.28918457 , 2161.61791992 ,
                    2233.05151367  , 2661.48388672 , 2397.94360352 , 2579.23974609 , 2125.00439453 , 2119.68359375 , 1631.13183594 ,
                    1546.22912598  , 1129.73095703 , 1030.00244141 , 725.08172607  , 641.80291748  , 438.453125    , 335.5932312   ,
                    225.46757507   , 175.30892944  , 115.84326172  , 101.24755859  , 64.85153198   , 56.26319504   , 33.94129181]
                }
    },
    
    "multiplicity" : {
        "data" : {
            20: [ 1.5         , 3.          , 6.   , 6.         , 9.07692308  , 8.76923077 , 11.69230769 ,
                9.92307692    , 12.         , 9.   , 10.        , 7.34615385  , 8.         , 5.          ,
                5.            , 3.          , 3.   , 1.92307692 , 2.          , 1.         , 0.          ,
                0.            , 0.          , 0.   , 0.         , 0.          , 0.         , 0.        ] ,
            100: [ 2.90909091 , 4.          , 9.   , 10.        , 15.63636364 , 16.        , 22.         ,
                21.           , 26.         , 23.  , 27.        , 23.27272727 , 27.        , 22.         ,
                24.           , 17.54545455 , 19.  , 14.        , 14.09090909 , 9.         , 7.63636364  ,
                4.54545455    , 3.54545455  , 2.   , 2.         , 1.          , 1.         , 0.        ] ,
            300: [ 3.         , 5.          , 11.  , 13.2       , 20.         , 22.8       , 32.         ,
                31.           , 37.8        , 36.  , 42.        , 39.         , 43.8       , 40.4        ,
                44.           , 35.2        , 40.2 , 33.        , 35.         , 26.6       , 25.         ,
                17.8          , 15.         , 10.  , 10.        , 6.          , 5.         , 2.6]
            },
        
        "mc" : {
            20: [ 1.   , 2.   , 5.  , 5.   , 8.2  , 8.   , 11.  ,
                9.8    , 11.2 , 9.  , 10.  , 8.   , 8.   , 5.6  ,
                5.     , 3.   , 3.  , 2.   , 1.   , 1.   , 0.   ,
                0.     , 0.   , 0.  , 0.   , 0.   , 0.   , 0. ] ,
            100: [ 1.5 , 3.   , 7.  , 9.   , 14.  , 17.  , 21.  ,
                23.    , 26.  , 26. , 28.  , 25.  , 27.  , 23.  ,
                24.    , 19.  , 19. , 14.  , 13.6 , 9.   , 7.   ,
                5.     , 3.2  , 2.  , 2.   , 1.   , 1.   , 0. ] ,
            300: [ 2.  , 4.   , 9.  , 13.  , 19.  , 25.6 , 31.  ,
                36.    , 40.  , 43. , 46.2 , 45.  , 48.  , 44.  ,
                46.    , 41.  , 42. , 35.  , 35.  , 28.  , 25.  ,
                18.2   , 16.  , 11. , 10.  , 6.   , 5.   , 3. ]
        }
    }
}
