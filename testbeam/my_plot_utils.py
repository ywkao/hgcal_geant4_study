#!/usr/bin/env python
import math
import copy
import array
import ROOT
ROOT.gROOT.SetBatch(True)

def load_single_histogram(fin, varName):
    v_hists = []
    histName = varName
    h = fin.Get(histName)
    v_hists.append(h)

    return v_hists

def load_histograms(fin, varName, processes):
    v_hists = []
    for p in processes:
        histName = varName + "_" + str(p)
        h = fin.Get(histName)
        v_hists.append(h)

    return v_hists

def set_graph(gr, ytitle, xtitle, color):
    gr.SetTitle("")
    gr.SetLineColor(color)
    gr.SetLineWidth(2)
    gr.SetMarkerStyle(20)
    gr.SetMarkerColor(color)
    gr.SetMinimum(0)
    gr.GetXaxis().SetTitle(xtitle)
    gr.GetYaxis().SetTitle(ytitle)
    gr.GetXaxis().SetTitleOffset(1.2)

sub_directory = {
    "hEta"                        : "EtaPhi",
    "hPhi"                        : "EtaPhi",

    "ADC_total_layer"             : "ADC",
    "ADC_120mum_layer"            : "ADC",
    "ADC_200mum_layer"            : "ADC",
    "ADC_300mum_layer"            : "ADC",
    "MIP_total_layer"             : "MIP",
    "MIP_120mum_layer"            : "MIP",
    "MIP_200mum_layer"            : "MIP",
    "MIP_300mum_layer"            : "MIP",
    "SIM_total_layer"             : "SIM",
    "SIM_120mum_layer"            : "SIM",
    "SIM_200mum_layer"            : "SIM",
    "SIM_300mum_layer"            : "SIM",
    "ADC_SimhitE_total_layer"     : "ADC_SimhitE",
    "ADC_SimhitE_120mum_layer"    : "ADC_SimhitE",
    "ADC_SimhitE_200mum_layer"    : "ADC_SimhitE",
    "ADC_SimhitE_300mum_layer"    : "ADC_SimhitE",
    "ADC_MIP_total_layer"         : "ADC_MIP",
    "ADC_MIP_120mum_layer"        : "ADC_MIP",
    "ADC_MIP_200mum_layer"        : "ADC_MIP",
    "ADC_MIP_300mum_layer"        : "ADC_MIP",
    "MIP_SimhitE_total_layer"     : "MIP_SimhitE",
    "MIP_SimhitE_120mum_layer"    : "MIP_SimhitE",
    "MIP_SimhitE_200mum_layer"    : "MIP_SimhitE",
    "MIP_SimhitE_300mum_layer"    : "MIP_SimhitE",

    "total_ADC_total"             : "total_ADC",
    "total_ADC_120mum"            : "total_ADC",
    "total_ADC_200mum"            : "total_ADC",
    "total_ADC_300mum"            : "total_ADC",
    "total_MIP_total"             : "total_MIP",
    "total_MIP_120mum"            : "total_MIP",
    "total_MIP_200mum"            : "total_MIP",
    "total_MIP_300mum"            : "total_MIP",
    "total_SIM_total"             : "total_SIM",
    "total_SIM_120mum"            : "total_SIM",
    "total_SIM_200mum"            : "total_SIM",
    "total_SIM_300mum"            : "total_SIM",
    "multiplicity_digis_total"    : "multiplicity_digis",
    "multiplicity_digis_120mum"   : "multiplicity_digis",
    "multiplicity_digis_200mum"   : "multiplicity_digis",
    "multiplicity_digis_300mum"   : "multiplicity_digis",
    "multiplicity_simhits_300mum" : "multiplicity_simhits",
    "multiplicity_simhits_total"  : "multiplicity_simhits",
    "multiplicity_simhits_120mum" : "multiplicity_simhits",
    "multiplicity_simhits_200mum" : "multiplicity_simhits",

    "rechit_amplitude_total" : "amplitude",
    "rechit_amplitudeHigh_total" : "amplitudeHigh",
    "rechit_amplitudeLow_total"  : "amplitudeLow",
}

ytitles = copy.deepcopy( sub_directory )
ytitles["hEta"] = "Eta"
ytitles["hPhi"] = "Phi"

x_D86 = [0.564,1.567,2.547,3.549,4.528,5.531,6.509,7.512,8.49,9.493,10.472,11.474,12.453,13.455,14.434,15.437,16.415,17.418,18.975,19.978,21.536,22.538,24.096,25.099,26.656,27.659]

def get_graph(varName, v_hists):
    #n = 26 
    n = 28 
    x, ex = array.array('d'), array.array('d')
    y, ey = array.array('d'), array.array('d')

    for i, h in enumerate(v_hists):
        #x.append(x_D86[i])
        x.append(float(i+1))
        y.append(h.GetMean())
        ex.append(0.)
        #ey.append(h.GetStdDev())
        ey.append(h.GetMeanError())

        gr = ROOT.TGraphErrors(n, x, y, ex, ey)
        #set_graph(gr, ytitles[varName], "Layer depth [ X_{0} ]", ROOT.kBlue)
        set_graph(gr, ytitles[varName], "Layer", ROOT.kBlue)

    return gr
