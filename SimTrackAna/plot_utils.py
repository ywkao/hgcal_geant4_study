#!/usr/bin/env python
import math
import array
import ROOT
ROOT.gROOT.SetBatch(True)

def load_histograms(fin, varName, processes):
    v_hists = []
    for p in processes:
        histName = "prodEE_DigiSim/" + varName + "_" + str(p)
        h = fin.Get(histName)
        v_hists.append(h)

    return v_hists

def set_graph(gr, ytitle, xtitle, color):
    gr.SetTitle("")
    gr.SetLineColor(color)
    gr.SetLineWidth(2)
    gr.SetMarkerStyle(20)
    gr.SetMarkerColor(color)
    gr.GetXaxis().SetTitle(xtitle)
    gr.GetYaxis().SetTitle(ytitle)
    gr.GetXaxis().SetTitleOffset(1.2)

sub_directory = {
    "ADC_120mum_layer"         : "ADC",
    "ADC_200mum_layer"         : "ADC",
    "ADC_300mum_layer"         : "ADC",
    "MIP_120mum_layer"         : "MIP",
    "MIP_200mum_layer"         : "MIP",
    "MIP_300mum_layer"         : "MIP",
    "SIM_120mum_layer"         : "SIM",
    "SIM_200mum_layer"         : "SIM",
    "SIM_300mum_layer"         : "SIM",
    "ADC_SimhitE_120mum_layer" : "ADC_SimhitE",
    "ADC_SimhitE_200mum_layer" : "ADC_SimhitE",
    "ADC_SimhitE_300mum_layer" : "ADC_SimhitE",
    "ADC_MIP_120mum_layer"     : "ADC_MIP",
    "ADC_MIP_200mum_layer"     : "ADC_MIP",
    "ADC_MIP_300mum_layer"     : "ADC_MIP",
    "MIP_SimhitE_120mum_layer" : "MIP_SimhitE",
    "MIP_SimhitE_200mum_layer" : "MIP_SimhitE",
    "MIP_SimhitE_300mum_layer" : "MIP_SimhitE",
}

ytitles = sub_directory

x_D86 = [0.564,1.567,2.547,3.549,4.528,5.531,6.509,7.512,8.49,9.493,10.472,11.474,12.453,13.455,14.434,15.437,16.415,17.418,18.975,19.978,21.536,22.538,24.096,25.099,26.656,27.659]

def get_graph(varName, v_hists):
    n = 26 
    x, ex = array.array('d'), array.array('d')
    y, ey = array.array('d'), array.array('d')

    for i, h in enumerate(v_hists):
        x.append(x_D86[i])
        y.append(h.GetMean())
        ex.append(0.)
        ey.append(h.GetStdDev())

        gr = ROOT.TGraphErrors(n, x, y, ex, ey)
        set_graph(gr, ytitles[varName], "Layer depth [ X_{0} ]", ROOT.kBlue)

    return gr
