#!/usr/bin/env python
import math
import copy
import array
import ROOT
ROOT.gROOT.SetBatch(True)

import MetaData as m

def load_single_histogram(fin, varName):
    histName = "prodEE_DigiSim/" + varName
    h = fin.Get(histName)
    return h

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
    gr.SetMinimum(0)
    gr.GetXaxis().SetTitle(xtitle)
    gr.GetYaxis().SetTitle(ytitle)
    gr.GetXaxis().SetTitleOffset(1.2)

sub_directory = {
        "hEta"                        : "EtaPhi",
        "hPhi"                        : "EtaPhi",
        "nt_hit_position"             : "",

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
        "total_MIP_coarse"             : "total_MIP",
        "total_MIP_fine"             : "total_MIP",
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
        "multiplicity_simhits_coarse"  : "multiplicity_simhits",
        "multiplicity_simhits_fine"  : "multiplicity_simhits",
        "multiplicity_simhits_120mum" : "multiplicity_simhits",
        "multiplicity_simhits_200mum" : "multiplicity_simhits",
        }

ytitles = copy.deepcopy( sub_directory )
ytitles["hEta"] = "Eta"
ytitles["hPhi"] = "Phi"
ytitles["nt_hit_position"] = "Rxy [cm]"

x_D86 = [0.564,1.567,2.547,3.549,4.528,5.531,6.509,7.512,8.49,9.493,10.472,11.474,12.453,13.455,14.434,15.437,16.415,17.418,18.975,19.978,21.536,22.538,24.096,25.099,26.656,27.659]

def get_graph(varName, v_hists, is_number_of_hits = False):
    lx = x_D86
    ly = [h.GetMean() for h in v_hists]
    ley = [h.GetMeanError() for h in v_hists]
    gr = get_graph_from_list(varName, lx, ly, ley, is_number_of_hits)
    return gr

def get_graph_from_list(varName, lx, ly, ley, normalize_to_unity = False):
    n = len(lx) 
    x, ex = array.array('d'), array.array('d')
    y, ey = array.array('d'), array.array('d')

    total = sum(ly) 
    scale = 1./total if normalize_to_unity else 1.
    ly = [ele*scale for ele in ly]
    ley = [ele*scale for ele in ley]

    for i, ele in enumerate(lx):
        x.append(ele)
        y.append(ly[i])
        ex.append(0.)
        ey.append(ley[i])

    gr = ROOT.TGraphErrors(n, x, y, ex, ey)
    set_graph(gr, ytitles[varName], "Layer depth [ X_{0} ]", ROOT.kBlue)
    if normalize_to_unity: gr.SetMaximum(0.12)

    #print ">>> get_graph_from_list::total =", total
    #print ">>> n, total, ly, ley =", n, total, ly, ley
    return gr


def record_fit_result(func):
    #global d_fit_const, d_fit_mean, d_fit_sigma

    fit_const = func.GetParameter(0)
    fit_mean  = func.GetParameter(1)
    fit_sigma = func.GetParameter(2)
    fitError_const = func.GetParError(0)
    fitError_mean  = func.GetParError(1)
    fitError_sigma = func.GetParError(2)
    print ">>> result:", fit_const, fit_mean, fit_sigma 
    print ">>> fit error:", fitError_const, fitError_mean, fitError_sigma 
    
    #d_fit_const["central"].append( func.GetParameter(0) )
    #d_fit_mean ["central"].append( func.GetParameter(1) )
    #d_fit_sigma["central"].append( func.GetParameter(2) )
    #d_fit_const[ "error" ].append( func.GetParError(0)  )
    #d_fit_mean [ "error" ].append( func.GetParError(1)  )
    #d_fit_sigma[ "error" ].append( func.GetParError(2)  )
