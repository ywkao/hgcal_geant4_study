#!/usr/bin/env python
import os, subprocess
import math
import copy
import array
import ROOT
ROOT.gROOT.SetBatch(True)

import MetaData as m

eos = "./eos"
def create_directory(dir_output):
    if not os.path.isdir(dir_output):
        subprocess.call("mkdir %s" % dir_output, shell=True)
        subprocess.call("cp -p %s/index.php %s" % (eos, dir_output), shell=True)

def annotate(rshift=0):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(43)
    latex.SetTextAlign(11)
    #latex.SetTextSize(24)
    latex.SetTextSize(20)
    latex.DrawLatex( 0.12, 0.912, "#bf{CMS} #it{Preliminary}" )
    latex.DrawLatex( 0.56+rshift, 0.912, "D86 Simulation with 1,000 events" )
    #latex.DrawLatex( 0.58+rshift, 0.912, "D86 Simulation with 1,000 events" )
    #latex.DrawLatex( 0.69+rshift, 0.912, "%s fb^{-1} (13 TeV)" % str(lumi["RunII"]) )

def draw_2D_ntuple(hname, v_hists, selection, color, is_first_plot=False):
    if is_first_plot: v_hists[0].Draw("r:z>>%s(2500, 300, 550, 3000, 0, 300)" % hname, selection)
    else:             v_hists[0].Draw("r:z>>%s(2500, 300, 550, 3000, 0, 300)" % hname, selection, "same")

    hnew = ROOT.gDirectory.Get(hname)
    hnew.SetTitle("")
    hnew.GetXaxis().SetTitle("Z [cm]")
    hnew.GetYaxis().SetTitle("Rxy [cm]")
    hnew.SetMarkerColor(color)
    if is_first_plot: hnew.Draw()
    else:             hnew.Draw("same")

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

def set_graph(gr, xtitle, ytitle, color):
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
        "total_MIP_coarse"            : "total_MIP",
        "total_MIP_fine"              : "total_MIP",
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
        "multiplicity_simhits_coarse" : "multiplicity_simhits",
        "multiplicity_simhits_fine"   : "multiplicity_simhits",
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
    lex = [0. for i in range(len(v_hists))]
    ley = [h.GetMeanError() for h in v_hists]
    gr = get_graph_from_list("Layer depth [ X_{0} ]", ytitles[varName], lx, ly, lex, ley, ROOT.kBlue, is_number_of_hits)
    return gr

def get_graph_from_list(xtitle, ytitle, lx, ly, lex, ley, color, normalize_to_unity = False):
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
        ex.append(lex[i])
        ey.append(ley[i])

    gr = ROOT.TGraphErrors(n, x, y, ex, ey)
    set_graph(gr, xtitle, ytitle, color)
    if normalize_to_unity: gr.SetMaximum(0.12)

    #print ">>> get_graph_from_list::total =", total
    #print ">>> n, total, ly, ley =", n, total, ly, ley
    return gr

#----------------------------------------------------------------------------------------------------

sigmaEoverE, error_sigmaEoverE = [], []
def reset_containers():
    global sigmaEoverE , error_sigmaEoverE
    sigmaEoverE, error_sigmaEoverE = [], []

fit_result = {}
def record_fit_result(myTags, func):
    print ">>> in record_fit_result...."
    global fit_result, sigmaEoverE, error_sigmaEoverE

    fit_const = func.GetParameter(0)
    fit_mean  = func.GetParameter(1)
    fit_sigma = func.GetParameter(2)
    fitError_const = func.GetParError(0)
    fitError_mean  = func.GetParError(1)
    fitError_sigma = func.GetParError(2)

    energyType = myTags[0] # MIP / SIM / ENE
    tag = myTags[1] # beam energy
    label = myTags[2] # odd / even / set0 / set1 / set2

    if not energyType in fit_result.keys():
        fit_result[energyType] = {}

    if not label in fit_result[energyType].keys():
        fit_result[energyType][label] = {}

    fit_result[energyType][label][tag] = {}
    fit_result[energyType][label][tag]["mean"] = fit_mean
    fit_result[energyType][label][tag]["error_mean"] = fitError_mean
    fit_result[energyType][label][tag]["sigma"] = fit_sigma
    fit_result[energyType][label][tag]["error_sigma"] = fitError_sigma

    ratio = fit_sigma/fit_mean
    uncertainty = ratio * math.sqrt( math.pow(fitError_mean/fit_mean, 2) + math.pow(fitError_sigma/fit_sigma, 2) )
    
    sigmaEoverE.append(ratio)
    error_sigmaEoverE.append(uncertainty)

    #return fit_mean, fit_sigma
    #print ">>> result:", fit_const, fit_mean, fit_sigma 
    #print ">>> fit error:", fitError_const, fitError_mean, fitError_sigma 
    
def set_stat_pad(stat, positions, color):
    if stat:
        #print ">>>>> check:", stat.GetName()
        stat.GetName()
        stat.SetX1NDC(positions[0])
        stat.SetY1NDC(positions[1])
        stat.SetX2NDC(positions[2])
        stat.SetY2NDC(positions[3])
        stat.SetTextSize(0.02)
        stat.SetTextColor(color)
        stat.SetLineColor(color)
    else:
        #print ">>>>> stat is null"
        return

def draw_and_fit_a_histogram(c1, h, myTags, xtitle, max_value, xRange, color, stat_position):
    maxbin = h.GetMaximumBin()
    bincenter = h.GetBinCenter(maxbin)
    #binwidth = h.GetBinWidth(maxbin)
    #fit_range_lower = bincenter - 6*binwidth
    #fit_range_upper = bincenter + 6*binwidth

    std_dev = h.GetStdDev()
    fit_range_lower = bincenter - 2*std_dev
    fit_range_upper = bincenter + 2*std_dev

    h.SetTitle("")
    h.SetMaximum(max_value*1.2)
    h.SetLineWidth(2)
    h.SetLineColor(color)
    h.GetXaxis().SetRangeUser(xRange[0], xRange[1])
    h.GetXaxis().SetTitleOffset(1.1)
    h.GetXaxis().SetTitle(xtitle)

    energyType = myTags[0] # MIP / SIM / ENE

    if energyType == "ENE":
        h.GetYaxis().SetTitle("Entries / 0.1 MeV")
    else:
        h.GetYaxis().SetTitle("Entries")

    h.Fit("gaus", "0", "", fit_range_lower, fit_range_upper)

    ## second fit
    #h.Draw()
    #func = h.GetListOfFunctions().FindObject("gaus")
    #fit_mean  = func.GetParameter(1)
    #fit_sigma = func.GetParameter(2)
    #fit_range_lower = fit_mean - 2*fit_sigma
    #fit_range_upper = fit_mean + 2*fit_sigma
    #h.Fit("gaus", "0", "", fit_range_lower, fit_range_upper)

    h.Draw()
    h.GetFunction("gaus").Draw("same")

    c1.Update()
    lof = h.GetListOfFunctions()
    record_fit_result( myTags, lof.FindObject("gaus") ) # record in pu.sigmaEoverE and pu.fit_result
    set_stat_pad( lof.FindObject("stats"), stat_position, color )

def get_strings_for_simple_plot(energyType, selection):
    if energyType == "ENE":
        xtitle = "Corrected energy [MeV]"
        outputName = "h_Edep_corrected_energy_%s_" % selection

        if selection == "set1_set2":
            hname_odd = "total_ENE_set1"
            hname_even = "total_ENE_set2"
            labels = ["set1", "set2"]

        if selection == "set0":
            hname_odd = "total_ENE_set0"
            hname_even = "total_ENE_set0"
            labels = ["set0", "set0"]

    if energyType == "MIP":
        xtitle = "Deposited energy [MIP]"
        outputName = "h_Edep_MIP_%s_" % selection

        if selection == "odd_even":
            hname_odd = "total_MIP_odd"
            hname_even = "total_MIP_even"
            labels = ["odd", "even"]

        if selection == "set1_set2":
            hname_odd = "total_MIP_set1"
            hname_even = "total_MIP_set2"
            labels = ["set1", "set2"]

        if selection == "set0":
            hname_odd = "total_MIP_set0"
            hname_even = "total_MIP_set0"
            labels = ["set0", "set0"]

    if energyType == "SIM":
        xtitle = "Deposited energy [MeV]"
        outputName = "h_Edep_SIM_%s_" % selection

        if selection == "odd_even":
            hname_odd = "total_SIM_odd"
            hname_even = "total_SIM_even"
            labels = ["odd", "even"]

        if selection == "set1_set2":
            hname_odd = "total_SIM_set1"
            hname_even = "total_SIM_set2"
            labels = ["set1", "set2"]

        if selection == "set0":
            hname_odd = "total_SIM_set0"
            hname_even = "total_SIM_set0"
            labels = ["set0", "set0"]

    return xtitle, outputName, hname_odd, hname_even, labels

#----------------------------------------------------------------------------------------------------

