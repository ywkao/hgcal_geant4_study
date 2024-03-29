#!/usr/bin/env python2
import math
import ROOT
ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetOptStat("emr")
ROOT.gStyle.SetOptStat(2210)
ROOT.gStyle.SetOptFit(1111)

my_stat_pos = [0.89, 0.87, 0.24, 0.08]
ROOT.gStyle.SetStatX(my_stat_pos[0])
ROOT.gStyle.SetStatY(my_stat_pos[1])
ROOT.gStyle.SetStatW(my_stat_pos[2])
ROOT.gStyle.SetStatH(my_stat_pos[3])
ROOT.gStyle.SetTextSize(1.2)

c1 = ROOT.TCanvas("c1", "", 800, 600)
c1.SetGrid()
c1.SetTicks(1,1)
c1.SetLeftMargin(0.12)
c1.SetRightMargin(0.08)

c2 = ROOT.TCanvas("c2", "", 5600, 2400)
c2.SetLeftMargin(0.)
c2.SetRightMargin(0.)
c2.Divide(7,4)

c3 = ROOT.TCanvas("c3", "", 1200, 600)
c3.SetGrid()
c3.SetTicks(1,1)
c3.SetLeftMargin(0.12)
c3.SetRightMargin(0.08)

import toolbox.plot_utils as pu
import toolbox.MetaData as m

#--------------------------------------------------
# parameters to set
#--------------------------------------------------
tags = []
label = {}
myRootfiles = []
flag_add_reference = False
fit_constraints = {}
colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kBlue-7, ROOT.kMagenta, ROOT.kRed-7]

#--------------------------------------------------
# plotting functions
#--------------------------------------------------
def make_plot(varName, dir_output="", bool_make_logitudinal_profile=False):
    global myRootfiles, flag_add_reference
    is_number_of_hits = "multiplicity" in varName
    
    #++++++++++++++++++++++++++++++
    # Initiate
    #++++++++++++++++++++++++++++++
    bool_ntuple = "nt_" in varName
    bool_this_is_eta_phi = "Eta" in varName or "Phi" in varName
    bool_single_figures = bool_ntuple or bool_this_is_eta_phi
    processes = [str(i) for i in range(1,27)]
    if len(dir_output) == 0:
        dir_output = m.specified_directory + "/" + pu.sub_directory[varName]

    #++++++++++++++++++++++++++++++
    # Load histograms
    #++++++++++++++++++++++++++++++
    vf = []
    v_v_hists = []
    for rootfile in myRootfiles:
        #print ">>> rootfile", rootfile
        fin = ROOT.TFile.Open(rootfile, "R")
        vf.append(fin)

        v_hists = []
        if bool_single_figures:
            v_hists.append( pu.load_single_histogram(fin, varName) )
            v_v_hists.append(v_hists)
        else:
            v_hists = pu.load_histograms(fin, varName, processes)
            v_v_hists.append(v_hists)

    #++++++++++++++++++++++++++++++
    # Eta or Phi or ntuple
    #++++++++++++++++++++++++++++++
    if bool_single_figures:
        c1.cd()
        for i, v_hists in enumerate(v_v_hists):
            if bool_this_is_eta_phi:
                if "Eta" in varName: v_hists[0].GetXaxis().SetRangeUser(1., 3.)
                if "Phi" in varName: v_hists[0].GetXaxis().SetRangeUser(-3., 3.)
                maximum = v_hists[0].GetMaximum()*1.3
                v_hists[0].SetMaximum(maximum)
                v_hists[0].SetMinimum(0.)
                v_hists[0].Draw()
            elif bool_ntuple:
                pu.draw_2D_ntuple("hnew"+"_"+tags[i]+"_w300", v_hists, "is_Silicon_w300>0", ROOT.kMagenta, True)
                pu.draw_2D_ntuple("hnew"+"_"+tags[i]+"_w200", v_hists, "is_Silicon_w200>0", ROOT.kGreen)
                pu.draw_2D_ntuple("hnew"+"_"+tags[i]+"_w120", v_hists, "is_Silicon_w120>0", ROOT.kRed)
            else: continue # no matched situation

            output = dir_output + "/" + varName + "_" + tags[i]
            c1.SaveAs(output + ".png")
            #c1.SaveAs(output + ".pdf")

    #++++++++++++++++++++++++++++++
    # Longitdinal profile
    #++++++++++++++++++++++++++++++
    if bool_make_logitudinal_profile:
        v_gr = []

        #------------------------------
        # loop over energy
        #------------------------------
        for i, v_hists in enumerate(v_v_hists):
            #gr = pu.get_graph(varName, v_hists, is_number_of_hits)
            gr = pu.get_graph(varName, v_hists, False)
            gr.SetLineStyle(2)
            if "efficiency" in varName:
                gr.SetMaximum(1.)

            if not flag_add_reference:
                gr.SetLineColor(colors[i])
                gr.SetMarkerColor(colors[i])
            else:
                gr.SetMarkerColor(ROOT.kBlue)
                gr.SetLineColor(ROOT.kBlue)
            v_gr.append(gr)

            # 26 plots in one
            c2.cd()
            for j, h in enumerate(v_hists):
                c2.cd(j+1)
                c2.SetGrid()
                c2.SetTicks(1,1)
                h.Draw()

            output = dir_output + "/" + varName + "_all_plots_" + tags[i]
            c2.SaveAs(output + ".png")
            c2.SaveAs(output + ".pdf")

            continue

            # individual plots
            c1.cd()
            pu.create_directory(dir_output)
            for j, h in enumerate(v_hists):
                h.Draw()
                output = dir_output + "/" + varName + "_" + processes[j] + "_" + tags[i]
                c1.SaveAs(output + ".png")
                c1.SaveAs(output + ".pdf")

        #------------------------------
        # ref from test beam
        #------------------------------
        xt, yt, tbr = "Layer depth [ X_{0} ]", pu.ytitles, m.test_beam_result
        v_gr_ref_mc, v_gr_ref_data, ref_ex, ref_ey = [], [], [0.]*28, [0.]*28
        for ene in [300, 100, 20]:
            ref_tag = "multiplicity" if 'multiplicity' in varName else "mips"
            gr = pu.get_graph_from_list(xt, yt[varName], tbr["layer_depth"], tbr[ref_tag]["mc"][ene], ref_ex, ref_ey, ROOT.kRed, is_number_of_hits)
            v_gr_ref_mc.append(gr)

            gr = pu.get_graph_from_list(xt, yt[varName], tbr["layer_depth"], tbr[ref_tag]["data"][ene], ref_ex, ref_ey, ROOT.kBlack, is_number_of_hits)
            v_gr_ref_data.append(gr)

        #------------------------------
        # logitudinal profile
        #------------------------------
        c1.cd()
        #legend = ROOT.TLegend(0.67, 0.65, 0.87, 0.85)
        #legend = ROOT.TLegend(0.65, 0.65, 0.85, 0.85)

        if not flag_add_reference:
            legend = ROOT.TLegend(0.70, 0.65, 0.85, 0.85)
            legend.SetLineColor(0)
            legend.SetTextSize(0.04)

            for i, gr in enumerate(v_gr):
                if i==0: gr.Draw("alp")
                else:    gr.Draw('lp;same')
                legend.AddEntry(gr, label[tags[i]], "lp")

                if i==0 and 'MIP' in varName:
                    gr.SetMaximum(1060.)
                    gr.SetMaximum(3000.)

                if i+1 == len(v_gr):
                    pu.annotate()
                    legend.Draw("same")

                    output = dir_output + "/" + varName
                    c1.SaveAs(output + ".png")
                    c1.SaveAs(output + ".pdf")

        else:
            legend = ROOT.TLegend(0.60, 0.65, 0.85, 0.85)
            legend.SetLineColor(0)
            legend.SetTextSize(0.025)
            for i, gr in enumerate(v_gr):
                gr.Draw("alp")
                v_gr_ref_mc[i].Draw('lp;same')
                v_gr_ref_data[i].Draw('lp;same')

                legend.Clear()
                legend.AddEntry(v_gr_ref_data[i], "2018 TB Data", "lp")
                legend.AddEntry(v_gr_ref_mc[i], "2018 TB MC", "lp")
                legend.AddEntry(gr, "D86 geometry", "lp")

                pu.annotate()
                legend.Draw("same")

                output = dir_output + "/" + varName + "_" + tags[i]
                c1.SaveAs(output + ".png")
                c1.SaveAs(output + ".pdf")

#----------------------------------------------------------------------------------------------------

def make_simple_plot(energyType, dir_output, selection):
    global myRootfiles, fit_constraints
    
    xtitle, outputName, hname_odd, hname_even, labels = pu.get_strings_for_simple_plot(energyType, selection)

    #++++++++++++++++++++++++++++++
    # Load histograms
    #++++++++++++++++++++++++++++++
    vf, v_v_hists = [], []
    for rootfile in myRootfiles:
        #print ">>> rootfile", rootfile
        fin = ROOT.TFile.Open(rootfile, "R")
        vf.append(fin)

        v_hists = []
        v_hists.append( pu.load_single_histogram(fin, hname_odd) )
        v_hists.append( pu.load_single_histogram(fin, hname_even) )
        v_v_hists.append(v_hists)

    c1.cd()
    for i, v_hists in enumerate(v_v_hists):
        #print "\n----------------------------------------------------------------------------------------------------\n"

        if energyType == "ENE" and "set" in selection: # and m.type_resolution == "resolution_unclustered":
            # type_resolution: resolution_clustered, resolution_unclustered
            # selection: set0, set1_set2, odd_even, etc.
            # tags = E300, E100, E20, etc.
            rebin_factor = m.rebin_factor[m.type_resolution][selection][tags[i]]
            v_hists[0].Rebin(rebin_factor)
            v_hists[1].Rebin(rebin_factor)

        pu.reset_containers()
        max_values = []
        max_values.append(v_hists[0].GetMaximum())
        max_values.append(v_hists[1].GetMaximum())
        max_value = max(max_values)

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextFont(43)
        latex.SetTextAlign(11)
        latex.SetTextSize(24)

        if not selection == "set0":
            #--------------------------------------------------
            # Edep odd/even layers
            #--------------------------------------------------
            xRanges = fit_constraints["set1set2"][energyType]["xRanges"]

            pu.draw_and_fit_a_histogram(c1, v_hists[0], [energyType, tags[i], labels[0]], xtitle, max_value, xRanges[i], ROOT.kBlue   , [0.60, 0.66, 0.88, 0.86] )
            pu.draw_and_fit_a_histogram(c1, v_hists[1], [energyType, tags[i], labels[1]], xtitle, max_value, xRanges[i], ROOT.kGreen+3, [0.60, 0.42, 0.88, 0.62] )

            v_hists[0].Draw("same")
            v_hists[0].GetFunction("gaus").SetLineStyle(1)
            v_hists[0].GetFunction("gaus").SetLineColorAlpha(ROOT.kRed, 0.40)
            v_hists[0].GetFunction("gaus").Draw("same")

            #--------------------------------------------------
            # result
            #--------------------------------------------------
            latex.SetTextColor(ROOT.kBlue)
            latex.DrawLatex( 0.47, 0.30, "#sigma#left(E_{%s}#right) / #bar{E}_{%s} = %.4f #pm %.4f" % ("def.", "def.", pu.sigmaEoverE[0], pu.error_sigmaEoverE[0]) )
            latex.SetTextColor(ROOT.kGreen+3)
            latex.DrawLatex( 0.47, 0.20, "#sigma#left(E_{%s}#right) / #bar{E}_{%s} = %.4f #pm %.4f" % ("alt.", "alt.", pu.sigmaEoverE[1], pu.error_sigmaEoverE[1]) )

        else:
            xRanges = fit_constraints["set0"][energyType]["xRanges"]

            pu.draw_and_fit_a_histogram(c1, v_hists[0], [energyType, tags[i], labels[0]] , xtitle, max_value, xRanges[i], ROOT.kBlue  , [0.60, 0.66, 0.88, 0.86] )

            latex.SetTextColor(ROOT.kBlue)
            latex.DrawLatex( 0.47, 0.30, "#sigma#left(E_{%s}#right) / #bar{E}_{%s} = %.4f #pm %.4f" % ("tot.", "tot.", pu.sigmaEoverE[0], pu.error_sigmaEoverE[0]) )

        c1.Update()
        pu.annotate()
        output = dir_output + "/" + outputName + tags[i]
        c1.SaveAs(output + ".png")
        c1.SaveAs(output + ".pdf")

#----------------------------------------------------------------------------------------------------

def run_linear_fit(dir_output, fit_result, label, dx, dy):
    x_range = m.linear_fit_parameter[label]["linear_fit_xrange"]
    y_range = m.draw_options_for_run_summary[m.type_resolution]["linear_fit_yrange"]


    Energy = ["E20", "E60", "E100", "E175", "E225", "E300"]
    lx  = [ dx[ene]["mean"]  for ene in Energy ]
    ly  = [ dy[ene]["mean"]  for ene in Energy ]
    lex = [ dx[ene]["sigma"] for ene in Energy ]
    ley = [ dy[ene]["sigma"] for ene in Energy ]

    gr = pu.get_graph_from_list("Energy (MIPs)", "Calibrated energy (GeV)", lx, ly, lex, ley, ROOT.kBlack)

    c1.cd()
    c1.Clear()
    gr.Draw("ap")
    gr.GetXaxis().SetLimits(x_range[0], x_range[1]) # along X
    gr.GetYaxis().SetRangeUser(y_range[0], y_range[1])

    f1 = ROOT.TF1('f1', "[0] + [1]*x", x_range[0], x_range[1])
    gr.Fit(f1, "", "", x_range[0], x_range[1])

    func = gr.GetListOfFunctions().FindObject("f1")
    fit_intercept       = func.GetParameter(0)
    fit_slope           = func.GetParameter(1)
    fit_error_intercept = func.GetParError(0)
    fit_error_slope     = func.GetParError(1)

    if not "linear_fit" in fit_result.keys(): fit_result["linear_fit"] = {}
    fit_result["linear_fit"][label] = {}
    fit_result["linear_fit"][label]["slope"] = fit_slope
    fit_result["linear_fit"][label]["intercept"] = fit_intercept
    fit_result["linear_fit"][label]["error_slope"] = fit_error_slope
    fit_result["linear_fit"][label]["error_intercept"] = fit_error_intercept

    print ">>> fit_intercept =", fit_intercept
    print ">>> fit_slope =", fit_slope
    print ">>> fit_error_intercept =", fit_error_intercept
    print ">>> fit_error_slope =", fit_error_slope

    my_stat_pos = [0.42, 0.87, 0.15, 0.15]
    ROOT.gStyle.SetStatX(my_stat_pos[0])
    ROOT.gStyle.SetStatY(my_stat_pos[1])
    ROOT.gStyle.SetStatW(my_stat_pos[2])
    ROOT.gStyle.SetStatH(my_stat_pos[3])

    pu.annotate()
    output = dir_output + "/correction_generatedShowerEnergy_MIPs_" + label
    c1.SaveAs(output + ".png")
    c1.SaveAs(output + ".pdf")

#----------------------------------------------------------------------------------------------------

def perform_consistency_check_v2(title, energy, unc_enegy, sigma, unc_sigma):
    e1 = energy[0]
    e2 = energy[1]
    s1 = sigma[0]
    s2 = sigma[1]

    errE1 = unc_enegy[0]
    errE2 = unc_enegy[1]
    errS1 = unc_sigma[0]
    errS2 = unc_sigma[1]

    ratio = e2/e1
    res_ratio = s2/s1
    test = math.sqrt(1./ratio) * res_ratio
    uncertainty = test * math.sqrt( 0.25*math.pow(errE1/e1,2) + 0.25*math.pow(errE2/e2,2) + math.pow(errS1/s1,2) + math.pow(errS2/s2,2) )

    print ">>> consitency check %5s: test = %.2f #pm %.2f (e2/e1 = %.2f)" % (title, test, uncertainty, ratio)

def perform_consistency_check_v1(title, energy, unc_enegy, resolution, unc_resolution):
    e1 = energy[0]
    e2 = energy[1]
    r1 = resolution[0]
    r2 = resolution[1]

    errE1 = unc_enegy[0]
    errE2 = unc_enegy[1]
    errR1 = unc_resolution[0]
    errR2 = unc_resolution[1]

    ratio = e2/e1
    res_ratio = r2/r1
    test = math.sqrt(ratio) * res_ratio
    uncertainty = test * math.sqrt( 0.25*math.pow(errE1/e1,2) + 0.25*math.pow(errE2/e2,2) + math.pow(errR1/r1,2) + math.pow(errR2/r2,2) )

    print ">>> consitency check %5s: test = %.2f #pm %.2f (e2/e1 = %.2f)" % (title, test, uncertainty, ratio)

def run_summary(title, dy1, dy2):
    dir_output = m.specified_directory
    options = m.draw_options_for_run_summary[title]

    c3.cd()
    c3.Clear()

    # title and range
    ytitle         = options["ytitle"]
    yrange         = options["yrange"]
    leg_pos        = options["leg_pos"]
    draw_goodness  = options["draw_goodness"]
    draw_lower_pad = options["draw_lower_pad"]
    reference_line = options["reference_line"]
    leg_option     = options["leg_option"]
    c3.SetLogy(options["useLog"])

    print title, "yrange = ", yrange

    # list of x and y
    Energy = ["E20", "E60", "E100", "E175", "E225", "E300"]
    lx  = [ float(ene.split('E')[1]) for ene in Energy ]
    lex = [0.]*6

    ly1, ly2, ley1, ley2 = [], [], [], []
    for ene in Energy:
        if title == "pvalue" :
            ly1.append( dy1[ene]["pvalue"] )
            ly2.append( dy2[ene]["pvalue"] )
            ley1.append(0.)
            ley2.append(0.)

        if title == "chi2ndf":
            if not dy1[ene]["ndf"]>0 or not dy2[ene]["ndf"]>0:
                ly1.append( dy1[ene]["chi2"] )
                ly2.append( dy2[ene]["chi2"] )
                print ">>>>>>>>>> [WARNING] NDF is zero."
            else:
                ly1.append( dy1[ene]["chi2"]/dy1[ene]["ndf"] )
                ly2.append( dy2[ene]["chi2"]/dy2[ene]["ndf"] )

            ley1.append(0.)
            ley2.append(0.)

        if title == "resolution_clustered" or title == "resolution_unclustered":
            fit_mean = dy1[ene]["mean"]
            fit_sigma = dy1[ene]["sigma"]
            fitError_mean  = dy1[ene]["error_mean"]
            fitError_sigma = dy1[ene]["error_sigma"]
            resolution = fit_sigma/fit_mean
            uncertainty = resolution * math.sqrt( math.pow(fitError_mean/fit_mean, 2) + math.pow(fitError_sigma/fit_sigma, 2) )
            ly1.append(resolution)
            ley1.append(uncertainty)
            m.resolution["set1"][m.energy_type][ene] = {}
            m.resolution["set1"][m.energy_type][ene]["mean"] = resolution
            m.resolution["set1"][m.energy_type][ene]["error"] = uncertainty

            fit_mean = dy2[ene]["mean"]
            fit_sigma = dy2[ene]["sigma"]
            fitError_mean  = dy2[ene]["error_mean"]
            fitError_sigma = dy2[ene]["error_sigma"]
            resolution = fit_sigma/fit_mean
            uncertainty = resolution * math.sqrt( math.pow(fitError_mean/fit_mean, 2) + math.pow(fitError_sigma/fit_sigma, 2) )
            ly2.append(resolution)
            ley2.append(uncertainty)
            m.resolution["set2"][m.energy_type][ene] = {}
            m.resolution["set2"][m.energy_type][ene]["mean"] = resolution
            m.resolution["set2"][m.energy_type][ene]["error"] = uncertainty

            # consistency check
            #if "MIP" in m.energy_type:

            my_energy = [dy1[ene]["mean"], dy2[ene]["mean"]]
            my_energy_unc = [dy1[ene]["error_mean"], dy2[ene]["error_mean"]]
            my_res = [ly1[-1], ly2[-1]]
            my_res_unc = [ley1[-1], ley2[-1]]
            my_sigma = [dy1[ene]["sigma"], dy2[ene]["sigma"]]
            my_sigma_unc = [dy1[ene]["error_sigma"], dy2[ene]["error_sigma"]]

            #perform_consistency_check_v1( ene, my_energy, my_energy_unc, my_res, my_res_unc )
            perform_consistency_check_v2( ene, my_energy, my_energy_unc, my_sigma, my_sigma_unc )

        if "bias" in title:
            if m.energy_type == "set1_set2_MIPs":
                ref = ( dy1[ene]["mean"] + dy2[ene]["mean"] ) / 2.
                rel_e1 = dy1[ene]["mean"] / ref
                rel_e2 = dy2[ene]["mean"] / ref
                unc_rel_e1 = dy1[ene]["error_mean"] / ref
                unc_rel_e2 = dy2[ene]["error_mean"] / ref

                bias_e1 = rel_e1 - 1.
                bias_e2 = rel_e2 - 1.

                ly1.append(bias_e1)
                ley1.append(unc_rel_e1)

                ly2.append(bias_e2)
                ley2.append(unc_rel_e2)

            if m.energy_type == "set1_set2_MeV":
                rel_e1 = dy1[ene]["mean"] / float(ene.split("E")[1])
                rel_e2 = dy2[ene]["mean"] / float(ene.split("E")[1])
                unc_rel_e1 = dy1[ene]["error_mean"] / float(ene.split("E")[1])
                unc_rel_e2 = dy2[ene]["error_mean"] / float(ene.split("E")[1])

                bias_e1 = rel_e1 - 1.
                bias_e2 = rel_e2 - 1.

                ly1.append(bias_e1)
                ley1.append(unc_rel_e1)

                ly2.append(bias_e2)
                ley2.append(unc_rel_e2)

        #if title == "linearity":
        #    mip = dy1[ene]["mean"] / fit_result["linear_fit"]["slope"]

        if title == "changes_in_resolution":
            res_mev = dy1["set1_set2_MeV"][ene]["mean"]
            res_mip = dy1["set1_set2_MIPs"][ene]["mean"]
            err_mev = dy1["set1_set2_MeV"][ene]["error"]
            err_mip = dy1["set1_set2_MIPs"][ene]["error"]
            difference = (res_mev - res_mip) / res_mip
            ratio = res_mev / res_mip
            uncertainty = ratio * math.sqrt( math.pow(err_mev/res_mev,2) + math.pow(err_mip/res_mip, 2) )
            ly1.append(difference)
            ley1.append(uncertainty)
            #print ">>>>> diff in set1 resolution: %.2f +/- %.2f (res_mev = %.5f, res_mip = %.5f)" % (difference, uncertainty, res_mev, res_mip)

            res_mev = dy2["set1_set2_MeV"][ene]["mean"]
            res_mip = dy2["set1_set2_MIPs"][ene]["mean"]
            err_mev = dy2["set1_set2_MeV"][ene]["error"]
            err_mip = dy2["set1_set2_MIPs"][ene]["error"]
            difference = (res_mev - res_mip) / res_mip
            ratio = res_mev / res_mip
            uncertainty = ratio * math.sqrt( math.pow(err_mev/res_mev,2) + math.pow(err_mip/res_mip, 2) )
            ly2.append(difference)
            ley2.append(uncertainty)
            #print ">>>>> diff in set2 resolution: %.2f +/- %.2f (res_mev = %.5f, res_mip = %.5f)" % (difference, uncertainty, res_mev, res_mip)

    # graphs
    gr1 = pu.get_graph_from_list("Positron energy (GeV)", ytitle, lx, ly1, lex, ley1,  ROOT.kBlue)
    gr2 = pu.get_graph_from_list("Positron energy (GeV)", ytitle, lx, ly2, lex, ley2,  ROOT.kGreen+3)

    #print gr1, gr2
    #print ly1, ly2
    #print ley1, ley2

    gr1.SetLineWidth(2)
    gr1.SetMarkerStyle(20)
    gr1.SetMarkerSize(1.25)
    gr1.GetYaxis().SetTitleSize(0.05)
    gr1.GetYaxis().SetTitleOffset(0.80)
    gr1.GetXaxis().SetLimits(0, 350) # along X
    gr1.GetYaxis().SetRangeUser(yrange[0], yrange[1]) # along Y

    gr2.SetLineWidth(2)
    gr2.SetMarkerStyle(20)
    gr2.SetMarkerSize(1.25)

    legend = ROOT.TLegend(leg_pos[0], leg_pos[1], leg_pos[2], leg_pos[3])
    legend.SetTextSize(0.04)
    legend.AddEntry(gr1, m.labels[0], leg_option)
    legend.AddEntry(gr2, m.labels[1], leg_option)

    # quantify difference
    if draw_lower_pad:
        gr1.GetXaxis().SetLabelSize(0.00)
        gr1.GetXaxis().SetTitleSize(0.00)

        ly, ley = [], []
        for i, resolution_E_set1 in enumerate(ly1):
            err1 = ley1[i]
            err2 = ley2[i]
            resolution_E_set2 = ly2[i]

            #difference = (resolution_E_set1 - resolution_E_set2) / resolution_E_set2  
            ratio = resolution_E_set2 / resolution_E_set1
            uncertainty = ratio * math.sqrt(pow(err1/resolution_E_set1, 2) + pow(err2/resolution_E_set2, 2))
            
            print ">>> check resolution:", i, ratio, uncertainty

            ly.append(ratio)
            ley.append(uncertainty)

        #ytitle = "#frac{#color[%d]{Blue} #minus #color[%d]{Green}}{#color[%d]{Green}}" % (ROOT.kBlue, ROOT.kGreen+3, ROOT.kGreen+3)
        #ytitle = "#frac{#color[%d]{Green}}{#color[%d]{Blue}}" % (ROOT.kGreen+3, ROOT.kBlue)
        ytitle = "#frac{#color[%d]{E (alternative)}}{#color[%d]{E (default)}}" % (ROOT.kGreen+3, ROOT.kBlue)
        gr_ratio = pu.get_graph_from_list("Positron energy (GeV)", ytitle , lx, ly, lex, ley,  ROOT.kBlack)
        gr_ratio.SetLineWidth(2)
        gr_ratio.SetMarkerStyle(20)
        gr_ratio.SetMarkerSize(1.25)

        gr_ratio.GetXaxis().SetLimits(0, 350) # along X
        gr_ratio.GetXaxis().SetTitleFont(42)
        gr_ratio.GetXaxis().SetTitleSize(0.12)
        gr_ratio.GetXaxis().SetLabelFont(42)
        gr_ratio.GetXaxis().SetLabelSize(0.10)
        gr_ratio.GetXaxis().SetLabelOffset(0.06)

        gr_ratio.GetYaxis().SetNdivisions(505)
        if m.energy_type == "set1_set2_MIPs":
            gr_ratio.GetYaxis().SetRangeUser(0.90, 1.70) # along Y
        if m.energy_type == "set1_set2_MeV":
            gr_ratio.GetYaxis().SetRangeUser(0.70, 1.30) # along Y
        gr_ratio.GetYaxis().SetTitleFont(42)
        gr_ratio.GetYaxis().SetTitleSize(0.10)
        gr_ratio.GetYaxis().SetTitleOffset(0.35)
        gr_ratio.GetYaxis().SetLabelFont(42)
        gr_ratio.GetYaxis().SetLabelSize(0.08)
        #gr_ratio.GetYaxis().SetLabelOffset(0.06)

        mainPad, ratPad = pu.init_pads()

        # main pad
        c3.cd()
        mainPad.Draw()
        mainPad.cd()
        ROOT.gPad.SetTicks()

        gr1.Draw("ap")
        gr2.Draw("p;same")
        legend.Draw("same")
        pu.annotate(0.12)

        # ratio pad
        c3.cd()
        ratPad.Draw()
        ratPad.cd()
        gr_ratio.Draw("ap")

        #if title == "resolution_clustered" or title == "resolution_unclustered":
        #    print ">>>>> perform fit"
        #    f1 = ROOT.TF1('f1', "[0]", 0, 320)
        #    gr_ratio.Fit(f1, "", "", 0, 320)

        if draw_goodness:
            #print ">>>>> check line positions", c3.GetUxmin(), reference_line, c3.GetUxmax(), reference_line
            line = ROOT.TLine( c3.GetUxmin(), reference_line, c3.GetUxmax(), reference_line )
            line.SetLineColor(ROOT.kRed)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw()


    else:
        gr1.GetXaxis().SetTitleSize(0.04)

        gr1.Draw("ap")
        gr2.Draw("p;same")
        legend.Draw("same")
        pu.annotate(0.12)

        if draw_goodness:
            #print ">>>>> check line positions", c3.GetUxmin(), reference_line, c3.GetUxmax(), reference_line
            line = ROOT.TLine( c3.GetUxmin(), reference_line, c3.GetUxmax(), reference_line )
            line.SetLineColor(ROOT.kRed)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw()

    token = title if not "resolution" in title else "resolution"
    output = dir_output + "/summary_" + token + "_" + m.energy_type
    c3.SaveAs(output + ".png")
    c3.SaveAs(output + ".pdf")

