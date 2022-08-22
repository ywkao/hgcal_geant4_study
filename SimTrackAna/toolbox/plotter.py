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
specified_directory = ""
flag_add_reference = False
fit_constraints = {}
colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kBlue-7, ROOT.kMagenta, ROOT.kRed-7]

#--------------------------------------------------
# plotting functions
#--------------------------------------------------
def make_plot(varName, dir_output="", bool_make_logitudinal_profile=False):
    global myRootfiles, specified_directory, flag_add_reference
    is_number_of_hits = "multiplicity" in varName
    
    #++++++++++++++++++++++++++++++
    # Initiate
    #++++++++++++++++++++++++++++++
    bool_ntuple = "nt_" in varName
    bool_this_is_eta_phi = "Eta" in varName or "Phi" in varName
    bool_single_figures = bool_ntuple or bool_this_is_eta_phi
    processes = [str(i) for i in range(1,27)]
    if len(dir_output) == 0:
        dir_output = specified_directory + "/" + pu.sub_directory[varName]

    #++++++++++++++++++++++++++++++
    # Load histograms
    #++++++++++++++++++++++++++++++
    vf = []
    v_v_hists = []
    for rootfile in myRootfiles:
        print ">>> rootfile", rootfile
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
    global myRootfiles, specified_directory, fit_constraints
    
    xtitle, outputName, hname_odd, hname_even, labels = pu.get_strings_for_simple_plot(energyType, selection)

    #++++++++++++++++++++++++++++++
    # Load histograms
    #++++++++++++++++++++++++++++++
    vf, v_v_hists = [], []
    for rootfile in myRootfiles:
        print ">>> rootfile", rootfile
        fin = ROOT.TFile.Open(rootfile, "R")
        vf.append(fin)

        v_hists = []
        v_hists.append( pu.load_single_histogram(fin, hname_odd) )
        v_hists.append( pu.load_single_histogram(fin, hname_even) )
        v_v_hists.append(v_hists)

    c1.cd()
    for i, v_hists in enumerate(v_v_hists):
        print "\n----------------------------------------------------------------------------------------------------\n"
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

            pu.draw_and_fit_a_histogram( c1, v_hists[0], [energyType, tags[i], labels[0]], xtitle, max_value, xRanges[i], ROOT.kBlue   , [0.60, 0.66, 0.88, 0.86] )
            pu.draw_and_fit_a_histogram( c1, v_hists[1], [energyType, tags[i], labels[1]], xtitle, max_value, xRanges[i], ROOT.kGreen+3, [0.60, 0.42, 0.88, 0.62] )

            v_hists[0].Draw("same")
            v_hists[0].GetFunction("gaus").Draw("same")

            #--------------------------------------------------
            # result
            #--------------------------------------------------
            latex.SetTextColor(ROOT.kBlue)
            latex.DrawLatex( 0.45, 0.30, "#sigma#left(E_{%s}#right) / #bar{E}_{%s} = %.4f #pm %.4f" % (labels[0], labels[0], pu.sigmaEoverE[0], pu.error_sigmaEoverE[0]) )
            latex.SetTextColor(ROOT.kGreen+3)
            latex.DrawLatex( 0.45, 0.20, "#sigma#left(E_{%s}#right) / #bar{E}_{%s} = %.4f #pm %.4f" % (labels[1], labels[1], pu.sigmaEoverE[1], pu.error_sigmaEoverE[1]) )

        else:
            xRanges = fit_constraints["set0"][energyType]["xRanges"]

            pu.draw_and_fit_a_histogram( c1, v_hists[0], [energyType, tags[i], labels[0]] , xtitle, max_value, xRanges[i], ROOT.kBlue  , [0.60, 0.66, 0.88, 0.86] )

            latex.SetTextColor(ROOT.kBlue)
            latex.DrawLatex( 0.45, 0.30, "#sigma#left(E_{%s}#right) / #bar{E}_{%s} = %.4f #pm %.4f" % (labels[0], labels[0], pu.sigmaEoverE[0], pu.error_sigmaEoverE[0]) )

        c1.Update()
        pu.annotate()
        output = dir_output + "/" + outputName + tags[i]
        c1.SaveAs(output + ".png")
        c1.SaveAs(output + ".pdf")

#----------------------------------------------------------------------------------------------------

def run_linear_fit(dir_output, label, dx, dy):
    #global specified_directory

    Energy = ["E20", "E60", "E100", "E175", "E225", "E300"]
    lx  = [ dx[ene]["mean"]  for ene in Energy ]
    ly  = [ dy[ene]["mean"]  for ene in Energy ]
    lex = [ dx[ene]["sigma"] for ene in Energy ]
    ley = [ dy[ene]["sigma"] for ene in Energy ]

    gr = pu.get_graph_from_list("Energy (MIPs)", "Generated shower energy (MeV)", lx, ly, lex, ley, ROOT.kBlack)

    c1.cd()
    c1.Clear()
    gr.Draw("ap")
    gr.GetXaxis().SetLimits(0, 30000) # alogn X
    gr.GetYaxis().SetRangeUser(0, 300) # v1p1, v1p2
    gr.GetYaxis().SetRangeUser(0, 100) # v2p1, v2p2

    f1 = ROOT.TF1('f1', "[0] + [1]*x", 0, 30000)
    gr.Fit(f1, "", "", 0, 30000)

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

def run_resolution_summary(dir_output, labels, dy1, dy2):
    #global specified_directory

    Energy = ["E20", "E60", "E100", "E175", "E225", "E300"]
    lx  = [ float(ene.split('E')[1]) for ene in Energy ]
    lex = [0.]*6

    ly1, ly2, ley1, ley2 = [], [], [], []
    for ene in Energy:
        fit_mean = dy1[ene]["mean"]
        fit_sigma = dy1[ene]["sigma"]
        fitError_mean  = dy1[ene]["error_mean"]
        fitError_sigma = dy1[ene]["error_sigma"]

        ratio = fit_sigma/fit_mean
        uncertainty = ratio * math.sqrt( math.pow(fitError_mean/fit_mean, 2) + math.pow(fitError_sigma/fit_sigma, 2) )

        ly1.append(ratio)
        ley1.append(uncertainty)


        fit_mean = dy2[ene]["mean"]
        fit_sigma = dy2[ene]["sigma"]
        fitError_mean  = dy2[ene]["error_mean"]
        fitError_sigma = dy2[ene]["error_sigma"]

        ratio = fit_sigma/fit_mean
        uncertainty = ratio * math.sqrt( math.pow(fitError_mean/fit_mean, 2) + math.pow(fitError_sigma/fit_sigma, 2) )

        ly2.append(ratio)
        ley2.append(uncertainty)


    gr1 = pu.get_graph_from_list("Positron energy (GeV)", "#sigma#left(E#right) / #bar{E}", lx, ly1, lex, ley1,  ROOT.kBlue)
    gr2 = pu.get_graph_from_list("Positron energy (GeV)", "#sigma#left(E#right) / #bar{E}", lx, ly2, lex, ley2,  ROOT.kGreen+3)


    #legend = ROOT.TLegend(0.52, 0.65, 0.88, 0.85)
    #legend.SetLineColor(0)
    legend = ROOT.TLegend(0.62, 0.65, 0.89, 0.85)
    legend.SetTextSize(0.04)
    legend.AddEntry(gr1, labels[0], "ep")
    legend.AddEntry(gr2, labels[1], "ep")

    c3.cd()
    c3.Clear()

    gr1.SetLineWidth(2)
    gr1.SetMarkerStyle(20)
    gr1.SetMarkerSize(1.25)
    gr1.GetXaxis().SetTitleSize(0.04)
    gr1.GetYaxis().SetTitleSize(0.04)
    gr1.Draw("ap")
    gr1.GetXaxis().SetLimits(0, 350) # alogn X
    gr1.GetYaxis().SetRangeUser(0.018, 0.043) # v1p1, v1p2
    gr1.GetYaxis().SetRangeUser(0.0, 0.040) # v2p1, v2p2

    gr2.SetLineWidth(2)
    gr2.SetMarkerStyle(20)
    gr2.SetMarkerSize(1.25)
    gr2.Draw("p;same")
    legend.Draw("same")

    pu.annotate(0.12)
    output = dir_output + "/summary_resolution"
    c3.SaveAs(output + ".png")
    c3.SaveAs(output + ".pdf")

