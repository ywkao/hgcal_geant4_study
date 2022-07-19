#!/usr/bin/env python
import os, subprocess
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

import toolbox.plot_utils as pu
import toolbox.MetaData as m

flag_add_reference = True
flag_add_reference = False

c1 = ROOT.TCanvas("c1", "", 800, 600)
c1.SetGrid()
c1.SetTicks(1,1)
c1.SetLeftMargin(0.12)
c1.SetRightMargin(0.08)

c2 = ROOT.TCanvas("c2", "", 5600, 2400)
c2.SetLeftMargin(0.)
c2.SetRightMargin(0.)
c2.Divide(7,4)

X0 = {
    "uniform" : [1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.],
    "alternative" : [0.564,1.567,2.547,3.549,4.528,5.531,6.509,7.512,8.49,9.493,10.472,11.474,12.453,13.455,14.434,15.437,16.415,17.418,18.975,19.978,21.536,22.538,24.096,25.099,26.656,27.659],
}

eos = "./eos"
def create_directory(dir_output):
    if not os.path.isdir(dir_output):
        subprocess.call("mkdir %s" % dir_output, shell=True)
        subprocess.call("cp -p %s/index.php %s" % (eos, dir_output), shell=True)

#----------------------------------------------------------------------------------------------------

def make_plot(varName, bool_make_logitudinal_profile):
    global myRootfiles, specified_directory, flag_add_reference
    is_number_of_hits = "multiplicity" in varName
    
    #++++++++++++++++++++++++++++++
    # Initiate
    #++++++++++++++++++++++++++++++
    bool_ntuple = "nt_" in varName
    bool_this_is_eta_phi = "Eta" in varName or "Phi" in varName
    bool_single_figures = bool_ntuple or bool_this_is_eta_phi
    dir_output = specified_directory + "/" + pu.sub_directory[varName]
    create_directory(dir_output)
    processes = [str(i) for i in range(1,27)]

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
                v_hists[0].GetXaxis().SetRangeUser(1., 3.) # Eta
                v_hists[0].Draw()
            elif bool_ntuple:
                #nt_hit_position = fs->make<TNtuple>("nt_hit_position","nt_hit_position", "r:z:is_Silicon_w120:is_Silicon_w200:is_Silicon_w300:is_Scintillator");
                pu.draw_2D_ntuple("hnew"+"_"+tags[i]+"_w300", v_hists, "is_Silicon_w300>0", ROOT.kMagenta, True)
                pu.draw_2D_ntuple("hnew"+"_"+tags[i]+"_w200", v_hists, "is_Silicon_w200>0", ROOT.kGreen)
                pu.draw_2D_ntuple("hnew"+"_"+tags[i]+"_w120", v_hists, "is_Silicon_w120>0", ROOT.kRed)
            else: continue # no matched situation

            output = dir_output + varName + "_" + tags[i]
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

            output = specified_directory + "/" + varName + "_all_plots_" + tags[i]
            c2.SaveAs(output + ".png")
            c2.SaveAs(output + ".pdf")

            continue

            # individual plots
            c1.cd()
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

                    output = specified_directory + "/" + varName
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

                output = specified_directory + "/" + varName + "_" + tags[i]
                c1.SaveAs(output + ".png")
                c1.SaveAs(output + ".pdf")

#----------------------------------------------------------------------------------------------------

def make_simple_plot(energyType, selection):
    global myRootfiles, specified_directory, fit_constraints
    
    xRanges = fit_constraints[energyType]["xRanges"]

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
            pu.draw_and_fit_a_histogram( c1, v_hists[0], [energyType, tags[i], labels[0]] , xtitle, max_value, xRanges[i], ROOT.kBlue  , [0.60, 0.66, 0.88, 0.86] )
            pu.draw_and_fit_a_histogram( c1, v_hists[1], [energyType, tags[i], labels[1]], xtitle, max_value, xRanges[i], ROOT.kGreen+3, [0.60, 0.42, 0.88, 0.62] )

            v_hists[0].Draw("same")
            v_hists[0].GetFunction("gaus").Draw("same")

            #--------------------------------------------------
            # result
            #--------------------------------------------------
            latex.SetTextColor(ROOT.kBlue)
            latex.DrawLatex( 0.55, 0.30, "#sigma#left(E_{%s}#right) / #bar{E}_{%s} = %.4f" % (labels[0], labels[0], pu.sigmaEoverE[0]) )
            latex.SetTextColor(ROOT.kGreen+3)
            latex.DrawLatex( 0.55, 0.20, "#sigma#left(E_{%s}#right) / #bar{E}_{%s} = %.4f" % (labels[1], labels[1], pu.sigmaEoverE[1]) )

        else:
            pu.draw_and_fit_a_histogram( c1, v_hists[0], [energyType, tags[i], labels[0]] , xtitle, max_value, xRanges[i], ROOT.kBlue  , [0.60, 0.66, 0.88, 0.86] )

            latex.SetTextColor(ROOT.kBlue)
            latex.DrawLatex( 0.55, 0.30, "#sigma#left(E_{%s}#right) / #bar{E}_{%s} = %.4f" % (labels[0], labels[0], pu.sigmaEoverE[0]) )

        c1.Update()
        pu.annotate()
        output = specified_directory + "/" + outputName + tags[i]
        c1.SaveAs(output + ".png")
        c1.SaveAs(output + ".pdf")

#----------------------------------------------------------------------------------------------------

def run(myfin, mydin):
    global flag_add_reference, myRootfiles, specified_directory
    myRootfiles = myfin
    specified_directory = mydin

    create_directory( specified_directory )
    make_plot( "hEta", False )
    #make_plot( "hPhi", False )

    #make_plot( "nt_hit_position", False )
    #make_simple_plot("MIP", "odd_even")
    #make_simple_plot("SIM", "odd_even")

    make_simple_plot("MIP", "set1_set2")
    make_simple_plot("MIP", "set0")

    make_simple_plot("SIM", "set1_set2")
    make_simple_plot("SIM", "set0")

    make_simple_plot("ENE", "set1_set2")
    make_simple_plot("ENE", "set0")

    return

    thickness = ["120mum", "200mum", "300mum", "total"]
    thickness = ["total", "coarse", "fine"] # consider 120, 200, 300 altogether
    thickness = ["total"]

    for t in thickness:
        make_plot( "multiplicity_simhits_%s" % t , True  )
        make_plot( "total_MIP_%s"            % t , True  )
        make_plot( "total_SIM_%s"            % t , True  )
        continue

        make_plot( "total_ADC_%s"            % t , True  )
        make_plot( "total_MIP_%s"            % t , True  )
        make_plot( "total_SIM_%s"            % t , True  )
        make_plot( "multiplicity_digis_%s"   % t , True  )
        make_plot( "multiplicity_simhits_%s" % t , True  )

        make_plot( "ADC_SimhitE_%s_layer"    % t , False )
        make_plot( "ADC_MIP_%s_layer"        % t , False )
        make_plot( "MIP_SimhitE_%s_layer"    % t , False )

        make_plot( "ADC_%s"                  % t , True  )
        make_plot( "MIP_%s"                  % t , True  )
        make_plot( "SIM_%s"                  % t , True  )

#----------------------------------------------------------------------------------------------------

def run_linear_fit(label, dx, dy):
    #global specified_directory

    Energy = ["E20", "E60", "E100", "E175", "E225", "E300"]
    lx  = [ dx[ene]["mean"]  for ene in Energy ]
    ly  = [ dy[ene]["mean"]  for ene in Energy ]
    lex = [ dx[ene]["sigma"] for ene in Energy ]
    ley = [ dy[ene]["sigma"] for ene in Energy ]

    gr = pu.get_graph_from_list("Energy (MIPs)", "Generated shower energy (keV)", lx, ly, lex, ley, ROOT.kBlack)

    c1.cd()
    c1.Clear()
    gr.Draw("ap")
    gr.GetXaxis().SetLimits(0, 30000) # alogn X
    gr.GetYaxis().SetRangeUser(0, 300000)

    f1 = ROOT.TF1('f1', "[0] + [1]*x", 0, 30000)
    gr.Fit(f1, "", "", 0, 30000)

    my_stat_pos = [0.42, 0.87, 0.15, 0.15]
    ROOT.gStyle.SetStatX(my_stat_pos[0])
    ROOT.gStyle.SetStatY(my_stat_pos[1])
    ROOT.gStyle.SetStatW(my_stat_pos[2])
    ROOT.gStyle.SetStatH(my_stat_pos[3])

    pu.annotate()
    directory = eos + "/R80To150_linearFit_v1p2/"
    output = directory + "correction_generatedShowerEnergy_MIPs_" + label
    create_directory(directory)
    c1.SaveAs(output + ".png")
    c1.SaveAs(output + ".pdf")

#----------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    myRootfiles, specified_directory, label = [], "", {}
    colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kBlue-7, ROOT.kMagenta, ROOT.kRed-7]

    #----------------------------------------------------------------------------------------------------

    tags = ["E300", "E100", "E20"]
    fit_constraints = m.fit_constraints_v1
    for tag in tags: label[tag] = tag.split("E")[1] + " GeV"
    run( m.input_files["R80To150"], eos + "/" + "R80To150_v4p1" )

    tags = ["E225", "E175", "E60"]
    fit_constraints = m.fit_constraints_v2
    for tag in tags: label[tag] = tag.split("E")[1] + " GeV"
    run( m.input_files["R80To150_v2"], eos + "/" + "R80To150_v4p1" )

    #----------------------------------------------------------------------------------------------------

    exit()

    fit_result = pu.fit_result

    for energyType in ["MIP", "SIM"]:
        for tag in ["E20", "E60", "E100", "E175", "E225", "E300"]:
            print energyType, tag
            print fit_result[energyType]["odd" ][tag]
            print fit_result[energyType]["even"][tag]
            print fit_result[energyType]["set0"][tag]
            print fit_result[energyType]["set1"][tag]
            print fit_result[energyType]["set2"][tag]
            print ""

    run_linear_fit("even", fit_result["MIP"]["odd" ], fit_result["SIM"]["odd" ])
    run_linear_fit("odd" , fit_result["MIP"]["even"], fit_result["SIM"]["even"])
    run_linear_fit("set0", fit_result["MIP"]["set0"], fit_result["SIM"]["set0"])
    run_linear_fit("set0_set1", fit_result["MIP"]["set1"], fit_result["SIM"]["set0"])
    run_linear_fit("set0_set2", fit_result["MIP"]["set2"], fit_result["SIM"]["set0"])

    #----------------------------------------------------------------------------------------------------

    exit()

# previous stdy {{{
    tags = ["nominal", "inverse_dEdx_X0", "dEdx_divided_X0", "inverse_X0", "applied_dEdx"]
    for tag in tags: label[tag] = tag
    run( m.input_files["X0_corrections"]    , eos + "/" + "R80To100_study_with_corrections_v2"    )
    exit()

    tags = ["nominal", "with_dEdx_weight"]
    for tag in tags: label[tag] = tag
    run( m.input_files["dEdxStudy"]    , eos + "/" + "R80To100_study_with_dEdx_weights"    )
    exit()

    #tags = ["nominal", "Turn_off_Compton", "Replace_PCB_by_Air", "TurnOffCompton_CutEle5mm", "TurnOffCompton_CutEle10mm"]
    #tags = ["nominal", "Turn_off_Compton", "TOC+airPCB", "TOC+airPCB+CutEle5mm", "TOC+airPCB+CutEle10mm"]
    tags = ["nominal", "Turn_off_Compton", "TOC+airPCB", "airPCB"]
    for tag in tags: label[tag] = tag
    run( m.input_files["pcbStudy"]    , eos + "/" + "R80To100_PCB_study_with_dEdx_weights_v2"    )
    exit()

    tags = ["nominal", "Turn_off_Compton", "TurnOffCompton_CutEle1mm", "TurnOffCompton_CutEle10mm", "TurnOffCompton_CutEle100mm", "TurnOffCompton_CutEle100mm"]
    for tag in tags: label[tag] = tag
    run( m.input_files["turnOffCompton_ProdCutElectron"]    , eos + "/" + "R80To100_turnOffCompton_ProdCutElectron_v3"    )
    exit()

    tags = ["nominal", "Turn_off_Compton", "TurnOffCompton_CutEle1mm", "TurnOffCompton_CutEle5mm", "TurnOffCompton_CutEle10mm", "TurnOffCompton_CutEle50mm", "TurnOffCompton_CutEle100mm"]
    for tag in tags: label[tag] = tag
    run( m.input_files["turnOffCompton_ProdCutElectron"]    , eos + "/" + "R80To100_turnOffCompton_ProdCutElectron_v2"    )
    exit()

    # extra
    tags = ["nominal", "Turn_off_Compton", "Turn_off_Conversion", "Turn_off_both"]
    for tag in tags: label[tag] = tag
    run( m.input_files["extraStudy"]    , eos + "/" + "R80To100_extra_study"    )
    exit()

    # compton
    tags = ["nominal", "Turn_off_Compton", "Turn_off_both", "ProdCut_electron"]
    for tag in tags: label[tag] = tag
    run( m.input_files["turnOffCompton"]    , eos + "/" + "R80To100_TurnOffCompton_comparison"    )
    exit()

    # prodCut study
    tags = ["nominal", "ProdCut_electron", "ProdCut_photon", "ProdCut_egamma"]
    for tag in tags: label[tag] = tag
    run( m.input_files["ProdCut1mm"]    , eos + "/" + "R80To100_ProdCut1mm"    )
    run( m.input_files["ProdCut100mm"]  , eos + "/" + "R80To100_ProdCut100mm"  )
    run( m.input_files["ProdCut1000mm"] , eos + "/" + "R80To100_ProdCut1000mm" )
    exit()

    #++++++++++++++++++++++++++++++
    # need to check root files
    #++++++++++++++++++++++++++++++

    # init study
    tags = ["E300", "E100", "E20"]
    for tag in tags: label[tag] = tag.split("E")[1] + " GeV"
    run( m.input_files["R35To60"], eos + "/" + "R35To60"  )
    run( m.input_files["R80To100"], eos + "/" + "R80To100" )

    # muon study
    tags = ["E100"]
    run( [input_files["muon"]], eos + "/" + "muon" )
#}}}
