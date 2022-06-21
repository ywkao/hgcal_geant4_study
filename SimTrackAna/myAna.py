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

import plot_utils as pu
import MetaData as m

flag_add_reference = True
flag_add_reference = False

eos = "./eos"
#rootfile = eos + "/" + "geantoutput_v3p1.root"
#rootfile = "geantoutput_v3p1.root"
#rootfile = "rootfiles/geantoutput_D86_R80To100_E100.root"
#fin = ROOT.TFile.Open(rootfile, "R")

c1 = ROOT.TCanvas("c1", "", 800, 600)
c1.SetGrid()
c1.SetTicks(1,1)
c1.SetLeftMargin(0.12)
c1.SetRightMargin(0.08)

c2 = ROOT.TCanvas("c2", "", 5600, 2400)
c2.SetLeftMargin(0.)
c2.SetRightMargin(0.)
c2.Divide(7,4)

#--------------------------------------------------
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
    latex.DrawLatex( 0.12, 0.912, "#bf{CMS} #it{work in progress}" )
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
            v_hists = pu.load_single_histogram(fin, varName)
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
            if bool_this_is_eta_phi: v_hists[0].Draw()
            elif bool_ntuple:
                #nt_hit_position = fs->make<TNtuple>("nt_hit_position","nt_hit_position", "r:z:is_Silicon_w120:is_Silicon_w200:is_Silicon_w300:is_Scintillator");
                draw_2D_ntuple("hnew"+"_"+tags[i]+"_w300", v_hists, "is_Silicon_w300>0", ROOT.kMagenta, True)
                draw_2D_ntuple("hnew"+"_"+tags[i]+"_w200", v_hists, "is_Silicon_w200>0", ROOT.kGreen)
                draw_2D_ntuple("hnew"+"_"+tags[i]+"_w120", v_hists, "is_Silicon_w120>0", ROOT.kRed)
            else: continue # no matched situation

            output = dir_output + "/" + varName + "_" + tags[i]
            c1.SaveAs(output + ".png")
            c1.SaveAs(output + ".pdf")

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
        v_gr_ref_mc, v_gr_ref_data, ref_ey = [], [], [0.]*28
        for ene in [300, 100, 20]:
            ref_tag = "multiplicity" if 'multiplicity' in varName else "mips"
            gr = pu.get_graph_from_list(varName, m.test_beam_result["layer_depth"], m.test_beam_result[ref_tag]["mc"][ene], ref_ey, is_number_of_hits)
            gr.SetMarkerColor(ROOT.kRed)
            gr.SetLineColor(ROOT.kRed)
            gr.SetLineStyle(2)
            v_gr_ref_mc.append(gr)

            gr = pu.get_graph_from_list(varName, m.test_beam_result["layer_depth"], m.test_beam_result[ref_tag]["data"][ene], ref_ey, is_number_of_hits)
            gr.SetMarkerColor(ROOT.kBlack)
            gr.SetLineColor(ROOT.kBlack)
            gr.SetLineStyle(2)
            v_gr_ref_data.append(gr)

        #------------------------------
        # logitudinal profile
        #------------------------------
        c1.cd()
        #legend = ROOT.TLegend(0.67, 0.65, 0.87, 0.85)
        #legend = ROOT.TLegend(0.65, 0.65, 0.85, 0.85)
        legend = ROOT.TLegend(0.60, 0.65, 0.85, 0.85)
        legend.SetLineColor(0)
        legend.SetTextSize(0.025)
        for i, gr in enumerate(v_gr):

            if not flag_add_reference:
                print i

                if i==0: gr.Draw("alp")
                else:    gr.Draw('lp;same')
                legend.AddEntry(gr, label[tags[i]], "lp")

                if i+1 == len(v_gr):
                    annotate()
                    legend.Draw("same")

                    output = specified_directory + "/" + varName
                    c1.SaveAs(output + ".png")
                    c1.SaveAs(output + ".pdf")

            else:
                gr.Draw("alp")
                v_gr_ref_mc[i].Draw('lp;same')
                v_gr_ref_data[i].Draw('lp;same')

                legend.Clear()
                legend.AddEntry(v_gr_ref_data[i], "2018 TB Data", "lp")
                legend.AddEntry(v_gr_ref_mc[i], "2018 TB MC", "lp")
                legend.AddEntry(gr, "D86 geometry", "lp")

                annotate()
                legend.Draw("same")

                output = specified_directory + "/" + varName + "_" + tags[i]
                c1.SaveAs(output + ".png")
                c1.SaveAs(output + ".pdf")

#--------------------------------------------------

def run(myfin, mydin):
    global flag_add_reference, myRootfiles, specified_directory
    myRootfiles = myfin
    specified_directory = mydin

    create_directory( specified_directory )
    #make_plot( "hEta", False )
    #make_plot( "hPhi", False )
    make_plot( "nt_hit_position", False )

    thickness = ["120mum", "200mum", "300mum", "total"]
    thickness = ["total", "coarse", "fine"] # consider 120, 200, 300 altogether
    thickness = ["total"]

    for t in thickness:
        make_plot( "multiplicity_simhits_%s" % t , True  )
        make_plot( "total_MIP_%s"            % t , True  )
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

#--------------------------------------------------

if __name__ == "__main__":
    myRootfiles, specified_directory, label = [], "", {}
    colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kGreen, ROOT.kRed]
    colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen]
    colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kBlue-7, ROOT.kRed-7]

    #tags = ["nominal", "Turn_off_Compton", "Replace_PCB_by_Air", "TurnOffCompton_CutEle5mm", "TurnOffCompton_CutEle10mm"]
    #tags = ["nominal", "Turn_off_Compton", "TOC+airPCB", "TOC+airPCB+CutEle5mm", "TOC+airPCB+CutEle10mm"]
    tags = ["nominal", "Turn_off_Compton", "TOC+airPCB", "airPCB"]
    for tag in tags: label[tag] = tag
    run( m.input_files["pcbStudy"]    , eos + "/" + "R80To100_PCB_study"    )
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

