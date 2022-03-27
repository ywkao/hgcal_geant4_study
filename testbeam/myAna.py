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

import my_plot_utils as pu

eos = "./eos_output"
rootfile = eos + "/" + "geantoutput_v3p1.root"
rootfile = "geantoutput_v3p1.root"
rootfile = "rootfiles/geantoutput_D86_R80To100_E100.root"
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
    latex.DrawLatex( 0.58+rshift, 0.912, "2018 Oct. test beam data v17" )
    #latex.DrawLatex( 0.58+rshift, 0.912, "D86 Simulation with 1,000 events" )
    #latex.DrawLatex( 0.69+rshift, 0.912, "%s fb^{-1} (13 TeV)" % str(lumi["RunII"]) )

def make_plot(varName, bool_make_logitudinal_profile):
    global myRootfiles, specified_directory
    
    # Initiate
    bool_this_is_eta_phi = "Eta" in varName or "Phi" in varName
    #dir_output = specified_directory + "/" + pu.sub_directory[varName]
    dir_output = specified_directory
    create_directory(dir_output)
    processes = [str(i) for i in range(1,27)]
    processes = [str(i) for i in range(1,29)]

    # Load histograms
    vf = []
    v_v_hists = []
    for rootfile in myRootfiles:
        fin = ROOT.TFile.Open(rootfile, "R")
        vf.append(fin)

        v_hists = []
        if bool_this_is_eta_phi:
            v_hists = pu.load_single_histogram(fin, varName)
            v_v_hists.append(v_hists)
        else:
            v_hists = pu.load_histograms(fin, varName, processes)
            v_v_hists.append(v_hists)

    # Eta or Phi
    if bool_this_is_eta_phi:
        c1.cd()
        for i, v_hists in enumerate(v_v_hists):
            v_hists[0].Draw()
            output = dir_output + "/" + varName + "_" + tags[i]
            c1.SaveAs(output + ".png")
            c1.SaveAs(output + ".pdf")

    # Longitdinal profile
    if bool_make_logitudinal_profile:
        v_gr = []

        # loop over energy
        for i, v_hists in enumerate(v_v_hists):
            gr = pu.get_graph(varName, v_hists)
            gr.SetLineColor(colors[i])
            gr.SetMarkerColor(colors[i])
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

        # logitudinal profile
        c1.cd()
        legend = ROOT.TLegend(0.67, 0.65, 0.87, 0.85)
        legend.SetLineColor(0)
        legend.SetTextSize(0.04)
        for i, gr in enumerate(v_gr):
            if i==0: gr.Draw("apc")
            else:    gr.Draw('pc;same')
            legend.AddEntry(gr, label[tags[i]], "l")

        annotate(0.02)
        legend.Draw("same")

        output = specified_directory + "/" + varName
        c1.SaveAs(output + ".png")
        c1.SaveAs(output + ".pdf")

#--------------------------------------------------

def run(myfin, mydin):
    global myRootfiles, specified_directory
    myRootfiles = myfin
    specified_directory = mydin

    create_directory( specified_directory )

    thickness = ["total"] # consider 120, 200, 300 altogether

    for t in thickness:
        make_plot( "rechit_amplitude_%s" % t , True  )
        #make_plot( "rechit_amplitudeHigh_%s" % t , True  )
        #make_plot( "rechit_amplitudeLow_%s"  % t , True  )

#--------------------------------------------------

if __name__ == "__main__":
    tags = ["E300", "E100", "E20"]
    label = {}
    label["E300"] = "300 GeV"
    label["E100"] = "100 GeV"
    label["E20"]  = "20 GeV"

    colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed]
    input_files = [
        "output/my_test_beam_data_e300.root",
        "output/my_test_beam_data_e100.root",
        "output/my_test_beam_data_e20.root",
    ]

    myRootfiles, specified_directory = [], ""
    run( input_files[0:3], eos + "/" + "data"  )

