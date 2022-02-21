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
fin = ROOT.TFile.Open(rootfile, "R")

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

def make_plot(myRootfiles, varName, bool_make_logitudinal_profile):

    dir_output = eos + "/" + pu.sub_directory[varName]
    if not os.path.isdir(dir_output):
        subprocess.call("mkdir %s" % dir_output, shell=True)
        subprocess.call("cp -p %s/index.php %s" % (eos, dir_output), shell=True)

    processes = ["1", "13", "26"]
    processes = [str(i) for i in range(1,27)]

    # loading
    vf = []
    v_v_hists = []
    for rootfile in myRootfiles:
        fin = ROOT.TFile.Open(rootfile, "R")
        vf.append(fin)

        v_hists = []
        if "Eta" in varName or "Phi" in varName:
            v_hists = pu.load_single_histogram(fin, varName)
            v_v_hists.append(v_hists)
        else:
            v_hists = pu.load_histograms(fin, varName, processes)
            v_v_hists.append(v_hists)

    # longitdinal profile
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

            output = eos + "/" + varName + "_all_plots_" + tags[i]
            c2.SaveAs(output + ".png")
            c2.SaveAs(output + ".pdf")

            # individual plots
            c1.cd()
            for j, h in enumerate(v_hists):
                h.Draw()
                output = dir_output + "/" + varName + "_" + processes[j] + "_" + tags[i]
                c1.SaveAs(output + ".png")
                c1.SaveAs(output + ".pdf")

        # logitudinal profile
        c1.cd()
        for i, gr in enumerate(v_gr):
            if i==0: gr.Draw("apc")
            else:    gr.Draw('pc;same')

        output = eos + "/" + varName
        c1.SaveAs(output + ".png")
        c1.SaveAs(output + ".pdf")

#--------------------------------------------------

def run(input_files):
    make_plot(input_files, "hEta", False )
    make_plot(input_files, "hPhi", False )

    thickness = ["120mum", "200mum", "300mum", "total"]
    thickness = ["total"] # consider 120, 200, 300 altogether

    for t in thickness:

        make_plot(input_files, "total_MIP_%s"            % t , True  )

        continue

        make_plot(input_files, "total_ADC_%s"            % t , True  )
        make_plot(input_files, "total_MIP_%s"            % t , True  )
        make_plot(input_files, "total_SIM_%s"            % t , True  )
        make_plot(input_files, "multiplicity_digis_%s"   % t , True  )
        make_plot(input_files, "multiplicity_simhits_%s" % t , True  )

        make_plot(input_files, "ADC_SimhitE_%s_layer"    % t , False )
        make_plot(input_files, "ADC_MIP_%s_layer"        % t , False )
        make_plot(input_files, "MIP_SimhitE_%s_layer"    % t , False )

        continue

        make_plot(input_files, "ADC_%s"            % t , True  )
        make_plot(input_files, "MIP_%s"            % t , True  )
        make_plot(input_files, "SIM_%s"            % t , True  )

if __name__ == "__main__":
    colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed]
    myRootfiles = [
        "rootfiles/geantoutput_D86_R35To60_E300.root",
        "rootfiles/geantoutput_D86_R35To60_E100.root",
        "rootfiles/geantoutput_D86_R35To60_E20.root",
        "rootfiles/geantoutput_D86_R80To100_E300.root",
        "rootfiles/geantoutput_D86_R80To100_E100.root",
        "rootfiles/geantoutput_D86_R80To100_E20.root",
    ]

    tags = ["E300", "E100", "E20"]

    run( myRootfiles[0:3] )
    #subprocess.call("ls -lhrt %s" % dir_output, shell=True)

