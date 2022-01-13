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

eos = "./eos_output"
rootfile = eos + "/" + "geantoutput_v3p1.root"
rootfile = "geantoutput_v3p1.root"
fin = ROOT.TFile.Open(rootfile, "R")

c1 = ROOT.TCanvas("c1", "", 800, 600)
c1.SetGrid()
c1.SetTicks(1,1)
c1.SetLeftMargin(0.12)
c1.SetRightMargin(0.08)

#--------------------------------------------------

def make_plot(varName, bool_make_profile):

    dir_output = eos + "/" + pu.sub_directory[varName]
    if not os.path.isdir(dir_output):
        subprocess.call("mkdir %s" % dir_output, shell=True)
        subprocess.call("cp -p %s/index.php %s" % (eos, dir_output), shell=True)

    processes = ["1", "13", "26"]
    processes = [str(i) for i in range(1,27)]
    v_hists = pu.load_histograms(fin, varName, processes)

    # longitdinal profile
    if bool_make_profile:
        gr = pu.get_graph(varName, v_hists)
        gr.Draw("apc")
        output = eos + "/" + varName
        c1.SaveAs(output + ".png")

    # individual plots
    for i, h in enumerate(v_hists):
        h.Draw()
        output = dir_output + "/" + varName + "_" + processes[i]
        c1.SaveAs(output + ".png")
        c1.SaveAs(output + ".pdf")

#--------------------------------------------------

def run():
    c1.cd()
    thickness = [200]
    thickness = [120, 200, 300]
    for t in thickness:
        make_plot("total_ADC_%dmum"            % t , True  )
        make_plot("total_MIP_%dmum"            % t , True  )
        make_plot("total_SIM_%dmum"            % t , True  )
        make_plot("multiplicity_digis_%dmum"   % t , True  )
        make_plot("multiplicity_simhits_%dmum" % t , True  )

        make_plot("ADC_%dmum_layer"            % t , True  )
        make_plot("MIP_%dmum_layer"            % t , True  )
        make_plot("SIM_%dmum_layer"            % t , True  )

        make_plot("ADC_SimhitE_%dmum_layer"    % t , False )
        make_plot("ADC_MIP_%dmum_layer"        % t , False )
        make_plot("MIP_SimhitE_%dmum_layer"    % t , False )

if __name__ == "__main__":
    run()
    #subprocess.call("ls -lhrt %s" % dir_output, shell=True)
