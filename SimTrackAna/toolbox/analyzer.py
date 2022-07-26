#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(2210)
ROOT.gStyle.SetOptFit(1111)

class HitAnalyzer:
    def __init__(self, input_files, directory, tags):
        self.input_files = input_files
        self.output_directory = directory
        self.tags = tags

        self.canvas = ROOT.TCanvas("canvas", "", 800, 600)
        self.canvas.SetGrid()
        self.canvas.SetTicks(1,1)
        self.canvas.SetLeftMargin(0.12)
        self.canvas.SetRightMargin(0.08)

    def check_name(self):
        print "input_files =", self.input_files
        print "output_directory =", self.output_directory

    #--------------------------------------------------

    def loop(self):
        vt, vf = [], []
        for i, rootfile in enumerate(self.input_files):
            print ">>> rootfile", rootfile
            fin = ROOT.TFile.Open(rootfile, "R")
            tree = fin.Get("prodEE_DigiSim/tr_hits")
            vf.append(fin)
            vt.append(tree)

            self.output_pdf_file = self.output_directory + "/output_hits_" + self.tags[i] + ".pdf"
            self.__make_plot(tree)

            #for i, event in enumerate(tree):
            #    print i

    #--------------------------------------------------

    def __make_plot(self, tree):
        # x, y, z, r, layerNo
        
        #------------------------------
        # Rxy-Z distribution
        #------------------------------
        legend = ROOT.TLegend(0.15, 0.65, 0.45, 0.85)
        legend.SetLineColor(0)
        legend.SetTextSize(0.04)

        vh = []
        colors = [ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta]
        wafer_types = ["is_Silicon_w120", "is_Silicon_w200", "is_Silicon_w300"]
        for i, wafer in enumerate(wafer_types):
            name, title = "h_hits_Rxy_Z_%s" % wafer, "Hits R_{xy} vs Z"
            h = ROOT.TH2D(name, title, 160, 300, 380, 600, 0, 300 )
            tree.Draw("r:z>>%s" % name, "%s==1" % wafer)
            htemp = ROOT.gPad.GetPrimitive(name)
            htemp.SetStats(0)
            htemp.SetLineWidth(2)
            htemp.SetLineColor(colors[i])
            htemp.SetMarkerColor(colors[i])
            htemp.GetXaxis().SetTitle("Z (cm)")
            htemp.GetYaxis().SetTitle("R_{xy} (cm)")
            htemp.GetXaxis().SetTitleOffset(1.1)
            htemp.GetYaxis().SetTitleOffset(1.2)
            vh.append(htemp)

            legend.AddEntry(htemp, wafer, "l")

        for i, h in enumerate(vh):
            if i==0: h.Draw()
            else:    h.Draw("same")

        legend.Draw("same")
        output_file = self.output_directory + "/hits_Rxy_Z"
        #self.canvas.SaveAs(output_file + ".png")
        self.canvas.Print(self.output_pdf_file + "(", "pdf")

        #------------------------------
        # X-Y distributions
        #------------------------------
        for layer in range(1,27):

            name = "h_2d_hits_layer%d" % layer
            title = "Layer-0%d" % layer if layer < 10 else "Layer-%d" % layer

            h = ROOT.TH2D(name, title, 680, -170, 170, 680, -170, 170 )
            tree.Draw("y:x>>%s" % name, "layerNo==%d" % layer)
            htemp = ROOT.gPad.GetPrimitive(name)
            htemp.SetStats(0)
            htemp.GetXaxis().SetTitle("X (cm)")
            htemp.GetYaxis().SetTitle("Y (cm)")
            htemp.GetXaxis().SetTitleOffset(1.1)
            htemp.GetYaxis().SetTitleOffset(1.0)

            tag = "0%d" % layer if layer < 10 else "%d" % layer
            output_file = self.output_directory + "/hits_x_y_layer" + tag
            #self.canvas.SaveAs(output_file + ".png")

            if layer == 26:
                self.canvas.Print(self.output_pdf_file + ")", "pdf")
            else:
                self.canvas.Print(self.output_pdf_file, "pdf")
