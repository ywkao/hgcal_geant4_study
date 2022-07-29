#!/usr/bin/env python
import copy
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

    #====================================================================================================

    def loop(self):
        vt, vf = [], []
        for i, rootfile in enumerate(self.input_files):
            #if not i==1: continue # use E100 for check

            print ">>> rootfile", rootfile
            fin = ROOT.TFile.Open(rootfile, "R")
            tree = fin.Get("prodEE_DigiSim/tr_hits")
            vf.append(fin)
            vt.append(tree)

            # add energy selection
            self.output_pdf_file = self.output_directory + "/output_hits_" + self.tags[i] + ".pdf"
            self.__make_plots(tree)

            break

            # total event
            self.check_individual_event = False 
            self.output_pdf_file = self.output_directory + "/output_hits_" + self.tags[i] + ".pdf"
            self.__retrieve_hit_plots(tree)

            break

            # specific event
            self.check_individual_event = True
            for evtNo in range(3, 5):
                self.evtNo = evtNo 
                self.output_pdf_file = self.output_directory + "/output_hits_" + self.tags[i] + "_evtNo%d.pdf" % self.evtNo
                self.__retrieve_hit_plots(tree)

    #====================================================================================================

    def __make_plots(self, tree):
        # branches: evtNo, layerNo, x, y, z, e, r, eta, phi, is_Silicon_w120, is_Silicon_w200, is_Silicon_w300, is_Scintillator
        
        #----------------------------------------
        # store total energy for each event
        #----------------------------------------
        total_energy = [0.]*1000
        for i, hit in enumerate(tree):
            total_energy[hit.evtNo] += hit.e

        #print "total_energy =", total_energy

        #----------------------------------------
        # register histograms
        #----------------------------------------
        legend = ROOT.TLegend(0.15, 0.65, 0.45, 0.85)
        legend.SetLineColor(0)
        legend.SetTextSize(0.04)

        vh_Rxy_z = []
        colors = [ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta]
        colors = [ROOT.kBlack, ROOT.kBlack, ROOT.kBlack]
        wafer_types = ["is_Silicon_w120", "is_Silicon_w200", "is_Silicon_w300"]
        for i, wafer in enumerate(wafer_types):
            name, title = "h_hits_Rxy_Z_%s" % wafer, "Hits R_{xy} vs Z"
            h = ROOT.TH2D(name, title, 160, 300, 380, 600, 0, 300 )
            h.SetStats(0)
            h.SetLineWidth(2)
            h.SetLineColor(colors[i])
            h.SetMarkerColor(colors[i])
            h.GetXaxis().SetTitle("Z (cm)")
            h.GetYaxis().SetTitle("R_{xy} (cm)")
            h.GetXaxis().SetTitleOffset(1.1)
            h.GetYaxis().SetTitleOffset(1.2)
            vh_Rxy_z.append(h)

            legend.AddEntry(h, wafer, "l")

        vh_x_y = []
        for layer in range(1,27):
            name = "h_2d_hits_layer%d" % layer
            title = "Layer-0%d" % layer if layer < 10 else "Layer-%d" % layer
            h = ROOT.TH2D(name, title, 680, -170, 170, 680, -170, 170 )
            h.SetStats(0)
            h.GetXaxis().SetTitle("X (cm)")
            h.GetYaxis().SetTitle("Y (cm)")
            h.GetXaxis().SetTitleOffset(1.1)
            h.GetYaxis().SetTitleOffset(1.0)
            vh_x_y.append(h)

        e_set0 = ROOT.TH1D("e_set0", ";E_{set0} (keV);Entries", 400, 0, 400000)
        e_set0.SetTitle("")
        e_set0.SetLineWidth(2)
        e_set0.SetLineColor(ROOT.kBlack)
        e_set0.GetXaxis().SetRangeUser(0, 150000)
        e_set0.GetXaxis().SetTitleOffset(1.1)

        h_Rxy_e = ROOT.TH2D("h_Rxy_e", ";E_{set0} (keV);R_{xy} (cm)", 60, 40000, 100000, 600, 0, 300)

        h_eta = ROOT.TH1D("h_eta", ";Eta;Entries", 20, 1., 3.)

        #----------------------------------------
        # deepcopy
        #----------------------------------------
        a, b = copy.deepcopy(vh_Rxy_z), copy.deepcopy(vh_x_y)
        c = e_set0.Clone()

        #----------------------------------------
        # fill histograms
        #----------------------------------------
        for evtNo in range(len(total_energy)):
            if total_energy[evtNo] > 65000.:
                e_set0.Fill(total_energy[evtNo])
            else:
                c.Fill(total_energy[evtNo])

        for i, hit in enumerate(tree):
            h_Rxy_e.Fill(total_energy[hit.evtNo], hit.r)
            h_eta.Fill(hit.eta)

            if total_energy[hit.evtNo] > 65000.:
                if hit.is_Silicon_w120: vh_Rxy_z[0].Fill(hit.z, hit.r)
                if hit.is_Silicon_w200: vh_Rxy_z[1].Fill(hit.z, hit.r)
                if hit.is_Silicon_w300: vh_Rxy_z[2].Fill(hit.z, hit.r)

                idx = hit.layerNo-1
                vh_x_y[idx].Fill(hit.x, hit.y)

            # to be observed (kRed)
            else:
                if hit.is_Silicon_w120: a[0].Fill(hit.z, hit.r)
                if hit.is_Silicon_w200: a[1].Fill(hit.z, hit.r)
                if hit.is_Silicon_w300: a[2].Fill(hit.z, hit.r)

                idx = hit.layerNo-1
                b[idx].Fill(hit.x, hit.y)

        #----------------------------------------
        # make plots
        #----------------------------------------
        for i, h in enumerate(vh_Rxy_z):
            if i==0: h.Draw()
            else:    h.Draw("same")

        for i, h in enumerate(a):
            h.SetMarkerColor(ROOT.kRed)
            h.Draw("same")

        #legend.Draw("same")
        output_file = self.output_directory + "/hits_Rxy_Z"
        self.canvas.SaveAs(output_file + ".png")
        self.canvas.Print(self.output_pdf_file + "(", "pdf")

        for i, h in enumerate(vh_x_y):
            h.Draw()
            b[i].SetMarkerColor(ROOT.kRed)
            b[i].Draw("same")

            layer = i+1
            tag = "0%d" % layer if layer < 10 else "%d" % layer
            output_file = self.output_directory + "/hits_x_y_layer" + tag
            self.canvas.SaveAs(output_file + ".png")

            if layer == 26:
                self.canvas.Print(self.output_pdf_file + ")", "pdf")
            else:
                self.canvas.Print(self.output_pdf_file, "pdf")

        c.SetFillColorAlpha(ROOT.kRed, 0.5)
        e_set0.Draw()
        c.Draw("same")
        self.canvas.SaveAs(self.output_directory + "/h_sim_E_set0.png")

        ROOT.TGaxis.SetExponentOffset(0.01, -0.05, "x")
        h_Rxy_e.Draw("colz")
        self.canvas.SaveAs(self.output_directory + "/h_Rxy_e.png")

        h_eta.Draw()
        self.canvas.SaveAs(self.output_directory + "/h_eta.png")


    def __retrieve_hit_plots(self, tree):
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
            selection = "evtNo==%d && %s==1" % (self.evtNo, wafer) if self.check_individual_event else "%s==1" % wafer
            h = ROOT.TH2D(name, title, 160, 300, 380, 600, 0, 300 )
            tree.Draw("r:z>>%s" % name, selection)
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
        self.canvas.SaveAs(output_file + ".png")
        self.canvas.Print(self.output_pdf_file + "(", "pdf")

        #------------------------------
        # X-Y distributions
        #------------------------------
        for layer in range(1,27):

            name = "h_2d_hits_layer%d" % layer
            title = "Layer-0%d" % layer if layer < 10 else "Layer-%d" % layer
            selection = "evtNo==%d && layerNo==%d" % (self.evtNo, layer) if self.check_individual_event else "layerNo==%d" % layer

            h = ROOT.TH2D(name, title, 680, -170, 170, 680, -170, 170 )
            tree.Draw("y:x>>%s" % name, selection)
            htemp = ROOT.gPad.GetPrimitive(name)
            htemp.SetStats(0)
            htemp.GetXaxis().SetTitle("X (cm)")
            htemp.GetYaxis().SetTitle("Y (cm)")
            htemp.GetXaxis().SetTitleOffset(1.1)
            htemp.GetYaxis().SetTitleOffset(1.0)

            tag = "0%d" % layer if layer < 10 else "%d" % layer
            output_file = self.output_directory + "/hits_x_y_layer" + tag
            self.canvas.SaveAs(output_file + ".png")

            if layer == 26:
                self.canvas.Print(self.output_pdf_file + ")", "pdf")
            else:
                self.canvas.Print(self.output_pdf_file, "pdf")
