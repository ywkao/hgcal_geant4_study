#!/usr/bin/env python
import copy
import math
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(2210)
ROOT.gStyle.SetOptFit(1111)

import toolbox.plot_utils as pu

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
            if not i==0: continue # use E100 for check
            self.tag = self.tags[i]

            print ">>> rootfile", rootfile
            self.fin = ROOT.TFile.Open(rootfile, "R")
            tree = self.fin.Get("prodEE_DigiSim/tr_hits")
            vf.append(self.fin)
            vt.append(tree)

            #--------------------------------------------------
            # total event
            #--------------------------------------------------
            self.flag_make_plots = True
            self.check_individual_event = False 
            self.output_pdf_file = self.output_directory + "/output_hits_" + self.tags[i] + ".pdf"
            self.__retrieve_hit_plots(tree)

            break

            #--------------------------------------------------
            # specific event
            #--------------------------------------------------
            self.check_individual_event = True

            self.flag_make_plots = True
            for evtNo in range(0, 5):
                self.evtNo = evtNo 
                self.output_pdf_file = self.output_directory + "/output_hits_" + self.tags[i] + "_evtNo%d.pdf" % self.evtNo
                self.__retrieve_hit_plots(tree)

            break

            # False for evaluating efficiency
            self.flag_make_plots = False

            # 30 min processing time
            self.efficiency = []
            for layer in range(1, 27):
                h = ROOT.TH1D("eff_layer-%d" % layer, ";Efficiency (%);Entries", 50, 0, 100)
                self.efficiency.append(h)

            for evtNo in range(0, 1000):
                self.evtNo = evtNo 
                self.output_pdf_file = self.output_directory + "/output_hits_" + self.tags[i] + "_evtNo%d.pdf" % self.evtNo
                self.__retrieve_hit_plots(tree)

            for layer in range(1, 27):
                self.efficiency[layer-1].Draw()
                # ploting
                str_number = "0%d" % layer if layer < 10 else "%d" % layer
                output_file = self.output_directory + "/hit_efficiency_x_y_layer" + str_number + "_" + self.tag
                self.canvas.SaveAs(output_file + ".png")

            break

            # add energy selection
            self.output_pdf_file = self.output_directory + "/output_hits_" + self.tags[i] + ".pdf"
            self.__make_plots(tree)

            break

    #====================================================================================================

    def __make_plots(self, tree): #{{{
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

        e_set0 = ROOT.TH1D("e_set0", ";E_{set0} (MeV);Entries", 400, 0, 400)
        e_set0.SetTitle("")
        e_set0.SetLineWidth(2)
        e_set0.SetLineColor(ROOT.kBlack)
        e_set0.GetXaxis().SetRangeUser(0, 150)
        e_set0.GetXaxis().SetTitleOffset(1.1)

        h_Rxy_e = ROOT.TH2D("h_Rxy_e", ";E_{set0} (MeV);R_{xy} (cm)", 60, 40, 100, 600, 0, 300)

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
            if total_energy[evtNo] > 65.:
                e_set0.Fill(total_energy[evtNo])
            else:
                c.Fill(total_energy[evtNo])

        for i, hit in enumerate(tree):
            h_Rxy_e.Fill(total_energy[hit.evtNo], hit.r)
            h_eta.Fill(hit.eta)

            if total_energy[hit.evtNo] > 65.:
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
        output_file = self.output_directory + "/hits_Rxy_Z_" + self.tag
        self.canvas.SaveAs(output_file + ".png")
        self.canvas.Print(self.output_pdf_file + "(", "pdf")

        for i, h in enumerate(vh_x_y):
            h.Draw()
            b[i].SetMarkerColor(ROOT.kRed)
            b[i].Draw("same")

            layer = i+1
            str_number = "0%d" % layer if layer < 10 else "%d" % layer
            output_file = self.output_directory + "/hits_x_y_layer" + str_number + "_" + self.tag
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

    #}}}

    def __retrieve_hit_plots(self, tree):

        # output pdf only when checking single events -> save disk space
        self.flag_make_png_file = not self.check_individual_event
        self.flag_make_pdf_file = self.check_individual_event

        # expected hit positions
        tree_max_cell = self.fin.Get("prodEE_DigiSim/tr_max_cell")
        tree_energy_weighted = self.fin.Get("prodEE_DigiSim/tr_energy_weighted")
        tree_linear_track = self.fin.Get("prodEE_DigiSim/tr_linear_trajectory")

        #------------------------------
        # Rxy-Z distribution
        #------------------------------
        if self.flag_make_plots:
            legend = self.get_TLegend(0.15, 0.65, 0.45, 0.85)

            # Rec hits
            self.tree = tree
            self.show_rec_hits = True
            vh = self.draw_RZ_hitstogram()

            vh[0].Draw()
            vh[1].Draw("same")
            vh[2].Draw("same")

            if self.check_individual_event:
                legend.AddEntry(vh[0], "Rec. hits", "l")
            else:
                legend.AddEntry(vh[0], self.wafer_types[0], "l")
                legend.AddEntry(vh[1], self.wafer_types[1], "l")
                legend.AddEntry(vh[2], self.wafer_types[2], "l")
                legend.Draw("same")

            # expected hit position
            if self.check_individual_event:
                self.show_rec_hits = False 

                self.algo_type = "max_cell"
                self.tree = tree_max_cell
                self.settings = (20, 0.3, ROOT.kBlue)
                h1 = self.draw_RZ_hitstogram()
                h1.Draw("same")

                self.algo_type = "energy_weighted"
                self.tree = tree_energy_weighted
                self.settings = (20, 0.3, ROOT.kGreen)
                h2 = self.draw_RZ_hitstogram()
                h2.Draw("same")

                self.algo_type = "linear_track"
                self.tree = tree_linear_track
                self.settings = (20, 0.3, ROOT.kRed)
                h3 = self.draw_RZ_hitstogram()
                h3.Draw("same")

                legend.AddEntry(h1, "Expected from max cell", "l")
                legend.AddEntry(h2, "Expected from energy weight", "l")
                legend.AddEntry(h3, "Expected from linear track", "l")
                legend.Draw("same")

            self.canvas.Update()

            output_file = self.output_directory + "/hits_Rxy_Z_" + self.tag
            if self.flag_make_pdf_file: self.canvas.Print(self.output_pdf_file + "(", "pdf")
            if self.flag_make_png_file: self.canvas.SaveAs(output_file + ".png")

        #------------------------------
        # X-Y distributions
        #------------------------------
        v_eff, v_eff_unc = [], []
        for layer in range(1,27):
            self.layer = layer
            self.title = "Layer-0%d" % layer if layer < 10 else "Layer-%d" % layer
            #self.selection = "evtNo==%d && layerNo==%d" % (self.evtNo, layer) if self.check_individual_event else "layerNo==%d" % layer
            pass_signal_region = "(signal_region_linear_track==1 || signal_region_linear_track==2)"

            # Rec hits
            self.tree = tree
            self.show_rec_hits = True

            self.name = "h_2d_hits_layer%d_outside_signal_region" % layer
            self.selection = "evtNo==%d && layerNo==%d && !%s"  % (self.evtNo, layer, pass_signal_region) if self.check_individual_event else "layerNo==%d" % layer
            self.settings = (20, 0.3, ROOT.kGray+2)
            hb = self.draw_XY_hitstogram()

            self.name = "h_2d_hits_layer%d_in_signal_region" % layer
            self.selection = "evtNo==%d && layerNo==%d && %s"  % (self.evtNo, layer, pass_signal_region) if self.check_individual_event else "layerNo==%d" % layer
            self.settings = (20, 0.3, ROOT.kRed-7)
            hs = self.draw_XY_hitstogram()

            ns, nb = hs.GetEntries(), hb.GetEntries()
            if ns+nb > 0.:
                #print ns, nb, ns/(nb+ns)
                efficiency = ns/(nb+ns)
                unc = math.sqrt( (1.-efficiency)*efficiency / (nb+ns) )
                v_eff.append(efficiency)
                v_eff_unc.append(unc)

                efficiency = 100 * ns/(nb+ns)
                content = "Efficiency = %.1f%%" % efficiency
            else:
                #print ns, nb, "N/A"
                efficiency = -999
                v_eff.append(-1)
                v_eff_unc.append(0)
                content = "Efficiency = N/A"

            if not self.flag_make_plots: 
                self.efficiency[layer-1].Fill(efficiency)
                continue

            hb.Draw()
            hs.Draw("same")

            latex = ROOT.TLatex()
            latex.SetNDC()
            latex.SetTextSize(0.04)
            latex.SetTextAlign(22)
            latex.SetTextFont(42)
            latex.DrawLatex(0.5, 0.2, content)

            # expected hit position
            if self.check_individual_event:
                self.show_rec_hits = False 
                self.name = "h_2d_hits_layer%d" % layer

                self.algo_type = "max_cell"
                self.tree = tree_max_cell
                self.settings = (29, 0.5, ROOT.kBlue)
                h1 = self.draw_XY_hitstogram()
                h1.Draw("same")

                self.algo_type = "energy_weighted"
                self.tree = tree_energy_weighted
                self.settings = (29, 0.5, ROOT.kGreen)
                h2 = self.draw_XY_hitstogram()
                h2.Draw("same")

                self.algo_type = "linear_track"
                self.tree = tree_linear_track
                self.settings = (29, 0.5, ROOT.kRed)
                h3 = self.draw_XY_hitstogram()
                h3.Draw("same")

                #legend = ROOT.TLegend(0.20, 0.20, 0.65, 0.35)
                #legend.SetTextSize(0.03)
                #legend.AddEntry(hs,  "Rec. hits in SR", "p")
                #legend.AddEntry(hb,  "Rec. hits not in SR", "p")
                #legend.AddEntry(h1, "Expected from max cell", "p")
                #legend.AddEntry(h2, "Expected from energy weight", "p")
                #legend.AddEntry(h3, "Expected from linear track", "p")
                #legend.Draw("same")

            # ploting
            str_number = "0%d" % layer if layer < 10 else "%d" % layer
            output_file = self.output_directory + "/hits_x_y_layer" + str_number + "_" + self.tag

            if self.flag_make_pdf_file:
                if layer == 26:
                    self.canvas.Print(self.output_pdf_file + ")", "pdf")
                else:
                    self.canvas.Print(self.output_pdf_file, "pdf")

            if self.flag_make_png_file:
                self.canvas.SaveAs(output_file + ".png")

        for layer in range(1,27):
            idx = layer-1
            print v_eff[idx], v_eff_unc[idx]

        output_file = self.output_directory + "/efficiency_profile" + self.tag + "_evtNo%d" % self.evtNo
        gr = pu.get_graph_from_list("Layer depth [ X_{0} ]", "Efficiency", pu.x_D86, v_eff, [0.]*26, v_eff_unc, ROOT.kBlack)
        gr.Draw("alp")
        self.canvas.SaveAs(output_file + ".png")

    #====================================================================================================

    def draw_XY_hitstogram(self):
        fullname = self.name if self.show_rec_hits else self.name + "_" + self.algo_type
        h = ROOT.TH2D(fullname, self.title, 680, -170, 170, 680, -170, 170 )

        if self.show_rec_hits:
            self.tree.Draw("y:x>>%s" % self.name, self.selection)
            h = ROOT.gPad.GetPrimitive(self.name)
            #print fullname, h.GetEntries(), h.GetMean(1), h.GetRMS(1)

        else: # expected position
            for i, hit in enumerate(self.tree):
                if not i == self.evtNo: continue
                x = hit.vx[self.layer-1]
                y = hit.vy[self.layer-1]
                h.Fill(x,y);
                break

        style, size, color = self.settings
        h.SetStats(0)
        h.SetMarkerStyle(style)
        h.SetMarkerSize(size)
        h.SetMarkerColor(color)
        h.GetXaxis().SetTitle("X (cm)")
        h.GetYaxis().SetTitle("Y (cm)")
        h.GetXaxis().SetTitleOffset(1.1)
        h.GetYaxis().SetTitleOffset(1.0)

        return h


    def draw_RZ_hitstogram(self):
        if self.show_rec_hits:
            vh, colors = [], []
            if self.check_individual_event:
                colors = [ROOT.kBlack, ROOT.kBlack, ROOT.kBlack]
            else:
                colors = [ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta]

            self.wafer_types = ["is_Silicon_w120", "is_Silicon_w200", "is_Silicon_w300"]
            for i, wafer in enumerate(self.wafer_types):
                name, title = "h_hits_Rxy_Z_%s" % wafer, "Hits R_{xy} vs Z"
                selection = "evtNo==%d && %s==1" % (self.evtNo, wafer) if self.check_individual_event else "%s==1" % wafer
                h = ROOT.TH2D(name, title, 160, 300, 380, 600, 0, 300 )
                self.tree.Draw("r:z>>%s" % name, selection)
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

            return vh

        else: # expected position
            name, title = "h_hits_Rxy_Z_%s" % self.algo_type, "Hits R_{xy} vs Z"
            h = ROOT.TH2D(name, title, 160, 300, 380, 600, 0, 300 )
            for i, hit in enumerate(self.tree):
                if not i == self.evtNo: continue
                for layer in range(1,27):
                    x = hit.vx[layer-1]
                    y = hit.vy[layer-1]
                    z = hit.vz[layer-1]
                    r = math.sqrt(pow(x,2)+pow(y,2))
                    h.Fill(z,r);
                break

            style, size, color = self.settings
            h.SetStats(0)
            h.SetLineWidth(2)
            h.SetLineColor(color)
            h.SetMarkerStyle(style)
            h.SetMarkerSize(size)
            h.SetMarkerColor(color)
            h.GetXaxis().SetTitle("Z (cm)")
            h.GetYaxis().SetTitle("R_{xy} (cm)")
            h.GetXaxis().SetTitleOffset(1.1)
            h.GetYaxis().SetTitleOffset(1.2)

            return h


    def get_TLegend(self, x0, y0, x1, y1):
        #legend = ROOT.TLegend(0.15, 0.65, 0.45, 0.85)
        legend = ROOT.TLegend(x0, y0, x1, y1)
        legend.SetLineColor(0)
        legend.SetTextSize(0.04)
        return legend

