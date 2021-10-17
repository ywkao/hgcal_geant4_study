#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <vector>

const int NLayers = 26;
bool flag_draw_plots = false;
TString global_tag = "dqm";
vector<TString> layer_tags = {"01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26"};

TFile *fin;
TString directory;
TCanvas *c = new TCanvas("c", "", 800, 600);

TGraphErrors* retrieve_histogram_info(TString stem, int color) //{{{
{
    vector<double> v_mean, v_stderr;

    // loop CEE layers
    for(int i=0; i<NLayers; ++i)
    {
        //TH1D* h = (TH1D*) fin->Get("DQMData/Run 1/HGCAL/Run summary/HGCalDigisV/HGCalEESensitive/ADC_layer_01");
        TString histname = directory + stem + "_layer_" + layer_tags[i];
        TH1D* h = (TH1D*) fin->Get(histname);
        
        double mean = h->GetMean();
        double stderr = h->GetStdDev();

        v_mean.push_back(mean);
        v_stderr.push_back(stderr);

        printf("mean = %.3f, ", mean);
        printf("stderr = %.3f\n", stderr);

        // skip making individual plots
        if(flag_draw_plots)
        {
            h->Draw();
            c->SaveAs("plots/dqm_D86_electron_E20_ADC_layer_01.png");
        }
    }

    // create graph
    double x[NLayers] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
    double ex[NLayers] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double *y = v_mean.data();
    double *ey = v_stderr.data();
    TGraphErrors *gr = new TGraphErrors(NLayers, x, y, ex, ey);
    gr->SetTitle("");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(color);
    gr->SetLineColor(color);

    return gr;
}
//}}}
void get_HGCAL_EE_info(TString type, TString stem) //{{{
{
    directory = "DQMData/Run 1/HGCAL/Run summary/" + type + "/HGCalEESensitive/";

    // retrieve information
    fin = TFile::Open("DQM_CloseByParticle_2026D86_Positron_E20_nEvents1000.root");
    TGraphErrors* gr_e20 = retrieve_histogram_info(stem, kBlack);

    fin = TFile::Open("DQM_CloseByParticle_2026D86_Positron_E100_nEvents1000.root");
    TGraphErrors* gr_e100 = retrieve_histogram_info(stem, kBlue);

    fin = TFile::Open("DQM_CloseByParticle_2026D86_Positron_E300_nEvents1000.root");
    TGraphErrors* gr_e300 = retrieve_histogram_info(stem, kRed);


    // make plots
    TString output = "plots/" + global_tag + "_" + type + "_" + stem + ".png";
    TLegend *legend = new TLegend(0.75, 0.70, 0.89, 0.86);

    legend->SetLineColor(0);
    legend->SetTextSize(0.03);
    legend->AddEntry(gr_e300 , "300 GeV" , "lp");
    legend->AddEntry(gr_e100 , "100 GeV" , "lp");
    legend->AddEntry(gr_e20  , "20 GeV"  , "lp");

    gr_e300->Draw("ALP");
    gr_e300->GetXaxis()->SetTitle("Layer");
    gr_e300->GetYaxis()->SetTitle(stem);

    gr_e100->Draw("LP;same");
    gr_e20->Draw("LP;same");
    legend->Draw("same");

    c->SaveAs(output);
}
//}}}
void read_DQM_files()
{
    c->SetLeftMargin(0.15);
    c->SetTicks(1,1);

    get_HGCAL_EE_info("HGCalDigisV", "ADC");
    get_HGCAL_EE_info("HGCalDigisV", "DigiOccupancy_Minus");
    get_HGCAL_EE_info("HGCalDigisV", "DigiOccupancy_Plus");

    get_HGCAL_EE_info("HGCalRecHitsV", "HitOccupancy_Minus");
    get_HGCAL_EE_info("HGCalRecHitsV", "HitOccupancy_Plus");
    get_HGCAL_EE_info("HGCalRecHitsV", "energy");

    get_HGCAL_EE_info("HGCalSimHitsV", "HitOccupancy_Minus");
    get_HGCAL_EE_info("HGCalSimHitsV", "HitOccupancy_Plus");
    get_HGCAL_EE_info("HGCalSimHitsV", "energy_time_0");
    get_HGCAL_EE_info("HGCalSimHitsV", "energy_time_1");
}
