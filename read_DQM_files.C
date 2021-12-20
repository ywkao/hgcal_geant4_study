#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <vector>

bool flag_draw_plots = false;
bool flag_debug = false;

TString global_tag = "dqm";
TString dir_source = "/eos/user/y/ykao/www/HGCAL_Geant4_project/electronBeam/";

TFile *fin;
vector<TString> input_files;
TString directory;
bool make_resolution;
TCanvas *c = new TCanvas("c", "", 800, 600);

TGraphErrors* retrieve_histogram_info(TString stem, int color) //{{{
{
    vector<double> v_mean, v_stderr;
    const int NLayers = global_tag.Contains("D86") ? 26 : 28;
    vector<TString> layer_tags_D86 = {"01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26"};
    vector<TString> layer_tags_D83 = {"01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28"};
    vector<TString> layer_tags = global_tag.Contains("D86") ? layer_tags_D86 : layer_tags_D83;

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

        if(flag_debug)
        {
            printf("mean = %.3f, ", mean);
            printf("stderr = %.3f\n", stderr);
        }

        // skip making individual plots
        if(flag_draw_plots)
        {
            h->Draw();
            c->SaveAs("plots/dqm_D86_electron_E20_ADC_layer_01.png");
        }
    }

    // create graph
    // /afs/cern.ch/work/a/apsallid/public/HGCal/Validation/D8[36].cumulative.xo
    double x_D86[26]  = {0.564,1.567,2.547,3.549,4.528,5.531,6.509,7.512,8.49,9.493,10.472,11.474,12.453,
                         13.455,14.434,15.437,16.415,17.418,18.975,19.978,21.536,22.538,24.096,25.099,26.656,27.659};
    double x_D83[28]  = {0.441,1.438,2.289,3.286,4.137,5.134,5.985,6.982,7.833,8.83,9.68,10.676,11.528,12.524,13.376,
                         14.373,15.224,16.221,17.072,18.069,18.92,19.917,20.768,21.765,22.616,23.613,24.464,25.461};
    double ex_D86[26] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double ex_D83[28] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    double *x = global_tag.Contains("D86") ? x_D86 : x_D83;
    double *ex = global_tag.Contains("D86") ? ex_D86 : ex_D83;
    double *y = v_mean.data();
    double *ey = v_stderr.data();

    if( !make_resolution )
        for(int i=0; i<NLayers; ++i) ey[i] = 0.;

    //printf("check ey[0] = %.2f\n", ey[0]);

    TGraphErrors *gr = new TGraphErrors(NLayers, x, y, ex, ey);
    gr->SetTitle("");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(color);
    gr->SetLineColor(color);

    return gr;
}
//}}}
void get_HGCAL_EE_info(TString type, TString stem, bool option_resolution) //{{{
{
    make_resolution = option_resolution;
    directory = "DQMData/Run 1/HGCAL/Run summary/" + type + "/HGCalEESensitive/";
    //printf("input_files[0] = %s\n", input_files[0].Data());
    //printf("input_files[1] = %s\n", input_files[1].Data());
    //printf("input_files[2] = %s\n", input_files[2].Data());

    // retrieve information
    fin = TFile::Open( dir_source + input_files[0] );
    TGraphErrors* gr_e20 = retrieve_histogram_info(stem, kBlack);

    fin = TFile::Open( dir_source + input_files[1] );
    TGraphErrors* gr_e100 = retrieve_histogram_info(stem, kBlue);

    fin = TFile::Open( dir_source + input_files[2] );
    TGraphErrors* gr_e300 = retrieve_histogram_info(stem, kRed);


    // make plots
    TString output = dir_source + "plots_v3/" + global_tag + "_" + type + "_" + stem + ".png";
    TLegend *legend = new TLegend(0.72, 0.70, 0.86, 0.86);

    legend->SetLineColor(0);
    legend->SetTextSize(0.04);
    legend->AddEntry(gr_e300 , "300 GeV" , "lp");
    legend->AddEntry(gr_e100 , "100 GeV" , "lp");
    legend->AddEntry(gr_e20  , "20 GeV"  , "lp");

    gr_e300->Draw("ALP");
    gr_e300->GetXaxis()->SetTitle("Layer depth [^{}X_{0} ]");
    gr_e300->GetXaxis()->SetTitleOffset(1.2);
    gr_e300->GetYaxis()->SetTitle(stem);

    gr_e100->Draw("LP;same");
    gr_e20->Draw("LP;same");
    legend->Draw("same");

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextFont(43);
    latex->SetTextAlign(11);
    latex->SetTextSize(22);
    latex->DrawLatex( 0.15, 0.92, global_tag );

    c->SaveAs(output);
}
//}}}
void make_plots() //{{{
{
    get_HGCAL_EE_info("HGCalSimHitsV", "HitOccupancy_Plus", true);
    get_HGCAL_EE_info("HGCalSimHitsV", "energy_time_0", false);

    /*
    get_HGCAL_EE_info("HGCalDigisV", "ADC", false);
    get_HGCAL_EE_info("HGCalDigisV", "DigiOccupancy_Minus", true);
    get_HGCAL_EE_info("HGCalDigisV", "DigiOccupancy_Plus", true);

    get_HGCAL_EE_info("HGCalRecHitsV", "HitOccupancy_Minus", true);
    get_HGCAL_EE_info("HGCalRecHitsV", "HitOccupancy_Plus", true);
    get_HGCAL_EE_info("HGCalRecHitsV", "energy", false);

    get_HGCAL_EE_info("HGCalSimHitsV", "HitOccupancy_Minus", true);
    get_HGCAL_EE_info("HGCalSimHitsV", "HitOccupancy_Plus", true);
    get_HGCAL_EE_info("HGCalSimHitsV", "energy_time_0", false);
    get_HGCAL_EE_info("HGCalSimHitsV", "energy_time_1", false);
    */
}
//}}}
void read_DQM_files() //{{{
{
    c->SetLeftMargin(0.15);
    c->SetTicks(1,1);

    global_tag = "dqm_D86_R80To100";  input_files = {"DQM_CloseByParticle_2026D86_Positron_E20_nEvents1000.root","DQM_CloseByParticle_2026D86_Positron_E100_nEvents1000.root","DQM_CloseByParticle_2026D86_Positron_E300_nEvents1000.root"};
    make_plots();

    global_tag = "dqm_D86_R35To60";   input_files = {"DQM_CloseByParticle_positron_D86_R35To60_E20.root","DQM_CloseByParticle_positron_D86_R35To60_E100.root","DQM_CloseByParticle_positron_D86_R35To60_E300.root"};
    make_plots();

    global_tag = "dqm_D86_R120To140"; input_files = {"DQM_CloseByParticle_positron_D86_R120To140_E20.root","DQM_CloseByParticle_positron_D86_R120To140_E100.root","DQM_CloseByParticle_positron_D86_R120To140_E300.root"};
    make_plots();

    global_tag = "dqm_D83_R35To60";   input_files = {"DQM_CloseByParticle_positron_D83_R35To60_E20.root","DQM_CloseByParticle_positron_D83_R35To60_E100.root","DQM_CloseByParticle_positron_D83_R35To60_E300.root"};
    make_plots();

    global_tag = "dqm_D83_R120To140"; input_files = {"DQM_CloseByParticle_positron_D83_R120To140_E20.root","DQM_CloseByParticle_positron_D83_R120To140_E100.root","DQM_CloseByParticle_positron_D83_R120To140_E300.root"};
    make_plots();
} //}}}
