#include "hits.C"

int main(int argc, char *argv[])
{
    TString input = argc<=1 ? "ntuple_496.root" : argv[1];
    TString tag   = argc<=2 ? "" : argv[2];

    //--------------------------------------------------
    // load & loop
    //--------------------------------------------------
    TTree *tree = 0;
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(input);
    if (!f || !f->IsOpen()) { f = new TFile(input); }
    TDirectory * dir = (TDirectory*)f->Get(input + ":/rechitntupler");
    dir->GetObject("hits",tree);

    hits mytree(tree);
    mytree.Loop();

    //--------------------------------------------------
    // make plots
    //--------------------------------------------------
    if(false) {
        TString output_heading = "./eos_output/testbeam_data_rechits_";
        mytree.MakePlot(output_heading);
    }

    //--------------------------------------------------
    // output
    //--------------------------------------------------
    TString path = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v17/";
    mytree.Report(input.ReplaceAll(path, ""));

    TFile *fout = new TFile("./output/output_" + tag + ".root", "RECREATE");
    fout->cd();
    mytree.Write();
    fout->Close();

    return 0;
}
