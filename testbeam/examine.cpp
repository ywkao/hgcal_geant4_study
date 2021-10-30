#include "hits.C"

int main(int argc, char *argv[])
{
    TString input = argc<=1 ? "ntuple_496.root" : argv[1];
    TString tag   = argc<=2 ? "" : argv[2];
    //printf(">>> %s: \n", input.Data());

    TTree *tree = 0;
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(input);
    if (!f || !f->IsOpen()) { f = new TFile(input); }
    TDirectory * dir = (TDirectory*)f->Get(input + ":/rechitntupler");
    dir->GetObject("hits",tree);

    hits mytree(tree);
    mytree.Loop();

    TString output = "/eos/user/y/ykao/www/HGCAL_Geant4_project/test/testbeam_data_rechits_layer" + tag + ".png";
    mytree.MakePlot(output);

    //TString path = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v17/";
    //mytree.Report(input.ReplaceAll(path, ""));

    return 0;
}
