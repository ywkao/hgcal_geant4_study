#include "hits.C"

int main(int argc, char *argv[])
{
    TString input = argc<=1 ? "ntuple_496.root" : argv[1];
    printf(">>> %s: ", input.Data());

    TTree *tree = 0;
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(input);
    if (!f || !f->IsOpen()) { f = new TFile(input); }
    TDirectory * dir = (TDirectory*)f->Get(input + ":/rechitntupler");
    dir->GetObject("hits",tree);

    hits mytree(tree);
    mytree.Loop();

    return 0;
}
