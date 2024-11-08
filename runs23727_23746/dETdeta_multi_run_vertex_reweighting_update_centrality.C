#include <iostream>
#include <TH2D.h>
#include <TH1D.h>
#include <TChain.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <sstream> //std::ostringstsream
#include <fstream> //std::ifstream
#include <iostream> //std::cout, std::endl
#include <cmath>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <string>
#include <set>
#include <TVector3.h>
#include <map>
#include <vector>
#include <TDatabasePDG.h>
#include <tuple>

using namespace std;

void dETdeta_multi_run_vertex_reweighting_update_centrality(const char* generator) {	

    TFile *file2 = TFile::Open("/sphenix/user/dlis/Projects/centrality/calib/mbdana_centrality_23696.root", "READ");
    if (!file2 || file2->IsZombie()) {
        std::cerr << "Error opening file2!" << std::endl;
        return 1;
    }

    TTree *tree2 = (TTree*)file2->Get("tn_centrality");
    if (!tree2) {
        std::cerr << "Error getting TTree from file2!" << std::endl;
        return 1;
    }

    float high;
    tree2->SetBranchAddress("high", &high);
	float new_data_centrality[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	Long64_t nentries2 = tree2->GetEntries();
    for (Long64_t i = 0; i < nentries2; i++) {
        tree2->GetEntry(i);
        if (i % 5 == 0 && i < 100) {
        	new_data_centrality[19-(i/5)] = high;
        }
    }
    new_data_centrality[0] = 0.0;

    file2->Close();

    TFile *out = new TFile(TString::Format("dETdeta_vertex_reweight_run23727_23745_%s_p015.root",generator),"UPDATE");
    TTree *tree = (TTree*)out->Get("T");
    if (!tree) {
        std::cerr << "Error getting TTree!" << std::endl;
        return 1;
    }

    TBranch *newBranch = tree->Branch("new_data_centrality", new_data_centrality, "new_data_centrality[20]/F");
    newBranch->Fill();

    out->cd();
    // Write the updated TTree to the file
    tree->Write("", TObject::kOverwrite);
	out->Close();
}
