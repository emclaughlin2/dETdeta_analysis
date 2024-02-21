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

std::set<std::tuple<int, int>> emcal_hot_dead_map;
std::set<std::tuple<int, int>> ihcal_hot_dead_map;
std::set<std::tuple<int, int>> ohcal_hot_dead_map;

float eta_spread_ihcal[24] = {0.0};
float eta_spread_ohcal[24] = {0.0};
int eta_spread_ihcal_count[24] = {0};
int eta_spread_ohcal_count[24] = {0};

void dETdeta_hot_dead_map(int runnumber = 23727, string optflag = "", int minus_z = -2, int plus_z = 2) {  

    if (optflag != "") optflag += "_";
    TFile *out = new TFile(Form("run%d_hotdeadmap_%s%d_%d.root", runnumber, optflag.c_str(), minus_z, plus_z),"RECREATE");

    vector<int> emcal_hot_dead_ieta;
    vector<int> emcal_hot_dead_iphi;
    vector<int> ihcal_hot_dead_ieta;
    vector<int> ihcal_hot_dead_iphi;
    vector<int> ohcal_hot_dead_ieta;
    vector<int> ohcal_hot_dead_iphi;
    vector<float> ihcal_eta_bin_centers;
    vector<float> ohcal_eta_bin_centers;

    TTree* T = new TTree("T", "");
    T->Branch("emcal_hot_dead_ieta", &emcal_hot_dead_ieta);
    T->Branch("emcal_hot_dead_iphi", &emcal_hot_dead_iphi);
    T->Branch("ihcal_hot_dead_ieta", &ihcal_hot_dead_ieta);
    T->Branch("ihcal_hot_dead_iphi", &ihcal_hot_dead_iphi);
    T->Branch("ohcal_hot_dead_ieta", &ohcal_hot_dead_ieta);
    T->Branch("ohcal_hot_dead_iphi", &ohcal_hot_dead_iphi);
    T->Branch("ihcal_eta_bin_centers", &ihcal_eta_bin_centers);
    T->Branch("ohcal_eta_bin_centers", &ohcal_eta_bin_centers);

    TH2F* h_2D_hot_dead_ihcal= new TH2F("h_2D_hot_dead_ihcal","",24,0.,24.,64,0.,64.);
    TH2F* h_2D_hot_dead_ohcal = new TH2F("h_2D_hot_dead_ohcal","",24,0.,24.,64,0.,64.);
    TH2F* h_2D_hot_dead_emcal = new TH2F("h_2D_hot_dead_emcal","",96,0.,96.,256,0.,256.);

    TChain chain("outt");
    TChain chain2("ttree");

    int hotmap[3][96][256];
    chain.SetBranchAddress("hotmap", hotmap);

    float m_vtx[3];
    int m_simtwrmult_ihcal;
    int m_simtwr_ihcal_ieta[1536];
    float m_simtwr_ihcal_eta[1536];
    int m_simtwrmult_ohcal;
    int m_simtwr_ohcal_ieta[1536];
    float m_simtwr_ohcal_eta[1536];
    chain2.SetBranchAddress("track_vtx", m_vtx);
    chain2.SetBranchAddress("sectorih", &m_simtwrmult_ihcal);
    chain2.SetBranchAddress("ihcaletabin", m_simtwr_ihcal_ieta);
    chain2.SetBranchAddress("ihetacor", m_simtwr_ihcal_eta);
    chain2.SetBranchAddress("sectoroh", &m_simtwrmult_ohcal);
    chain2.SetBranchAddress("ohcaletabin", m_simtwr_ohcal_ieta);
    chain2.SetBranchAddress("ohetacor", m_simtwr_ohcal_eta);
    
    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
    for (int i = 0; i < 20; i++) {
        TString wildcardPath = TString::Format("%sevents_p008_zs_%d_%sdata_cor_%d.root", inputDirectory, runnumber, optflag.c_str(), i);
        if (i == 0) std::cout << wildcardPath << std::endl;
        chain.Add(wildcardPath);
        chain2.Add(wildcardPath);
    }

    Long64_t nEntries = chain.GetEntries();
    std::cout << "hot tower tree entries " << nEntries << std::endl;

    for (Long64_t entry = 0; entry < nEntries; entry++) {
        chain.GetEntry(entry);
        
        for (int j = 0; j < 96; j++) {
            for (int k = 0; k < 256; k++) {
                if (hotmap[0][j][k] == 1) {
                    tuple<int, int> new_hot_tower = make_tuple(j,k);
                    emcal_hot_dead_map.insert(new_hot_tower);
                    h_2D_hot_dead_emcal->Fill(j,k);
                }
                if (hotmap[1][j][k] == 1 && j < 24 && k < 64) {
                    tuple<int, int> new_hot_tower = make_tuple(j,k);
                    ihcal_hot_dead_map.insert(new_hot_tower);
                    h_2D_hot_dead_ihcal->Fill(j,k);
                }
                if (hotmap[2][j][k] == 1 && j < 24 && k < 64) {
                    tuple<int, int> new_hot_tower = make_tuple(j,k);
                    ohcal_hot_dead_map.insert(new_hot_tower);
                    h_2D_hot_dead_ohcal->Fill(j,k);
                }
            }
        }
    }

    for (const auto& tuple : emcal_hot_dead_map) {
        h_2D_hot_dead_emcal->Fill(get<0>(tuple), get<1>(tuple));
        emcal_hot_dead_ieta.push_back(get<0>(tuple));
        emcal_hot_dead_iphi.push_back(get<1>(tuple));
    }

    for (const auto& tuple : ihcal_hot_dead_map) {
        h_2D_hot_dead_ihcal->Fill(get<0>(tuple), get<1>(tuple));
        ihcal_hot_dead_ieta.push_back(get<0>(tuple));
        ihcal_hot_dead_iphi.push_back(get<1>(tuple));
    }

    for (const auto& tuple : ohcal_hot_dead_map) {
        h_2D_hot_dead_ohcal->Fill(get<0>(tuple), get<1>(tuple));
        ohcal_hot_dead_ieta.push_back(get<0>(tuple));
        ohcal_hot_dead_iphi.push_back(get<1>(tuple));
    }


    nEntries = chain2.GetEntries();
    std::cout << "tower tree entries " << nEntries << std::endl;

    for (Long64_t entry = 0; entry < nEntries; entry++) {
        chain2.GetEntry(entry);

        if (isnan(m_vtx[2])) { continue; }
        if (m_vtx[2] < minus_z || m_vtx[2] > plus_z) { continue; }
        
        for (int i = 0; i < m_simtwrmult_ihcal; i++) {
            eta_spread_ihcal[m_simtwr_ihcal_ieta[i]] += m_simtwr_ihcal_eta[i];
            eta_spread_ihcal_count[m_simtwr_ihcal_ieta[i]] += 1;
        }

        for (int i = 0; i < m_simtwrmult_ohcal; i++) {
            eta_spread_ohcal[m_simtwr_ohcal_ieta[i]] += m_simtwr_ohcal_eta[i];
            eta_spread_ohcal_count[m_simtwr_ohcal_ieta[i]] += 1;
        }

    }

    for (int i = 0; i < 24; i++) {
        std::cout << "ihcal " << eta_spread_ihcal[i]/eta_spread_ihcal_count[i] << std::endl;
        std::cout << "ohcal " << eta_spread_ohcal[i]/eta_spread_ohcal_count[i] << std::endl;
    }

    for (int i = 0; i < 24; i++) {
        ihcal_eta_bin_centers.push_back(eta_spread_ihcal[i]/eta_spread_ihcal_count[i]);
        ohcal_eta_bin_centers.push_back(eta_spread_ohcal[i]/eta_spread_ohcal_count[i]);
    }


    T->Fill();
    
    out->Write();
    out->Close();
}
