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

double emcal_eta_bin_centers[24];
double ihcal_eta_bin_centers[24];
double ohcal_eta_bin_centers[24];
double calo_eta_bin_centers[24];
double hcal_eta_bin_centers[24];

double vertex_reweight[200] = {1.0};
std::vector<float> centrality_bin;

const double eta_bin_centers[24] = {-1.05417,-0.9625,-0.870833,-0.779167,-0.6875,-0.595833,-0.504167,-0.4125,-0.320833,-0.229167,
	-0.1375,-0.0458333,0.0458333,0.1375,0.229167,0.320833,0.4125,0.504167,0.595833,0.6875,0.779167,0.870833,0.9625,1.05417};

float em_zs_calib[96][256];
float ih_zs_calib[24][64];
float oh_zs_calib[24][64];

const int cemcSize = 24576;
const int ihcalSize = 1536;
const int ohcalSize = 1536;
const int mbdSize = 256;
const int g4Size = 100000;
const int vtxSize = 3;

void fill_hot_dead_map_eta_bin_centers() {

	// hot dead maps 
    vector<int> *emcal_hot_dead_ieta = nullptr;
    vector<int> *emcal_hot_dead_iphi = nullptr;
    vector<int> *ihcal_hot_dead_ieta = nullptr;
    vector<int> *ihcal_hot_dead_iphi = nullptr;
    vector<int> *ohcal_hot_dead_ieta = nullptr;
    vector<int> *ohcal_hot_dead_iphi = nullptr;
    vector<float> *ihcal_eta_bins = nullptr;
    vector<float> *ohcal_eta_bins = nullptr;

	TFile *hotdeadfile = new TFile("/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/run54912_hotdeadmap_z_-10_10_ana450_2024p009_fixed_build.root", "READ"); 
	TTree *hotdeadtree = dynamic_cast<TTree*>(hotdeadfile->Get("T"));
	hotdeadtree->SetBranchAddress("emcal_hot_dead_ieta", &emcal_hot_dead_ieta);
    hotdeadtree->SetBranchAddress("emcal_hot_dead_iphi", &emcal_hot_dead_iphi);
    hotdeadtree->SetBranchAddress("ihcal_hot_dead_ieta", &ihcal_hot_dead_ieta);
    hotdeadtree->SetBranchAddress("ihcal_hot_dead_iphi", &ihcal_hot_dead_iphi);
    hotdeadtree->SetBranchAddress("ohcal_hot_dead_ieta", &ohcal_hot_dead_ieta);
    hotdeadtree->SetBranchAddress("ohcal_hot_dead_iphi", &ohcal_hot_dead_iphi);
    hotdeadtree->SetBranchAddress("ihcal_eta_bin_centers", &ihcal_eta_bins);
    hotdeadtree->SetBranchAddress("ohcal_eta_bin_centers", &ohcal_eta_bins);
	hotdeadtree->GetEntry(0);

	for (int i = 0; i < emcal_hot_dead_ieta->size(); i++) {
		tuple<int, int> new_hot_tower = make_tuple((*emcal_hot_dead_ieta)[i], (*emcal_hot_dead_iphi)[i]); 
        emcal_hot_dead_map.insert(new_hot_tower);
	}

	for (int i = 0; i < ihcal_hot_dead_ieta->size(); i++) {
		tuple<int, int> new_hot_tower = make_tuple((*ihcal_hot_dead_ieta)[i], (*ihcal_hot_dead_iphi)[i]); 
        ihcal_hot_dead_map.insert(new_hot_tower);
	}

	for (int i = 0; i < ohcal_hot_dead_ieta->size(); i++) {
		tuple<int, int> new_hot_tower = make_tuple((*ohcal_hot_dead_ieta)[i], (*ohcal_hot_dead_iphi)[i]); 
        ohcal_hot_dead_map.insert(new_hot_tower);
	}

	for (int i = 0; i < ihcal_eta_bins->size(); i++) {
		ihcal_eta_bin_centers[i] = (*ihcal_eta_bins)[i];
		emcal_eta_bin_centers[i] = (*ihcal_eta_bins)[i];
		calo_eta_bin_centers[i] = (*ihcal_eta_bins)[i];
		hcal_eta_bin_centers[i] = (*ihcal_eta_bins)[i];
	}
	for (int i = 0; i < ohcal_eta_bins->size(); i++) {
		ohcal_eta_bin_centers[i] = (*ohcal_eta_bins)[i];
	}

	hotdeadfile->Close();

}

void dETdeta_analysis_pedestal_subtraction(int runnumber = 23727, int reweighted_bins = 1, const char* opt_tag = "") {

	string filename = "dETdeta_analysis_pedestal_subtraction_final_calib";
	filename += "_" + to_string(runnumber);
	string opttag = opt_tag;
	if (strcmp(opt_tag,"")) { filename += "_" + opttag; }
  	string weighttag = (reweighted_bins?"reweighted_bins":"default_bins");
  	filename += "_" + weighttag + ".root";

	fill_hot_dead_map_eta_bin_centers();

	TFile *out = new TFile(filename.c_str(),"RECREATE");
	
	TH2F* h_2D_ihcal_calib = new TH2F("h_2D_ihcal_calib","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_ohcal_calib = new TH2F("h_2D_ohcal_calib","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_emcal_calib = new TH2F("h_2D_emcal_calib","",96,0.,96.,256,0.,256.);

  	TH2F* h_2D_ihcal_calibT = new TH2F("h_2D_ihcal_calibT","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_ohcal_calibT = new TH2F("h_2D_ohcal_calibT","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_emcal_calibT = new TH2F("h_2D_emcal_calibT","",96,0.,96.,256,0.,256.);
	
	TH1F* h_event_energy = new TH1F("h_event_energy","", 6000,-1000,5000);
	TH1F* h_event_hcal_energy = new TH1F("h_event_hcal_energy","", 6000,-1000,5000);
	TH1F* h_event_emcal_energy = new TH1F("h_event_emcal_energy","", 6000,-1000,5000);
	TH1F* h_event_ihcal_energy = new TH1F("h_event_ihcal_energy","", 2000,-1000,1000);
	TH1F* h_event_ohcal_energy = new TH1F("h_event_ohcal_energy","", 2000,-1000,1000);

	TH1F* h_emcal = new TH1F("h_emcal","",1500,-5,10);
	TH1F* h_ihcal = new TH1F("h_ihcal","",1500,-5,10);
	TH1F* h_ohcal = new TH1F("h_ohcal","",1500,-5,10);
	
	int emcal_num_bins = 24;
	double emcal_bin_edges[25];
    if (reweighted_bins) {
    	emcal_bin_edges[0] = emcal_eta_bin_centers[0] - 0.5 * (emcal_eta_bin_centers[1] - emcal_eta_bin_centers[0]);
	    for (int i = 1; i < emcal_num_bins; ++i) { emcal_bin_edges[i] = (emcal_eta_bin_centers[i] + emcal_eta_bin_centers[i - 1]) / 2.0; }
	    emcal_bin_edges[emcal_num_bins] = emcal_eta_bin_centers[emcal_num_bins - 1] + 0.5 * (emcal_eta_bin_centers[emcal_num_bins - 1] - emcal_eta_bin_centers[emcal_num_bins - 2]);
    } else {
    	emcal_bin_edges[0] = eta_bin_centers[0] - 0.5 * (eta_bin_centers[1] - eta_bin_centers[0]);
    	for (int i = 1; i < emcal_num_bins; ++i) { emcal_bin_edges[i] = (eta_bin_centers[i] + eta_bin_centers[i - 1]) / 2.0; }
    	emcal_bin_edges[emcal_num_bins] = eta_bin_centers[emcal_num_bins - 1] + 0.5 * (eta_bin_centers[emcal_num_bins - 1] - eta_bin_centers[emcal_num_bins - 2]);
    }
    TProfile* h_eT_eta_emcal_profile = new TProfile("h_eT_eta_emcal_profile","",emcal_num_bins, emcal_bin_edges);
	TH1F* h_eT_eta_emcal_profile_hist = new TH1F("h_eT_eta_emcal_profile_hist","",emcal_num_bins, emcal_bin_edges);
	TProfile* h_eta_emcal_profile = new TProfile("h_eta_emcal_profile","",emcal_num_bins, emcal_bin_edges);

	int ihcal_num_bins = 24;
    double ihcal_bin_edges[25];
    if (reweighted_bins) {
    	ihcal_bin_edges[0] = ihcal_eta_bin_centers[0] - 0.5 * (ihcal_eta_bin_centers[1] - ihcal_eta_bin_centers[0]);
	    for (int i = 1; i < ihcal_num_bins; ++i) { ihcal_bin_edges[i] = (ihcal_eta_bin_centers[i] + ihcal_eta_bin_centers[i - 1]) / 2.0; }
	    ihcal_bin_edges[ihcal_num_bins] = ihcal_eta_bin_centers[ihcal_num_bins - 1] + 0.5 * (ihcal_eta_bin_centers[ihcal_num_bins - 1] - ihcal_eta_bin_centers[ihcal_num_bins - 2]);
    } else {
    	ihcal_bin_edges[0] = eta_bin_centers[0] - 0.5 * (eta_bin_centers[1] - eta_bin_centers[0]);
    	for (int i = 1; i < ihcal_num_bins; ++i) { ihcal_bin_edges[i] = (eta_bin_centers[i] + eta_bin_centers[i - 1]) / 2.0; }
    	ihcal_bin_edges[ihcal_num_bins] = eta_bin_centers[ihcal_num_bins - 1] + 0.5 * (eta_bin_centers[ihcal_num_bins - 1] - eta_bin_centers[ihcal_num_bins - 2]);
    }
    TProfile* h_eT_eta_ihcal_profile = new TProfile("h_eT_eta_ihcal_profile","",ihcal_num_bins, ihcal_bin_edges);
	TH1F* h_eT_eta_ihcal_profile_hist = new TH1F("h_eT_eta_ihcal_profile_hist","",ihcal_num_bins, ihcal_bin_edges);
	TProfile* h_eta_ihcal_profile = new TProfile("h_eta_ihcal_profile","",ihcal_num_bins, ihcal_bin_edges);

	int ohcal_num_bins = 24;
    double ohcal_bin_edges[25];
    if (reweighted_bins) {
	    ohcal_bin_edges[0] = ohcal_eta_bin_centers[0] - 0.5 * (ohcal_eta_bin_centers[1] - ohcal_eta_bin_centers[0]);
	    for (int i = 1; i < ohcal_num_bins; ++i) { ohcal_bin_edges[i] = (ohcal_eta_bin_centers[i] + ohcal_eta_bin_centers[i - 1]) / 2.0; }
	    ohcal_bin_edges[ohcal_num_bins] = ohcal_eta_bin_centers[ohcal_num_bins - 1] + 0.5 * (ohcal_eta_bin_centers[ohcal_num_bins - 1] - ohcal_eta_bin_centers[ohcal_num_bins - 2]);
	} else {
		ohcal_bin_edges[0] = eta_bin_centers[0] - 0.5 * (eta_bin_centers[1] - eta_bin_centers[0]);
	    for (int i = 1; i < ohcal_num_bins; ++i) { ohcal_bin_edges[i] = (eta_bin_centers[i] + eta_bin_centers[i - 1]) / 2.0; }
	    ohcal_bin_edges[ohcal_num_bins] = eta_bin_centers[ohcal_num_bins - 1] + 0.5 * (eta_bin_centers[ohcal_num_bins - 1] - eta_bin_centers[ohcal_num_bins - 2]);
	}
	TProfile* h_eT_eta_ohcal_profile = new TProfile("h_eT_eta_ohcal_profile","",ohcal_num_bins, ohcal_bin_edges);
	TH1F* h_eT_eta_ohcal_profile_hist = new TH1F("h_eT_eta_ohcal_profile_hist","",ohcal_num_bins, ohcal_bin_edges);
	TProfile* h_eta_ohcal_profile = new TProfile("h_eta_ohcal_profile","",ohcal_num_bins, ohcal_bin_edges);

	int calo_num_bins = 24;
    double calo_bin_edges[25];
    if (reweighted_bins) {
    	calo_bin_edges[0] = calo_eta_bin_centers[0] - 0.5 * (calo_eta_bin_centers[1] - calo_eta_bin_centers[0]);
	    for (int i = 1; i < calo_num_bins; ++i) { calo_bin_edges[i] = (calo_eta_bin_centers[i] + calo_eta_bin_centers[i - 1]) / 2.0; }
	    calo_bin_edges[calo_num_bins] = calo_eta_bin_centers[calo_num_bins - 1] + 0.5 * (calo_eta_bin_centers[calo_num_bins - 1] - calo_eta_bin_centers[calo_num_bins - 2]);
    } else {
    	calo_bin_edges[0] = eta_bin_centers[0] - 0.5 * (eta_bin_centers[1] - eta_bin_centers[0]);
    	for (int i = 1; i < calo_num_bins; ++i) { calo_bin_edges[i] = (eta_bin_centers[i] + eta_bin_centers[i - 1]) / 2.0; }
    	calo_bin_edges[calo_num_bins] = eta_bin_centers[calo_num_bins - 1] + 0.5 * (eta_bin_centers[calo_num_bins - 1] - eta_bin_centers[calo_num_bins - 2]);
    }
    TProfile* h_eT_eta_calo_profile = new TProfile("h_eT_eta_calo_profile","",calo_num_bins, calo_bin_edges);
	TH1F* h_eT_eta_calo_profile_hist = new TH1F("h_eT_eta_calo_profile_hist","",calo_num_bins, calo_bin_edges);
	TProfile* h_eta_calo_profile = new TProfile("h_eta_calo_profile","",calo_num_bins, calo_bin_edges);

	int hcal_num_bins = 24;
    double hcal_bin_edges[25];
    if (reweighted_bins) {
    	hcal_bin_edges[0] = hcal_eta_bin_centers[0] - 0.5 * (hcal_eta_bin_centers[1] - hcal_eta_bin_centers[0]);
	    for (int i = 1; i < hcal_num_bins; ++i) { hcal_bin_edges[i] = (hcal_eta_bin_centers[i] + hcal_eta_bin_centers[i - 1]) / 2.0; }
	    hcal_bin_edges[hcal_num_bins] = hcal_eta_bin_centers[hcal_num_bins - 1] + 0.5 * (hcal_eta_bin_centers[hcal_num_bins - 1] - hcal_eta_bin_centers[hcal_num_bins - 2]);
    } else {
    	hcal_bin_edges[0] = eta_bin_centers[0] - 0.5 * (eta_bin_centers[1] - eta_bin_centers[0]);
    	for (int i = 1; i < hcal_num_bins; ++i) { hcal_bin_edges[i] = (eta_bin_centers[i] + eta_bin_centers[i - 1]) / 2.0; }
    	hcal_bin_edges[hcal_num_bins] = eta_bin_centers[hcal_num_bins - 1] + 0.5 * (eta_bin_centers[hcal_num_bins - 1] - eta_bin_centers[hcal_num_bins - 2]);
    }
    TProfile* h_eT_eta_hcal_profile = new TProfile("h_eT_eta_hcal_profile","",hcal_num_bins, hcal_bin_edges);
	TH1F* h_eT_eta_hcal_profile_hist = new TH1F("h_eT_eta_hcal_profile_hist","",hcal_num_bins, hcal_bin_edges);
	TProfile* h_eta_hcal_profile = new TProfile("h_eta_hcal_profile","",hcal_num_bins, hcal_bin_edges);
    
    TChain chain("ttree");

    const char* inputDirectory = "/sphenix/tg/tg01/commissioning/CaloCalibWG/egm2153/detdeta_run24auau/";
    if (!strcmp("zs_200_100_100_fixed_build",opt_tag)) {
    	for (int i = 0; i < 3000; i++) {
			TString wildcardPath = TString::Format("%sevents_%d_ZS_200_100_100_pedestal_unc_%d.root", inputDirectory, runnumber, i);
			chain.Add(wildcardPath);
		}
    } else if (!strcmp("zs_60_30_30_fixed_build",opt_tag)) {
    	for (int i = 0; i < 3000; i++) {
			TString wildcardPath = TString::Format("%sevents_%d_ZS_60_30_30_pedestal_unc_%d.root", inputDirectory, runnumber, i);
			chain.Add(wildcardPath);
		}
    } else {
    	for (int i = 0; i < 3000; i++) {
			TString wildcardPath = TString::Format("%sevents_%d_pedestal_unc_%d.root", inputDirectory, runnumber, i);
			chain.Add(wildcardPath);
		}
    }

    int m_simtwrmult_cemc;
    float m_simtwr_cemc_e[cemcSize];
    int m_simtwr_cemc_ieta[cemcSize];
    int m_simtwr_cemc_iphi[cemcSize];
    float m_simtwr_cemc_eta[cemcSize];

    int m_simtwrmult_ihcal;
    float m_simtwr_ihcal_e[ihcalSize];
    int m_simtwr_ihcal_ieta[ihcalSize];
    int m_simtwr_ihcal_iphi[ihcalSize];
    float m_simtwr_ihcal_eta[ihcalSize];
    
    int m_simtwrmult_ohcal;
    float m_simtwr_ohcal_e[ohcalSize];
    int m_simtwr_ohcal_ieta[ohcalSize];
    int m_simtwr_ohcal_iphi[ohcalSize];
    float m_simtwr_ohcal_eta[ohcalSize];

    chain.SetBranchStatus("*", 0);

     // Set branch addresses
    int use_emcal = 1;
    int use_hcal = 1;
    if (use_emcal) {
    	chain.SetBranchStatus("sectorem", 1);
	    chain.SetBranchStatus("emcalen", 1);
	    chain.SetBranchStatus("emcaletabin", 1);
	    chain.SetBranchStatus("emcalphibin", 1);
	    chain.SetBranchStatus("emetacor", 1);
	    chain.SetBranchAddress("sectorem", &m_simtwrmult_cemc);
	    chain.SetBranchAddress("emcalen", m_simtwr_cemc_e);
	    chain.SetBranchAddress("emcaletabin", m_simtwr_cemc_ieta);
	    chain.SetBranchAddress("emcalphibin", m_simtwr_cemc_iphi);
	    chain.SetBranchAddress("emetacor", m_simtwr_cemc_eta);
	}
	if (use_hcal) {
		chain.SetBranchStatus("sectorih", 1);
	    chain.SetBranchStatus("ihcalen", 1);
	    chain.SetBranchStatus("ihcaletabin", 1);
	    chain.SetBranchStatus("ihcalphibin", 1);
	    chain.SetBranchStatus("ihetacor", 1);
	    chain.SetBranchAddress("sectorih", &m_simtwrmult_ihcal);
	    chain.SetBranchAddress("ihcalen", m_simtwr_ihcal_e);
	    chain.SetBranchAddress("ihcaletabin", m_simtwr_ihcal_ieta);
	    chain.SetBranchAddress("ihcalphibin", m_simtwr_ihcal_iphi);
	    chain.SetBranchAddress("ihetacor", m_simtwr_ihcal_eta);

	    chain.SetBranchStatus("sectoroh", 1);
	    chain.SetBranchStatus("ohcalen", 1);
	    chain.SetBranchStatus("ohcaletabin", 1);
	    chain.SetBranchStatus("ohcalphibin", 1);
	    chain.SetBranchStatus("ohetacor", 1);
	    chain.SetBranchAddress("sectoroh", &m_simtwrmult_ohcal);
	    chain.SetBranchAddress("ohcalen", m_simtwr_ohcal_e);
	    chain.SetBranchAddress("ohcaletabin", m_simtwr_ohcal_ieta);
	    chain.SetBranchAddress("ohcalphibin", m_simtwr_ohcal_iphi);
	    chain.SetBranchAddress("ohetacor", m_simtwr_ohcal_eta);
	}

	int eventnumber = 0;
	float totalweights = 0.0;

    Long64_t nEntries = chain.GetEntries();
    std::cout << nEntries << std::endl;
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 1000; ++entry) {
        chain.GetEntry(entry);
    	if (eventnumber % 1000 == 0) cout << "event " << eventnumber << endl;

		float emcale = 0;
		float ihcale = 0;
		float ohcale = 0;
		float totale = 0;
		
		float E_emcalbin[96] = {0};
		float E_ihcal[24] = {0};
  		float E_ohcal[24] = {0};
  		float eta_emcal[96] = {0};
  		float eta_ihcal[24] = {0};
  		float eta_ohcal[24] = {0};

  		float E_emcal[24] = {0};
  		float E_calo[24] = {0};
  		float E_hcal[24] = {0};
  		float vz_weight = 1.0;

  		int emcaleta[24] = {0};
  		int hcaleta[24] = {0};
  		int caloeta[24] = {0};

  		eventnumber++;
  		totalweights += vz_weight;
		if (use_emcal) {
			for (int i = 0; i < m_simtwrmult_cemc; i++) {
				std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i]);
			    auto it = emcal_hot_dead_map.find(hot_tower);
			    if (it != emcal_hot_dead_map.end()) { continue; }
				h_2D_emcal_calib->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight);
				h_2D_emcal_calibT->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
				emcale += m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]); 
				h_emcal->Fill(m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]), vz_weight);
				E_emcalbin[m_simtwr_cemc_ieta[i]] += m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]);
				eta_emcal[m_simtwr_cemc_ieta[i]] = m_simtwr_cemc_eta[i];
			}
		}
		if (use_hcal) {
			for (int i = 0; i < m_simtwrmult_ihcal; i++) {
				std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i]);
			    auto it = ihcal_hot_dead_map.find(hot_tower);
			    if (it != ihcal_hot_dead_map.end()) { continue; }
				double simtwr_ihcal_e = m_simtwr_ihcal_e[i]*1.031443343;
				h_2D_ihcal_calib->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], simtwr_ihcal_e*vz_weight);
				h_2D_ihcal_calibT->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], simtwr_ihcal_e*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
				ihcale += simtwr_ihcal_e/cosh(m_simtwr_ihcal_eta[i]); 
				h_ihcal->Fill(simtwr_ihcal_e/cosh(m_simtwr_ihcal_eta[i]), vz_weight);
				E_ihcal[m_simtwr_ihcal_ieta[i]] += simtwr_ihcal_e/cosh(m_simtwr_ihcal_eta[i]);
				eta_ihcal[m_simtwr_ihcal_ieta[i]] = m_simtwr_ihcal_eta[i];
			}

			for (int i = 0; i < m_simtwrmult_ohcal; i++) {
				std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i]);
			    auto it = ohcal_hot_dead_map.find(hot_tower);
			    if (it != ohcal_hot_dead_map.end()) { continue; }
			    double simtwr_ohcal_e = m_simtwr_ohcal_e[i]*1.029804356;
				h_2D_ohcal_calib->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], simtwr_ohcal_e*vz_weight);
				h_2D_ohcal_calibT->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], simtwr_ohcal_e*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
				ohcale += simtwr_ohcal_e/cosh(m_simtwr_ohcal_eta[i]); 
				h_ohcal->Fill(simtwr_ohcal_e/cosh(m_simtwr_ohcal_eta[i]), vz_weight);
				E_ohcal[m_simtwr_ohcal_ieta[i]] += simtwr_ohcal_e/cosh(m_simtwr_ohcal_eta[i]);
				eta_ohcal[m_simtwr_ohcal_ieta[i]] = m_simtwr_ohcal_eta[i];
			}
		}

		for (int i = 0; i < 96; ++i) {
    		double eta = eta_emcal[i];
		    for (int bin = 0; bin < emcal_num_bins; ++bin) {
		        if (eta >= emcal_bin_edges[bin] && eta < emcal_bin_edges[bin + 1]) {
		            E_emcal[bin] += E_emcalbin[i];
		            emcaleta[bin] += 1;
		            break;
		        }
		    }
		    for (int bin = 0; bin < calo_num_bins; ++bin) {
		        if (eta >= calo_bin_edges[bin] && eta < calo_bin_edges[bin + 1]) {
		            E_calo[bin] += E_emcalbin[i];
		            caloeta[bin] += 1;
		            break;
		        }
		    }
		}
		for (int i = 0; i < 24; i++) {
			double ih_eta = eta_ihcal[i];
			double oh_eta = eta_ohcal[i];
			for (int bin = 0; bin < ihcal_num_bins; ++bin) {
				if (ih_eta >= calo_bin_edges[bin] && ih_eta < calo_bin_edges[bin + 1]) {
					E_calo[bin] += E_ihcal[i];
					E_hcal[bin] += E_ihcal[i];
					caloeta[bin] += 1;
					hcaleta[bin] += 1;
					break;
				}
			}
			for (int bin = 0; bin < ohcal_num_bins; ++bin) {
				if (oh_eta >= calo_bin_edges[bin] && oh_eta < calo_bin_edges[bin + 1]) {
					E_calo[bin] += E_ohcal[i];
					E_hcal[bin] += E_ohcal[i];
					caloeta[bin] += 1;
					hcaleta[bin] += 1;
					break;
				}
			}
		}

		for (int i = 0; i < 24; i++) {
			h_eT_eta_ihcal_profile->Fill(eta_ihcal[i],E_ihcal[i],vz_weight);
			h_eT_eta_ohcal_profile->Fill(eta_ohcal[i],E_ohcal[i],vz_weight);
			h_eta_ihcal_profile->Fill(eta_ihcal[i],1);
			h_eta_ohcal_profile->Fill(eta_ohcal[i],1);
			if (E_emcal[i] != 0) {
				h_eT_eta_emcal_profile->Fill(eta_bin_centers[i],E_emcal[i],vz_weight);
				h_eta_emcal_profile->Fill(eta_bin_centers[i],emcaleta[i]);
			}
			if (E_calo[i] != 0) {
				h_eT_eta_calo_profile->Fill(eta_bin_centers[i],E_calo[i],vz_weight);
				h_eta_calo_profile->Fill(eta_bin_centers[i],caloeta[i]);
			}
			if (E_hcal[i] != 0) {
				h_eT_eta_hcal_profile->Fill(eta_bin_centers[i],E_hcal[i],vz_weight);
				h_eta_hcal_profile->Fill(eta_bin_centers[i],hcaleta[i]);
			}

		}
		
		totale = ihcale + ohcale;
		if (use_emcal && use_hcal) h_event_energy->Fill(totale + emcale, vz_weight);
		if (use_hcal) h_event_hcal_energy->Fill(totale, vz_weight);
		if (use_emcal) h_event_emcal_energy->Fill(emcale, vz_weight);
		if (use_hcal) h_event_ihcal_energy->Fill(ihcale, vz_weight);
		if (use_hcal) h_event_ohcal_energy->Fill(ohcale, vz_weight);
		
	}

	if (use_emcal) {
		h_2D_emcal_calib->Scale(1.0/totalweights);
		h_2D_emcal_calibT->Scale(1.0/totalweights);
		h_emcal->Scale(1.0/totalweights);
		for (int i = 1; i <= emcal_num_bins; ++i) {
	        double bin_width = emcal_bin_edges[i] - emcal_bin_edges[i - 1];
	        h_eT_eta_emcal_profile_hist->SetBinContent(i, h_eT_eta_emcal_profile->GetBinContent(i) / bin_width);
	        h_eT_eta_emcal_profile_hist->SetBinError(i, h_eT_eta_emcal_profile->GetBinError(i) / bin_width);
	    }
	}
	if (use_hcal) {
		h_2D_ihcal_calib->Scale(1.0/totalweights);
		h_2D_ihcal_calibT->Scale(1.0/totalweights);
		h_ihcal->Scale(1.0/totalweights);
		for (int i = 1; i <= ihcal_num_bins; ++i) {
	        double bin_width = ihcal_bin_edges[i] - ihcal_bin_edges[i - 1];
	        h_eT_eta_ihcal_profile_hist->SetBinContent(i, h_eT_eta_ihcal_profile->GetBinContent(i) / bin_width);
	        h_eT_eta_ihcal_profile_hist->SetBinError(i, h_eT_eta_ihcal_profile->GetBinError(i) / bin_width);
	    }

		h_2D_ohcal_calib->Scale(1.0/totalweights);
		h_2D_ohcal_calibT->Scale(1.0/totalweights);
		h_ohcal->Scale(1.0/totalweights);
		for (int i = 1; i <= ohcal_num_bins; ++i) {
	        double bin_width = ohcal_bin_edges[i] - ohcal_bin_edges[i - 1];
	        h_eT_eta_ohcal_profile_hist->SetBinContent(i, h_eT_eta_ohcal_profile->GetBinContent(i) / bin_width);
	        h_eT_eta_ohcal_profile_hist->SetBinError(i, h_eT_eta_ohcal_profile->GetBinError(i) / bin_width);
	    }

		for (int i = 1; i <= hcal_num_bins; ++i) {
	        double bin_width = hcal_bin_edges[i] - hcal_bin_edges[i - 1];
	        h_eT_eta_hcal_profile_hist->SetBinContent(i, h_eT_eta_hcal_profile->GetBinContent(i) / bin_width);
	        h_eT_eta_hcal_profile_hist->SetBinError(i, h_eT_eta_hcal_profile->GetBinError(i) / bin_width);
	    }
	}

	if (use_emcal && use_hcal) {
		for (int i = 1; i <= calo_num_bins; ++i) {
	        double bin_width = calo_bin_edges[i] - calo_bin_edges[i - 1];
	        h_eT_eta_calo_profile_hist->SetBinContent(i, h_eT_eta_calo_profile->GetBinContent(i) / bin_width);
	        h_eT_eta_calo_profile_hist->SetBinError(i, h_eT_eta_calo_profile->GetBinError(i) / bin_width);
	    }
	}

	if (use_emcal && use_hcal) h_event_energy->Scale(1.0/totalweights);
	if (use_hcal) h_event_hcal_energy->Scale(1.0/totalweights);
	if (use_emcal) h_event_emcal_energy->Scale(1.0/totalweights);
	if (use_hcal) h_event_ihcal_energy->Scale(1.0/totalweights);
	if (use_hcal) h_event_ohcal_energy->Scale(1.0/totalweights);
	
	out->Write();
	out->Close();


}
