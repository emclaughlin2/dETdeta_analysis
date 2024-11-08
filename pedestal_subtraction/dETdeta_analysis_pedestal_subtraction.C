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

std::set<std::tuple<int, int>> emcal_hot_dead_map = {{13,232},{47,138},{48,7}};
std::set<std::tuple<int, int>> ihcal_hot_dead_map = {{8,32},{7,51}};
std::set<std::tuple<int, int>> ohcal_hot_dead_map;

double emcal_eta_bin_centers[24];
double ihcal_eta_bin_centers[24];
double ohcal_eta_bin_centers[24];

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

	TFile *hotdeadfile = new TFile("/sphenix/user/egm2153/calib_study/detdeta/analysis/runs23727_23746/run23727_hotdeadmap_z_-20_20_new_status.root", "READ"); 
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
	}
	for (int i = 0; i < ohcal_eta_bins->size(); i++) {
		ohcal_eta_bin_centers[i] = (*ohcal_eta_bins)[i];
	}

	hotdeadfile->Close();

}

void fill_zs_cross_calib() {

	float emcal_zs_calib[96][256];
	float ihcal_zs_calib[24][64];
	float ohcal_zs_calib[24][64];
	TFile *zscalibfile = new TFile("/sphenix/user/egm2153/calib_study/zs_testing/zs_testing_run23727_p008_new.root", "READ");
	TTree *zscalibtree = dynamic_cast<TTree*>(zscalibfile->Get("zs_calib_tree"));
	zscalibtree->SetBranchAddress("emcal_zs_calib", emcal_zs_calib);
	zscalibtree->SetBranchAddress("ihcal_zs_calib", ihcal_zs_calib);
	zscalibtree->SetBranchAddress("ohcal_zs_calib", ohcal_zs_calib);
	zscalibtree->GetEntry(0);
	for (int i = 0; i < 96; i++) {
		for (int j = 0; j < 256; j++) {
			em_zs_calib[i][j] = emcal_zs_calib[i][j];
			if (i < 24 && j < 64) {
				ih_zs_calib[i][j] = ihcal_zs_calib[i][j];
				oh_zs_calib[i][j] = ohcal_zs_calib[i][j];
			}
		}
	}
	zscalibfile->Close();
}

void dETdeta_analysis_pedestal_subtraction(int runnumber = 23727, int reweighted_bins = 1, int zs = 0, int zs_value = 10, int cross_calib = 1, int time = 0, const char* opt_tag = "") {

	string filename = "dETdeta_analysis_pedestal_subtraction";
	filename += "_" + to_string(runnumber);
	string opttag = opt_tag;
	if (strcmp(opt_tag,"")) { filename += "_" + opttag; }
	string zstag; 
	if (zs == 0) { zstag = "nozs"; }
	if (zs == 1) { zstag = "zs_" + to_string(zs_value) + "ADC"; }
	string crosscalibtag = (cross_calib?"":"no_cross_calib_");
	string timetag = (time?"emcal_timecut_":"");
  	string weighttag = (reweighted_bins?"reweighted_bins":"default_bins");
  	filename += "_" + zstag + "_" + crosscalibtag + timetag + weighttag + ".root";

	fill_hot_dead_map_eta_bin_centers();
	fill_zs_cross_calib();

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

	TH2F* h_em_zero_zscross = new TH2F("h_em_zero_zscross","",96,0,96,256,0,256);
	TH2F* h_ih_zero_zscross = new TH2F("h_ih_zero_zscross","",24,0,24,64,0,64);
	TH2F* h_oh_zero_zscross = new TH2F("h_oh_zero_zscross","",24,0,24,64,0,64);

	TH1I* h_em_zs_rate = new TH1I("h_em_zs_rate","",2,0,2);
	TH1I* h_ih_zs_rate = new TH1I("h_ih_zs_rate","",2,0,2);
	TH1I* h_oh_zs_rate = new TH1I("h_oh_zs_rate","",2,0,2);
	
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
    TH1F* h_eT_emcal = new TH1F("h_eT_emcal","",emcal_num_bins,emcal_bin_edges);

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
    TH1F* h_eT_ihcal = new TH1F("h_eT_ihcal","",ihcal_num_bins,ihcal_bin_edges);

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
	TH1F* h_eT_ohcal = new TH1F("h_eT_ohcal","",ohcal_num_bins,ohcal_bin_edges);
    
    TChain chain("ttree");

    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
	for (int i = 0; i < 74; i++) {
		TString wildcardPath = TString::Format("%sevents_ped_sub_%d__unc_%d.root", inputDirectory, runnumber, i);
		chain.Add(wildcardPath);
	}

    int m_simtwrmult_cemc;
    float m_simtwr_cemc_e[cemcSize];
    float m_simtwr_cemc_zs_e[cemcSize];
    int m_simtwr_cemc_ieta[cemcSize];
    int m_simtwr_cemc_iphi[cemcSize];
    float m_simtwr_cemc_eta[cemcSize];
    int m_simtwr_cemc_adc[cemcSize];
    int m_simtwr_cemc_zs_adc[cemcSize];
    float m_simtwr_cemc_time[cemcSize];
    int m_simtwrmult_ihcal;
    float m_simtwr_ihcal_e[ihcalSize];
    float m_simtwr_ihcal_zs_e[ihcalSize];
    int m_simtwr_ihcal_ieta[ihcalSize];
    int m_simtwr_ihcal_iphi[ihcalSize];
    float m_simtwr_ihcal_eta[ihcalSize];
    int m_simtwr_ihcal_adc[ihcalSize];
    int m_simtwr_ihcal_zs_adc[ihcalSize];
    float m_simtwr_ihcal_time[ihcalSize];
    int m_simtwrmult_ohcal;
    float m_simtwr_ohcal_e[ohcalSize];
    float m_simtwr_ohcal_zs_e[ohcalSize];
    int m_simtwr_ohcal_ieta[ohcalSize];
    int m_simtwr_ohcal_iphi[ohcalSize];
    float m_simtwr_ohcal_eta[ohcalSize];
    int m_simtwr_ohcal_adc[ohcalSize];
    int m_simtwr_ohcal_zs_adc[ohcalSize];
    float m_simtwr_ohcal_time[ohcalSize];

     // Set branch addresses
    int use_emcal = 1;
    int use_hcal = 0;
    if (use_emcal) {
	    chain.SetBranchAddress("sectorem", &m_simtwrmult_cemc);
	    chain.SetBranchAddress("emcalen", m_simtwr_cemc_e);
	    chain.SetBranchAddress("emcaletabin", m_simtwr_cemc_ieta);
	    chain.SetBranchAddress("emcalphibin", m_simtwr_cemc_iphi);
	    chain.SetBranchAddress("emetacor", m_simtwr_cemc_eta);
	    chain.SetBranchAddress("emcalzs", m_simtwr_cemc_zs_e);
		chain.SetBranchAddress("emcaladc", m_simtwr_cemc_adc);
		chain.SetBranchAddress("emcalzsadc", m_simtwr_cemc_zs_adc);
		chain.SetBranchAddress("emcalt", m_simtwr_cemc_time);
	}
	if (use_hcal) {
	    chain.SetBranchAddress("sectorih", &m_simtwrmult_ihcal);
	    chain.SetBranchAddress("ihcalen", m_simtwr_ihcal_e);
	    chain.SetBranchAddress("ihcaletabin", m_simtwr_ihcal_ieta);
	    chain.SetBranchAddress("ihcalphibin", m_simtwr_ihcal_iphi);
	    chain.SetBranchAddress("ihetacor", m_simtwr_ihcal_eta);
	    chain.SetBranchAddress("ihcalzs", m_simtwr_ihcal_zs_e);
	    chain.SetBranchAddress("ihcaladc", m_simtwr_ihcal_adc);
	    chain.SetBranchAddress("ihcalzsadc", m_simtwr_ihcal_zs_adc);
	    chain.SetBranchAddress("ihcalt", m_simtwr_ihcal_time);

	    chain.SetBranchAddress("sectoroh", &m_simtwrmult_ohcal);
	    chain.SetBranchAddress("ohcalen", m_simtwr_ohcal_e);
	    chain.SetBranchAddress("ohcaletabin", m_simtwr_ohcal_ieta);
	    chain.SetBranchAddress("ohcalphibin", m_simtwr_ohcal_iphi);
	    chain.SetBranchAddress("ohetacor", m_simtwr_ohcal_eta);
		chain.SetBranchAddress("ohcalzs", m_simtwr_ohcal_zs_e);
	    chain.SetBranchAddress("ohcaladc", m_simtwr_ohcal_adc);
	    chain.SetBranchAddress("ohcalzsadc", m_simtwr_ohcal_zs_adc);
	    chain.SetBranchAddress("ohcalt", m_simtwr_ohcal_time);
	}

	int eventnumber = 0;
	float totalweights = 0.0;

    Long64_t nEntries = chain.GetEntries();
    std::cout << nEntries << std::endl;
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 20000; ++entry) {
        chain.GetEntry(entry);
    	if (eventnumber % 1000 == 0) cout << "event " << eventnumber << endl;

		float emcale = 0;
		float ihcale = 0;
		float ohcale = 0;
		float totale = 0;
		float truthe = 0;

		float E_emcal[96] = {0};
		float E_ihcal[24] = {0};
  		float E_ohcal[24] = {0};

  		float vz_weight = 1.0;

  		eventnumber++;
  		totalweights += vz_weight;
		if (use_emcal) {
			for (int i = 0; i < m_simtwrmult_cemc; i++) {
				if (m_simtwr_cemc_ieta[i] < 8) { continue; }
				int ieta = m_simtwr_cemc_ieta[i]; // needed to manually change mapping due to ADC board switching just previous to pedestal data taking
				if (m_simtwr_cemc_ieta[i] >= 56 && m_simtwr_cemc_ieta[i] < 64 && m_simtwr_cemc_iphi[i] >= 32 && m_simtwr_cemc_iphi[i] < 40) {
					ieta += 8;
				}
				if (m_simtwr_cemc_ieta[i] >= 32 && m_simtwr_cemc_ieta[i] < 40 && m_simtwr_cemc_iphi[i] >= 144 && m_simtwr_cemc_iphi[i] < 152) {
					ieta -= 8;
				}
				std::tuple<int, int> hot_tower = std::make_tuple(ieta, m_simtwr_cemc_iphi[i]);
			    auto it = emcal_hot_dead_map.find(hot_tower);
			    if (it != emcal_hot_dead_map.end()) { continue; }
				if (zs == 1 && m_simtwr_cemc_zs_adc[i] < zs_value) {
					float zs_e = 0; // = m_simtwr_cemc_zs_e[i]; 
					if (em_zs_calib[ieta][m_simtwr_cemc_iphi[i]] == 0) {
						h_em_zero_zscross->Fill(ieta, m_simtwr_cemc_iphi[i]);
					} else {
						if (cross_calib) zs_e = m_simtwr_cemc_zs_e[i]/em_zs_calib[ieta][m_simtwr_cemc_iphi[i]];
						else zs_e = m_simtwr_cemc_zs_e[i];
					}
					h_2D_emcal_calib->Fill(ieta, m_simtwr_cemc_iphi[i], zs_e*vz_weight);
					h_2D_emcal_calibT->Fill(ieta, m_simtwr_cemc_iphi[i], zs_e*vz_weight/cosh(m_simtwr_cemc_eta[i]));
					emcale += zs_e/cosh(m_simtwr_cemc_eta[i]); 
					h_emcal->Fill(zs_e/cosh(m_simtwr_cemc_eta[i]), vz_weight);
					h_eT_emcal->Fill(m_simtwr_cemc_eta[i],zs_e*vz_weight/cosh(m_simtwr_cemc_eta[i]));
					h_em_zs_rate->Fill(0);
				} else {
					if (time && (m_simtwr_cemc_time[i] < -2 || m_simtwr_cemc_time[i] > 2)) { continue; }
					h_2D_emcal_calib->Fill(ieta, m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight);
					h_2D_emcal_calibT->Fill(ieta, m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
					emcale += m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]); 
					h_emcal->Fill(m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]), vz_weight);
					h_eT_emcal->Fill(m_simtwr_cemc_eta[i],m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
					h_em_zs_rate->Fill(1);
				}
			}
		}
		if (use_hcal) {
			for (int i = 0; i < m_simtwrmult_ihcal; i++) {
				std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i]);
			    auto it = ihcal_hot_dead_map.find(hot_tower);
			    if (it != ihcal_hot_dead_map.end()) { continue; }
			    if (zs == 1 && m_simtwr_ihcal_zs_adc[i] < zs_value) { 
			    	float zs_e = 0; // m_simtwr_ihcal_zs_e[i]; 
					if (ih_zs_calib[m_simtwr_ihcal_ieta[i]][m_simtwr_ihcal_iphi[i]] == 0) {
						h_ih_zero_zscross->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i]);
					} else {
						if (cross_calib) zs_e = m_simtwr_ihcal_zs_e[i]/ih_zs_calib[m_simtwr_ihcal_ieta[i]][m_simtwr_ihcal_iphi[i]];
						else zs_e = m_simtwr_ihcal_zs_e[i];
					}
					h_2D_ihcal_calib->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], zs_e*vz_weight);
					h_2D_ihcal_calibT->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], zs_e*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
					ihcale += zs_e/cosh(m_simtwr_ihcal_eta[i]); 
					h_ihcal->Fill(zs_e/cosh(m_simtwr_ihcal_eta[i]), vz_weight);
					h_eT_ihcal->Fill(m_simtwr_ihcal_eta[i],zs_e*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
					h_ih_zs_rate->Fill(0);
			    } else {
					h_2D_ihcal_calib->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]*vz_weight);
					h_2D_ihcal_calibT->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
					ihcale += m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]); 
					h_ihcal->Fill(m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]), vz_weight);
					h_eT_ihcal->Fill(m_simtwr_ihcal_eta[i],m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
					h_ih_zs_rate->Fill(1);
				}
			}

			for (int i = 0; i < m_simtwrmult_ohcal; i++) {
				std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i]);
			    auto it = ohcal_hot_dead_map.find(hot_tower);
			    if (it != ohcal_hot_dead_map.end()) { continue; }
			    if (zs == 1 && m_simtwr_ohcal_zs_adc[i] < zs_value) {
			    	float zs_e = 0; // m_simtwr_ihcal_zs_e[i]; 
					if (oh_zs_calib[m_simtwr_ohcal_ieta[i]][m_simtwr_ohcal_iphi[i]] == 0) {
						h_oh_zero_zscross->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i]);
					} else {
						if (cross_calib) zs_e = m_simtwr_ohcal_zs_e[i]/oh_zs_calib[m_simtwr_ohcal_ieta[i]][m_simtwr_ohcal_iphi[i]];
						else zs_e = m_simtwr_ohcal_zs_e[i];
					}
			    	h_2D_ohcal_calib->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], zs_e*vz_weight);
					h_2D_ohcal_calibT->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], zs_e*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
					ohcale += zs_e/cosh(m_simtwr_ohcal_eta[i]); 
					h_ohcal->Fill(zs_e/cosh(m_simtwr_ohcal_eta[i]), vz_weight);
					h_eT_ohcal->Fill(m_simtwr_ohcal_eta[i],zs_e*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
					h_oh_zs_rate->Fill(0);
			    } else {
					h_2D_ohcal_calib->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]*vz_weight);
					h_2D_ohcal_calibT->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
					ohcale += m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]); 
					h_ohcal->Fill(m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]), vz_weight);
					h_eT_ohcal->Fill(m_simtwr_ohcal_eta[i],m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
					h_oh_zs_rate->Fill(1);
				}
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
		h_eT_emcal->Scale(1.0/totalweights);
		//h_em_zs_rate->Scale(1.0/h_em_zs_rate->GetEntries());
		for (int i = 1; i <= emcal_num_bins; ++i) {
	        double bin_width = emcal_bin_edges[i] - emcal_bin_edges[i - 1];
	        h_eT_emcal->SetBinContent(i, h_eT_emcal->GetBinContent(i) / bin_width);
	        h_eT_emcal->SetBinError(i, h_eT_emcal->GetBinError(i) / bin_width);
	    }
	}
	if (use_hcal) {
		h_2D_ihcal_calib->Scale(1.0/totalweights);
		h_2D_ihcal_calibT->Scale(1.0/totalweights);
		h_ihcal->Scale(1.0/totalweights);
		h_eT_ihcal->Scale(1.0/totalweights);
		//h_ih_zs_rate->Scale(1.0/h_ih_zs_rate->GetEntries());
		for (int i = 1; i <= ihcal_num_bins; ++i) {
	        double bin_width = ihcal_bin_edges[i] - ihcal_bin_edges[i - 1];
	        h_eT_ihcal->SetBinContent(i, h_eT_ihcal->GetBinContent(i) / bin_width);
	        h_eT_ihcal->SetBinError(i, h_eT_ihcal->GetBinError(i) / bin_width);
	    }

		h_2D_ohcal_calib->Scale(1.0/totalweights);
		h_2D_ohcal_calibT->Scale(1.0/totalweights);
		h_ohcal->Scale(1.0/totalweights);
		h_eT_ohcal->Scale(1.0/totalweights);
		//h_oh_zs_rate->Scale(1.0/h_oh_zs_rate->GetEntries());
		for (int i = 1; i <= ohcal_num_bins; ++i) {
	        double bin_width = ohcal_bin_edges[i] - ohcal_bin_edges[i - 1];
	        h_eT_ohcal->SetBinContent(i, h_eT_ohcal->GetBinContent(i) / bin_width);
	        h_eT_ohcal->SetBinError(i, h_eT_ohcal->GetBinError(i) / bin_width);
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
