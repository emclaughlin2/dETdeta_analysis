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

void fill_hot_dead_map_eta_bin_centers(int runnumber, float minus_z, float plus_z) {

	// hot dead maps 
    vector<int> *emcal_hot_dead_ieta = nullptr;
    vector<int> *emcal_hot_dead_iphi = nullptr;
    vector<int> *ihcal_hot_dead_ieta = nullptr;
    vector<int> *ihcal_hot_dead_iphi = nullptr;
    vector<int> *ohcal_hot_dead_ieta = nullptr;
    vector<int> *ohcal_hot_dead_iphi = nullptr;
    vector<float> *ihcal_eta_bins = nullptr;
    vector<float> *ohcal_eta_bins = nullptr;

	TFile *hotdeadfile = new TFile(TString::Format("run%d_hotdeadmap_z_%d_%d.root", runnumber, int(floor(minus_z)), int(floor(plus_z))), "READ");
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

	//std::cout << "emcal hot dead map input size " << emcal_hot_dead_ieta->size() << std::endl;
	for (int i = 0; i < emcal_hot_dead_ieta->size(); i++) {
		tuple<int, int> new_hot_tower = make_tuple((*emcal_hot_dead_ieta)[i], (*emcal_hot_dead_iphi)[i]); 
        emcal_hot_dead_map.insert(new_hot_tower);
	}
	//std::cout << "emcal hot dead map set size " << emcal_hot_dead_map.size() << std::endl;
	
	//std::cout << "ihcal hot dead map input size " << ihcal_hot_dead_ieta->size() << std::endl;
	for (int i = 0; i < ihcal_hot_dead_ieta->size(); i++) {
		tuple<int, int> new_hot_tower = make_tuple((*ihcal_hot_dead_ieta)[i], (*ihcal_hot_dead_iphi)[i]); 
        ihcal_hot_dead_map.insert(new_hot_tower);
	}
	//std::cout << "ihcal hot dead map set size " << ihcal_hot_dead_map.size() << std::endl;

	//std::cout << "ohcal hot dead map input size " << ohcal_hot_dead_ieta->size() << std::endl;
	for (int i = 0; i < ohcal_hot_dead_ieta->size(); i++) {
		tuple<int, int> new_hot_tower = make_tuple((*ohcal_hot_dead_ieta)[i], (*ohcal_hot_dead_iphi)[i]); 
        ohcal_hot_dead_map.insert(new_hot_tower);
	}
	//std::cout << "ohcal hot dead map set size " << ohcal_hot_dead_map.size() << std::endl;

	for (int i = 0; i < ihcal_eta_bins->size(); i++) {
		ihcal_eta_bin_centers[i] = (*ihcal_eta_bins)[i];
		emcal_eta_bin_centers[i] = (*ihcal_eta_bins)[i];
	}
	for (int i = 0; i < ohcal_eta_bins->size(); i++) {
		ohcal_eta_bin_centers[i] = (*ohcal_eta_bins)[i];
	}

	hotdeadfile->Close();

}

void fill_zvertex_centrality(int dataormc, int runnumber, const char* generator) {
	
	TFile *zvertexfile;
	float vz_reweight[200]; 
	float mc_cent[20];
	float data_cent[20];
	if (dataormc) {
		zvertexfile = new TFile(TString::Format("dETdeta_vertex_reweight_run%d_%s.root", runnumber, generator), "READ");
		TTree *vertextree = dynamic_cast<TTree*>(zvertexfile->Get("T"));
		vertextree->SetBranchAddress("vertex_reweight", vz_reweight);
		vertextree->SetBranchAddress("mc_centrality", mc_cent);
		vertextree->GetEntry(0);
		for (int i = 0; i < 200; i++) {
			vertex_reweight[i] = vz_reweight[i];
		}
		centrality_bin.assign(mc_cent, mc_cent+20);
	} else {
		zvertexfile = new TFile(TString::Format("dETdeta_vertex_reweight_run%d_epos.root",runnumber), "READ");
		TTree *vertextree = dynamic_cast<TTree*>(zvertexfile->Get("T"));
		vertextree->SetBranchAddress("data_centrality", data_cent);
		vertextree->GetEntry(0);
		centrality_bin.assign(data_cent, data_cent+20);
	}
	zvertexfile->Close();

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

void dETdeta_analysis(int runnumber = 23727, const char* generator = "", float minus_z = -2, float plus_z = 2, int dataormc = 0, int reweighting = 1, int central = 0, int zs = 0, int zs_value = 10, int time = 0, const char* opt_tag = "") {

	string filename = "dETdeta_analysis_";
	filename += to_string(runnumber) + "_z=0";
	string opttag = opt_tag;
	if (strcmp(opt_tag,"")) { filename += "_" + opttag; }
	string zstag; 
	if (zs == 0) { zstag = "nozs"; }
	if (zs == 1) { zstag = "zs_" + to_string(zs_value) + "ADC"; }
	if (zs == 2) { zstag = "zs_abs" + to_string(zs_value) + "ADC"; }
	string timetag = (time?"emcal_timecut_":"");
  	string dattag = (dataormc?"mc":"data");
  	string weighttag = (reweighting?"reweight":"noweight");
  	string centtag = (central?"0-10":"0-90"); 
  	string gentag = generator;
  	if (!strcmp(generator,"")) filename += "_" + zstag + "_" + timetag + dattag + "_" + weighttag + "_" + centtag + ".root";
  	else filename += "_" + zstag + "_" + timetag + dattag + "_" + weighttag + "_" + centtag + "_" + gentag + ".root";

	fill_hot_dead_map_eta_bin_centers(runnumber, minus_z, plus_z);
	fill_zvertex_centrality(dataormc, runnumber, generator);
	fill_zs_cross_calib();
	assert(centrality_bin.size() == 20);

	std::cout << "centrality_bins" << std::endl;
	for (int i = 0; i < 20; i++) {
		std::cout << centrality_bin[i] << " ";
	}
	std::cout << std::endl;

	TFile *out = new TFile(filename.c_str(),"RECREATE");
	TH1F* h_vz = new TH1F("h_vz","",400, -100, 100);
	TH1F* h_vz_reweight = new TH1F("h_vz_reweight","",400, -100, 100);
	TH1F* h_mbd = new TH1F("h_mbd","",250000,0,250000);

	TH1F* h_event_truth_energy = new TH1F("h_event_truth_energy","",10000,0,10000);
	TH1F* hetdeta = new TH1F("hetdeta","",120,-6,6);
	TH1F* hetdeta_zoom = new TH1F("hetdeta_zoom","",220,-1.1,1.1);
	
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

	TH1F* h_emcal = new TH1F("h_emcal","",1000,0,10);
	TH1F* h_ihcal = new TH1F("h_ihcal","",1000,0,10);
	TH1F* h_ohcal = new TH1F("h_ohcal","",1000,0,10);

	TH2F* h_em_zero_zscross = new TH2F("h_em_zero_zscross","",96,0,96,256,0,256);
	TH2F* h_ih_zero_zscross = new TH2F("h_ih_zero_zscross","",24,0,24,64,0,64);
	TH2F* h_oh_zero_zscross = new TH2F("h_oh_zero_zscross","",24,0,24,64,0,64);
	
	int emcal_num_bins = 24;
	double emcal_bin_edges[25];
    if (!dataormc || reweighting) {
    	emcal_bin_edges[0] = emcal_eta_bin_centers[0] - 0.5 * (emcal_eta_bin_centers[1] - emcal_eta_bin_centers[0]);
	    for (int i = 1; i < emcal_num_bins; ++i) { emcal_bin_edges[i] = (emcal_eta_bin_centers[i] + emcal_eta_bin_centers[i - 1]) / 2.0; }
	    emcal_bin_edges[emcal_num_bins] = emcal_eta_bin_centers[emcal_num_bins - 1] + 0.5 * (emcal_eta_bin_centers[emcal_num_bins - 1] - emcal_eta_bin_centers[emcal_num_bins - 2]);
    } else {
    	emcal_bin_edges[0] = eta_bin_centers[0] - 0.5 * (eta_bin_centers[1] - eta_bin_centers[0]);
    	for (int i = 1; i < emcal_num_bins; ++i) { emcal_bin_edges[i] = (eta_bin_centers[i] + eta_bin_centers[i - 1]) / 2.0; }
    	emcal_bin_edges[emcal_num_bins] = eta_bin_centers[emcal_num_bins - 1] + 0.5 * (eta_bin_centers[emcal_num_bins - 1] - eta_bin_centers[emcal_num_bins - 2]);
    }
    TH1F* h_eT_emcal = new TH1F("h_eT_emcal","",emcal_num_bins,emcal_bin_edges);
	TH1F* hetdeta_emcalbin = new TH1F("hetdeta_emcalbin","",emcal_num_bins, emcal_bin_edges);

	int ihcal_num_bins = 24;
    double ihcal_bin_edges[25];
    if (!dataormc || reweighting) {
    	ihcal_bin_edges[0] = ihcal_eta_bin_centers[0] - 0.5 * (ihcal_eta_bin_centers[1] - ihcal_eta_bin_centers[0]);
	    for (int i = 1; i < ihcal_num_bins; ++i) { ihcal_bin_edges[i] = (ihcal_eta_bin_centers[i] + ihcal_eta_bin_centers[i - 1]) / 2.0; }
	    ihcal_bin_edges[ihcal_num_bins] = ihcal_eta_bin_centers[ihcal_num_bins - 1] + 0.5 * (ihcal_eta_bin_centers[ihcal_num_bins - 1] - ihcal_eta_bin_centers[ihcal_num_bins - 2]);
    } else {
    	ihcal_bin_edges[0] = eta_bin_centers[0] - 0.5 * (eta_bin_centers[1] - eta_bin_centers[0]);
    	for (int i = 1; i < ihcal_num_bins; ++i) { ihcal_bin_edges[i] = (eta_bin_centers[i] + eta_bin_centers[i - 1]) / 2.0; }
    	ihcal_bin_edges[ihcal_num_bins] = eta_bin_centers[ihcal_num_bins - 1] + 0.5 * (eta_bin_centers[ihcal_num_bins - 1] - eta_bin_centers[ihcal_num_bins - 2]);
    }
    TH1F* h_eT_ihcal = new TH1F("h_eT_ihcal","",ihcal_num_bins,ihcal_bin_edges);
	TH1F* hetdeta_ihcalbin = new TH1F("hetdeta_ihcalbin","",ihcal_num_bins, ihcal_bin_edges);

	int ohcal_num_bins = 24;
    double ohcal_bin_edges[25];
    if (!dataormc || reweighting) {
	    ohcal_bin_edges[0] = ohcal_eta_bin_centers[0] - 0.5 * (ohcal_eta_bin_centers[1] - ohcal_eta_bin_centers[0]);
	    for (int i = 1; i < ohcal_num_bins; ++i) { ohcal_bin_edges[i] = (ohcal_eta_bin_centers[i] + ohcal_eta_bin_centers[i - 1]) / 2.0; }
	    ohcal_bin_edges[ohcal_num_bins] = ohcal_eta_bin_centers[ohcal_num_bins - 1] + 0.5 * (ohcal_eta_bin_centers[ohcal_num_bins - 1] - ohcal_eta_bin_centers[ohcal_num_bins - 2]);
	} else {
		ohcal_bin_edges[0] = eta_bin_centers[0] - 0.5 * (eta_bin_centers[1] - eta_bin_centers[0]);
	    for (int i = 1; i < ohcal_num_bins; ++i) { ohcal_bin_edges[i] = (eta_bin_centers[i] + eta_bin_centers[i - 1]) / 2.0; }
	    ohcal_bin_edges[ohcal_num_bins] = eta_bin_centers[ohcal_num_bins - 1] + 0.5 * (eta_bin_centers[ohcal_num_bins - 1] - eta_bin_centers[ohcal_num_bins - 2]);
	}
	TH1F* h_eT_ohcal = new TH1F("h_eT_ohcal","",ohcal_num_bins,ohcal_bin_edges);
	TH1F* hetdeta_ohcalbin = new TH1F("hetdeta_ohcalbin","",ohcal_num_bins, ohcal_bin_edges);
    
    TChain chain("ttree");

    if (dataormc && !strcmp(generator, "epos")) {
    	// location of EPOS files 
	    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/eposrun/plugdoor_condor/";
	    for (int i = 0; i < 1000; i++) {
	    	TString wildcardPath = TString::Format("%sOutDir%d/epos_sim_output.root", inputDirectory, i);
	    	chain.Add(wildcardPath);
	    }
    } else if (dataormc && !strcmp(generator, "hijing")) {
    	// location of HIJING files 
	    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
	    for (int i = 0; i < 2500; i++) {
	    	TString wildcardPath = TString::Format("%sevents_20240110_run101_nopileup_mc_cor_%d.root", inputDirectory, i);
	    	chain.Add(wildcardPath);
	    }
    } else if (dataormc && !strcmp(generator, "ampt")) {
    	// location of AMPT files 
	    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/amptrun/plugdoor_condor/";
	    for (int i = 0; i < 1000; i++) {
	    	TString wildcardPath = TString::Format("%sOutDir%d/ampt_sim_output.root", inputDirectory , i);
	    	chain.Add(wildcardPath);
	    }
    } else if (!dataormc) {
    	const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
		//TString wildcardPath = TString::Format("%sevents_p008_zs_%d_data_cor*.root", inputDirectory, runnumber); 
		//TString wildcardPath = TString::Format("%sevents_p008_zs_%d_z-2cm_data_cor_*.root", inputDirectory, runnumber);
		for (int i = 0; i < 230; i++) {
			TString wildcardPath = TString::Format("%sevents_p008_zs_%d_data_cor_%d.root", inputDirectory, runnumber, i);
			chain.Add(wildcardPath);
		}
		//TString wildcardPath = TString::Format("%sevents_20240112_p007_23696_data_cor*.root", inputDirectory); 
		//TString wildcardPath = TString::Format("%sevents_pA_%d_data_cor_0.root", inputDirectory, runnumber); 
    	
    } else {
    	std::cout << "generator/data not found" << std::endl;
    	return;
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
    int m_sectormb;
    float m_mbenergy[mbdSize];
    int m_g4;
    float m_g4_e[g4Size];
    float m_g4_eta[g4Size];
    float m_g4_pt[g4Size];
    float m_g4_pz[g4Size];
    float m_vtx[vtxSize];

     // Set branch addresses
    chain.SetBranchAddress("sectorem", &m_simtwrmult_cemc);
    chain.SetBranchAddress("emcalen", m_simtwr_cemc_e);
    chain.SetBranchAddress("emcaletabin", m_simtwr_cemc_ieta);
    chain.SetBranchAddress("emcalphibin", m_simtwr_cemc_iphi);
    chain.SetBranchAddress("emetacor", m_simtwr_cemc_eta);
    if (!dataormc) {
    	chain.SetBranchAddress("emcalzs", m_simtwr_cemc_zs_e);
	    chain.SetBranchAddress("emcaladc", m_simtwr_cemc_adc);
	    chain.SetBranchAddress("emcalzsadc", m_simtwr_cemc_zs_adc);
	    chain.SetBranchAddress("emcalt", m_simtwr_cemc_time);
	}

    chain.SetBranchAddress("sectorih", &m_simtwrmult_ihcal);
    chain.SetBranchAddress("ihcalen", m_simtwr_ihcal_e);
    chain.SetBranchAddress("ihcaletabin", m_simtwr_ihcal_ieta);
    chain.SetBranchAddress("ihcalphibin", m_simtwr_ihcal_iphi);
    chain.SetBranchAddress("ihetacor", m_simtwr_ihcal_eta);
    if (!dataormc) {
    	chain.SetBranchAddress("ihcalzs", m_simtwr_ihcal_zs_e);
	    chain.SetBranchAddress("ihcaladc", m_simtwr_ihcal_adc);
	    chain.SetBranchAddress("ihcalzsadc", m_simtwr_ihcal_zs_adc);
	    chain.SetBranchAddress("ihcalt", m_simtwr_ihcal_time);
	}

    chain.SetBranchAddress("sectoroh", &m_simtwrmult_ohcal);
    chain.SetBranchAddress("ohcalen", m_simtwr_ohcal_e);
    chain.SetBranchAddress("ohcaletabin", m_simtwr_ohcal_ieta);
    chain.SetBranchAddress("ohcalphibin", m_simtwr_ohcal_iphi);
    chain.SetBranchAddress("ohetacor", m_simtwr_ohcal_eta);
    if (!dataormc) {
    	chain.SetBranchAddress("ohcalzs", m_simtwr_ohcal_zs_e);
	    chain.SetBranchAddress("ohcaladc", m_simtwr_ohcal_adc);
	    chain.SetBranchAddress("ohcalzsadc", m_simtwr_ohcal_zs_adc);
	    chain.SetBranchAddress("ohcalt", m_simtwr_ohcal_time);
	}
	
    chain.SetBranchAddress("sectormb", &m_sectormb);
    chain.SetBranchAddress("mbenrgy", &m_mbenergy);

    if (dataormc) {
    	chain.SetBranchAddress("truthpar_n", &m_g4);
	    chain.SetBranchAddress("truthpar_e", m_g4_e);
	    chain.SetBranchAddress("truthpar_eta", m_g4_eta);
	    chain.SetBranchAddress("truthpar_pt", m_g4_pt);
	    chain.SetBranchAddress("truthpar_pz", m_g4_pz);
    }

    chain.SetBranchAddress("track_vtx", m_vtx);

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
		float truthe = 0;

		float E_emcal[96] = {0};
		float E_ihcal[24] = {0};
  		float E_ohcal[24] = {0};

  		float vz_weight = 1.0;

  		eventnumber++;
  		// require that simulation could reconstruct a vertex for the event
  		if(isnan(m_vtx[2])) { continue; }
  		if (m_vtx[2] < minus_z || m_vtx[2] > plus_z) { continue; }

  		float totalcharge = 0.0;
  		for (int i = 0; i < m_sectormb; i++) {
  			totalcharge += m_mbenergy[i];
  		}
  		if (central && totalcharge < centrality_bin[18]) { continue; }
  		h_mbd->Fill(totalcharge);

  		h_vz->Fill(m_vtx[2]);
  		if (reweighting && dataormc) { vz_weight = vertex_reweight[int(floor(m_vtx[2]*2)+100)]; }
  		h_vz_reweight->Fill(m_vtx[2],vz_weight);
  		totalweights += vz_weight;

  		if (dataormc) {
			for (int i = 0; i < m_g4; i++) {
	    		float theta = atan(m_g4_pt[i] / m_g4_pz[i]);
	    		float ET = m_g4_e[i] * abs(sin(theta));
	    		hetdeta->Fill(m_g4_eta[i], ET*vz_weight);
	    		hetdeta_emcalbin->Fill(m_g4_eta[i], ET*vz_weight);
	    		hetdeta_ihcalbin->Fill(m_g4_eta[i], ET*vz_weight);
	    		hetdeta_ohcalbin->Fill(m_g4_eta[i], ET*vz_weight);
	    		if (fabs(m_g4_eta[i]) <= 1.1) {
	    			truthe += ET;
	    			hetdeta_zoom->Fill(m_g4_eta[i], ET*vz_weight);
	    		} 
			}
			h_event_truth_energy->Fill(truthe,vz_weight);
		}
		
		for (int i = 0; i < m_simtwrmult_cemc; i++) {
			if (m_simtwr_cemc_ieta[i] < 8) { continue; }
			std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i]);
		    auto it = emcal_hot_dead_map.find(hot_tower);
		    if (it != emcal_hot_dead_map.end()) { continue; }
			if (!dataormc && zs == 1 && m_simtwr_cemc_zs_adc[i] < zs_value) {
				float zs_e = 0; // = m_simtwr_cemc_zs_e[i]; 
				if (em_zs_calib[m_simtwr_cemc_ieta[i]][m_simtwr_cemc_iphi[i]] == 0) {
					h_em_zero_zscross->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i]);
				} else {
					zs_e = m_simtwr_cemc_zs_e[i]/em_zs_calib[m_simtwr_cemc_ieta[i]][m_simtwr_cemc_iphi[i]];
				}
				h_2D_emcal_calib->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], zs_e*vz_weight);
				h_2D_emcal_calibT->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], zs_e*vz_weight/cosh(m_simtwr_cemc_eta[i]));
				emcale += zs_e/cosh(m_simtwr_cemc_eta[i]); 
				h_emcal->Fill(zs_e/cosh(m_simtwr_cemc_eta[i]), vz_weight);
				h_eT_emcal->Fill(m_simtwr_cemc_eta[i],zs_e*vz_weight/cosh(m_simtwr_cemc_eta[i]));
			} else if (!dataormc && zs == 2 && m_simtwr_cemc_zs_adc[i] < zs_value && m_simtwr_cemc_zs_adc[i] > -1.0*zs_value) {
				//float zs_e = m_simtwr_cemc_zs_e[i]; 
				float zs_e = m_simtwr_cemc_zs_e[i]*em_zs_calib[m_simtwr_cemc_ieta[i]][m_simtwr_cemc_iphi[i]];
				h_2D_emcal_calib->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], zs_e*vz_weight);
				h_2D_emcal_calibT->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], zs_e*vz_weight/cosh(m_simtwr_cemc_eta[i]));
				emcale += zs_e/cosh(m_simtwr_cemc_eta[i]); 
				h_emcal->Fill(zs_e/cosh(m_simtwr_cemc_eta[i]), vz_weight);
				h_eT_emcal->Fill(m_simtwr_cemc_eta[i],zs_e*vz_weight/cosh(m_simtwr_cemc_eta[i]));
			} else {
				if (!dataormc && time && (m_simtwr_cemc_time[i] < -2 || m_simtwr_cemc_time[i] > 2)) { continue; }
				h_2D_emcal_calib->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight);
				h_2D_emcal_calibT->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
				emcale += m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]); 
				h_emcal->Fill(m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]), vz_weight);
				h_eT_emcal->Fill(m_simtwr_cemc_eta[i],m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
			}
		}

		for (int i = 0; i < m_simtwrmult_ihcal; i++) {
			std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i]);
		    auto it = ihcal_hot_dead_map.find(hot_tower);
		    if (it != ihcal_hot_dead_map.end()) { continue; }
		    if (!dataormc && zs == 1 && m_simtwr_ihcal_zs_adc[i] < zs_value) { 
		    	float zs_e = 0; // m_simtwr_ihcal_zs_e[i]; 
				if (ih_zs_calib[m_simtwr_ihcal_ieta[i]][m_simtwr_ihcal_iphi[i]] == 0) {
					h_ih_zero_zscross->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i]);
				} else {
					zs_e = m_simtwr_ihcal_zs_e[i]/ih_zs_calib[m_simtwr_ihcal_ieta[i]][m_simtwr_ihcal_iphi[i]];
				}
				h_2D_ihcal_calib->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], zs_e*vz_weight);
				h_2D_ihcal_calibT->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], zs_e*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
				ihcale += zs_e/cosh(m_simtwr_ihcal_eta[i]); 
				h_ihcal->Fill(zs_e/cosh(m_simtwr_ihcal_eta[i]), vz_weight);
				h_eT_ihcal->Fill(m_simtwr_ihcal_eta[i],zs_e*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
		    } else if (!dataormc && zs == 2 && m_simtwr_ihcal_zs_adc[i] < zs_value && m_simtwr_ihcal_zs_adc[i] > -1.0*zs_value) {
		    	//float zs_e = m_simtwr_ihcal_zs_e[i]; 
				float zs_e = m_simtwr_ihcal_zs_e[i]*ih_zs_calib[m_simtwr_ihcal_ieta[i]][m_simtwr_ihcal_iphi[i]];
				h_2D_ihcal_calib->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], zs_e*vz_weight);
				h_2D_ihcal_calibT->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], zs_e*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
				ihcale += zs_e/cosh(m_simtwr_ihcal_eta[i]); 
				h_ihcal->Fill(zs_e/cosh(m_simtwr_ihcal_eta[i]), vz_weight);
				h_eT_ihcal->Fill(m_simtwr_ihcal_eta[i],zs_e*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
			} else {
				h_2D_ihcal_calib->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]*vz_weight);
				h_2D_ihcal_calibT->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
				ihcale += m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]); 
				h_ihcal->Fill(m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]), vz_weight);
				h_eT_ihcal->Fill(m_simtwr_ihcal_eta[i],m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
			}
		}

		for (int i = 0; i < m_simtwrmult_ohcal; i++) {
			std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i]);
		    auto it = ohcal_hot_dead_map.find(hot_tower);
		    if (it != ohcal_hot_dead_map.end()) { continue; }
		    if (!dataormc && zs == 1 && m_simtwr_ohcal_zs_adc[i] < zs_value) {
		    	float zs_e = 0; // m_simtwr_ihcal_zs_e[i]; 
				if (oh_zs_calib[m_simtwr_ohcal_ieta[i]][m_simtwr_ohcal_iphi[i]] == 0) {
					h_oh_zero_zscross->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i]);
				} else {
					zs_e = m_simtwr_ohcal_zs_e[i]/oh_zs_calib[m_simtwr_ohcal_ieta[i]][m_simtwr_ohcal_iphi[i]];
				}
		    	h_2D_ohcal_calib->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], zs_e*vz_weight);
				h_2D_ohcal_calibT->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], zs_e*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
				ohcale += zs_e/cosh(m_simtwr_ohcal_eta[i]); 
				h_ohcal->Fill(zs_e/cosh(m_simtwr_ohcal_eta[i]), vz_weight);
				h_eT_ohcal->Fill(m_simtwr_ohcal_eta[i],zs_e*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
		    } else if (!dataormc && zs == 2 && m_simtwr_ohcal_zs_adc[i] < zs_value && m_simtwr_ohcal_zs_adc[i] > 1.0*zs_value) {
		    			    	//float zs_e = m_simtwr_ihcal_zs_e[i]; 
				float zs_e = m_simtwr_ohcal_zs_e[i]*oh_zs_calib[m_simtwr_ohcal_ieta[i]][m_simtwr_ohcal_iphi[i]];
		    	h_2D_ohcal_calib->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], zs_e*vz_weight);
				h_2D_ohcal_calibT->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], zs_e*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
				ohcale += zs_e/cosh(m_simtwr_ohcal_eta[i]); 
				h_ohcal->Fill(zs_e/cosh(m_simtwr_ohcal_eta[i]), vz_weight);
				h_eT_ohcal->Fill(m_simtwr_ohcal_eta[i],zs_e*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
			} else {
				h_2D_ohcal_calib->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]*vz_weight);
				h_2D_ohcal_calibT->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
				ohcale += m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]); 
				h_ohcal->Fill(m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]), vz_weight);
				h_eT_ohcal->Fill(m_simtwr_ohcal_eta[i],m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
			}
		}
		
		totale = ihcale + ohcale;
		h_event_energy->Fill(totale + emcale, vz_weight);
		h_event_hcal_energy->Fill(totale, vz_weight);
		h_event_emcal_energy->Fill(emcale, vz_weight);
		h_event_ihcal_energy->Fill(ihcale, vz_weight);
		h_event_ohcal_energy->Fill(ohcale, vz_weight);
		
	}

	if (dataormc) {
		hetdeta->Scale(1.0/totalweights);
		hetdeta_zoom->Scale(1.0/totalweights);
		hetdeta->Scale(1.0/0.1);
		hetdeta_zoom->Scale(1.0/0.01);
		h_event_truth_energy->Scale(1.0/totalweights);
		hetdeta_ihcalbin->Scale(1.0/totalweights);
		hetdeta_ohcalbin->Scale(1.0/totalweights);
	}

	h_2D_emcal_calib->Scale(1.0/totalweights);
	h_2D_emcal_calibT->Scale(1.0/totalweights);
	h_emcal->Scale(1.0/totalweights);
	h_eT_emcal->Scale(1.0/totalweights);
	for (int i = 1; i <= emcal_num_bins; ++i) {
        double bin_width = emcal_bin_edges[i] - emcal_bin_edges[i - 1];
        h_eT_emcal->SetBinContent(i, h_eT_emcal->GetBinContent(i) / bin_width);
        if (dataormc) hetdeta_emcalbin->SetBinContent(i, hetdeta_emcalbin->GetBinContent(i)/ bin_width);
    }

	h_2D_ihcal_calib->Scale(1.0/totalweights);
	h_2D_ihcal_calibT->Scale(1.0/totalweights);
	h_ihcal->Scale(1.0/totalweights);
	h_eT_ihcal->Scale(1.0/totalweights);
	for (int i = 1; i <= ihcal_num_bins; ++i) {
        double bin_width = ihcal_bin_edges[i] - ihcal_bin_edges[i - 1];
        h_eT_ihcal->SetBinContent(i, h_eT_ihcal->GetBinContent(i) / bin_width);
        if (dataormc) hetdeta_ihcalbin->SetBinContent(i, hetdeta_ihcalbin->GetBinContent(i)/ bin_width);
    }

	h_2D_ohcal_calib->Scale(1.0/totalweights);
	h_2D_ohcal_calibT->Scale(1.0/totalweights);
	h_ohcal->Scale(1.0/totalweights);
	h_eT_ohcal->Scale(1.0/totalweights);
	for (int i = 1; i <= ohcal_num_bins; ++i) {
        double bin_width = ohcal_bin_edges[i] - ohcal_bin_edges[i - 1];
        h_eT_ohcal->SetBinContent(i, h_eT_ohcal->GetBinContent(i) / bin_width);
        if (dataormc) hetdeta_ohcalbin->SetBinContent(i, hetdeta_ohcalbin->GetBinContent(i) / bin_width);
    }

	h_event_energy->Scale(1.0/totalweights);
	h_event_hcal_energy->Scale(1.0/totalweights);
	h_event_emcal_energy->Scale(1.0/totalweights);
	h_event_ihcal_energy->Scale(1.0/totalweights);
	h_event_ohcal_energy->Scale(1.0/totalweights);
	
	out->Write();
	out->Close();


}
