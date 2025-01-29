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

std::set<std::tuple<int, int>> emcal_hot_dead_map = {{7,27}};
std::set<std::tuple<int, int>> ihcal_hot_dead_map = {{8,32},{7,51},{10,18},{10,19},{15,61}};
std::set<std::tuple<int, int>> ohcal_hot_dead_map = {{6,42},{7,33},{15,20},{15,23},{16,20}};

double emcal_eta_bin_centers[24];
double ihcal_eta_bin_centers[24];
double ohcal_eta_bin_centers[24];
double calo_eta_bin_centers[24];
double hcal_eta_bin_centers[24];

double vertex_reweight[200] = {1.0};
std::vector<float> centrality_bin;
std::vector<float> tight_centrality_bin;

const double eta_bin_centers[24] = {-1.05417,-0.9625,-0.870833,-0.779167,-0.6875,-0.595833,-0.504167,-0.4125,-0.320833,-0.229167,
	-0.1375,-0.0458333,0.0458333,0.1375,0.229167,0.320833,0.4125,0.504167,0.595833,0.6875,0.779167,0.870833,0.9625,1.05417};

const int cemcSize = 24576;
const int ihcalSize = 1536;
const int ohcalSize = 1536;
const int mbdSize = 256;
const int g4Size = 100000;
const int vtxSize = 3;

const int nruns = 1;
int good_runs[nruns] = {54911}; // [54911, 54914]
int good_run_length[nruns] = {511}; // [511, 447]

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

	TFile *hotdeadfile = new TFile("/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/run_by_run_syst/run54911_hotdeadmap_z_-10_10_ana450_2024p009_fixed_build.root", "READ"); 
	TTree *hotdeadtree = dynamic_cast<TTree*>(hotdeadfile->Get("T"));
	hotdeadtree->SetBranchAddress("emcal_hot_dead_ieta", &emcal_hot_dead_ieta);
    hotdeadtree->SetBranchAddress("emcal_hot_dead_iphi", &emcal_hot_dead_iphi);
    hotdeadtree->SetBranchAddress("ihcal_hot_dead_ieta", &ihcal_hot_dead_ieta);
    hotdeadtree->SetBranchAddress("ihcal_hot_dead_iphi", &ihcal_hot_dead_iphi);
    hotdeadtree->SetBranchAddress("ohcal_hot_dead_ieta", &ohcal_hot_dead_ieta);
    hotdeadtree->SetBranchAddress("ohcal_hot_dead_iphi", &ohcal_hot_dead_iphi);

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
	hotdeadfile->Close();

	TFile *etabinfile = new TFile("/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/run_by_run_syst/run54911_hotdeadmap_z_-10_10_ana450_2024p009_fixed_build.root", "READ"); 
	TTree *etabintree = dynamic_cast<TTree*>(etabinfile->Get("T"));
	etabintree->SetBranchAddress("ihcal_eta_bin_centers", &ihcal_eta_bins);
    etabintree->SetBranchAddress("ohcal_eta_bin_centers", &ohcal_eta_bins);
	etabintree->GetEntry(0);
	for (int i = 0; i < ihcal_eta_bins->size(); i++) {
		ihcal_eta_bin_centers[i] = (*ihcal_eta_bins)[i];
		emcal_eta_bin_centers[i] = (*ihcal_eta_bins)[i];
		calo_eta_bin_centers[i] = (*ihcal_eta_bins)[i];
		hcal_eta_bin_centers[i] = (*ihcal_eta_bins)[i];
	}
	for (int i = 0; i < ohcal_eta_bins->size(); i++) {
		ohcal_eta_bin_centers[i] = (*ihcal_eta_bins)[i]; // EDITED TO USE THE SAME BINNING FOR ALL CALOS
	}

	etabinfile->Close();

}

float mc_centrality(int centbin) {
	int centrality_index = centbin/5;
	return centrality_bin[19-centrality_index];
}

float mc_tight_centrality(int centbin) {
	return tight_centrality_bin[99-centbin];
}

void fill_zvertex_centrality(int dataormc, const char* generator) {
	
	TFile *zvertexfile;
	float vz_reweight[200]; 
	float mc_cent[20];
	float data_cent[20];
	float mc_tight_cent[100]; 
	if (dataormc) {
		//zvertexfile = new TFile(TString::Format("/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/dETdeta_vertex_reweight_run54912_%s_new_2024p007.root", generator), "READ");
		zvertexfile = new TFile(TString::Format("/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/run_by_run_syst/dETdeta_vertex_reweight_run54911_%s_ana450_2024p009_fixed_build.root", generator), "READ");
		TTree *vertextree = dynamic_cast<TTree*>(zvertexfile->Get("T"));
		vertextree->SetBranchAddress("vertex_reweight", vz_reweight);
		vertextree->SetBranchAddress("mc_centrality", mc_cent);
		vertextree->SetBranchAddress("mc_centrality_tight", mc_tight_cent); 
		vertextree->GetEntry(0);
		for (int i = 0; i < 200; i++) {
			vertex_reweight[i] = vz_reweight[i];
		}
		centrality_bin.assign(mc_cent, mc_cent+20);
		tight_centrality_bin.assign(mc_tight_cent, mc_tight_cent+100); 
	} else { 
		//zvertexfile = new TFile("/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/dETdeta_vertex_reweight_run54912_reweight_hijing_new_2024p007.root", "READ");
		zvertexfile = new TFile("/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/run_by_run_syst/dETdeta_vertex_reweight_run54911_reweight_epos_2024_ana450_2024p009_fixed_build.root", "READ");
		TTree *vertextree = dynamic_cast<TTree*>(zvertexfile->Get("T"));
		//vertextree->SetBranchAddress("data_centrality", data_cent); // changed from data_centrality to use new centrality bins
		//vertextree->GetEntry(0);
		float data_central[] = {0.0,0.0,3.4412263484727412,9.090884292111504,16.905063169580032,28.022382830087857,44.265968300738564,67.84341228546377,100.2493262652662,143.30257633270628,198.03245290421535,265.35366129472857,346.1799984953793,441.36798663537985,552.532564809136,682.210081562082,833.0489753725108,1011.1971134894196,1223.2265691250534,1484.3956409201726};
		centrality_bin.assign(data_central, data_central+20);
	}
	zvertexfile->Close();

}

void dETdeta_multi_run_analysis(const char* generator = "", float minus_z = -2, float plus_z = 2, int dataormc = 0, int reweighting = 1, int central = 0, int min_cent = 0, int max_cent = 5, int zs = 0, int zs_value = 10, int time = 0, const char* emsyst = "", const char* ihsyst = "", const char* ohsyst = "", const char* opt_tag = "", const char* subdirectory = "", int closure = 0) {

	string filename = "/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/run_by_run_syst/";
	if (!strcmp(subdirectory,"")) {
		filename += "dETdeta_analysis_allruns"; 
	} else {
		string dirtag = subdirectory;
		filename += dirtag + "/dETdeta_analysis_allruns";
	}
	string opttag = opt_tag;
	if (strcmp(opt_tag,"")) { filename += "_" + opttag; }
	string emsysttag = emsyst;
	if (strcmp(emsyst,"")) { filename += "_emsyst" + emsysttag; }
	string ihsysttag = ihsyst;
	if (strcmp(ihsyst,"")) { filename += "_ihsyst" + ihsysttag; }
	string ohsysttag = ohsyst;
	if (strcmp(ohsyst,"")) { filename += "_ohsyst" + ohsysttag; }
	
	string zstag; 
	if (zs == 0) { zstag = "nozs"; }
	if (zs == 1) { zstag = "HCal_zs_" + to_string(zs_value) + "ADC_EMCal_zs_" + to_string(zs_value + 10) + "ADC"; }
	if (zs == 2) { zstag = "zs_ecut" + to_string(zs_value) + "ADC"; }
	string timetag = (time?"emcal_timecut_":"");
	string closuretag;
	if (!dataormc || closure == 0) { closuretag = ""; }
	if (dataormc && closure == 1) { closuretag = "half_closure_correction_"; }
	if (dataormc && closure == 2) { closuretag = "half_closure_dataset_"; }
  	string dattag = (dataormc?"mc":"data");
  	string weighttag = (reweighting?"reweight":"noweight");
  	string centtag = (central?to_string(min_cent)+"-"+to_string(max_cent):"0-90"); 
  	string gentag = generator;
  	if (!strcmp(generator,"")) filename += "_" + zstag + "_" + closuretag + timetag + dattag + "_" + weighttag + "_" + centtag + ".root";
  	else filename += "_" + zstag + "_" + closuretag + timetag + dattag + "_" + weighttag + "_" + centtag + "_" + gentag + ".root";

	fill_hot_dead_map_eta_bin_centers();
	fill_zvertex_centrality(dataormc, generator);
	assert(centrality_bin.size() == 20);
	if (dataormc) assert(tight_centrality_bin.size() == 100);

	std::cout << filename << std::endl;

	std::cout << "centrality_bins" << std::endl;
	for (int i = 0; i < 20; i++) {
		std::cout << centrality_bin[i] << " ";
	}
	std::cout << std::endl;
	if (dataormc) {
		std::cout << "tight_centrality_bins" << std::endl;
		for (int i = 0; i < 100; i++) {
			std::cout << tight_centrality_bin[i] << " ";
		}
		std::cout << std::endl;
	}

	TFile *out = new TFile(filename.c_str(),"RECREATE");
	TH1F* h_vz = new TH1F("h_vz","",400, -100, 100);
	TH1F* h_vz_reweight = new TH1F("h_vz_reweight","",400, -100, 100);
	TH1F* h_mbd = new TH1F("h_mbd","",250000,0,250000);
	TH1F* h_npart = new TH1F("h_npart","",500,0,500);
	TH1F* h_cent = new TH1F("h_cent","",100,-0.5,99.5);

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
    //TH1F* h_eT_emcal = new TH1F("h_eT_emcal","",emcal_num_bins,emcal_bin_edges);
	TH1F* hetdeta_emcalbin = new TH1F("hetdeta_emcalbin","",emcal_num_bins, emcal_bin_edges);
	//TH1F* h_eT_eta_emcal = new TH1F("h_eT_eta_emcal","",emcal_num_bins, emcal_bin_edges);
	TProfile* h_eT_eta_emcal_profile = new TProfile("h_eT_eta_emcal_profile","",emcal_num_bins, emcal_bin_edges);
	TH1F* h_eT_eta_emcal_profile_hist = new TH1F("h_eT_eta_emcal_profile_hist","",emcal_num_bins, emcal_bin_edges);
	TProfile* h_eta_emcal_profile = new TProfile("h_eta_emcal_profile","",emcal_num_bins, emcal_bin_edges);

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
    //TH1F* h_eT_ihcal = new TH1F("h_eT_ihcal","",ihcal_num_bins,ihcal_bin_edges);
	TH1F* hetdeta_ihcalbin = new TH1F("hetdeta_ihcalbin","",ihcal_num_bins, ihcal_bin_edges);
	//TH1F* h_eT_eta_ihcal = new TH1F("h_eT_eta_ihcal","",ihcal_num_bins, ihcal_bin_edges);
	TProfile* h_eT_eta_ihcal_profile = new TProfile("h_eT_eta_ihcal_profile","",ihcal_num_bins, ihcal_bin_edges);
	TH1F* h_eT_eta_ihcal_profile_hist = new TH1F("h_eT_eta_ihcal_profile_hist","",ihcal_num_bins, ihcal_bin_edges);
	TProfile* h_eta_ihcal_profile = new TProfile("h_eta_ihcal_profile","",ihcal_num_bins, ihcal_bin_edges);
	
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
	//TH1F* h_eT_ohcal = new TH1F("h_eT_ohcal","",ohcal_num_bins,ohcal_bin_edges);
	TH1F* hetdeta_ohcalbin = new TH1F("hetdeta_ohcalbin","",ohcal_num_bins, ohcal_bin_edges);
	//TH1F* h_eT_eta_ohcal = new TH1F("h_eT_eta_ohcal","",ohcal_num_bins, ohcal_bin_edges);
	TProfile* h_eT_eta_ohcal_profile = new TProfile("h_eT_eta_ohcal_profile","",ohcal_num_bins, ohcal_bin_edges);
	TH1F* h_eT_eta_ohcal_profile_hist = new TH1F("h_eT_eta_ohcal_profile_hist","",ohcal_num_bins, ohcal_bin_edges);
	TProfile* h_eta_ohcal_profile = new TProfile("h_eta_ohcal_profile","",ohcal_num_bins, ohcal_bin_edges);
	
	int calo_num_bins = 24;
    double calo_bin_edges[25];
    if (!dataormc || reweighting) {
    	calo_bin_edges[0] = calo_eta_bin_centers[0] - 0.5 * (calo_eta_bin_centers[1] - calo_eta_bin_centers[0]);
	    for (int i = 1; i < calo_num_bins; ++i) { calo_bin_edges[i] = (calo_eta_bin_centers[i] + calo_eta_bin_centers[i - 1]) / 2.0; }
	    calo_bin_edges[calo_num_bins] = calo_eta_bin_centers[calo_num_bins - 1] + 0.5 * (calo_eta_bin_centers[calo_num_bins - 1] - calo_eta_bin_centers[calo_num_bins - 2]);
    } else {
    	calo_bin_edges[0] = eta_bin_centers[0] - 0.5 * (eta_bin_centers[1] - eta_bin_centers[0]);
    	for (int i = 1; i < calo_num_bins; ++i) { calo_bin_edges[i] = (eta_bin_centers[i] + eta_bin_centers[i - 1]) / 2.0; }
    	calo_bin_edges[calo_num_bins] = eta_bin_centers[calo_num_bins - 1] + 0.5 * (eta_bin_centers[calo_num_bins - 1] - eta_bin_centers[calo_num_bins - 2]);
    }
    //TH1F* h_eT_calo = new TH1F("h_eT_calo","",calo_num_bins,calo_bin_edges);
	TH1F* hetdeta_calobin = new TH1F("hetdeta_calobin","",calo_num_bins, calo_bin_edges);
	//TH1F* h_eT_eta_calo = new TH1F("h_eT_eta_calo","",calo_num_bins, calo_bin_edges);
	TProfile* h_eT_eta_calo_profile = new TProfile("h_eT_eta_calo_profile","",calo_num_bins, calo_bin_edges);
	TH1F* h_eT_eta_calo_profile_hist = new TH1F("h_eT_eta_calo_profile_hist","",calo_num_bins, calo_bin_edges);
	TProfile* h_eta_calo_profile = new TProfile("h_eta_calo_profile","",calo_num_bins, calo_bin_edges);

	int hcal_num_bins = 24;
    double hcal_bin_edges[25];
    if (!dataormc || reweighting) {
    	hcal_bin_edges[0] = hcal_eta_bin_centers[0] - 0.5 * (hcal_eta_bin_centers[1] - hcal_eta_bin_centers[0]);
	    for (int i = 1; i < hcal_num_bins; ++i) { hcal_bin_edges[i] = (hcal_eta_bin_centers[i] + hcal_eta_bin_centers[i - 1]) / 2.0; }
	    hcal_bin_edges[hcal_num_bins] = hcal_eta_bin_centers[hcal_num_bins - 1] + 0.5 * (hcal_eta_bin_centers[hcal_num_bins - 1] - hcal_eta_bin_centers[hcal_num_bins - 2]);
    } else {
    	hcal_bin_edges[0] = eta_bin_centers[0] - 0.5 * (eta_bin_centers[1] - eta_bin_centers[0]);
    	for (int i = 1; i < hcal_num_bins; ++i) { hcal_bin_edges[i] = (eta_bin_centers[i] + eta_bin_centers[i - 1]) / 2.0; }
    	hcal_bin_edges[hcal_num_bins] = eta_bin_centers[hcal_num_bins - 1] + 0.5 * (eta_bin_centers[hcal_num_bins - 1] - eta_bin_centers[hcal_num_bins - 2]);
    }
    //TH1F* h_eT_hcal = new TH1F("h_eT_hcal","",hcal_num_bins,hcal_bin_edges);
	TH1F* hetdeta_hcalbin = new TH1F("hetdeta_hcalbin","",hcal_num_bins, hcal_bin_edges);
	//TH1F* h_eT_eta_hcal = new TH1F("h_eT_eta_hcal","",hcal_num_bins, hcal_bin_edges);
	TProfile* h_eT_eta_hcal_profile = new TProfile("h_eT_eta_hcal_profile","",hcal_num_bins, hcal_bin_edges);
	TH1F* h_eT_eta_hcal_profile_hist = new TH1F("h_eT_eta_hcal_profile_hist","",hcal_num_bins, hcal_bin_edges);
	TProfile* h_eta_hcal_profile = new TProfile("h_eta_hcal_profile","",hcal_num_bins, hcal_bin_edges);

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
    } else if (dataormc && !strcmp(generator, "reweight_hijing")) {
	    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
	    for (int i = 0; i < 2500; i++) { 
    		TString wildcardPath = TString::Format("%sevents_hijing_reweighted_phenix_run14_mc_cor_%d.root", inputDirectory, i); 
	    	//TString wildcardPath = TString::Format("%sevents_hijing_reweighted_brahms_run14_mc_cor_%d.root", inputDirectory, i); 
	    	//TString wildcardPath = TString::Format("%sevents_hijing_reweighted_brahms_run14_eta_-6_6_mc_cor_%d.root", inputDirectory, i);
	    	//TString wildcardPath = TString::Format("%sevents_hijing_reweighted_phenix_run14_eta_-6_6_mc_cor_%d.root", inputDirectory, i);
	    	chain.Add(wildcardPath);
	    }
	} else if (dataormc && !strcmp(generator, "ampt")) {
    	// location of AMPT files 
	    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/amptrun/plugdoor_condor/";
	    for (int i = 0; i < 1000; i++) {
	    	TString wildcardPath = TString::Format("%sOutDir%d/ampt_sim_output.root", inputDirectory , i);
	    	chain.Add(wildcardPath);
	    }
    } else if (dataormc && !strcmp(generator, "reweight_ampt")) {
	    const char* inputDirectory = "/sphenix/tg/tg01/commissioning/CaloCalibWG/egm2153/detdeta_datasets/"; 
		for (int i = 0; i < 3000; i++) { 
    		TString wildcardPath = TString::Format("%sevents_ampt_reweighted_run14_ampt_cor_%d.root", inputDirectory, i);
	    	chain.Add(wildcardPath);
	    }
	} else if (dataormc && !strcmp(generator, "reweight_epos")) {
	    const char* inputDirectory = "/sphenix/tg/tg01/commissioning/CaloCalibWG/egm2153/detdeta_datasets/"; 
		for (int i = 0; i < 3000; i++) { 
    		TString wildcardPath = TString::Format("%sevents_epos_reweighted_run14_epos_cor_%d.root", inputDirectory, i);
	    	chain.Add(wildcardPath);
	    }  
	} else if (dataormc && !strcmp(generator, "reweight_ampt_2024")) {
	    const char* inputDirectory = "/sphenix/tg/tg01/commissioning/CaloCalibWG/egm2153/detdeta_run24auau/"; 
	    int data_start = 0; int data_end = 5000;
	    if (closure == 1) { data_end = 2500; }
	    if (closure == 2) { data_start = 2500; }
		for (int i = data_start; i < data_end; i++) { 
    		TString wildcardPath = TString::Format("%sevents_ampt_reweighted_run14_fixed_build_ampt_cor_%d.root", inputDirectory, i);
	    	chain.Add(wildcardPath);
	    }
	} else if (dataormc && !strcmp(generator, "reweight_epos_2024")) {
	    const char* inputDirectory = "/sphenix/tg/tg01/commissioning/CaloCalibWG/egm2153/detdeta_run24auau/"; 
	    int data_start = 0; int data_end = 5000;
	    if (closure == 1) { data_end = 2500; }
	    if (closure == 2) { data_start = 2500; }
		for (int i = data_start; i < data_end; i++) { 
    		TString wildcardPath = TString::Format("%sevents_epos_reweighted_run14_fixed_build_epos_cor_%d.root", inputDirectory, i);
	    	chain.Add(wildcardPath);
	    }
	} else if (dataormc && !strcmp(generator, "reweight_hijing_2024")) {
	    const char* inputDirectory = "/sphenix/tg/tg01/commissioning/CaloCalibWG/egm2153/detdeta_run24auau/";
	    int data_start = 0; int data_end = 5000;
	    if (closure == 1) { data_end = 2500; }
	    if (closure == 2) { data_start = 2500; }
	    for (int i = data_start; i < data_end; i++) { 
    		TString wildcardPath = TString::Format("%sevents_hijing_reweighted_run14_fixed_build_mc_cor_%d.root", inputDirectory, i); 
	    	chain.Add(wildcardPath);
	    }
	} else if (!dataormc) { 
    	const char* inputDirectory = "/sphenix/tg/tg01/commissioning/CaloCalibWG/egm2153/detdeta_run24auau/";
    	for (int r = 0; r < nruns; r++) {
            for (int s = 0; s < good_run_length[r]; s++) {
            	//TString wildcardPath = TString::Format("%sevents_new_calib_12_12_24_ana450_2024p009_%d_100_50_50_zs_data_cor_%d.root", inputDirectory, good_runs[r], s); 
                //TString wildcardPath = TString::Format("%sevents_ana450_2024p009_%d_100_50_50_zs_hcal_tsc_data_cor_%d.root", inputDirectory, good_runs[r], s);
                TString wildcardPath = TString::Format("%sevents_ana450_2024p009_%d_fixed_build_data_cor_%d.root", inputDirectory, good_runs[r], s);
                //TString wildcardPath = TString::Format("%snew_calib_12_12_24_trig10_events_withzscc_ana450_2024p009_%d_data_cor_%d.root", inputDirectory, good_runs[r], s); 
                chain.Add(wildcardPath); 
            }
        }
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
    float m_mbtime[mbdSize];
    float m_mbdtotalq;
    int m_g4;
    float m_g4_e[g4Size];
    float m_g4_eta[g4Size];
    float m_g4_pt[g4Size];
    float m_g4_pz[g4Size];
    float m_vtx[vtxSize];
    int m_npart;
    bool m_isMinBias;
    int m_centbin;
    float m_simtwr_cemc_chi2[cemcSize];
    float m_simtwr_ihcal_chi2[ihcalSize];
    float m_simtwr_ohcal_chi2[ohcalSize];

     // Set branch addresses
    chain.SetBranchAddress("sectorem", &m_simtwrmult_cemc);
    if (dataormc || !strcmp(emsyst,"")) {
    	chain.SetBranchAddress("emcalen", m_simtwr_cemc_e);
    } else {
    	string emsystbranch = "emcalsyst" + emsysttag;
    	chain.SetBranchAddress(emsystbranch.c_str(), m_simtwr_cemc_e);
    }
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
    if (dataormc || !strcmp(ihsyst,"")) {
    	chain.SetBranchAddress("ihcalen", m_simtwr_ihcal_e);
    } else {
    	string ihsystbranch = "ihcalsyst" + ihsysttag;
    	chain.SetBranchAddress(ihsystbranch.c_str(), m_simtwr_ihcal_e);
    }
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
    if (dataormc || !strcmp(ohsyst,"")) {
    	chain.SetBranchAddress("ohcalen", m_simtwr_ohcal_e);
    } else {
    	string ohsystbranch = "ohcalsyst" + ohsysttag;
    	chain.SetBranchAddress(ohsystbranch.c_str(), m_simtwr_ohcal_e);
    }
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
    chain.SetBranchAddress("mbenrgy", m_mbenergy);
    chain.SetBranchAddress("mbd_total_q", &m_mbdtotalq); // moved to use old MC
   	chain.SetBranchAddress("mbtime", m_mbtime); // moved to use old MC
    if (!dataormc) {
    	chain.SetBranchAddress("isMinBias", &m_isMinBias);
	    chain.SetBranchAddress("centbin", &m_centbin);
    }
    if (dataormc) {
    	chain.SetBranchAddress("npart",&m_npart);
    	chain.SetBranchAddress("truthpar_n", &m_g4);
	    chain.SetBranchAddress("truthpar_e", m_g4_e);
	    chain.SetBranchAddress("truthpar_eta", m_g4_eta);
	    chain.SetBranchAddress("truthpar_pt", m_g4_pt);
	    chain.SetBranchAddress("truthpar_pz", m_g4_pz);
    }
    chain.SetBranchAddress("track_vtx", m_vtx);
    
    if (!dataormc) {
    	chain.SetBranchAddress("emchi2", m_simtwr_cemc_chi2);
    	chain.SetBranchAddress("ihchi2", m_simtwr_ihcal_chi2);
    	chain.SetBranchAddress("ohchi2", m_simtwr_ohcal_chi2);
    }

	int eventnumber = 0;
	float totalweights = 0.0;

    Long64_t nEntries = chain.GetEntries();
    std::cout << nEntries << std::endl;
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 10000; ++entry) {
        chain.GetEntry(entry);
    	if (eventnumber % 1000 == 0) cout << "event " << eventnumber << endl;

		float emcale = 0;
		float ihcale = 0;
		float ohcale = 0;
		float totale = 0;
		float truthe = 0;
		
		float E_emcalbin[96] = {0};
		float E_ihcal[24] = {0};
  		float E_ohcal[24] = {0};
  		float eta_emcal[96] = {0};
  		float eta_ihcal[24] = {0};
  		float eta_ohcal[24] = {0};

  		float E_emcal[24] = {0};
  		float E_calo[24] = {0};
  		float E_hcal[24] = {0};

  		int emcaleta[24] = {0};
  		int hcaleta[24] = {0};
  		int caloeta[24] = {0};

  		float vz_weight = 1.0;

  		eventnumber++;

  		// require that simulation could reconstruct a vertex for the event
  		if(isnan(m_vtx[2])) { continue; }
  		if (m_vtx[2] < minus_z || m_vtx[2] > plus_z) { continue; }
  		if (!dataormc && !m_isMinBias) { continue; }

  		/*
  		// testing for intermittant hot towers 
		bool good_event = true;
		if (!dataormc) {
			for (int i = 0; i < m_simtwrmult_cemc; i++) {
				if (m_simtwr_cemc_chi2[i] > 10000) { 
					good_event = false;
					break;
				}
				if (i < 1536) {
					if (m_simtwr_ihcal_chi2[i] > 10000) {
						good_event = false;
						break;
					}
					if (m_simtwr_ohcal_chi2[i] > 10000) {
						good_event = false;
						break;
					}

				}
			}
		}
		if (!good_event) { continue; }
		// testing for intermittant hot towers 
		*/

  		int mbd_nhits1 = 0;
    	int mbd_nhits2 = 0;
    	float totalcharge = 0.0;

  		if (dataormc) {
	    	//for (int i = 0; i < 64; i++) {
	    		//if (m_mbenergy[i] > 0.5) { mbd_nhits1 += 1; }
    			//if (m_mbenergy[i+64] > 0.5) { mbd_nhits2 += 1; }
	    	//	if (m_mbenergy[i] > 0.5 && m_mbtime[i] < 25.0) { mbd_nhits1 += 1; } // moved to use old MC
    		//	if (m_mbenergy[i+64] > 0.5 && m_mbtime[i] < 25.0) { mbd_nhits2 += 1; } // moved to use old MC
    		//}
    		//if (mbd_nhits1 < 2 || mbd_nhits2 < 2) { continue; } // MC minimum bias definition 
  			
  			for (int i = 0; i < m_sectormb; i++) {
  				//if (m_mbenergy[i] > 0.5) { totalcharge += m_mbenergy[i]; }
  				if (m_mbenergy[i] > 0.5 && m_mbtime[i] < 25.0) { totalcharge += m_mbenergy[i]; } // moved to use old MC
  			}
  		}

  		if (dataormc && central) {
  			float max_charge = mc_tight_centrality(max_cent);
  			float min_charge = mc_tight_centrality(min_cent);
  			if (totalcharge > mc_tight_centrality(min_cent) || totalcharge < mc_tight_centrality(max_cent)) {
  				continue; 
  			}
  		}

  		if (!dataormc && central && (m_centbin > max_cent || m_centbin < (min_cent + 1))) { continue; } 

  		if (!dataormc) { totalcharge = m_mbdtotalq; }
  		h_mbd->Fill(totalcharge);
  		if (!dataormc) { h_cent->Fill(m_centbin); }

  		h_vz->Fill(m_vtx[2]);
  		if (reweighting && dataormc) { 
  			vz_weight = vertex_reweight[int(floor(m_vtx[2]*2)+100)]; 
  			if (vz_weight == 0) { std::cout << "VZ WEIGHT OF ZERO FOR VZ = " << m_vtx[2] << std::endl; }
  		}
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
	    		hetdeta_calobin->Fill(m_g4_eta[i], ET*vz_weight);
	    		hetdeta_hcalbin->Fill(m_g4_eta[i], ET*vz_weight);
	    		if (fabs(m_g4_eta[i]) <= 1.1) {
	    			truthe += ET;
	    			hetdeta_zoom->Fill(m_g4_eta[i], ET*vz_weight);
	    		} 
			}
			h_event_truth_energy->Fill(truthe,vz_weight);
			h_npart->Fill(m_npart,vz_weight);
		}
		
		for (int i = 0; i < m_simtwrmult_cemc; i++) {
			std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i]);
		    auto it = emcal_hot_dead_map.find(hot_tower);
		    if (it != emcal_hot_dead_map.end()) { continue; }
		    if (!dataormc && m_simtwr_cemc_chi2[i] > 10000) { 
		    	std::cout << "EMCal high chi2 channel " << m_simtwr_cemc_ieta[i] << " " << m_simtwr_cemc_iphi[i] << " energy " << m_simtwr_cemc_e[i] << " chi2 " << m_simtwr_cemc_chi2[i] << std::endl;
		    	continue; 
		    }
		    h_2D_emcal_calib->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight);
			h_2D_emcal_calibT->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
			emcale += m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]); 
			h_emcal->Fill(m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]), vz_weight);
			//h_eT_emcal->Fill(m_simtwr_cemc_eta[i],m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
			//h_eT_calo->Fill(m_simtwr_cemc_eta[i],m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
			E_emcalbin[m_simtwr_cemc_ieta[i]] += m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]);
			eta_emcal[m_simtwr_cemc_ieta[i]] = m_simtwr_cemc_eta[i];
		}
		
		for (int i = 0; i < m_simtwrmult_ihcal; i++) {
			std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i]);
		    auto it = ihcal_hot_dead_map.find(hot_tower);
		    if (it != ihcal_hot_dead_map.end()) { continue; }
		    if (isnan(m_simtwr_ihcal_e[i])) { continue; }
		    if (!dataormc && m_simtwr_ihcal_chi2[i] > 10000) { 
		    	std::cout << "IHCal high chi2 channel " << m_simtwr_ihcal_ieta[i] << " " << m_simtwr_ihcal_iphi[i] << " energy " << m_simtwr_ihcal_e[i] << " chi2 " << m_simtwr_ihcal_chi2[i] << std::endl;
		    	continue; 
		    }
		    h_2D_ihcal_calib->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]*vz_weight);
			h_2D_ihcal_calibT->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
			ihcale += m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]); 
			h_ihcal->Fill(m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]), vz_weight);
			//h_eT_ihcal->Fill(m_simtwr_ihcal_eta[i],m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
			//h_eT_calo->Fill(m_simtwr_ihcal_eta[i],m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
			E_ihcal[m_simtwr_ihcal_ieta[i]] += m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]);
			eta_ihcal[m_simtwr_ihcal_ieta[i]] = m_simtwr_ihcal_eta[i];
		}

		for (int i = 0; i < m_simtwrmult_ohcal; i++) {
			std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i]);
		    auto it = ohcal_hot_dead_map.find(hot_tower);
		    if (it != ohcal_hot_dead_map.end()) { continue; }
		    if (!dataormc && m_simtwr_ohcal_chi2[i] > 10000) { 
		    	std::cout << "OHCal high chi2 channel " << m_simtwr_ohcal_ieta[i] << " " << m_simtwr_ohcal_iphi[i] << " energy " << m_simtwr_ohcal_e[i] << " chi2 " << m_simtwr_ohcal_chi2[i] << std::endl;
		    	continue; 
		    }
		    h_2D_ohcal_calib->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]*vz_weight);
			h_2D_ohcal_calibT->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
			ohcale += m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]); 
			h_ohcal->Fill(m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]), vz_weight);
			//h_eT_ohcal->Fill(m_simtwr_ohcal_eta[i],m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
			//h_eT_calo->Fill(m_simtwr_ohcal_eta[i],m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
			E_ohcal[m_simtwr_ohcal_ieta[i]] += m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]);
			eta_ohcal[m_simtwr_ohcal_ieta[i]] = m_simtwr_ohcal_eta[i];
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
			//h_eT_eta_ihcal->Fill(eta_ihcal[i],E_ihcal[i]*vz_weight);
			//h_eT_eta_ohcal->Fill(eta_ohcal[i],E_ohcal[i]*vz_weight);
			h_eT_eta_ihcal_profile->Fill(eta_ihcal[i],E_ihcal[i],vz_weight);
			h_eT_eta_ohcal_profile->Fill(eta_ohcal[i],E_ohcal[i],vz_weight);
			h_eta_ihcal_profile->Fill(eta_ihcal[i],1);
			h_eta_ohcal_profile->Fill(eta_ohcal[i],1);
			if (E_emcal[i] != 0) {
				//h_eT_eta_emcal->Fill(emcal_eta_bin_centers[i],E_emcal[i]*vz_weight);
				h_eT_eta_emcal_profile->Fill(emcal_eta_bin_centers[i],E_emcal[i],vz_weight);
				h_eta_emcal_profile->Fill(emcal_eta_bin_centers[i],emcaleta[i]);
			}
			if (E_calo[i] != 0) {
				//h_eT_eta_calo->Fill(calo_eta_bin_centers[i],E_calo[i]*vz_weight);
				h_eT_eta_calo_profile->Fill(calo_eta_bin_centers[i],E_calo[i],vz_weight);
				h_eta_calo_profile->Fill(calo_eta_bin_centers[i],caloeta[i]);
			}
			if (E_hcal[i] != 0) {
				//h_eT_eta_hcal->Fill(hcal_eta_bin_centers[i],E_hcal[i]*vz_weight);
				h_eT_eta_hcal_profile->Fill(hcal_eta_bin_centers[i],E_hcal[i],vz_weight);
				h_eta_hcal_profile->Fill(hcal_eta_bin_centers[i],hcaleta[i]);
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
		hetdeta_emcalbin->Scale(1.0/totalweights);
		hetdeta_calobin->Scale(1.0/totalweights);
		hetdeta_hcalbin->Scale(1.0/totalweights);
	}

	h_2D_emcal_calib->Scale(1.0/totalweights);
	h_2D_emcal_calibT->Scale(1.0/totalweights);
	h_emcal->Scale(1.0/totalweights);
	//h_eT_emcal->Scale(1.0/totalweights);
	//h_eT_eta_emcal->Scale(1.0/totalweights);
	for (int i = 1; i <= emcal_num_bins; ++i) {
        double bin_width = emcal_bin_edges[i] - emcal_bin_edges[i - 1];
        //h_eT_emcal->SetBinContent(i, h_eT_emcal->GetBinContent(i) / bin_width);
        //h_eT_emcal->SetBinError(i, h_eT_emcal->GetBinError(i) / bin_width);
        //h_eT_eta_emcal->SetBinContent(i, h_eT_eta_emcal->GetBinContent(i) / bin_width);
        //h_eT_eta_emcal->SetBinError(i, h_eT_eta_emcal->GetBinError(i) / bin_width);
        h_eT_eta_emcal_profile_hist->SetBinContent(i, h_eT_eta_emcal_profile->GetBinContent(i) / bin_width);
        h_eT_eta_emcal_profile_hist->SetBinError(i, h_eT_eta_emcal_profile->GetBinError(i) / bin_width);
        if (dataormc) {
        	hetdeta_emcalbin->SetBinContent(i, hetdeta_emcalbin->GetBinContent(i)/ bin_width);
        	hetdeta_emcalbin->SetBinError(i, hetdeta_emcalbin->GetBinError(i)/ bin_width);
        }
    }

	h_2D_ihcal_calib->Scale(1.0/totalweights);
	h_2D_ihcal_calibT->Scale(1.0/totalweights);
	h_ihcal->Scale(1.0/totalweights);
	//h_eT_ihcal->Scale(1.0/totalweights);
	//h_eT_eta_ihcal->Scale(1.0/totalweights);
	for (int i = 1; i <= ihcal_num_bins; ++i) {
        double bin_width = ihcal_bin_edges[i] - ihcal_bin_edges[i - 1];
        //h_eT_ihcal->SetBinContent(i, h_eT_ihcal->GetBinContent(i) / bin_width);
        //h_eT_ihcal->SetBinError(i, h_eT_ihcal->GetBinError(i) / bin_width);
		//h_eT_eta_ihcal->SetBinContent(i, h_eT_eta_ihcal->GetBinContent(i) / bin_width);
        //h_eT_eta_ihcal->SetBinError(i, h_eT_eta_ihcal->GetBinError(i) / bin_width);
        h_eT_eta_ihcal_profile_hist->SetBinContent(i, h_eT_eta_ihcal_profile->GetBinContent(i) / bin_width);
        h_eT_eta_ihcal_profile_hist->SetBinError(i, h_eT_eta_ihcal_profile->GetBinError(i) / bin_width);
        if (dataormc) {
        	hetdeta_ihcalbin->SetBinContent(i, hetdeta_ihcalbin->GetBinContent(i)/ bin_width);
        	hetdeta_ihcalbin->SetBinError(i, hetdeta_ihcalbin->GetBinError(i)/ bin_width);
        }
    }

	h_2D_ohcal_calib->Scale(1.0/totalweights);
	h_2D_ohcal_calibT->Scale(1.0/totalweights);
	h_ohcal->Scale(1.0/totalweights);
	//h_eT_ohcal->Scale(1.0/totalweights);
	//h_eT_eta_ohcal->Scale(1.0/totalweights);
	for (int i = 1; i <= ohcal_num_bins; ++i) {
        double bin_width = ohcal_bin_edges[i] - ohcal_bin_edges[i - 1];
        //h_eT_ohcal->SetBinContent(i, h_eT_ohcal->GetBinContent(i) / bin_width);
        //h_eT_ohcal->SetBinError(i, h_eT_ohcal->GetBinError(i) / bin_width);
        //h_eT_eta_ohcal->SetBinContent(i, h_eT_eta_ohcal->GetBinContent(i) / bin_width);
        //h_eT_eta_ohcal->SetBinError(i, h_eT_eta_ohcal->GetBinError(i) / bin_width);
        h_eT_eta_ohcal_profile_hist->SetBinContent(i, h_eT_eta_ohcal_profile->GetBinContent(i) / bin_width);
        h_eT_eta_ohcal_profile_hist->SetBinError(i, h_eT_eta_ohcal_profile->GetBinError(i) / bin_width);
        if (dataormc) {
        	hetdeta_ohcalbin->SetBinContent(i, hetdeta_ohcalbin->GetBinContent(i) / bin_width);
        	hetdeta_ohcalbin->SetBinError(i, hetdeta_ohcalbin->GetBinError(i) / bin_width);
        }
    }

    //h_eT_calo->Scale(1.0/totalweights);
    //h_eT_eta_calo->Scale(1.0/totalweights);
	for (int i = 1; i <= calo_num_bins; ++i) {
        double bin_width = calo_bin_edges[i] - calo_bin_edges[i - 1];
        //h_eT_calo->SetBinContent(i, h_eT_calo->GetBinContent(i) / bin_width);
        //h_eT_calo->SetBinError(i, h_eT_calo->GetBinError(i) / bin_width);
        //h_eT_eta_calo->SetBinContent(i, h_eT_eta_calo->GetBinContent(i) / bin_width);
        //h_eT_eta_calo->SetBinError(i, h_eT_eta_calo->GetBinError(i) / bin_width);
        h_eT_eta_calo_profile_hist->SetBinContent(i, h_eT_eta_calo_profile->GetBinContent(i) / bin_width);
        h_eT_eta_calo_profile_hist->SetBinError(i, h_eT_eta_calo_profile->GetBinError(i) / bin_width);
        if (dataormc) {
        	hetdeta_calobin->SetBinContent(i, hetdeta_calobin->GetBinContent(i) / bin_width);
        	hetdeta_calobin->SetBinError(i, hetdeta_calobin->GetBinError(i) / bin_width);
        }
    }

    //h_eT_hcal->Scale(1.0/totalweights);
    //h_eT_eta_hcal->Scale(1.0/totalweights);
	for (int i = 1; i <= hcal_num_bins; ++i) {
        double bin_width = hcal_bin_edges[i] - hcal_bin_edges[i - 1];
        //h_eT_hcal->SetBinContent(i, h_eT_hcal->GetBinContent(i) / bin_width);
        //h_eT_hcal->SetBinError(i, h_eT_hcal->GetBinError(i) / bin_width);
        //h_eT_eta_hcal->SetBinContent(i, h_eT_eta_hcal->GetBinContent(i) / bin_width);
        //h_eT_eta_hcal->SetBinError(i, h_eT_eta_hcal->GetBinError(i) / bin_width);
        h_eT_eta_hcal_profile_hist->SetBinContent(i, h_eT_eta_hcal_profile->GetBinContent(i) / bin_width);
        h_eT_eta_hcal_profile_hist->SetBinError(i, h_eT_eta_hcal_profile->GetBinError(i) / bin_width);
        if (dataormc) {
        	hetdeta_hcalbin->SetBinContent(i, hetdeta_hcalbin->GetBinContent(i) / bin_width);
        	hetdeta_hcalbin->SetBinError(i, hetdeta_hcalbin->GetBinError(i) / bin_width);
        }
    }

	h_event_energy->Scale(1.0/totalweights);
	h_event_hcal_energy->Scale(1.0/totalweights);
	h_event_emcal_energy->Scale(1.0/totalweights);
	h_event_ihcal_energy->Scale(1.0/totalweights);
	h_event_ohcal_energy->Scale(1.0/totalweights);

	TH1F* h_emcal_correction;
	TH1F* h_ihcal_correction;
	TH1F* h_ohcal_correction;
	TH1F* h_calo_correction;
	TH1F* h_hcal_correction;

	if (dataormc) {
		h_emcal_correction = dynamic_cast<TH1F *>(h_eT_eta_emcal_profile_hist->Clone("h_emcal_correction"));
		h_emcal_correction->Divide(hetdeta_emcalbin);
		h_ihcal_correction = dynamic_cast<TH1F *>(h_eT_eta_ihcal_profile_hist->Clone("h_ihcal_correction"));
		h_ihcal_correction->Divide(hetdeta_ihcalbin);
		h_ohcal_correction = dynamic_cast<TH1F *>(h_eT_eta_ohcal_profile_hist->Clone("h_ohcal_correction"));
		h_ohcal_correction->Divide(hetdeta_ohcalbin);
		h_calo_correction = dynamic_cast<TH1F *>(h_eT_eta_calo_profile_hist->Clone("h_calo_correction"));
		h_calo_correction->Divide(hetdeta_calobin);
		h_hcal_correction = dynamic_cast<TH1F *>(h_eT_eta_hcal_profile_hist->Clone("h_hcal_correction"));
		h_hcal_correction->Divide(hetdeta_hcalbin);
	}
	
	out->Write();
	out->Close();


}
