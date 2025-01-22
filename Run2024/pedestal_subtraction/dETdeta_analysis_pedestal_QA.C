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
#include <TProfile2D.h>

using namespace std;

std::set<std::tuple<int, int>> emcal_hot_dead_map;
std::set<std::tuple<int, int>> ihcal_hot_dead_map;
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

void fill_hot_dead_map() {

	// hot dead maps 
    vector<int> *emcal_hot_dead_ieta = nullptr;
    vector<int> *emcal_hot_dead_iphi = nullptr;
    vector<int> *ihcal_hot_dead_ieta = nullptr;
    vector<int> *ihcal_hot_dead_iphi = nullptr;
    vector<int> *ohcal_hot_dead_ieta = nullptr;
    vector<int> *ohcal_hot_dead_iphi = nullptr;
    vector<float> *ihcal_eta_bins = nullptr;
    vector<float> *ohcal_eta_bins = nullptr;

	TFile *hotdeadfile = new TFile("/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/run54912_hotdeadmap_z_-20_20_new_2024p007.root", "READ");
	TTree *hotdeadtree = dynamic_cast<TTree*>(hotdeadfile->Get("T"));
	hotdeadtree->SetBranchAddress("emcal_hot_dead_ieta", &emcal_hot_dead_ieta);
    hotdeadtree->SetBranchAddress("emcal_hot_dead_iphi", &emcal_hot_dead_iphi);
    hotdeadtree->SetBranchAddress("ihcal_hot_dead_ieta", &ihcal_hot_dead_ieta);
    hotdeadtree->SetBranchAddress("ihcal_hot_dead_iphi", &ihcal_hot_dead_iphi);
    hotdeadtree->SetBranchAddress("ohcal_hot_dead_ieta", &ohcal_hot_dead_ieta);
    hotdeadtree->SetBranchAddress("ohcal_hot_dead_iphi", &ohcal_hot_dead_iphi);

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
	hotdeadfile->Close();

}

void dETdeta_analysis_pedestal_QA(int runnumber = 23727, const char* opt_tag = "") {

	string filename = "dETdeta_analysis_QA";
	filename += "_" + to_string(runnumber);
	string opttag = opt_tag;
	if (strcmp(opt_tag,"")) { filename += "_" + opttag; }
  	filename += ".root";

  	fill_hot_dead_map();

	TFile *out = new TFile(filename.c_str(),"RECREATE");
	
	TH2F* h_2D_ihcal_calib = new TH2F("h_2D_ihcal_calib","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_ohcal_calib = new TH2F("h_2D_ohcal_calib","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_emcal_calib = new TH2F("h_2D_emcal_calib","",96,0.,96.,256,0.,256.);

  	TH2F* h_2D_ihcal_calibT = new TH2F("h_2D_ihcal_calibT","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_ohcal_calibT = new TH2F("h_2D_ohcal_calibT","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_emcal_calibT = new TH2F("h_2D_emcal_calibT","",96,0.,96.,256,0.,256.);

  	TProfile2D* h_2D_emcal_chi2 = new TProfile2D("h_2D_emcal_chi2","",96,0.,96.,256,0.,256.);
  	TProfile2D* h_2D_ihcal_chi2 = new TProfile2D("h_2D_ihcal_chi2","",24,0.,24.,64,0.,64.);
  	TProfile2D* h_2D_ohcal_chi2 = new TProfile2D("h_2D_ohcal_chi2","",24,0.,24.,64,0.,64.);
	  
    TChain chain("ttree");

    const char* inputDirectory = "/sphenix/tg/tg01/commissioning/CaloCalibWG/egm2153/detdeta_run24auau/";
	for (int i = 0; i < 1; i++) {
		TString wildcardPath = TString::Format("%sevents_%d_pedestal_cor_%d.root", inputDirectory, runnumber, i);
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
    float m_simtwr_cemc_chi2[cemcSize];
    int m_simtwrmult_ihcal;
    float m_simtwr_ihcal_e[ihcalSize];
    float m_simtwr_ihcal_zs_e[ihcalSize];
    int m_simtwr_ihcal_ieta[ihcalSize];
    int m_simtwr_ihcal_iphi[ihcalSize];
    float m_simtwr_ihcal_eta[ihcalSize];
    int m_simtwr_ihcal_adc[ihcalSize];
    int m_simtwr_ihcal_zs_adc[ihcalSize];
    float m_simtwr_ihcal_time[ihcalSize];
    float m_simtwr_ihcal_chi2[ihcalSize];
    int m_simtwrmult_ohcal;
    float m_simtwr_ohcal_e[ohcalSize];
    float m_simtwr_ohcal_zs_e[ohcalSize];
    int m_simtwr_ohcal_ieta[ohcalSize];
    int m_simtwr_ohcal_iphi[ohcalSize];
    float m_simtwr_ohcal_eta[ohcalSize];
    int m_simtwr_ohcal_adc[ohcalSize];
    int m_simtwr_ohcal_zs_adc[ohcalSize];
    float m_simtwr_ohcal_time[ohcalSize];
    float m_simtwr_ohcal_chi2[ohcalSize];

     // Set branch addresses
    int use_emcal = 1;
    int use_hcal = 1;
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
		chain.SetBranchAddress("emchi2", m_simtwr_cemc_chi2);
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
	    chain.SetBranchAddress("ihchi2", m_simtwr_ihcal_chi2);

	    chain.SetBranchAddress("sectoroh", &m_simtwrmult_ohcal);
	    chain.SetBranchAddress("ohcalen", m_simtwr_ohcal_e);
	    chain.SetBranchAddress("ohcaletabin", m_simtwr_ohcal_ieta);
	    chain.SetBranchAddress("ohcalphibin", m_simtwr_ohcal_iphi);
	    chain.SetBranchAddress("ohetacor", m_simtwr_ohcal_eta);
		chain.SetBranchAddress("ohcalzs", m_simtwr_ohcal_zs_e);
	    chain.SetBranchAddress("ohcaladc", m_simtwr_ohcal_adc);
	    chain.SetBranchAddress("ohcalzsadc", m_simtwr_ohcal_zs_adc);
	    chain.SetBranchAddress("ohcalt", m_simtwr_ohcal_time);
	    chain.SetBranchAddress("ohchi2", m_simtwr_ohcal_chi2);
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
				std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i]);
			    auto it = emcal_hot_dead_map.find(hot_tower);
			    if (it != emcal_hot_dead_map.end()) { continue; }
				h_2D_emcal_calib->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight);
				h_2D_emcal_calibT->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
				h_2D_emcal_chi2->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_chi2[i]);
				emcale += m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]); 
			}
		}
		if (use_hcal) {
			for (int i = 0; i < m_simtwrmult_ihcal; i++) {
				std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i]);
			    auto it = ihcal_hot_dead_map.find(hot_tower);
			    if (it != ihcal_hot_dead_map.end()) { continue; }
				h_2D_ihcal_calib->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]*vz_weight);
				h_2D_ihcal_calibT->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
				h_2D_ihcal_chi2->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_chi2[i]);
			}

			for (int i = 0; i < m_simtwrmult_ohcal; i++) {
				std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i]);
			    auto it = ohcal_hot_dead_map.find(hot_tower);
			    if (it != ohcal_hot_dead_map.end()) { continue; }
		    	h_2D_ohcal_calib->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]*vz_weight);
				h_2D_ohcal_calibT->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
				h_2D_ohcal_chi2->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_chi2[i]);
			}
		}
		
	}

	if (use_emcal) {
		h_2D_emcal_calib->Scale(1.0/totalweights);
		h_2D_emcal_calibT->Scale(1.0/totalweights);
	}
	if (use_hcal) {
		h_2D_ihcal_calib->Scale(1.0/totalweights);
		h_2D_ihcal_calibT->Scale(1.0/totalweights);
		h_2D_ohcal_calib->Scale(1.0/totalweights);
		h_2D_ohcal_calibT->Scale(1.0/totalweights);
	}
	if (use_hcal) {
		for (int i = 0; i < 24; i++) {
			for (int j = 0; j < 64; j++) {
				if (h_2D_ihcal_chi2->GetBinContent(h_2D_ihcal_chi2->FindBin(i,j)) > 10) { 
					std::cout << "IHCal ieta " << i << " iphi " << j << " mean chi2 " << h_2D_ihcal_chi2->GetBinContent(h_2D_ihcal_chi2->FindBin(i,j)) << std::endl;
				}
				if (h_2D_ohcal_chi2->GetBinContent(h_2D_ohcal_chi2->FindBin(i,j)) > 10) { 
					std::cout << "OHCal ieta " << i << " iphi " << j << " mean chi2 " << h_2D_ohcal_chi2->GetBinContent(h_2D_ohcal_chi2->FindBin(i,j)) << std::endl;
				}
			}
		}
	}
	if (use_emcal) {
		for (int i = 0; i < 96; i++) {
			for (int j = 0; j < 256; j++) {
				if (h_2D_emcal_chi2->GetBinContent(h_2D_emcal_chi2->FindBin(i,j)) > 10) { 
					std::cout << "EMCal ieta " << i << " iphi " << j << " mean chi2 " << h_2D_emcal_chi2->GetBinContent(h_2D_emcal_chi2->FindBin(i,j)) << std::endl;
				}	
			}
		}
	}
	
	out->Write();
	out->Close();


}
