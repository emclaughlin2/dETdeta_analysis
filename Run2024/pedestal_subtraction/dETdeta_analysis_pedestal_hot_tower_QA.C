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

const int cemcSize = 24576;
const int ihcalSize = 1536;
const int ohcalSize = 1536;
const int mbdSize = 256;
const int g4Size = 100000;
const int vtxSize = 3;

void dETdeta_analysis_pedestal_hot_tower_QA(int segment) {

	string outfilename = "/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/pedestal_subtraction/hot_tower_output/hot_tower_QA_run_54256_segment_" + to_string(segment) + ".root";
	TFile *out = new TFile(outfilename.c_str(),"RECREATE");

	TH1I* h_emcal = new TH1I("h_emcal","",cemcSize,0,cemcSize);
	TH1I* h_ihcal = new TH1I("h_ihcal","",ihcalSize,0,ihcalSize);
	TH1I* h_ohcal = new TH1I("h_ohcal","",ohcalSize,0,ohcalSize);

	TH2F* h_emcal_2D = new TH2F("h_emcal_2D","",96,0,96,256,0,256);
	TH2F* h_ihcal_2D = new TH2F("h_ihcal_2D","",24,0,24,64,0,64);
	TH2F* h_ohcal_2D = new TH2F("h_ohcal_2D","",24,0,24,64,0,64);
	  
    TChain chain("ttree");

    string infilename = "/sphenix/tg/tg01/commissioning/CaloCalibWG/egm2153/detdeta_run24auau/events_54256_pedestal_unc_" + to_string(segment) + ".root";
	chain.Add(infilename.c_str());

    int m_simtwrmult_cemc;
    float m_simtwr_cemc_e[cemcSize];
    int m_simtwr_cemc_ieta[cemcSize];
    int m_simtwr_cemc_iphi[cemcSize];
    float m_simtwr_cemc_chi2[cemcSize];
    int m_simtwrmult_ihcal;
    float m_simtwr_ihcal_e[ihcalSize];
    int m_simtwr_ihcal_ieta[ihcalSize];
    int m_simtwr_ihcal_iphi[ihcalSize];
    float m_simtwr_ihcal_chi2[ihcalSize];
    int m_simtwrmult_ohcal;
    float m_simtwr_ohcal_e[ohcalSize];
    int m_simtwr_ohcal_ieta[ohcalSize];
    int m_simtwr_ohcal_iphi[ohcalSize];
    float m_simtwr_ohcal_chi2[ohcalSize];

     // Set branch addresses
    int use_emcal = 1;
    int use_hcal = 1;
    if (use_emcal) {
	    chain.SetBranchAddress("sectorem", &m_simtwrmult_cemc);
	    chain.SetBranchAddress("emcalen", m_simtwr_cemc_e);
	    chain.SetBranchAddress("emcaletabin", m_simtwr_cemc_ieta);
	    chain.SetBranchAddress("emcalphibin", m_simtwr_cemc_iphi);
		chain.SetBranchAddress("emchi2", m_simtwr_cemc_chi2);
	}
	if (use_hcal) {
	    chain.SetBranchAddress("sectorih", &m_simtwrmult_ihcal);
	    chain.SetBranchAddress("ihcalen", m_simtwr_ihcal_e);
	    chain.SetBranchAddress("ihcaletabin", m_simtwr_ihcal_ieta);
	    chain.SetBranchAddress("ihcalphibin", m_simtwr_ihcal_iphi);
	    chain.SetBranchAddress("ihchi2", m_simtwr_ihcal_chi2);

	    chain.SetBranchAddress("sectoroh", &m_simtwrmult_ohcal);
	    chain.SetBranchAddress("ohcalen", m_simtwr_ohcal_e);
	    chain.SetBranchAddress("ohcaletabin", m_simtwr_ohcal_ieta);
	    chain.SetBranchAddress("ohcalphibin", m_simtwr_ohcal_iphi);
	    chain.SetBranchAddress("ohchi2", m_simtwr_ohcal_chi2);
	}

	int eventnumber = 0;
	float totalweights = 0.0;
	bool found_hot_tower = false;

    Long64_t nEntries = chain.GetEntries();
    std::cout << nEntries << std::endl;
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        chain.GetEntry(entry);
    	if (eventnumber % 1000 == 0) cout << "event " << eventnumber << endl;

  		eventnumber++;
		if (use_emcal) {
			for (int i = 0; i < m_simtwrmult_cemc; i++) {
				if (m_simtwr_cemc_e[i] < -20) {
					found_hot_tower = true;
					h_emcal->Fill(i);
					h_emcal_2D->Fill(m_simtwr_cemc_ieta[i],m_simtwr_cemc_iphi[i]);
				}
			}
		}
		if (use_hcal) {
			for (int i = 0; i < m_simtwrmult_ihcal; i++) {
				if (m_simtwr_ihcal_e[i] < -10) {
					found_hot_tower = true;
					h_ihcal->Fill(i);
					h_ihcal_2D->Fill(m_simtwr_ihcal_ieta[i],m_simtwr_ihcal_iphi[i]);
				}
			}

			for (int i = 0; i < m_simtwrmult_ohcal; i++) {
		    	if (m_simtwr_ohcal_e[i] < -10) {
		    		found_hot_tower = true;
					h_ohcal->Fill(i);
					h_ohcal_2D->Fill(m_simtwr_ohcal_ieta[i],m_simtwr_ohcal_iphi[i]);
				}
			}
		}
		
	}
	if (found_hot_tower) {
		out->Write();
	}
	out->Close();

}
