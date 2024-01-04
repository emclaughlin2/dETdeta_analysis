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
#include <string.h>

using namespace std;

// hijing centrality binning
const float centrality_binning[11] = {0.0,3.0,2155.0,5270.0,11058.0,20421.0,34539.0,54158.0,80942.0,117224.0,250000.0};

void dETdeta_vertex_mc(const char* generator = "hijing") {	

	string outfile = "hjing_vertex_centrality.root";
	if (!strcmp(generator, "epos")) {
		outfile = "epos_vertex_centrality.root";
	} else if (!strcmp(generator, "ampt")) {
		outfile = "ampt_vertex_centrality.root";
	}

	TDatabasePDG *_pdg = new TDatabasePDG();

	TFile *out = new TFile(outfile.c_str(),"RECREATE");
	TH1F* h_vz = new TH1F("h_vz","",400, -100, 100);

	TH2F* h_2D_ihcal_calib = new TH2F("h_2D_ihcal_calib","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_ohcal_calib = new TH2F("h_2D_ohcal_calib","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_emcal_calib = new TH2F("h_2D_emcal_calib","",96,0.,96.,256,0.,256.);

  	TH2F* h_2D_ihcal_hits = new TH2F("h_2D_ihcal_hits","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_ohcal_hits = new TH2F("h_2D_ohcal_hits","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_emcal_hits = new TH2F("h_2D_emcal_hits","",96,0.,96.,256,0.,256.);
	
	TH1F* h_event_energy = new TH1F("h_event_energy","", 5000,0,5000);
	TH1F* h_event_hcal_energy = new TH1F("h_event_hcal_energy","", 5000,0,5000);
	TH1F* h_event_emcal_energy = new TH1F("h_event_emcal_energy","", 5000,0,5000);
	TH1F* h_event_ihcal_energy = new TH1F("h_event_ihcal_energy","", 1000,0,1000);
	TH1F* h_event_ohcal_energy = new TH1F("h_event_ohcal_energy","", 1000,0,1000);

	TH1F* h_emcal = new TH1F("h_emcal","",1000,0,10);
	TH1F* h_ihcal = new TH1F("h_ihcal","",1000,0,10);
	TH1F* h_ohcal = new TH1F("h_ohcal","",1000,0,10);
    
    TChain chain("ttree");

    if (!strcmp(generator, "epos")) {
    	// location of EPOS files 
	    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/eposrun/condor/";
	    for (int i = 0; i < 900; i++) {
	    	TString wildcardPath = TString::Format("%sOutDir%d/epos_sim_output.root", inputDirectory, i);
	    	chain.Add(wildcardPath);
	    }
    } else if (!strcmp(generator, "hijing")) {
    	// location of HIJING files 
	    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
	    for (int i = 0; i < 555; i++) {
	    	TString wildcardPath = TString::Format("%sevents_20231122_nopileup_mc_cor_%d.root", inputDirectory, i);
	    	chain.Add(wildcardPath);
	    }
    } else if (!strcmp(generator, "ampt")) {
    	// location of AMPT files 
	    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/amptrun/condor/";
	    for (int i = 0; i < 1000; i++) {
	    	TString wildcardPath = TString::Format("%sOutDir%d/ampt_sim_output.root", inputDirectory , i);
	    	chain.Add(wildcardPath);
	    }
    } else {
    	std::cout << "generator not found" << std::endl;
    	return;
    }

	TTreeReader reader(&chain);
	
	TTreeReaderValue<int> m_simtwrmult_cemc(reader, "sectorem");
	TTreeReaderArray<float> m_simtwr_cemc_e(reader, "emcalen");
	TTreeReaderArray<int> m_simtwr_cemc_ieta(reader, "emcaletabin");
	TTreeReaderArray<int> m_simtwr_cemc_iphi(reader, "emcalphibin");
	TTreeReaderArray<float> m_simtwr_cemc_eta(reader, "emetacor");

	TTreeReaderValue<int> m_simtwrmult_ihcal(reader, "sectorih");
	TTreeReaderArray<float> m_simtwr_ihcal_e(reader, "ihcalen");
	TTreeReaderArray<int> m_simtwr_ihcal_ieta(reader, "ihcaletabin");
	TTreeReaderArray<int> m_simtwr_ihcal_iphi(reader, "ihcalphibin");
	TTreeReaderArray<float> m_simtwr_ihcal_eta(reader, "ihetacor");

	TTreeReaderValue<int> m_simtwrmult_ohcal(reader, "sectoroh");
	TTreeReaderArray<float> m_simtwr_ohcal_e(reader, "ohcalen");
	TTreeReaderArray<int> m_simtwr_ohcal_ieta(reader, "ohcaletabin");
	TTreeReaderArray<int> m_simtwr_ohcal_iphi(reader, "ohcalphibin");
	TTreeReaderArray<float> m_simtwr_ohcal_eta(reader, "ohetacor");

	TTreeReaderValue<int> m_sectormb(reader, "sectormb");
	TTreeReaderArray<float> m_mbenergy(reader, "mbenrgy");

	TTreeReaderArray<float> m_vtx(reader, "track_vtx");
	TTreeReaderArray<float> mbd_vtx(reader, "mbd_vtx");

	TH1F* h_mbd = new TH1F("h_mbd","",250000,0,250000);

	int eventnumber = 0;
	float delta_eta = 0.09167;
	float delta_em_eta = 0.022918;

    while (reader.Next()) {
    	if (eventnumber % 1000 == 0) cout << "event " << eventnumber << endl;

		float emcale = 0;
		float ihcale = 0;
		float ohcale = 0;
		float totale = 0;
		float truthe = 0;

		float E_emcal[96] = {0};
		float E_ihcal[24] = {0};
  		float E_ohcal[24] = {0};

  		float mbde = 0;

  		// require that simulation could reconstruct a vertex for the event
  		eventnumber++;
  		h_vz->Fill(m_vtx[2]);

  		for (int i = 0; i < *m_sectormb; i++) {
  			mbde += m_mbenergy[i];
  		}
  		h_mbd->Fill(mbde);
		
		for (int i = 0; i < *m_simtwrmult_cemc; i++) {
			h_2D_emcal_calib->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]);
			if (m_simtwr_cemc_e[i] > 0.2) h_2D_emcal_hits->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i]);
			emcale += m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]); 
			h_emcal->Fill(m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]));
		}

		for (int i = 0; i < *m_simtwrmult_ihcal; i++) {
			h_2D_ihcal_calib->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]);
			if (m_simtwr_ihcal_e[i] > 0.05) h_2D_ihcal_hits->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i]);
			ihcale += m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]); 
			h_ihcal->Fill(m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]));
		}
		for (int i = 0; i < *m_simtwrmult_ohcal; i++) {
			h_2D_ohcal_calib->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]);
			if (m_simtwr_ohcal_e[i] > 0.1) h_2D_ohcal_hits->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i]);
			ohcale += m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]); 
			h_ohcal->Fill(m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]));
		}
		
		totale = ihcale + ohcale;
		h_event_energy->Fill(totale + emcale);
		h_event_hcal_energy->Fill(totale);
		h_event_emcal_energy->Fill(emcale);
		h_event_ihcal_energy->Fill(ihcale);
		h_event_ohcal_energy->Fill(ohcale);
	}
	
	out->Write();
	out->Close();
}
