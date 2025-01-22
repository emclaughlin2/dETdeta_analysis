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

void dETdeta_em_had_energy_fraction(const char* generator) {	

	TFile *out = new TFile(TString::Format("dETdeta_em_had_energy_fraction_%s.root", generator),"RECREATE");
	
	TH1F* h_had_frac_emcal = new TH1F("h_had_frac_emcal","",100,0,1);
	TH1F* h_had_frac_ihcal = new TH1F("h_had_frac_ihcal","",100,0,1);
	TH1F* h_had_frac_ohcal = new TH1F("h_had_frac_ohcal","",100,0,1);
	TH1F* h_had_frac_calo = new TH1F("h_had_frac_calo","",100,0,1);
	TH1F* h_had_frac_hcal = new TH1F("h_had_frac_hcal","",100,0,1);

	TH1F* h_had_energy_emcal = new TH1F("h_had_energy_emcal","",1000,0,1000);
	TH1F* h_had_energy_ihcal = new TH1F("h_had_energy_ihcal","",1000,0,1000);
	TH1F* h_had_energy_ohcal = new TH1F("h_had_energy_ohcal","",1000,0,1000);
	TH1F* h_had_energy_calo = new TH1F("h_had_energy_calo","",1000,0,1000);
	TH1F* h_had_energy_hcal = new TH1F("h_had_energy_hcal","",1000,0,1000);

	TH1F* h_em_energy_emcal = new TH1F("h_em_energy_emcal","",1000,0,1000);
	TH1F* h_em_energy_ihcal = new TH1F("h_em_energy_ihcal","",1000,0,1000);
	TH1F* h_em_energy_ohcal = new TH1F("h_em_energy_ohcal","",1000,0,1000);
	TH1F* h_em_energy_calo = new TH1F("h_em_energy_calo","",1000,0,1000);
	TH1F* h_em_energy_hcal = new TH1F("h_em_energy_hcal","",1000,0,1000);

	TChain emcalchain("EnergyTree");
	TChain ihcalchain("EnergyTree");
	TChain ohcalchain("EnergyTree");

	const char* mcInputDirectory = "/sphenix/tg/tg01/commissioning/CaloCalibWG/egm2153/detdeta_run24auau/";
	if (!strcmp(generator, "reweight_hijing")) {
		for (int i = 0; i < 5000; i++) {
		    TString emcalWildcardPath = TString::Format("%sEMCal_energy_hijing_reweighted_run14_%d_rw_hijing.root", mcInputDirectory, i); 
		    emcalchain.Add(emcalWildcardPath);
		    TString ihcalWildcardPath = TString::Format("%sIHCal_energy_hijing_reweighted_run14_%d_rw_hijing.root", mcInputDirectory, i); 
		    ihcalchain.Add(ihcalWildcardPath);
		    TString ohcalWildcardPath = TString::Format("%sOHCal_energy_hijing_reweighted_run14_%d_rw_hijing.root", mcInputDirectory, i); 
		    ohcalchain.Add(ohcalWildcardPath);
		}
	} else if (!strcmp(generator, "reweight_ampt")) {
		for (int i = 0; i < 5000; i++) {
		    TString emcalWildcardPath = TString::Format("%sEMCal_energy_ampt_reweighted_run14_%d_rw_ampt.root", mcInputDirectory, i); 
		    emcalchain.Add(emcalWildcardPath);
		    TString ihcalWildcardPath = TString::Format("%sIHCal_energy_ampt_reweighted_run14_%d_rw_ampt.root", mcInputDirectory, i); 
		    ihcalchain.Add(ihcalWildcardPath);
		    TString ohcalWildcardPath = TString::Format("%sOHCal_energy_ampt_reweighted_run14_%d_rw_ampt.root", mcInputDirectory, i); 
		    ohcalchain.Add(ohcalWildcardPath);
		}
	} else if (!strcmp(generator, "reweight_epos")) {
		for (int i = 0; i < 5000; i++) {
		    TString emcalWildcardPath = TString::Format("%sEMCal_energy_epos_reweighted_run14_%d_rw_epos.root", mcInputDirectory, i); 
		    emcalchain.Add(emcalWildcardPath);
		    TString ihcalWildcardPath = TString::Format("%sIHCal_energy_epos_reweighted_run14_%d_rw_epos.root", mcInputDirectory, i); 
		    ihcalchain.Add(ihcalWildcardPath);
		    TString ohcalWildcardPath = TString::Format("%sOHCal_energy_epos_reweighted_run14_%d_rw_epos.root", mcInputDirectory, i); 
		    ohcalchain.Add(ohcalWildcardPath);
		}
	}

	float em_emcal, had_emcal;
	float em_ihcal, had_ihcal;
	float em_ohcal, had_ohcal;
	emcalchain.SetBranchAddress("em_energy", &em_emcal);
	emcalchain.SetBranchAddress("had_energy", &had_emcal);
	ihcalchain.SetBranchAddress("em_energy", &em_ihcal);
	ihcalchain.SetBranchAddress("had_energy", &had_ihcal);
	ohcalchain.SetBranchAddress("em_energy", &em_ohcal);
	ohcalchain.SetBranchAddress("had_energy", &had_ohcal);

	const Long64_t nEntriesEmcal = emcalchain.GetEntries();
	const Long64_t nEntriesIhcal = ihcalchain.GetEntries();
	const Long64_t nEntriesOhcal = ohcalchain.GetEntries();

	std::cout << "emcal chain entries " << nEntriesEmcal << std::endl;
	std::cout << "ihcal chain entries " << nEntriesIhcal << std::endl;
	std::cout << "ohcal chain entries " << nEntriesOhcal << std::endl;

	float emcal_sf = 2e-02;
	float ihcal_sf = 0.162166;
	float ohcal_sf = 3.38021e-02;

	if (nEntriesEmcal != nEntriesIhcal || nEntriesEmcal != nEntriesOhcal) {
	    std::cout << "Error: TTrees have different numbers of entries!" << std::endl;
	    return; 
	}

	const Long64_t nEntries = nEntriesEmcal; 
	for (Long64_t entry = 0; entry < nEntries; ++entry) {
	    emcalchain.GetEntry(entry);
	    ihcalchain.GetEntry(entry);
	    ohcalchain.GetEntry(entry);
	   	float em_hcal = em_ihcal/ihcal_sf + em_ohcal/ohcal_sf;
	    float had_hcal = had_ihcal/ihcal_sf + had_ohcal/ohcal_sf;
	    float em_calo = em_emcal/emcal_sf + em_ihcal/ihcal_sf + em_ohcal/ohcal_sf;
	    float had_calo = had_emcal/emcal_sf + had_ihcal/ihcal_sf + had_ohcal/ohcal_sf;

	    h_em_energy_emcal->Fill(em_emcal/emcal_sf);
	    h_had_energy_emcal->Fill(had_emcal/emcal_sf);
	    h_em_energy_ihcal->Fill(em_ihcal/ihcal_sf);
	    h_had_energy_ihcal->Fill(had_ihcal/ihcal_sf);
	    h_em_energy_ohcal->Fill(em_ohcal/ohcal_sf);
	    h_had_energy_ohcal->Fill(had_ohcal/ohcal_sf);
	    h_em_energy_hcal->Fill(em_hcal);
	    h_had_energy_hcal->Fill(had_hcal);
	    h_em_energy_calo->Fill(em_calo);
	    h_had_energy_calo->Fill(had_calo);

	    h_had_frac_emcal->Fill((had_emcal)/(em_emcal+had_emcal));
	    h_had_frac_ihcal->Fill((had_ihcal)/(em_ihcal+had_ihcal));
	    h_had_frac_ohcal->Fill((had_ohcal)/(em_ohcal+had_ohcal));
	    h_had_frac_hcal->Fill((had_hcal)/(em_hcal+had_hcal));
	    h_had_frac_calo->Fill((had_calo)/(em_calo+had_calo));
	}

	out->Write();
	out->Close();
}
