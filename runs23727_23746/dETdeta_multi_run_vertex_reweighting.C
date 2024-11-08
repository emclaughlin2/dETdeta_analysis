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

int good_runs[8] = {23727, 23735, 23737, 23738, 23739, 23740, 23743, 23745};
//int good_runs[5] = {23727, 23735, 23737, 23743, 23745};

void dETdeta_multi_run_vertex_reweighting(const char* generator) {	

	TFile *out = new TFile(TString::Format("dETdeta_vertex_reweight_run23727_23745_%s_p015.root", generator),"RECREATE");
	
	TH1F* h_vz_data = new TH1F("h_vz_data","",200,-50,50);
	TH1F* h_vz_data_zoom = new TH1F("h_vz_data_zoom","",200,-50,50);
	TH1F* h_vz_mc = new TH1F("h_vz_mc","",200,-50,50);

	TH1F* h_mbd_data = new TH1F("h_mbd_data","",12000,0,6000);
	TH1F* h_mbd_mc = new TH1F("h_mbd_mc","",250000,0,250000);

	TChain mcchain("ttree");
	TChain datachain("ttree");
    

    const char* dataInputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
    for (int r = 0; r < 8; r++) { 
    	TString dataWildcardPath = TString::Format("%sevents_p015_%d_data_cor_*.root", dataInputDirectory, good_runs[r]);
    	//TString dataWildcardPath = TString::Format("%sevents_p011_zs_%d_data_cor_*.root", dataInputDirectory, good_runs[r]);
    	datachain.Add(dataWildcardPath);
    }
	TTreeReader datareader(&datachain);

	if (!strcmp(generator, "hijing")) {
		const char* mcInputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
		for (int i = 0; i < 2500; i++) {
		    	TString mcWildcardPath = TString::Format("%sevents_20240110_run101_nopileup_mc_cor_%d.root", mcInputDirectory, i);
		    	mcchain.Add(mcWildcardPath);
		}
	} else if (!strcmp(generator, "reweight_hijing")) {
	    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
	    for (int i = 0; i < 2500; i++) { 
    		//TString wildcardPath = TString::Format("%sevents_hijing_run101_reweighted_mc_cor_%d.root", inputDirectory, i);
	    	TString wildcardPath = TString::Format("%sevents_20240316_run101_nopileup_mc_cor_%d.root", inputDirectory, i); //edited 3.15.24
	    	mcchain.Add(wildcardPath);
	    }
	} else if (!strcmp(generator, "epos")) {
		const char* mcInputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/eposrun/plugdoor_condor/";
	    for (int i = 0; i < 1000; i++) {
	    	TString mcWildcardPath = TString::Format("%sOutDir%d/epos_sim_output.root", mcInputDirectory, i);
	    	mcchain.Add(mcWildcardPath);
	    }
	} else if (!strcmp(generator, "ampt")) {
		const char* mcInputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/amptrun/plugdoor_condor/";
	    for (int i = 0; i < 1000; i++) {
	    	TString mcWildcardPath = TString::Format("%sOutDir%d/ampt_sim_output.root", mcInputDirectory, i);
	    	mcchain.Add(mcWildcardPath);
	    }
	} else if (!strcmp(generator, "reweight_ampt")) {
	    const char* inputDirectory = "/sphenix/tg/tg01/commissioning/CaloCalibWG/egm2153/detdeta_datasets/"; 
		for (int i = 0; i < 3000; i++) { 
    		TString mcWildcardPath = TString::Format("%sevents_ampt_reweighted_run14_ampt_cor_%d.root", inputDirectory, i);
	    	mcchain.Add(mcWildcardPath);
	    }
	} else if (!strcmp(generator, "reweight_epos")) {
	    const char* inputDirectory = "/sphenix/tg/tg01/commissioning/CaloCalibWG/egm2153/detdeta_datasets/"; 
		for (int i = 0; i < 3000; i++) { 
    		TString mcWildcardPath = TString::Format("%sevents_epos_reweighted_run14_epos_cor_%d.root", inputDirectory, i);
	    	mcchain.Add(mcWildcardPath);
	    }
	}
	TTreeReader mcreader(&mcchain);


	TTreeReaderArray<float> m_vtx(datareader, "track_vtx");
	TTreeReaderValue<int> m_sectormb(datareader, "sectormb");
	TTreeReaderArray<float> m_mbenergy(datareader, "mbenrgy");
	TTreeReaderValue<bool> m_isMinBias(datareader, "isMinBias");

	TTreeReaderArray<float> track_vtx(mcreader, "track_vtx");
	TTreeReaderValue<int> sectormb(mcreader, "sectormb");
	TTreeReaderArray<float> mbenergy(mcreader, "mbenrgy");

	// get data vertex distribution and mbd charge distribution
	int eventnumber = 0;
    while (datareader.Next()) {
    	if (eventnumber % 10000 == 0) cout << "data event " << eventnumber << endl;
  		float mbde = 0;
  		eventnumber++;
  		if (!*m_isMinBias) { continue; }
  		if (isnan(m_vtx[2])) { continue; }
  		h_vz_data->Fill(m_vtx[2]);
  		if (fabs(m_vtx[2]) < 20.0) h_vz_data_zoom->Fill(m_vtx[2]);

  		for (int i = 0; i < *m_sectormb; i++) {
  			mbde += m_mbenergy[i];
  		}
  		h_mbd_data->Fill(mbde);	
	}
	eventnumber = 0;
	// get mc vertex distribution and mbd charge distribution
	while (mcreader.Next()) {
    	if (eventnumber % 10000 == 0) cout << "mc event " << eventnumber << endl;
  		float mbde = 0;
  		eventnumber++;
  		h_vz_mc->Fill(track_vtx[2]);

  		for (int i = 0; i < *sectormb; i++) {
  			mbde += mbenergy[i];
  		}
  		h_mbd_mc->Fill(mbde);	
	}

	h_vz_data->Scale(1.0/h_vz_data->GetEntries());
	h_vz_mc->Scale(1.0/h_vz_mc->GetEntries());

	TH1F* vz_ratio = static_cast<TH1F*>(h_vz_data->Clone("vz_ratio"));
	vz_ratio->Divide(h_vz_mc);

	TTree* T = new TTree("T", "");
	float mc_cent[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	float data_cent[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	float vertex_reweight[200] = {0.};
    T->Branch("data_centrality", data_cent, "data_cent[20]/F");
    T->Branch("mc_centrality", mc_cent, "mc_cent[20]/F");
    T->Branch("vertex_reweight", vertex_reweight, "vertex_reweight[200]/F");

	int mc_nq = 21;
	double mc_q[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	double mc_xq[] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
	h_mbd_mc->GetQuantiles(mc_nq,mc_q,mc_xq);
	std::cout << "mc centrality" << std::endl;
	for (int i = 0; i < 21; i++) {
		std::cout << mc_q[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "mc cent array" << std::endl;
	std::cout << mc_cent[0] << " ";
	for (int i = 1; i < 20; i++) {
		mc_cent[i] = float(mc_q[i+1]);
		std::cout << mc_cent[i] << " " << std::endl;
	}
	std::cout << std::endl;

	int data_nq = 19;
	double data_q[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	double data_xq[] = {0.0,0.0555,0.11111,0.16666,0.22222,0.27777,0.33333,0.38888,0.44444,0.5,0.55555,0.61111,0.66666,0.72222,0.77777,0.83333,0.88888,0.94444,1.0};
	h_mbd_data->GetQuantiles(data_nq,data_q,data_xq);
	std::cout << "data centrality" << std::endl;
	for (int i = 0; i < 19; i++) {
		std::cout << data_q[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "data cent array" << std::endl;
	std::cout << data_cent[0] << " ";
	for (int i = 1; i < 20; i++) {
		data_cent[i] = float(data_q[i-1]);
		std::cout << data_cent[i] << " " << std::endl;
	}
	std::cout << std::endl;

	// for data just add a zero to the front of the array (index 0)
	// for mc get rid of the value at array index 1

	for (int i = 1; i < vz_ratio->GetNbinsX() + 1; i++) {
		vertex_reweight[i-1] = vz_ratio->GetBinContent(i);
	}

	T->Fill();
	T->Write();
	h_vz_data->Write();
	h_vz_data_zoom->Write();
	h_vz_mc->Write();
	vz_ratio->Write();
	h_mbd_data->Write();
	h_mbd_mc->Write();
	
	out->Close();
}
