#include <iostream>
#include <TH2F.h>
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

void dETdeta_centrality_npart_calc(int runnumber, const char* generator) {	

	TFile *out = new TFile(TString::Format("dETdeta_centrality_npart_run%d_%s_p011.root",runnumber, generator),"RECREATE");

	TH1F* h_mbd_data = new TH1F("h_mbd_data","",12000,0,6000);
	TH1F* h_mbd_mc = new TH1F("h_mbd_mc","",12000,0,6000);
	TH1F* h_npart = new TH1F("h_npart","",500,0,500);
	TH2F* npart_vs_mbd = new TH2F("npart_vs_mbd","",500,0,500,12000,0,6000);

	TChain mcchain("ttree");
	TChain datachain("ttree");
    
    const char* dataInputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
    for (int i = 0; i < 230; i++) {  
    	TString dataWildcardPath = TString::Format("%sevents_p011_zs_%d_data_cor_%d.root", dataInputDirectory, runnumber, i);
    	datachain.Add(dataWildcardPath);
    }
	TTreeReader datareader(&datachain);

	if (!strcmp(generator, "hijing")) {
		const char* mcInputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
		for (int i = 0; i < 2500; i++) {
		    	TString mcWildcardPath = TString::Format("%sevents_hijing_run101_reweighted_mc_cor_%d.root", mcInputDirectory, i);
		    	mcchain.Add(mcWildcardPath);
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
	}
	TTreeReader mcreader(&mcchain);

	TTreeReaderArray<float> m_vtx(datareader, "track_vtx");
	TTreeReaderValue<int> m_sectormb(datareader, "sectormb");
	TTreeReaderArray<float> m_mbenergy(datareader, "mbenrgy");
	TTreeReaderValue<bool> m_isMinBias(datareader, "isMinBias");

	TTreeReaderArray<float> track_vtx(mcreader, "track_vtx");
	TTreeReaderValue<int> sectormb(mcreader, "sectormb");
	TTreeReaderArray<float> mbenergy(mcreader, "mbenrgy");
	TTreeReaderValue<int> npart(mcreader, "npart");

	// get data vertex distribution and mbd charge distribution
	int eventnumber = 0;
    while (datareader.Next()) {
    	if (eventnumber % 1000 == 0) cout << "data event " << eventnumber << endl;
  		eventnumber++;
  		if (!*m_isMinBias) { continue; }
  		if (isnan(m_vtx[2])) { continue; }
        if (m_vtx[2] < -50 || m_vtx[2] > 50) { continue; }
  		float mbde = 0;
  		for (int i = 0; i < *m_sectormb; i++) {
  			mbde += m_mbenergy[i];
  		}
  		h_mbd_data->Fill(mbde);	
	}
	eventnumber = 0;
	// get mc vertex distribution and mbd charge distribution
	while (mcreader.Next()) {
    	if (eventnumber % 1000 == 0) cout << "mc event " << eventnumber << endl;
  		float mbde = 0;
  		eventnumber++;
  		h_npart->Fill(*npart);
  		for (int i = 0; i < *sectormb; i++) {
  			mbde += mbenergy[i];
  		}
  		h_mbd_mc->Fill(mbde);	
  		npart_vs_mbd->Fill(*npart, mbde);
	}

	h_mbd_data->Write();
	h_mbd_mc->Write();
	h_npart->Write();
	npart_vs_mbd->Write();
	
	out->Close();
}
