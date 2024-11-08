#include <iostream>
#include <TH2F.h>
#include <TH1F.h>
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

void dETdeta_multi_run_centrality_check() {	

	TFile *out = new TFile("dETdeta_centrality_check_run23727_23745_p015.root","RECREATE");
	TH1F* h_cent_ismb = new TH1F("h_cent_ismb","",100,0,100);
	TH1F* h_cent_ismb_notnan = new TH1F("h_cent_ismb_notnan","",100,0,100);
	TH1F* h_cent_ismb_vz_lt_20 = new TH1F("h_cent_ismb_vz_lt_20","",100,0,100);

	TChain datachain("ttree");

    const char* dataInputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
    for (int r = 0; r < 8; r++) { 
    	TString dataWildcardPath = TString::Format("%sevents_p015_%d_data_cor_*.root", dataInputDirectory, good_runs[r]);
    	//TString dataWildcardPath = TString::Format("%sevents_p011_zs_%d_data_cor_*.root", dataInputDirectory, good_runs[r]);
    	datachain.Add(dataWildcardPath);
    }
	TTreeReader datareader(&datachain);

	TTreeReaderArray<float> m_vtx(datareader, "track_vtx");
	TTreeReaderValue<int> m_sectormb(datareader, "sectormb");
	TTreeReaderArray<float> m_mbenergy(datareader, "mbenrgy");
	TTreeReaderValue<bool> m_isMinBias(datareader, "isMinBias");
	TTreeReaderValue<int> m_centbin(datareader, "centbin");

	// get data vertex distribution and mbd charge distribution
	int eventnumber = 0;
    while (datareader.Next()) {
    	if (eventnumber % 10000 == 0) cout << "data event " << eventnumber << endl;
  		float mbde = 0;
  		eventnumber++;
  		if (!*m_isMinBias) { continue; }
  		h_cent_ismb->Fill(*m_centbin);
  		if (isnan(m_vtx[2])) { continue; }
  		h_cent_ismb_notnan->Fill(*m_centbin);
  		if (m_vtx[2] > 20 || m_vtx[2] < -20) { continue; }
  		h_cent_ismb_vz_lt_20->Fill(*m_centbin);
	}

	h_cent_ismb->Write();
	h_cent_ismb_notnan->Write();
	h_cent_ismb_vz_lt_20->Write();
	out->Close();
}
