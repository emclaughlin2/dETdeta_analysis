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

void dETdeta_vertex_reweighting(int runnumber) {	

	TFile *out = new TFile(TString::Format("dETdeta_vertex_reweight_run%d.root",runnumber),"RECREATE");
	
	TH1F* h_vz_data = new TH1F("h_vz_data","",200,-50,50);
	TH1F* h_vz_mc = new TH1F("h_vz_mc","",200,-50,50);

	TH1F* h_mbd_data = new TH1F("h_mbd_data","",12000,0,6000);
	TH1F* h_mbd_mc = new TH1F("h_mbd_mc","",250000,0,250000);

	TChain mcchain("ttree");
	TChain datachain("ttree");
    
    const char* dataInputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
    for (int i = 0; i < 5; i++) {
    	TString dataWildcardPath = TString::Format("%sevents_20231211_p004_23696_data_cor_%d.root", dataInputDirectory, i);
    	datachain.Add(dataWildcardPath);
    }
	TTreeReader datareader(&datachain);

    const char* mcInputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
	for (int i = 0; i < 5; i++) {
	    	TString mcWildcardPath = TString::Format("%sevents_20231122_nopileup_mc_cor_%d.root", mcInputDirectory, i);
	    	mcchain.Add(mcWildcardPath);
	}
	TTreeReader mcreader(&mcchain);

	TTreeReaderArray<float> m_vtx(datareader, "track_vtx");
	TTreeReaderValue<int> m_sectormb(datareader, "sectormb");
	TTreeReaderArray<float> m_mbenergy(datareader, "mbenrgy");

	TTreeReaderArray<float> track_vtx(mcreader, "track_vtx");
	TTreeReaderValue<int> sectormb(mcreader, "sectormb");
	TTreeReaderArray<float> mbenergy(mcreader, "mbenrgy");

	// get data vertex distribution and mbd charge distribution
	int eventnumber = 0;
    while (datareader.Next()) {
    	if (eventnumber % 1000 == 0) cout << "event " << eventnumber << endl;
  		float mbde = 0;
  		eventnumber++;
  		h_vz_data->Fill(m_vtx[2]);

  		for (int i = 0; i < *m_sectormb; i++) {
  			mbde += m_mbenergy[i];
  		}
  		h_mbd_data->Fill(mbde);	
	}

	// get mc vertex distribution and mbd charge distribution
	while (mcreader.Next()) {
    	if (eventnumber % 1000 == 0) cout << "event " << eventnumber << endl;
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
	double mc_q[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	double data_q[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	float vertex_reweight[200] = {0.};
    T->Branch("data_centrality", data_q, "data_q[18]/F");
    T->Branch("mc_centrality", mc_q, "mc_q[20]/F");
    T->Branch("vertex_reweight", vertex_reweight, "vertex_reweight[200]/F");

	int mc_nq = 21;
	double mc_xq[] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
	h_mbd_mc->GetQuantiles(mc_nq,mc_q,mc_xq);

	int data_nq = 19;
	double data_xq[] = {0.0,0.0555,0.11111,0.16666,0.22222,0.27777,0.33333,0.38888,0.44444,0.5,0.55555,0.61111,0.66666,0.72222,0.77777,0.83333,0.88888,0.94444,1.0};
	h_mbd_data->GetQuantiles(data_nq,data_q,data_xq);


	for (int i = 1; i < vz_ratio->GetNbinsX() + 1; i++) {
		vertex_reweight[i-1] = vz_ratio->GetBinContent(i);
	}

	T->Fill();
	T->Write();
	h_vz_data->Write();
	h_vz_mc->Write();
	vz_ratio->Write();
	h_mbd_data->Write();
	h_mbd_mc->Write();
	
	out->Close();
}
