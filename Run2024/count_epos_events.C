#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TH1.h>

void count_epos_events(string outfile = "epos_reweighted_minbias_info.root") {
    // Input directory and file pattern
    const char* mcInputDirectory = "/sphenix/tg/tg01/commissioning/CaloCalibWG/egm2153/detdeta_run24auau/";
    TChain mcchain("ttree");  // Adjust "ttree" to your actual TTree name if different

    // Add files to TChain
    for (int i = 0; i < 5000; i++) {
        TString mcWildcardPath = TString::Format("%sevents_epos_reweighted_run14_waveformsim_nonoise_nozs_fixed_build_epos_cor_%d.root", mcInputDirectory, i);
        mcchain.Add(mcWildcardPath);
    }

    mcchain.SetBranchStatus("*", 0);

    // Variables to read branches
    float mbenrgy[128];  // Assuming max size of 128 (64+64)
    float mbtime[128];
    float bimp;
    int npart;

    // Enable only the required branches
    mcchain.SetBranchStatus("mbenrgy", 1);
    mcchain.SetBranchStatus("mbtime", 1);
    mcchain.SetBranchStatus("bimp", 1);
    mcchain.SetBranchStatus("npart", 1);

    // Set branch addresses
    mcchain.SetBranchAddress("mbenrgy", mbenrgy);
    mcchain.SetBranchAddress("mbtime", mbtime);
    mcchain.SetBranchAddress("bimp", &bimp);
    mcchain.SetBranchAddress("npart", &npart);

    // Counters
    int total_events = 0;
    int passed_events = 0;
    int notime_passed_events = 0;
    int low_passed_events = 0;

    TH1F* h_total_events = new TH1F("h_total_events","",1,0,1);
    TH1F* h_passed_events = new TH1F("h_passed_events","",1,0,1);
    TH1F* h_total_bimp = new TH1F("h_total_bimp","",300,0,30);
    TH1F* h_passed_bimp = new TH1F("h_passed_bimp","",300,0,30);
    
    TH1F* h_notime_passed_events = new TH1F("h_notime_passed_events","",1,0,1);
    TH1F* h_low_passed_events = new TH1F("h_low_passed_events","",1,0,1);
    TH1F* h_notime_passed_bimp = new TH1F("h_notime_passed_bimp","",300,0,30);
    TH1F* h_low_passed_bimp = new TH1F("h_low_passed_bimp","",300,0,30);

    TH1F* h_total_npart = new TH1F("h_total_npart","",400,0,400);
    TH1F* h_low_passed_npart = new TH1F("h_low_passed_npart","",400,0,400);

    // Loop over events
    Long64_t nEntries = mcchain.GetEntries();
    std::cout << "Total events in chain: " << nEntries << std::endl;

    for (Long64_t entry = 0; entry < nEntries; entry++) {
        mcchain.GetEntry(entry);
        if (entry % 10000 == 0) { std::cout << entry << std::endl;}

        int mbd_nhits1 = 0;
        int mbd_nhits2 = 0;

	int mbd_low_nhits1 = 0;
	int mbd_low_nhits2 = 0; 

        int mbd_notime_nhits1 = 0;
	int mbd_notime_nhits2 = 0;
	
	h_total_bimp->Fill(bimp);
	h_total_npart->Fill(npart);

        for (int i = 0; i < 64; i++) {
            if (mbenrgy[i] > 0.5 && mbtime[i] < 25.0) { mbd_nhits1 += 1; }
            if (mbenrgy[i+64] > 0.5 && mbtime[i] < 25.0) { mbd_nhits2 += 1; }
       	    if (mbenrgy[i] > 0.4) { mbd_low_nhits1 += 1; }
	    if (mbenrgy[i+64] > 0.4) {mbd_low_nhits2 += 1; }
	    if (mbenrgy[i] > 0.5) { mbd_notime_nhits1 += 1; }
	    if (mbenrgy[i+64] > 0.5) { mbd_notime_nhits2 += 1; }    
	}

        // Apply event selection
        if (mbd_nhits1 > 2 && mbd_nhits2 > 2) {
            h_passed_bimp->Fill(bimp);
            passed_events++;
        }

	if (mbd_notime_nhits1 > 2 && mbd_notime_nhits2 > 2) {
		h_notime_passed_bimp->Fill(bimp);
		notime_passed_events++;
	}

	if (mbd_low_nhits1 > 2 && mbd_low_nhits2 > 2) {
		h_low_passed_bimp->Fill(bimp);
		h_low_passed_npart->Fill(npart);
		low_passed_events++;
	}
        total_events++;
    }

    // Print results
    std::cout << "Total events processed: " << total_events << std::endl;
    std::cout << "Events passing cut (mbd_nhits1 > 2 && mbd_nhits2 > 2): " << passed_events << std::endl;
    std::cout << "Events passing no time cut: " << notime_passed_events << std::endl;
    std::cout << "Events passing low cut: " << low_passed_events << std::endl;

    h_total_events->SetBinContent(1, total_events);
    h_passed_events->SetBinContent(1, passed_events);
    h_notime_passed_events->SetBinContent(1, notime_passed_events);
    h_low_passed_events->SetBinContent(1, low_passed_events);

    TFile *out = new TFile(outfile.c_str(),"RECREATE");
    out->cd();
    h_total_events->Write();
    h_passed_events->Write();
    h_notime_passed_events->Write();
    h_low_passed_events->Write();
    h_total_bimp->Write();
    h_passed_bimp->Write();
    h_notime_passed_bimp->Write();
    h_low_passed_bimp->Write();
    h_total_npart->Write();
    h_low_passed_npart->Write();
    out->Close();

}

