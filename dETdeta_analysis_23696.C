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

std::set<std::tuple<int, int>> emcal_hot_dead_map_23696 = {{28,83},{39,115},{80, 228},{80, 236},{80, 237},{81, 236},{64, 164},{59,183},{72,188},{76,104},{88,104},{89, 104},{90, 104},{91, 104},
									{88,105},{89, 105},{90, 105},{91, 105},{88,108},{89, 108},{90, 107},{91, 107},{91,106},{88,111},{89,110},{89,111},{48,7},{58,7},{75,27},{9,200},{13,232},{22,157},
									{24,78},{24,208},{24,209},{24,210},{24,211},{24,212},{24,213},{24,214},{24,215},{25,208},{25,209},{25,210},{25,211},{25,212},{25,213},{25,214},{25,215},{26,208},
									{26,209},{26,210},{26,211},{26,212},{26,213},{26,214},{26,215},{27,208},{27,209},{27,210},{27,211},{27,212},{27,213},{27,214},{27,215},{28,83},{28,208},{28,209},
									{28,210},{28,211},{28,212},{28,213},{28,214},{28,215},{29,208},{29,209},{29,210},{29,211},{29,212},{29,213},{29,214},{29,215},{30,208},{30,209},{30,210},{30,211},
									{30,212},{30,213},{30,214},{30,215},{31,104},{31,208},{31,209},{31,210},{31,211},{31,212},{31,213},{31,214},{31,215},{32,144},{32,145},{32,146},{32,147},{32,148},
									{32,149},{32,150},{32,151},{33,144},{33,145},{33,146},{33,147},{33,148},{33,149},{33,150},{33,151},{34,144},{34,145},{34,146},{34,147},{34,148},{34,149},{34,150},
									{34,151},{35,144},{35,145},{35,146},{35,147},{35,148},{35,149},{35,150},{35,151},{36,144},{36,145},{36,146},{36,147},{36,148},{36,149},{36,150},{36,151},{37,144},
									{37,145},{37,146},{37,147},{37,148},{37,149},{37,150},{37,151},{38,144},{38,145},{38,146},{38,147},{38,148},{38,149},{38,150},{38,151},{38,219},{39,144},{39,145},
									{39,146},{39,147},{39,148},{39,149},{39,150},{39,151},{47,80},{47,96},{47,138},{48,215},{48,231},{48,253},{48,255},{63,77},{76,186},{76,187},{77,186},{77,187},
									{78,186},{78,187},{79,186},{79,187},{80,149},{88,40},{88,41},{88,42},{88,43},{88,168},{88,169},{88,170},{88,171},{88,224},{88,225},{88,226},{88,227},{89,40},
									{89,41},{89,42},{89,43},{89,168},{89,169},{89,170},{89,171},{89,224},{89,225},{89,226},{89,227},{90,40},{90,41},{90,42},{90,43},{90,111},{90,168},{90,169},
									{90,170},{90,171},{90,224},{90,225},{90,226},{90,227},{91,40},{91,41},{91,42},{91,43},{91,168},{91,169},{91,170},{91,171},{91,224},{91,225},{91,226},{91,227},
									{92,0},{92,40},{92,41},{92,42},{92,43},{92,168},{92,169},{92,170},{92,171},{92,224},{92,225},{92,226},{92,227},{93,40},{93,41},{93,42},{93,43},{93,168},
									{93,169},{93,170},{93,171},{93,224},{93,225},{93,226},{93,227},{94,40},{94,41},{94,42},{94,43},{94,96},{94,106},{94,157},{94,168},{94,169},{94,170},{94,171},
									{94,197},{94,219},{94,224},{94,225},{94,226},{94,227}};

const double vz_hijing_reweight[40] = {2.232659, 2.097269, 2.058908, 1.995611, 1.875155, 1.844663, 1.720625, 1.695003, 1.606846, 1.580020, 1.483648, 1.519997, 1.448841, 1.322138, 
					1.324515, 1.260287, 1.211295, 1.208241, 1.138972, 1.091909, 1.057621, 1.047813, 0.967096, 0.939584, 0.903100, 0.878061, 0.818851, 0.814870, 0.767065, 0.764742, 
					0.710518, 0.687330, 0.647735, 0.615111, 0.604018, 0.569581, 0.538452, 0.506414, 0.503946, 0.490534};

const double vz_epos_reweight[40] = {2.176031, 2.183098, 2.067832, 1.894137, 1.801472, 1.771395, 1.775465, 1.630111, 1.569169, 1.706506, 1.582617, 1.449648, 1.499269, 1.335510, 
					1.394603, 1.318539, 1.214925, 1.164297, 1.112682, 1.131543, 1.070854, 1.034719, 0.985470, 0.944210, 0.918652, 0.884022, 0.833071, 0.825844, 0.771094, 0.716402, 
					0.696809, 0.680178, 0.670844, 0.620484, 0.622772, 0.586456, 0.562244, 0.529274, 0.504112, 0.468999};

const double vz_ampt_reweight[40] = {2.322869, 2.339063, 1.993545, 1.884416, 1.747539, 1.754988, 1.738317, 1.611920, 1.573347, 1.562594, 1.451963, 1.471386, 1.376459, 1.321487, 
					1.328429, 1.247453, 1.190785, 1.210302, 1.118670, 1.097008, 1.110321, 1.032524, 0.969075, 0.960457, 0.924667, 0.855760, 0.810511, 0.790480, 0.755335, 0.752490, 
					0.714910, 0.706133, 0.659971, 0.580305, 0.610281, 0.573972, 0.538400, 0.533292, 0.502327, 0.496472};

const double ihcal_eta_bin_centers[24] = {-1.0535162174898431, -0.9618021382744392, -0.8700893703958079, -0.7783717834698298,
										 -0.6866574781378039, -0.59494737316278, -0.5032464095864394, -0.4115406126317369,
										 -0.3198493439076178, -0.22815189383493045, -0.13647855478639742, -0.044804882018267274,
										 0.046855637524394804, 0.13851525981103677, 0.2301533056227234, 0.3217972649933214, 
										 0.41343376653551833, 0.5050601463927868, 0.5966824430063085, 0.6883066285326807, 
										 0.7799266452410983, 0.8715489274547986, 0.9631703160231737, 1.0547930824432712};

const double ohcal_eta_bin_centers[24] = {-1.0538018547715944, -0.9621090096478763, -0.8704159228831849, -0.7787217241822206,
										 -0.6870278330109779, -0.59533641448814, -0.5036475993078456, -0.41196332668527996,
										 -0.3202788771117685, -0.22859494086986867, -0.13691896519812005, -0.04525182544603272,
										 0.04641353992004345, 0.1380760971648327, 0.22973335216638652, 0.3213823957243743, 
										 0.413029965511396, 0.5046704423603338, 0.5963193647641828, 0.6879607567415511, 
										 0.7796011271114296, 0.8712416652855527, 0.9628809430983388, 1.0545228311327837};

const double eta_bin_centers[24] = {-1.05417,-0.9625,-0.870833,-0.779167,-0.6875,-0.595833,-0.504167,-0.4125,-0.320833,-0.229167,
	-0.1375,-0.0458333,0.0458333,0.1375,0.229167,0.320833,0.4125,0.504167,0.595833,0.6875,0.779167,0.870833,0.9625,1.05417};

const int cemcSize = 24576;
const int ihcalSize = 1536;
const int ohcalSize = 1536;
const int mbdSize = 256;
const int g4Size = 100000;
const int vtxSize = 3;

void dETdeta_analysis_23696(int dataormc = 0, int reweighting = 1, int central = 0, const char* generator = "") {	

	TDatabasePDG *_pdg = new TDatabasePDG();

	string filename = "dETdeta_analysis_23696_z=0";
  	string dattag = (dataormc?"mc":"data");
  	string weighttag = (reweighting?"reweight":"noweight");
  	string centtag = (central?"0-10":"0-90"); 
  	string gentag = generator;
  	if (!strcmp(generator,"")) filename += "_" + dattag + "_" + weighttag + "_" + centtag + ".root";
  	else filename += "_" + dattag + "_" + weighttag + "_" + centtag + "_" + gentag + ".root";

	TFile *out = new TFile(filename.c_str(),"RECREATE");
	TH1F* h_vz = new TH1F("h_vz","",400, -100, 100);
	TH1F* h_vz_reweight = new TH1F("h_vz_reweight","",400, -100, 100);
	TH1F* h_mbd = new TH1F("h_mbd","",250000,0,250000);

	TH1F* h_event_truth_energy = new TH1F("h_event_truth_energy","",10000,0,10000);
	TH1F* hetdeta = new TH1F("hetdeta","",120,-6,6);
	TH1F* hetdeta_zoom = new TH1F("hetdeta_zoom","",220,-1.1,1.1);
	
	TH2F* h_2D_ihcal_calib = new TH2F("h_2D_ihcal_calib","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_ohcal_calib = new TH2F("h_2D_ohcal_calib","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_emcal_calib = new TH2F("h_2D_emcal_calib","",96,0.,96.,256,0.,256.);

  	TH2F* h_2D_ihcal_calibT = new TH2F("h_2D_ihcal_calibT","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_ohcal_calibT = new TH2F("h_2D_ohcal_calibT","",24,0.,24.,64,0.,64.);
  	TH2F* h_2D_emcal_calibT = new TH2F("h_2D_emcal_calibT","",96,0.,96.,256,0.,256.);
	
	TH1F* h_event_energy = new TH1F("h_event_energy","", 5000,0,5000);
	TH1F* h_event_hcal_energy = new TH1F("h_event_hcal_energy","", 5000,0,5000);
	TH1F* h_event_emcal_energy = new TH1F("h_event_emcal_energy","", 5000,0,5000);
	TH1F* h_event_ihcal_energy = new TH1F("h_event_ihcal_energy","", 1000,0,1000);
	TH1F* h_event_ohcal_energy = new TH1F("h_event_ohcal_energy","", 1000,0,1000);

	TH1F* h_emcal = new TH1F("h_emcal","",1000,0,10);
	TH1F* h_ihcal = new TH1F("h_ihcal","",1000,0,10);
	TH1F* h_ohcal = new TH1F("h_ohcal","",1000,0,10);

	TH1F* h_eT_emcal = new TH1F("h_eT_emcal","",220,-1.1,1.1);

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
    TH1F* h_eT_ihcal = new TH1F("h_eT_ihcal","",ihcal_num_bins,ihcal_bin_edges);
	TH1F* h_eT_ihcal_noend = new TH1F("h_eT_ihcal_noend","",ihcal_num_bins,ihcal_bin_edges);
	TH1F* hetdeta_ihcalbin = new TH1F("hetdeta_ihcalbin","",ihcal_num_bins, ihcal_bin_edges);

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
	TH1F* h_eT_ohcal = new TH1F("h_eT_ohcal","",ohcal_num_bins,ohcal_bin_edges);
	TH1F* h_eT_ohcal_noend = new TH1F("h_eT_ohcal_noend","",ohcal_num_bins,ohcal_bin_edges);
	TH1F* hetdeta_ohcalbin = new TH1F("hetdeta_ohcalbin","",ohcal_num_bins, ohcal_bin_edges);

	TH1F* h_eta_spread_ihcal[24];
	TH1F* h_eta_spread_ohcal[24];
	for (int i = 0; i < 24; i++) {
		h_eta_spread_ihcal[i] = new TH1F(TString::Format("h_eta_spread_ihcal_%d",i),"",400,-2,2);
		h_eta_spread_ohcal[i] = new TH1F(TString::Format("h_eta_spread_ohcal_%d",i),"",400,-2,2);
	}
    
    TChain chain("ttree");

    if (dataormc && !strcmp(generator, "epos")) {
    	// location of EPOS files 
	    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/eposrun/condor/";
	    for (int i = 0; i < 1000; i++) {
	    	TString wildcardPath = TString::Format("%sOutDir%d/epos_sim_output.root", inputDirectory, i);
	    	chain.Add(wildcardPath);
	    }
    } else if (dataormc && !strcmp(generator, "plugdoor_epos")) { 
    	const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/eposrun/plugdoor_condor/";
	    for (int i = 0; i < 1000; i++) {
	    	TString wildcardPath = TString::Format("%sOutDir%d/epos_sim_output.root", inputDirectory, i);
	    	chain.Add(wildcardPath);
	    }
    } else if (dataormc && !strcmp(generator, "hijing")) {
    	// location of HIJING files 
	    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
	    for (int i = 0; i < 555; i++) {
	    	TString wildcardPath = TString::Format("%sevents_20231122_nopileup_mc_cor_%d.root", inputDirectory, i);
	    	chain.Add(wildcardPath);
	    }
	    //TString wildcardPath = "/sphenix/user/jocl/projects/sandbox/datatemp/merged_dEdeta_20231129_21615_mc_cor_555.root";
    	//chain.Add(wildcardPath);
    } else if (dataormc && !strcmp(generator, "ampt")) {
    	// location of AMPT files 
	    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/amptrun/condor/";
	    for (int i = 0; i < 1000; i++) {
	    	TString wildcardPath = TString::Format("%sOutDir%d/ampt_sim_output.root", inputDirectory , i);
	    	chain.Add(wildcardPath);
	    }
    } else if (!dataormc) {
    	const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
		TString wildcardPath = TString::Format("%sevents_20231211_p004_23696_data_cor*.root", inputDirectory); 
    	chain.Add(wildcardPath);
    } else {
    	std::cout << "generator/data not found" << std::endl;
    	return;
    }

    int m_simtwrmult_cemc;
    float m_simtwr_cemc_e[cemcSize];
    int m_simtwr_cemc_ieta[cemcSize];
    int m_simtwr_cemc_iphi[cemcSize];
    float m_simtwr_cemc_eta[cemcSize];
    int m_simtwrmult_ihcal;
    float m_simtwr_ihcal_e[ihcalSize];
    int m_simtwr_ihcal_ieta[ihcalSize];
    int m_simtwr_ihcal_iphi[ihcalSize];
    float m_simtwr_ihcal_eta[ihcalSize];
    int m_simtwrmult_ohcal;
    float m_simtwr_ohcal_e[ohcalSize];
    int m_simtwr_ohcal_ieta[ohcalSize];
    int m_simtwr_ohcal_iphi[ohcalSize];
    float m_simtwr_ohcal_eta[ohcalSize];
    int m_sectormb;
    float m_mbenergy[mbdSize];
    int m_g4;
    float m_g4_e[g4Size];
    float m_g4_eta[g4Size];
    float m_g4_pt[g4Size];
    float m_g4_pz[g4Size];
    float m_vtx[vtxSize];
    float mbd_vtx[vtxSize];

     // Set branch addresses
    chain.SetBranchAddress("sectorem", &m_simtwrmult_cemc);
    chain.SetBranchAddress("emcalen", m_simtwr_cemc_e);
    chain.SetBranchAddress("emcaletabin", m_simtwr_cemc_ieta);
    chain.SetBranchAddress("emcalphibin", m_simtwr_cemc_iphi);
    chain.SetBranchAddress("emetacor", m_simtwr_cemc_eta);

    chain.SetBranchAddress("sectorih", &m_simtwrmult_ihcal);
    chain.SetBranchAddress("ihcalen", m_simtwr_ihcal_e);
    chain.SetBranchAddress("ihcaletabin", m_simtwr_ihcal_ieta);
    chain.SetBranchAddress("ihcalphibin", m_simtwr_ihcal_iphi);
    chain.SetBranchAddress("ihetacor", m_simtwr_ihcal_eta);

    chain.SetBranchAddress("sectoroh", &m_simtwrmult_ohcal);
    chain.SetBranchAddress("ohcalen", m_simtwr_ohcal_e);
    chain.SetBranchAddress("ohcaletabin", m_simtwr_ohcal_ieta);
    chain.SetBranchAddress("ohcalphibin", m_simtwr_ohcal_iphi);
    chain.SetBranchAddress("ohetacor", m_simtwr_ohcal_eta);

    chain.SetBranchAddress("sectormb", &m_sectormb);
    chain.SetBranchAddress("mbenrgy", &m_mbenergy);

    if (dataormc) {
    	chain.SetBranchAddress("truthpar_n", &m_g4);
	    chain.SetBranchAddress("truthpar_e", m_g4_e);
	    chain.SetBranchAddress("truthpar_eta", m_g4_eta);
	    chain.SetBranchAddress("truthpar_pt", m_g4_pt);
	    chain.SetBranchAddress("truthpar_pz", m_g4_pz);
    }

    chain.SetBranchAddress("track_vtx", m_vtx);
    chain.SetBranchAddress("mbd_vtx", mbd_vtx);

	int eventnumber = 0;
	float totalweights = 0.0;
	float delta_eta = 0.09167;
	float delta_em_eta = 0.022918;

	float hijing_cent[] = {0.0,1941.0,3020.0,4478.0,6418.0,8983.0,12221.0,16172.0,21010.0,26857.0,33720.0,41673.0,50856.0,61295.0,73376.0,87385.0,103226.0,121633.0,143473.0,250000.0};
	float epos_cent[] = {0.0,2122.0,3341.0,5071.0,7311.0,9977.0,13405.0,17700.0,22914.0,29130.0,36602.0,45471.0,55881.0,67429.0,81229.0,97175.0,115926.0,138395.0,166584.0,250000.0};
	float ampt_cent[] = {0.0,2228.0,3562.0,5316.0,7688.0,10525.0,14400.0,19223.0,24427.0,30656.0,38090.0,47382.0,57843.0,69679.0,83789.0,100250.0,121187.0,142123.0,165105.0,250000.0};
	float data_cent[] = {0.0,2.0,18.0,31.0,48.0,71.0,101.0,140.0,189.0,247.0,316.0,396.0,490.0,598.0,724.0,868.0,1035.0,1227.0,1452.0,6000.0};
    std::vector<float> centrality_bin;
    if (dataormc && !strcmp(generator, "hijing")) centrality_bin.assign(hijing_cent, hijing_cent+20);
    if (dataormc && !strcmp(generator, "epos")) centrality_bin.assign(epos_cent, epos_cent+20);
    if (dataormc && !strcmp(generator, "ampt")) centrality_bin.assign(ampt_cent, ampt_cent+20);
    if (dataormc && !strcmp(generator, "plugdoor_epos")) centrality_bin.assign(epos_cent, epos_cent+20);
    if (!dataormc) centrality_bin.assign(data_cent, data_cent+20);

    assert(centrality_bin.size() == 20);

    Long64_t nEntries = chain.GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 100000; ++entry) {
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
  		// require that simulation could reconstruct a vertex for the event
  		if (m_vtx[2] <= -2.0 || m_vtx[2] >= 2.0) { continue; }

  		float totalcharge = 0.0;
  		for (int i = 0; i < m_sectormb; i++) {
  			totalcharge += m_mbenergy[i];
  		}
  		if (central && totalcharge < centrality_bin[18]) { continue; }
  		h_mbd->Fill(totalcharge);

  		h_vz->Fill(m_vtx[2]);
  		if (reweighting && dataormc && !strcmp(generator, "hijing")) { vz_weight = vz_hijing_reweight[int(floor(m_vtx[2]*2)+28)]; }
  		if (reweighting && dataormc && !strcmp(generator, "epos")) { vz_weight = vz_epos_reweight[int(floor(m_vtx[2]*2)+28)]; }
  		if (reweighting && dataormc && !strcmp(generator, "plugdoor_epos")) { vz_weight = vz_epos_reweight[int(floor(m_vtx[2]*2)+28)]; }
  		if (reweighting && dataormc && !strcmp(generator, "ampt")) { vz_weight = vz_ampt_reweight[int(floor(m_vtx[2]*2)+28)]; }
  		h_vz_reweight->Fill(m_vtx[2],vz_weight);
  		totalweights += vz_weight;

  		if (dataormc) {
			for (int i = 0; i < m_g4; i++) {
	    		float theta = atan(m_g4_pt[i] / m_g4_pz[i]);
	    		float ET = m_g4_e[i] * abs(sin(theta));
	    		hetdeta->Fill(m_g4_eta[i], ET*vz_weight);
	    		hetdeta_ihcalbin->Fill(m_g4_eta[i], ET*vz_weight);
	    		hetdeta_ohcalbin->Fill(m_g4_eta[i], ET*vz_weight);
	    		if (fabs(m_g4_eta[i]) <= 1.1) {
	    			truthe += ET;
	    			hetdeta_zoom->Fill(m_g4_eta[i], ET*vz_weight);
	    		} 
			}
			h_event_truth_energy->Fill(truthe,vz_weight);
		}
		
		for (int i = 0; i < m_simtwrmult_cemc; i++) {
			if (m_simtwr_cemc_ieta[i] < 9) { continue; }
			if (m_simtwr_cemc_ieta[i] >= 48 && m_simtwr_cemc_ieta[i] <= 95 && m_simtwr_cemc_iphi[i] >= 208 && m_simtwr_cemc_iphi[i] <= 215) { continue; }
			if (m_simtwr_cemc_ieta[i] >= 48 && m_simtwr_cemc_ieta[i] <= 95 && m_simtwr_cemc_iphi[i] >= 240 && m_simtwr_cemc_iphi[i] <= 247) { continue; }
			if (m_simtwr_cemc_ieta[i] >= 64 && m_simtwr_cemc_ieta[i] <= 71 && m_simtwr_cemc_iphi[i] >= 64 && m_simtwr_cemc_iphi[i] <= 71) { continue; }
			if (m_simtwr_cemc_ieta[i] < 48 && m_simtwr_cemc_iphi[i] < 64) { continue; }
			if (m_simtwr_cemc_ieta[i] > 94) { continue; }
			std::tuple<int, int> hot_tower = std::make_tuple(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i]);
		    auto it = emcal_hot_dead_map_23696.find(hot_tower);
		    if (it != emcal_hot_dead_map_23696.end()) { continue; }
			h_2D_emcal_calib->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight);
			h_2D_emcal_calibT->Fill(m_simtwr_cemc_ieta[i], m_simtwr_cemc_iphi[i], m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
			emcale += m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]); 
			h_emcal->Fill(m_simtwr_cemc_e[i]/cosh(m_simtwr_cemc_eta[i]), vz_weight);
			h_eT_emcal->Fill(m_simtwr_cemc_eta[i],m_simtwr_cemc_e[i]*vz_weight/cosh(m_simtwr_cemc_eta[i]));
		}


		for (int i = 0; i < m_simtwrmult_ihcal; i++) {
			if (m_simtwr_ihcal_ieta[i] == 8 && m_simtwr_ihcal_iphi[i] == 32) { continue; }
			if (m_simtwr_ihcal_ieta[i] == 7 && m_simtwr_ihcal_iphi[i] == 51) { continue; }
			h_2D_ihcal_calib->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]*vz_weight);
			h_2D_ihcal_calibT->Fill(m_simtwr_ihcal_ieta[i], m_simtwr_ihcal_iphi[i], m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
			ihcale += m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]); 
			h_ihcal->Fill(m_simtwr_ihcal_e[i]/cosh(m_simtwr_ihcal_eta[i]), vz_weight);
			h_eT_ihcal->Fill(m_simtwr_ihcal_eta[i],m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
			if (m_simtwr_ihcal_ieta[i] >= 2 && m_simtwr_ihcal_ieta[i] <= 21) h_eT_ihcal_noend->Fill(m_simtwr_ihcal_eta[i],m_simtwr_ihcal_e[i]*vz_weight/cosh(m_simtwr_ihcal_eta[i]));
			h_eta_spread_ihcal[m_simtwr_ihcal_ieta[i]]->Fill(m_simtwr_ihcal_eta[i], vz_weight);
		}


		for (int i = 0; i < m_simtwrmult_ohcal; i++) {
			h_2D_ohcal_calib->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]*vz_weight);
			h_2D_ohcal_calibT->Fill(m_simtwr_ohcal_ieta[i], m_simtwr_ohcal_iphi[i], m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
			ohcale += m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]); 
			h_ohcal->Fill(m_simtwr_ohcal_e[i]/cosh(m_simtwr_ohcal_eta[i]), vz_weight);
			h_eT_ohcal->Fill(m_simtwr_ohcal_eta[i],m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
			if (m_simtwr_ohcal_ieta[i] >= 2 && m_simtwr_ohcal_ieta[i] <= 21) h_eT_ohcal_noend->Fill(m_simtwr_ohcal_eta[i],m_simtwr_ohcal_e[i]*vz_weight/cosh(m_simtwr_ohcal_eta[i]));
			h_eta_spread_ohcal[m_simtwr_ohcal_ieta[i]]->Fill(m_simtwr_ohcal_eta[i], vz_weight);
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
	}

	h_2D_emcal_calib->Scale(1.0/totalweights);
	h_2D_emcal_calibT->Scale(1.0/totalweights);
	h_emcal->Scale(1.0/totalweights);
	h_eT_emcal->Scale(1.0/totalweights);
	h_eT_emcal->Scale(1.0/0.01);

	h_2D_ihcal_calib->Scale(1.0/totalweights);
	h_2D_ihcal_calibT->Scale(1.0/totalweights);
	h_ihcal->Scale(1.0/totalweights);
	h_eT_ihcal->Scale(1.0/totalweights);
	h_eT_ihcal_noend->Scale(1.0/totalweights);
	for (int i = 1; i <= ihcal_num_bins; ++i) {
        double bin_width = ihcal_bin_edges[i] - ihcal_bin_edges[i - 1];
        h_eT_ihcal->SetBinContent(i, h_eT_ihcal->GetBinContent(i) / bin_width);
        h_eT_ihcal_noend->SetBinContent(i, h_eT_ihcal_noend->GetBinContent(i) / bin_width);
        if (dataormc) hetdeta_ihcalbin->SetBinContent(i, hetdeta_ihcalbin->GetBinContent(i)/ bin_width);
    }

	h_2D_ohcal_calib->Scale(1.0/totalweights);
	h_2D_ohcal_calibT->Scale(1.0/totalweights);
	h_ohcal->Scale(1.0/totalweights);
	h_eT_ohcal->Scale(1.0/totalweights);
	h_eT_ohcal_noend->Scale(1.0/totalweights);
	for (int i = 1; i <= ohcal_num_bins; ++i) {
        double bin_width = ohcal_bin_edges[i] - ohcal_bin_edges[i - 1];
        h_eT_ohcal->SetBinContent(i, h_eT_ohcal->GetBinContent(i) / bin_width);
        h_eT_ohcal_noend->SetBinContent(i, h_eT_ohcal_noend->GetBinContent(i) / bin_width);
        if (dataormc) hetdeta_ohcalbin->SetBinContent(i, hetdeta_ohcalbin->GetBinContent(i) / bin_width);
    }

	for (int i = 0; i < 24; i++) {
		h_eta_spread_ihcal[i]->Scale(1.0/totalweights);
		h_eta_spread_ohcal[i]->Scale(1.0/totalweights);
		h_eta_spread_ihcal[i]->Scale(1.0/0.01);
		h_eta_spread_ohcal[i]->Scale(1.0/0.01);
	}

	h_event_energy->Scale(1.0/totalweights);
	h_event_hcal_energy->Scale(1.0/totalweights);
	h_event_emcal_energy->Scale(1.0/totalweights);
	h_event_ihcal_energy->Scale(1.0/totalweights);
	h_event_ohcal_energy->Scale(1.0/totalweights);
	
	out->Write();
	out->Close();


}
