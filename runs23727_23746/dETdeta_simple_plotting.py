import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import numpy as np
import pdb

tag = 'p015_w_hcal'
cents = ['0-5','5-10','10-20','20-30','30-40','40-50','50-60']
for cent in cents:
	hijingfile = 'MC/dETdeta_analysis_allruns_hijing_run14_nozs_mc_reweight_'+cent+'_reweight_hijing.root'
	datafile = 'updated_centrality/dETdeta_analysis_allruns_test_cent_HCal_zs_30ADC_EMCal_zs_40ADC_data_noweight_'+cent+'.root'
	outfile = 'dETdeta_plots_w_hijing_reweight_updated_centrality_'+tag+'_'+cent+'.root'

	f1 = ROOT.TFile.Open(hijingfile)
	h_eT_truth = TH1F(f1.Get("hetdeta_ihcalbin"))
	h_eT_truth_ihcalbin = TH1F(f1.Get("hetdeta_ihcalbin"))
	h_eT_truth_ohcalbin = TH1F(f1.Get("hetdeta_ohcalbin"))
	h_eT_truth_calobin = TH1F(f1.Get("hetdeta_calobin"))
	h_eT_sim_emcal = TH1F(f1.Get("h_eT_emcal"))
	h_eT_sim_ihcal = TH1F(f1.Get("h_eT_ihcal"))
	h_eT_sim_ohcal = TH1F(f1.Get("h_eT_ohcal"))
	h_eT_sim_calo = TH1F(f1.Get("h_eT_calo"))
	h_eT_truth.SetDirectory(0)
	h_eT_truth_ihcalbin.SetDirectory(0)
	h_eT_truth_ohcalbin.SetDirectory(0)
	h_eT_truth_calobin.SetDirectory(0)
	h_eT_sim_emcal.SetDirectory(0)
	h_eT_sim_ihcal.SetDirectory(0)
	h_eT_sim_ohcal.SetDirectory(0)
	h_eT_sim_calo.SetDirectory(0)
	f1.Close()
	f2 = ROOT.TFile.Open(datafile)
	h_eT_data_emcal = TH1F(f2.Get("h_eT_emcal"))
	h_eT_data_ihcal = TH1F(f2.Get("h_eT_ihcal"))
	h_eT_data_ohcal = TH1F(f2.Get("h_eT_ohcal"))
	h_eT_data_calo = TH1F(f2.Get("h_eT_calo"))
	h_eT_data_emcal.SetDirectory(0)
	h_eT_data_ihcal.SetDirectory(0)
	h_eT_data_ohcal.SetDirectory(0)
	h_eT_data_calo.SetDirectory(0)
	f2.Close()

	file = ROOT.TFile(outfile, "RECREATE")

	emcal_ratio_hijing = TH1F(h_eT_sim_emcal.Clone("emcal_ratio_hijing"))
	emcal_ratio_hijing.Divide(h_eT_truth)
	emcal_detdeta_hijing = TH1F(h_eT_data_emcal.Clone("emcal_detdeta_hijing"))
	emcal_detdeta_hijing.Divide(emcal_ratio_hijing)
	ihcal_ratio_hijing = TH1F(h_eT_sim_ihcal.Clone("ihcal_ratio_hijing"))
	ihcal_ratio_hijing.Divide(h_eT_truth_ihcalbin)
	ihcal_detdeta_hijing = TH1F(h_eT_data_ihcal.Clone("ihcal_detdeta_hijing"))
	ihcal_detdeta_hijing.Divide(ihcal_ratio_hijing)
	ohcal_ratio_hijing = TH1F(h_eT_sim_ohcal.Clone("ohcal_ratio_hijing"))
	ohcal_ratio_hijing.Divide(h_eT_truth_ohcalbin)
	ohcal_detdeta_hijing = TH1F(h_eT_data_ohcal.Clone("ohcal_detdeta_hijing"))
	ohcal_detdeta_hijing.Divide(ohcal_ratio_hijing)
	calo_ratio_hijing = TH1F(h_eT_sim_calo.Clone("calo_ratio_hijing"))
	calo_ratio_hijing.Divide(h_eT_truth_calobin)
	calo_detdeta_hijing = TH1F(h_eT_data_calo.Clone("calo_detdeta_hijing"))
	calo_detdeta_hijing.Divide(calo_ratio_hijing)

	h_eT_data_hcal = TH1F(h_eT_data_ihcal.Clone("h_eT_data_hcal"))
	for i in range(1, h_eT_data_hcal.GetNbinsX() + 1):
		h_eT_data_hcal.SetBinContent(i, h_eT_data_ihcal.GetBinContent(i) + h_eT_data_ohcal.GetBinContent(i))
		h_eT_data_hcal.SetBinError(i, np.sqrt(h_eT_data_ihcal.GetBinError(i)**2 + h_eT_data_ohcal.GetBinError(i)**2))
	h_eT_sim_hcal = TH1F(h_eT_sim_ihcal.Clone("h_eT_sim_hcal"))
	for i in range(1, h_eT_sim_hcal.GetNbinsX() + 1):
		h_eT_sim_hcal.SetBinContent(i, h_eT_sim_ihcal.GetBinContent(i) + h_eT_sim_ohcal.GetBinContent(i))
		h_eT_sim_hcal.SetBinError(i, np.sqrt(h_eT_sim_ihcal.GetBinError(i)**2 + h_eT_sim_ohcal.GetBinError(i)**2))
	hcal_ratio_hijing = TH1F(h_eT_sim_hcal.Clone("hcal_ratio_hijing"))
	hcal_ratio_hijing.Divide(h_eT_truth_ihcalbin)
	hcal_detdeta_hijing = TH1F(h_eT_data_hcal.Clone("hcal_detdeta_hijing"))
	hcal_detdeta_hijing.Divide(hcal_ratio_hijing)

	emcal_detdeta_hijing.Write()
	ihcal_detdeta_hijing.Write()
	ohcal_detdeta_hijing.Write()
	calo_detdeta_hijing.Write()
	h_eT_sim_hcal.Write()
	h_eT_data_hcal.Write()
	hcal_detdeta_hijing.Write()

	file.Close()
