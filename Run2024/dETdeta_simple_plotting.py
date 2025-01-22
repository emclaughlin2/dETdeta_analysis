import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import numpy as np
import pdb

tag = 'ana450_2024p009_w_hcal_tsc'
cents = ['0-5','5-10','10-20','20-30','30-40','40-50','50-60']

pedestalfile = 'pedestal_subtraction/dETdeta_analysis_pedestal_subtraction_54256_ana450_default_bins.root'
f3 = ROOT.TFile.Open(pedestalfile)
h_eT_pedestal_emcal = TH1F(f3.Get("h_eT_eta_emcal_profile_hist"))
h_eT_pedestal_ihcal = TH1F(f3.Get("h_eT_eta_ihcal_profile_hist"))
h_eT_pedestal_ohcal = TH1F(f3.Get("h_eT_eta_ohcal_profile_hist"))
h_eT_pedestal_calo = TH1F(f3.Get("h_eT_eta_calo_profile_hist"))
h_eT_pedestal_emcal.SetDirectory(0)
h_eT_pedestal_ihcal.SetDirectory(0)
h_eT_pedestal_ohcal.SetDirectory(0)
h_eT_pedestal_calo.SetDirectory(0)
f3.Close()

for cent in cents:
	hijingfile = 'tprofile/dETdeta_analysis_allruns_run14_with_npart_nozs_mc_reweight_'+cent+'_reweight_hijing.root'
	datafile = 'tprofile/dETdeta_analysis_allruns_ana450_2024p009_100_50_50_ZS_hcal_tsc_nozs_data_noweight_'+cent+'.root'
	outfile = 'tprofile/dETdeta_plots_'+tag+'_'+cent+'.root'

	f1 = ROOT.TFile.Open(hijingfile)
	h_eT_truth = TH1F(f1.Get("hetdeta_ihcalbin"))
	h_eT_truth_ihcalbin = TH1F(f1.Get("hetdeta_ihcalbin"))
	h_eT_truth_ohcalbin = TH1F(f1.Get("hetdeta_ohcalbin"))
	h_eT_truth_calobin = TH1F(f1.Get("hetdeta_calobin"))
	h_eT_sim_emcal = TH1F(f1.Get("h_eT_eta_emcal_profile_hist"))
	h_eT_sim_ihcal = TH1F(f1.Get("h_eT_eta_ihcal_profile_hist"))
	h_eT_sim_ohcal = TH1F(f1.Get("h_eT_eta_ohcal_profile_hist"))
	h_eT_sim_calo = TH1F(f1.Get("h_eT_eta_calo_profile_hist"))
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
	h_eT_data_emcal = TH1F(f2.Get("h_eT_eta_emcal_profile_hist"))
	h_eT_data_ihcal = TH1F(f2.Get("h_eT_eta_ihcal_profile_hist"))
	h_eT_data_ohcal = TH1F(f2.Get("h_eT_eta_ohcal_profile_hist"))
	h_eT_data_calo = TH1F(f2.Get("h_eT_eta_calo_profile_hist"))
	h_eT_data_emcal.SetDirectory(0)
	h_eT_data_ihcal.SetDirectory(0)
	h_eT_data_ohcal.SetDirectory(0)
	h_eT_data_calo.SetDirectory(0)
	f2.Close()

	file = ROOT.TFile(outfile, "RECREATE")

	for i in range(1, h_eT_data_emcal.GetNbinsX() + 1):
		h_eT_data_emcal.SetBinContent(i, h_eT_data_emcal.GetBinContent(i) - h_eT_pedestal_emcal.GetBinContent(i))
	for i in range(1, h_eT_data_ihcal.GetNbinsX() + 1):
		h_eT_data_ihcal.SetBinContent(i, h_eT_data_ihcal.GetBinContent(i) - h_eT_pedestal_ihcal.GetBinContent(i))
	for i in range(1, h_eT_data_ohcal.GetNbinsX() + 1):
		h_eT_data_ohcal.SetBinContent(i, h_eT_data_ohcal.GetBinContent(i) - h_eT_pedestal_ohcal.GetBinContent(i))
	for i in range(1, h_eT_data_calo.GetNbinsX() + 1):
		h_eT_data_calo.SetBinContent(i, h_eT_data_calo.GetBinContent(i) - h_eT_pedestal_calo.GetBinContent(i))

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
