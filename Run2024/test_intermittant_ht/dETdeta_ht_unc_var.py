import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import numpy as np
import pdb
import sys

cents = ['0-1','1-2','2-3','3-4','4-5','55-56','56-57','57-58','58-59','59-60']
for cent in cents:
	datafile = 'dETdeta_analysis_allruns_ana450_2024p009_100_50_50_ZS_hcal_tsc_intermittant_ht_events_nozs_data_noweight_'+cent+'.root'
	vardatafile = 'dETdeta_analysis_allruns_ana450_2024p009_100_50_50_ZS_hcal_tsc_no_intermittant_ht_events_nozs_data_noweight_'+cent+'.root'
	outfile = 'dETdeta_variation_tprofile_ht_'+cent+'.root'
	f2 = ROOT.TFile.Open(datafile)
	h_eT_data_emcal = TH1F(f2.Get("h_eT_eta_emcal_profile_hist"))
	h_eT_data_ihcal = TH1F(f2.Get("h_eT_eta_ihcal_profile_hist"))
	h_eT_data_ohcal = TH1F(f2.Get("h_eT_eta_ohcal_profile_hist"))
	h_eT_data_calo = TH1F(f2.Get("h_eT_eta_calo_profile_hist"))
	h_eT_data_hcal = TH1F(f2.Get("h_eT_eta_hcal_profile_hist"))
	h_eT_data_emcal.SetDirectory(0)
	h_eT_data_ihcal.SetDirectory(0)
	h_eT_data_ohcal.SetDirectory(0)
	h_eT_data_calo.SetDirectory(0)
	h_eT_data_hcal.SetDirectory(0)
	f2.Close()
	f5 = ROOT.TFile.Open(vardatafile)
	h_eT_var_emcal = TH1F(f5.Get("h_eT_eta_emcal_profile_hist"))
	h_eT_var_ihcal = TH1F(f5.Get("h_eT_eta_ihcal_profile_hist"))
	h_eT_var_ohcal = TH1F(f5.Get("h_eT_eta_ohcal_profile_hist"))
	h_eT_var_calo = TH1F(f5.Get("h_eT_eta_calo_profile_hist"))
	h_eT_var_hcal = TH1F(f5.Get("h_eT_eta_hcal_profile_hist"))
	h_eT_var_emcal.SetDirectory(0)
	h_eT_var_ihcal.SetDirectory(0)
	h_eT_var_ohcal.SetDirectory(0)
	h_eT_var_calo.SetDirectory(0)
	h_eT_var_hcal.SetDirectory(0)
	f5.Close()

	file = ROOT.TFile(outfile, "RECREATE")

	# find the deviation from nominal 
	emcal_reco_dev = TH1F(h_eT_data_emcal.Clone("emcal_reco_dev"))
	emcal_reco_dev.SetXTitle("#eta")
	emcal_reco_dev.SetYTitle("Reco dE_{T}/d#eta dev. [GeV]")
	for i in range(1, emcal_reco_dev.GetNbinsX() + 1):
		emcal_reco_dev.SetBinContent(i, np.abs(h_eT_data_emcal.GetBinContent(i) - h_eT_var_emcal.GetBinContent(i)))
		emcal_reco_dev.SetBinError(i, np.sqrt(h_eT_data_emcal.GetBinError(i)**2 + h_eT_var_emcal.GetBinError(i)**2))
	emcal_reco_dev.Divide(h_eT_data_emcal)
	ihcal_reco_dev = TH1F(h_eT_data_ihcal.Clone("ihcal_reco_dev"))
	ihcal_reco_dev.SetXTitle("#eta")
	ihcal_reco_dev.SetYTitle("Reco dE_{T}/d#eta dev. [GeV]")
	for i in range(1, ihcal_reco_dev.GetNbinsX() + 1):
		ihcal_reco_dev.SetBinContent(i, np.abs(h_eT_data_ihcal.GetBinContent(i) - h_eT_var_ihcal.GetBinContent(i)))
		ihcal_reco_dev.SetBinError(i, np.sqrt(h_eT_data_ihcal.GetBinError(i)**2 + h_eT_var_ihcal.GetBinError(i)**2))
	ihcal_reco_dev.Divide(h_eT_data_ihcal)
	ohcal_reco_dev = TH1F(h_eT_data_ohcal.Clone("ohcal_reco_dev"))
	ohcal_reco_dev.SetXTitle("#eta")
	ohcal_reco_dev.SetYTitle("Reco dE_{T}/d#eta dev. [GeV]")
	for i in range(1, ohcal_reco_dev.GetNbinsX() + 1):
		ohcal_reco_dev.SetBinContent(i, np.abs(h_eT_data_ohcal.GetBinContent(i) - h_eT_var_ohcal.GetBinContent(i)))
		ohcal_reco_dev.SetBinError(i, np.sqrt(h_eT_data_ohcal.GetBinError(i)**2 + h_eT_var_ohcal.GetBinError(i)**2))
	ohcal_reco_dev.Divide(h_eT_data_ohcal)
	
	calo_reco_dev = TH1F(h_eT_data_calo.Clone("calo_reco_dev"))
	calo_reco_dev.SetXTitle("#eta")
	calo_reco_dev.SetYTitle("Reco dE_{T}/d#eta dev. [GeV]")
	for i in range(1, calo_reco_dev.GetNbinsX() + 1):
		calo_reco_dev.SetBinContent(i, np.abs(h_eT_data_calo.GetBinContent(i) - h_eT_var_calo.GetBinContent(i)))
		calo_reco_dev.SetBinError(i, np.sqrt(h_eT_data_calo.GetBinError(i)**2 + h_eT_var_calo.GetBinError(i)**2))
	calo_reco_dev.Divide(h_eT_data_calo)
		
	hcal_reco_dev = TH1F(h_eT_data_hcal.Clone("hcal_reco_dev"))
	hcal_reco_dev.SetXTitle("#eta")
	hcal_reco_dev.SetYTitle("Reco dE_{T}/d#eta dev. [GeV]")
	for i in range(1, hcal_reco_dev.GetNbinsX() + 1):
		hcal_reco_dev.SetBinContent(i, np.abs(h_eT_data_hcal.GetBinContent(i) - h_eT_var_hcal.GetBinContent(i)))
		hcal_reco_dev.SetBinError(i, np.sqrt(h_eT_data_hcal.GetBinError(i)**2 + h_eT_var_hcal.GetBinError(i)**2))
	hcal_reco_dev.Divide(h_eT_data_hcal)

	emcal_reco_dev.Write()
	ihcal_reco_dev.Write()
	ohcal_reco_dev.Write()
	calo_reco_dev.Write()
	hcal_reco_dev.Write()
	file.Close()
