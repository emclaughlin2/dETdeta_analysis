import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import numpy as np
import pdb
import sys

cents = ['0-5','5-10','10-20','20-30','30-40','40-50','50-60']
tags = ['ihsyst1','ihsyst2','ihsyst3','ihsyst4']
for cent in cents:
	for tag in tags:
		mcfile = 'MC/dETdeta_analysis_allruns_hijing_run14_nozs_mc_reweight_'+cent+'_reweight_hijing.root'
		datafile = 'new_cent_bins/dETdeta_analysis_allruns_HCal_zs_30ADC_EMCal_zs_40ADC_data_noweight_'+cent+'.root'
		vardatafile = 'hcal_syst/dETdeta_analysis_allruns_p015_'+tag+'_HCal_zs_30ADC_EMCal_zs_40ADC_data_noweight_'+cent+'.root'
		outfile = 'hcal_syst/dETdeta_variation_'+tag+'_w_hcal_'+cent+'.root'

		f1 = ROOT.TFile.Open(mcfile)
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
		f5 = ROOT.TFile.Open(vardatafile)
		h_eT_var_emcal = TH1F(f5.Get("h_eT_emcal"))
		h_eT_var_ihcal = TH1F(f5.Get("h_eT_ihcal"))
		h_eT_var_ohcal = TH1F(f5.Get("h_eT_ohcal"))
		h_eT_var_calo = TH1F(f5.Get("h_eT_calo"))
		h_eT_var_emcal.SetDirectory(0)
		h_eT_var_ihcal.SetDirectory(0)
		h_eT_var_ohcal.SetDirectory(0)
		h_eT_var_calo.SetDirectory(0)
		f5.Close()

		file = ROOT.TFile(outfile, "RECREATE")

		# find the deviation from nominal 
		emcal_reco_dev = TH1F(h_eT_data_emcal.Clone("emcal_reco_dev"))
		emcal_reco_dev.SetXTitle("#eta")
		emcal_reco_dev.SetYTitle("Reco dE_{T}/d#eta dev. [GeV]")
		for i in range(1, emcal_reco_dev.GetNbinsX() + 1):
			emcal_reco_dev.SetBinContent(i, np.abs(h_eT_data_emcal.GetBinContent(i) - h_eT_var_emcal.GetBinContent(i)))
			emcal_reco_dev.SetBinError(i, np.sqrt(h_eT_data_emcal.GetBinError(i)**2 + h_eT_var_emcal.GetBinError(i)**2))
		ihcal_reco_dev = TH1F(h_eT_data_ihcal.Clone("ihcal_reco_dev"))
		ihcal_reco_dev.SetXTitle("#eta")
		ihcal_reco_dev.SetYTitle("Reco dE_{T}/d#eta dev. [GeV]")
		for i in range(1, ihcal_reco_dev.GetNbinsX() + 1):
			ihcal_reco_dev.SetBinContent(i, np.abs(h_eT_data_ihcal.GetBinContent(i) - h_eT_var_ihcal.GetBinContent(i)))
			ihcal_reco_dev.SetBinError(i, np.sqrt(h_eT_data_ihcal.GetBinError(i)**2 + h_eT_var_ihcal.GetBinError(i)**2))
		ohcal_reco_dev = TH1F(h_eT_data_ohcal.Clone("ohcal_reco_dev"))
		ohcal_reco_dev.SetXTitle("#eta")
		ohcal_reco_dev.SetYTitle("Reco dE_{T}/d#eta dev. [GeV]")
		for i in range(1, ohcal_reco_dev.GetNbinsX() + 1):
			ohcal_reco_dev.SetBinContent(i, np.abs(h_eT_data_ohcal.GetBinContent(i) - h_eT_var_ohcal.GetBinContent(i)))
			ohcal_reco_dev.SetBinError(i, np.sqrt(h_eT_data_ohcal.GetBinError(i)**2 + h_eT_var_ohcal.GetBinError(i)**2))
		calo_reco_dev = TH1F(h_eT_data_calo.Clone("calo_reco_dev"))
		calo_reco_dev.SetXTitle("#eta")
		calo_reco_dev.SetYTitle("Reco dE_{T}/d#eta dev. [GeV]")
		for i in range(1, calo_reco_dev.GetNbinsX() + 1):
			calo_reco_dev.SetBinContent(i, np.abs(h_eT_data_calo.GetBinContent(i) - h_eT_var_calo.GetBinContent(i)))
			calo_reco_dev.SetBinError(i, np.sqrt(h_eT_data_calo.GetBinError(i)**2 + h_eT_var_calo.GetBinError(i)**2))
		

		h_eT_data_hcal = TH1F(h_eT_data_ihcal.Clone("h_eT_data_hcal"))
		for i in range(1, h_eT_data_hcal.GetNbinsX() + 1):
			h_eT_data_hcal.SetBinContent(i, h_eT_data_ihcal.GetBinContent(i) + h_eT_data_ohcal.GetBinContent(i))
			h_eT_data_hcal.SetBinError(i, np.sqrt(h_eT_data_ihcal.GetBinError(i)**2 + h_eT_data_ohcal.GetBinError(i)**2))
		h_eT_sim_hcal = TH1F(h_eT_sim_ihcal.Clone("h_eT_sim_hcal"))
		for i in range(1, h_eT_sim_hcal.GetNbinsX() + 1):
			h_eT_sim_hcal.SetBinContent(i, h_eT_sim_ihcal.GetBinContent(i) + h_eT_sim_ohcal.GetBinContent(i))
			h_eT_sim_hcal.SetBinError(i, np.sqrt(h_eT_sim_ihcal.GetBinError(i)**2 + h_eT_sim_ohcal.GetBinError(i)**2))
		h_eT_var_hcal = TH1F(h_eT_var_ihcal.Clone("h_eT_var_hcal"))
		for i in range(1, h_eT_var_hcal.GetNbinsX() + 1):
			h_eT_var_hcal.SetBinContent(i, h_eT_var_ihcal.GetBinContent(i) + h_eT_var_ohcal.GetBinContent(i))
			h_eT_var_hcal.SetBinError(i, np.sqrt(h_eT_var_ihcal.GetBinError(i)**2 + h_eT_var_ohcal.GetBinError(i)**2))

		# correction factors from MC
		emcal_ratio_mc = TH1F(h_eT_sim_emcal.Clone("emcal_ratio_mc"))
		emcal_ratio_mc.Divide(h_eT_truth)
		ihcal_ratio_mc = TH1F(h_eT_sim_ihcal.Clone("ihcal_ratio_mc"))
		ihcal_ratio_mc.Divide(h_eT_truth_ihcalbin)
		ohcal_ratio_mc = TH1F(h_eT_sim_ohcal.Clone("ohcal_ratio_mc"))
		ohcal_ratio_mc.Divide(h_eT_truth_ohcalbin)
		calo_ratio_mc = TH1F(h_eT_sim_calo.Clone("calo_ratio_mc"))
		calo_ratio_mc.Divide(h_eT_truth_calobin)
		hcal_ratio_mc = TH1F(h_eT_sim_hcal.Clone("hcal_ratio_mc"))
		hcal_ratio_mc.Divide(h_eT_truth_ihcalbin)

		# fully corrected histograms for nominal data
		emcal_detdeta_mc = TH1F(h_eT_data_emcal.Clone("emcal_detdeta_mc"))
		emcal_detdeta_mc.Divide(emcal_ratio_mc)
		ihcal_detdeta_mc = TH1F(h_eT_data_ihcal.Clone("ihcal_detdeta_mc"))
		ihcal_detdeta_mc.Divide(ihcal_ratio_mc)
		ohcal_detdeta_mc = TH1F(h_eT_data_ohcal.Clone("ohcal_detdeta_mc"))
		ohcal_detdeta_mc.Divide(ohcal_ratio_mc)
		calo_detdeta_mc = TH1F(h_eT_data_calo.Clone("calo_detdeta_mc"))
		calo_detdeta_mc.Divide(calo_ratio_mc)
		hcal_detdeta_mc = TH1F(h_eT_data_hcal.Clone("hcal_detdeta_mc"))
		hcal_detdeta_mc.Divide(hcal_ratio_mc)	

		# fully corrected histograms for variation data
		emcal_var_detdeta_mc = TH1F(h_eT_var_emcal.Clone("emcal_var_detdeta_mc"))
		emcal_var_detdeta_mc.Divide(emcal_ratio_mc)
		ihcal_var_detdeta_mc = TH1F(h_eT_var_ihcal.Clone("ihcal_var_detdeta_mc"))
		ihcal_var_detdeta_mc.Divide(ihcal_ratio_mc)
		ohcal_var_detdeta_mc = TH1F(h_eT_var_ohcal.Clone("ohcal_var_detdeta_mc"))
		ohcal_var_detdeta_mc.Divide(ohcal_ratio_mc)
		calo_var_detdeta_mc = TH1F(h_eT_var_calo.Clone("calo_var_detdeta_mc"))
		calo_var_detdeta_mc.Divide(calo_ratio_mc)
		hcal_var_detdeta_mc = TH1F(h_eT_var_hcal.Clone("hcal_var_detdeta_mc"))
		hcal_var_detdeta_mc.Divide(hcal_ratio_mc)

		emcal_detdeta_dev = TH1F(emcal_detdeta_mc.Clone("emcal_detdeta_dev"))
		emcal_detdeta_dev.SetXTitle("#eta")
		emcal_detdeta_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
		for i in range(1, emcal_detdeta_dev.GetNbinsX() + 1):
			emcal_detdeta_dev.SetBinContent(i, np.abs(emcal_detdeta_mc.GetBinContent(i) - emcal_var_detdeta_mc.GetBinContent(i)))
			emcal_detdeta_dev.SetBinError(i, np.sqrt(emcal_detdeta_mc.GetBinError(i)**2 + emcal_var_detdeta_mc.GetBinError(i)**2))

		ihcal_detdeta_dev = TH1F(ihcal_detdeta_mc.Clone("ihcal_detdeta_dev"))
		ihcal_detdeta_dev.SetXTitle("#eta")
		ihcal_detdeta_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
		for i in range(1, ihcal_detdeta_dev.GetNbinsX() + 1):
			ihcal_detdeta_dev.SetBinContent(i, np.abs(ihcal_detdeta_mc.GetBinContent(i) - ihcal_var_detdeta_mc.GetBinContent(i)))
			ihcal_detdeta_dev.SetBinError(i, np.sqrt(ihcal_detdeta_mc.GetBinError(i)**2 + ihcal_var_detdeta_mc.GetBinError(i)**2))

		ohcal_detdeta_dev = TH1F(ohcal_detdeta_mc.Clone("ohcal_detdeta_dev"))
		ohcal_detdeta_dev.SetXTitle("#eta")
		ohcal_detdeta_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
		for i in range(1, ohcal_detdeta_dev.GetNbinsX() + 1):
			ohcal_detdeta_dev.SetBinContent(i, np.abs(ohcal_detdeta_mc.GetBinContent(i) - ohcal_var_detdeta_mc.GetBinContent(i)))
			ohcal_detdeta_dev.SetBinError(i, np.sqrt(ohcal_detdeta_mc.GetBinError(i)**2 + ohcal_var_detdeta_mc.GetBinError(i)**2))

		calo_detdeta_dev = TH1F(calo_detdeta_mc.Clone("calo_detdeta_dev"))
		calo_detdeta_dev.SetXTitle("#eta")
		calo_detdeta_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
		for i in range(1, calo_detdeta_dev.GetNbinsX() + 1):
			calo_detdeta_dev.SetBinContent(i, np.abs(calo_detdeta_mc.GetBinContent(i) - calo_var_detdeta_mc.GetBinContent(i)))
			calo_detdeta_dev.SetBinError(i, np.sqrt(calo_detdeta_mc.GetBinError(i)**2 + calo_var_detdeta_mc.GetBinError(i)**2))

		hcal_detdeta_dev = TH1F(hcal_detdeta_mc.Clone("hcal_detdeta_dev"))
		hcal_detdeta_dev.SetXTitle("#eta")
		hcal_detdeta_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
		for i in range(1, hcal_detdeta_dev.GetNbinsX() + 1):
			hcal_detdeta_dev.SetBinContent(i, np.abs(hcal_detdeta_mc.GetBinContent(i) - hcal_var_detdeta_mc.GetBinContent(i)))
			hcal_detdeta_dev.SetBinError(i, np.sqrt(hcal_detdeta_mc.GetBinError(i)**2 + hcal_var_detdeta_mc.GetBinError(i)**2))

		emcal_reco_dev.Write()
		ihcal_reco_dev.Write()
		ohcal_reco_dev.Write()
		calo_reco_dev.Write()
		emcal_detdeta_dev.Write()
		ihcal_detdeta_dev.Write()
		ohcal_detdeta_dev.Write()
		calo_detdeta_dev.Write()
		hcal_detdeta_dev.Write()
		file.Close()
