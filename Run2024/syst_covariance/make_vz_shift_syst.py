import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import numpy as np
import pdb
import sys

cents = ['0-5','5-10','10-20','20-30','30-40','40-50','50-60']
tag = 'vz_-3cm'
outfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/syst_covariance/dETdeta_signed_variation_'+tag+'.root'
emcal_detdeta_dev = TH1F("emcal_detdeta_dev","",7,0,7)
hcal_detdeta_dev = TH1F("hcal_detdeta_dev","",7,0,7)
calo_detdeta_dev = TH1F("calo_detdeta_dev","",7,0,7)
for c, cent in enumerate(cents):
	mcfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_run14_with_npart_nozs_mc_reweight_'+cent+'_reweight_epos_2024.root'
	datafile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_ana450_2024p009_100_50_50_ZS_hcal_tsc_emcal_calib_iter15_nozs_data_noweight_'+cent+'.root'
	vardatafile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_ana450_2024p009_100_50_50_ZS_hcal_tsc_emcal_calib_iter15,_nozs_data_noweight_'+cent+'_'+tag+'.root'

	f1 = ROOT.TFile.Open(mcfile)
	emcal_correction = TH1F(f1.Get("h_emcal_correction"))
	calo_correction = TH1F(f1.Get("h_calo_correction"))
	hcal_correction = TH1F(f1.Get("h_hcal_correction"))
	emcal_correction.SetDirectory(0)
	calo_correction.SetDirectory(0)
	hcal_correction.SetDirectory(0)
	f1.Close()
	f2 = ROOT.TFile.Open(datafile)
	h_eT_data_emcal = TH1F(f2.Get("h_eT_eta_emcal_profile_hist"))
	h_eT_data_calo = TH1F(f2.Get("h_eT_eta_calo_profile_hist"))
	h_eT_data_hcal = TH1F(f2.Get("h_eT_eta_hcal_profile_hist"))
	h_eT_data_emcal.SetDirectory(0)
	h_eT_data_calo.SetDirectory(0)
	h_eT_data_hcal.SetDirectory(0)
	f2.Close()
	f5 = ROOT.TFile.Open(vardatafile)
	h_eT_var_emcal = TH1F(f5.Get("h_eT_eta_emcal_profile_hist"))
	h_eT_var_calo = TH1F(f5.Get("h_eT_eta_calo_profile_hist"))
	h_eT_var_hcal = TH1F(f5.Get("h_eT_eta_hcal_profile_hist"))
	h_eT_var_emcal.SetDirectory(0)
	h_eT_var_calo.SetDirectory(0)
	h_eT_var_hcal.SetDirectory(0)
	f5.Close()

	# fully corrected histograms for nominal data
	emcal_detdeta_mc = TH1F(h_eT_data_emcal.Clone("emcal_detdeta_mc"))
	emcal_detdeta_mc.Divide(emcal_correction)
	calo_detdeta_mc = TH1F(h_eT_data_calo.Clone("calo_detdeta_mc"))
	calo_detdeta_mc.Divide(calo_correction)
	hcal_detdeta_mc = TH1F(h_eT_data_hcal.Clone("hcal_detdeta_mc"))
	hcal_detdeta_mc.Divide(hcal_correction) 

	# fully corrected histograms for variation data
	emcal_var_detdeta_mc = TH1F(h_eT_var_emcal.Clone("emcal_var_detdeta_mc"))
	emcal_var_detdeta_mc.Divide(emcal_correction)
	calo_var_detdeta_mc = TH1F(h_eT_var_calo.Clone("calo_var_detdeta_mc"))
	calo_var_detdeta_mc.Divide(calo_correction)
	hcal_var_detdeta_mc = TH1F(h_eT_var_hcal.Clone("hcal_var_detdeta_mc"))
	hcal_var_detdeta_mc.Divide(hcal_correction)

	emcal_detdeta_dev12 = TH1F(emcal_detdeta_mc.Clone("emcal_detdeta_dev12"))
	emcal_detdeta_dev12.SetXTitle("#eta")
	emcal_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
	for i in range(1, emcal_detdeta_dev12.GetNbinsX() + 1):
		emcal_detdeta_dev12.SetBinContent(i, emcal_detdeta_mc.GetBinContent(i) - emcal_var_detdeta_mc.GetBinContent(i))
		emcal_detdeta_dev12.SetBinError(i, np.sqrt(emcal_detdeta_mc.GetBinError(i)**2 + emcal_var_detdeta_mc.GetBinError(i)**2))

	calo_detdeta_dev12 = TH1F(calo_detdeta_mc.Clone("calo_detdeta_dev12"))
	calo_detdeta_dev12.SetXTitle("#eta")
	calo_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
	for i in range(1, calo_detdeta_dev12.GetNbinsX() + 1):
		calo_detdeta_dev12.SetBinContent(i, calo_detdeta_mc.GetBinContent(i) - calo_var_detdeta_mc.GetBinContent(i))
		calo_detdeta_dev12.SetBinError(i, np.sqrt(calo_detdeta_mc.GetBinError(i)**2 + calo_var_detdeta_mc.GetBinError(i)**2))

	hcal_detdeta_dev12 = TH1F(hcal_detdeta_mc.Clone("hcal_detdeta_dev12"))
	hcal_detdeta_dev12.SetXTitle("#eta")
	hcal_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
	for i in range(1, hcal_detdeta_dev12.GetNbinsX() + 1):
		hcal_detdeta_dev12.SetBinContent(i, hcal_detdeta_mc.GetBinContent(i) - hcal_var_detdeta_mc.GetBinContent(i))
		hcal_detdeta_dev12.SetBinError(i, np.sqrt(hcal_detdeta_mc.GetBinError(i)**2 + hcal_var_detdeta_mc.GetBinError(i)**2))

	bins = emcal_detdeta_dev12.GetNbinsX()
	emcal_detdeta_dev12.Rebin(bins)
	emcal_detdeta_dev12.Scale(1.0/bins)
	emcal_detdeta_dev.SetBinContent(c+1, emcal_detdeta_dev12.GetBinContent(1))
	
	bins = calo_detdeta_dev12.GetNbinsX()
	calo_detdeta_dev12.Rebin(bins)
	calo_detdeta_dev12.Scale(1.0/bins)
	calo_detdeta_dev.SetBinContent(c+1, calo_detdeta_dev12.GetBinContent(1))

	bins = hcal_detdeta_dev12.GetNbinsX()
	hcal_detdeta_dev12.Rebin(bins)
	hcal_detdeta_dev12.Scale(1.0/bins)
	hcal_detdeta_dev.SetBinContent(c+1, hcal_detdeta_dev12.GetBinContent(1))
	
file = ROOT.TFile(outfile, "RECREATE")
file.cd()
emcal_detdeta_dev.Write()
calo_detdeta_dev.Write()
hcal_detdeta_dev.Write()
file.Write()
file.Close()