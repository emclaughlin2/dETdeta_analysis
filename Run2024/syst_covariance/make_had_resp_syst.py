import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import numpy as np
import pdb
import sys

had_resp_unc = 0.117
emcal_had_frac = 0.61
ihcal_had_frac = 0.961
ohcal_had_frac = 0.972
calo_had_frac = 0.694
hcal_had_frac = 0.97

cents = ['0-5','5-10','10-20','20-30','30-40','40-50','50-60']
tag = 'had_resp'
for cent in cents:
    mcfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_run14_with_npart_nozs_mc_reweight_'+cent+'_reweight_epos_2024.root'
    datafile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_ana450_2024p009_100_50_50_ZS_hcal_tsc_emcal_calib_iter15_nozs_data_noweight_'+cent+'.root'
    outfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_variation_'+tag+'_'+cent+'.root'

    f1 = ROOT.TFile.Open(mcfile)
    h_emcal_correction = TH1F(f1.Get("h_emcal_correction"))
    h_ihcal_correction = TH1F(f1.Get("h_ihcal_correction"))
    h_ohcal_correction = TH1F(f1.Get("h_ohcal_correction"))
    h_calo_correction = TH1F(f1.Get("h_calo_correction"))
    h_hcal_correction = TH1F(f1.Get("h_hcal_correction"))
    h_emcal_correction.SetDirectory(0)
    h_ihcal_correction.SetDirectory(0)
    h_ohcal_correction.SetDirectory(0)
    h_calo_correction.SetDirectory(0)
    h_hcal_correction.SetDirectory(0)
    f1.Close()
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
    
    file = ROOT.TFile(outfile, "RECREATE")

    emcal_detdeta = TH1F(h_eT_data_emcal.Clone("emcal_detdeta"))
    emcal_detdeta.Divide(h_emcal_correction)
    ihcal_detdeta = TH1F(h_eT_data_ihcal.Clone("ihcal_detdeta"))
    ihcal_detdeta.Divide(h_ihcal_correction)
    ohcal_detdeta = TH1F(h_eT_data_ohcal.Clone("ohcal_detdeta"))
    ohcal_detdeta.Divide(h_ohcal_correction)
    calo_detdeta = TH1F(h_eT_data_calo.Clone("calo_detdeta"))
    calo_detdeta.Divide(h_calo_correction)
    hcal_detdeta = TH1F(h_eT_data_hcal.Clone("hcal_detdeta"))
    hcal_detdeta.Divide(h_hcal_correction) 

    emcal_detdeta_mean = emcal_detdeta.Clone("emcal_detdeta_mean")
    bins = emcal_detdeta_mean.GetNbinsX()
    emcal_detdeta_mean.Rebin(bins)
    emcal_detdeta_mean.Scale(1.0/bins)

    ihcal_detdeta_mean = ihcal_detdeta.Clone("ihcal_detdeta_mean")
    bins = ihcal_detdeta_mean.GetNbinsX()
    ihcal_detdeta_mean.Rebin(bins)
    ihcal_detdeta_mean.Scale(1.0/bins)

    ohcal_detdeta_mean = ohcal_detdeta.Clone("ohcal_detdeta_mean")
    bins = ohcal_detdeta_mean.GetNbinsX()
    ohcal_detdeta_mean.Rebin(bins)
    ohcal_detdeta_mean.Scale(1.0/bins)

    calo_detdeta_mean = calo_detdeta.Clone("calo_detdeta_mean")
    bins = calo_detdeta_mean.GetNbinsX()
    calo_detdeta_mean.Rebin(bins)
    calo_detdeta_mean.Scale(1.0/bins)

    hcal_detdeta_mean = hcal_detdeta.Clone("hcal_detdeta_mean")
    bins = hcal_detdeta_mean.GetNbinsX()
    hcal_detdeta_mean.Rebin(bins)
    hcal_detdeta_mean.Scale(1.0/bins)

    emcal_detdeta_dev = TH1F(emcal_detdeta.Clone("emcal_detdeta_dev"))
    emcal_detdeta_dev.SetXTitle("#eta")
    emcal_detdeta_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, emcal_detdeta_dev.GetNbinsX() + 1):
        emcal_detdeta_dev.SetBinContent(i, (emcal_detdeta_mean.GetBinContent(1)*2.0*had_resp_unc*emcal_had_frac)/(np.sqrt(12))) # RMS of uniform distribution
        emcal_detdeta_dev.SetBinError(i, (emcal_detdeta_mean.GetBinError(1)*2.0*had_resp_unc*emcal_had_frac)/(np.sqrt(12)))

    ihcal_detdeta_dev = TH1F(ihcal_detdeta.Clone("ihcal_detdeta_dev"))
    ihcal_detdeta_dev.SetXTitle("#eta")
    ihcal_detdeta_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, ihcal_detdeta_dev.GetNbinsX() + 1):
        ihcal_detdeta_dev.SetBinContent(i, (ihcal_detdeta_mean.GetBinContent(1)*2.0*had_resp_unc*ihcal_had_frac)/(np.sqrt(12)))
        ihcal_detdeta_dev.SetBinError(i, (ihcal_detdeta_mean.GetBinError(1)*2.0*had_resp_unc*ihcal_had_frac)/(np.sqrt(12)))

    ohcal_detdeta_dev = TH1F(ohcal_detdeta.Clone("ohcal_detdeta_dev"))
    ohcal_detdeta_dev.SetXTitle("#eta")
    ohcal_detdeta_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, ohcal_detdeta_dev.GetNbinsX() + 1):
        ohcal_detdeta_dev.SetBinContent(i, (ohcal_detdeta_mean.GetBinContent(1)*2.0*had_resp_unc*ohcal_had_frac)/(np.sqrt(12)))
        ohcal_detdeta_dev.SetBinError(i, (ohcal_detdeta_mean.GetBinError(1)*2.0*had_resp_unc*ohcal_had_frac)/(np.sqrt(12)))

    calo_detdeta_dev = TH1F(calo_detdeta.Clone("calo_detdeta_dev"))
    calo_detdeta_dev.SetXTitle("#eta")
    calo_detdeta_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, calo_detdeta_dev.GetNbinsX() + 1):
        calo_detdeta_dev.SetBinContent(i, (calo_detdeta_mean.GetBinContent(1)*2.0*had_resp_unc*calo_had_frac)/(np.sqrt(12)))
        calo_detdeta_dev.SetBinError(i, (calo_detdeta_mean.GetBinError(1)*2.0*had_resp_unc*calo_had_frac)/(np.sqrt(12)))

    hcal_detdeta_dev = TH1F(hcal_detdeta.Clone("hcal_detdeta_dev"))
    hcal_detdeta_dev.SetXTitle("#eta")
    hcal_detdeta_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, hcal_detdeta_dev.GetNbinsX() + 1):
        hcal_detdeta_dev.SetBinContent(i, (hcal_detdeta_mean.GetBinContent(1)*2.0*had_resp_unc*hcal_had_frac)/(np.sqrt(12)))
        hcal_detdeta_dev.SetBinError(i, (hcal_detdeta_mean.GetBinError(1)*2.0*had_resp_unc*hcal_had_frac)/(np.sqrt(12)))

    emcal_detdeta_dev.Write()
    ihcal_detdeta_dev.Write()
    ohcal_detdeta_dev.Write()
    calo_detdeta_dev.Write()
    hcal_detdeta_dev.Write()
    file.Close()