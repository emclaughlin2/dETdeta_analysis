import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import numpy as np
import pdb
import sys

cents = ['0-5','5-10','10-20','20-30','30-40','40-50','50-60']
tag = 'MC'
for cent in cents:
    mcfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_run14_with_npart_mb_nozs_mc_reweight_'+cent+'_reweight_epos_2024.root'
    mcfile1 = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_run14_with_npart_mb_nozs_mc_reweight_'+cent+'_reweight_ampt_2024.root'
    mcfile2 = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_run14_with_npart_mb_nozs_mc_reweight_'+cent+'_reweight_hijing_2024.root'
    datafile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_ana450_2024p009_100_50_50_ZS_hcal_scaled_emcal_calib_iter26_nozs_data_noweight_'+cent+'.root'
    outfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_variation1_'+tag+'_'+cent+'.root'

    f1 = ROOT.TFile.Open(mcfile)
    h_emcal_correction1 = TH1F(f1.Get("h_emcal_correction"))
    h_ihcal_correction1 = TH1F(f1.Get("h_ihcal_correction"))
    h_ohcal_correction1 = TH1F(f1.Get("h_ohcal_correction"))
    h_calo_correction1 = TH1F(f1.Get("h_calo_correction"))
    h_hcal_correction1 = TH1F(f1.Get("h_hcal_correction"))
    h_emcal_correction1.SetDirectory(0)
    h_ihcal_correction1.SetDirectory(0)
    h_ohcal_correction1.SetDirectory(0)
    h_calo_correction1.SetDirectory(0)
    h_hcal_correction1.SetDirectory(0)
    f1.Close()
    f1 = ROOT.TFile.Open(mcfile1)
    h_emcal_correction2 = TH1F(f1.Get("h_emcal_correction"))
    h_ihcal_correction2 = TH1F(f1.Get("h_ihcal_correction"))
    h_ohcal_correction2 = TH1F(f1.Get("h_ohcal_correction"))
    h_calo_correction2 = TH1F(f1.Get("h_calo_correction"))
    h_hcal_correction2 = TH1F(f1.Get("h_hcal_correction"))
    h_emcal_correction2.SetDirectory(0)
    h_ihcal_correction2.SetDirectory(0)
    h_ohcal_correction2.SetDirectory(0)
    h_calo_correction2.SetDirectory(0)
    h_hcal_correction2.SetDirectory(0)
    f1.Close()
    f1 = ROOT.TFile.Open(mcfile2)
    h_emcal_correction3 = TH1F(f1.Get("h_emcal_correction"))
    h_ihcal_correction3 = TH1F(f1.Get("h_ihcal_correction"))
    h_ohcal_correction3 = TH1F(f1.Get("h_ohcal_correction"))
    h_calo_correction3 = TH1F(f1.Get("h_calo_correction"))
    h_hcal_correction3 = TH1F(f1.Get("h_hcal_correction"))
    h_emcal_correction3.SetDirectory(0)
    h_ihcal_correction3.SetDirectory(0)
    h_ohcal_correction3.SetDirectory(0)
    h_calo_correction3.SetDirectory(0)
    h_hcal_correction3.SetDirectory(0)
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

    emcal_detdeta_mc1 = TH1F(h_eT_data_emcal.Clone("emcal_detdeta_mc1"))
    emcal_detdeta_mc1.Divide(h_emcal_correction1)
    ihcal_detdeta_mc1 = TH1F(h_eT_data_ihcal.Clone("ihcal_detdeta_mc1"))
    ihcal_detdeta_mc1.Divide(h_ihcal_correction1)
    ohcal_detdeta_mc1 = TH1F(h_eT_data_ohcal.Clone("ohcal_detdeta_mc1"))
    ohcal_detdeta_mc1.Divide(h_ohcal_correction1)
    calo_detdeta_mc1 = TH1F(h_eT_data_calo.Clone("calo_detdeta_mc1"))
    calo_detdeta_mc1.Divide(h_calo_correction1)
    hcal_detdeta_mc1 = TH1F(h_eT_data_hcal.Clone("hcal_detdeta_mc1"))
    hcal_detdeta_mc1.Divide(h_hcal_correction1) 

    emcal_detdeta_mc2 = TH1F(h_eT_data_emcal.Clone("emcal_detdeta_mc2"))
    emcal_detdeta_mc2.Divide(h_emcal_correction2)
    ihcal_detdeta_mc2 = TH1F(h_eT_data_ihcal.Clone("ihcal_detdeta_mc2"))
    ihcal_detdeta_mc2.Divide(h_ihcal_correction2)
    ohcal_detdeta_mc2 = TH1F(h_eT_data_ohcal.Clone("ohcal_detdeta_mc2"))
    ohcal_detdeta_mc2.Divide(h_ohcal_correction2)
    calo_detdeta_mc2 = TH1F(h_eT_data_calo.Clone("calo_detdeta_mc2"))
    calo_detdeta_mc2.Divide(h_calo_correction2)
    hcal_detdeta_mc2 = TH1F(h_eT_data_hcal.Clone("hcal_detdeta_mc2"))
    hcal_detdeta_mc2.Divide(h_hcal_correction2) 

    emcal_detdeta_mc3 = TH1F(h_eT_data_emcal.Clone("emcal_detdeta_mc3"))
    emcal_detdeta_mc3.Divide(h_emcal_correction3)
    ihcal_detdeta_mc3 = TH1F(h_eT_data_ihcal.Clone("ihcal_detdeta_mc3"))
    ihcal_detdeta_mc3.Divide(h_ihcal_correction3)
    ohcal_detdeta_mc3 = TH1F(h_eT_data_ohcal.Clone("ohcal_detdeta_mc3"))
    ohcal_detdeta_mc3.Divide(h_ohcal_correction3)
    calo_detdeta_mc3 = TH1F(h_eT_data_calo.Clone("calo_detdeta_mc3"))
    calo_detdeta_mc3.Divide(h_calo_correction3)
    hcal_detdeta_mc3 = TH1F(h_eT_data_hcal.Clone("hcal_detdeta_mc3"))
    hcal_detdeta_mc3.Divide(h_hcal_correction3)

    emcal_detdeta_dev12 = TH1F(emcal_detdeta_mc1.Clone("emcal_detdeta_dev12"))
    emcal_detdeta_dev12.SetXTitle("#eta")
    emcal_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, emcal_detdeta_dev12.GetNbinsX() + 1):
        emcal_detdeta_dev12.SetBinContent(i, np.abs(emcal_detdeta_mc1.GetBinContent(i) - emcal_detdeta_mc2.GetBinContent(i)))
        emcal_detdeta_dev12.SetBinError(i, np.sqrt(emcal_detdeta_mc1.GetBinError(i)**2 + emcal_detdeta_mc2.GetBinError(i)**2))

    ihcal_detdeta_dev12 = TH1F(ihcal_detdeta_mc1.Clone("ihcal_detdeta_dev12"))
    ihcal_detdeta_dev12.SetXTitle("#eta")
    ihcal_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, ihcal_detdeta_dev12.GetNbinsX() + 1):
        ihcal_detdeta_dev12.SetBinContent(i, np.abs(ihcal_detdeta_mc1.GetBinContent(i) - ihcal_detdeta_mc2.GetBinContent(i)))
        ihcal_detdeta_dev12.SetBinError(i, np.sqrt(ihcal_detdeta_mc1.GetBinError(i)**2 + ihcal_detdeta_mc2.GetBinError(i)**2))

    ohcal_detdeta_dev12 = TH1F(ohcal_detdeta_mc1.Clone("ohcal_detdeta_dev12"))
    ohcal_detdeta_dev12.SetXTitle("#eta")
    ohcal_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, ohcal_detdeta_dev12.GetNbinsX() + 1):
        ohcal_detdeta_dev12.SetBinContent(i, np.abs(ohcal_detdeta_mc1.GetBinContent(i) - ohcal_detdeta_mc2.GetBinContent(i)))
        ohcal_detdeta_dev12.SetBinError(i, np.sqrt(ohcal_detdeta_mc1.GetBinError(i)**2 + ohcal_detdeta_mc2.GetBinError(i)**2))

    calo_detdeta_dev12 = TH1F(calo_detdeta_mc1.Clone("calo_detdeta_dev12"))
    calo_detdeta_dev12.SetXTitle("#eta")
    calo_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, calo_detdeta_dev12.GetNbinsX() + 1):
        calo_detdeta_dev12.SetBinContent(i, np.abs(calo_detdeta_mc1.GetBinContent(i) - calo_detdeta_mc2.GetBinContent(i)))
        calo_detdeta_dev12.SetBinError(i, np.sqrt(calo_detdeta_mc1.GetBinError(i)**2 + calo_detdeta_mc2.GetBinError(i)**2))

    hcal_detdeta_dev12 = TH1F(hcal_detdeta_mc1.Clone("hcal_detdeta_dev12"))
    hcal_detdeta_dev12.SetXTitle("#eta")
    hcal_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, hcal_detdeta_dev12.GetNbinsX() + 1):
        hcal_detdeta_dev12.SetBinContent(i, np.abs(hcal_detdeta_mc1.GetBinContent(i) - hcal_detdeta_mc2.GetBinContent(i)))
        hcal_detdeta_dev12.SetBinError(i, np.sqrt(hcal_detdeta_mc1.GetBinError(i)**2 + hcal_detdeta_mc2.GetBinError(i)**2))

    emcal_detdeta_dev13 = TH1F(emcal_detdeta_mc1.Clone("emcal_detdeta_dev13"))
    emcal_detdeta_dev13.SetXTitle("#eta")
    emcal_detdeta_dev13.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, emcal_detdeta_dev13.GetNbinsX() + 1):
        emcal_detdeta_dev13.SetBinContent(i, np.abs(emcal_detdeta_mc1.GetBinContent(i) - emcal_detdeta_mc3.GetBinContent(i)))
        emcal_detdeta_dev13.SetBinError(i, np.sqrt(emcal_detdeta_mc1.GetBinError(i)**2 + emcal_detdeta_mc3.GetBinError(i)**2))

    ihcal_detdeta_dev13 = TH1F(ihcal_detdeta_mc1.Clone("ihcal_detdeta_dev13"))
    ihcal_detdeta_dev13.SetXTitle("#eta")
    ihcal_detdeta_dev13.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, ihcal_detdeta_dev13.GetNbinsX() + 1):
        ihcal_detdeta_dev13.SetBinContent(i, np.abs(ihcal_detdeta_mc1.GetBinContent(i) - ihcal_detdeta_mc3.GetBinContent(i)))
        ihcal_detdeta_dev13.SetBinError(i, np.sqrt(ihcal_detdeta_mc1.GetBinError(i)**2 + ihcal_detdeta_mc3.GetBinError(i)**2))

    ohcal_detdeta_dev13 = TH1F(ohcal_detdeta_mc1.Clone("ohcal_detdeta_dev13"))
    ohcal_detdeta_dev13.SetXTitle("#eta")
    ohcal_detdeta_dev13.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, ohcal_detdeta_dev13.GetNbinsX() + 1):
        ohcal_detdeta_dev13.SetBinContent(i, np.abs(ohcal_detdeta_mc1.GetBinContent(i) - ohcal_detdeta_mc3.GetBinContent(i)))
        ohcal_detdeta_dev13.SetBinError(i, np.sqrt(ohcal_detdeta_mc1.GetBinError(i)**2 + ohcal_detdeta_mc3.GetBinError(i)**2))

    calo_detdeta_dev13 = TH1F(calo_detdeta_mc1.Clone("calo_detdeta_dev13"))
    calo_detdeta_dev13.SetXTitle("#eta")
    calo_detdeta_dev13.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, calo_detdeta_dev13.GetNbinsX() + 1):
        calo_detdeta_dev13.SetBinContent(i, np.abs(calo_detdeta_mc1.GetBinContent(i) - calo_detdeta_mc3.GetBinContent(i)))
        calo_detdeta_dev13.SetBinError(i, np.sqrt(calo_detdeta_mc1.GetBinError(i)**2 + calo_detdeta_mc3.GetBinError(i)**2))

    hcal_detdeta_dev13 = TH1F(hcal_detdeta_mc1.Clone("hcal_detdeta_dev13"))
    hcal_detdeta_dev13.SetXTitle("#eta")
    hcal_detdeta_dev13.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, hcal_detdeta_dev13.GetNbinsX() + 1):
        hcal_detdeta_dev13.SetBinContent(i, np.abs(hcal_detdeta_mc1.GetBinContent(i) - hcal_detdeta_mc3.GetBinContent(i)))
        hcal_detdeta_dev13.SetBinError(i, np.sqrt(hcal_detdeta_mc1.GetBinError(i)**2 + hcal_detdeta_mc3.GetBinError(i)**2))

    emcal_detdeta_dev = emcal_detdeta_dev12.Clone("emcal_detdeta_dev")
    for i in range(1, emcal_detdeta_dev12.GetNbinsX() + 1):
        emcal_detdeta_dev.SetBinContent(i, np.average([emcal_detdeta_dev12.GetBinContent(i), emcal_detdeta_dev13.GetBinContent(i)]))
        emcal_detdeta_dev.SetBinError(i, 0)
    
    ihcal_detdeta_dev = ihcal_detdeta_dev12.Clone("ihcal_detdeta_dev")
    for i in range(1, ihcal_detdeta_dev12.GetNbinsX() + 1):
        ihcal_detdeta_dev.SetBinContent(i, np.average([ihcal_detdeta_dev12.GetBinContent(i), ihcal_detdeta_dev13.GetBinContent(i)]))
        ihcal_detdeta_dev.SetBinError(i, 0)
    
    ohcal_detdeta_dev = ohcal_detdeta_dev12.Clone("ohcal_detdeta_dev")
    for i in range(1, ohcal_detdeta_dev12.GetNbinsX() + 1):
        ohcal_detdeta_dev.SetBinContent(i, np.average([ohcal_detdeta_dev12.GetBinContent(i), ohcal_detdeta_dev13.GetBinContent(i)]))
        ohcal_detdeta_dev.SetBinError(i, 0)
    
    calo_detdeta_dev = calo_detdeta_dev12.Clone("calo_detdeta_dev")
    for i in range(1, calo_detdeta_dev12.GetNbinsX() + 1):
        calo_detdeta_dev.SetBinContent(i, np.average([calo_detdeta_dev12.GetBinContent(i), calo_detdeta_dev13.GetBinContent(i)]))
        calo_detdeta_dev.SetBinError(i, 0)
        
    hcal_detdeta_dev = hcal_detdeta_dev12.Clone("hcal_detdeta_dev")
    for i in range(1, hcal_detdeta_dev12.GetNbinsX() + 1):
        hcal_detdeta_dev.SetBinContent(i, np.average([hcal_detdeta_dev12.GetBinContent(i), hcal_detdeta_dev13.GetBinContent(i)]))
        hcal_detdeta_dev.SetBinError(i, 0)
    
    emcal_detdeta_dev.Write()
    ihcal_detdeta_dev.Write()
    ohcal_detdeta_dev.Write()
    calo_detdeta_dev.Write()
    hcal_detdeta_dev.Write()
    file.Close()