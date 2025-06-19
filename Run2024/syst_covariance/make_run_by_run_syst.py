import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import numpy as np
import pdb
import sys

cents = ['0-5','5-10','10-20','20-30','30-40','40-50','50-60']
tag = 'run_by_run'
outfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/syst_covariance/dETdeta_signed_variation_run_by_run.root'
emcal_detdeta_dev = TH1F("emcal_detdeta_dev","",7,0,7)
hcal_detdeta_dev = TH1F("hcal_detdeta_dev","",7,0,7)
calo_detdeta_dev = TH1F("calo_detdeta_dev","",7,0,7)
pos_eta_emcal_detdeta_dev = TH2F("pos_eta_emcal_detdeta_dev","",3,0,3,7,0,7)
pos_eta_hcal_detdeta_dev = TH2F("pos_eta_hcal_detdeta_dev","",3,0,3,7,0,7)
pos_eta_calo_detdeta_dev = TH2F("pos_eta_calo_detdeta_dev","",3,0,3,7,0,7)
neg_eta_emcal_detdeta_dev = TH2F("neg_eta_emcal_detdeta_dev","",3,0,3,7,0,7)
neg_eta_hcal_detdeta_dev = TH2F("neg_eta_hcal_detdeta_dev","",3,0,3,7,0,7)
neg_eta_calo_detdeta_dev = TH2F("neg_eta_calo_detdeta_dev","",3,0,3,7,0,7)
for c, cent in enumerate(cents):
    mcfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_run14_with_npart_nozs_mc_reweight_'+cent+'_reweight_epos_2024.root'
    mcfile1 = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/run_by_run_syst/dETdeta_analysis_allruns_run54911_acceptance_nozs_mc_reweight_'+cent+'_reweight_epos_2024.root'
    datafile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_ana450_2024p009_100_50_50_ZS_hcal_tsc_emcal_calib_iter15_nozs_data_noweight_'+cent+'.root'
    datafile1 = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/run_by_run_syst/dETdeta_analysis_allruns_ana450_2024p009_54911_nozs_data_noweight_'+cent+'.root'

    f1 = ROOT.TFile.Open(mcfile)
    h_emcal_correction1 = TH1F(f1.Get("h_emcal_correction"))
    h_calo_correction1 = TH1F(f1.Get("h_calo_correction"))
    h_hcal_correction1 = TH1F(f1.Get("h_hcal_correction"))
    h_emcal_correction1.SetDirectory(0)
    h_calo_correction1.SetDirectory(0)
    h_hcal_correction1.SetDirectory(0)
    f1.Close()

    f1 = ROOT.TFile.Open(mcfile1)
    h_emcal_correction2 = TH1F(f1.Get("h_emcal_correction"))
    h_calo_correction2 = TH1F(f1.Get("h_calo_correction"))
    h_hcal_correction2 = TH1F(f1.Get("h_hcal_correction"))
    h_emcal_correction2.SetDirectory(0)
    h_calo_correction2.SetDirectory(0)
    h_hcal_correction2.SetDirectory(0)
    f1.Close()
    
    f2 = ROOT.TFile.Open(datafile)
    h_eT_data_emcal1 = TH1F(f2.Get("h_eT_eta_emcal_profile_hist"))
    h_eT_data_calo1 = TH1F(f2.Get("h_eT_eta_calo_profile_hist"))
    h_eT_data_hcal1 = TH1F(f2.Get("h_eT_eta_hcal_profile_hist"))
    h_eT_data_emcal1.SetDirectory(0)
    h_eT_data_calo1.SetDirectory(0)
    h_eT_data_hcal1.SetDirectory(0)
    f2.Close()

    f3 = ROOT.TFile.Open(datafile1)
    h_eT_data_emcal2 = TH1F(f3.Get("h_eT_eta_emcal_profile_hist"))
    h_eT_data_calo2 = TH1F(f3.Get("h_eT_eta_calo_profile_hist"))
    h_eT_data_hcal2 = TH1F(f3.Get("h_eT_eta_hcal_profile_hist"))
    h_eT_data_emcal2.SetDirectory(0)
    h_eT_data_calo2.SetDirectory(0)
    h_eT_data_hcal2.SetDirectory(0)
    f3.Close()

    emcal_detdeta_data1 = TH1F(h_eT_data_emcal1.Clone("emcal_detdeta_data1"))
    emcal_detdeta_data1.Divide(h_emcal_correction1)
    calo_detdeta_data1 = TH1F(h_eT_data_calo1.Clone("calo_detdeta_data1"))
    calo_detdeta_data1.Divide(h_calo_correction1)
    hcal_detdeta_data1 = TH1F(h_eT_data_hcal1.Clone("hcal_detdeta_data1"))
    hcal_detdeta_data1.Divide(h_hcal_correction1) 

    emcal_detdeta_data2 = TH1F(h_eT_data_emcal2.Clone("emcal_detdeta_data2"))
    emcal_detdeta_data2.Divide(h_emcal_correction2)
    calo_detdeta_data2 = TH1F(h_eT_data_calo2.Clone("calo_detdeta_data2"))
    calo_detdeta_data2.Divide(h_calo_correction2)
    hcal_detdeta_data2 = TH1F(h_eT_data_hcal2.Clone("hcal_detdeta_data2"))
    hcal_detdeta_data2.Divide(h_hcal_correction2) 

    emcal_detdeta_dev12 = TH1F(emcal_detdeta_data1.Clone("emcal_detdeta_dev12"))
    emcal_detdeta_dev12.SetXTitle("#eta")
    emcal_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, emcal_detdeta_dev12.GetNbinsX() + 1):
        emcal_detdeta_dev12.SetBinContent(i, emcal_detdeta_data1.GetBinContent(i) - emcal_detdeta_data2.GetBinContent(i))
        emcal_detdeta_dev12.SetBinError(i, np.sqrt(emcal_detdeta_data1.GetBinError(i)**2 + emcal_detdeta_data2.GetBinError(i)**2))

    calo_detdeta_dev12 = TH1F(calo_detdeta_data1.Clone("calo_detdeta_dev12"))
    calo_detdeta_dev12.SetXTitle("#eta")
    calo_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, calo_detdeta_dev12.GetNbinsX() + 1):
        calo_detdeta_dev12.SetBinContent(i, calo_detdeta_data1.GetBinContent(i) - calo_detdeta_data2.GetBinContent(i))
        calo_detdeta_dev12.SetBinError(i, np.sqrt(calo_detdeta_data1.GetBinError(i)**2 + calo_detdeta_data2.GetBinError(i)**2))

    hcal_detdeta_dev12 = TH1F(hcal_detdeta_data1.Clone("hcal_detdeta_dev12"))
    hcal_detdeta_dev12.SetXTitle("#eta")
    hcal_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
    for i in range(1, hcal_detdeta_dev12.GetNbinsX() + 1):
        hcal_detdeta_dev12.SetBinContent(i, hcal_detdeta_data1.GetBinContent(i) - hcal_detdeta_data2.GetBinContent(i))
        hcal_detdeta_dev12.SetBinError(i, np.sqrt(hcal_detdeta_data1.GetBinError(i)**2 + hcal_detdeta_data2.GetBinError(i)**2))

    bins = emcal_detdeta_dev12.GetNbinsX()
    neg_eta_emcal_detdeta_dev.SetBinContent(neg_eta_emcal_detdeta_dev.GetBin(1, c+1), emcal_detdeta_dev12.GetBinContent(1))
    pos_eta_emcal_detdeta_dev.SetBinContent(pos_eta_emcal_detdeta_dev.GetBin(1, c+1), emcal_detdeta_dev12.GetBinContent(6))
    neg_eta_emcal_detdeta_dev.SetBinContent(neg_eta_emcal_detdeta_dev.GetBin(2, c+1), emcal_detdeta_dev12.GetBinContent(2))
    pos_eta_emcal_detdeta_dev.SetBinContent(pos_eta_emcal_detdeta_dev.GetBin(2, c+1), emcal_detdeta_dev12.GetBinContent(5))
    neg_eta_emcal_detdeta_dev.SetBinContent(neg_eta_emcal_detdeta_dev.GetBin(3, c+1), emcal_detdeta_dev12.GetBinContent(3))
    pos_eta_emcal_detdeta_dev.SetBinContent(pos_eta_emcal_detdeta_dev.GetBin(3, c+1), emcal_detdeta_dev12.GetBinContent(4))
    emcal_detdeta_dev12.Rebin(bins)
    emcal_detdeta_dev12.Scale(1.0/bins)
    emcal_detdeta_dev.SetBinContent(c+1, emcal_detdeta_dev.GetBinContent(c+1) + emcal_detdeta_dev12.GetBinContent(1))

    bins = calo_detdeta_dev12.GetNbinsX()
    neg_eta_calo_detdeta_dev.SetBinContent(neg_eta_calo_detdeta_dev.GetBin(1, c+1), calo_detdeta_dev12.GetBinContent(1))
    pos_eta_calo_detdeta_dev.SetBinContent(pos_eta_calo_detdeta_dev.GetBin(1, c+1), calo_detdeta_dev12.GetBinContent(6))
    neg_eta_calo_detdeta_dev.SetBinContent(neg_eta_calo_detdeta_dev.GetBin(2, c+1), calo_detdeta_dev12.GetBinContent(2))
    pos_eta_calo_detdeta_dev.SetBinContent(pos_eta_calo_detdeta_dev.GetBin(2, c+1), calo_detdeta_dev12.GetBinContent(5))
    neg_eta_calo_detdeta_dev.SetBinContent(neg_eta_calo_detdeta_dev.GetBin(3, c+1), calo_detdeta_dev12.GetBinContent(3))
    pos_eta_calo_detdeta_dev.SetBinContent(pos_eta_calo_detdeta_dev.GetBin(3, c+1), calo_detdeta_dev12.GetBinContent(4))
    calo_detdeta_dev12.Rebin(bins)
    calo_detdeta_dev12.Scale(1.0/bins)
    calo_detdeta_dev.SetBinContent(c+1, calo_detdeta_dev.GetBinContent(c+1) + calo_detdeta_dev12.GetBinContent(1))

    bins = hcal_detdeta_dev12.GetNbinsX()
    neg_eta_hcal_detdeta_dev.SetBinContent(neg_eta_hcal_detdeta_dev.GetBin(1, c+1), hcal_detdeta_dev12.GetBinContent(1))
    pos_eta_hcal_detdeta_dev.SetBinContent(pos_eta_hcal_detdeta_dev.GetBin(1, c+1), hcal_detdeta_dev12.GetBinContent(6))
    neg_eta_hcal_detdeta_dev.SetBinContent(neg_eta_hcal_detdeta_dev.GetBin(2, c+1), hcal_detdeta_dev12.GetBinContent(2))
    pos_eta_hcal_detdeta_dev.SetBinContent(pos_eta_hcal_detdeta_dev.GetBin(2, c+1), hcal_detdeta_dev12.GetBinContent(5))
    neg_eta_hcal_detdeta_dev.SetBinContent(neg_eta_hcal_detdeta_dev.GetBin(3, c+1), hcal_detdeta_dev12.GetBinContent(3))
    pos_eta_hcal_detdeta_dev.SetBinContent(pos_eta_hcal_detdeta_dev.GetBin(3, c+1), hcal_detdeta_dev12.GetBinContent(4))
    hcal_detdeta_dev12.Rebin(bins)
    hcal_detdeta_dev12.Scale(1.0/bins)
    hcal_detdeta_dev.SetBinContent(c+1, hcal_detdeta_dev.GetBinContent(c+1) + hcal_detdeta_dev12.GetBinContent(1))
    
file = ROOT.TFile(outfile, "RECREATE")
file.cd()
pos_eta_emcal_detdeta_dev.Write()
pos_eta_calo_detdeta_dev.Write()
pos_eta_hcal_detdeta_dev.Write()
neg_eta_emcal_detdeta_dev.Write()
neg_eta_calo_detdeta_dev.Write()
neg_eta_hcal_detdeta_dev.Write()
emcal_detdeta_dev.Write()
calo_detdeta_dev.Write()
hcal_detdeta_dev.Write()
file.Write()
file.Close()