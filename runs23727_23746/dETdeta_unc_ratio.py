import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import numpy as np
import pdb

hijingfile = 'dETdeta_analysis_23727_z=0_p011_nozs_mc_reweight_0-5_hijing.root'
datafile = 'dETdeta_analysis_23727_z=0_emcal_nosyst_nozs_data_noweight_0-5.root'
eposfile = 'dETdeta_analysis_23727_z=0_p011_nozs_mc_reweight_0-5_epos.root'
amptfile = 'dETdeta_analysis_23727_z=0_p011_nozs_mc_reweight_0-5_ampt.root'
vardatafile = 'dETdeta_analysis_23727_z=0_emcal_syst2_nozs_data_noweight_0-5.root'
outfile = 'dETdeta_variation_EMCAL_SYST2_ratio.root'

f1 = ROOT.TFile.Open(hijingfile)
h_eT_truth = TH1F(f1.Get("hetdeta_ihcalbin"))
h_eT_truth_ihcalbin = TH1F(f1.Get("hetdeta_ihcalbin"))
h_eT_truth_ohcalbin = TH1F(f1.Get("hetdeta_ohcalbin"))
h_eT_sim_emcal = TH1F(f1.Get("h_eT_emcal"))
h_eT_sim_ihcal = TH1F(f1.Get("h_eT_ihcal"))
h_eT_sim_ohcal = TH1F(f1.Get("h_eT_ohcal"))
h_eT_truth.SetDirectory(0)
h_eT_truth_ihcalbin.SetDirectory(0)
h_eT_truth_ohcalbin.SetDirectory(0)
h_eT_sim_emcal.SetDirectory(0)
h_eT_sim_ihcal.SetDirectory(0)
h_eT_sim_ohcal.SetDirectory(0)
f1.Close()
f2 = ROOT.TFile.Open(datafile)
h_eT_data_emcal = TH1F(f2.Get("h_eT_emcal"))
h_eT_data_ihcal = TH1F(f2.Get("h_eT_ihcal"))
h_eT_data_ohcal = TH1F(f2.Get("h_eT_ohcal"))
h_eT_data_emcal.SetDirectory(0)
h_eT_data_ihcal.SetDirectory(0)
h_eT_data_ohcal.SetDirectory(0)
f2.Close()
f3 = ROOT.TFile.Open(eposfile)
h_eT_truth_epos = TH1F(f3.Get("hetdeta_ihcalbin"))
h_eT_truth_ihcalbin_epos = TH1F(f3.Get("hetdeta_ihcalbin"))
h_eT_truth_ohcalbin_epos = TH1F(f3.Get("hetdeta_ohcalbin"))
h_eT_epos_emcal = TH1F(f3.Get("h_eT_emcal"))
h_eT_epos_ihcal = TH1F(f3.Get("h_eT_ihcal"))
h_eT_epos_ohcal = TH1F(f3.Get("h_eT_ohcal"))
h_eT_truth_epos.SetDirectory(0)
h_eT_truth_ihcalbin_epos.SetDirectory(0)
h_eT_truth_ohcalbin_epos.SetDirectory(0)
h_eT_epos_emcal.SetDirectory(0)
h_eT_epos_ihcal.SetDirectory(0)
h_eT_epos_ohcal.SetDirectory(0)
f3.Close()
f4 = ROOT.TFile.Open(amptfile)
h_eT_truth_ampt = TH1F(f4.Get("hetdeta_ihcalbin"))
h_eT_truth_ihcalbin_ampt = TH1F(f4.Get("hetdeta_ihcalbin"))
h_eT_truth_ohcalbin_ampt = TH1F(f4.Get("hetdeta_ohcalbin"))
h_eT_ampt_emcal = TH1F(f4.Get("h_eT_emcal"))
h_eT_ampt_ihcal = TH1F(f4.Get("h_eT_ihcal"))
h_eT_ampt_ohcal = TH1F(f4.Get("h_eT_ohcal"))
h_eT_truth_ampt.SetDirectory(0)
h_eT_truth_ihcalbin_ampt.SetDirectory(0)
h_eT_truth_ohcalbin_ampt.SetDirectory(0)
h_eT_ampt_emcal.SetDirectory(0)
h_eT_ampt_ihcal.SetDirectory(0)
h_eT_ampt_ohcal.SetDirectory(0)
f4.Close()
f5 = ROOT.TFile.Open(vardatafile)
h_eT_var_emcal = TH1F(f5.Get("h_eT_emcal"))
h_eT_var_ihcal = TH1F(f5.Get("h_eT_ihcal"))
h_eT_var_ohcal = TH1F(f5.Get("h_eT_ohcal"))
h_eT_var_emcal.SetDirectory(0)
h_eT_var_ihcal.SetDirectory(0)
h_eT_var_ohcal.SetDirectory(0)
f5.Close()

file = ROOT.TFile(outfile, "RECREATE")

# find the deviation from nominal 
emcal_reco_dev = TH1F(h_eT_data_emcal.Clone("emcal_reco_dev"))
emcal_reco_dev.SetXTitle("#eta")
emcal_reco_dev.SetYTitle("Reco dE_{T}/d#eta dev. [GeV]")
for i in range(1, emcal_reco_dev.GetNbinsX() + 1):
	emcal_reco_dev.SetBinContent(i, h_eT_data_emcal.GetBinContent(i) - h_eT_var_emcal.GetBinContent(i))
	emcal_reco_dev.SetBinError(i, np.sqrt(h_eT_data_emcal.GetBinError(i)**2 + h_eT_var_emcal.GetBinError(i)**2))
ihcal_reco_dev = TH1F(h_eT_data_ihcal.Clone("ihcal_reco_dev"))
ihcal_reco_dev.SetXTitle("#eta")
ihcal_reco_dev.SetYTitle("Reco dE_{T}/d#eta dev. [GeV]")
for i in range(1, ihcal_reco_dev.GetNbinsX() + 1):
	ihcal_reco_dev.SetBinContent(i, h_eT_data_ihcal.GetBinContent(i) - h_eT_var_ihcal.GetBinContent(i))
	ihcal_reco_dev.SetBinError(i, np.sqrt(h_eT_data_ihcal.GetBinError(i)**2 + h_eT_var_ihcal.GetBinError(i)**2))
ohcal_reco_dev = TH1F(h_eT_data_ohcal.Clone("ohcal_reco_dev"))
ohcal_reco_dev.SetXTitle("#eta")
ohcal_reco_dev.SetYTitle("Reco dE_{T}/d#eta dev. [GeV]")
for i in range(1, ohcal_reco_dev.GetNbinsX() + 1):
	ohcal_reco_dev.SetBinContent(i, h_eT_data_ohcal.GetBinContent(i) - h_eT_var_ohcal.GetBinContent(i))
	ohcal_reco_dev.SetBinError(i, np.sqrt(h_eT_data_ohcal.GetBinError(i)**2 + h_eT_var_ohcal.GetBinError(i)**2))

# correction factors from MC
emcal_ratio_hijing = TH1F(h_eT_sim_emcal.Clone("emcal_ratio_hijing"))
emcal_ratio_hijing.Divide(h_eT_truth)
emcal_ratio_epos = TH1F(h_eT_epos_emcal.Clone("emcal_ratio_epos"))
emcal_ratio_epos.Divide(h_eT_truth_epos)
emcal_ratio_ampt = TH1F(h_eT_ampt_emcal.Clone("emcal_ratio_ampt"))
emcal_ratio_ampt.Divide(h_eT_truth_ampt)

ihcal_ratio_hijing = TH1F(h_eT_sim_ihcal.Clone("ihcal_ratio_hijing"))
ihcal_ratio_hijing.Divide(h_eT_truth_ihcalbin)
ihcal_ratio_epos = TH1F(h_eT_epos_ihcal.Clone("ihcal_ratio_epos"))
ihcal_ratio_epos.Divide(h_eT_truth_ihcalbin_epos)
ihcal_ratio_ampt = TH1F(h_eT_ampt_ihcal.Clone("ihcal_ratio_ampt"))
ihcal_ratio_ampt.Divide(h_eT_truth_ihcalbin_ampt)

ohcal_ratio_hijing = TH1F(h_eT_sim_ohcal.Clone("ohcal_ratio_hijing"))
ohcal_ratio_hijing.Divide(h_eT_truth_ohcalbin)
ohcal_ratio_epos = TH1F(h_eT_epos_ohcal.Clone("ohcal_ratio_epos"))
ohcal_ratio_epos.Divide(h_eT_truth_ohcalbin_epos)
ohcal_ratio_ampt = TH1F(h_eT_ampt_ohcal.Clone("ohcal_ratio_ampt"))
ohcal_ratio_ampt.Divide(h_eT_truth_ohcalbin_ampt)

# fully corrected histograms for nominal data
emcal_detdeta_hijing = TH1F(h_eT_data_emcal.Clone("emcal_detdeta_hijing"))
emcal_detdeta_hijing.Divide(emcal_ratio_hijing)
emcal_detdeta_epos = TH1F(h_eT_data_emcal.Clone("emcal_detdeta_epos"))
emcal_detdeta_epos.Divide(emcal_ratio_epos)
emcal_detdeta_ampt = TH1F(h_eT_data_emcal.Clone("emcal_detdeta_ampt"))
emcal_detdeta_ampt.Divide(emcal_ratio_ampt)

ihcal_detdeta_hijing = TH1F(h_eT_data_ihcal.Clone("ihcal_detdeta_hijing"))
ihcal_detdeta_hijing.Divide(ihcal_ratio_hijing)
ihcal_detdeta_epos = TH1F(h_eT_data_ihcal.Clone("ihcal_detdeta_epos"))
ihcal_detdeta_epos.Divide(ihcal_ratio_epos)
ihcal_detdeta_ampt = TH1F(h_eT_data_ihcal.Clone("ihcal_detdeta_ampt"))
ihcal_detdeta_ampt.Divide(ihcal_ratio_ampt)

ohcal_detdeta_hijing = TH1F(h_eT_data_ohcal.Clone("ohcal_detdeta_hijing"))
ohcal_detdeta_hijing.Divide(ohcal_ratio_hijing)
ohcal_detdeta_epos = TH1F(h_eT_data_ohcal.Clone("ohcal_detdeta_epos"))
ohcal_detdeta_epos.Divide(ohcal_ratio_epos)
ohcal_detdeta_ampt = TH1F(h_eT_data_ohcal.Clone("ohcal_detdeta_ampt"))
ohcal_detdeta_ampt.Divide(ohcal_ratio_ampt)

# fully corrected histograms for variation data
emcal_var_detdeta_hijing = TH1F(h_eT_var_emcal.Clone("emcal_var_detdeta_hijing"))
emcal_var_detdeta_hijing.Divide(emcal_ratio_hijing)
emcal_var_detdeta_epos = TH1F(h_eT_var_emcal.Clone("emcal_var_detdeta_epos"))
emcal_var_detdeta_epos.Divide(emcal_ratio_epos)
emcal_var_detdeta_ampt = TH1F(h_eT_var_emcal.Clone("emcal_var_detdeta_ampt"))
emcal_var_detdeta_ampt.Divide(emcal_ratio_ampt)

ihcal_var_detdeta_hijing = TH1F(h_eT_var_ihcal.Clone("ihcal_var_detdeta_hijing"))
ihcal_var_detdeta_hijing.Divide(ihcal_ratio_hijing)
ihcal_var_detdeta_epos = TH1F(h_eT_var_ihcal.Clone("ihcal_var_detdeta_epos"))
ihcal_var_detdeta_epos.Divide(ihcal_ratio_epos)
ihcal_var_detdeta_ampt = TH1F(h_eT_var_ihcal.Clone("ihcal_var_detdeta_ampt"))
ihcal_var_detdeta_ampt.Divide(ihcal_ratio_ampt)

ohcal_var_detdeta_hijing = TH1F(h_eT_var_ohcal.Clone("ohcal_var_detdeta_hijing"))
ohcal_var_detdeta_hijing.Divide(ohcal_ratio_hijing)
ohcal_var_detdeta_epos = TH1F(h_eT_var_ohcal.Clone("ohcal_var_detdeta_epos"))
ohcal_var_detdeta_epos.Divide(ohcal_ratio_epos)
ohcal_var_detdeta_ampt = TH1F(h_eT_var_ohcal.Clone("ohcal_var_detdeta_ampt"))
ohcal_var_detdeta_ampt.Divide(ohcal_ratio_ampt)

emcal_hijing_dev = TH1F(emcal_detdeta_hijing.Clone("emcal_hijing_dev"))
emcal_hijing_dev.SetXTitle("#eta")
emcal_hijing_dev.SetYTitle("dE_{T}/d#eta dev. ratio")
for i in range(1, emcal_hijing_dev.GetNbinsX() + 1):
	if (emcal_detdeta_hijing.GetBinContent(i) != 0):
		emcal_hijing_dev.SetBinContent(i, emcal_var_detdeta_hijing.GetBinContent(i)/emcal_detdeta_hijing.GetBinContent(i))
	#emcal_hijing_dev.SetBinError(i, np.sqrt(emcal_detdeta_hijing.GetBinError(i)**2 + emcal_var_detdeta_hijing.GetBinError(i)**2))
emcal_epos_dev = TH1F(emcal_detdeta_epos.Clone("emcal_epos_dev"))
emcal_epos_dev.SetXTitle("#eta")
emcal_epos_dev.SetYTitle("dE_{T}/d#eta dev. ratio")
for i in range(1, emcal_epos_dev.GetNbinsX() + 1):
	if (emcal_detdeta_epos.GetBinContent(i) != 0):
		emcal_epos_dev.SetBinContent(i, emcal_var_detdeta_epos.GetBinContent(i)/emcal_detdeta_epos.GetBinContent(i))
	#emcal_epos_dev.SetBinError(i, np.sqrt(emcal_detdeta_epos.GetBinError(i)**2 + emcal_var_detdeta_epos.GetBinError(i)**2))
emcal_ampt_dev = TH1F(emcal_detdeta_ampt.Clone("emcal_ampt_dev"))
emcal_ampt_dev.SetXTitle("#eta")
emcal_ampt_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
for i in range(1, emcal_ampt_dev.GetNbinsX() + 1):
	emcal_ampt_dev.SetBinContent(i, emcal_detdeta_ampt.GetBinContent(i) - emcal_var_detdeta_ampt.GetBinContent(i))
	emcal_ampt_dev.SetBinError(i, np.sqrt(emcal_detdeta_ampt.GetBinError(i)**2 + emcal_var_detdeta_ampt.GetBinError(i)**2))

ihcal_hijing_dev = TH1F(ihcal_detdeta_hijing.Clone("ihcal_hijing_dev"))
ihcal_hijing_dev.SetXTitle("#eta")
ihcal_hijing_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
for i in range(1, ihcal_hijing_dev.GetNbinsX() + 1):
	ihcal_hijing_dev.SetBinContent(i, ihcal_detdeta_hijing.GetBinContent(i) - ihcal_var_detdeta_hijing.GetBinContent(i))
	ihcal_hijing_dev.SetBinError(i, np.sqrt(ihcal_detdeta_hijing.GetBinError(i)**2 + ihcal_var_detdeta_hijing.GetBinError(i)**2))
ihcal_epos_dev = TH1F(ihcal_detdeta_epos.Clone("ihcal_epos_dev"))
ihcal_epos_dev.SetXTitle("#eta")
ihcal_epos_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
for i in range(1, ihcal_epos_dev.GetNbinsX() + 1):
	ihcal_epos_dev.SetBinContent(i, ihcal_detdeta_epos.GetBinContent(i) - ihcal_var_detdeta_epos.GetBinContent(i))
	ihcal_epos_dev.SetBinError(i, np.sqrt(ihcal_detdeta_epos.GetBinError(i)**2 + ihcal_var_detdeta_epos.GetBinError(i)**2))
ihcal_ampt_dev = TH1F(ihcal_detdeta_ampt.Clone("ihcal_ampt_dev"))
ihcal_ampt_dev.SetXTitle("#eta")
ihcal_ampt_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
for i in range(1, ihcal_ampt_dev.GetNbinsX() + 1):
	ihcal_ampt_dev.SetBinContent(i, ihcal_detdeta_ampt.GetBinContent(i) - ihcal_var_detdeta_ampt.GetBinContent(i))
	ihcal_ampt_dev.SetBinError(i, np.sqrt(ihcal_detdeta_ampt.GetBinError(i)**2 + ihcal_var_detdeta_ampt.GetBinError(i)**2))

ohcal_hijing_dev = TH1F(ohcal_detdeta_hijing.Clone("ohcal_hijing_dev"))
ohcal_hijing_dev.SetXTitle("#eta")
ohcal_hijing_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
for i in range(1, ohcal_hijing_dev.GetNbinsX() + 1):
	ohcal_hijing_dev.SetBinContent(i, ohcal_detdeta_hijing.GetBinContent(i) - ohcal_var_detdeta_hijing.GetBinContent(i))
	ohcal_hijing_dev.SetBinError(i, np.sqrt(ohcal_detdeta_hijing.GetBinError(i)**2 + ohcal_var_detdeta_hijing.GetBinError(i)**2))
ohcal_epos_dev = TH1F(ohcal_detdeta_epos.Clone("ohcal_epos_dev"))
ohcal_epos_dev.SetXTitle("#eta")
ohcal_epos_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
for i in range(1, ohcal_epos_dev.GetNbinsX() + 1):
	ohcal_epos_dev.SetBinContent(i, ohcal_detdeta_epos.GetBinContent(i) - ohcal_var_detdeta_epos.GetBinContent(i))
	ohcal_epos_dev.SetBinError(i, np.sqrt(ohcal_detdeta_epos.GetBinError(i)**2 + ohcal_var_detdeta_epos.GetBinError(i)**2))
ohcal_ampt_dev = TH1F(ohcal_detdeta_ampt.Clone("ohcal_ampt_dev"))
ohcal_ampt_dev.SetXTitle("#eta")
ohcal_ampt_dev.SetYTitle("dE_{T}/d#eta dev. [GeV]")
for i in range(1, ohcal_ampt_dev.GetNbinsX() + 1):
	ohcal_ampt_dev.SetBinContent(i, ohcal_detdeta_ampt.GetBinContent(i) - ohcal_var_detdeta_ampt.GetBinContent(i))
	ohcal_ampt_dev.SetBinError(i, np.sqrt(ohcal_detdeta_ampt.GetBinError(i)**2 + ohcal_var_detdeta_ampt.GetBinError(i)**2))

emcal_reco_dev.Write()
ihcal_reco_dev.Write()
ohcal_reco_dev.Write()
emcal_hijing_dev.Write()
emcal_epos_dev.Write()
emcal_ampt_dev.Write()
ihcal_hijing_dev.Write()
ihcal_epos_dev.Write()
ihcal_ampt_dev.Write()
ohcal_hijing_dev.Write()
ohcal_epos_dev.Write()
ohcal_ampt_dev.Write()
file.Close()
