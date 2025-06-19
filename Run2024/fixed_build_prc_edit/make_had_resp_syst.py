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

rgb = [[0, 0, 0], [230, 25, 75], [60, 180, 75], [255, 225, 25], [0, 130, 200], [245, 130, 48], [145, 30, 180], [70, 240, 240], [240, 50, 230], [210, 245, 60], [250, 190, 212], [0, 128, 128], [220, 190, 255], [170, 110, 40], [255, 250, 200], [128, 0, 0], [170, 255, 195], [128, 128, 0], [255, 215, 180], [0, 0, 128], [128, 128, 128], [34, 139, 34]]
colors = [TColor.GetColor(rgb[i][0],rgb[i][1],rgb[i][2]) for i in range(len(rgb))]

hadfracfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/dETdeta_em_had_energy_fraction_reweight_epos.root'
f = ROOT.TFile.Open(hadfracfile)
h_had_frac_emcal = TH1F(f.Get("h_had_frac_emcal"))
h_had_frac_ihcal = TH1F(f.Get("h_had_frac_ihcal"))
h_had_frac_ohcal = TH1F(f.Get("h_had_frac_ohcal"))
h_had_frac_calo = TH1F(f.Get("h_had_frac_calo"))
h_had_frac_hcal = TH1F(f.Get("h_had_frac_hcal"))
h_had_frac_emcal.SetDirectory(0)  
h_had_frac_ihcal.SetDirectory(0)  
h_had_frac_ohcal.SetDirectory(0)  
h_had_frac_calo.SetDirectory(0)   
h_had_frac_hcal.SetDirectory(0)
f.Close()

emcal_canvas = TCanvas("emcal_canvas","",600,500)
leg = ROOT.TLegend(.2,.75,.5,.89)
leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
leg.AddEntry("","Reweighted EPOS","")
emcal_canvas.SetLeftMargin(0.18)
leg.SetTextSize(0.035)
leg.SetBorderSize(0)
emcal_mean = h_had_frac_emcal.GetMean()
leg.AddEntry(h_had_frac_emcal,f"EMCal mean = {emcal_mean:.2f}","lep")
h_had_frac_emcal.SetStats(0)
h_had_frac_emcal.SetLineColor(colors[0])
h_had_frac_emcal.SetMarkerColor(colors[0])
h_had_frac_emcal.SetXTitle("Had Energy / Total Energy")
h_had_frac_emcal.SetYTitle("Events")
h_had_frac_emcal.GetYaxis().SetTitleOffset(1.8)
h_had_frac_emcal.Draw()
leg.Draw()
emcal_canvas.SaveAs('syst_plots/had_frac_emcal.png')

ihcal_canvas = TCanvas("ihcal_canvas","",600,500)
leg = ROOT.TLegend(.2,.75,.5,.89)
leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
leg.AddEntry("","Reweighted EPOS","")
ihcal_canvas.SetLeftMargin(0.18)
leg.SetTextSize(0.035)
leg.SetBorderSize(0)
ihcal_mean = h_had_frac_ihcal.GetMean()
leg.AddEntry(h_had_frac_ihcal,f"IHCal mean = {ihcal_mean:.2f}","lep")
h_had_frac_ihcal.SetStats(0)
h_had_frac_ihcal.SetLineColor(colors[0])
h_had_frac_ihcal.SetMarkerColor(colors[0])
h_had_frac_ihcal.SetXTitle("Had Energy / Total Energy")
h_had_frac_ihcal.SetYTitle("Events")
h_had_frac_ihcal.GetYaxis().SetTitleOffset(1.8)
h_had_frac_ihcal.Draw()
leg.Draw()
ihcal_canvas.SaveAs('syst_plots/had_frac_ihcal.png')

ohcal_canvas = TCanvas("ohcal_canvas","",600,500)
leg = ROOT.TLegend(.2,.75,.5,.89)
leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
leg.AddEntry("","Reweighted EPOS","")
ohcal_canvas.SetLeftMargin(0.18)
leg.SetTextSize(0.035)
leg.SetBorderSize(0)
ohcal_mean = h_had_frac_ohcal.GetMean()
leg.AddEntry(h_had_frac_ohcal,f"OHCal mean = {ohcal_mean:.2f}","lep")
h_had_frac_ohcal.SetStats(0)
h_had_frac_ohcal.SetLineColor(colors[0])
h_had_frac_ohcal.SetMarkerColor(colors[0])
h_had_frac_ohcal.SetXTitle("Had Energy / Total Energy")
h_had_frac_ohcal.SetYTitle("Events")
h_had_frac_ohcal.GetYaxis().SetTitleOffset(1.8)
h_had_frac_ohcal.Draw()
leg.Draw()
ohcal_canvas.SaveAs('syst_plots/had_frac_ohcal.png')

calo_canvas = TCanvas("calo_canvas","",600,500)
leg = ROOT.TLegend(.2,.75,.5,.89)
leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
leg.AddEntry("","Reweighted EPOS","")
calo_canvas.SetLeftMargin(0.18)
leg.SetTextSize(0.035)
leg.SetBorderSize(0)
calo_mean = h_had_frac_calo.GetMean()
leg.AddEntry(h_had_frac_calo,f"Full Calo mean = {calo_mean:.2f}","lep")
h_had_frac_calo.SetStats(0)
h_had_frac_calo.SetLineColor(colors[0])
h_had_frac_calo.SetMarkerColor(colors[0])
h_had_frac_calo.SetXTitle("Had Energy / Total Energy")
h_had_frac_calo.SetYTitle("Events")
h_had_frac_calo.GetYaxis().SetTitleOffset(1.8)
h_had_frac_calo.Draw()
leg.Draw()
calo_canvas.SaveAs('syst_plots/had_frac_calo.png')

hcal_canvas = TCanvas("hcal_canvas","",600,500)
leg = ROOT.TLegend(.2,.75,.5,.89)
leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
leg.AddEntry("","Reweighted EPOS","")
hcal_canvas.SetLeftMargin(0.18)
leg.SetTextSize(0.035)
leg.SetBorderSize(0)
hcal_mean = h_had_frac_hcal.GetMean()
leg.AddEntry(h_had_frac_hcal,f"HCal mean = {hcal_mean:.2f}","lep")
h_had_frac_hcal.SetStats(0)
h_had_frac_hcal.SetLineColor(colors[0])
h_had_frac_hcal.SetMarkerColor(colors[0])
h_had_frac_hcal.SetXTitle("Had Energy / Total Energy")
h_had_frac_hcal.SetYTitle("Events")
h_had_frac_hcal.GetYaxis().SetTitleOffset(1.8)
h_had_frac_hcal.Draw()
leg.Draw()
hcal_canvas.SaveAs('syst_plots/had_frac_hcal.png')


cents = ['0-5','5-10','10-20','20-30','30-40','40-50','50-60','60-70']
tag = 'had_resp'
for cent in cents:
    mcfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build_prc_edit/dETdeta_analysis_allruns_run14_with_centbin_nozs_mc_reweight_'+cent+'_reweight_epos_2024.root'
    datafile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build_prc_edit/dETdeta_analysis_allruns_ana450_2024p009_100_50_50_ZS_hcal_scaled_emcal_calib_iter26_nozs_data_noweight_'+cent+'.root'
    outfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build_prc_edit/dETdeta_variation_'+tag+'_'+cent+'.root'

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
