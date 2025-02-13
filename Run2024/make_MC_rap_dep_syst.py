import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import numpy as np
import pdb
import sys

gROOT.LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
gROOT.ProcessLine("SetsPhenixStyle()")

rgb = [[230, 25, 75], [60, 180, 75], [255, 225, 25], [0, 130, 200], [245, 130, 48], [145, 30, 180], [70, 240, 240], [240, 50, 230], [210, 245, 60], [250, 190, 212], [0, 128, 128], [220, 190, 255], [170, 110, 40], [255, 250, 200], [128, 0, 0], [170, 255, 195], [128, 128, 0], [255, 215, 180], [0, 0, 128], [128, 128, 128], [34, 139, 34], [0, 0, 0]]
colors = [TColor.GetColor(rgb[i][0],rgb[i][1],rgb[i][2]) for i in range(len(rgb))]

mcsystfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_run14_with_npart_nozs_mc_reweight_0-5_reweight_epos_2024.root'
mcsystfile1 = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_run14_with_npart_nozs_mc_reweight_0-5_reweight_brahms_epos_2024.root'
datasystfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_ana450_2024p009_100_50_50_ZS_hcal_tsc_emcal_calib_iter15_nozs_data_noweight_0-5.root'

f1 = ROOT.TFile.Open(mcsystfile)
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
f1 = ROOT.TFile.Open(mcsystfile1)
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
f2 = ROOT.TFile.Open(datasystfile)
h_eT_data_emcal1 = TH1F(f2.Get("h_eT_eta_emcal_profile_hist"))
h_eT_data_ihcal1 = TH1F(f2.Get("h_eT_eta_ihcal_profile_hist"))
h_eT_data_ohcal1 = TH1F(f2.Get("h_eT_eta_ohcal_profile_hist"))
h_eT_data_calo1 = TH1F(f2.Get("h_eT_eta_calo_profile_hist"))
h_eT_data_hcal1 = TH1F(f2.Get("h_eT_eta_hcal_profile_hist"))
h_eT_data_emcal1.SetDirectory(0)
h_eT_data_ihcal1.SetDirectory(0)
h_eT_data_ohcal1.SetDirectory(0)
h_eT_data_calo1.SetDirectory(0)
h_eT_data_hcal1.SetDirectory(0)
f2.Close()

emcal_detdeta_mc1 = TH1F(h_eT_data_emcal1.Clone("emcal_detdeta_mc1"))
emcal_detdeta_mc1.Divide(h_emcal_correction1)
ihcal_detdeta_mc1 = TH1F(h_eT_data_ihcal1.Clone("ihcal_detdeta_mc1"))
ihcal_detdeta_mc1.Divide(h_ihcal_correction1)
ohcal_detdeta_mc1 = TH1F(h_eT_data_ohcal1.Clone("ohcal_detdeta_mc1"))
ohcal_detdeta_mc1.Divide(h_ohcal_correction1)
calo_detdeta_mc1 = TH1F(h_eT_data_calo1.Clone("calo_detdeta_mc1"))
calo_detdeta_mc1.Divide(h_calo_correction1)
hcal_detdeta_mc1 = TH1F(h_eT_data_hcal1.Clone("hcal_detdeta_mc1"))
hcal_detdeta_mc1.Divide(h_hcal_correction1) 

emcal_detdeta_mc2 = TH1F(h_eT_data_emcal1.Clone("emcal_detdeta_mc2"))
emcal_detdeta_mc2.Divide(h_emcal_correction2)
ihcal_detdeta_mc2 = TH1F(h_eT_data_ihcal1.Clone("ihcal_detdeta_mc2"))
ihcal_detdeta_mc2.Divide(h_ihcal_correction2)
ohcal_detdeta_mc2 = TH1F(h_eT_data_ohcal1.Clone("ohcal_detdeta_mc2"))
ohcal_detdeta_mc2.Divide(h_ohcal_correction2)
calo_detdeta_mc2 = TH1F(h_eT_data_calo1.Clone("calo_detdeta_mc2"))
calo_detdeta_mc2.Divide(h_calo_correction2)
hcal_detdeta_mc2 = TH1F(h_eT_data_hcal1.Clone("hcal_detdeta_mc2"))
hcal_detdeta_mc2.Divide(h_hcal_correction2) 

emcal_detdeta_dev12 = TH1F(emcal_detdeta_mc1.Clone("emcal_detdeta_dev_percent"))
emcal_detdeta_dev12.SetDirectory(0)
emcal_detdeta_dev12.SetXTitle("#eta")
emcal_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
for i in range(1, emcal_detdeta_dev12.GetNbinsX() + 1):
    emcal_detdeta_dev12.SetBinContent(i, np.abs(emcal_detdeta_mc1.GetBinContent(i) - emcal_detdeta_mc2.GetBinContent(i)))
    emcal_detdeta_dev12.SetBinError(i, np.sqrt(emcal_detdeta_mc1.GetBinError(i)**2 + emcal_detdeta_mc2.GetBinError(i)**2))
emcal_detdeta_dev12.Divide(emcal_detdeta_mc1)

ihcal_detdeta_dev12 = TH1F(ihcal_detdeta_mc1.Clone("ihcal_detdeta_dev_percent"))
ihcal_detdeta_dev12.SetDirectory(0)
ihcal_detdeta_dev12.SetXTitle("#eta")
ihcal_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
for i in range(1, ihcal_detdeta_dev12.GetNbinsX() + 1):
    ihcal_detdeta_dev12.SetBinContent(i, np.abs(ihcal_detdeta_mc1.GetBinContent(i) - ihcal_detdeta_mc2.GetBinContent(i)))
    ihcal_detdeta_dev12.SetBinError(i, np.sqrt(ihcal_detdeta_mc1.GetBinError(i)**2 + ihcal_detdeta_mc2.GetBinError(i)**2))
ihcal_detdeta_dev12.Divide(ihcal_detdeta_mc1)

ohcal_detdeta_dev12 = TH1F(ohcal_detdeta_mc1.Clone("ohcal_detdeta_dev_percent"))
ohcal_detdeta_dev12.SetDirectory(0)
ohcal_detdeta_dev12.SetXTitle("#eta")
ohcal_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
for i in range(1, ohcal_detdeta_dev12.GetNbinsX() + 1):
    ohcal_detdeta_dev12.SetBinContent(i, np.abs(ohcal_detdeta_mc1.GetBinContent(i) - ohcal_detdeta_mc2.GetBinContent(i)))
    ohcal_detdeta_dev12.SetBinError(i, np.sqrt(ohcal_detdeta_mc1.GetBinError(i)**2 + ohcal_detdeta_mc2.GetBinError(i)**2))
ohcal_detdeta_dev12.Divide(ohcal_detdeta_mc1)

calo_detdeta_dev12 = TH1F(calo_detdeta_mc1.Clone("calo_detdeta_dev_percent"))
calo_detdeta_dev12.SetDirectory(0)
calo_detdeta_dev12.SetXTitle("#eta")
calo_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
for i in range(1, calo_detdeta_dev12.GetNbinsX() + 1):
    calo_detdeta_dev12.SetBinContent(i, np.abs(calo_detdeta_mc1.GetBinContent(i) - calo_detdeta_mc2.GetBinContent(i)))
    calo_detdeta_dev12.SetBinError(i, np.sqrt(calo_detdeta_mc1.GetBinError(i)**2 + calo_detdeta_mc2.GetBinError(i)**2))
calo_detdeta_dev12.Divide(calo_detdeta_mc1)

hcal_detdeta_dev12 = TH1F(hcal_detdeta_mc1.Clone("hcal_detdeta_dev_percent"))
hcal_detdeta_dev12.SetDirectory(0)
hcal_detdeta_dev12.SetXTitle("#eta")
hcal_detdeta_dev12.SetYTitle("dE_{T}/d#eta dev. [GeV]")
for i in range(1, hcal_detdeta_dev12.GetNbinsX() + 1):
    hcal_detdeta_dev12.SetBinContent(i, np.abs(hcal_detdeta_mc1.GetBinContent(i) - hcal_detdeta_mc2.GetBinContent(i)))
    hcal_detdeta_dev12.SetBinError(i, np.sqrt(hcal_detdeta_mc1.GetBinError(i)**2 + hcal_detdeta_mc2.GetBinError(i)**2))
hcal_detdeta_dev12.Divide(hcal_detdeta_mc1)

canvas = TCanvas("canvas","",600,500)
leg = ROOT.TLegend(.45,.55,.8,.89)
leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
leg.AddEntry("","0-5% cent.","")
canvas.SetLeftMargin(0.18)
leg.SetBorderSize(0)
leg.AddEntry(emcal_detdeta_dev12,"EMCal","lep")
leg.AddEntry(hcal_detdeta_dev12,"HCal","lep")
leg.AddEntry(calo_detdeta_dev12,"Full Calo","lep")
leg.AddEntry(ihcal_detdeta_dev12,"IHCal","lep")
leg.AddEntry(ohcal_detdeta_dev12,"OHCal","lep")
emcal_detdeta_dev12.SetStats(0)
emcal_detdeta_dev12.SetLineColor(colors[0])
emcal_detdeta_dev12.SetMarkerColor(colors[0])
ihcal_detdeta_dev12.SetStats(0)
ihcal_detdeta_dev12.SetLineColor(colors[1])
ihcal_detdeta_dev12.SetMarkerColor(colors[1])
ohcal_detdeta_dev12.SetStats(0)
ohcal_detdeta_dev12.SetLineColor(colors[3])
ohcal_detdeta_dev12.SetMarkerColor(colors[3])
calo_detdeta_dev12.SetStats(0)
calo_detdeta_dev12.SetLineColor(colors[4])
calo_detdeta_dev12.SetMarkerColor(colors[4])
hcal_detdeta_dev12.SetStats(0)
hcal_detdeta_dev12.SetLineColor(colors[5])
hcal_detdeta_dev12.SetMarkerColor(colors[5])
emcal_detdeta_dev12.GetYaxis().SetRangeUser(0,0.1)
emcal_detdeta_dev12.SetXTitle("#eta")
emcal_detdeta_dev12.SetYTitle("MC Rap. Dep. dE_{T}/d#eta Uncertainty")
emcal_detdeta_dev12.GetYaxis().SetTitleOffset(1.8)
emcal_detdeta_dev12.Draw()
ihcal_detdeta_dev12.Draw('same')
ohcal_detdeta_dev12.Draw('same')
calo_detdeta_dev12.Draw('same')
hcal_detdeta_dev12.Draw('same')
leg.Draw()
canvas.SaveAs('syst_plots/MC_rap_dep_syst.png')

cents = ['0-5','5-10','10-20','20-30','30-40','40-50','50-60']
tag = 'MC_rap_dep'
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

    emcal_detdeta_mc = TH1F(h_eT_data_emcal.Clone("emcal_detdeta_mc"))
    emcal_detdeta_mc.Divide(h_emcal_correction)
    ihcal_detdeta_mc = TH1F(h_eT_data_ihcal.Clone("ihcal_detdeta_mc"))
    ihcal_detdeta_mc.Divide(h_ihcal_correction)
    ohcal_detdeta_mc = TH1F(h_eT_data_ohcal.Clone("ohcal_detdeta_mc"))
    ohcal_detdeta_mc.Divide(h_ohcal_correction)
    calo_detdeta_mc = TH1F(h_eT_data_calo.Clone("calo_detdeta_mc"))
    calo_detdeta_mc.Divide(h_calo_correction)
    hcal_detdeta_mc = TH1F(h_eT_data_hcal.Clone("hcal_detdeta_mc"))
    hcal_detdeta_mc.Divide(h_hcal_correction) 

    emcal_detdeta_dev = emcal_detdeta_mc.Clone("emcal_detdeta_dev")
    for i in range(1, emcal_detdeta_mc.GetNbinsX() + 1):
        emcal_detdeta_dev.SetBinContent(i, emcal_detdeta_dev12.GetBinContent(i)*emcal_detdeta_mc.GetBinContent(i))
        emcal_detdeta_dev.SetBinError(i, 0)
    
    ihcal_detdeta_dev = ihcal_detdeta_mc.Clone("ihcal_detdeta_dev")
    for i in range(1, ihcal_detdeta_mc.GetNbinsX() + 1):
        ihcal_detdeta_dev.SetBinContent(i, ihcal_detdeta_dev12.GetBinContent(i)*ihcal_detdeta_mc.GetBinContent(i))
        ihcal_detdeta_dev.SetBinError(i, 0)
    
    ohcal_detdeta_dev = ohcal_detdeta_mc.Clone("ohcal_detdeta_dev")
    for i in range(1, ohcal_detdeta_mc.GetNbinsX() + 1):
        ohcal_detdeta_dev.SetBinContent(i, ohcal_detdeta_dev12.GetBinContent(i)*ohcal_detdeta_mc.GetBinContent(i))
        ohcal_detdeta_dev.SetBinError(i, 0)
    
    calo_detdeta_dev = calo_detdeta_mc.Clone("calo_detdeta_dev")
    for i in range(1, calo_detdeta_mc.GetNbinsX() + 1):
        calo_detdeta_dev.SetBinContent(i, calo_detdeta_dev12.GetBinContent(i)*calo_detdeta_mc.GetBinContent(i))
        calo_detdeta_dev.SetBinError(i, 0)
        
    hcal_detdeta_dev = hcal_detdeta_mc.Clone("hcal_detdeta_dev")
    for i in range(1, hcal_detdeta_mc.GetNbinsX() + 1):
        hcal_detdeta_dev.SetBinContent(i, hcal_detdeta_dev12.GetBinContent(i)*hcal_detdeta_mc.GetBinContent(i))
        hcal_detdeta_dev.SetBinError(i, 0)
    
    emcal_detdeta_dev12.Write()
    ihcal_detdeta_dev12.Write()
    ohcal_detdeta_dev12.Write()
    calo_detdeta_dev12.Write()
    hcal_detdeta_dev12.Write()
    emcal_detdeta_dev.Write()
    ihcal_detdeta_dev.Write()
    ohcal_detdeta_dev.Write()
    calo_detdeta_dev.Write()
    hcal_detdeta_dev.Write()
    file.Close()