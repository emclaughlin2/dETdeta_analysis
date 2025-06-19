import ROOT
import numpy as np
from array import array
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import csv

gROOT.LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
gROOT.ProcessLine("SetsPhenixStyle()")

# systematics 
# MC 
# MC rap dep
# run by run 
# had resp
# emcal 1
# emcal 2
# emcal 3
# ihcal 1
# ihcal 2
# ihcal 3
# ohcal 1
# ohcal 2
# ohcal 3
# zs
# z vertex

cents = ['0-5','5-10','10-20','20-30','30-40','40-50','50-60']
rho_corr = {'0-5':0.7755, '5-10':0.7658, '10-20':0.7518, '20-30':0.7306, '30-40':0.7001, '40-50':0.6443, '50-60':0.5859}
for cent in cents:
	y_max = {'0-5': 170, '5-10': 140, '10-20': 110, '20-30': 80, '30-40': 50, '40-50': 40, '50-60': 30}
	tag = ['MC','MC_rap_dep','run_by_run','had_resp','1_emsyst1','2_emsyst2','3_emsyst3','1_ihsyst1',
				 '2_ihsyst2','3_ihsyst3','1_ohsyst1','2_ohsyst2','3_ohsyst3','zs_60_30_30ADC','vz_-3cm','total']
	taglabels = ['MC','MC Rapidity Dep.','Acceptance','Had. Resp.','EMsyst1','EMsyst2','EMsyst3',
				 'IHsyst1','IHsyst2','IHsyst3','OHsyst1','OHsyst2','OHsyst3','ZS','Vz Res.','Total']
	rgb = [[230, 25, 75], [60, 180, 75], [255, 225, 25], [0, 130, 200], [245, 130, 48], [145, 30, 180], [70, 240, 240], [240, 50, 230], [210, 245, 60], [250, 190, 212], [0, 128, 128], [220, 190, 255], [170, 110, 40], [128, 128, 128], [128, 0, 0], [0, 0, 0], [128, 128, 0], [255, 215, 180], [0, 0, 128], [34, 139, 34]]
	colors = [TColor.GetColor(rgb[i][0],rgb[i][1],rgb[i][2]) for i in range(len(rgb))]

	emcal_dev = []
	hcal_dev = []

	hcal_emcal_ratio_dev = []

	for i in range(len(tag)-1):
		filename = '../fixed_build/dETdeta_variation_'+tag[i]+'_'+cent+'.root'
		f = ROOT.TFile.Open(filename)
		print(i, filename)
		emcal_dev.append(TH1F(f.Get("emcal_detdeta_dev")))
		hcal_dev.append(TH1F(f.Get("hcal_detdeta_dev")))
		if i == 0: # for use in final nominal plots 
			emcal_detdeta_dev = TH1F(f.Get("emcal_detdeta_dev"))
			hcal_detdeta_dev = TH1F(f.Get("hcal_detdeta_dev"))
			emcal_detdeta_dev.SetDirectory(0)
			hcal_detdeta_dev.SetDirectory(0)
		emcal_dev[i].Rebin(4)
		emcal_dev[i].Scale(1.0/4)
		hcal_dev[i].Rebin(4)
		hcal_dev[i].Scale(1.0/4)
		emcal_dev[i].SetDirectory(0)
		hcal_dev[i].SetDirectory(0)
		f.Close()

	for i in range(len(emcal_dev)):
		for j in range(1, emcal_dev[i].GetNbinsX() + 1):
			emcal_dev[i].SetBinError(j,0)
			hcal_dev[i].SetBinError(j,0)

	mcfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_run14_with_centbin_nozs_mc_reweight_'+cent+'_reweight_epos_2024.root'
	datafile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_ana450_2024p009_100_50_50_ZS_hcal_scaled_emcal_calib_iter26_nozs_data_noweight_'+cent+'.root'
	f1 = ROOT.TFile.Open(mcfile)
	h_emcal_correction = TH1F(f1.Get("h_emcal_correction"))
	h_hcal_correction = TH1F(f1.Get("h_hcal_correction"))
	h_emcal_correction.SetDirectory(0)
	h_hcal_correction.SetDirectory(0)
	f1.Close()
	f2 = ROOT.TFile.Open(datafile)
	h_eT_data_emcal = TH1F(f2.Get("h_eT_eta_emcal_profile_hist"))
	h_eT_data_hcal = TH1F(f2.Get("h_eT_eta_hcal_profile_hist"))
	h_eT_data_emcal.SetDirectory(0)
	h_eT_data_hcal.SetDirectory(0)
	f2.Close()

	emcal_detdeta = TH1F(h_eT_data_emcal.Clone("emcal_detdeta"))
	emcal_detdeta.Divide(h_emcal_correction)
	hcal_detdeta = TH1F(h_eT_data_hcal.Clone("hcal_detdeta"))
	hcal_detdeta.Divide(h_hcal_correction) 

	emcal_detdeta.Rebin(4)
	emcal_detdeta.Scale(1.0/4)
	hcal_detdeta.Rebin(4)
	hcal_detdeta.Scale(1.0/4)

	for i in range(len(emcal_dev)):
		print(f"variation {i}",end=': ')
		for j in range(1, emcal_dev[i].GetNbinsX() + 1):
			emcal_dev[i].SetBinError(j,0)
			hcal_dev[i].SetBinError(j,0)
			print(emcal_dev[i].GetBinContent(j),end=',')
		print()

	hcal_emcal_detdeta_ratio = TH1F(emcal_detdeta.Clone("emcal_hcal_detdeta_ratio"))

	print("nominal ratio:",end='')
	for i in range(1, hcal_emcal_detdeta_ratio.GetNbinsX() + 1):
		hcal_emcal_detdeta_ratio.SetBinContent(i, hcal_detdeta.GetBinContent(i)/emcal_detdeta.GetBinContent(i))
		print(hcal_emcal_detdeta_ratio.GetBinContent(i),end=', ')
	print()

	for i in range(len(emcal_dev)):
		hcal_emcal_ratio_dev.append(emcal_dev[i].Clone(f"hcal_emcal_ratio_{i}"))
	
	for i in range(len(emcal_dev)):
		print(f"variation ratio {i}",end=': ')
		for j in range(1, emcal_dev[i].GetNbinsX() + 1):
			error = np.sqrt((hcal_dev[i].GetBinContent(j)/hcal_detdeta.GetBinContent(j))**2 + (emcal_dev[i].GetBinContent(j)/emcal_detdeta.GetBinContent(j))**2 - 2*(rho_corr[cent]*hcal_dev[i].GetBinContent(j)*emcal_dev[i].GetBinContent(j))/(hcal_detdeta.GetBinContent(j)*emcal_detdeta.GetBinContent(j)))
			hcal_emcal_ratio_dev[i].SetBinContent(j, error*hcal_emcal_detdeta_ratio.GetBinContent(j))
			print(hcal_emcal_ratio_dev[i].GetBinContent(j),end=', ')
		print()

	hcal_emcal_total_dev = TH1F(hcal_emcal_ratio_dev[0].Clone("hcal_emcal_total_dev"))
	hcal_emcal_total = np.zeros(6)

	for i in range(len(hcal_emcal_ratio_dev)):
		for j in range(1, hcal_emcal_ratio_dev[i].GetNbinsX() + 1):
			hcal_emcal_total[j-1] += hcal_emcal_ratio_dev[i].GetBinContent(j)**2
			
	for i in range(1, hcal_emcal_total_dev.GetNbinsX() + 1):
			hcal_emcal_total_dev.SetBinContent(i, np.sqrt(hcal_emcal_total[i-1]))

	hcal_emcal_ratio_dev.append(hcal_emcal_total_dev)

	for i in range(1, hcal_emcal_detdeta_ratio.GetNbinsX() + 1):
		hcal_emcal_detdeta_ratio.SetBinError(i, hcal_emcal_ratio_dev[-1].GetBinContent(i))

	emcal_canvas = TCanvas("emcal_canvas","",500,600)
	leg = ROOT.TLegend(.45,.8,.8,.89)
	leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
	leg.AddEntry("",f"{cent} % cent.","")
	emcal_canvas.SetLeftMargin(0.15)
	leg.SetBorderSize(0)
	hcal_emcal_detdeta_ratio.SetLineColor(1)
	hcal_emcal_detdeta_ratio.SetMarkerColor(1)
	hcal_emcal_detdeta_ratio.SetYTitle("dE_{T}/d#eta_{EMCal}/dE_{T}/d#eta_{HCal}")
	hcal_emcal_detdeta_ratio.SetXTitle("#eta")
	hcal_emcal_detdeta_ratio.Draw()
	leg.Draw()
	emcal_canvas.SaveAs('hcal_emcal_ratio_syst_'+cent+'.png')