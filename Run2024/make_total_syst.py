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
for cent in cents:
	y_max = {'0-5': 170, '5-10': 140, '10-20': 110, '20-30': 80, '30-40': 50, '40-50': 40, '50-60': 30}
	tag = ['MC','MC_rap_dep','run_by_run','had_resp','1_emsyst1','2_emsyst2','3_emsyst3','1_ihsyst1',
				 '2_ihsyst2','3_ihsyst3','1_ohsyst1','2_ohsyst2','3_ohsyst3','zs_60_30_30ADC','vz_-3cm','total']
	taglabels = ['MC','MC Rapidity Dep.','Acceptance','Had. Resp.','EMsyst1','EMsyst2','EMsyst3',
				 'IHsyst1','IHsyst2','IHsyst3','OHsyst1','OHsyst2','OHsyst3','ZS','Vz Res.','Total']
	rgb = [[230, 25, 75], [60, 180, 75], [255, 225, 25], [0, 130, 200], [245, 130, 48], [145, 30, 180], [70, 240, 240], [240, 50, 230], [210, 245, 60], [250, 190, 212], [0, 128, 128], [220, 190, 255], [170, 110, 40], [128, 128, 128], [128, 0, 0], [0, 0, 0], [128, 128, 0], [255, 215, 180], [0, 0, 128], [34, 139, 34]]
	colors = [TColor.GetColor(rgb[i][0],rgb[i][1],rgb[i][2]) for i in range(len(rgb))]

	emcal_dev = []
	ihcal_dev = []
	ohcal_dev = []
	calo_dev = []
	hcal_dev = []

	for i in range(len(tag)-1):
		filename = 'fixed_build/dETdeta_variation_'+tag[i]+'_'+cent+'.root'
		f = ROOT.TFile.Open(filename)
		print(i, filename)
		emcal_dev.append(TH1F(f.Get("emcal_detdeta_dev")))
		ihcal_dev.append(TH1F(f.Get("ihcal_detdeta_dev")))
		ohcal_dev.append(TH1F(f.Get("ohcal_detdeta_dev")))
		calo_dev.append(TH1F(f.Get("calo_detdeta_dev")))
		hcal_dev.append(TH1F(f.Get("hcal_detdeta_dev")))
		if i == 0: # for use in final nominal plots 
			emcal_detdeta_dev = TH1F(f.Get("emcal_detdeta_dev"))
			ihcal_detdeta_dev = TH1F(f.Get("ihcal_detdeta_dev"))
			ohcal_detdeta_dev = TH1F(f.Get("ohcal_detdeta_dev"))
			calo_detdeta_dev = TH1F(f.Get("calo_detdeta_dev"))
			hcal_detdeta_dev = TH1F(f.Get("hcal_detdeta_dev"))
			emcal_detdeta_dev.SetDirectory(0)
			ihcal_detdeta_dev.SetDirectory(0)
			ohcal_detdeta_dev.SetDirectory(0)
			calo_detdeta_dev.SetDirectory(0)
			hcal_detdeta_dev.SetDirectory(0)
		emcal_dev[i].Rebin(2)
		emcal_dev[i].Scale(1.0/2)
		ihcal_dev[i].Rebin(2)
		ihcal_dev[i].Scale(1.0/2)
		ohcal_dev[i].Rebin(2)
		ohcal_dev[i].Scale(1.0/2)
		calo_dev[i].Rebin(2)
		calo_dev[i].Scale(1.0/2)
		hcal_dev[i].Rebin(2)
		hcal_dev[i].Scale(1.0/2)
		emcal_dev[i].SetDirectory(0)
		ihcal_dev[i].SetDirectory(0)
		ohcal_dev[i].SetDirectory(0)
		calo_dev[i].SetDirectory(0)
		hcal_dev[i].SetDirectory(0)
		f.Close()

	for i in range(len(emcal_dev)):
		for j in range(1, emcal_dev[i].GetNbinsX() + 1):
			emcal_dev[i].SetBinError(j,0)
			ihcal_dev[i].SetBinError(j,0)
			ohcal_dev[i].SetBinError(j,0)
			calo_dev[i].SetBinError(j,0)
			hcal_dev[i].SetBinError(j,0)

	emcal_total_dev = TH1F(emcal_dev[0].Clone("emcal_total_dev"))
	ihcal_total_dev = TH1F(ihcal_dev[0].Clone("ihcal_total_dev"))
	ohcal_total_dev = TH1F(ohcal_dev[0].Clone("ohcal_total_dev"))
	calo_total_dev = TH1F(calo_dev[0].Clone("calo_total_dev"))
	hcal_total_dev = TH1F(hcal_dev[0].Clone("hcal_total_dev"))

	emcal_total = np.zeros(12)
	ihcal_total = np.zeros(12)
	ohcal_total = np.zeros(12)
	calo_total = np.zeros(12)
	hcal_total = np.zeros(12)

	for i in range(len(emcal_dev)):
		for j in range(1, emcal_dev[i].GetNbinsX() + 1):
			emcal_total[j-1] += emcal_dev[i].GetBinContent(j)**2
	for i in range(len(ihcal_dev)):
		for j in range(1, ihcal_dev[i].GetNbinsX() + 1):
			ihcal_total[j-1] += ihcal_dev[i].GetBinContent(j)**2
	for i in range(len(ohcal_dev)):
		for j in range(1, ohcal_dev[i].GetNbinsX() + 1):
			ohcal_total[j-1] += ohcal_dev[i].GetBinContent(j)**2
	for i in range(len(calo_dev)):
		for j in range(1, calo_dev[i].GetNbinsX() + 1):
			calo_total[j-1] += calo_dev[i].GetBinContent(j)**2
	for i in range(len(hcal_dev)):
		for j in range(1, hcal_dev[i].GetNbinsX() + 1):
			hcal_total[j-1] += hcal_dev[i].GetBinContent(j)**2
			
	for i in range(1, emcal_total_dev.GetNbinsX() + 1):
			emcal_total_dev.SetBinContent(i, np.sqrt(emcal_total[i-1]))
	for i in range(1, ihcal_total_dev.GetNbinsX() + 1):
			ihcal_total_dev.SetBinContent(i, np.sqrt(ihcal_total[i-1]))
	for i in range(1, ohcal_total_dev.GetNbinsX() + 1):
			ohcal_total_dev.SetBinContent(i, np.sqrt(ohcal_total[i-1]))
	for i in range(1, calo_total_dev.GetNbinsX() + 1):
			calo_total_dev.SetBinContent(i, np.sqrt(calo_total[i-1]))
	for i in range(1, hcal_total_dev.GetNbinsX() + 1):
			hcal_total_dev.SetBinContent(i, np.sqrt(hcal_total[i-1]))

	emcal_dev.append(emcal_total_dev)
	ihcal_dev.append(ihcal_total_dev)
	ohcal_dev.append(ohcal_total_dev)
	calo_dev.append(calo_total_dev)
	hcal_dev.append(hcal_total_dev)

	testfile = "fixed_build/dETdeta_analysis_allruns_ana450_2024p009_100_50_50_ZS_hcal_scaled_emcal_calib_iter26_nozs_data_noweight_"+cent+".root"
	f = ROOT.TFile.Open(testfile)
	emcal_detdeta_dev = TH1F(f.Get("h_eT_eta_emcal_profile_hist").Clone("emcal_detdeta_dev"))
	ihcal_detdeta_dev = TH1F(f.Get("h_eT_eta_ihcal_profile_hist").Clone("ihcal_detdeta_dev"))
	ohcal_detdeta_dev = TH1F(f.Get("h_eT_eta_ohcal_profile_hist").Clone("ohcal_detdeta_dev"))
	calo_detdeta_dev = TH1F(f.Get("h_eT_eta_calo_profile_hist").Clone("calo_detdeta_dev"))
	hcal_detdeta_dev = TH1F(f.Get("h_eT_eta_hcal_profile_hist").Clone("hcal_detdeta_dev"))
	emcal_detdeta_dev.SetDirectory(0)
	ihcal_detdeta_dev.SetDirectory(0)
	ohcal_detdeta_dev.SetDirectory(0)
	calo_detdeta_dev.SetDirectory(0)
	hcal_detdeta_dev.SetDirectory(0)
	f.Close()

	for i in range(1, emcal_detdeta_dev.GetNbinsX() + 1):
		emcal_detdeta_dev.SetBinContent(i, emcal_dev[-1].GetBinContent((i-1)//2+1))
		emcal_detdeta_dev.SetBinError(i, 0)
	for i in range(1, ihcal_detdeta_dev.GetNbinsX() + 1):
		ihcal_detdeta_dev.SetBinContent(i, ihcal_dev[-1].GetBinContent((i-1)//2+1))
		ihcal_detdeta_dev.SetBinError(i, 0)
	for i in range(1, ohcal_detdeta_dev.GetNbinsX() + 1):
		ohcal_detdeta_dev.SetBinContent(i, ohcal_dev[-1].GetBinContent((i-1)//2+1))
		ohcal_detdeta_dev.SetBinError(i, 0)
	for i in range(1, calo_detdeta_dev.GetNbinsX() + 1):
		calo_detdeta_dev.SetBinContent(i, calo_dev[-1].GetBinContent((i-1)//2+1))
		calo_detdeta_dev.SetBinError(i, 0)
	for i in range(1, hcal_detdeta_dev.GetNbinsX() + 1):
		hcal_detdeta_dev.SetBinContent(i, hcal_dev[-1].GetBinContent((i-1)//2+1))
		hcal_detdeta_dev.SetBinError(i, 0)
		
	outfile = ROOT.TFile.Open('fixed_build/dETdeta_total_variation_'+cent+'.root',"RECREATE")
	emcal_detdeta_dev.Write()
	ihcal_detdeta_dev.Write()
	ohcal_detdeta_dev.Write()
	calo_detdeta_dev.Write()
	hcal_detdeta_dev.Write()

# systematics 
# 0 MC 
# 1 MC rap dep
# 2 run by run 
# 3 had resp
# 4 emcal 1
# 5 emcal 2
# 6 emcal 3
# 7 ihcal 1
# 8 ihcal 2
# 9 ihcal 3
# 10 ohcal 1
# 11 ohcal 2
# 12 ohcal 3
# 13 zs
# 14 z vertex

	emcal_canvas = TCanvas("emcal_canvas","",500,600)
	leg = ROOT.TLegend(.45,.4,.8,.89)
	leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
	leg.AddEntry("",f"{cent} % cent.","")
	emcal_canvas.SetLeftMargin(0.15)
	leg.SetBorderSize(0)
	for i in range(len(emcal_dev)):
		if i >= 7 and i <= 12:
			continue
		leg.AddEntry(emcal_dev[i],taglabels[i],"lep")
		emcal_dev[i].SetStats(0)
		emcal_dev[i].SetMarkerStyle(20)
		emcal_dev[i].SetLineColor(colors[i])
		emcal_dev[i].SetMarkerColor(colors[i])
		if i == 0:
			emcal_dev[i].GetYaxis().SetRangeUser(0,y_max[cent])
			emcal_dev[i].SetXTitle("#eta")
			emcal_dev[i].SetYTitle("EMCal dE_{T}/d#eta Uncertainty  [GeV]")
			emcal_dev[i].Draw()
		else:
			emcal_dev[i].Draw('same')
	leg.Draw()
	emcal_canvas.SaveAs('syst_plots/emcal_total_syst_'+cent+'.png')

	ihcal_canvas = TCanvas("ihcal_canvas","",500,600)
	leg = ROOT.TLegend(.45,.43,.8,.89)
	leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
	leg.AddEntry("",f"{cent} % cent.","")
	ihcal_canvas.SetLeftMargin(0.15)
	leg.SetBorderSize(0)
	for i in range(len(ihcal_dev)):
		if i >= 4 and i <= 6: continue
		if i >= 10 and i <= 12: continue
		leg.AddEntry(ihcal_dev[i],taglabels[i],"lep")
		ihcal_dev[i].SetStats(0)
		ihcal_dev[i].SetMarkerStyle(20)
		ihcal_dev[i].SetLineColor(colors[i])
		ihcal_dev[i].SetMarkerColor(colors[i])
		if i == 0:
			ihcal_dev[i].GetYaxis().SetRangeUser(0,y_max[cent])
			ihcal_dev[i].SetXTitle("#eta")
			ihcal_dev[i].SetYTitle("IHCal dE_{T}/d#eta Uncertainty [GeV]")
			ihcal_dev[i].Draw()
		else:
			ihcal_dev[i].Draw('same')
	leg.Draw()
	ihcal_canvas.SaveAs('syst_plots/ihcal_total_syst_'+cent+'.png')

	ohcal_canvas = TCanvas("ohcal_canvas","",500,600)
	leg = ROOT.TLegend(.45,.43,.7,.89)
	leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
	leg.AddEntry("",f"{cent} % cent.","")
	ohcal_canvas.SetLeftMargin(0.15)
	leg.SetBorderSize(0)
	for i in range(len(ohcal_dev)):
		if i >= 4 and i <= 9: continue
		leg.AddEntry(ohcal_dev[i],taglabels[i],"lep")
		ohcal_dev[i].SetStats(0)
		ohcal_dev[i].SetMarkerStyle(20)
		ohcal_dev[i].SetLineColor(colors[i])
		ohcal_dev[i].SetMarkerColor(colors[i])
		if i == 0:
			ohcal_dev[i].GetYaxis().SetRangeUser(0,y_max[cent])
			ohcal_dev[i].SetXTitle("#eta")
			ohcal_dev[i].SetYTitle("OHCal dE_{T}/d#eta Uncertainty [GeV]")
			ohcal_dev[i].Draw()
		else:
			ohcal_dev[i].Draw('same')
	leg.Draw()
	ohcal_canvas.SaveAs('syst_plots/ohcal_total_syst_'+cent+'.png')

	calo_canvas = TCanvas("calo_canvas","",500,600)
	leg = ROOT.TLegend(.18,.45,.9,.89)
	leg.SetNColumns(2)
	leg.SetTextSize(0.045)
	leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
	leg.AddEntry("",f"{cent} % cent.","")
	calo_canvas.SetLeftMargin(0.15)
	leg.SetBorderSize(0)
	for i in range(len(calo_dev)):
		leg.AddEntry(calo_dev[i],taglabels[i],"lep")
		calo_dev[i].SetStats(0)
		calo_dev[i].SetMarkerStyle(20)
		calo_dev[i].SetLineColor(colors[i])
		calo_dev[i].SetMarkerColor(colors[i])
		if i == 0:
			calo_dev[i].GetYaxis().SetRangeUser(0,y_max[cent])
			calo_dev[i].SetXTitle("#eta")
			calo_dev[i].SetYTitle("dE_{T}/d#eta Uncertainty  [GeV]")
			calo_dev[i].Draw()
		else:
			calo_dev[i].Draw('same')
	leg.Draw()
	calo_canvas.SaveAs('syst_plots/calo_total_syst_'+cent+'.png')

	hcal_canvas = TCanvas("hcal_canvas","",500,600)
	leg = ROOT.TLegend(.18,.55,.9,.89)
	leg.SetNColumns(2)
	leg.SetTextSize(0.045)
	leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
	leg.AddEntry("",f"{cent} % cent.","")
	hcal_canvas.SetLeftMargin(0.15)
	leg.SetBorderSize(0)
	for i in range(len(hcal_dev)):
		if i >= 4 and i <= 6: continue
		leg.AddEntry(hcal_dev[i],taglabels[i],"lep")
		hcal_dev[i].SetStats(0)
		hcal_dev[i].SetMarkerStyle(20)
		hcal_dev[i].SetLineColor(colors[i])
		hcal_dev[i].SetMarkerColor(colors[i])
		if i == 0:
			hcal_dev[i].GetYaxis().SetRangeUser(0,y_max[cent])
			hcal_dev[i].SetXTitle("#eta")
			hcal_dev[i].SetYTitle("HCal dE_{T}/d#eta Uncertainty [GeV]")
			hcal_dev[i].Draw()
		else:
			hcal_dev[i].Draw('same')
	leg.Draw()
	hcal_canvas.SaveAs('syst_plots/hcal_total_syst_'+cent+'.png')

	emcal_canvas.Write()
	ihcal_canvas.Write()
	ohcal_canvas.Write()
	calo_canvas.Write()
	hcal_canvas.Write()
	outfile.Write()
	outfile.Close()

	mcfile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_run14_with_centbin_nozs_mc_reweight_'+cent+'_reweight_epos_2024.root'
	datafile = '/sphenix/user/egm2153/calib_study/detdeta/analysis/Run2024/fixed_build/dETdeta_analysis_allruns_ana450_2024p009_100_50_50_ZS_hcal_scaled_emcal_calib_iter26_nozs_data_noweight_'+cent+'.root'
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

	emcal_avg = emcal_detdeta_mean.GetBinContent(1)
	ihcal_avg = ihcal_detdeta_mean.GetBinContent(1)
	ohcal_avg = ohcal_detdeta_mean.GetBinContent(1)
	calo_avg = calo_detdeta_mean.GetBinContent(1)
	hcal_avg = hcal_detdeta_mean.GetBinContent(1)

	# systematics 
	# 0 MC 
	# 1 MC rap dep
	# 2 run by run 
	# 3 had resp
	# 4 emcal 1
	# 5 emcal 2
	# 6 emcal 3
	# 7 ihcal 1
	# 8 ihcal 2
	# 9 ihcal 3
	# 10 ohcal 1
	# 11 ohcal 2
	# 12 ohcal 3
	# 13 zs
	# 14 z vertex

	si = [0,0,1,2,3,3,3,3,3,3,3,3,3,4,5,6]
	emcal_dev_avgt = np.zeros(20)
	ihcal_dev_avgt = np.zeros(20)
	ohcal_dev_avgt = np.zeros(20)
	calo_dev_avgt = np.zeros(20)
	hcal_dev_avgt = np.zeros(20)

	for i in range(len(emcal_dev)):
		for j in range(1,13):
			emcal_dev_avgt[i] += emcal_dev[i].GetBinContent(j)/12.0
			ihcal_dev_avgt[i] += ihcal_dev[i].GetBinContent(j)/12.0
			ohcal_dev_avgt[i] += ohcal_dev[i].GetBinContent(j)/12.0
			calo_dev_avgt[i] += calo_dev[i].GetBinContent(j)/12.0
			hcal_dev_avgt[i] += hcal_dev[i].GetBinContent(j)/12.0
		emcal_dev_avgt[i] = emcal_dev_avgt[i]**2
		ihcal_dev_avgt[i] = ihcal_dev_avgt[i]**2
		ohcal_dev_avgt[i] = ohcal_dev_avgt[i]**2
		calo_dev_avgt[i] = calo_dev_avgt[i]**2
		hcal_dev_avgt[i] = hcal_dev_avgt[i]**2

	emcal_dev_avg = np.zeros(7)
	ihcal_dev_avg = np.zeros(7)
	ohcal_dev_avg = np.zeros(7)
	calo_dev_avg = np.zeros(7)
	hcal_dev_avg = np.zeros(7)

	for i in range(len(emcal_dev)):
		emcal_dev_avg[si[i]] += emcal_dev_avgt[i]
		ihcal_dev_avg[si[i]] += ihcal_dev_avgt[i]
		ohcal_dev_avg[si[i]] += ohcal_dev_avgt[i]
		calo_dev_avg[si[i]] += calo_dev_avgt[i]
		hcal_dev_avg[si[i]] += hcal_dev_avgt[i]

	for i in range(7):
		emcal_dev_avg[i] = np.sqrt(emcal_dev_avg[i])
		ihcal_dev_avg[i] = np.sqrt(ihcal_dev_avg[i])
		ohcal_dev_avg[i] = np.sqrt(ohcal_dev_avg[i])
		calo_dev_avg[i] = np.sqrt(calo_dev_avg[i])
		hcal_dev_avg[i] = np.sqrt(hcal_dev_avg[i])

	for i in range(7):
		emcal_dev_avg[i] /= emcal_avg
		ihcal_dev_avg[i] /= ihcal_avg
		ohcal_dev_avg[i] /= ohcal_avg
		calo_dev_avg[i] /= calo_avg
		hcal_dev_avg[i] /= hcal_avg

	emcalfilename = 'syst_plots/avg_syst_unc.csv'
	ihcalfilename = 'syst_plots/avg_syst_unc.csv'
	ohcalfilename = 'syst_plots/avg_syst_unc.csv'
	calofilename = 'syst_plots/avg_syst_unc.csv'
	hcalfilename = 'syst_plots/avg_syst_unc.csv'

	with open(emcalfilename, mode='a', newline='') as emcalfile:
		writer = csv.writer(emcalfile)
		writer.writerow(emcal_dev_avg)

	with open(hcalfilename, mode='a', newline='') as hcalfile:
		writer = csv.writer(hcalfile)
		writer.writerow(hcal_dev_avg)
		
	with open(calofilename, mode='a', newline='') as calofile:
		writer = csv.writer(calofile)
		writer.writerow(calo_dev_avg)

	with open(ihcalfilename, mode='a', newline='') as ihcalfile:
		writer = csv.writer(ihcalfile)
		writer.writerow(ihcal_dev_avg)
		
	with open(ohcalfilename, mode='a', newline='') as ohcalfile:
		writer = csv.writer(ohcalfile)
		writer.writerow(ohcal_dev_avg)
