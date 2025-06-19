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
emcal_rho_corr = [[ 0.9197, 0.9421, 0.9521, 0.9396, 0.9582, 0.9157, 0.9758], [0.9787, 0.9833, 0.9812, 0.9790, 0.9799, 0.8332, 0.9852], [0.9901, 0.9940, 0.9961, 0.9963, 0.9970, 0.8674, 0.9605]]
hcal_rho_corr = [[0.8406, 0.8408, 0.8338, 0.8662, 0.8576, 0.5406, 0.9458], [0.9619, 0.9713, 0.9671, 0.9647, 0.9576, 0.5061, 0.9578], [0.9932, 0.9921, 0.9892, 0.9900, 0.9849, 0.4699, 0.9153]]
calo_rho_corr = [[0.8484, 0.8723, 0.8948, 0.8924, 0.9146, 0.7243, 0.9624], [0.9742, 0.9812, 0.9807, 0.9773, 0.9764, 0.6318, 0.9759], [0.9896, 0.9924, 0.9938, 0.9943, 0.9938, 0.6403, 0.9483]]
for c, cent in enumerate(cents):
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

	emcal_ratio_dev = []
	ihcal_ratio_dev = []
	ohcal_ratio_dev = []
	calo_ratio_dev = []
	hcal_ratio_dev = []

	for i in range(len(tag)-1):
		filename = '../fixed_build/dETdeta_variation_'+tag[i]+'_'+cent+'.root'
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
		emcal_dev[i].Rebin(4)
		emcal_dev[i].Scale(1.0/4)
		ihcal_dev[i].Rebin(4)
		ihcal_dev[i].Scale(1.0/4)
		ohcal_dev[i].Rebin(4)
		ohcal_dev[i].Scale(1.0/4)
		calo_dev[i].Rebin(4)
		calo_dev[i].Scale(1.0/4)
		hcal_dev[i].Rebin(4)
		hcal_dev[i].Scale(1.0/4)
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

	emcal_detdeta.Rebin(4)
	emcal_detdeta.Scale(1.0/4)
	ihcal_detdeta.Rebin(4)
	ihcal_detdeta.Scale(1.0/4)
	ohcal_detdeta.Rebin(4)
	ohcal_detdeta.Scale(1.0/4)
	calo_detdeta.Rebin(4)
	calo_detdeta.Scale(1.0/4)
	hcal_detdeta.Rebin(4)
	hcal_detdeta.Scale(1.0/4)

	for i in range(len(emcal_dev)):
		print(f"variation {i}",end=': ')
		for j in range(1, emcal_dev[i].GetNbinsX() + 1):
			emcal_dev[i].SetBinError(j,0)
			ihcal_dev[i].SetBinError(j,0)
			ohcal_dev[i].SetBinError(j,0)
			calo_dev[i].SetBinError(j,0)
			hcal_dev[i].SetBinError(j,0)
			print(emcal_dev[i].GetBinContent(j),end=',')
		print()

	emcal_detdeta_ratio = TH1F(emcal_detdeta.Clone("emcal_detdeta_ratio"))
	ihcal_detdeta_ratio = TH1F(ihcal_detdeta.Clone("ihcal_detdeta_ratio"))
	ohcal_detdeta_ratio = TH1F(ohcal_detdeta.Clone("ohcal_detdeta_ratio"))
	calo_detdeta_ratio = TH1F(calo_detdeta.Clone("calo_detdeta_ratio"))
	hcal_detdeta_ratio = TH1F(hcal_detdeta.Clone("hcal_detdeta_ratio"))

	print("nominal ratio:",end='')
	for i in range(1, emcal_detdeta_ratio.GetNbinsX() + 1):
		emcal_detdeta_ratio.SetBinContent(i, emcal_detdeta.GetBinContent(i)/emcal_detdeta.GetBinContent(emcal_detdeta.GetNbinsX() + 1 - i))
		ihcal_detdeta_ratio.SetBinContent(i, ihcal_detdeta.GetBinContent(i)/ihcal_detdeta.GetBinContent(ihcal_detdeta.GetNbinsX() + 1 - i))
		ohcal_detdeta_ratio.SetBinContent(i, ohcal_detdeta.GetBinContent(i)/ohcal_detdeta.GetBinContent(ohcal_detdeta.GetNbinsX() + 1 - i))
		calo_detdeta_ratio.SetBinContent(i, calo_detdeta.GetBinContent(i)/calo_detdeta.GetBinContent(calo_detdeta.GetNbinsX() + 1 - i))
		hcal_detdeta_ratio.SetBinContent(i, hcal_detdeta.GetBinContent(i)/hcal_detdeta.GetBinContent(hcal_detdeta.GetNbinsX() + 1 - i))
		print(emcal_detdeta_ratio.GetBinContent(i),end=', ')
	print()

	for i in range(len(emcal_dev)):
		emcal_ratio_dev.append(emcal_dev[i].Clone(f"emcal_ratio_{i}"))
		ihcal_ratio_dev.append(ihcal_dev[i].Clone(f"ihcal_ratio_{i}"))
		ohcal_ratio_dev.append(ohcal_dev[i].Clone(f"ohcal_ratio_{i}"))
		calo_ratio_dev.append(calo_dev[i].Clone(f"calo_ratio_{i}"))
		hcal_ratio_dev.append(hcal_dev[i].Clone(f"hcal_ratio_{i}"))

	for i in range(len(emcal_dev)):
		print(f"variation ratio {i}",end=': ')
		for j in range(1, emcal_dev[i].GetNbinsX() + 1):
			if j <= 3: ie = j - 1
			else: ie = 6 - j
			print(c, ie)
			error = np.sqrt((emcal_dev[i].GetBinContent(j)/emcal_detdeta.GetBinContent(j))**2 + (emcal_dev[i].GetBinContent(emcal_dev[i].GetNbinsX() + 1 - j)/emcal_detdeta.GetBinContent(emcal_detdeta.GetNbinsX() + 1 - j))**2 - 2*(emcal_rho_corr[ie][c]*emcal_dev[i].GetBinContent(j)*emcal_dev[i].GetBinContent(emcal_dev[i].GetNbinsX() + 1 - j))/(emcal_detdeta.GetBinContent(j)*emcal_detdeta.GetBinContent(emcal_detdeta.GetNbinsX() + 1 - j)))
			emcal_ratio_dev[i].SetBinContent(j, error*emcal_detdeta_ratio.GetBinContent(j))
			error = np.sqrt((ihcal_dev[i].GetBinContent(j)/ihcal_detdeta.GetBinContent(j))**2 + (ihcal_dev[i].GetBinContent(ihcal_dev[i].GetNbinsX() + 1 - j)/ihcal_detdeta.GetBinContent(ihcal_detdeta.GetNbinsX() + 1 - j))**2 - 2*(hcal_rho_corr[ie][c]*ihcal_dev[i].GetBinContent(j)*ihcal_dev[i].GetBinContent(ihcal_dev[i].GetNbinsX() + 1 - j))/(ihcal_detdeta.GetBinContent(j)*ihcal_detdeta.GetBinContent(ihcal_detdeta.GetNbinsX() + 1 - j)))
			ihcal_ratio_dev[i].SetBinContent(j, error*ihcal_detdeta_ratio.GetBinContent(j))
			error = np.sqrt((ohcal_dev[i].GetBinContent(j)/ohcal_detdeta.GetBinContent(j))**2 + (ohcal_dev[i].GetBinContent(ohcal_dev[i].GetNbinsX() + 1 - j)/ohcal_detdeta.GetBinContent(ohcal_detdeta.GetNbinsX() + 1 - j))**2 - 2*(hcal_rho_corr[ie][c]*ohcal_dev[i].GetBinContent(j)*ohcal_dev[i].GetBinContent(ohcal_dev[i].GetNbinsX() + 1 - j))/(ohcal_detdeta.GetBinContent(j)*ohcal_detdeta.GetBinContent(ohcal_detdeta.GetNbinsX() + 1 - j)))
			ohcal_ratio_dev[i].SetBinContent(j, error*ohcal_detdeta_ratio.GetBinContent(j))
			error = np.sqrt((calo_dev[i].GetBinContent(j)/calo_detdeta.GetBinContent(j))**2 + (calo_dev[i].GetBinContent(calo_dev[i].GetNbinsX() + 1 - j)/calo_detdeta.GetBinContent(calo_detdeta.GetNbinsX() + 1 - j))**2 - 2*(calo_rho_corr[ie][c]*calo_dev[i].GetBinContent(j)*calo_dev[i].GetBinContent(calo_dev[i].GetNbinsX() + 1 - j))/(calo_detdeta.GetBinContent(j)*calo_detdeta.GetBinContent(calo_detdeta.GetNbinsX() + 1 - j)))
			calo_ratio_dev[i].SetBinContent(j, error*calo_detdeta_ratio.GetBinContent(j))
			error = np.sqrt((hcal_dev[i].GetBinContent(j)/hcal_detdeta.GetBinContent(j))**2 + (hcal_dev[i].GetBinContent(hcal_dev[i].GetNbinsX() + 1 - j)/hcal_detdeta.GetBinContent(hcal_detdeta.GetNbinsX() + 1 - j))**2 - 2*(hcal_rho_corr[ie][c]*hcal_dev[i].GetBinContent(j)*hcal_dev[i].GetBinContent(hcal_dev[i].GetNbinsX() + 1 - j))/(hcal_detdeta.GetBinContent(j)*hcal_detdeta.GetBinContent(hcal_detdeta.GetNbinsX() + 1 - j)))
			hcal_ratio_dev[i].SetBinContent(j, error*hcal_detdeta_ratio.GetBinContent(j))
			print(hcal_ratio_dev[i].GetBinContent(j),end=', ')
		print()

	emcal_total_dev = TH1F(emcal_ratio_dev[0].Clone("emcal_total_dev"))
	ihcal_total_dev = TH1F(ihcal_ratio_dev[0].Clone("ihcal_total_dev"))
	ohcal_total_dev = TH1F(ohcal_ratio_dev[0].Clone("ohcal_total_dev"))
	calo_total_dev = TH1F(calo_ratio_dev[0].Clone("calo_total_dev"))
	hcal_total_dev = TH1F(hcal_ratio_dev[0].Clone("hcal_total_dev"))

	emcal_total = np.zeros(6)
	ihcal_total = np.zeros(6)
	ohcal_total = np.zeros(6)
	calo_total = np.zeros(6)
	hcal_total = np.zeros(6)

	for i in range(len(emcal_ratio_dev)):
		for j in range(1, emcal_ratio_dev[i].GetNbinsX() + 1):
			emcal_total[j-1] += emcal_ratio_dev[i].GetBinContent(j)**2
	for i in range(len(ihcal_ratio_dev)):
		for j in range(1, ihcal_ratio_dev[i].GetNbinsX() + 1):
			ihcal_total[j-1] += ihcal_ratio_dev[i].GetBinContent(j)**2
	for i in range(len(ohcal_ratio_dev)):
		for j in range(1, ohcal_ratio_dev[i].GetNbinsX() + 1):
			ohcal_total[j-1] += ohcal_ratio_dev[i].GetBinContent(j)**2
	for i in range(len(calo_ratio_dev)):
		for j in range(1, calo_ratio_dev[i].GetNbinsX() + 1):
			calo_total[j-1] += calo_ratio_dev[i].GetBinContent(j)**2
	for i in range(len(hcal_ratio_dev)):
		for j in range(1, hcal_ratio_dev[i].GetNbinsX() + 1):
			hcal_total[j-1] += hcal_ratio_dev[i].GetBinContent(j)**2
			
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

	emcal_ratio_dev.append(emcal_total_dev)
	ihcal_ratio_dev.append(ihcal_total_dev)
	ohcal_ratio_dev.append(ohcal_total_dev)
	calo_ratio_dev.append(calo_total_dev)
	hcal_ratio_dev.append(hcal_total_dev)

	for i in range(1, emcal_detdeta_ratio.GetNbinsX() + 1):
		emcal_detdeta_ratio.SetBinError(i, emcal_ratio_dev[-1].GetBinContent(i))
		ihcal_detdeta_ratio.SetBinError(i, ihcal_ratio_dev[-1].GetBinContent(i))
		ohcal_detdeta_ratio.SetBinError(i, ohcal_ratio_dev[-1].GetBinContent(i))
		calo_detdeta_ratio.SetBinError(i, calo_ratio_dev[-1].GetBinContent(i))
		hcal_detdeta_ratio.SetBinError(i, hcal_ratio_dev[-1].GetBinContent(i))

	emcal_canvas = TCanvas("emcal_canvas","",500,600)
	leg = ROOT.TLegend(.45,.8,.8,.89)
	leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
	leg.AddEntry("",f"{cent} % cent.","")
	emcal_canvas.SetLeftMargin(0.15)
	leg.SetBorderSize(0)
	emcal_detdeta_ratio.SetLineColor(1)
	emcal_detdeta_ratio.SetMarkerColor(1)
	emcal_detdeta_ratio.SetYTitle("EMCal dE_{T}/d#eta(+#eta)/dE_{T}/d#eta(-#eta)")
	emcal_detdeta_ratio.SetXTitle("#eta")
	emcal_detdeta_ratio.GetXaxis().SetRangeUser(0,1)
	emcal_detdeta_ratio.GetYaxis().SetRangeUser(0.9,1.1)
	line = ROOT.TLine(0,1,1.1,1)
	line.SetLineColor(2)
	line.SetLineWidth(2)
	emcal_detdeta_ratio.Draw()
	line.Draw('same')
	emcal_detdeta_ratio.Draw('same')
	leg.Draw()
	emcal_canvas.SaveAs('emcal_eta_ratio_syst_'+cent+'.png')

	hcal_canvas = TCanvas("hcal_canvas","",500,600)
	leg = ROOT.TLegend(.45,.8,.8,.89)
	leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
	leg.AddEntry("",f"{cent} % cent.","")
	hcal_canvas.SetLeftMargin(0.15)
	leg.SetBorderSize(0)
	hcal_detdeta_ratio.SetLineColor(1)
	hcal_detdeta_ratio.SetMarkerColor(1)
	hcal_detdeta_ratio.SetYTitle("HCal dE_{T}/d#eta(+#eta)/dE_{T}/d#eta(-#eta)")
	hcal_detdeta_ratio.SetXTitle("#eta")
	hcal_detdeta_ratio.GetXaxis().SetRangeUser(0,1)
	hcal_detdeta_ratio.GetYaxis().SetRangeUser(0.9,1.1)
	line = ROOT.TLine(0,1,1.1,1)
	line.SetLineColor(2)
	line.SetLineWidth(2)
	hcal_detdeta_ratio.Draw()
	line.Draw('same')
	hcal_detdeta_ratio.Draw('same')
	leg.Draw()
	hcal_canvas.SaveAs('hcal_eta_ratio_syst_'+cent+'.png')

	calo_canvas = TCanvas("calo_canvas","",500,600)
	leg = ROOT.TLegend(.45,.8,.8,.89)
	leg.AddEntry("","#bf{#it{sPHENIX}} Internal","")
	leg.AddEntry("",f"{cent} % cent.","")
	calo_canvas.SetLeftMargin(0.15)
	leg.SetBorderSize(0)
	calo_detdeta_ratio.SetLineColor(1)
	calo_detdeta_ratio.SetMarkerColor(1)
	calo_detdeta_ratio.SetYTitle("Full Calo dE_{T}/d#eta(+#eta)/dE_{T}/d#eta(-#eta)")
	calo_detdeta_ratio.SetXTitle("#eta")
	calo_detdeta_ratio.GetXaxis().SetRangeUser(0,1)
	calo_detdeta_ratio.GetYaxis().SetRangeUser(0.9,1.1)
	line = ROOT.TLine(0,1,1.1,1)
	line.SetLineColor(2)
	line.SetLineWidth(2)
	calo_detdeta_ratio.Draw()
	line.Draw('same')
	calo_detdeta_ratio.Draw('same')
	leg.Draw()
	calo_canvas.SaveAs('calo_eta_ratio_syst_'+cent+'.png')
