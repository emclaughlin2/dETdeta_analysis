import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import numpy as np
from array import array
import pdb

tag = 'p011'
cents = ['0-5','5-15','15-25']
leg_tags = ['emcal hijing','ihcal hijing','ohcal hijing','emcal epos','ihcal epos','ohcal epos','emcal ampt','ihcal ampt','ohcal ampt']
emcal_detdeta_hijing = []
emcal_detdeta_epos = []
emcal_detdeta_ampt = []
ihcal_detdeta_hijing = []
ihcal_detdeta_epos = []
ihcal_detdeta_ampt = []
ohcal_detdeta_hijing = []
ohcal_detdeta_epos = []
ohcal_detdeta_ampt = []
for i, cent in enumerate(cents):
	infile = 'dETdeta_plots_w_hijing_reweight_p011_new_'+cent+'.root'
	f1 = ROOT.TFile.Open(infile)
	emcal_detdeta_hijing.append(TH1F(f1.Get("emcal_detdeta_hijing")))
	emcal_detdeta_epos.append(TH1F(f1.Get("emcal_detdeta_epos")))
	emcal_detdeta_ampt.append(TH1F(f1.Get("emcal_detdeta_ampt")))
	ihcal_detdeta_hijing.append(TH1F(f1.Get("ihcal_detdeta_hijing")))
	ihcal_detdeta_epos.append(TH1F(f1.Get("ihcal_detdeta_epos")))
	ihcal_detdeta_ampt.append(TH1F(f1.Get("ihcal_detdeta_ampt")))
	ohcal_detdeta_hijing.append(TH1F(f1.Get("ohcal_detdeta_hijing")))
	ohcal_detdeta_epos.append(TH1F(f1.Get("ohcal_detdeta_epos")))
	ohcal_detdeta_ampt.append(TH1F(f1.Get("ohcal_detdeta_ampt")))
	emcal_detdeta_hijing[i].SetDirectory(0)
	emcal_detdeta_epos[i].SetDirectory(0)
	emcal_detdeta_ampt[i].SetDirectory(0)
	ihcal_detdeta_hijing[i].SetDirectory(0)
	ihcal_detdeta_epos[i].SetDirectory(0)
	ihcal_detdeta_ampt[i].SetDirectory(0)
	ohcal_detdeta_hijing[i].SetDirectory(0)
	ohcal_detdeta_epos[i].SetDirectory(0)
	ohcal_detdeta_ampt[i].SetDirectory(0)
	f1.Close()

outfile = ROOT.TFile('dETdeta_cent_plots_w_hijing_reweight_p011_new.root','RECREATE')

emcal_detdeta_hijing_mean = np.zeros(3)
emcal_detdeta_epos_mean = np.zeros(3)
emcal_detdeta_ampt_mean = np.zeros(3)
ihcal_detdeta_hijing_mean = np.zeros(3)
ihcal_detdeta_epos_mean = np.zeros(3)
ihcal_detdeta_ampt_mean = np.zeros(3)
ohcal_detdeta_hijing_mean = np.zeros(3)
ohcal_detdeta_epos_mean = np.zeros(3)
ohcal_detdeta_ampt_mean = np.zeros(3)
    
for i in range(len(cents)):
    for j in range(6, 20):
        emcal_detdeta_hijing_mean[i] += emcal_detdeta_hijing[i].GetBinContent(j)/14.0
        emcal_detdeta_epos_mean[i] += emcal_detdeta_epos[i].GetBinContent(j)/14.0
        emcal_detdeta_ampt_mean[i] += emcal_detdeta_ampt[i].GetBinContent(j)/14.0
        ihcal_detdeta_hijing_mean[i] += ihcal_detdeta_hijing[i].GetBinContent(j)/14.0
        ihcal_detdeta_epos_mean[i] += ihcal_detdeta_epos[i].GetBinContent(j)/14.0
        ihcal_detdeta_ampt_mean[i] += ihcal_detdeta_ampt[i].GetBinContent(j)/14.0
        ohcal_detdeta_hijing_mean[i] += ohcal_detdeta_hijing[i].GetBinContent(j)/14.0
        ohcal_detdeta_epos_mean[i] += ohcal_detdeta_epos[i].GetBinContent(j)/14.0
        ohcal_detdeta_ampt_mean[i] += ohcal_detdeta_ampt[i].GetBinContent(j)/14.0

#for i in range(len(cents)):
#    print(emcal_detdeta_hijing_mean[i],emcal_detdeta_epos_mean[i],emcal_detdeta_ampt_mean[i],end=' ')
#    print(ihcal_detdeta_hijing_mean[i],ihcal_detdeta_epos_mean[i],ihcal_detdeta_ampt_mean[i],end=' ')
#    print(ohcal_detdeta_hijing_mean[i],ohcal_detdeta_epos_mean[i],ohcal_detdeta_ampt_mean[i])

emcal_detdeta_hijing_mean = emcal_detdeta_hijing_mean.tolist()
emcal_detdeta_epos_mean = emcal_detdeta_epos_mean.tolist()
emcal_detdeta_ampt_mean = emcal_detdeta_ampt_mean.tolist()
ihcal_detdeta_hijing_mean = ihcal_detdeta_hijing_mean.tolist()
ihcal_detdeta_epos_mean = ihcal_detdeta_epos_mean.tolist()
ihcal_detdeta_ampt_mean = ihcal_detdeta_ampt_mean.tolist()
ohcal_detdeta_hijing_mean = ohcal_detdeta_hijing_mean.tolist()
ohcal_detdeta_epos_mean = ohcal_detdeta_epos_mean.tolist()
ohcal_detdeta_ampt_mean = ohcal_detdeta_ampt_mean.tolist()

eh = array('f', emcal_detdeta_hijing_mean)
ih = array('f', ihcal_detdeta_hijing_mean)
oh = array('f', ohcal_detdeta_hijing_mean)
ee = array('f', emcal_detdeta_epos_mean)
ie = array('f', ihcal_detdeta_epos_mean)
oe = array('f', ohcal_detdeta_epos_mean)
ea = array('f', emcal_detdeta_ampt_mean)
ia = array('f', ihcal_detdeta_ampt_mean)
oa = array('f', ohcal_detdeta_ampt_mean)
x = [0,1,2]
xarray = array('f', x)

arrays = [eh, ih, oh, ee, ie, oe, ea, ia, oa]

# Create a canvas
canvas = ROOT.TCanvas("dETdeta_vs_cent", "", 800, 600)

# Create and plot TGraphs for each array
graphs = []
colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kOrange, ROOT.kCyan, ROOT.kMagenta, ROOT.kYellow, ROOT.kBlack, ROOT.kGray]
for i, arr in enumerate(arrays):
    graphs.append(ROOT.TGraph(3, xarray, arr))
    graphs[i].SetMarkerStyle(20)
    graphs[i].SetMarkerSize(1)
    graphs[i].SetMarkerColor(colors[i])
    graphs[i].GetYaxis().SetRangeUser(300,800)
    if i == 0:
        graphs[i].Draw("AP")
    else:
        graphs[i].Draw("P,same")

# Set the legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
for i, graph in enumerate(graphs):
    legend.AddEntry(graph,leg_tags[i], "p")
legend.Draw()

# Draw the canvas
canvas.Draw()

canvas.Write()
outfile.Close()
