import ROOT
import sys

cent = sys.argv[1]

input_file = ROOT.TFile("dETdeta_analysis_pedestal_subtraction_54256_nozs_no_cross_calib_default_bins.root", "READ")
h_eT_emcal = input_file.Get("h_eT_emcal")
h_eT_ihcal = input_file.Get("h_eT_ihcal")
h_eT_ohcal = input_file.Get("h_eT_ohcal")
h_eT_calo = input_file.Get("h_eT_calo")
h_eT_hcal = input_file.Get("h_eT_hcal")

emcal_detdeta_dev = h_eT_emcal.Clone("emcal_detdeta_dev")
ihcal_detdeta_dev = h_eT_ihcal.Clone("ihcal_detdeta_dev")
ohcal_detdeta_dev = h_eT_ohcal.Clone("ohcal_detdeta_dev")
calo_detdeta_dev = h_eT_calo.Clone("calo_detdeta_dev")
hcal_detdeta_dev = h_eT_hcal.Clone("hcal_detdeta_dev")

output_file = ROOT.TFile("dETdeta_variation_ZS_w_hcal_"+cent+".root", "RECREATE")
emcal_detdeta_dev.Write()
ihcal_detdeta_dev.Write()
ohcal_detdeta_dev.Write()
calo_detdeta_dev.Write()
hcal_detdeta_dev.Write()
output_file.Close()
input_file.Close()

