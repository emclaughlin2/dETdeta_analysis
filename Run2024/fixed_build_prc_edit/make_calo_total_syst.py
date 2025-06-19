import ROOT
import numpy as np

def compute_uncorrelated_uncertainty(emcal_hist, hcal_hist, rho):
    """Computes the uncorrelated uncertainties for EMCal and HCal separately."""
    uncorrelated_emcal = emcal_hist.Clone("uncorrelated_emcal_dev")
    uncorrelated_hcal = hcal_hist.Clone("uncorrelated_hcal_dev")

    uncorrelated_emcal.Reset()
    uncorrelated_hcal.Reset()

    for bin_idx in range(1, emcal_hist.GetNbinsX() + 1):
        sigma_emcal = emcal_hist.GetBinContent(bin_idx)
        sigma_hcal = hcal_hist.GetBinContent(bin_idx)

        # Compute uncorrelated uncertainties for EMCal and HCal separately
        sigma_uncorrelated_emcal = np.sqrt((1 - rho) * sigma_emcal**2)
        sigma_uncorrelated_hcal = np.sqrt((1 - rho) * sigma_hcal**2)

        uncorrelated_emcal.SetBinContent(bin_idx, sigma_uncorrelated_emcal)
        uncorrelated_hcal.SetBinContent(bin_idx, sigma_uncorrelated_hcal)

    return uncorrelated_emcal, uncorrelated_hcal

tags = [
    "1_emsyst1", "2_emsyst2", "3_emsyst3",
    "1_ihsyst1", "2_ihsyst2", "3_ihsyst3",
    "1_ohsyst1", "2_ohsyst2", "3_ohsyst3",
    "MC_rap_dep"
]

cents = ["0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60"]

hist_names = ["emcal_detdeta_dev", "hcal_detdeta_dev"]

# Correlation coefficient for MC_rap_dep only
# Assume a correlation coefficient (adjust as needed)
rho_correlated = [0.7755, 0.7658, 0.7518, 0.7306, 0.7001, 0.6443, 0.5859]

for c, cent in enumerate(cents):
    sum_uncorrelated_emcal = None
    sum_uncorrelated_hcal = None

    # Find a valid reference histogram for binning
    reference_hist = None
    for tag in tags:
        file_name = f"fixed_build/dETdeta_variation_{tag}_{cent}.root"
        f = ROOT.TFile.Open(file_name, "READ")
        if not f or f.IsZombie():
            continue

        for hist_name in hist_names:
            hist = f.Get(hist_name)
            if hist and hist.InheritsFrom("TH1"):
                reference_hist = hist.Clone()
                reference_hist.SetDirectory(0)
                break
        f.Close()
        if reference_hist:
            break

    if not reference_hist:
        print(f"Error: No valid reference histogram found for cent {cent}. Skipping.")
        continue

    # Initialize empty histograms for quadratic sums
    sum_uncorrelated_emcal = reference_hist.Clone("sum_uncorrelated_emcal")
    sum_uncorrelated_hcal = reference_hist.Clone("sum_uncorrelated_hcal")
    sum_uncorrelated_emcal.Reset()
    sum_uncorrelated_hcal.Reset()
    sum_uncorrelated_emcal.SetDirectory(0)
    sum_uncorrelated_hcal.SetDirectory(0)

    # Loop over tag files and accumulate quadratic sums
    for tag in tags:
        file_name = f"fixed_build/dETdeta_variation_{tag}_{cent}.root"
        f = ROOT.TFile.Open(file_name, "READ")

        if not f or f.IsZombie():
            print(f"Warning: Could not open {file_name}")
            continue

        emcal_hist = f.Get("emcal_detdeta_dev")
        hcal_hist = f.Get("hcal_detdeta_dev")

        if not emcal_hist or not hcal_hist or not emcal_hist.InheritsFrom("TH1") or not hcal_hist.InheritsFrom("TH1"):
            print(f"Warning: Missing histograms in {file_name}")
            f.Close()
            continue

        emcal_hist.SetDirectory(0)
        hcal_hist.SetDirectory(0)
        f.Close()

        # Determine rho value
        rho = rho_correlated[c] if tag == "MC_rap_dep" else 0  # Fully uncorrelated except for MC_rap_dep

        # Compute uncorrelated uncertainties
        uncorrelated_emcal, uncorrelated_hcal = compute_uncorrelated_uncertainty(emcal_hist, hcal_hist, rho)

        # Sum the uncorrelated uncertainties in quadrature
        for bin_idx in range(1, sum_uncorrelated_emcal.GetNbinsX() + 1):
            current_val_emcal = sum_uncorrelated_emcal.GetBinContent(bin_idx)
            new_val_emcal = uncorrelated_emcal.GetBinContent(bin_idx)
            sum_uncorrelated_emcal.SetBinContent(bin_idx, np.sqrt(current_val_emcal**2 + new_val_emcal**2))

            current_val_hcal = sum_uncorrelated_hcal.GetBinContent(bin_idx)
            new_val_hcal = uncorrelated_hcal.GetBinContent(bin_idx)
            sum_uncorrelated_hcal.SetBinContent(bin_idx, np.sqrt(current_val_hcal**2 + new_val_hcal**2))

    # Save the results
    output_file = ROOT.TFile.Open(f"dETdeta_calo_variation_{cent}.root", "RECREATE")
    sum_uncorrelated_emcal.Write("uncorrelated_total_emcal_dev")
    sum_uncorrelated_hcal.Write("uncorrelated_total_hcal_dev")
    output_file.Close()
    print(f"Saved {output_file.GetName()}")