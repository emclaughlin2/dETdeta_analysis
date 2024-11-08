#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <cdbobjects/CDBTTree.h> 
#include "TowerInfoDefs.h"

R__LOAD_LIBRARY(libcdbobjects)
R__LOAD_LIBRARY(libcalo_io.so)

void mc_decalib_factors(){

    float em_mc_decalib[96][256] = {{0.0}};
    float ih_mc_decalib[24][64] = {{0.0}};
    float oh_mc_decalib[24][64] = {{0.0}};

    TFile *emfile = TFile::Open("/sphenix/u/bseidlitz/work/macros/calibrations/calo/emcal_calib_year1/fin23714_23746/calib_emcal_23726_23746.root");
    TTree *emtree = (TTree*)emfile->Get("Multiple");

    Int_t IID;
    Float_t FFemc_datadriven_qm1_correction;
    TBranch *b_IID = emtree->GetBranch("IID");
    TBranch *b_FFemc_datadriven_qm1_correction = emtree->GetBranch("FFemc_datadriven_qm1_correction");
    b_IID->SetAddress(&IID);
    b_FFemc_datadriven_qm1_correction->SetAddress(&FFemc_datadriven_qm1_correction);

    Long64_t nEntries = emtree->GetEntries();
    std::cout << nEntries << std::endl;
    for (Long64_t i = 0; i < nEntries; ++i) {
      emtree->GetEntry(i);
      int ieta = TowerInfoDefs::getCaloTowerEtaBin(IID);
      int iphi = TowerInfoDefs::getCaloTowerPhiBin(IID);
      em_mc_decalib[ieta][iphi] = FFemc_datadriven_qm1_correction;
      if (i == 1388) std::cout << ieta << " " << iphi << " " << FFemc_datadriven_qm1_correction << std::endl;
    }
    emfile->Close();

    TFile *ihfile = TFile::Open("/sphenix/u/bseidlitz/work/forChris/calibHCal_apr7/ihcal_cdb_calib.root");
    TTree *ihtree = (TTree*)ihfile->Get("Multiple");

    Int_t ihIID;
    Float_t Fihcal_abscalib_mip;
    TBranch *b_ihIID = ihtree->GetBranch("IID");
    TBranch *b_Fihcal_abscalib_mip = ihtree->GetBranch("Fihcal_abscalib_mip");
    b_ihIID->SetAddress(&ihIID);
    b_Fihcal_abscalib_mip->SetAddress(&Fihcal_abscalib_mip);

    nEntries = ihtree->GetEntries();
    std::cout << nEntries << std::endl;
    for (Long64_t i = 0; i < nEntries; ++i) {
      ihtree->GetEntry(i);
      int ieta = TowerInfoDefs::getCaloTowerEtaBin(ihIID);
      int iphi = TowerInfoDefs::getCaloTowerPhiBin(ihIID);
      ih_mc_decalib[ieta][iphi] = Fihcal_abscalib_mip;
      if (i == 1388) std::cout << ieta << " " << iphi << " " << Fihcal_abscalib_mip << std::endl;
    }
    ihfile->Close();

    TFile *ohfile = TFile::Open("/sphenix/u/bseidlitz/work/forChris/calibHCal_apr7/ohcal_cdb_calib.root");
    TTree *ohtree = (TTree*)ohfile->Get("Multiple");

    Int_t ohIID;
    Float_t Fohcal_abscalib_mip;
    TBranch *b_ohIID = ohtree->GetBranch("IID");
    TBranch *b_Fohcal_abscalib_mip = ohtree->GetBranch("Fohcal_abscalib_mip");
    b_ohIID->SetAddress(&ohIID);
    b_Fohcal_abscalib_mip->SetAddress(&Fohcal_abscalib_mip);

    nEntries = ohtree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
      ohtree->GetEntry(i);
      int ieta = TowerInfoDefs::getCaloTowerEtaBin(ohIID);
      int iphi = TowerInfoDefs::getCaloTowerPhiBin(ohIID);
      oh_mc_decalib[ieta][iphi] = Fohcal_abscalib_mip;
      if (i == 1388) std::cout << ieta << " " << iphi << " " << Fohcal_abscalib_mip << std::endl;
    }
    ohfile->Close();

    for (int i = 0; i < 24; i++) {
      for (int j = 0; j < 64; j++) {
        std::cout << oh_mc_decalib[i][j] << " ";
      }
    }
    std::cout << std::endl;

    TFile *outfile = new TFile("/sphenix/user/egm2153/calib_study/detdeta/analysis/mc_decalib_factors.root", "RECREATE");
    TTree *tree = new TTree("T", "");
    tree->Branch("em_mc_decalib", em_mc_decalib, "em_mc_decalib[96][256]/F");
    tree->Branch("ih_mc_decalib", ih_mc_decalib, "ih_mc_decalib[24][64]/F");
    tree->Branch("oh_mc_decalib", oh_mc_decalib, "oh_mc_decalib[24][64]/F");
    tree->Fill();
    Int_t write_status = tree->Write();
    if (write_status < 0) {
        std::cerr << "Error writing tree to file: " << write_status << std::endl;
    }

    // Close the output file
    outfile->Close();
    delete outfile;
}
