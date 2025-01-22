#include "cdbHistConv.C"


void genStatusCDBTrees(){

    /*
    const char* inputFile = "files2.txt";
    std::ifstream fileList(inputFile);
    if (!fileList.is_open()) {
        std::cerr << "Error opening file list" << std::endl;
        return;
    }    

    std::string line;
    while (std::getline(fileList, line)) {
      string infile = Form("mergedQA/%s",line.c_str());
      if (!fin) {cout << "no file" << endl; return;}

      // Parse the tag and run number from the filename
      std::string tag;
      std::string runNumber;
      std::size_t pos1 = line.find("HIST_CALO_") + 10;
      std::size_t pos2 = line.find("-", pos1);
      std::size_t pos3 = line.find(".root");

      if (pos1 != std::string::npos && pos2 != std::string::npos && pos3 != std::string::npos) {
          tag = line.substr(pos1, pos2 - pos1);
          runNumber = line.substr(pos2 + 1, pos3 - pos2 - 1);
      } else {
          std::cerr << "Error parsing filename: " << line << std::endl;
          continue;
      }
      */

      TFile* fin = new TFile("hist_new_calib_12_12_24_trig10_events_ana450_2024p009_54912.root");
      std::string tag = "ana450_2024p009";
      std::string runNumber = "54912";

      cout << "doing " << tag.c_str() << "  " << runNumber.c_str() << endl;
      string payloadName;

      //string path = "cdbStatusFiles/";

      // CEMC
      string detector = "CEMC";
      TH2F* h1 = (TH2F*) fin->Get("h_CaloFittingQA_cemc_etaphi_ZScrosscalib");
      if (h1) {
        string payloadName = detector + "_ZSCrossCalib"+ "_" +tag + "_" + runNumber +".root"; 
        histToCaloCDBTree(payloadName, "ratio", 0, h1);
      }
      

      // HCALIN 
      detector = "HCALIN";
      h1 = (TH2F*) fin->Get("h_CaloFittingQA_ihcal_etaphi_ZScrosscalib");
      if (h1) {
        payloadName = detector + "_ZSCrossCalib"+ "_" + tag + "_" + runNumber +".root";
        histToCaloCDBTree(payloadName, "ratio", 1, h1);
      }

      // HCALOUT 
      detector = "HCALOUT";
      h1 = (TH2F*) fin->Get("h_CaloFittingQA_ohcal_etaphi_ZScrosscalib");
      if (h1) {
        payloadName = detector + "_ZSCrossCalib"+ "_" + tag + "_" + runNumber +".root";
        histToCaloCDBTree(payloadName, "ratio", 1, h1);
      } 

      fin->Close();

}
