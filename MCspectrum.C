#include <iostream>
#include <string>
#include <vector>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TDatabasePDG.h>
#include <TMath.h>

#include <cmath>
#include <math.h>

const double PI_M = 3.14159265358979323846;

// make TChain
TChain *chain = new TChain("ttree");

/*
// HIJING/AMPT Npart centrality bounds
int bound10 = 275;
int bound20 = 194;
int bound40 = 89;
int bound60 = 32;
int bound92 = 3;
*/

/*
// EPOS Npart centrality bounds
int bound10 = 278;
int bound20 = 201;
int bound40 = 97;
int bound60 = 37;
int bound92 = 3;
*/

// AMPT MBD charge centrality bounds
int bound10 = 142123;
int bound20 = 100250;
int bound40 = 47382;
int bound60 = 19223;
int bound92 = 2228.0;

const int ncentbins = 5;

double ptbins[27] = {0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 4.5, 5, 5.5, 6};
double kaonptbins[23] = {0.5, 0.6, 0.7, 0.8, 0.9, 1,
                         1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
                         1.7, 1.8, 1.9, 2, 2.2, 2.4,
                         2.6, 2.8, 3, 3.5, 4};
int findcentbin(int ipart)
{
    if (ipart > bound10)
    {
        return 0;
    }
    else if (ipart > bound20 && ipart <= bound10)
    {
        return 1;
    }
    else if (ipart > bound40 && ipart <= bound20)
    {
        return 2;
    }
    else if (ipart > bound60 && ipart <= bound40)
    {
        return 3;
    }
    else if (ipart > bound92 && ipart <= bound60)
    {
        return 4;
    }
    else
    {
        return -1;
    }
}

TH1F *h_pimi0010_ratio;
TH1F *h_pimi1020_ratio;
TH1F *h_pimi2040_ratio;
TH1F *h_pimi4060_ratio;
TH1F *h_pimi6092_ratio;
TH1F *h_pip0010_ratio;
TH1F *h_pip1020_ratio;
TH1F *h_pip2040_ratio;
TH1F *h_pip4060_ratio;
TH1F *h_pip6092_ratio;
TH1F *h_p0010_ratio;
TH1F *h_p1020_ratio;
TH1F *h_p2040_ratio;
TH1F *h_p4060_ratio;
TH1F *h_p6092_ratio;
TH1F *h_pbar0010_ratio;
TH1F *h_pbar1020_ratio;
TH1F *h_pbar2040_ratio;
TH1F *h_pbar4060_ratio;
TH1F *h_pbar6092_ratio;
TH1F *h_kp0010_ratio;
TH1F *h_kp1020_ratio;
TH1F *h_kp2040_ratio;
TH1F *h_kp4060_ratio;
TH1F *h_kp6092_ratio;
TH1F *h_kmi0010_ratio;
TH1F *h_kmi1020_ratio;
TH1F *h_kmi2040_ratio;
TH1F *h_kmi4060_ratio;
TH1F *h_kmi6092_ratio;

float avgcent[ncentbins] = {325.8, 236.1, 141.5, 61.6, 14.7};

float findcorrection(int npart, int pid, float pt)
{

    // if not proton anti-proton neutron anti-neutron pi+ pi- K+ K- K0 K0bar K long K short omega return 1
    // if (abs(pid) != 2212 && abs(pid) != 2112 && abs(pid) != 211 && abs(pid) != 321 abs(pid) != 311 && pid != 111 && pid != 130 && pid != 310 && abs(pid) > 4000 && abs(pid) < 2000)
    // {
    //     return 1;
    // }
    // if(pid==321) std::cout<<"kplus"<<std::endl;
    float weight[ncentbins] = {0};
    float scale = 0;

    if (npart > avgcent[0] || npart < avgcent[ncentbins - 1])
    {
        if (npart > avgcent[0])
            weight[0] = 1;
        if (npart < avgcent[4])
            weight[4] = 1;
    }
    else
    {
        // use interpolation here
        // first find which two bins the npart falls in between
        int lowerBin = -1;
        int upperBin = -1;

        for (int i = 0; i < ncentbins - 1; i++)
        {
            if (npart <= avgcent[i] && npart >= avgcent[i + 1])
            {
                lowerBin = i;
                upperBin = i + 1;
                break;
            }
        }
        // interpolate
        weight[upperBin] = (avgcent[lowerBin] - npart) / (avgcent[lowerBin] - avgcent[upperBin]);
        weight[lowerBin] = (npart - avgcent[upperBin]) / (avgcent[lowerBin] - avgcent[upperBin]);
        // print all weights and the sum
    }

    if (pid == 211)
    {
        // loop over cent bins
        for (int i = 0; i < ncentbins; i++)
        {
            if (weight[i] == 0)
                continue;
            if (i == 0)
            {
                scale += weight[i] * h_pip0010_ratio->Interpolate(pt);
            }
            else if (i == 1)
            {
                scale += weight[i] * h_pip1020_ratio->Interpolate(pt);
            }
            else if (i == 2)
            {
                scale += weight[i] * h_pip2040_ratio->Interpolate(pt);
            }
            else if (i == 3)
            {
                scale += weight[i] * h_pip4060_ratio->Interpolate(pt);
            }
            else if (i == 4)
            {
                scale += weight[i] * h_pip6092_ratio->Interpolate(pt);
            }
        }
    }
    else if (pid == -211 || pid == 111)
    {
        // loop over cent bins
        for (int i = 0; i < ncentbins; i++)
        {
            if (weight[i] == 0)
                continue;
            if (i == 0)
            {
                scale += weight[i] * h_pimi0010_ratio->Interpolate(pt);
            }

            else if (i == 1)
            {
                scale += weight[i] * h_pimi1020_ratio->Interpolate(pt);
            }
            else if (i == 2)
            {
                scale += weight[i] * h_pimi2040_ratio->Interpolate(pt);
            }
            else if (i == 3)
            {
                scale += weight[i] * h_pimi4060_ratio->Interpolate(pt);
            }
            else if (i == 4)
            {
                scale += weight[i] * h_pimi6092_ratio->Interpolate(pt);
            }
        }
    }
    else if (pid == 2212 || pid == 2112)
    {
        
        // loop over cent bins
        for (int i = 0; i < ncentbins; i++)
        {
            if (weight[i] == 0)
                continue;
            if (i == 0)
            {
                scale += weight[i] * h_p0010_ratio->Interpolate(pt);
            }
            else if (i == 1)
            {
                scale += weight[i] * h_p1020_ratio->Interpolate(pt);
            }
            else if (i == 2)
            {
                scale += weight[i] * h_p2040_ratio->Interpolate(pt);
            }
            else if (i == 3)
            {
                scale += weight[i] * h_p4060_ratio->Interpolate(pt);
            }
            else if (i == 4)
            {
                scale += weight[i] * h_p6092_ratio->Interpolate(pt);
            }
        }
    }
    else if (abs(pid) > 2000 && abs(pid) < 4000)
    {
        
        // loop over cent bins
        for (int i = 0; i < ncentbins; i++)
        {
            if (weight[i] == 0)
                continue;
            if (i == 0)
            {
                scale += weight[i] * h_pbar0010_ratio->Interpolate(pt);
            }
            else if (i == 1)
            {
                scale += weight[i] * h_pbar1020_ratio->Interpolate(pt);
            }
            else if (i == 2)
            {
                scale += weight[i] * h_pbar2040_ratio->Interpolate(pt);
            }
            else if (i == 3)
            {
                scale += weight[i] * h_pbar4060_ratio->Interpolate(pt);
            }
            else if (i == 4)
            {
                scale += weight[i] * h_pbar6092_ratio->Interpolate(pt);
            }
        }
    }
    // k+
    else if (pid == 321 || pid == 130 || pid == 310 || pid == 311)
    {
        // loop over cent bins
        for (int i = 0; i < ncentbins; i++)
        {
            if (weight[i] == 0)
                continue;
            if (i == 0)
            {
                scale += weight[i] * h_kp0010_ratio->Interpolate(pt);
            }
            else if (i == 1)
            {
                scale += weight[i] * h_kp1020_ratio->Interpolate(pt);
            }
            else if (i == 2)
            {
                scale += weight[i] * h_kp2040_ratio->Interpolate(pt);
            }
            else if (i == 3)
            {
                scale += weight[i] * h_kp4060_ratio->Interpolate(pt);
            }
            else if (i == 4)
            {
                scale += weight[i] * h_kp6092_ratio->Interpolate(pt);
            }
        }
    }
    // k-
    else if (pid == -321)
    {
        // loop over cent bins
        for (int i = 0; i < ncentbins; i++)
        {
            if (weight[i] == 0)
                continue;
            if (i == 0)
            {
                scale += weight[i] * h_kmi0010_ratio->Interpolate(pt);
            }
            else if (i == 1)
            {
                scale += weight[i] * h_kmi1020_ratio->Interpolate(pt);
            }
            else if (i == 2)
            {
                scale += weight[i] * h_kmi2040_ratio->Interpolate(pt);
            }
            else if (i == 3)
            {
                scale += weight[i] * h_kmi4060_ratio->Interpolate(pt);
            }
            else if (i == 4)
            {
                scale += weight[i] * h_kmi6092_ratio->Interpolate(pt);
            }
        }
    }
    else
    {
        scale = 1;
    }
    // if k+ print pt npart and scale
    if(pid==22 && scale!=1)
        std::cout <<pid<<" "<< pt << " " << npart << " " << scale << std::endl;
    return scale;
}

void addFilesToChain(TChain *chain, const std::string &path)
{
    void *dirp = gSystem->OpenDirectory(path.c_str());
    const char *dirname;
    while ((dirname = gSystem->GetDirEntry(dirp)))
    {
        if (TString(dirname).Contains("OutDir"))
        {
            std::string fullPath = path + "/" + dirname + "/ampt_sim_output.root";
            chain->Add(fullPath.c_str());
        }
    }
}

void MCspectrum()
{

    addFilesToChain(chain, "/sphenix/user/egm2153/calib_study/detdeta/amptrun/condor");
    /*
    const char* inputDirectory = "/sphenix/user/egm2153/calib_study/detdeta/runsimana0/output/evt/";
        for (int i = 0; i < 555; i++) {
            TString wildcardPath = TString::Format("%sevents_20231122_nopileup_mc_cor_%d.root", inputDirectory, i);
            chain->Add(wildcardPath);
    }
    */

    // read data from tree
    int m_hepmc = -999;
    int m_npart_proj = -999;
    int m_npart_targ = 0;
    int m_ncoll = -999;
    int m_ncoll_hard = -999;
    int m_hepmc_pid[20000];
    float m_hepmc_pt[20000];
    float m_hepmc_eta[20000];
    float m_hepmc_phi[20000];
    float m_hepmc_pz[20000];
    float m_hepmc_e[20000];
    int m_sectormb = 0;
    float m_mbenergy[256];


    chain->SetBranchAddress("truthpar_nh", &m_hepmc);
    chain->SetBranchAddress("npart", &m_npart_proj);
    //chain->SetBranchAddress("m_npart_targ", &m_npart_targ);
    chain->SetBranchAddress("ncoll", &m_ncoll);
    //chain->SetBranchAddress("m_ncoll_hard", &m_ncoll_hard);
    chain->SetBranchAddress("truthparh_id", m_hepmc_pid);
    chain->SetBranchAddress("truthparh_pt", m_hepmc_pt);
    chain->SetBranchAddress("truthparh_eta", m_hepmc_eta);
    chain->SetBranchAddress("truthparh_pz", m_hepmc_pz);
    chain->SetBranchAddress("truthparh_phi", m_hepmc_phi);
    chain->SetBranchAddress("truthparh_e", m_hepmc_e);
    chain->SetBranchAddress("sectormb", &m_sectormb);
    chain->SetBranchAddress("mbenrgy", m_mbenergy);

    // read in the upweight factor
    TFile *f_upweight = new TFile("/sphenix/user/shuhangli/dETdeta/macro/scalecompare.root");
    h_pimi0010_ratio = (TH1F *)f_upweight->Get("h_pimi0010_ratio");
    h_pimi1020_ratio = (TH1F *)f_upweight->Get("h_pimi1020_ratio");
    h_pimi2040_ratio = (TH1F *)f_upweight->Get("h_pimi2040_ratio");
    h_pimi4060_ratio = (TH1F *)f_upweight->Get("h_pimi4060_ratio");
    h_pimi6092_ratio = (TH1F *)f_upweight->Get("h_pimi6092_ratio");

    h_pip0010_ratio = (TH1F *)f_upweight->Get("h_pip0010_ratio");
    h_pip1020_ratio = (TH1F *)f_upweight->Get("h_pip1020_ratio");
    h_pip2040_ratio = (TH1F *)f_upweight->Get("h_pip2040_ratio");
    h_pip4060_ratio = (TH1F *)f_upweight->Get("h_pip4060_ratio");
    h_pip6092_ratio = (TH1F *)f_upweight->Get("h_pip6092_ratio");

    h_p0010_ratio = (TH1F *)f_upweight->Get("h_p0010_ratio");
    h_p1020_ratio = (TH1F *)f_upweight->Get("h_p1020_ratio");
    h_p2040_ratio = (TH1F *)f_upweight->Get("h_p2040_ratio");
    h_p4060_ratio = (TH1F *)f_upweight->Get("h_p4060_ratio");
    h_p6092_ratio = (TH1F *)f_upweight->Get("h_p6092_ratio");

    h_pbar0010_ratio = (TH1F *)f_upweight->Get("h_pbar0010_ratio");
    h_pbar1020_ratio = (TH1F *)f_upweight->Get("h_pbar1020_ratio");
    h_pbar2040_ratio = (TH1F *)f_upweight->Get("h_pbar2040_ratio");
    h_pbar4060_ratio = (TH1F *)f_upweight->Get("h_pbar4060_ratio");
    h_pbar6092_ratio = (TH1F *)f_upweight->Get("h_pbar6092_ratio");

    h_kp0010_ratio = (TH1F *)f_upweight->Get("h_kp0010_ratio");
    h_kp1020_ratio = (TH1F *)f_upweight->Get("h_kp1020_ratio");
    h_kp2040_ratio = (TH1F *)f_upweight->Get("h_kp2040_ratio");
    h_kp4060_ratio = (TH1F *)f_upweight->Get("h_kp4060_ratio");
    h_kp6092_ratio = (TH1F *)f_upweight->Get("h_kp6092_ratio");

    h_kmi0010_ratio = (TH1F *)f_upweight->Get("h_kmi0010_ratio");
    h_kmi1020_ratio = (TH1F *)f_upweight->Get("h_kmi1020_ratio");
    h_kmi2040_ratio = (TH1F *)f_upweight->Get("h_kmi2040_ratio");
    h_kmi4060_ratio = (TH1F *)f_upweight->Get("h_kmi4060_ratio");
    h_kmi6092_ratio = (TH1F *)f_upweight->Get("h_kmi6092_ratio");

    TFile *fout = new TFile("AMPTspectrum.root", "RECREATE");
    // make histograms
    TH1F *hnpart = new TH1F("hnpart", "hnpart", 400, 0, 400);
    TH1F *hprotonpt[ncentbins];
    TH1F *hprotonbarpt[ncentbins];
    TH1F *hpipluspt[ncentbins];
    TH1F *hpiminuspt[ncentbins];
    TH1F *hkpluspt[ncentbins];
    TH1F *hkminuspt[ncentbins];
    // calculate the rapidity range for each pt bin for eta<abs(0.35)
    float protonweight[26];
    for (int i = 0; i < 26; i++)
    {
        float pt = (ptbins[i] + ptbins[i + 1]) / 2;
        float eta = 0.35;
        float m = 0.938272;
        float y = log((sqrt(pt * pt * cosh(eta) * cosh(eta) + m * m) + pt * sinh(eta)) / sqrt(pt * pt + m * m));
        protonweight[i] = 2 * y * (ptbins[i + 1] - ptbins[i]);
        std::cout << eta << " " << pt << " " << y << " " << protonweight[i] << std::endl;
    }
    float pionweight[26];
    for (int i = 0; i < 26; i++)
    {
        float pt = (ptbins[i] + ptbins[i + 1]) / 2;
        float eta = 0.35;
        float m = 0.139570;
        float y = log((sqrt(pt * pt * cosh(eta) * cosh(eta) + m * m) + pt * sinh(eta)) / sqrt(pt * pt + m * m));
        pionweight[i] = 2 * y * (ptbins[i + 1] - ptbins[i]);
        std::cout << eta << " " << pt << " " << y << " " << pionweight[i] << std::endl;
    }
    float kaonweight[22];
    for (int i = 0; i < 22; i++)
    {
        float pt = (ptbins[i] + ptbins[i + 1]) / 2;
        float eta = 0.35;
        float m = 0.493677;
        float y = log((sqrt(pt * pt * cosh(eta) * cosh(eta) + m * m) + pt * sinh(eta)) / sqrt(pt * pt + m * m));
        kaonweight[i] = 2 * y * (ptbins[i + 1] - ptbins[i]);
        std::cout << eta << " " << pt << " " << y << " " << kaonweight[i] << std::endl;
    }
    TH1F *hET[ncentbins];
    for (int i = 0; i < ncentbins; i++)
    {
        hprotonpt[i] = new TH1F(Form("hprotonpt_%d", i), Form("hprotonpt_%d", i), 26, ptbins);
        hprotonbarpt[i] = new TH1F(Form("hprotonbarpt_%d", i), Form("hprotonbarpt_%d", i), 26, ptbins);
        hpipluspt[i] = new TH1F(Form("hpipluspt_%d", i), Form("hpipluspt_%d", i), 26, ptbins);
        hpiminuspt[i] = new TH1F(Form("hpiminuspt_%d", i), Form("hpiminuspt_%d", i), 26, ptbins);
        hkpluspt[i] = new TH1F(Form("hkpluspt_%d", i), Form("hkpluspt_%d", i), 26, ptbins);
        hkminuspt[i] = new TH1F(Form("hkminuspt_%d", i), Form("hkminuspt_%d", i), 26, ptbins);
        //hkpluspt[i] = new TH1F(Form("hkpluspt_%d", i), Form("hkpluspt_%d", i), 22, kaonptbins);
        //hkminuspt[i] = new TH1F(Form("hkminuspt_%d", i), Form("hkminuspt_%d", i), 22, kaonptbins);
        hET[i] = new TH1F(Form("hET_%d", i), Form("hET_%d", i), 24, -1.1, 1.1);
    }

    TProfile *hdetdeta = new TProfile("hdetdeta", "hdetdeta", 400, 0, 400);
    TProfile *hdetdeta_noscale = new TProfile("hdetdeta_noscale", "hdetdeta_noscale", 400, 0, 400);

    // loop over events
    int nentries = chain->GetEntries();
    int nevents[ncentbins] = {0};
    for (int i = 0; i < nentries; i++)
    {
        if (i % 10000 == 0)
            std::cout << "Event " << i << " / " << nentries << std::endl;
        chain->GetEntry(i);
        hnpart->Fill(m_npart_proj + m_npart_targ);
        float mbd_total = 0;
        for (int j = 0; j < m_sectormb; j++) {
            mbd_total += m_mbenergy[j];
        }

        // find centrality bin
        //int centbin = findcentbin(m_npart_proj + m_npart_targ); // for cent from Npart
        int centbin = findcentbin(mbd_total); // for cent from MBD charge 
        if (centbin == -1)
        {
            continue;
        }
        nevents[centbin]++;
        float ET = 0;
        float ET_noscale = 0;
        // loop over particles
        for (int j = 0; j < m_hepmc; j++)
        {
            // if proton protonbar pi+ pi-
            // if (m_hepmc_pid[j] == 2212 || m_hepmc_pid[j] == -2212 || m_hepmc_pid[j] == 211 || m_hepmc_pid[j] == -211)
            {

                float pz = m_hepmc_pz[j];
                // calcualte eta using cmath
                float eta = m_hepmc_eta[j];

                float pt = m_hepmc_pt[j];
                float e = m_hepmc_e[j];
                float e_noscale = m_hepmc_e[j];
                float theta = atan2(pt, pz);
                int pid = m_hepmc_pid[j];
                float scale = 1;
                

                // if (pid == 211 || pid == -211 || pid == 2212 || pid == -2212 || pid == 111)
                
                {
                    if (pid == 111)
                    {
                        if (gRandom->Rndm() > 0.5)
                        {
                            scale = findcorrection(m_npart_proj + m_npart_targ, 211, pt);
                        }
                        else
                        {
                            scale = findcorrection(m_npart_proj + m_npart_targ, -211, pt);
                        }
                    }
                    else
                    {
                        scale = findcorrection(m_npart_proj + m_npart_targ, pid, pt);
                    }
                    //use 1 as scale
                    scale = 1;
                    e = e * scale;
                }
                hET[centbin]->Fill(eta, e * abs(sin(theta)));
                if (eta > 0 && eta < 1)
                {
                    ET += e * sin(theta);
                    ET_noscale += e_noscale * sin(theta);
                }

                if (abs(eta) > 0.35)
                {
                    continue;
                }
                if (m_hepmc_pid[j] == 2212)
                {
                    hprotonpt[centbin]->Fill(pt, scale / (2 * PI_M * pt));
                    // std::cout<<"at event "<<i<<" "<<px <<" "<<py<<" "<<pz<<" "<<pt<<" "<<eta<<std::endl;
                    //  hprotonpt[centbin]->Fill(pt);
                    //  std::cout<<pt<< " "<<centbin<<std::endl;
                }
                else if (m_hepmc_pid[j] == -2212)
                {
                    hprotonbarpt[centbin]->Fill(pt, scale / (2 * PI_M * pt));
                }
                else if (m_hepmc_pid[j] == 211)
                {
                    hpipluspt[centbin]->Fill(pt, scale / (2 * PI_M * pt));
                }
                else if (m_hepmc_pid[j] == -211)
                {
                    hpiminuspt[centbin]->Fill(pt, scale / (2 * PI_M * pt));
                }
                else if (m_hepmc_pid[j] == 321)
                {
                    //std::cout << pt << " " << scale << std::endl;
                    hkpluspt[centbin]->Fill(pt, scale / (2 * PI_M * pt));
                }
                else if (m_hepmc_pid[j] == -321)
                {
                    hkminuspt[centbin]->Fill(pt, scale / (2 * PI_M * pt));
                }
            }
        }
        hdetdeta->Fill(m_npart_proj + m_npart_targ, ET / (0.5 * (m_npart_proj + m_npart_targ)));
        hdetdeta_noscale->Fill(m_npart_proj + m_npart_targ, ET_noscale / (0.5 * (m_npart_proj + m_npart_targ)));
    }
    // scale the pt spectra
    for (int i = 0; i < ncentbins; i++)
    {

        //
        std::cout << nevents[i] << std::endl;
        hprotonpt[i]->Scale(1. / nevents[i]);
        hprotonbarpt[i]->Scale(1. / nevents[i]);
        hpipluspt[i]->Scale(1. / nevents[i]);
        hpiminuspt[i]->Scale(1. / nevents[i]);
        hkpluspt[i]->Scale(1. / nevents[i]);
        hkminuspt[i]->Scale(1. / nevents[i]);
        hET[i]->Scale(1. / nevents[i]);
        hET[i]->Scale(1. / hET[i]->GetBinWidth(1));
        
        for (int j = 1; j <= 26; j++)
        {
            hprotonpt[i]->SetBinContent(j, hprotonpt[i]->GetBinContent(j) / protonweight[j - 1]);
            hprotonbarpt[i]->SetBinContent(j, hprotonbarpt[i]->GetBinContent(j) / protonweight[j - 1]);
            hpipluspt[i]->SetBinContent(j, hpipluspt[i]->GetBinContent(j) / pionweight[j - 1]);
            hpiminuspt[i]->SetBinContent(j, hpiminuspt[i]->GetBinContent(j) / pionweight[j - 1]);

            // set bin error
            hprotonpt[i]->SetBinError(j, hprotonpt[i]->GetBinError(j) / protonweight[j - 1]);
            hprotonbarpt[i]->SetBinError(j, hprotonbarpt[i]->GetBinError(j) / protonweight[j - 1]);
            hpipluspt[i]->SetBinError(j, hpipluspt[i]->GetBinError(j) / pionweight[j - 1]);
            hpiminuspt[i]->SetBinError(j, hpiminuspt[i]->GetBinError(j) / pionweight[j - 1]);
        }
        for (int j = 1; j <= 22; j++)
        {
            hkpluspt[i]->SetBinContent(j, hkpluspt[i]->GetBinContent(j) / kaonweight[j - 1]);
            hkminuspt[i]->SetBinContent(j, hkminuspt[i]->GetBinContent(j) / kaonweight[j - 1]);
            hkpluspt[i]->SetBinError(j, hkpluspt[i]->GetBinError(j) / kaonweight[j - 1]);
            hkminuspt[i]->SetBinError(j, hkminuspt[i]->GetBinError(j) / kaonweight[j - 1]);
        }
        
    }

    fout->Write();
    fout->Close();
}