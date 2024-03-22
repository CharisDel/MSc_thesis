// For more details see https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/ZZjj2l2vRun2#Data_and_MC_samples
// Script for plotting 15/06/20 Dimitrii Krasnopevtsev
#include <TRandom3.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include "TProfile.h"
#include "TMath.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TLatex.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
using namespace std;

// DoReco function to calculate events from trees

std::vector<float> DoReco(TFile *file, int MC, string directory, string fileName, TH1F *Histogram, int Harry)
{

    TH1::SetDefaultSumw2(kTRUE);

    TTree *tree = file->Get<TTree>("tree");

    Int_t nentries = (Int_t)tree->GetEntries();

    Double_t event_3CR = 0.;
    Double_t weight = 1.;

    Double_t M2Lep = 0.;
    Double_t met_tst = 0.;
    Double_t met_signif = 0.;
    Double_t dMetZPhi = 0.;
    Double_t MetOHT = 0.;
    Double_t dLepR = 0.;
    Double_t n_bjets = 0.;
    Double_t n_jets = 0.;
    Double_t leading_pT_lepton = 0;
    Double_t subleading_pT_lepton = 0;
    Double_t detajj = 0.;
    Double_t mjj = 0.;
    Double_t leading_jet_pt = 0.;
    Double_t second_jet_pt = 0.;
    Double_t event_type = 0.;
    Double_t Z_pT = 0.;
    Double_t ST = 0.;
    Double_t STjj = 0.;

    int counter = 0;
    //  Double_t mTW = 0.;

    std::vector<float> events;
    events.clear();

    vector<TString> branches = {
        //"M2Lep",
        "met_tst",
        "dMetZPhi",
        "MetOHT",
        //"dLepR",
        //"leading_pT_lepton",
        //"subleading_pT_lepton",
        //"Z_pT",
        //"met_signif",
        "n_jets",
        //"n_bjets",
        //"detajj",
        //"mjj",
        //"leading_jet_pt",
        //"second_jet_pt",
        //"event_3CR",
        //"event_type",
        "global_weight"};

    tree->SetBranchStatus("*", 0);

    for (const auto &branch : branches)
    {
        tree->SetBranchStatus(branch, 1);
    }

    // tree->SetBranchAddress("event_3CR", &event_3CR);
    // tree->SetBranchAddress("M2Lep", &M2Lep);
    tree->SetBranchAddress("met_tst", &met_tst);
    // tree->SetBranchAddress("met_signif", &met_signif);
    tree->SetBranchAddress("dMetZPhi", &dMetZPhi);
    tree->SetBranchAddress("MetOHT", &MetOHT);
    // tree->SetBranchAddress("dLepR", &dLepR);
    // tree->SetBranchAddress("leading_pT_lepton", &leading_pT_lepton);
    // tree->SetBranchAddress("subleading_pT_lepton", &subleading_pT_lepton);
    tree->SetBranchAddress("n_jets", &n_jets);
    // tree->SetBranchAddress("n_bjets", &n_bjets);
    // tree->SetBranchAddress("detajj", &detajj);
    // tree->SetBranchAddress("mjj", &mjj);
    // tree->SetBranchAddress("leading_jet_pt", &leading_jet_pt);
    // tree->SetBranchAddress("second_jet_pt", &second_jet_pt);
    // tree->SetBranchAddress("event_type", &event_type);
    // tree->SetBranchAddress("Z_pT", &Z_pT);

    double signal = 0.;
    double signaler = 0.;
    double signal_njets_lt_1 = 0.;
    double signal_njets_eq_1 = 0.;
    double signal_njets_gt_1 = 0.;
    double signaler_njets_lt_1 = 0.;
    double signaler_njets_eq_1 = 0.;
    double signaler_njets_gt_1 = 0.;
    double signal_njets_less_1 = 0.;
    double signal_njets_just_1 = 0.;
    double signal_njets_big_1 = 0.;
    double signaler_njets_less_1 = 0.;
    double signaler_njets_just_1 = 0.;
    double signaler_njets_big_1 = 0.;

    // give each event the appropriate weight (data events should have weight = 1)
    // global_weight=weight_pileup*weight_gen*weight_exp*weight_trig*weight_jets*weight_jvt + normalization to the MC cross section and L=139 fb^-1
    if (MC == 1)
    {
        tree->SetBranchAddress("global_weight", &weight);
    }

    // Loop over events
    for (Int_t i = 0; i < nentries; i++)
    {

        tree->GetEntry(i);

        if (directory == "/3lCR/" || directory == "/emCR_B/") /*&& (event_3CR == 0 || event_3CR == 1 || event_3CR == 2 || event_3CR == 3 || event_3CR == 4))*/ // an thelw na parw apla gia to 3lCR ta μμμ,μμe,μee,eee.
        {
            signal += weight;            // signal yield is sum of weights
            signaler += weight * weight; // keep Sum(weights^2) for calculating the error on the signal yield
        }

        if (directory == "/Zjets0/" || directory == "/Zjets1/" || directory == "/Zjets2/")
        {
            if ((met_tst < 80) || (met_tst > 80 && MetOHT < 0.45))
            {
                signal += weight;
                signaler += weight * weight;
            }
        }

        if (directory == "/emCR_A/")
        {
            signal += weight;
            signaler += weight * weight;
            if (fileName == "/Z_jets_ee.root")
            {
                if (n_jets < 1 && Harry == 1)
                {
                    signal_njets_lt_1 += weight;
                    signaler_njets_lt_1 += weight * weight;
                }
                if (n_jets == 1 && Harry == 2)
                {
                    signal_njets_eq_1 += weight;
                    signaler_njets_eq_1 += weight * weight;
                }
                if (n_jets > 1 && Harry == 3)
                {
                    signal_njets_gt_1 += weight;
                    signaler_njets_gt_1 += weight * weight;
                }
            }
            if (fileName == "/Z_jets_mumu.root")
            {
                if (n_jets < 1 && Harry == 1)
                {
                    signal_njets_less_1 += weight;
                    signaler_njets_less_1 += weight * weight;
                }

                if (n_jets == 1 && Harry == 2)
                {
                    signal_njets_just_1 += weight;
                    signaler_njets_just_1 += weight * weight;
                }
                if (n_jets > 1 && Harry == 3)
                {
                    signal_njets_big_1 += weight;
                    signaler_njets_big_1 += weight * weight;
                }
            }
        }

        if (directory == "/SR/")
        {
            // Define the cuts using a boolean variable (Reco-cuts)
            bool cuts = (dMetZPhi > 2.536255 && met_tst > 101.130124 &&
                         MetOHT > 0.737131);

            /*bool cuts = (dMetZPhi > 2.5 && met_tst > 101 &&
                         MetOHT > 0.74);*/

            if (cuts)
            {
                signal += weight;
                signaler += weight * weight;
                counter += 1;
                if (fileName == "/Z_jets_ee.root")
                {
                    if (n_jets < 1 && Harry == 1)
                    {
                        // ST = Z_pT + met_tst;
                        // STjj = Z_pT + met_tst + leading_jet_pt + second_jet_pt;
                        signal_njets_lt_1 += weight;
                        signaler_njets_lt_1 += weight * weight;
                        Histogram->Fill(dMetZPhi, weight);
                    }
                    if (n_jets == 1 && Harry == 2)
                    {
                        // ST = Z_pT + met_tst;
                        // STjj = Z_pT + met_tst + leading_jet_pt + second_jet_pt;
                        signal_njets_eq_1 += weight;
                        signaler_njets_eq_1 += weight * weight;
                        Histogram->Fill(dMetZPhi, weight);
                    }
                    if (n_jets > 1 && Harry == 3)
                    {
                        // ST = Z_pT + met_tst;
                        // STjj = Z_pT + met_tst + leading_jet_pt + second_jet_pt;
                        signal_njets_gt_1 += weight;
                        signaler_njets_gt_1 += weight * weight;
                        Histogram->Fill(dMetZPhi, weight);
                    }
                }
                if (fileName == "/Z_jets_mumu.root")
                {
                    if (n_jets < 1 && Harry == 1)
                    {
                        // ST = Z_pT + met_tst;
                        // STjj = Z_pT + met_tst + leading_jet_pt + second_jet_pt;
                        signal_njets_less_1 += weight;
                        signaler_njets_less_1 += weight * weight;
                        Histogram->Fill(dMetZPhi, weight);
                    }

                    if (n_jets == 1 && Harry == 2)
                    {
                        // ST = Z_pT + met_tst;
                        // STjj = Z_pT + met_tst + leading_jet_pt + second_jet_pt;
                        signal_njets_just_1 += weight;
                        signaler_njets_just_1 += weight * weight;
                        Histogram->Fill(dMetZPhi, weight);
                    }
                    if (n_jets > 1 && Harry == 3)
                    {
                        // ST = Z_pT + met_tst;
                        // STjj = Z_pT + met_tst + leading_jet_pt + second_jet_pt;
                        signal_njets_big_1 += weight;
                        signaler_njets_big_1 += weight * weight;
                        Histogram->Fill(dMetZPhi, weight);
                    }
                }
                // ST = Z_pT + met_tst;
                // STjj = Z_pT + met_tst + leading_jet_pt + second_jet_pt;
                if (fileName != "/Z_jets_ee.root" && fileName != "/Z_jets_mumu.root")
                {
                    Histogram->Fill(dMetZPhi, weight);
                }
                // Histogram->Fill(ST,weight);
            }
        }
    }
    cout << "N=" << signal << "+-" << sqrt(signaler) << endl; // signal yield and error on this yield
    printf("ENTRIES FOR EACH HISTOGRAM ARE: %f as well as times got into is : %d \n", Histogram->GetEntries(), counter);

    events.push_back(signal);
    events.push_back(sqrt(signaler));

    events.push_back(signal_njets_lt_1);
    events.push_back(sqrt(signaler_njets_lt_1));
    events.push_back(signal_njets_eq_1);
    events.push_back(sqrt(signaler_njets_eq_1));
    events.push_back(signal_njets_gt_1);
    events.push_back(sqrt(signaler_njets_gt_1));

    events.push_back(signal_njets_less_1);
    events.push_back(sqrt(signaler_njets_less_1));
    events.push_back(signal_njets_just_1);
    events.push_back(sqrt(signaler_njets_just_1));
    events.push_back(signal_njets_big_1);
    events.push_back(sqrt(signaler_njets_big_1));

    return events;
}

std::vector<float> DoReco_FID(TFile *file, string dir_FID)
{

    TH1::SetDefaultSumw2(kTRUE);

    TTree *tree = file->Get<TTree>("tree");

    Int_t nentries = (Int_t)tree->GetEntries();

    std::vector<float> events_FID;
    events_FID.clear();

    double signal_FID = 0.;
    double signaler_FID = 0.;
    Double_t weight_FID = 1.;

    // give each event the appropriate weight (data events should have weight = 1)
    // global_weight=weight_pileup*weight_gen*weight_exp*weight_trig*weight_jets*weight_jvt + normalization to the MC cross section and L=139 fb^-1

    tree->SetBranchAddress("global_weight", &weight_FID);

    // Loop over events
    for (Int_t i = 0; i < nentries; i++)
    {

        tree->GetEntry(i);

        if (dir_FID == "/FIDUCIAL/")
        {
            signal_FID += weight_FID;                // signal yield is sum of weights
            signaler_FID += weight_FID * weight_FID; // keep Sum(weights^2) for calculating the error on the signal yield
        }
    }
    cout << "N=" << signal_FID << "+-" << sqrt(signaler_FID) << endl; // signal yield and error on this yield

    events_FID.push_back(signal_FID);
    events_FID.push_back(sqrt(signaler_FID));

    return events_FID;
}

// Change of directory//
std::vector<std::string> directories = {"/emCR_B/", "/3lCR/", "/Zjets2/", "/Zjets1/", "/Zjets0/", "/emCR_A/", "/SR/"};
std::vector<std::string> fileNames = {"/DATA.root", "/llll.root", "/lllljj.root", "/llqq.root", "/llvv.root", "/llvvjj.root", "/llvvjj_WW.root", "/top.root", "/ttbarV_ttbarVV.root", "/VVV.root", "/W_jets.root", "/Wt.root", "/WW.root", "/WZ.root", "/WZ_jj.root", "/Z_jets_ee.root", "/Z_jets_mumu.root", "/Ztt.root"};

std::vector<std::string> dir = {"/FIDUCIAL/"};
std::vector<std::string> root = {"/EWK_SAMPLE.root", "/QCD_qq_SAMPLE.root", "/QCD_gg_SAMPLE.root"};

// the main method

int vol0()
{
    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    Double_t N_SR;
    Double_t N_SR_error;
    Double_t scale_Signal;
    Double_t scale_Signal_error;

    Double_t scaleWZ = 0.0;
    Double_t scaleWZ_error = 0.0;
    Double_t scaleAllTop = 0.0;
    Double_t scaleAllTop_error = 0.0;
    Double_t scaleWW = 0.0;
    Double_t scaleWW_error = 0.0;
    Double_t scaleZjets0 = 0.0;
    Double_t scaleZjets0_error = 0.0;
    Double_t scaleZjets1 = 0.0;
    Double_t scaleZjets1_error = 0.0;
    Double_t scaleZjets2 = 0.0;
    Double_t scaleZjets2_error = 0.0;
    Double_t sfsignal = 0.0;
    Double_t sfsignal_error = 0.0;
    Double_t asa;

    // Float_t xbins[61] = {-0.5, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 530, 560, 590, 620, 650, 680, 710, 750, 790, 840, 890, 940, 1000, 1100, 1500, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000}; ///// Histograms for the desired variables
    TH1F *HDATA = new TH1F("HDATA", "", 10, 2.2, 3.2);
    TH1F *Hllll = new TH1F("Hllll", "", 10, 2.2, 3.2);
    TH1F *Hlllljj = new TH1F("Hlllljj", "", 10, 2.2, 3.2);
    TH1F *Hllqq = new TH1F("Hllqq", "", 10, 2.2, 3.2);
    TH1F *Hllvv = new TH1F("Hllvv", "", 10, 2.2, 3.2);
    TH1F *Hllvvjj = new TH1F("Hllvvjj", "", 10, 2.2, 3.2);
    TH1F *Hllvvjj_WW = new TH1F("Hllvvjj_WW", "", 10, 2.2, 3.2);
    TH1F *HTop = new TH1F("HTop", "", 10, 2.2, 3.2);
    TH1F *HttbarV_ttbarVV = new TH1F("HttbarV_ttbarVV", "", 10, 2.2, 3.2);
    TH1F *HVVV = new TH1F("HVVV", "", 10, 2.2, 3.2);
    TH1F *HW_jets = new TH1F("HW_jets", "", 10, 2.2, 3.2);
    TH1F *HWt = new TH1F("HWt", "", 10, 2.2, 3.2);
    TH1F *HWW = new TH1F("HWW", "", 10, 2.2, 3.2);
    TH1F *HWZ = new TH1F("HWZ", "", 10, 2.2, 3.2);
    TH1F *HWZ_jj = new TH1F("HWZ_jj", "", 10, 2.2, 3.2);
    TH1F *HZ_jets_ee = new TH1F("HZ_jets_ee", "", 10, 2.2, 3.2);
    TH1F *HZ_jets_ee_0 = new TH1F("HZ_jets_ee_0", "", 10, 2.2, 3.2);
    TH1F *HZ_jets_ee_1 = new TH1F("HZ_jets_ee_1", "", 10, 2.2, 3.2);
    TH1F *HZ_jets_ee_2 = new TH1F("HZ_jets_ee_2", "", 10, 2.2, 3.2);
    TH1F *HZ_jets_mumu = new TH1F("HZ_jets_mumu", "", 10, 2.2, 3.2);
    TH1F *HZ_jets_mumu_0 = new TH1F("HZ_jets_mumu_0", "", 10, 2.2, 3.2);
    TH1F *HZ_jets_mumu_1 = new TH1F("HZ_jets_mumu_1", "", 10, 2.2, 3.2);
    TH1F *HZ_jets_mumu_2 = new TH1F("HZ_jets_mumu_2", "", 10, 2.2, 3.2);
    TH1F *HZtt = new TH1F("HZtt", "", 10, 2.2, 3.2);

    TH1F *HWZscaled = new TH1F("HWZCscaled", "", 10, 2.2, 3.2);
    TH1F *HWWscaled = new TH1F("HWWscaled", "", 10, 2.2, 3.2);
    TH1F *HZjetsscaled = new TH1F("HZjetsscaled", "", 10, 2.2, 3.2);
    TH1F *HTopsscaled = new TH1F("HTopsscaled", "", 10, 2.2, 3.2);

    TH1F *HZjets0 = new TH1F("HZjets0", "", 10, 2.2, 3.2);
    TH1F *HZjets1 = new TH1F("HZjets1", "", 10, 2.2, 3.2);
    TH1F *HZjets2 = new TH1F("HZjets2", "", 10, 2.2, 3.2);
    TH1F *HZjets0scaled = new TH1F("HZjets0scaled", "", 10, 2.2, 3.2);
    TH1F *HZjets1scaled = new TH1F("HZjets1scaled", "", 10, 2.2, 3.2);
    TH1F *HZjets2scaled = new TH1F("HZjets2scaled", "", 10, 2.2, 3.2);

    // Create a THStack object
    // THStack *stack = new THStack("stack", "Stacked Histogram");

    TH1F *HOther = new TH1F("HOther", "", 10, 2.2, 3.2);
    TH1F *HZjets = new TH1F("HZjets", "", 10, 2.2, 3.2);
    TH1F *HSignal = new TH1F("HSignal", "", 10, 2.2, 3.2);
    TH1F *HTops = new TH1F("HTops", "", 10, 2.2, 3.2);
    TH1F *HBKG = new TH1F("HBKG", "", 10, 2.2, 3.2);
    TH1F *HSignif = new TH1F("HSignif", "", 10, 2.2, 3.2);
    TH1F *HSignif_norm = new TH1F("HSignif_norm", "", 10, 2.2, 3.2);

    TH1F *HMC = new TH1F("HMC", "", 10, 2.2, 3.2);

    // for significance /////
    double fourl, fourl_er, fourljj, fourljj_er, doublelqq, doublelqq_er, doublelvv, doublelvv_er, doublelvvjj, doublelvvjj_er, doublelvvjjWW, doublelvvjjWW_er, TOP, TOP_er, ttbar, ttbar_er, tripleV, tripleV_er, WJ, WJ_er, WT, WT_er, ww, ww_er, wz, wz_er, WZJJ, WZJJ_er, ZJEE, ZJEE_er, ZJEE_1, ZJEE_1_er, ZJEE_2, ZJEE_2_er, ZJMM, ZJMM_er, ZJMM_1, ZJMM_1_er, ZJMM_2, ZJMM_2_er, ZTT, ZTT_er;
    for (const auto &directory : directories)
    {

        Double_t Data = 0.0, Data_error = 0.0;
        Double_t llll = 0.0, llll_error = 0.0;
        Double_t lllljj = 0.0, lllljj_error = 0.0;
        Double_t llqq = 0.0, llqq_error = 0.0;
        Double_t llvv = 0.0, llvv_error = 0.0;
        Double_t llvvjj = 0.0, llvvjj_error = 0.0;
        Double_t llvvjjWW = 0.0, llvvjjWW_error = 0.0;
        Double_t top = 0.0, top_error = 0.0;
        Double_t ttbarVttbarVV = 0.0, ttbarVttbarVV_error = 0.0;
        Double_t VVV = 0.0, VVV_error = 0.0;
        Double_t Wjets = 0.0, Wjets_error = 0.0;
        Double_t Wt = 0.0, Wt_error = 0.0;
        Double_t WW = 0.0, WW_error = 0.0;
        Double_t WZ = 0.0, WZ_error = 0.0;
        Double_t WZjj = 0.0, WZjj_error = 0.0;
        Double_t Zjetsee = 0.0, Zjetsee_error = 0.0;
        Double_t Zjetsmumu = 0.0, Zjetsmumu_error = 0.0;
        Double_t Ztt = 0.0, Ztt_error = 0.0;

        Double_t Zjetsee_0 = 0.0, Zjetsee_error_0 = 0.0;
        Double_t Zjetsee_1 = 0.0, Zjetsee_error_1 = 0.0;
        Double_t Zjetsee_2 = 0.0, Zjetsee_error_2 = 0.0;

        Double_t Zjetsmumu_0 = 0.0, Zjetsmumu_error_0 = 0.0;
        Double_t Zjetsmumu_1 = 0.0, Zjetsmumu_error_1 = 0.0;
        Double_t Zjetsmumu_2 = 0.0, Zjetsmumu_error_2 = 0.0;

        cout << "=============" << directory << "=============" << endl;

        for (const auto &fileName : fileNames)
        {

            std::string filepath = "./SAMPLES" + directory + fileName;

            if (fileName == "/DATA.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_data;
                events_data = DoReco(file, 0, directory, fileName, HDATA, 0);
                Data += events_data.at(0);
                Data_error += events_data.at(1);
                file->Close();
            }
            if (fileName == "/llll.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_llll;
                events_llll = DoReco(file, 1, directory, fileName, Hllll, 0);
                llll += events_llll.at(0);
                llll_error += events_llll.at(1);
                if (directory == "/SR/")
                {
                    fourl = llll;
                    fourl_er = llll_error;
                }
                file->Close();
            }
            if (fileName == "/lllljj.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_lllljj;
                events_lllljj = DoReco(file, 1, directory, fileName, Hlllljj, 0);
                lllljj += events_lllljj.at(0);
                lllljj_error += events_lllljj.at(1);
                if (directory == "/SR/")
                {
                    fourljj = lllljj;
                    fourljj_er = lllljj_error;
                }
                file->Close();
            }
            if (fileName == "/llqq.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_llqq;
                events_llqq = DoReco(file, 1, directory, fileName, Hllqq, 0);
                llqq += events_llqq.at(0);
                llqq_error += events_llqq.at(1);
                if (directory == "/SR/")
                {
                    doublelqq = llqq;
                    doublelqq_er = llqq_error;
                }
                file->Close();
            }
            if (fileName == "/llvv.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_llvv;
                events_llvv = DoReco(file, 1, directory, fileName, Hllvv, 0);
                llvv += events_llvv.at(0);
                llvv_error += events_llvv.at(1);
                if (directory == "/SR/")
                {
                    doublelvv = llvv;
                    doublelvv_er = llvv_error;
                }
                file->Close();
            }
            if (fileName == "/llvvjj.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_llvvjj;
                events_llvvjj = DoReco(file, 1, directory, fileName, Hllvvjj, 0);
                llvvjj += events_llvvjj.at(0);
                llvvjj_error += events_llvvjj.at(1);
                if (directory == "/SR/")
                {
                    doublelvvjj = llvvjj;
                    doublelvvjj_er = llvvjj_error;
                }
                file->Close();
            }
            if (fileName == "/llvvjj_WW.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_llvvjj_WW;
                events_llvvjj_WW = DoReco(file, 1, directory, fileName, Hllvvjj_WW, 0);
                llvvjjWW += events_llvvjj_WW.at(0);
                llvvjjWW_error += events_llvvjj_WW.at(1);
                if (directory == "/SR/")
                {
                    doublelvvjjWW = llvvjjWW;
                    doublelvvjjWW_er = llvvjjWW_error;
                }
                file->Close();
            }
            if (fileName == "/top.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_top;
                events_top = DoReco(file, 1, directory, fileName, HTop, 0);
                top += events_top.at(0);
                top_error += events_top.at(1);
                if (directory == "/SR/")
                {
                    TOP = top;
                    TOP_er = top_error;
                }
                file->Close();
            }
            if (fileName == "/ttbarV_ttbarVV.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_ttbarV_ttbarVV;
                events_ttbarV_ttbarVV = DoReco(file, 1, directory, fileName, HttbarV_ttbarVV, 0);
                ttbarVttbarVV += events_ttbarV_ttbarVV.at(0);
                ttbarVttbarVV_error += events_ttbarV_ttbarVV.at(1);
                if (directory == "/SR/")
                {
                    ttbar = ttbarVttbarVV;
                    ttbar_er = ttbarVttbarVV_error;
                }
                file->Close();
            }
            if (fileName == "/VVV.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_VVV;
                events_VVV = DoReco(file, 1, directory, fileName, HVVV, 0);
                VVV += events_VVV.at(0);
                VVV_error += events_VVV.at(1);
                if (directory == "/SR/")
                {
                    tripleV = VVV;
                    tripleV_er = VVV_error;
                }
                file->Close();
            }
            if (fileName == "/W_jets.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_W_jets;
                events_W_jets = DoReco(file, 1, directory, fileName, HW_jets, 0);
                Wjets += events_W_jets.at(0);
                Wjets_error += events_W_jets.at(1);
                if (directory == "/SR/")
                {
                    WJ = Wjets;
                    WJ_er = Wjets_error;
                }
                file->Close();
            }
            if (fileName == "/Wt.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_Wt;
                events_Wt = DoReco(file, 1, directory, fileName, HWt, 0);
                Wt += events_Wt.at(0);
                Wt_error += events_Wt.at(1);
                if (directory == "/SR/")
                {
                    WT = Wt;
                    WT_er = Wt_error;
                }
                file->Close();
            }
            if (fileName == "/WW.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_WW;
                events_WW = DoReco(file, 1, directory, fileName, HWW, 0);
                WW += events_WW.at(0);
                WW_error += events_WW.at(1);
                if (directory == "/SR/")
                {
                    ww = WW;
                    ww_er = WW_error;
                }
                file->Close();
            }
            if (fileName == "/WZ.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_WZ;
                events_WZ = DoReco(file, 1, directory, fileName, HWZ, 0);
                WZ += events_WZ.at(0);
                WZ_error += events_WZ.at(1);
                if (directory == "/SR/")
                {
                    wz = WZ;
                    wz_er = WZ_error;
                }
                file->Close();
            }
            if (fileName == "/WZ_jj.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_WZ_jj;
                events_WZ_jj = DoReco(file, 1, directory, fileName, HWZ_jj, 0);
                WZjj += events_WZ_jj.at(0);
                WZjj_error += events_WZ_jj.at(1);
                if (directory == "/SR/")
                {
                    WZJJ = WZjj;
                    WZJJ_er = WZjj_error;
                }
                file->Close();
            }
            if (fileName == "/Z_jets_ee.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_Z_jets_ee;
                events_Z_jets_ee = DoReco(file, 1, directory, fileName, HZ_jets_ee, 0);
                Zjetsee += events_Z_jets_ee.at(0);
                Zjetsee_error += events_Z_jets_ee.at(1);
                file->Close();

                if (directory == "/SR/" || directory == "/emCR_A/")
                {
                    TFile *file = new TFile(filepath.c_str());
                    std::vector<float> events_Z_jets_ee_0;
                    events_Z_jets_ee_0 = DoReco(file, 1, directory, fileName, HZ_jets_ee_0, 1);
                    Zjetsee_0 += events_Z_jets_ee_0.at(2);
                    Zjetsee_error_0 += events_Z_jets_ee_0.at(3);
                    if (directory == "/SR/")
                    {
                        ZJEE = Zjetsee_0;
                        ZJEE_er = Zjetsee_error_0;
                    }

                    std::vector<float> events_Z_jets_ee_1;
                    events_Z_jets_ee_1 = DoReco(file, 1, directory, fileName, HZ_jets_ee_1, 2);
                    Zjetsee_1 += events_Z_jets_ee_1.at(4);
                    Zjetsee_error_1 += events_Z_jets_ee_1.at(5);
                    if (directory == "/SR/")
                    {
                        ZJEE_1 = Zjetsee_1;
                        ZJEE_1_er = Zjetsee_error_1;
                    }

                    std::vector<float> events_Z_jets_ee_2;
                    events_Z_jets_ee_2 = DoReco(file, 1, directory, fileName, HZ_jets_ee_2, 3);
                    Zjetsee_2 += events_Z_jets_ee_2.at(6);
                    Zjetsee_error_2 += events_Z_jets_ee_2.at(7);
                    if (directory == "/SR/")
                    {
                        ZJEE_2 = Zjetsee_2;
                        ZJEE_2_er = Zjetsee_error_2;
                    }

                    file->Close();
                }
            }
            if (fileName == "/Z_jets_mumu.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_Z_jets_mumu;
                events_Z_jets_mumu = DoReco(file, 1, directory, fileName, HZ_jets_mumu, 0);
                Zjetsmumu += events_Z_jets_mumu.at(0);
                Zjetsmumu_error += events_Z_jets_mumu.at(1);
                file->Close();

                if (directory == "/SR/" || directory == "/emCR_A/")
                {
                    TFile *file = new TFile(filepath.c_str());
                    std::vector<float> events_Z_jets_mumu_0;
                    events_Z_jets_mumu_0 = DoReco(file, 1, directory, fileName, HZ_jets_mumu_0, 1);
                    Zjetsmumu_0 += events_Z_jets_mumu_0.at(8);
                    Zjetsmumu_error_0 += events_Z_jets_mumu_0.at(9);
                    if (directory == "/SR/")
                    {
                        ZJMM = Zjetsmumu_0;
                        ZJMM_er = Zjetsmumu_error_0;
                    }

                    std::vector<float> events_Z_jets_mumu_1;
                    events_Z_jets_mumu_1 = DoReco(file, 1, directory, fileName, HZ_jets_mumu_1, 2);
                    Zjetsmumu_1 += events_Z_jets_mumu_1.at(10);
                    Zjetsmumu_error_1 += events_Z_jets_mumu_1.at(11);
                    if (directory == "/SR/")
                    {
                        ZJMM_1 = Zjetsmumu_1;
                        ZJMM_1_er = Zjetsmumu_error_1;
                    }

                    std::vector<float> events_Z_jets_mumu_2;
                    events_Z_jets_mumu_2 = DoReco(file, 1, directory, fileName, HZ_jets_mumu_2, 3);
                    Zjetsmumu_2 += events_Z_jets_mumu_2.at(12);
                    Zjetsmumu_error_2 += events_Z_jets_mumu_2.at(13);
                    if (directory == "/SR/")
                    {
                        ZJMM_2 = Zjetsmumu_2;
                        ZJMM_2_er = Zjetsmumu_error_2;
                    }

                    file->Close();
                }
            }
            if (fileName == "/Ztt.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_Ztt;
                events_Ztt = DoReco(file, 1, directory, fileName, HZtt, 0);
                Ztt += events_Ztt.at(0);
                Ztt_error += events_Ztt.at(1);
                if (directory == "/SR/")
                {
                    ZTT = Ztt;
                    ZTT_er = Ztt_error;
                }
                file->Close();
            }
        }

        //////// these are for all of 'em (3lCR,emCR_B,emCR_A,Zjets)//////////////
        Double_t SignalZZ = llvvjj + llvv;
        Double_t SignalZZ_error = sqrt(pow(llvvjj_error, 2) + pow(llvv_error, 2));
        Double_t EWKZZ = llvvjj;
        Double_t QCDZZ = llvv;
        Double_t Zjets = Zjetsmumu + Zjetsee;
        Double_t Zjets_error = sqrt(pow(Zjetsee_error, 2) + pow(Zjetsmumu_error, 2));
        Double_t AllTop = top + ttbarVttbarVV + Wt;
        Double_t AllTop_error = sqrt(pow(top_error, 2) + pow(ttbarVttbarVV_error, 2) + pow(Wt_error, 2));
        Double_t llklp = llll + llqq + VVV + Wjets + Ztt + lllljj + llvvjjWW + WZjj;
        Double_t llklp_error = sqrt(pow(llll_error, 2) + pow(llqq_error, 2) + pow(VVV_error, 2) + pow(Wjets_error, 2) + pow(Ztt_error, 2) + pow(llvvjj_error, 2) + pow(llvvjjWW_error, 2) + pow(WZjj_error, 2));
        Double_t BKG = WZ + AllTop + WW + llklp + Zjets;
        Double_t BKG_error = sqrt(pow(WZ_error, 2) + pow(AllTop_error, 2) + pow(WW_error, 2) + pow(llklp_error, 2) + pow(Zjets_error, 2));
        Double_t dataMCerror = (Data / BKG) * sqrt(pow(Data_error / Data, 2) + pow(BKG_error / BKG, 2));

        //////// these are for emCR_B ////////////
        Double_t WZemCR_B = WZ;
        Double_t WZemCR_B_error = WZ_error;
        Double_t BKG_emCR_B = WZemCR_B + AllTop + WW + llklp + Zjets;
        Double_t BKG_emCR_B_error = sqrt(pow(WZemCR_B_error, 2) + pow(AllTop_error, 2) + pow(WW_error, 2) + pow(llklp_error, 2) + pow(Zjets_error, 2));
        Double_t dataMCerror_emCR_B = (Data / BKG_emCR_B) * sqrt(pow(Data_error / Data, 2) + pow(BKG_emCR_B_error / BKG_emCR_B, 2));
        Double_t NR = AllTop;
        Double_t NR_error = sqrt(pow(AllTop_error, 2));
        Double_t nonNR = SignalZZ + llklp + WW + WZemCR_B + Zjets;
        Double_t nonNR_error = sqrt(pow(SignalZZ_error, 2) + pow(llklp_error, 2) + pow(WW_error, 2) + pow(WZemCR_B_error, 2) + pow(Zjets_error, 2));
        Double_t Numerator1 = Data - nonNR;
        Double_t Numerator1_error = sqrt(pow(Data_error, 2) + pow(nonNR_error, 2));
        Double_t DatanonNRoverNR = (Numerator1 / NR) * sqrt(pow(Numerator1_error / Numerator1, 2) + pow(NR_error / NR, 2));
        Double_t SignalZZNRerror = (SignalZZ / NR) * sqrt(pow(SignalZZ_error / SignalZZ, 2) + pow(NR_error / NR, 2));

        if (directory == "/emCR_B/")
        {
            printf("|                              |all                    |\n"
                   "|DATA                          |%.2f +- %.2f       |\n"
                   "|EWKZZ                         |%.2f +- %.2f           |\n"
                   "|QCDZZ                         |%.2f +- %.2f           |\n"
                   "|WZ                            |%.2f +- %.2f           |\n"
                   "|Zjets                         |%.2f +- %.2f           |\n"
                   "|top-ttbarV-ttbarVV-Wt         |%.2f +- %.2f       |\n"
                   "|WW                            |%.2f +- %.2f          |\n"
                   "|llll-llqq-VVV-Wjets-Ztt       |%.2f +- %.2f           |\n"
                   "|Signal(ZZ)                    |%.2f +- %.2f           |\n"
                   "|Background                    |%.2f +- %.2f       |\n"
                   "|Data/MC                       |%.2f +- %.2f           |\n"
                   "|Data - nonNR/NR               |%.2f +- %.2f           |\n"
                   "|Signal/NR(percent)            |%.2f +- %.2f           |\n",
                   Data, Data_error, EWKZZ, llvvjj_error, QCDZZ, llvv_error, WZemCR_B, WZemCR_B_error, Zjets, Zjets_error, AllTop, AllTop_error, WW, WW_error, llklp, llklp_error, SignalZZ, SignalZZ_error, BKG_emCR_B, BKG_emCR_B_error, (Data / BKG_emCR_B), dataMCerror_emCR_B, ((Data - nonNR) / NR), DatanonNRoverNR, 100 * (SignalZZ / NR), SignalZZNRerror);
            //////// Scale Factor for AllTop ///////////////////
            scaleAllTop = Numerator1 / NR;
            scaleAllTop_error = DatanonNRoverNR;
            cout << " The purity of emCR_B is (%): " << 100 * AllTop / (AllTop + SignalZZ + WZ + WW + Zjets + llklp) << endl;
        }

        //////// these are for 3lCR with scaled AllTop ///////////
        Double_t AllTop_3lCR = scaleAllTop * AllTop;
        Double_t AllTop_3lCR_error = abs(AllTop_3lCR) * sqrt(pow(scaleAllTop_error / scaleAllTop, 2) + pow(AllTop_error / AllTop, 2));
        Double_t nonWZ = SignalZZ + llklp + WW + AllTop_3lCR + Zjets;
        Double_t nonWZ_error = sqrt(pow(SignalZZ_error, 2) + pow(llklp_error, 2) + pow(WW_error, 2) + pow(AllTop_3lCR_error, 2) + pow(Zjets_error, 2));
        Double_t SignalZZWZerror = (SignalZZ / WZ) * sqrt(pow(SignalZZ_error / SignalZZ, 2) + pow(WZ_error / WZ, 2));
        Double_t Numerator = Data - nonWZ;
        Double_t Numeratorerror = sqrt(pow(Data_error, 2) + pow(nonWZ_error, 2));
        Double_t DatanonWZoverWZ = (Numerator / WZ) * sqrt(pow(Numeratorerror / Numerator, 2) + pow(WZ_error / WZ, 2));

        if (directory == "/3lCR/")
        {
            printf("|                              |all                    |\n"
                   "|DATA                          |%.2f +- %.2f       |\n"
                   "|EWKZZ                         |%.2f +- %.2f           |\n"
                   "|QCDZZ                         |%.2f +- %.2f           |\n"
                   "|WZ                            |%.2f +- %.2f       |\n"
                   "|Zjets                         |%.2f +- %.2f         |\n"
                   "|top-ttbarV-ttbarVV-Wt         |%.2f +- %.2f          |\n"
                   "|WW                            |%.2f +- %.2f           |\n"
                   "|llll-llqq-VVV-Wjets-Ztt       |%.2f +- %.2f          |\n"
                   "|Signal(ZZ)                    |%.2f +- %.2f           |\n"
                   "|Background                    |%.2f +- %.2f       |\n"
                   "|Data/MC                       |%.2f +- %.2f           |\n"
                   "|Data - nonWZ/WZ               |%.2f +- %.2f           |\n"
                   "|Signal/WZ(percent)            |%.2f +- %.2f           |\n",
                   Data, Data_error, EWKZZ, llvvjj_error, QCDZZ, llvv_error, WZ, WZ_error, Zjets, Zjets_error, AllTop_3lCR, AllTop_3lCR_error, WW, WW_error, llklp, llklp_error, SignalZZ, SignalZZ_error, BKG, BKG_error, (Data / BKG), dataMCerror, ((Data - nonWZ) / WZ), DatanonWZoverWZ, 100 * (SignalZZ / WZ), SignalZZWZerror);

            //////// Scale Factor for WZ ///////////////////
            scaleWZ = Numerator / WZ;
            scaleWZ_error = DatanonWZoverWZ;
            cout << " The purity of 3lCR is (%): " << 100 * WZ / (SignalZZ + WZ + AllTop + WW + Zjets + llklp) << endl;
        }

        //////// these are for Zjets with scaled WZ, scaled AllTop ////////////
        Double_t WZ_Zjets = scaleWZ * WZ;
        Double_t WZ_Zjets_error = abs(WZ_Zjets) * sqrt(pow(scaleWZ_error / scaleWZ, 2) + pow(WZ_error / WZ, 2));
        Double_t AllTop_Zjets = scaleAllTop * AllTop;
        Double_t AllTop_Zjets_error = abs(AllTop_Zjets) * sqrt(pow(scaleAllTop_error / scaleAllTop, 2) + pow(AllTop_error / AllTop, 2));
        Double_t WW_Zjets = WW;
        Double_t WW_Zjets_error = WW_error;
        Double_t BKG_Zjets = WZ_Zjets + AllTop_Zjets + WW_Zjets + llklp + Zjets;
        Double_t BKG_Zjets_error = sqrt(pow(WZ_Zjets_error, 2) + pow(AllTop_Zjets_error, 2) + pow(WW_Zjets_error, 2) + pow(llklp_error, 2) + pow(Zjets_error, 2));
        Double_t dataMCerror_Zjets = (Data / BKG_Zjets) * sqrt(pow(Data_error / Data, 2) + pow(BKG_Zjets_error / BKG_Zjets, 2));
        Double_t nonZjets = SignalZZ + llklp + WW_Zjets + WZ_Zjets + AllTop_Zjets;
        Double_t nonZjets_error = sqrt(pow(SignalZZ_error, 2) + pow(llklp_error, 2) + pow(WW_Zjets_error, 2) + pow(WZ_Zjets_error, 2) + pow(AllTop_Zjets_error, 2));
        Double_t Numerator2 = Data - nonZjets;
        Double_t Numerator2_error = sqrt(pow(Data_error, 2) + pow(nonZjets_error, 2));
        Double_t DatanonZjetsoverZjets = (Numerator2 / Zjets) * sqrt(pow(Numerator2_error / Numerator2, 2) + pow(Zjets_error / Zjets, 2));
        Double_t SignalZZZjetserror = (SignalZZ / Zjets) * sqrt(pow(SignalZZ_error / SignalZZ, 2) + pow(Zjets_error / Zjets, 2));

        //////// these are for Zjets with n_jets = 2 with scaled WZ, scaled AllTop ////////////
        if (directory == "/Zjets2/")
        {
            printf("|                              |all                     |\n"
                   "|DATA                          |%.2f +- %.2f      |\n"
                   "|EWKZZ                         |%.2f +- %.2f           |\n"
                   "|QCDZZ                         |%.2f +- %.2f        |\n"
                   "|WZ                            |%.2f +- %.2f        |\n"
                   "|Zjets                         |%.2f +- %.2f      |\n"
                   "|top-ttbarV-ttbarVV-Wt         |%.2f +- %.2f        |\n"
                   "|WW                            |%.2f +- %.2f         |\n"
                   "|llll-llqq-VVV-Wjets-Ztt       |%.2f +- %.2f          |\n"
                   "|Signal(ZZ)                    |%.2f +- %.2f        |\n"
                   "|Background                    |%.2f +- %.2f      |\n"
                   "|Data/MC                       |%.2f +- %.2f            |\n"
                   "|Data - nonZjets/Zjets         |%.2f +- %.2f            |\n"
                   "|Signal/Zjets(percent)         |%.2f +- %.2f           |\n",
                   Data, Data_error, EWKZZ, llvvjj_error, QCDZZ, llvv_error, WZ_Zjets, WZ_Zjets_error, Zjets, Zjets_error, AllTop_Zjets, AllTop_Zjets_error, WW_Zjets, WW_Zjets_error, llklp, llklp_error, SignalZZ, SignalZZ_error, BKG_Zjets, BKG_Zjets_error, (Data / BKG_Zjets), dataMCerror_Zjets, ((Data - nonZjets) / Zjets), DatanonZjetsoverZjets, 100 * (SignalZZ / Zjets), SignalZZZjetserror);
            //////// Scale Factor for Zjets2 ///////////////////
            scaleZjets2 = Numerator2 / Zjets;
            scaleZjets2_error = DatanonZjetsoverZjets;
            cout << " The purity of Zjets_2 CR is (%): " << 100 * Zjets / (Zjets + SignalZZ + WZ + WW + AllTop + llklp) << endl;
        }

        //////// these are for Zjets with n_jets = 1 with scaled WZ, scaled AllTop ////////////
        if (directory == "/Zjets1/")
        {
            printf("|                              |all                     |\n"
                   "|DATA                          |%.2f +- %.2f      |\n"
                   "|EWKZZ                         |%.2f +- %.2f           |\n"
                   "|QCDZZ                         |%.2f +- %.2f        |\n"
                   "|WZ                            |%.2f +- %.2f        |\n"
                   "|Zjets                         |%.2f +- %.2f      |\n"
                   "|top-ttbarV-ttbarVV-Wt         |%.2f +- %.2f        |\n"
                   "|WW                            |%.2f +- %.2f         |\n"
                   "|llll-llqq-VVV-Wjets-Ztt       |%.2f +- %.2f          |\n"
                   "|Signal(ZZ)                    |%.2f +- %.2f        |\n"
                   "|Background                    |%.2f +- %.2f      |\n"
                   "|Data/MC                       |%.2f +- %.2f            |\n"
                   "|Data - nonZjets/Zjets         |%.2f +- %.2f            |\n"
                   "|Signal/Zjets(percent)         |%.2f +- %.2f           |\n",
                   Data, Data_error, EWKZZ, llvvjj_error, QCDZZ, llvv_error, WZ_Zjets, WZ_Zjets_error, Zjets, Zjets_error, AllTop_Zjets, AllTop_Zjets_error, WW_Zjets, WW_Zjets_error, llklp, llklp_error, SignalZZ, SignalZZ_error, BKG_Zjets, BKG_Zjets_error, (Data / BKG_Zjets), dataMCerror_Zjets, ((Data - nonZjets) / Zjets), DatanonZjetsoverZjets, 100 * (SignalZZ / Zjets), SignalZZZjetserror);
            //////// Scale Factor for Zjets1 ///////////////////
            scaleZjets1 = Numerator2 / Zjets;
            scaleZjets1_error = DatanonZjetsoverZjets;
            cout << " The purity of Zjets_1 CR is: " << 100 * Zjets / (Zjets + SignalZZ + WZ + WW + AllTop + llklp) << endl;
        }

        if (directory == "/Zjets0/")
        {
            printf("|                              |all                     |\n"
                   "|DATA                          |%.2f +- %.2f      |\n"
                   "|EWKZZ                         |%.2f +- %.2f           |\n"
                   "|QCDZZ                         |%.2f +- %.2f        |\n"
                   "|WZ                            |%.2f +- %.2f        |\n"
                   "|Zjets                         |%.2f +- %.2f      |\n"
                   "|top-ttbarV-ttbarVV-Wt         |%.2f +- %.2f        |\n"
                   "|WW                            |%.2f +- %.2f         |\n"
                   "|llll-llqq-VVV-Wjets-Ztt       |%.2f +- %.2f          |\n"
                   "|Signal(ZZ)                    |%.2f +- %.2f        |\n"
                   "|Background                    |%.2f +- %.2f      |\n"
                   "|Data/MC                       |%.2f +- %.2f            |\n"
                   "|Data - nonZjets/Zjets         |%.2f +- %.2f            |\n"
                   "|Signal/Zjets(percent)         |%.2f +- %.2f           |\n",
                   Data, Data_error, EWKZZ, llvvjj_error, QCDZZ, llvv_error, WZ_Zjets, WZ_Zjets_error, Zjets, Zjets_error, AllTop_Zjets, AllTop_Zjets_error, WW_Zjets, WW_Zjets_error, llklp, llklp_error, SignalZZ, SignalZZ_error, BKG_Zjets, BKG_Zjets_error, (Data / BKG_Zjets), dataMCerror_Zjets, ((Data - nonZjets) / Zjets), DatanonZjetsoverZjets, 100 * (SignalZZ / Zjets), SignalZZZjetserror);
            //////// Scale Factor for Zjets0 ///////////////////
            scaleZjets0 = Numerator2 / Zjets;
            scaleZjets0_error = DatanonZjetsoverZjets;
            cout << " The purity of Zjets_0 CR is: " << 100 * Zjets / (SignalZZ + Zjets + WZ + WW + AllTop + llklp) << endl;
        }

        //////// these are for emCR_A with scaled WZ, scaled AllTop and scaled Zjets0,1,2 ////////////
        Double_t WZemCR_A = scaleWZ * WZ;
        Double_t WZemCR_A_error = abs(WZemCR_A) * sqrt(pow(scaleWZ_error / scaleWZ, 2) + pow(WZ_error / WZ, 2));
        Double_t AllTop_A = scaleAllTop * AllTop;
        Double_t AllTop_A_error = abs(AllTop_A) * sqrt(pow(scaleAllTop_error / scaleAllTop, 2) + pow(AllTop_error / AllTop, 2));
        Double_t ZjetsemCR_A = scaleZjets0 * (Zjetsee_0 + Zjetsmumu_0) + scaleZjets1 * (Zjetsee_1 + Zjetsmumu_1) + scaleZjets2 * (Zjetsee_2 + Zjetsmumu_2);
        Double_t ZjetsemCR_A_error = sqrt(
            (scaleZjets0 * Zjetsee_error_0) * (scaleZjets0 * Zjetsee_error_0) +
            (scaleZjets1 * Zjetsee_error_1) * (scaleZjets1 * Zjetsee_error_1) +
            (scaleZjets2 * Zjetsee_error_2) * (scaleZjets2 * Zjetsee_error_2) +
            (scaleZjets0 * Zjetsmumu_error_0) * (scaleZjets0 * Zjetsmumu_error_0) +
            (scaleZjets1 * Zjetsmumu_error_1) * (scaleZjets1 * Zjetsmumu_error_1) +
            (scaleZjets2 * Zjetsmumu_error_2) * (scaleZjets2 * Zjetsmumu_error_2));

        Double_t BKG_emCR_A = WZemCR_A + AllTop_A + WW + llklp + ZjetsemCR_A;
        Double_t BKG_emCR_A_error = sqrt(pow(WZemCR_A_error, 2) + pow(AllTop_A_error, 2) + pow(WW_error, 2) + pow(llklp_error, 2) + pow(ZjetsemCR_A_error, 2));
        Double_t dataMCerror_emCR_A = (Data / BKG_emCR_A) * sqrt(pow(Data_error / Data, 2) + pow(BKG_emCR_A_error / BKG_emCR_A, 2));
        Double_t NR_A = WW;
        Double_t NR_A_error = sqrt(pow(WW_error, 2));
        Double_t nonNR_A = SignalZZ + llklp + WZemCR_A + ZjetsemCR_A + AllTop_A;
        Double_t nonNR_A_error = sqrt(pow(SignalZZ_error, 2) + pow(llklp_error, 2) + pow(WZemCR_A_error, 2) + pow(ZjetsemCR_A_error, 2) + pow(AllTop_A_error, 2));
        Double_t Numerator3 = Data - nonNR_A;
        Double_t Numerator3_error = sqrt(pow(Data_error, 2) + pow(nonNR_A_error, 2));
        Double_t DatanonNRoverNR_A = (Numerator3 / NR_A) * sqrt(pow(Numerator3_error / Numerator3, 2) + pow(NR_A_error / NR_A, 2));
        Double_t SignalZZNRerror_A = (SignalZZ / NR_A) * sqrt(pow(SignalZZ_error / SignalZZ, 2) + pow(NR_A_error / NR_A, 2));
        if (directory == "/emCR_A/")
        {
            printf("|                              |all                    |\n"
                   "|DATA                          |%.2f +- %.2f       |\n"
                   "|EWKZZ                         |%.2f +- %.2f           |\n"
                   "|QCDZZ                         |%.2f +- %.2f           |\n"
                   "|WZ                            |%.2f +- %.2f          |\n"
                   "|Zjets                         |%.2f +- %.2f           |\n"
                   "|top-ttbarV-ttbarVV-Wt         |%.2f +- %.2f       |\n"
                   "|WW                            |%.2f +- %.2f         |\n"
                   "|llll-llqq-VVV-Wjets-Ztt       |%.2f +- %.2f         |\n"
                   "|Signal(ZZ)                    |%.2f +- %.2f           |\n"
                   "|Background                    |%.2f +- %.2f       |\n"
                   "|Data/MC                       |%.2f +- %.2f           |\n"
                   "|Data - nonNR/NR               |%.2f +- %.2f           |\n"
                   "|Signal/NR(percent)            |%.2f +- %.2f           |\n",
                   Data, Data_error, EWKZZ, llvvjj_error, QCDZZ, llvv_error, WZemCR_A, WZemCR_A_error, ZjetsemCR_A, ZjetsemCR_A_error, AllTop_A, AllTop_A_error, WW, WW_error, llklp, llklp_error, SignalZZ, SignalZZ_error, BKG_emCR_A, BKG_emCR_A_error, (Data / BKG_emCR_A), dataMCerror_emCR_A, ((Data - nonNR_A) / NR_A), DatanonNRoverNR_A, 100 * (SignalZZ / NR_A), SignalZZNRerror_A);
            //////// Scale Factor for WW ///////////////////
            scaleWW = Numerator3 / NR_A;
            scaleWW_error = DatanonNRoverNR_A;
            cout << " The purity of emCR_A is (%): " << 100 * WW / (WW + WZ + SignalZZ + AllTop + Zjets + llklp) << endl;
        }

        if (directory == "/SR/")
        {
            Double_t AllTop_SR = scaleAllTop * AllTop;
            Double_t AllTop_SR_error = abs(AllTop_SR) * sqrt(pow(scaleAllTop_error / scaleAllTop, 2) + pow(AllTop_error / AllTop, 2));
            Double_t WZ_SR = scaleWZ * WZ;
            Double_t WZ_SR_error = abs(WZ_SR) * sqrt(pow(scaleWZ_error / scaleWZ, 2) + pow(WZ_error / WZ, 2));
            Double_t WW_SR = scaleWW * WW;
            Double_t WW_SR_error = abs(WW_SR) * sqrt(pow(scaleWW_error / scaleWW, 2) + pow(WW_error / WW, 2));

            Double_t ZjetsSR = scaleZjets0 * (Zjetsee_0 + Zjetsmumu_0) + scaleZjets1 * (Zjetsee_1 + Zjetsmumu_1) + scaleZjets2 * (Zjetsee_2 + Zjetsmumu_2);
            Double_t ZjetsSR_unscaled = (Zjetsee_0 + Zjetsmumu_0) + (Zjetsee_1 + Zjetsmumu_1) + (Zjetsee_2 + Zjetsmumu_2);
            Double_t ZjetsSR_un_error = sqrt(pow(Zjetsee_0, 2) + pow(Zjetsmumu_0, 2) + pow(Zjetsee_1, 2) + pow(Zjetsmumu_1, 2) + pow(Zjetsee_2, 2) + pow(Zjetsmumu_2, 2));
            Double_t ZjetsSR_error = sqrt(
                (scaleZjets0 * Zjetsee_error_0) * (scaleZjets0 * Zjetsee_error_0) +
                (scaleZjets1 * Zjetsee_error_1) * (scaleZjets1 * Zjetsee_error_1) +
                (scaleZjets2 * Zjetsee_error_2) * (scaleZjets2 * Zjetsee_error_2) +
                (scaleZjets0 * Zjetsmumu_error_0) * (scaleZjets0 * Zjetsmumu_error_0) +
                (scaleZjets1 * Zjetsmumu_error_1) * (scaleZjets1 * Zjetsmumu_error_1) +
                (scaleZjets2 * Zjetsmumu_error_2) * (scaleZjets2 * Zjetsmumu_error_2));

            Double_t non = WZ_SR + AllTop_SR + WW_SR + ZjetsSR + llklp;
            Double_t non_err = pow(WZ_SR_error, 2) + pow(AllTop_SR_error, 2) + pow(WW_SR_error, 2) + pow(ZjetsSR_error, 2) + pow(llklp_error, 2);
            Double_t panw = Data - non;
            Double_t panw_error = sqrt(pow(Data_error, 2) + pow(non_err, 2));
            Double_t DatanonoverNr = (panw / SignalZZ) * sqrt(pow((panw_error / panw), 2) + pow((SignalZZ_error / SignalZZ), 2));

            Double_t Z0 = Zjetsee_0 + Zjetsmumu_0;
            Double_t Z0_error = sqrt(pow(Zjetsee_error_0, 2) + pow(Zjetsmumu_error_0, 2));
            Double_t scaledZ0 = scaleZjets0 * Z0;
            Double_t scaledZ0_error = abs(scaledZ0) * sqrt(pow(scaleZjets0_error / scaleZjets0, 2) + pow(Z0_error / Z0, 2));
            Double_t Z1 = Zjetsee_1 + Zjetsmumu_1;
            Double_t Z1_error = sqrt(pow(Zjetsee_error_1, 2) + pow(Zjetsmumu_error_1, 2));
            Double_t scaledZ1 = scaleZjets1 * Z1;
            Double_t scaledZ1_error = abs(scaledZ1) * sqrt(pow(scaleZjets1_error / scaleZjets1, 2) + pow(Z1_error / Z1, 2));
            Double_t Z2 = Zjetsee_2 + Zjetsmumu_2;
            Double_t Z2_error = sqrt(pow(Zjetsee_error_2, 2) + pow(Zjetsmumu_error_2, 2));
            Double_t scaledZ2 = scaleZjets2 * Z2;
            Double_t scaledZ2_error = abs(scaledZ2) * sqrt(pow(scaleZjets2_error / scaleZjets2, 2) + pow(Z2_error / Z2, 2));

            Double_t BKG_SR = WZ_SR + AllTop_SR + WW_SR + llklp + ZjetsSR;
            Double_t BKG_SR_error = sqrt(pow(WZ_SR_error, 2) + pow(AllTop_SR_error, 2) + pow(WW_SR_error, 2) + pow(llklp_error, 2) + pow(ZjetsSR_error, 2));
            Double_t MC_all = BKG_SR + SignalZZ;
            Double_t MC_all_error = sqrt(pow(BKG_SR_error, 2) + pow(SignalZZ_error, 2));
            Double_t MC = WZ + AllTop + WW + Z0 + Z1 + Z2 + llklp + SignalZZ;
            Double_t MC_error = sqrt(pow(WZ_error, 2) + pow(SignalZZ_error, 2) + pow(AllTop_error, 2) + pow(WW_error, 2) + pow(llklp_error, 2) + pow(Z0_error, 2) + pow(Z1_error, 2) + pow(Z2_error, 2));
            Double_t dataMCerror_SR = (MC_all / MC_all) * sqrt(pow(MC_all_error / MC_all, 2) + pow(MC_all_error / MC_all, 2));
            Double_t SignalZZNRerror = (SignalZZ / SignalZZ) * sqrt(pow(SignalZZ_error / SignalZZ, 2) + pow(SignalZZ_error / SignalZZ, 2));
            Double_t DATAMC = (MC / MC) * sqrt(pow(MC_error / MC, 2) + pow(MC_error / MC, 2));
            Double_t errorMCallminussignal = sqrt(pow(MC_all_error, 2) + pow(SignalZZ_error, 2));
            Double_t errorpanw = sqrt(pow(MC_all_error, 2) + pow(errorMCallminussignal, 2));
            Double_t erroroliko = 1 * sqrt(pow(errorpanw / (SignalZZ), 2) + pow(SignalZZ_error / SignalZZ, 2));

            printf("|                              |all                    |\n"
                   "|Data                          |%.2f +- %.2f      |\n"
                   "|DATA(fake)-MC                 |%.2f +- %.2f      |\n"
                   "|EWKZZ                         |%.2f +- %.2f           |\n"
                   "|QCDZZ                         |%.2f +- %.2f       |\n"
                   "|WZ                            |%.2f +- %.2f        |\n"
                   "|Zjets                         |%.2f +- %.2f       |\n"
                   "|Zjets0                        |%.2f +- %.2f         |\n"
                   "|Zjets1                        |%.2f +- %.2f         |\n"
                   "|Zjets2                        |%.2f +- %.2f         |\n"
                   "|top-ttbarV-ttbarVV-Wt         |%.2f +- %.2f          |\n"
                   "|WW                            |%.2f +- %.2f          |\n"
                   "|llll-llqq-VVV-Wjets-Ztt       |%.2f +- %.2f          |\n"
                   "|Signal(ZZ)                    |%.2f +- %.2f       |\n"
                   "|Background                    |%.2f +- %.2f      |\n"
                   "|Data/MC                       |%.2f +- %.2f           |\n"
                   "|Data - nonNR/NR               |%.2f +- %.2f           |\n"
                   "|Signal/NR(percent)            |%.2f +- %.2f        |\n",
                   Data, Data_error, MC_all, sqrt(MC_all), EWKZZ, llvvjj_error, QCDZZ, llvv_error, WZ_SR, WZ_SR_error, ZjetsSR, ZjetsSR_error, scaledZ0, scaledZ0_error, scaledZ1, scaledZ1_error, scaledZ2, scaledZ2_error, AllTop_SR, AllTop_SR_error, WW_SR, WW_SR_error, llklp, llklp_error, SignalZZ, SignalZZ_error, BKG_SR, BKG_SR_error, (MC_all / MC_all), dataMCerror_SR, ((MC_all - (MC_all - SignalZZ)) / SignalZZ), erroroliko, 100 * (SignalZZ / SignalZZ), SignalZZNRerror);

            printf("|                              |all                    |\n"
                   "|Data                          |%.2f +- %.2f      |\n"
                   "|DATA (fake)-MC                |%.2f +- %.2f      |\n"
                   "|EWKZZ                         |%.2f +- %.2f           |\n"
                   "|QCDZZ                         |%.2f +- %.2f       |\n"
                   "|WZ                            |%.2f +- %.2f        |\n"
                   "|Zjets                         |%.2f +- %.2f       |\n"
                   "|Zjets0                        |%.2f +- %.2f         |\n"
                   "|Zjets1                        |%.2f +- %.2f         |\n"
                   "|Zjets2                        |%.2f +- %.2f         |\n"
                   "|top-ttbarV-ttbarVV-Wt         |%.2f +- %.2f          |\n"
                   "|WW                            |%.2f +- %.2f          |\n"
                   "|llll-llqq-VVV-Wjets-Ztt       |%.2f +- %.2f          |\n"
                   "|Signal(ZZ)                    |%.2f +- %.2f       |\n"
                   "|Background                    |%.2f +- %.2f      |\n"
                   "|Data/MC                       |%.2f +- %.2f           |\n"
                   "|Data - nonNR/NR               |%.2f +- %.2f           |\n"
                   "|Signal/NR(percent)            |%.2f +- %.2f        |\n",
                   Data, Data_error, MC, sqrt(MC), EWKZZ, llvvjj_error, QCDZZ, llvv_error, WZ, WZ_error, Zjets, Zjets_error, Z0, Z0_error, Z1, Z1_error, Z2, Z2_error, AllTop, AllTop_error, WW, WW_error, llklp, llklp_error, SignalZZ, SignalZZ_error, BKG, BKG_error, (MC / MC), DATAMC, ((Data - non) / SignalZZ), DatanonoverNr, 100 * (SignalZZ / NR), SignalZZNRerror);

            N_SR = Data - WZ_SR - AllTop_SR - WW_SR - ZjetsSR - llklp;
            Double_t signalunscaled_error = sqrt(pow(Data_error, 2) + pow(WZ_error, 2) + pow(AllTop_error, 2) + pow(WW_error, 2) + pow(ZjetsSR_un_error, 2));
            Double_t signalunscaled = Data - WZ - AllTop - WW - (Zjetsee_0 + Zjetsmumu_0 + Zjetsee_1 + Zjetsmumu_1 + Zjetsee_2 + Zjetsmumu_2) - llklp;
            N_SR_error = sqrt(pow(Data_error, 2) + pow(WZ_SR_error, 2) + pow(WW_SR_error, 2) + pow(AllTop_SR_error, 2) + pow(ZjetsSR_error, 2) + pow(llklp_error, 2));
            scale_Signal = N_SR / SignalZZ;
            scale_Signal_error = abs(scale_Signal) * sqrt(pow(N_SR_error / N_SR, 2) + pow(SignalZZ_error / SignalZZ, 2));
            sfsignal = (MC_all - BKG_SR) / SignalZZ;
            Double_t er = sqrt(pow(MC_all_error, 2) + pow(BKG_SR_error, 2));
            cout << " Mc summation is (scaled all) +- " << MC_all << " +- " << MC_all_error << endl;
            sfsignal_error = abs(sfsignal) * sqrt(pow(er / (MC_all - BKG_SR), 2) + pow(SignalZZ_error / SignalZZ, 2));
            cout << "eisai konta sto 0.39???  " << sfsignal_error << endl;   
            printf("The N_SR signal after the corrections-scales is: %f +- %f\n", N_SR, N_SR_error);
            printf("THe signal for SR region without corrections is: %f +- %f\n", signalunscaled, signalunscaled_error);
            printf("Signal/BKG (unscaled) : %f\n", (signalunscaled / (WZ + AllTop + WW + llklp + (Zjetsee_0 + Zjetsmumu_0 + Zjetsee_1 + Zjetsmumu_1 + Zjetsee_2 + Zjetsmumu_2))));
            printf("Signal/BKG (scaled) : %f\n", (N_SR / (WZ_SR + AllTop_SR + WW_SR + llklp + ZjetsSR)));
            printf("The number that Zjets events found in SR region is %f\n", ZjetsSR);
            printf("The number for Zjets0 is %f +- %f\n", (Zjetsee_0 + Zjetsmumu_0), sqrt(pow(Zjetsee_error_0, 2) + pow(Zjetsmumu_error_0, 2)));
            printf("The number for Zjets1 is %f +- %f\n", (Zjetsee_1 + Zjetsmumu_1), sqrt(pow(Zjetsee_error_1, 2) + pow(Zjetsmumu_error_1, 2)));
            printf("The number for Zjets2 is %f +- %f\n", (Zjetsee_2 + Zjetsmumu_2), sqrt(pow(Zjetsee_error_2, 2) + pow(Zjetsmumu_error_2, 2)));
            asa = erroroliko;

            ////////// ssytimatiko sflma kapws ///////////
            TCanvas *c200 = new TCanvas("c200", "Systematic", 0., 0., 700, 700);
            TH1F *h_mu_s = new TH1F("H. Systematic", "", 200, 0.8, 1.2);
            TRandom3 rnd;
            for (int i = 1; i <= 50000; i++)
            {
                Double_t WZ_SYST = rnd.Gaus(704.05, 18.00);
                Double_t TOP_SYST = rnd.Gaus(102.46, 3.00);
                Double_t WW_SYST = rnd.Gaus(56.25, 4.17);
                Double_t ZJETS2_SYST = rnd.Gaus(1.54, 0.81);
                Double_t ZJETS1_SYST = rnd.Gaus(56.63, 13.78);
                Double_t ZJETS0_SYST = rnd.Gaus(147.21, 34.19);
                Double_t OTHER_SYST = rnd.Gaus(45.97, 2.06);
                Double_t SIGNAL_SYST = rnd.Gaus(1607.30, 13.93);
                //Double_t PSEUDO_DATA = WZ_SYST + TOP_SYST + WW_SYST + ZJETS0_SYST + ZJETS1_SYST + ZJETS2_SYST + OTHER_SYST + SIGNAL_SYST;
                //Double_t PSEUDO_DATA_SYST = rnd.Gaus(PSEUDO_DATA, sqrt(PSEUDO_DATA));
                Double_t PSEUDO_DATA_SYST = 2721.40;
                Double_t scaling = (PSEUDO_DATA_SYST - WZ_SYST - TOP_SYST - WW_SYST - ZJETS2_SYST - ZJETS0_SYST - ZJETS1_SYST - OTHER_SYST) / (SIGNAL_SYST);
                c200->cd();
                h_mu_s->Fill(scaling);
                /*if ((i % 10000) == 0)
                {
                    cout << " WZ_SYST = " << (WZ_SYST-704.05)/18.00 << " sigma " << endl;
                    cout << " TOP_SYST = " << (TOP_SYST-102.46)/3 << " sigma " << endl;
                }*/
            }
            h_mu_s->GetXaxis()->SetTitle("#mu_{s}"); // Here needs adjustment
            h_mu_s->GetYaxis()->SetTitle("Entries");
            h_mu_s->SetFillColor(41); // Signal
            h_mu_s->Draw("hist");
            Double_t systematic = h_mu_s->GetRMS();

            cout << " The systematic ancertainty is: " << systematic << endl;
        }
    }

    printf("scaling factors are: sf_WZ = %f +- %f\n"
           "sf_top = %f +- %f\n"
           "sf_WW = %f +- %f\n"
           "sf_zjets0 = %f +- %f\n"
           "sf_zjets1 = %f +- %f\n"
           "sf_zjets2 = %f +- %f\n"
           "sf_signal = %f +- %f\n"
           "mu_s = %f +- %f\n",
           scaleWZ, scaleWZ_error, scaleAllTop, scaleAllTop_error, scaleWW, scaleWW_error, scaleZjets0, scaleZjets0_error, scaleZjets1, scaleZjets1_error, scaleZjets2, scaleZjets2_error, scale_Signal, scale_Signal_error, sfsignal, asa);

    // Some useful outputs
    double totalBKG = fourl + fourljj + doublelqq + doublelvvjjWW + TOP + ttbar + tripleV + WJ + WT + ww + wz + WZJJ + ZJEE + ZJEE_1 + ZJEE_2 + ZJMM + ZJMM_1 + ZJMM_2 + ZTT;
    double totalBKG_er = sqrt(pow(fourl_er, 2) + pow(fourljj_er, 2) + pow(doublelqq_er, 2) + pow(doublelvvjjWW_er, 2) + pow(TOP_er, 2) + pow(ttbar_er, 2) + pow(tripleV_er, 2) + pow(WJ_er, 2) + pow(WT_er, 2) + pow(ww_er, 2) + pow(wz_er, 2) + pow(WZJJ_er, 2) + pow(ZJEE_er, 2) + pow(ZJEE_1_er, 2) + pow(ZJEE_2_er, 2) + pow(ZJMM_er, 2) + pow(ZJMM_1_er, 2) + pow(ZJMM_2_er, 2) + pow(ZTT_er, 2));
    cout << "total Bkg=" << totalBKG << "+-" << totalBKG_er << endl;
    cout << "S/B=" << (doublelvv + doublelvvjj) / totalBKG << endl;
    Double_t S = doublelvv + doublelvvjj;
    Double_t B = totalBKG;
    printf("S and B (signal and background) are %f and %f\n", S, B);
    Double_t Z = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
    cout << " Significance = " << Z << endl;

    TCanvas *c1 = new TCanvas("c1", "Canvadi", 0., 0., 700, 700);
    TPad *pad_1 = new TPad("pad_1", "This is pad_1", 0.01, 0.30, 1., 1.);
    TPad *pad_2 = new TPad("pad_2", "This is pad_2", 0.01, 0.01, 1., 0.30);

    pad_1->SetBorderSize(0);
    pad_1->SetBottomMargin(0.02);
    pad_1->Draw();
    pad_2->SetBottomMargin(0.35);
    pad_2->SetTopMargin(0.0);
    pad_2->SetBorderSize(0);
    pad_2->Draw();
    pad_1->cd();
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Adding Zjetsee0 + Zjetsmumu0
    for (int i = 1; i <= HZ_jets_ee_0->GetNbinsX(); i++)
    {
        Double_t binsofZ0 = HZ_jets_ee_0->GetBinContent(i) + HZ_jets_mumu_0->GetBinContent(i);
        HZjets0->SetBinContent(i, binsofZ0);
    }
    TH1F *h02 = (TH1F *)HZjets0->Clone("h02");
    // Scaling HZjets0
    for (int i = 1; i <= HZjets0->GetNbinsX(); i++)
    {
        Double_t binscaledZ0 = HZjets0->GetBinContent(i);
        Double_t scaledTimesContentZjets0 = binscaledZ0 * scaleZjets0;
        HZjets0scaled->SetBinContent(i, scaledTimesContentZjets0);
    }
    // Adding Zjetsee1 + Zjetsmumu1
    for (int i = 1; i <= HZ_jets_ee_1->GetNbinsX(); i++)
    {
        Double_t binsofZ1 = HZ_jets_ee_1->GetBinContent(i) + HZ_jets_mumu_1->GetBinContent(i);
        HZjets1->SetBinContent(i, binsofZ1);
    }
    TH1F *h03 = (TH1F *)HZjets1->Clone("h03");
    // Scaling HZjets1
    for (int i = 1; i <= HZjets1->GetNbinsX(); i++)
    {
        Double_t binscaledZ1 = HZjets1->GetBinContent(i);
        Double_t scaledTimesContentZjets1 = binscaledZ1 * scaleZjets1;
        HZjets1scaled->SetBinContent(i, scaledTimesContentZjets1);
    }
    // Adding Zjetsee2 + Zjetsmumu2
    for (int i = 1; i <= HZ_jets_ee_2->GetNbinsX(); i++)
    {
        Double_t binsofZ2 = HZ_jets_ee_2->GetBinContent(i) + HZ_jets_mumu_2->GetBinContent(i);
        HZjets2->SetBinContent(i, binsofZ2);
    }
    TH1F *h04 = (TH1F *)HZjets2->Clone("h04");
    // Scaling HZjets2
    for (int i = 1; i <= HZjets2->GetNbinsX(); i++)
    {
        Double_t binscaledZ2 = HZjets2->GetBinContent(i);
        Double_t scaledTimesContentZjets2 = binscaledZ2 * scaleZjets2;
        HZjets2scaled->SetBinContent(i, scaledTimesContentZjets2);
    }
    // Correction for Zjets histogram with scaling factors
    for (int i = 1; i <= HZjetsscaled->GetNbinsX(); i++)
    {
        Double_t binsofZjets = HZjets0scaled->GetBinContent(i) + HZjets1scaled->GetBinContent(i) + HZjets2scaled->GetBinContent(i);
        HZjetsscaled->SetBinContent(i, binsofZjets);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Adding all top together
    for (int i = 1; i <= HTop->GetNbinsX(); i++)
    {
        Double_t binContentT = HTop->GetBinContent(i) + HttbarV_ttbarVV->GetBinContent(i) + HWt->GetBinContent(i);
        HTops->SetBinContent(i, binContentT);
    }
    // Correction for Tops histogram with scaling factor scaleAllTop
    for (int i = 1; i <= HTops->GetNbinsX(); i++)
    {
        Double_t binContentTops = HTops->GetBinContent(i);
        Double_t scaledTimesContentTops = binContentTops * scaleAllTop;
        HTopsscaled->SetBinContent(i, scaledTimesContentTops);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Correction for WZ histogram with scaling factor scaleWZ
    for (int i = 1; i <= HWZ->GetNbinsX(); i++)
    {
        Double_t binContentWZ = HWZ->GetBinContent(i);
        Double_t scaledTimesContentWZ = binContentWZ * scaleWZ;
        HWZscaled->SetBinContent(i, scaledTimesContentWZ);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Correction for WW histogram with scaling factor scaleWZ
    for (int i = 1; i <= HWW->GetNbinsX(); i++)
    {
        Double_t binContentWW = HWW->GetBinContent(i);
        Double_t scaledTimesContentWW = binContentWW * scaleWW;
        HWWscaled->SetBinContent(i, scaledTimesContentWW);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // adding all histogramms to create the backgroung Histogram
    /*for (int i = 1; i <= HBKG->GetNbinsX(); i++)
    {
        Double_t bkgcontent = HWWscaled->GetBinContent(i) + HWZscaled->GetBinContent(i) + HTopsscaled->GetBinContent(i) + HZjetsscaled->GetBinContent(i) + Hllll->GetBinContent(i) + Hllqq->GetBinContent(i) + HVVV->GetBinContent(i) + HW_jets->GetBinContent(i) + HZtt->GetBinContent(i) + Hlllljj->GetBinContent(i) + Hllvvjj_WW->GetBinContent(i) + HWZ_jj->GetBinContent(i);
        HBKG->SetBinContent(i, bkgcontent);
    }*/

    /*HBKG->Add(HWWscaled);
    HBKG->Add(HWZscaled);
    HBKG->Add(HTopsscaled);
    HBKG->Add(HZjetsscaled);
    HBKG->Add(Hllll);
    HBKG->Add(Hllqq);
    HBKG->Add(HVVV);
    HBKG->Add(HW_jets);
    HBKG->Add(HZtt);
    HBKG->Add(Hlllljj);
    HBKG->Add(Hllvvjj_WW);
    HBKG->Add(HWZ_jj);
    TH1F *hBKG = (TH1F *)HBKG->Clone("hBKG");
    TH1F *h01 = (TH1F *)HBKG->Clone("h01");*/

    HBKG->Add(HWW);
    HBKG->Add(HWZ);
    HBKG->Add(HTops);
    HBKG->Add(HZjets0);
    HBKG->Add(HZjets1);
    HBKG->Add(HZjets2);
    HBKG->Add(Hllll);
    HBKG->Add(Hllqq);
    HBKG->Add(HVVV);
    HBKG->Add(HW_jets);
    HBKG->Add(HZtt);
    HBKG->Add(Hlllljj);
    HBKG->Add(Hllvvjj_WW);
    HBKG->Add(HWZ_jj);
    TH1F *hBKG = (TH1F *)HBKG->Clone("hBKG");
    TH1F *h01 = (TH1F *)HBKG->Clone("h01");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // One method to produce correct HSignal with correct Integral = calculated Signal
    // HSignal->Add(HDATA);
    // HSignal->Add(HBKG, -1);
    // One method to produce correct HSignal with correct Integral = calculated Signal
    // Hllvv->Add(Hllvvjj);

    // unscaled signal
    // HSignal->Add(Hllvv);

    // scaled signal ->1653.25 (cuts)
    // HSignal->Add(Hllvv, scale_Signal); // in order to be equal with the above HSignal calculation method (returns false entries... ~517k)

    HSignal->Add(Hllvv);
    HSignal->Add(Hllvvjj);
    // TH1F *hSignal_Sig = (TH1F *)HSignal->Clone("hSignal_Sig");
    // HSignal->Scale(scale_Signal);
    TH1F *h00 = (TH1F *)HSignal->Clone("h00");
    TH1F *hSignal = (TH1F *)HSignal->Clone("hSignal");

    // Draw the first histogram (HSignal) with a transparent red color
    HSignal->SetFillColorAlpha(38, 0.5); // 0.5 sets the transparency (0: fully transparent, 1: opaque)
    HSignal->SetTitle("");
    // Double_t scaling_HSignal = 1 / HSignal->Integral();
    // HSignal->Scale(scaling_HSignal);
    HSignal->Draw("hist");
    HSignal->SetStats(0);
    HSignal->GetYaxis()->SetTitle("Events");
    HSignal->GetYaxis()->SetTitleSize(0.045);
    // HSignal->SetMinimum(0);
    HSignal->SetStats(1111);
    HSignal->GetXaxis()->SetTitle("#DeltaR_{ll}");

    // Draw the second histogram (HBKG) on the same canvas with a transparent blue color
    HBKG->SetFillColorAlpha(46, 0.5); // Set the fill color to blue with 50% transparency
    HBKG->SetFillStyle(3344);
    // Double_t scaling_BKG = 1 / HBKG->Integral();
    // HBKG->Scale(scaling_BKG);
    HBKG->Draw("hist SAMES");
    HBKG->SetStats(1111);

    // Add a legend to distinguish between the histograms
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the position of the legend as needed
    legend->AddEntry(HSignal, "Signal", "f");
    legend->AddEntry(HBKG, "Background", "f");
    legend->Draw();

    TLatex *tex = new TLatex(0.1, 0.75, "#intL dt = 138.9 fb^{-1}, #sqrt{s} = 13 TeV");
    tex->SetNDC();
    tex->SetTextFont(12);
    tex->SetTextSize(0.03); // Adjust the value to make the text smaller
    tex->SetLineWidth(2);
    tex->Draw();

    TLatex *tex2 = new TLatex(0.1, 0.70, "Signal Region (SR)");
    tex2->SetNDC();
    tex2->SetTextFont(12);
    // tex->SetTextSize(0.03); // Adjust the value to make the text smaller
    tex2->SetLineWidth(2);
    tex2->Draw();

    // gStyle->SetOptStat(111);

    // Show the canvas
    pad_1->RedrawAxis();
    // c1->Update();
    //////////rearrange HSignal and HBKG without scale to get Significance ///////////////////////////////////////////
    /*TH1F *h00 = (TH1F *)HSignal->Clone("h00");
    TH1F *h01 = (TH1F *)HBKG->Clone("h01");
    TH1F *hBKG_Sig = (TH1F *)HBKG->Clone("hBKG_Sig");*/
    // h01->Scale(1 / scaling_BKG);
    // h00->Scale(1 / scaling_HSignal);

    Double_t Integ = 0.0;
    for (int i = 1; i <= HSignal->GetNbinsX(); i++)
    {
        Integ += HSignal->GetBinContent(i);
        printf("The values that h01 got are for bin : %d , content is : %f\n", i, HSignal->GetBinContent(i));
    }

    Double_t Integbkg = 0.0;
    for (int i = 1; i <= HBKG->GetNbinsX(); i++)
    {
        Integbkg += HBKG->GetBinContent(i);
        printf("The values that h01 got are for bin : %d , content is : %f\n", i, HBKG->GetBinContent(i));
    }

    Double_t integral = HSignal->Integral();
    cout << "Total Weighted Event Count (signal): " << integral << endl;

    Double_t integralbkg = HBKG->Integral();
    cout << "Total Weighted Event Count (Bkg): " << integralbkg << endl;

    printf("The integral of HSignal or equivalently ??? the number of entries are: %f\n", Integ);

    printf("The integral of HBKG or equivalently ??? the number of entries are: %f\n", Integbkg);
    double entriesSignal = HSignal->GetEntries();
    double entriesBkg = HBKG->GetEntries();
    printf("The real entries for the HSignal are %f and for HBKG are %f\n", entriesSignal, entriesBkg);

    pad_2->cd();

    // to plot significance

    // Double_t Z_bin=0;
    for (int bin = 1; bin <= HSignal->GetNbinsX(); ++bin)
    {
        // count from this bin and up AUTO XRISIMOPOIOUSA KANONIKA!!!!!!!!!!
        Double_t B = h01->Integral(bin, h01->GetSize());
        Double_t S = h00->Integral(bin, h00->GetSize());
        // bin per bin significance
        /*Double_t S = h00->GetBinContent(bin);
        Double_t B = h01->GetBinContent(bin);*/

        if (B > 0 && S > 0)
        {
            Double_t Z_bin = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
            // cout << "B=" << B << ", S=" << S << ", Z_bin=" << Z_bin << endl;
            HSignif->SetBinContent(bin, Z_bin);
        }
    }
    Double_t max = HSignif->GetMaximum();
    printf("The Signal Significance max was found : %f\n", max);

    HSignif->SetFillColorAlpha(8, 0.5);
    HSignif->SetFillStyle(3144);
    HSignif->SetMarkerColor(8);
    HSignif->SetLineColor(8);
    HSignif->Draw("");
    HSignif->SetStats(0);

    HSignif->GetXaxis()->SetTitleSize(0.1);
    HSignif->GetXaxis()->SetTitle("#DeltaR_{ll}");
    HSignif->GetXaxis()->SetTitleOffset(1);
    HSignif->GetXaxis()->SetLabelSize(0.1);

    HSignif->GetYaxis()->SetTitleOffset(0.35);
    HSignif->GetYaxis()->SetTitleSize(0.12);
    HSignif->GetYaxis()->SetLabelSize(0.1);
    //  HSignif->GetYaxis()->SetTitle("#frac{Data}{MC}");
    HSignif->GetYaxis()->SetTitle("Sig. Significance");

    pad_2->SetGrid();
    pad_2->RedrawAxis();
    c1->Update();
    // here we create copies of the above to scale them to achieve integral = 1
    TCanvas *c3 = new TCanvas("c3", "Canvadaki", 0., 0., 700, 700);
    TPad *pad_3 = new TPad("pad_3", "This is pad_3", 0.01, 0.30, 1., 1.);
    TPad *pad_4 = new TPad("pad_4", "This is pad_4", 0.01, 0.01, 1., 0.30);

    pad_3->SetBorderSize(0);
    pad_3->SetBottomMargin(0.02);
    pad_3->Draw();
    pad_4->SetBottomMargin(0.35);
    pad_4->SetTopMargin(0.0);
    pad_4->SetBorderSize(0);
    pad_4->Draw();
    pad_3->cd();
    Double_t scaling_HSignal = 1 / hSignal->Integral();
    hSignal->Scale(scaling_HSignal);

    hSignal->Draw("hist");

    hSignal->SetStats(1111);
    hSignal->GetYaxis()->SetTitle("Events");
    hSignal->GetXaxis()->SetTitle("M_{ll}");
    TH1F *hSignal_Sig = (TH1F *)HSignal->Clone("hSignal_Sig");
    hSignal->SetFillColorAlpha(38, 0.5);
    hSignal->SetMinimum(0);
    hSignal->SetTitle("");

    Double_t scaling_BKG = 1 / hBKG->Integral();
    hBKG->Scale(scaling_BKG);
    hBKG->Draw("hist SAMES");
    hBKG->SetFillColorAlpha(46, 0.5); // Set the fill color to blue with 50% transparency
    hBKG->SetFillStyle(3344);
    hBKG->SetStats(1111);
    TH1F *hBKG_Sig = (TH1F *)HBKG->Clone("hBKG_Sig");

    // Add a legend to distinguish between the histograms
    TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the position of the legend as needed
    legend1->AddEntry(hSignal, "Signal", "f");
    legend1->AddEntry(hBKG, "Background", "f");
    legend1->Draw();

    TLatex *tex1 = new TLatex(0.1, 0.75, "#intL dt = 138.9 fb^{-1}, #sqrt{s} = 13 TeV");
    tex1->SetNDC();
    tex1->SetTextFont(12);
    tex1->SetTextSize(0.03); // Adjust the value to make the text smaller
    tex1->SetLineWidth(2);
    tex1->Draw();

    TLatex *tex3 = new TLatex(0.1, 0.70, "(N-1) Optimal Signal Region (OSR)");
    tex3->SetNDC();
    tex3->SetTextFont(12);
    tex3->SetTextSize(0.03); // Adjust the value to make the text smaller
    tex3->SetLineWidth(2);
    tex3->Draw();

    // hSignal->GetXaxis()->SetLabelOffset(1);
    hSignal->GetYaxis()->SetTitleSize(0.045);
    // Show the canvas
    pad_3->RedrawAxis();

    pad_4->cd();

    // to plot significance

    Double_t Z_bin = 0;
    for (int bin_norm = 1; bin_norm <= hSignal_Sig->GetNbinsX(); ++bin_norm)
    {
        // count from this bin and up AUTO XRISIMOPOIOUSA KANONIKA!!!!!!!!!!
        Double_t B_norm = hBKG_Sig->Integral(bin_norm, hBKG_Sig->GetSize());
        Double_t S_norm = hSignal_Sig->Integral(bin_norm, hSignal_Sig->GetSize());
        // bin per bin significance
        /*Double_t S_norm = hSignal_Sig->GetBinContent(bin_norm);
        Double_t B_norm = hBKG_Sig->GetBinContent(bin_norm);*/

        if (B_norm > 0 && S_norm > 0)
        {
            Double_t Z_bin_norm = sqrt(2 * ((S_norm + B_norm) * log(1 + (S_norm / B_norm)) - S_norm));
            // cout << "B=" << B << ", S=" << S << ", Z_bin=" << Z_bin << endl;
            HSignif_norm->SetBinContent(bin_norm, Z_bin_norm);
        }
    }
    Double_t max_norm = HSignif_norm->GetMaximum();
    printf("The Normalized Signal Significance max was found : %f\n", max_norm);

    HSignif_norm->SetFillColorAlpha(8, 0.5);
    HSignif_norm->SetFillStyle(3144);
    HSignif_norm->SetMarkerColor(8);
    HSignif_norm->SetLineColor(8);
    HSignif_norm->Draw("");
    HSignif_norm->SetStats(0);

    HSignif_norm->GetXaxis()->SetTitleSize(0.1);
    HSignif_norm->GetXaxis()->SetTitle("E^{miss}_{T}/H_{T}");
    HSignif_norm->GetXaxis()->SetTitleOffset(1);
    HSignif_norm->GetXaxis()->SetLabelSize(0.1);

    HSignif_norm->GetYaxis()->SetTitleOffset(0.35);
    HSignif_norm->GetYaxis()->SetTitleSize(0.12);
    HSignif_norm->GetYaxis()->SetLabelSize(0.1);
    // HSignif->GetYaxis()->SetTitle("#frac{Data}{MC}");
    HSignif_norm->GetYaxis()->SetTitle("Sig. Significance");

    HSignif_norm->GetYaxis()->SetLabelSize(0.09);

    pad_4->SetGrid();
    pad_4->RedrawAxis();
    c3->Update();

    ///////////////////////////////////////////////// FIDUCIAL /////////////////////////////////////////////////////////////////////////////////
    Double_t Data_EWK_FID = 0.0, Data_EWK_FID_error = 0.0;
    Double_t Data_QCD_qq_FID = 0.0, Data_QCD_qq_FID_error = 0.0;
    Double_t Data_QCD_gg_FID = 0.0, Data_QCD_gg_FID_error = 0.0;
    for (const auto &dir_FID : dir)
    {
        cout << "=============" << dir_FID << "=============" << endl;

        for (const auto &root_FID : root)
        {
            std::string filepath_FID = "./SAMPLES" + dir_FID + root_FID;

            if (root_FID == "/EWK_SAMPLE.root")
            {
                cout << "===" << root_FID << "===" << endl;
                TFile *file = new TFile(filepath_FID.c_str());
                std::vector<float> events_data_EWK_FID;
                events_data_EWK_FID = DoReco_FID(file, dir_FID);
                Data_EWK_FID += events_data_EWK_FID.at(0);
                Data_EWK_FID_error += events_data_EWK_FID.at(1);
                file->Close();
            }

            if (root_FID == "/QCD_qq_SAMPLE.root")
            {
                cout << "===" << root_FID << "===" << endl;
                TFile *file = new TFile(filepath_FID.c_str());
                std::vector<float> events_data_QCD_qq_FID;
                events_data_QCD_qq_FID = DoReco_FID(file, dir_FID);
                Data_QCD_qq_FID += events_data_QCD_qq_FID.at(0);
                Data_QCD_qq_FID_error += events_data_QCD_qq_FID.at(1);
                file->Close();
            }

            if (root_FID == "/QCD_gg_SAMPLE.root")
            {
                cout << "===" << root_FID << "===" << endl;
                TFile *file = new TFile(filepath_FID.c_str());
                std::vector<float> events_data_QCD_gg_FID;
                events_data_QCD_gg_FID = DoReco_FID(file, dir_FID);
                Data_QCD_gg_FID += events_data_QCD_gg_FID.at(0);
                Data_QCD_gg_FID_error += events_data_QCD_gg_FID.at(1);
                file->Close();
            }
        }
    }
    printf("The number of Data for EWK_SAMPLE.root is %f +- %f\n", Data_EWK_FID, Data_EWK_FID_error);
    printf("The number of Data for QCD_qq_SAMPLE is %f +- %f\n", Data_QCD_qq_FID, Data_QCD_qq_FID_error);
    printf("The number of Data for QCD_gg_SAMPLE.root is %f +- %f\n", Data_QCD_gg_FID, Data_QCD_gg_FID_error);

    ///////////////////////////////////////////////////////// Cross-sections /////////////////////////////////////////////////////////////////////////
    Double_t Lumi_error = 0.017 * 139;
    Double_t C_ZZ_qq = (N_SR / Data_QCD_qq_FID);
    Double_t Cross_qq = (N_SR) / (139 * C_ZZ_qq); // units of fbarn

    Double_t C_ZZ_gg = (N_SR / Data_QCD_gg_FID);
    Double_t Cross_gg = (N_SR) / (139 * C_ZZ_gg); // units of fbarn

    Double_t C_ZZ_EWK = (N_SR / Data_EWK_FID);
    Double_t Cross_EWK = (N_SR) / (139 * C_ZZ_EWK); // units of fbarn

    Double_t summation = Data_EWK_FID + Data_QCD_gg_FID + Data_QCD_qq_FID;
    Double_t summation_er = sqrt(pow(Data_EWK_FID_error, 2) + pow(Data_QCD_gg_FID_error, 2) + pow(Data_QCD_qq_FID_error, 2));
    cout << " the total MC signal for EWK + QCD_gg and QCD_qq is: " << summation << "+-" << summation_er << endl;
    Double_t C_ZZ_all = (N_SR / summation);
    Double_t C_ZZ_all_er = sqrt(pow(N_SR_error / summation, 2) + pow((N_SR * summation_er / (summation * summation)), 2));
    Double_t Cross_all = (N_SR) / (139 * C_ZZ_all);
    Double_t partial_N_SR = 1.0 / (139 * C_ZZ_all);
    Double_t partial_Lumi = -N_SR / (139 * 139 * C_ZZ_all);
    Double_t partial_C_ZZ_all = -N_SR / (139 * C_ZZ_all * C_ZZ_all);

    Double_t Cross_all_error = sqrt(pow(partial_N_SR * N_SR_error, 2) + pow(partial_Lumi * Lumi_error, 2) + pow(partial_C_ZZ_all * C_ZZ_all_er, 2));

    printf("Efficiency for qq is %f\n", C_ZZ_qq);
    printf("Efficiency for gg is %f\n", C_ZZ_gg);
    printf("Efficiency for EWK is %f\n", C_ZZ_EWK);
    printf("Efficiency for all is %f +- %f\n", C_ZZ_all, C_ZZ_all_er);
    printf("The calculated Cross-section for QCD_qq is %f fb\n", Cross_qq);
    printf("The calculated Cross-section for QCD_gg is %f fb\n", Cross_gg);
    printf("The calculated Cross-section for EWK is %f fb\n", Cross_EWK);
    printf("The calculated overall Cross-section for pp->ZZ->llvv is (the summation over the above cross-sections) %f fb\n", Cross_EWK + Cross_gg + Cross_qq);
    printf("The calculated overall Cross-section for pp->ZZ->llvv is (the summation over the data) %f +- %f fb\n", Cross_all, Cross_all_error);

    Double_t sf_signal = 1.132;
    Double_t sf_signal_er = 0.0315;
    Double_t Lum = 139;
    Double_t Lum_er = 2.363;
    Double_t final = (sf_signal * summation) / Lum;
    Double_t final_er = final * sqrt(pow(sf_signal_er / sf_signal, 2) + pow(summation_er / summation, 2) + pow(Lum_er / Lum, 2));
    cout << "Cross section including only statistical uncertainties is: " << final << " +- " << final_er << endl;

    TCanvas *c2 = new TCanvas("c2", "Canvadaki", 0., 0., 700, 700);
    c2->Divide(1, 3); // Divides the canvas into 1 column and 3 rows
    c2->cd(1);        // Select the first pad
    h02->SetLineColor(kRed);
    h02->SetFillColor(kRed);
    h02->Draw(); // Plot histogram h01 in the first pad
    c2->cd(2);   // Select the second pad
    h03->SetLineColor(kBlue);
    h03->SetFillColor(kBlue);
    h03->Draw(); // Plot histogram h
    c2->cd(3);
    h04->SetLineColor(kGreen);
    h04->SetFillColor(kGreen);
    h04->Draw();

    TCanvas *c4 = new TCanvas("c4", "Can", 0., 0., 700, 700);
    TPad *pad_5 = new TPad("pad_5", "This is pad_met", 0.01, 0.30, 1., 1.);
    TPad *pad_6 = new TPad("pad_6", "This is pad_met2", 0.01, 0.01, 1., 0.30);

    pad_5->SetBorderSize(0);
    pad_5->SetBottomMargin(0.02);
    pad_5->Draw();
    pad_6->SetBottomMargin(0.35);
    pad_6->SetTopMargin(0.0);
    pad_6->SetBorderSize(0);
    pad_6->Draw();
    pad_5->cd();

    // Prosthesi olwn twn others
    for (int i = 1; i <= Hllll->GetNbinsX(); i++)
    {
        double binContent = Hllll->GetBinContent(i) + /*Hllvvjj->GetBinContent(i)*/ +Hllqq->GetBinContent(i) + Hllvvjj_WW->GetBinContent(i) + HVVV->GetBinContent(i) + HW_jets->GetBinContent(i) + HWZ_jj->GetBinContent(i) + HZtt->GetBinContent(i) + Hlllljj->GetBinContent(i);
        HOther->SetBinContent(i, binContent);
    }
    TH1F *Hother = (TH1F *)HOther->Clone("Hother");
    TH1F *Hzjets0 = (TH1F *)HZjets0->Clone("Hzjets2");
    TH1F *Hzjets1 = (TH1F *)HZjets1->Clone("Hzjets1");
    TH1F *Hzjets2 = (TH1F *)HZjets2->Clone("Hzjets2");
    TH1F *Hww = (TH1F *)HWW->Clone("Hww");
    TH1F *Hwz = (TH1F *)HWZ->Clone("Hwz");
    TH1F *Htop = (TH1F *)HTops->Clone("Htop");

    Hzjets0->Scale(scaleZjets0);
    Hzjets1->Scale(scaleZjets1);
    Hzjets2->Scale(scaleZjets2);
    Hww->Scale(scaleWW);
    Hwz->Scale(scaleWZ);
    Htop->Scale(scaleAllTop);

    /*Hwz->Scale(1.012);
    Htop->Scale(1.01);
    Hww->Scale(1.464);
    Hzjets0->Scale(1.352);
    Hzjets1->Scale(1.322);
    Hzjets2->Scale(1.066);
    HSignal->Scale(1.006);*/

    HMC->Add(Hother);
    HMC->Add(Hzjets0);
    HMC->Add(Hzjets1);
    HMC->Add(Hzjets2);
    HMC->Add(Hww);
    HMC->Add(Hwz);
    HMC->Add(Htop);
    HMC->Add(HSignal);

    HOther->SetFillColor(32);
    // HZjets->SetFillColor(46);
    HZjets0->SetFillColor(46);
    HZjets1->SetFillColor(7);
    HZjets2->SetFillColor(8);
    HSignal->SetFillColor(41);
    HTops->SetFillColor(38);
    HWW->SetFillColor(28);
    HWZscaled->SetFillColor(93);
    HWZ->SetFillColor(93);

    HWZ->Scale(scaleWZ);
    HTops->Scale(scaleAllTop);
    HWW->Scale(scaleWW);
    HZjets0->Scale(scaleZjets0);
    HZjets1->Scale(scaleZjets1);
    HZjets2->Scale(scaleZjets2);
    // HSignal->Scale(scale_Signal);

    /*HWZ->Scale(1.012);
    HTops->Scale(1.01);
    HWW->Scale(1.464);
    HZjets0->Scale(1.352);
    HZjets1->Scale(1.322);
    HZjets2->Scale(1.066);
    HSignal->Scale(1.006);*/

    ////// apo edw kai katw allazw an thelw na nai scaled or unscaled
    // HZjets->Add(HOther);
    HZjets0->Add(HOther);
    HZjets1->Add(HZjets0);
    HZjets2->Add(HZjets1);
    HTops->Add(HZjets2);
    HWW->Add(HTops);
    HWZ->Add(HWW);
    HSignal->Add(HWZ);
    TH1F *h000 = (TH1F *)HSignal->Clone("h000");

    HSignal->Draw("hist");
    HWZ->Draw("hist same");
    HWW->Draw("hist same");
    HTops->Draw("hist same");
    HZjets2->Draw("hist same");
    HZjets1->Draw("hist same");
    HZjets0->Draw("hist same");
    HOther->Draw("hist same");
    /*HDATA->SetMarkerStyle(20);
    HDATA->SetMarkerSize(0.8);
    HDATA->SetLineColor(1);
    HDATA->Draw("sameE0X0");*/
    HMC->SetMarkerStyle(20);
    HMC->SetMarkerSize(0.8);
    HMC->SetLineColor(1);
    HMC->Draw("sameE0X0");

    // HSignal->GetXaxis()->SetTitle("M_jj");
    HSignal->GetYaxis()->SetTitle("Events");
    HSignal->GetYaxis()->SetTitleSize(0.045);
    HSignal->SetStats(0);
    HSignal->GetXaxis()->SetLabelOffset(-999);

    TLatex *tex5 = new TLatex(0.1, 0.75, "#intL dt = 138.9 fb^{-1}, #sqrt{s} = 13 TeV");
    tex5->SetNDC();
    tex5->SetTextFont(12);
    tex5->SetTextSize(0.03); // Adjust the value to make the text smaller
    tex5->SetLineWidth(2);
    tex5->Draw();

    TLatex *tex4 = new TLatex(0.1, 0.70, "Optimal Signal Region");
    tex4->SetNDC();
    tex4->SetTextFont(12);
    tex4->SetTextSize(0.03); // Adjust the value to make the text smaller
    tex4->SetLineWidth(2);
    tex4->Draw();

    TLegend *leg3 = new TLegend(0.65, 0.65, 0.90, 0.90, NULL, "brNDC");
    TLegendEntry *leg_entry3;

    leg_entry3 = leg3->AddEntry(HDATA, "'Data'", "lp");
    leg_entry3 = leg3->AddEntry(HSignal, "Signal", "f");
    leg_entry3 = leg3->AddEntry(HWZ, "WZ scaled", "f");
    leg_entry3 = leg3->AddEntry(HWW, "WW scaled", "f");
    leg_entry3 = leg3->AddEntry(HZjets0, "Zjets0 scaled", "f");
    leg_entry3 = leg3->AddEntry(HZjets1, "Zjets1 scaled", "f");
    leg_entry3 = leg3->AddEntry(HZjets2, "Zjets2 scaled", "f");
    leg_entry3 = leg3->AddEntry(HTops, "top scaled", "f");
    leg_entry3 = leg3->AddEntry(HOther, "Other", "f");

    leg3->SetLineColor(0);
    leg3->SetBorderSize(0); // Set the border size to 0 to remove the border
    leg3->SetMargin(0.5);   // Adjust the margin value to change the size of the legend box
    leg3->Draw();

    pad_6->cd();

    TH1F *h001 = (TH1F *)HMC->Clone("h001");
    h000->Sumw2(1);
    h001->Sumw2(1);

    h001->Divide(h000);
    h001->Draw("E0X0");
    h001->SetName("Fit");
    TF1 *pol0 = new TF1("pol0", "pol0", 2.5, 3.2);
    h001->Fit("pol0", "R", "", 2.5, 3.2);
    pol0->Draw("same");
    h001->GetYaxis()->SetRangeUser(0.0, 2.5);
    h001->GetYaxis()->SetTitleSize(0.10);
    h001->GetYaxis()->SetLabelSize(0.07);
    h001->GetXaxis()->SetTitleSize(0.10);
    h001->GetXaxis()->SetLabelSize(0.08);
    h001->GetYaxis()->SetTitleOffset(0.45);

    Double_t Integral = 0.0;
    for (int i = 1; i <= h001->GetNbinsX(); i++)
    {
        Integral += h001->GetBinContent(i);
        printf("The values that h01 got are for bin : %d , content is : %f\n", i, h001->GetBinContent(i));
        printf("The ratio Data/MC is : %f\n", HMC->GetBinContent(i) / HSignal->GetBinContent(i));
    }
    printf("The integral of Data/MC is: %f\n", Integral);

    TLine *line = new TLine(2.5, 1, 3.2, 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(kRed);
    line->Draw("same");

    h001->GetXaxis()->SetTitle("#Delta#phi(E^{miss}_{T}, Z)");
    h001->GetYaxis()->SetTitle("#frac{Data}{MC}");

    double chi2 = 0.0;
    int ndof = 0; // Degrees of freedom

    for (int i = 1; i <= h001->GetNbinsX(); i++)
    {
        double data_over_MC = h001->GetBinContent(i);
        double error = h001->GetBinError(i);

        // Check for zero errors to avoid division by zero
        if (error > 0.0)
        {
            double term = (data_over_MC - (pol0->GetParameter(0))) / error;
            chi2 += term * term;
            ndof++;
        }
    }

    double reduced_chi2 = chi2 / (ndof - 1);

    // Create a TLatex object to display the reduced chi-squared value
    /*TLatex *chi2_label = new TLatex(0.7, 0.85, Form("#chi^{2}/ndf = %.2f/%d = %.2f", chi2, ndof, reduced_chi2));
    chi2_label->SetNDC();
    chi2_label->SetTextFont(12);
    chi2_label->SetTextSize(0.1);
    chi2_label->SetLineWidth(2);
    chi2_label->Draw();*/

    cout << "Reduced chi-squared value is : " << chi2 << " / " << (ndof - 1) << " or " << reduced_chi2 << endl;

    c4->Update();

    // Stop the timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    std::cout << "Time taken for script to run was: " << duration / pow(10, 6) << "seconds" << std::endl;

    return 0;
}