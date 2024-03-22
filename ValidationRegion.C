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
#include "TLatex.h"
#include "TPad.h"
#include "TLine.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
using namespace std;

// DoReco function to calculate events from trees

std::vector<float> DoReco(TFile *file, TH1F *Histogram, int MC, int Harry, string fileName)
{

    TH1::SetDefaultSumw2(kTRUE);

    TTree *tree = file->Get<TTree>("tree");

    Int_t nentries = (Int_t)tree->GetEntries();
    // cout << "number of entries" <<nentries <<endl;

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
    Double_t event_3CR = 0.;
    Double_t event_type = 0.;
    Double_t mTW = 0.;
    int counter = 0;
    Double_t Z_pT = 0.;

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
        "Z_pT",
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
    tree->SetBranchAddress("Z_pT", &Z_pT);

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

        if ((met_tst > 85 && MetOHT > 0.45 && MetOHT < 0.50) || (met_tst > 80 && met_tst < 85 && MetOHT > 0.45))
        {
            counter += 1;
            signal += weight;            // signal yield is sum of weights
            signaler += weight * weight; // keep Sum(weights^2) for calculating the error on the signal yield

            if (fileName == "/Z_jets_ee.root")
            {
                if (n_jets < 1 && Harry == 1)
                {
                    signal_njets_lt_1 += weight;
                    signaler_njets_lt_1 += weight * weight;
                    Histogram->Fill(Z_pT, weight);
                }
                if (n_jets == 1 && Harry == 2)
                {
                    signal_njets_eq_1 += weight;
                    signaler_njets_eq_1 += weight * weight;
                    Histogram->Fill(Z_pT, weight);
                }
                if (n_jets > 1 && Harry == 3)
                {
                    signal_njets_gt_1 += weight;
                    signaler_njets_gt_1 += weight * weight;
                    Histogram->Fill(Z_pT, weight);
                }
            }
            if (fileName == "/Z_jets_mumu.root")
            {
                if (n_jets < 1 && Harry == 1)
                {
                    signal_njets_less_1 += weight;
                    signaler_njets_less_1 += weight * weight;
                    Histogram->Fill(Z_pT, weight);
                }

                if (n_jets == 1 && Harry == 2)
                {
                    signal_njets_just_1 += weight;
                    signaler_njets_just_1 += weight * weight;
                    Histogram->Fill(Z_pT, weight);
                }
                if (n_jets > 1 && Harry == 3)
                {
                    signal_njets_big_1 += weight;
                    signaler_njets_big_1 += weight * weight;
                    Histogram->Fill(Z_pT, weight);
                }
            }
            if (fileName != "/Z_jets_ee.root" && fileName != "/Z_jets_mumu.root")
            {
                Histogram->Fill(Z_pT, weight);
            }
        }
    }

    cout << "N=" << signal << "+-" << sqrt(signaler) << endl; // signal yield and error on this yield
    printf("ENTRIES FOR EACH HISTOGRAM ARE: %f as well as times got into if : %d \n", Histogram->GetEntries(), counter);

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

// Change of directory//
std::vector<std::string> directories = {/*"/3lCR/",  "/emCR_B/" , "/emCR_A/" ,*/ "/Zjets/" /*, "/SR/"*/};
std::vector<std::string> fileNames = {"/DATA.root", "/llll.root", "/lllljj.root", "/llqq.root", "/llvv.root", "/llvvjj.root", "/llvvjj_WW.root", "/top.root", "/ttbarV_ttbarVV.root", "/VVV.root", "/W_jets.root", "/Wt.root", "/WW.root", "/WZ.root", "/WZ_jj.root", "/Z_jets_ee.root", "/Z_jets_mumu.root", "/Ztt.root"};

// the main method

int ValidationRegion()
{

    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    //////////// Simultaneous fit scaling factors ////////////////

    Double_t scaleWZ = 1.012;
    Double_t scaleWZ_error = 0.0225;
    Double_t scaleAllTop = 1.012;
    Double_t scaleAllTop_error = 0.0125;
    Double_t scaleZjets0 = 1.352;
    Double_t scaleZjets0_error = 0.0305;
    Double_t scaleZjets1 = 1.322;
    Double_t scaleZjets1_error = 0.0205; // na balw tis alles times twn scaling factor. pio aples. tis prwtes pou briskw apo to simultaneous!!!
    Double_t scaleZjets2 = 1.066;
    Double_t scaleZjets2_error = 0.017;
    Double_t scaleWW = 1.438;
    Double_t scaleWW_error = 0.087;
    Double_t scaleSignal = 1.006;
    Double_t scaleSignal_error = 0.0325;

    //////////// Per region scaling factors ////////////////
    // Double_t scaleWZ = 1.023562;
    // Double_t scaleWZ_error = 0.024493;
    // Double_t scaleAllTop = 1.011241;
    // Double_t scaleAllTop_error = 0.013766;
    // Double_t scaleZjets0 = 1.364025;
    // Double_t scaleZjets0_error = 0.084158;
    // Double_t scaleZjets1 = 1.327656;
    // Double_t scaleZjets1_error = 0.037128; // na balw tis alles times twn scaling factor. pio aples. tis prwtes pou briskw apo to simultaneous!!!
    // Double_t scaleZjets2 = 1.066263;
    // Double_t scaleZjets2_error = 0.030255;
    // Double_t scaleWW = 1.460575;
    // Double_t scaleWW_error = 0.099886;
    // Double_t scaleSignal = 1.00;
    // Double_t scaleSignal_error = 0.037566;

    // Float_t xbins[9] = {70, 95, 120, 145, 170, 195, 250, 350, 1000}; // met_tst
    Float_t xbins[14] = {75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 500}; // Z_pT
    // Float_t xbins[9] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2};
    ///// Histograms for the desired variables
    TH1F *HDATA = new TH1F("HDATA", "", 13, xbins);
    TH1F *Hllll = new TH1F("Hllll", "", 13, xbins);
    TH1F *Hlllljj = new TH1F("Hlllljj", "", 13, xbins);
    TH1F *Hllqq = new TH1F("Hllqq", "", 13, xbins);
    TH1F *Hllvv = new TH1F("Hllvv", "", 13, xbins);
    TH1F *Hllvvjj = new TH1F("Hllvvjj", "", 13, xbins);
    TH1F *Hllvvjj_WW = new TH1F("Hllvvjj_WW", "", 13, xbins);
    TH1F *HTop = new TH1F("HTop", "", 13, xbins);
    TH1F *HttbarV_ttbarVV = new TH1F("HttbarV_ttbarVV", "", 13, xbins);
    TH1F *HVVV = new TH1F("HVVV", "", 13, xbins);
    TH1F *HW_jets = new TH1F("HW_jets", "", 13, xbins);
    TH1F *HWt = new TH1F("HWt", "", 13, xbins);
    TH1F *HWW = new TH1F("HWW", "", 13, xbins);
    TH1F *HWZ = new TH1F("HWZ", "", 13, xbins);
    TH1F *HWZ_jj = new TH1F("HWZ_jj", "", 13, xbins);
    TH1F *HZ_jets_ee = new TH1F("HZ_jets_ee", "", 13, xbins);
    TH1F *HZ_jets_mumu = new TH1F("HZ_jets_mumu", "", 13, xbins);
    TH1F *HZtt = new TH1F("HZtt", "", 13, xbins);
    TH1F *HZ_jets_ee_0 = new TH1F("HZ_jets_ee_0", "", 13, xbins);
    TH1F *HZ_jets_ee_1 = new TH1F("HZ_jets_ee_1", "", 13, xbins);
    TH1F *HZ_jets_ee_2 = new TH1F("HZ_jets_ee_2", "", 13, xbins);
    TH1F *HZ_jets_mumu_0 = new TH1F("HZ_jets_mumu_0", "", 13, xbins);
    TH1F *HZ_jets_mumu_1 = new TH1F("HZ_jets_mumu_1", "", 13, xbins);
    TH1F *HZ_jets_mumu_2 = new TH1F("HZ_jets_mumu_2", "", 13, xbins);
    TH1F *HZjets0 = new TH1F("HZjets0", "", 13, xbins);
    TH1F *HZjets1 = new TH1F("HZjets1", "", 13, xbins);
    TH1F *HZjets2 = new TH1F("HZjets2", "", 13, xbins);

    TH1F *HWZscaled = new TH1F("HWZCscaled", "", 13, xbins);
    TH1F *HTopsscaled = new TH1F("HTopsscaled", "", 13, xbins);
    // TH1F *HWWscaled = new TH1F("HTopsscaled", "", 13, xbins);

    TH1F *HOther = new TH1F("HOther", "", 13, xbins);
    TH1F *HZjets = new TH1F("HZjets", "", 13, xbins);
    TH1F *HSignal = new TH1F("HSignal", "", 13, xbins);
    TH1F *HTops = new TH1F("HTops", "", 13, xbins);

    for (const auto &directory : directories)
    {

        Int_t i = 1;
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
                events_data = DoReco(file, HDATA, 0, 0, fileName);
                Data += events_data.at(0);
                Data_error += events_data.at(1);
                file->Close();
            }
            if (fileName == "/llll.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_llll;
                events_llll = DoReco(file, Hllll, 1, 0, fileName);
                llll += events_llll.at(0);
                llll_error += events_llll.at(1);
                file->Close();
            }
            if (fileName == "/lllljj.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_lllljj;
                events_lllljj = DoReco(file, Hlllljj, 1, 0, fileName);
                lllljj += events_lllljj.at(0);
                lllljj_error += events_lllljj.at(1);
                file->Close();
            }
            if (fileName == "/llqq.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_llqq;
                events_llqq = DoReco(file, Hllqq, 1, 0, fileName);
                llqq += events_llqq.at(0);
                llqq_error += events_llqq.at(1);
                file->Close();
            }
            if (fileName == "/llvv.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_llvv;
                events_llvv = DoReco(file, Hllvv, 1, 0, fileName);
                llvv += events_llvv.at(0);
                llvv_error += events_llvv.at(1);
                file->Close();
            }
            if (fileName == "/llvvjj.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_llvvjj;
                events_llvvjj = DoReco(file, Hllvvjj, 1, 0, fileName);
                llvvjj += events_llvvjj.at(0);
                llvvjj_error += events_llvvjj.at(1);
                file->Close();
            }
            if (fileName == "/llvvjj_WW.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_llvvjj_WW;
                events_llvvjj_WW = DoReco(file, Hllvvjj_WW, 1, 0, fileName);
                llvvjjWW += events_llvvjj_WW.at(0);
                llvvjjWW_error += events_llvvjj_WW.at(1);
                file->Close();
            }
            if (fileName == "/top.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_top;
                events_top = DoReco(file, HTop, 1, 0, fileName);
                top += events_top.at(0);
                top_error += events_top.at(1);
                file->Close();
            }
            if (fileName == "/ttbarV_ttbarVV.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_ttbarV_ttbarVV;
                events_ttbarV_ttbarVV = DoReco(file, HttbarV_ttbarVV, 1, 0, fileName);
                ttbarVttbarVV += events_ttbarV_ttbarVV.at(0);
                ttbarVttbarVV_error += events_ttbarV_ttbarVV.at(1);
                file->Close();
            }
            if (fileName == "/VVV.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_VVV;
                events_VVV = DoReco(file, HVVV, 1, 0, fileName);
                VVV += events_VVV.at(0);
                VVV_error += events_VVV.at(1);
                file->Close();
            }
            if (fileName == "/W_jets.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_W_jets;
                events_W_jets = DoReco(file, HW_jets, 1, 0, fileName);
                Wjets += events_W_jets.at(0);
                Wjets_error += events_W_jets.at(1);
                file->Close();
            }
            if (fileName == "/Wt.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_Wt;
                events_Wt = DoReco(file, HWt, 1, 0, fileName);
                Wt += events_Wt.at(0);
                Wt_error += events_Wt.at(1);
                file->Close();
            }
            if (fileName == "/WW.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_WW;
                events_WW = DoReco(file, HWW, 1, 0, fileName);
                WW += events_WW.at(0);
                WW_error += events_WW.at(1);
                file->Close();
            }
            if (fileName == "/WZ.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_WZ;
                events_WZ = DoReco(file, HWZ, 1, 0, fileName);
                WZ += events_WZ.at(0);
                WZ_error += events_WZ.at(1);
                file->Close();
            }
            if (fileName == "/WZ_jj.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_WZ_jj;
                events_WZ_jj = DoReco(file, HWZ_jj, 1, 0, fileName);
                WZjj += events_WZ_jj.at(0);
                WZjj_error += events_WZ_jj.at(1);
                file->Close();
            }
            if (fileName == "/Z_jets_ee.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_Z_jets_ee;
                events_Z_jets_ee = DoReco(file, HZ_jets_ee, 1, 0, fileName);
                Zjetsee += events_Z_jets_ee.at(0);
                Zjetsee_error += events_Z_jets_ee.at(1);
                file->Close();

                if (directory == "/Zjets/")
                {
                    TFile *file = new TFile(filepath.c_str());
                    std::vector<float> events_Z_jets_ee_0;
                    events_Z_jets_ee_0 = DoReco(file, HZ_jets_ee_0, 1, 1, fileName);
                    Zjetsee_0 += events_Z_jets_ee_0.at(2);
                    Zjetsee_error_0 += events_Z_jets_ee_0.at(3);

                    std::vector<float> events_Z_jets_ee_1;
                    events_Z_jets_ee_1 = DoReco(file, HZ_jets_ee_1, 1, 2, fileName);
                    Zjetsee_1 += events_Z_jets_ee_1.at(4);
                    Zjetsee_error_1 += events_Z_jets_ee_1.at(5);

                    std::vector<float> events_Z_jets_ee_2;
                    events_Z_jets_ee_2 = DoReco(file, HZ_jets_ee_2, 1, 3, fileName);
                    Zjetsee_2 += events_Z_jets_ee_2.at(6);
                    Zjetsee_error_2 += events_Z_jets_ee_2.at(7);

                    file->Close();
                }
            }
            if (fileName == "/Z_jets_mumu.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_Z_jets_mumu;
                events_Z_jets_mumu = DoReco(file, HZ_jets_mumu, 1, 0, fileName);
                Zjetsmumu += events_Z_jets_mumu.at(0);
                Zjetsmumu_error += events_Z_jets_mumu.at(1);
                file->Close();
                if (directory == "/Zjets/")
                {
                    TFile *file = new TFile(filepath.c_str());
                    std::vector<float> events_Z_jets_mumu_0;
                    events_Z_jets_mumu_0 = DoReco(file, HZ_jets_mumu_0, 1, 1, fileName);
                    Zjetsmumu_0 += events_Z_jets_mumu_0.at(8);
                    Zjetsmumu_error_0 += events_Z_jets_mumu_0.at(9);

                    std::vector<float> events_Z_jets_mumu_1;
                    events_Z_jets_mumu_1 = DoReco(file, HZ_jets_mumu_1, 1, 2, fileName);
                    Zjetsmumu_1 += events_Z_jets_mumu_1.at(10);
                    Zjetsmumu_error_1 += events_Z_jets_mumu_1.at(11);

                    std::vector<float> events_Z_jets_mumu_2;
                    events_Z_jets_mumu_2 = DoReco(file, HZ_jets_mumu_2, 1, 3, fileName);
                    Zjetsmumu_2 += events_Z_jets_mumu_2.at(12);
                    Zjetsmumu_error_2 += events_Z_jets_mumu_2.at(13);

                    file->Close();
                }
            }
            if (fileName == "/Ztt.root")
            {
                cout << "===" << fileName << "===" << endl;
                TFile *file = new TFile(filepath.c_str());
                std::vector<float> events_Ztt;
                events_Ztt = DoReco(file, HZtt, 1, 0, fileName);
                Ztt += events_Ztt.at(0);
                Ztt_error += events_Ztt.at(1);
                file->Close();
            }
        }

        //////// these are for all of 'em (3lCR,emCR_B,emCR_A,Zjets)//////////////
        Double_t SignalZZ = llvvjj + llvv;
        Double_t SignalZZ_error = sqrt(pow(llvvjj_error, 2) + pow(llvv_error, 2));
        Double_t SignalZZ_A = scaleSignal * SignalZZ;
        Double_t SignalZZ_A_error = abs(SignalZZ_A) * sqrt(pow(scaleSignal_error / scaleSignal, 2) + pow(SignalZZ_error / SignalZZ, 2));
        Double_t EWKZZ = llvvjj;
        Double_t QCDZZ = llvv;
        Double_t AllTop = top + ttbarVttbarVV + Wt;
        Double_t AllTop_error = sqrt(pow(top_error, 2) + pow(ttbarVttbarVV_error, 2) + pow(Wt_error, 2));
        Double_t llklp = llll + llqq + VVV + Wjets + Ztt + lllljj + llvvjjWW + WZjj;
        Double_t llklp_error = sqrt(pow(llll_error, 2) + pow(llqq_error, 2) + pow(VVV_error, 2) + pow(Wjets_error, 2) + pow(Ztt_error, 2) + pow(llvvjj_error, 2) + pow(llvvjjWW_error, 2) + pow(WZjj_error, 2));

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

        Double_t WW_A = scaleWW * WW;
        Double_t WW_A_error = abs(WW_A) * sqrt(pow(scaleWW_error / scaleWW, 2) + pow(WW_error / WW, 2));
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
        Double_t Zjets = Z0 + Z1 + Z2;
        Double_t Zjets_error = sqrt(pow(Z0_error, 2) + pow(Z1_error, 2) + pow(Z2_error, 2));
        Double_t BKG = WZ + AllTop + WW + llklp + Z0 + Z1 + Z2;
        Double_t BKG_error = sqrt(pow(WZ_error, 2) + pow(AllTop_error, 2) + pow(WW_error, 2) + pow(llklp_error, 2) + pow(Z0_error, 2) + pow(Z1_error, 2) + pow(Z2_error, 2));
        Double_t BKG_emCR_A = WZemCR_A + AllTop_A + WW_A + llklp + ZjetsemCR_A;
        Double_t BKG_emCR_A_error = sqrt(pow(WZemCR_A_error, 2) + pow(AllTop_A_error, 2) + pow(WW_A_error, 2) + pow(llklp_error, 2) + pow(ZjetsemCR_A_error, 2));
        Double_t MC = WZ + AllTop + WW + llklp + Z0 + Z1 + Z2 + SignalZZ;
        Double_t MC_error = sqrt(pow(WZ_error, 2) + pow(AllTop_error, 2) + pow(WW_error, 2) + pow(llklp_error, 2) + pow(Z0_error, 2) + pow(Z1_error, 2) + pow(Z2_error, 2) + pow(SignalZZ_error, 2));
        Double_t MC_A = WZemCR_A + AllTop_A + WW_A + llklp + scaledZ0 + scaledZ1 + scaledZ2 + SignalZZ_A;
        Double_t MC_A_error = sqrt(pow(WZemCR_A_error, 2) + pow(AllTop_A_error, 2) + pow(WW_A_error, 2) + pow(llklp_error, 2) + pow(scaledZ0_error, 2) + pow(scaledZ1_error, 2) + pow(scaledZ2_error, 2) + pow(SignalZZ_A_error, 2));
        Double_t DATAMC = Data / MC;
        Double_t dataMCerror_emCR_A = DATAMC * sqrt(pow(Data_error / Data, 2) + pow(MC_error / MC, 2));
        Double_t Contamination = SignalZZ / MC;
        Double_t Contamination_error = Contamination * sqrt(pow(SignalZZ_error / SignalZZ, 2) + pow(MC_error / MC, 2));
        if (directory == "/Zjets/")
        {
            printf("|                              |all                    |\n"
                   "|DATA                          |%.2f +- %.2f       |\n"
                   "|EWKZZ                         |%.2f +- %.2f           |\n"
                   "|QCDZZ                         |%.2f +- %.2f           |\n"
                   "|WZ                            |%.2f +- %.2f          |\n"
                   "|Zjets                         |%.2f +- %.2f           |\n"
                   "|Zjets0                        |%.2f +- %.2f         |\n"
                   "|Zjets1                        |%.2f +- %.2f         |\n"
                   "|Zjets2                        |%.2f +- %.2f         |\n"
                   "|top-ttbarV-ttbarVV-Wt         |%.2f +- %.2f       |\n"
                   "|WW                            |%.2f +- %.2f         |\n"
                   "|llll-llqq-VVV-Wjets-Ztt       |%.2f +- %.2f         |\n"
                   "|Signal(ZZ)                    |%.2f +- %.2f           |\n"
                   "|Background                    |%.2f +- %.2f       |\n"
                   "|Data/MC                       |%.2f +- %.2f           |\n"
                   "|Contamination (percent)       |%.2f +- %.2f           |\n",
                   Data, Data_error, EWKZZ, llvvjj_error, QCDZZ, llvv_error, WZ, WZ_error, Zjets, Zjets_error, Z0, Z0_error, Z1, Z1_error, Z2, Z2_error, AllTop, AllTop_error, WW, WW_error, llklp, llklp_error, SignalZZ, SignalZZ_error, BKG, BKG_error, DATAMC, dataMCerror_emCR_A, 100 * Contamination, Contamination_error);
        }
    }
    /////////////////////////////// WITHOUT SCALING FACTORS //////////////////////////////////////
    TCanvas *c1 = new TCanvas("c1", "Canvadi", 0., 0., 700, 700);

    TPad *pad_0 = new TPad("pad_0", "This is pad_met", 0.01, 0.30, 1., 1.);
    TPad *pad_1 = new TPad("pad_1", "This is pad_met2", 0.01, 0.01, 1., 0.30);

    pad_0->SetBorderSize(0);
    pad_0->SetBottomMargin(0.02);
    pad_0->Draw();
    pad_1->SetBottomMargin(0.35);
    pad_1->SetTopMargin(0.0);
    pad_1->SetBorderSize(0);
    pad_1->Draw();
    pad_0->cd();

    // Prosthesi olwn twn others
    for (int i = 1; i <= Hllll->GetNbinsX(); i++)
    {
        double binContent = Hllll->GetBinContent(i) + /*Hllvvjj->GetBinContent(i)*/ +Hllqq->GetBinContent(i) + Hllvvjj_WW->GetBinContent(i) + HVVV->GetBinContent(i) + HW_jets->GetBinContent(i) + HWZ_jj->GetBinContent(i) + HZtt->GetBinContent(i) + Hlllljj->GetBinContent(i);
        HOther->SetBinContent(i, binContent);
    }

    // Prosthesi olwn twn Zjets
    /*for (int i = 1; i <= HZ_jets_ee->GetNbinsX(); i++)
    {
        double binContentZ = HZ_jets_ee->GetBinContent(i) + HZ_jets_mumu->GetBinContent(i);
        HZjets->SetBinContent(i, binContentZ);
    }*/
    // Prothsesi twn Zjets0
    for (int i = 1; i <= HZ_jets_ee_0->GetNbinsX(); i++)
    {
        double binContentZ0 = HZ_jets_ee_0->GetBinContent(i) + HZ_jets_mumu_0->GetBinContent(i);
        HZjets0->SetBinContent(i, binContentZ0);
    }
    // Prothsesi twn Zjets1
    for (int i = 1; i <= HZ_jets_ee_1->GetNbinsX(); i++)
    {
        double binContentZ1 = HZ_jets_ee_1->GetBinContent(i) + HZ_jets_mumu_1->GetBinContent(i);
        HZjets1->SetBinContent(i, binContentZ1);
    }
    // Prothsesi twn Zjets2
    for (int i = 1; i <= HZ_jets_ee_2->GetNbinsX(); i++)
    {
        double binContentZ2 = HZ_jets_ee_2->GetBinContent(i) + HZ_jets_mumu_2->GetBinContent(i);
        HZjets2->SetBinContent(i, binContentZ2);
    }

    // Prosthesi gia to signal
    for (int i = 1; i <= Hllvvjj->GetNbinsX(); i++)
    {
        double binContentS = Hllvvjj->GetBinContent(i) + Hllvv->GetBinContent(i);
        HSignal->SetBinContent(i, binContentS);
    }

    // Prosthesi gia ta top events
    for (int i = 1; i <= HTop->GetNbinsX(); i++)
    {
        double binContentT = HTop->GetBinContent(i) + HttbarV_ttbarVV->GetBinContent(i) + HWt->GetBinContent(i);
        HTops->SetBinContent(i, binContentT);
    }

    // Correction for WZ histogram with scaling factor scaleWZ
    // HWZscaled->Scale(scaleWZ);

    // Correction for Top histogram with scaling factor scaleAllTop
    // HTopsscaled->Scale(scaleAllTop);

    // scaling for WW
    TH1F *HWWscaled = (TH1F *)HWW->Clone("HWWscaled");
    HWWscaled->Scale(scaleWW);

    // HWZ->Scale(scaleWZ);
    // HTops->Scale(scaleAllTop);
    // HZjets2->Scale(scaleZjets2);
    // HZjets1->Scale(scaleZjets1); // PR sfs or SF sfs are declared in the upper body od code
    // HZjets0->Scale(scaleZjets0);
    // HSignal->Scale(scaleSignal);
    // HWW->Scale(scaleWW);

    HOther->SetFillColor(32);
    HZjets0->SetFillColor(46);
    HZjets1->SetFillColor(7);
    HZjets2->SetFillColor(8);
    HSignal->SetFillColor(41);
    HTops->SetFillColor(38);
    HWW->SetFillColor(28);
    HWZ->SetFillColor(93);

    ////// apo edw kai katw allazw an thelw na nai scaled or unscaled
    HZjets0->Add(HOther);
    HZjets1->Add(HZjets0);
    HZjets2->Add(HZjets1);
    HTops->Add(HZjets2);
    HWW->Add(HTops);
    HWZ->Add(HWW);
    HSignal->Add(HWZ);

    HSignal->Draw("hist");
    HWZ->Draw("hist same");
    HWW->Draw("hist same");
    HTops->Draw("hist same");
    HZjets2->Draw("hist same");
    HZjets1->Draw("hist same");
    HZjets0->Draw("hist same");
    HOther->Draw("hist same");
    HDATA->SetMarkerStyle(20);
    HDATA->SetMarkerSize(0.8);
    HDATA->SetLineColor(1);
    HDATA->Draw("sameE0X0");

    // HSignal->GetXaxis()->SetTitle("M_jj");
    HSignal->GetYaxis()->SetTitle("Events");
    HSignal->GetYaxis()->SetTitleSize(0.045);
    HSignal->SetStats(0);
    HSignal->GetXaxis()->SetLabelOffset(-999);

    TLatex *tex = new TLatex(0.1, 0.75, "#intL dt = 138.9 fb^{-1}, #sqrt{s} = 13 TeV");
    tex->SetNDC();
    tex->SetTextFont(12);
    tex->SetTextSize(0.03); // Adjust the value to make the text smaller
    tex->SetLineWidth(2);
    tex->Draw();

    TLatex *tex2 = new TLatex(0.1, 0.70, "Validation Region (VR)");
    tex2->SetNDC();
    tex2->SetTextFont(12);
    tex->SetTextSize(0.03); // Adjust the value to make the text smaller
    tex2->SetLineWidth(2);
    tex2->Draw();

    TLegend *leg = new TLegend(0.65, 0.65, 0.90, 0.90, NULL, "brNDC");
    TLegendEntry *leg_entry;

    leg_entry = leg->AddEntry(HDATA, "Data", "lp");
    leg_entry = leg->AddEntry(HSignal, "Signal", "f");
    leg_entry = leg->AddEntry(HWZ, "WZ", "f");
    leg_entry = leg->AddEntry(HWW, "WW", "f");
    leg_entry = leg->AddEntry(HZjets0, "Zjets0", "f");
    leg_entry = leg->AddEntry(HZjets1, "Zjets1", "f");
    leg_entry = leg->AddEntry(HZjets2, "Zjets2", "f");
    leg_entry = leg->AddEntry(HTops, "top", "f");
    leg_entry = leg->AddEntry(HOther, "Other", "f");

    // leg_entry = leg->AddEntry(HDATA, "Data", "lp");
    // leg_entry = leg->AddEntry(HSignal, "Signal scaled", "f");
    // leg_entry = leg->AddEntry(HWZ, "WZ scaled", "f");
    // leg_entry = leg->AddEntry(HWW, "WW scaled", "f");
    // leg_entry = leg->AddEntry(HZjets0, "Zjets0 scaled", "f");
    // leg_entry = leg->AddEntry(HZjets1, "Zjets1 scaled", "f");
    // leg_entry = leg->AddEntry(HZjets2, "Zjets2 scaled", "f");
    // leg_entry = leg->AddEntry(HTops, "top scaled", "f");
    // leg_entry = leg->AddEntry(HOther, "Other", "f");

    leg->SetLineColor(0);
    leg->SetBorderSize(0); // Set the border size to 0 to remove the border
    leg->SetMargin(0.5);   // Adjust the margin value to change the size of the legend box
    leg->Draw();

    pad_1->cd();

    TH1F *h00 = (TH1F *)HSignal->Clone("h00");
    TH1F *h01 = (TH1F *)HDATA->Clone("h01");
    h00->Sumw2(1);
    h01->Sumw2(1);

    h01->Divide(h00);
    h01->Draw("E0X0");
    h01->SetName("Fit");
    TF1 *pol0 = new TF1("pol0", "pol0", 75, 500);
    h01->Fit("pol0", "R", "", 75, 500);
    pol0->Draw("same");
    h01->SetStats(1111);
    h01->GetYaxis()->SetRangeUser(0.0, 2.5);
    h01->GetYaxis()->SetTitleSize(0.10);
    h01->GetYaxis()->SetLabelSize(0.07);
    h01->GetXaxis()->SetTitleSize(0.10);
    h01->GetXaxis()->SetLabelSize(0.08);
    h01->GetYaxis()->SetTitleOffset(0.45);

    Double_t Integ = 0.0;
    for (int i = 1; i <= h01->GetNbinsX(); i++)
    {
        Integ += h01->GetBinContent(i);
        printf("The values that h01 got are for bin : %d , content is : %f\n", i, h01->GetBinContent(i));
        printf("The ratio Data/MC is : %f\n", HDATA->GetBinContent(i) / HSignal->GetBinContent(i));
    }
    printf("The integral of Data/MC is: %f\n", Integ);

    TLine *line = new TLine(75, 1, 500, 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(kRed);
    line->Draw("same");

    h01->GetXaxis()->SetTitle("p_{T}^{Z} [GeV]");
    h01->GetYaxis()->SetTitle("#frac{Data}{MC}");

    double chi2 = 0.0;
    int ndof = 0; // Degrees of freedom

    for (int i = 1; i <= h01->GetNbinsX(); i++)
    {
        double data_over_MC = h01->GetBinContent(i);
        double error = h01->GetBinError(i);

        // Check for zero errors to avoid division by zero
        if (error > 0.0)
        {
            double term = (data_over_MC - (pol0->GetParameter(0))) / error;
            chi2 += term * term;
            ndof++;
        }
        //cout << "the error is " << error << endl;
    }

    double reduced_chi2 = chi2 / (ndof-1);
    Double_t redchi2c = pol0->GetChisquare();

    // Create a TLatex object to display the reduced chi-squared value
    /*TLatex *chi2_label = new TLatex(0.7, 0.85, Form("#chi^{2}/ndf = %.2f/%d = %.2f", chi2, ndof, reduced_chi2));
    chi2_label->SetNDC();
    chi2_label->SetTextFont(12);
    chi2_label->SetTextSize(0.1);
    chi2_label->SetLineWidth(2);
    chi2_label->Draw();*/

    cout << "Reduced chi-squared value is : " << chi2 << " / " << (ndof-1) << " or " << reduced_chi2 << endl;
    cout << "chis square = " << redchi2c << endl;

    c1->Update();

    // Stop the timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    std::cout << "Time taken for script to run was: " << duration / pow(10, 6) << "seconds" << std::endl;

    return 0;
}