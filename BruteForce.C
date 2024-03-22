// For more details see https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/ZZjj2l2vRun2#Data_and_MC_samples
// Script for plotting 15/06/20 Dimitrii Krasnopevtsev
#include <TRandom3.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include "TProfile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLegendEntry.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
using namespace std;

// gStyle->SetLegendFont(42);
// gStyle->SetPalette(1);

// Method to fill hitograms (MET and MET_significance) with events properly weighted
//

std::vector<float> DoReco(TFile *file, TH1F *Histogram, int MC, double testdMetZPhi, double testmet_tst, double testMetOHT, string directory, string fileName, int Harry)
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

        if (directory == "/SR/")
        {
            // Define the cuts using a boolean variable (Reco-cuts)
            bool cuts = (dMetZPhi > testdMetZPhi && met_tst > testmet_tst &&
                         MetOHT > testMetOHT);

            if (cuts)
            {
                signal += weight;
                signaler += weight * weight;
                counter += 1;
                if (fileName == "/Z_jets_ee.root")
                {
                    if (n_jets < 1 && Harry == 1)
                    {
                        signal_njets_lt_1 += weight;
                        signaler_njets_lt_1 += weight * weight;
                        Histogram->Fill(met_tst, weight);
                    }
                    if (n_jets == 1 && Harry == 2)
                    {
                        signal_njets_eq_1 += weight;
                        signaler_njets_eq_1 += weight * weight;
                        Histogram->Fill(met_tst, weight);
                    }
                    if (n_jets > 1 && Harry == 3)
                    {
                        signal_njets_gt_1 += weight;
                        signaler_njets_gt_1 += weight * weight;
                        Histogram->Fill(met_tst, weight);
                    }
                }
                if (fileName == "/Z_jets_mumu.root")
                {
                    if (n_jets < 1 && Harry == 1)
                    {
                        signal_njets_less_1 += weight;
                        signaler_njets_less_1 += weight * weight;
                        Histogram->Fill(met_tst, weight);
                    }

                    if (n_jets == 1 && Harry == 2)
                    {
                        signal_njets_just_1 += weight;
                        signaler_njets_just_1 += weight * weight;
                        Histogram->Fill(met_tst, weight);
                    }
                    if (n_jets > 1 && Harry == 3)
                    {
                        signal_njets_big_1 += weight;
                        signaler_njets_big_1 += weight * weight;
                        Histogram->Fill(met_tst, weight);
                    }
                }

                if (fileName != "/Z_jets_ee.root" && fileName != "/Z_jets_mumu.root")
                {
                    Histogram->Fill(met_tst, weight);
                }
            }
        }
    }

    /*cout << "N=" << signal << "+-" << sqrt(signaler) << endl; // signal yield and error on this yield
    printf("ENTRIES FOR EACH HISTOGRAM ARE: %f as well as times got into is : %d \n", Histogram->GetEntries(), counter);*/

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
std::vector<std::string> directories = {"/SR/"};
std::vector<std::string> fileNames = {"/DATA.root", "/llll.root", "/lllljj.root", "/llqq.root", "/llvv.root", "/llvvjj.root", "/llvvjj_WW.root", "/top.root", "/ttbarV_ttbarVV.root", "/VVV.root", "/W_jets.root", "/Wt.root", "/WW.root", "/WZ.root", "/WZ_jj.root", "/Z_jets_ee.root", "/Z_jets_mumu.root", "/Ztt.root"};

// the main method

int BruteForce()
{
    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    Float_t Z;
    Float_t maxsig = -1;
    Int_t a = -1;
    Int_t b = -1;
    Int_t c = -1;
    Int_t z = -1;
    Int_t metritiri = 0;

    double testdMetZPhi[19] = {2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0, 3.05, 3.1};
    double testmet_tst[1] = {70};
    double testMetOHT[16] = {0.3, 0.4, 0.43, 0.46, 0.49, 0.52, 0.55, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.75, 0.8};

    for (Int_t j = 0; j < 19; j++)
    {
        for (Int_t k = 0; k < 1; k++)
        {
            for (Int_t l = 0; l < 16; l++)
            {
                Float_t xbins[32] = {70, 75, 80, 85, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 115, 120, 130, 140, 160, 180, 200, 250, 300, 400, 500, 700, 900, 1100, 1300, 1500, 1700};

                TH1F *HDATA = new TH1F("HDATA", "", 31, xbins);
                TH1F *Hllll = new TH1F("Hllll", "", 31, xbins);
                TH1F *Hlllljj = new TH1F("Hlllljj", "", 31, xbins);
                TH1F *Hllqq = new TH1F("Hllqq", "", 31, xbins);
                TH1F *Hllvv = new TH1F("Hllvv", "", 31, xbins);
                TH1F *Hllvvjj = new TH1F("Hllvvjj", "", 31, xbins);
                TH1F *Hllvvjj_WW = new TH1F("Hllvvjj_WW", "", 31, xbins);
                TH1F *HTop = new TH1F("HTop", "", 31, xbins);
                TH1F *HttbarV_ttbarVV = new TH1F("HttbarV_ttbarVV", "", 31, xbins);
                TH1F *HVVV = new TH1F("HVVV", "", 31, xbins);
                TH1F *HW_jets = new TH1F("HW_jets", "", 31, xbins);
                TH1F *HWt = new TH1F("HWt", "", 31, xbins);
                TH1F *HWW = new TH1F("HWW", "", 31, xbins);
                TH1F *HWZ = new TH1F("HWZ", "", 31, xbins);
                TH1F *HWZ_jj = new TH1F("HWZ_jj", "", 31, xbins);
                TH1F *HZ_jets_ee = new TH1F("HZ_jets_ee", "", 31, xbins);
                TH1F *HZ_jets_ee_0 = new TH1F("HZ_jets_ee_0", "", 31, xbins);
                TH1F *HZ_jets_ee_1 = new TH1F("HZ_jets_ee_1", "", 31, xbins);
                TH1F *HZ_jets_ee_2 = new TH1F("HZ_jets_ee_2", "", 31, xbins);
                TH1F *HZ_jets_mumu = new TH1F("HZ_jets_mumu", "", 31, xbins);
                TH1F *HZ_jets_mumu_0 = new TH1F("HZ_jets_mumu_0", "", 31, xbins);
                TH1F *HZ_jets_mumu_1 = new TH1F("HZ_jets_mumu_1", "", 31, xbins);
                TH1F *HZ_jets_mumu_2 = new TH1F("HZ_jets_mumu_2", "", 31, xbins);
                TH1F *HZtt = new TH1F("HZtt", "", 31, xbins);
                TH1F *HZjets0 = new TH1F("HZjets0", "", 31, xbins);
                TH1F *HZjets1 = new TH1F("HZjets1", "", 31, xbins);
                TH1F *HZjets2 = new TH1F("HZjets2", "", 31, xbins);
                TH1F *HOther = new TH1F("HOther", "", 31, xbins);
                TH1F *HSignal = new TH1F("HSignal", "", 31, xbins);
                TH1F *HTops = new TH1F("HTops", "", 31, xbins);
                TH1F *HBKG = new TH1F("HBKG", "", 31, xbins);

                printf("------------ Combination: j = %d, k= %d, l = %d-----------------\n", j, k, l);

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
                            events_data = DoReco(file, HDATA, 0, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            Data += events_data.at(0);
                            Data_error += events_data.at(1);
                            file->Close();
                        }
                        if (fileName == "/llll.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_llll;
                            events_llll = DoReco(file, Hllll, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            llll += events_llll.at(0);
                            llll_error += events_llll.at(1);
                            fourl = llll;
                            fourl_er = llll_error;
                            file->Close();
                        }
                        if (fileName == "/lllljj.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_lllljj;
                            events_lllljj = DoReco(file, Hlllljj, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            lllljj += events_lllljj.at(0);
                            lllljj_error += events_lllljj.at(1);
                            fourljj = lllljj;
                            fourljj_er = lllljj_error;
                            file->Close();
                        }
                        if (fileName == "/llqq.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_llqq;
                            events_llqq = DoReco(file, Hllqq, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            llqq += events_llqq.at(0);
                            llqq_error += events_llqq.at(1);
                            doublelqq = llqq;
                            doublelqq_er = llqq_error;
                            file->Close();
                        }
                        if (fileName == "/llvv.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_llvv;
                            events_llvv = DoReco(file, Hllvv, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            llvv += events_llvv.at(0);
                            llvv_error += events_llvv.at(1);
                            doublelvv = llvv;
                            doublelvv_er = llvv_error;
                            file->Close();
                        }
                        if (fileName == "/llvvjj.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_llvvjj;
                            events_llvvjj = DoReco(file, Hllvvjj, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            llvvjj += events_llvvjj.at(0);
                            llvvjj_error += events_llvvjj.at(1);
                            doublelvvjj = llvvjj;
                            doublelvvjj_er = llvvjj_error;
                            file->Close();
                        }
                        if (fileName == "/llvvjj_WW.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_llvvjj_WW;
                            events_llvvjj_WW = DoReco(file, Hllvvjj_WW, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            llvvjjWW += events_llvvjj_WW.at(0);
                            llvvjjWW_error += events_llvvjj_WW.at(1);
                            doublelvvjjWW = llvvjjWW;
                            doublelvvjjWW_er = llvvjjWW_error;
                            file->Close();
                        }
                        if (fileName == "/top.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_top;
                            events_top = DoReco(file, HTop, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            top += events_top.at(0);
                            top_error += events_top.at(1);
                            TOP = top;
                            TOP_er = top_error;
                            file->Close();
                        }
                        if (fileName == "/ttbarV_ttbarVV.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_ttbarV_ttbarVV;
                            events_ttbarV_ttbarVV = DoReco(file, HttbarV_ttbarVV, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            ttbarVttbarVV += events_ttbarV_ttbarVV.at(0);
                            ttbarVttbarVV_error += events_ttbarV_ttbarVV.at(1);
                            ttbar = ttbarVttbarVV;
                            ttbar_er = ttbarVttbarVV_error;
                            file->Close();
                        }
                        if (fileName == "/VVV.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_VVV;
                            events_VVV = DoReco(file, HVVV, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            VVV += events_VVV.at(0);
                            VVV_error += events_VVV.at(1);
                            tripleV = VVV;
                            tripleV_er = VVV_error;
                            file->Close();
                        }
                        if (fileName == "/W_jets.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_W_jets;
                            events_W_jets = DoReco(file, HW_jets, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            Wjets += events_W_jets.at(0);
                            Wjets_error += events_W_jets.at(1);
                            WJ = Wjets;
                            WJ_er = Wjets_error;
                            file->Close();
                        }
                        if (fileName == "/Wt.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_Wt;
                            events_Wt = DoReco(file, HWt, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            Wt += events_Wt.at(0);
                            Wt_error += events_Wt.at(1);
                            WT = Wt;
                            WT_er = Wt_error;
                            file->Close();
                        }
                        if (fileName == "/WW.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_WW;
                            events_WW = DoReco(file, HWW, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            WW += events_WW.at(0);
                            WW_error += events_WW.at(1);
                            ww = WW;
                            ww_er = WW_error;
                            file->Close();
                        }
                        if (fileName == "/WZ.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_WZ;
                            events_WZ = DoReco(file, HWZ, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            WZ += events_WZ.at(0);
                            WZ_error += events_WZ.at(1);
                            wz = WZ;
                            wz_er = WZ_error;
                            file->Close();
                        }
                        if (fileName == "/WZ_jj.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_WZ_jj;
                            events_WZ_jj = DoReco(file, HWZ_jj, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            WZjj += events_WZ_jj.at(0);
                            WZjj_error += events_WZ_jj.at(1);
                            WZJJ = WZjj;
                            WZJJ_er = WZjj_error;
                            file->Close();
                        }
                        if (fileName == "/Z_jets_ee.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_Z_jets_ee;
                            events_Z_jets_ee = DoReco(file, HZ_jets_ee, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            Zjetsee += events_Z_jets_ee.at(0);
                            Zjetsee_error += events_Z_jets_ee.at(1);
                            file->Close();

                            if (directory == "/SR/")
                            {
                                TFile *file = new TFile(filepath.c_str());
                                std::vector<float> events_Z_jets_ee_0;
                                events_Z_jets_ee_0 = DoReco(file, HZ_jets_ee_0, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 1);
                                Zjetsee_0 += events_Z_jets_ee_0.at(2);
                                Zjetsee_error_0 += events_Z_jets_ee_0.at(3);
                                ZJEE = Zjetsee_0;
                                ZJEE_er = Zjetsee_error_0;

                                std::vector<float> events_Z_jets_ee_1;
                                events_Z_jets_ee_1 = DoReco(file, HZ_jets_ee_1, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 2);
                                Zjetsee_1 += events_Z_jets_ee_1.at(4);
                                Zjetsee_error_1 += events_Z_jets_ee_1.at(5);
                                ZJEE_1 = Zjetsee_1;
                                ZJEE_1_er = Zjetsee_error_1;

                                std::vector<float> events_Z_jets_ee_2;
                                events_Z_jets_ee_2 = DoReco(file, HZ_jets_ee_2, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 3);
                                Zjetsee_2 += events_Z_jets_ee_2.at(6);
                                Zjetsee_error_2 += events_Z_jets_ee_2.at(7);
                                ZJEE_2 = Zjetsee_2;
                                ZJEE_2_er = Zjetsee_error_2;

                                file->Close();
                            }
                        }
                        if (fileName == "/Z_jets_mumu.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_Z_jets_mumu;
                            events_Z_jets_mumu = DoReco(file, HZ_jets_mumu, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            Zjetsmumu += events_Z_jets_mumu.at(0);
                            Zjetsmumu_error += events_Z_jets_mumu.at(1);
                            file->Close();

                            if (directory == "/SR/")
                            {
                                TFile *file = new TFile(filepath.c_str());
                                std::vector<float> events_Z_jets_mumu_0;
                                events_Z_jets_mumu_0 = DoReco(file, HZ_jets_mumu_0, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 1);
                                Zjetsmumu_0 += events_Z_jets_mumu_0.at(8);
                                Zjetsmumu_error_0 += events_Z_jets_mumu_0.at(9);
                                ZJMM = Zjetsmumu_0;
                                ZJMM_er = Zjetsmumu_error_0;

                                std::vector<float> events_Z_jets_mumu_1;
                                events_Z_jets_mumu_1 = DoReco(file, HZ_jets_mumu_1, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 2);
                                Zjetsmumu_1 += events_Z_jets_mumu_1.at(10);
                                Zjetsmumu_error_1 += events_Z_jets_mumu_1.at(11);
                                ZJMM_1 = Zjetsmumu_1;
                                ZJMM_1_er = Zjetsmumu_error_1;

                                std::vector<float> events_Z_jets_mumu_2;
                                events_Z_jets_mumu_2 = DoReco(file, HZ_jets_mumu_2, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 3);
                                Zjetsmumu_2 += events_Z_jets_mumu_2.at(12);
                                Zjetsmumu_error_2 += events_Z_jets_mumu_2.at(13);
                                ZJMM_2 = Zjetsmumu_2;
                                ZJMM_2_er = Zjetsmumu_error_2;

                                file->Close();
                            }
                        }
                        if (fileName == "/Ztt.root")
                        {
                            cout << "===" << fileName << "===" << endl;
                            TFile *file = new TFile(filepath.c_str());
                            std::vector<float> events_Ztt;
                            events_Ztt = DoReco(file, HZtt, 1, testdMetZPhi[j], testmet_tst[k], testMetOHT[l], directory, fileName, 0);
                            Ztt += events_Ztt.at(0);
                            Ztt_error += events_Ztt.at(1);
                            ZTT = Ztt;
                            ZTT_er = Ztt_error;
                            file->Close();
                        }
                    }
                }

                // Some useful outputs
                double totalBKG = fourl + fourljj + doublelqq + doublelvvjjWW + TOP + ttbar + tripleV + WJ + WT + ww + wz + WZJJ + ZJEE + ZJEE_1 + ZJEE_2 + ZJMM + ZJMM_1 + ZJMM_2 + ZTT;
                double totalBKG_er = sqrt(pow(fourl_er, 2) + pow(fourljj_er, 2) + pow(doublelqq_er, 2) + pow(doublelvvjjWW_er, 2) + pow(TOP_er, 2) + pow(ttbar_er, 2) + pow(tripleV_er, 2) + pow(WJ_er, 2) + pow(WT_er, 2) + pow(ww_er, 2) + pow(wz_er, 2) + pow(WZJJ_er, 2) + pow(ZJEE_er, 2) + pow(ZJEE_1_er, 2) + pow(ZJEE_2_er, 2) + pow(ZJMM_er, 2) + pow(ZJMM_1_er, 2) + pow(ZJMM_2_er, 2) + pow(ZTT_er, 2));
                cout << "total Bkg=" << totalBKG << "+-" << totalBKG_er << endl;
                cout << "S/B=" << (doublelvv + doublelvvjj) / totalBKG << endl;
                Double_t S = doublelvv + doublelvvjj;
                Double_t B = totalBKG;
                Double_t Z = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
                cout << " Significance = " << Z << endl;

                gROOT->SetBatch(kTRUE); // That way nothing is gonna appear in terms of histograms
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

                // Adding Zjetsee1 + Zjetsmumu1
                for (int i = 1; i <= HZ_jets_ee_1->GetNbinsX(); i++)
                {
                    Double_t binsofZ1 = HZ_jets_ee_1->GetBinContent(i) + HZ_jets_mumu_1->GetBinContent(i);
                    HZjets1->SetBinContent(i, binsofZ1);
                }
                TH1F *h03 = (TH1F *)HZjets1->Clone("h03");

                // Adding Zjetsee2 + Zjetsmumu2
                for (int i = 1; i <= HZ_jets_ee_2->GetNbinsX(); i++)
                {
                    Double_t binsofZ2 = HZ_jets_ee_2->GetBinContent(i) + HZ_jets_mumu_2->GetBinContent(i);
                    HZjets2->SetBinContent(i, binsofZ2);
                }
                TH1F *h04 = (TH1F *)HZjets2->Clone("h04");

                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // Adding all top together
                for (int i = 1; i <= HTop->GetNbinsX(); i++)
                {
                    Double_t binContentT = HTop->GetBinContent(i) + HttbarV_ttbarVV->GetBinContent(i) + HWt->GetBinContent(i);
                    HTops->SetBinContent(i, binContentT);
                }

                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

                HSignal->Add(Hllvv);
                HSignal->Add(Hllvvjj);

                // Double_t Z_bin=0;
                for (int bin = 1; bin < HSignal->GetSize(); ++bin)
                {
                    // bin per bin significance
                    // Double_t B = Hmet_pt_Zjets->GetBinContent(bin);
                    // Double_t S = Hmet_pt_sigQCD_new->GetBinContent(bin)+Hmet_pt_sigEWK_new->GetBinContent(bin);

                    // count from this bin and up
                    Double_t B = HBKG->Integral(bin, HBKG->GetSize());
                    Double_t S = HSignal->Integral(bin, HSignal->GetSize());

                    if (B > 0 && S > 0)
                    {
                        Double_t Z_bin = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));

                        if (Z_bin > maxsig)
                        {
                            b = j;
                            c = k;
                            z = l;
                            maxsig = Z_bin;
                            cout << "The maximum Significance is: " << maxsig << endl;
                            // printf("Found for bin = %f\n",bin);
                            printf("So met_tst was: %f\n", xbins[bin - 1]);
                            printf("Combination j = %d, k = %d, l = %d\n", j, k, l);
                        }
                    }
                }

                delete c1;

                delete HDATA;
                delete Hllll;
                delete Hlllljj;
                delete Hllqq;
                delete Hllvv;
                delete Hllvvjj;
                delete Hllvvjj_WW;
                delete HTop;
                delete HttbarV_ttbarVV;
                delete HVVV;
                delete HW_jets;
                delete HWt;
                delete HWZ;
                delete HWW;
                delete HWZ_jj;
                delete HZ_jets_ee;
                delete HZ_jets_ee_0;
                delete HZ_jets_ee_1;
                delete HZ_jets_ee_2;
                delete HZ_jets_mumu;
                delete HZ_jets_mumu_0;
                delete HZ_jets_mumu_1;
                delete HZ_jets_mumu_2;
                delete HZtt;
                delete HZjets0;
                delete HZjets1;
                delete HZjets2;
                delete HOther;
                delete HSignal;
                delete HTops;
                delete HBKG;

                metritiri += 1;
                Printf("TO METRITIRIIIIII # %d #\n", metritiri);
            }
        }
    }

    printf("Maximum significance is %f and found for combination of j = %d, k = %d, l = %d \n", maxsig, b, c, z);
    printf("The best cuts that can be applied are: DMetZPhi > %f, met_tst > %f, MetOHT > %f \n", testdMetZPhi[b], testmet_tst[c], testMetOHT[z]);//Need to search the specific combination and find met_tst because it is derived by the histogram 
    // Stop the timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    std::cout << "Time taken by script to run: " << duration / pow(10, 6) << " seconds"
              << " or " << duration / (60 * pow(10, 6)) << " minutes" << std::endl;

    return 0;
}