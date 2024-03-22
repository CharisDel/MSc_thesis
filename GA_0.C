#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TMath.h"
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <random>
#include <fstream>
using namespace std;

// DoReco function to calculate events from trees
std::vector<float> DoReco(TFile *file, const std::vector<double> &solution)
{
    TTree *tree = file->Get<TTree>("tree");

    Int_t nentries = (Int_t)tree->GetEntries();

    Double_t met_tst = 0.;
    Double_t dMetZPhi = 0.;
    Double_t MetOHT = 0.;
    //Double_t dLepR = 0.;
    Double_t weight = 0.0;

    std::vector<float> events;
    events.clear();

    vector<TString> branches = {
        //"M2Lep",
        "met_tst",
        "dMetZPhi",
        "MetOHT",
        //"dLepR",
        // "leading_pT_lepton",
        // "subleading_pT_lepton",
        // "Z_pT",
        // "n_jets",
        // "n_bjets",
        // "detajj",
        // "mjj",
        // "leading_jet_pt",
        // "second_jet_pt",
        //"event_3CR",
        //"event_type",
        "global_weight"};

    tree->SetBranchStatus("*", 0);

    for (const auto &branch : branches)
    {
        tree->SetBranchStatus(branch, 1);
    }

    tree->SetBranchAddress("met_tst", &met_tst);
    tree->SetBranchAddress("dMetZPhi", &dMetZPhi);
    tree->SetBranchAddress("MetOHT", &MetOHT);
    //tree->SetBranchAddress("dLepR", &dLepR);
    tree->SetBranchAddress("global_weight", &weight); // I don't use Data, so i don't need to set it independently as 1.

    double signal = 0.;
    double signaler = 0.;

    // Extract input values from the vector
    /*met_tst = inputValues[0];
    dLepR = inputValues[1];
    dMetZPhi = inputValues[2];
    MetOHT = inputValues[3];*/

    // Loop over events
    for (Int_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);
        if (met_tst > solution[0] && dMetZPhi > solution[1] /*&& dLepR < solution[2]*/ && MetOHT > solution[2])
        {
            signal += weight;            // signal yield is sum of weights
            signaler += weight * weight; // keep Sum(weights^2) for calculating the error on the signal yield
        }
    }
    events.push_back(signal);
    events.push_back(sqrt(signaler));

    return events;
}

// Significance calculation function
double Significance(const std::vector<double> &solution)
{
    // Directory to read from
    std::string directory = "/SR/";
    std::vector<std::string> fileNames = {"/DATA.root", "/llll.root", "/lllljj.root", "/llqq.root", "/llvv.root", "/llvvjj.root", "/llvvjj_WW.root", "/top.root", "/ttbarV_ttbarVV.root", "/VVV.root", "/W_jets.root", "/Wt.root", "/WW.root", "/WZ.root", "/WZ_jj.root", "/Z_jets_ee.root", "/Z_jets_mumu.root", "/Ztt.root"};

    // for significance /////
    double fourl, fourl_er, fourljj, fourljj_er, doublelqq, doublelqq_er, doublelvv, doublelvv_er, doublelvvjj, doublelvvjj_er, doublelvvjjWW, doublelvvjjWW_er, TOP, TOP_er, ttbar, ttbar_er, tripleV, tripleV_er, WJ, WJ_er, WT, WT_er, ww, ww_er, wz, wz_er, WZJJ, WZJJ_er, ZJEE, ZJEE_er, ZJEE_1, ZJEE_1_er, ZJEE_2, ZJEE_2_er, ZJMM, ZJMM_er, ZJMM_1, ZJMM_1_er, ZJMM_2, ZJMM_2_er, ZTT, ZTT_er, ZJETSEE, ZJETSEE_er, ZJETSMM, ZJETSMM_er;

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
            events_data = DoReco(file, solution);
            Data += events_data.at(0);
            Data_error += events_data.at(1);
            file->Close();
        }
        if (fileName == "/llll.root")
        {
            cout << "===" << fileName << "===" << endl;
            TFile *file = new TFile(filepath.c_str());
            std::vector<float> events_llll;
            events_llll = DoReco(file, solution);
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
            events_lllljj = DoReco(file, solution);
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
            events_llqq = DoReco(file, solution);
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
            events_llvv = DoReco(file, solution);
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
            events_llvvjj = DoReco(file, solution);
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
            events_llvvjj_WW = DoReco(file, solution);
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
            events_top = DoReco(file, solution);
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
            events_ttbarV_ttbarVV = DoReco(file, solution);
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
            events_VVV = DoReco(file, solution);
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
            events_W_jets = DoReco(file, solution);
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
            events_Wt = DoReco(file, solution);
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
            events_WW = DoReco(file, solution);
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
            events_WZ = DoReco(file, solution);
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
            events_WZ_jj = DoReco(file, solution);
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
            events_Z_jets_ee = DoReco(file, solution);
            Zjetsee += events_Z_jets_ee.at(0);
            Zjetsee_error += events_Z_jets_ee.at(1);
            ZJETSEE = Zjetsee;
            ZJETSEE_er = Zjetsee_error;
            file->Close();
        }
        if (fileName == "/Z_jets_mumu.root")
        {
            cout << "===" << fileName << "===" << endl;
            TFile *file = new TFile(filepath.c_str());
            std::vector<float> events_Z_jets_mumu;
            events_Z_jets_mumu = DoReco(file, solution);
            Zjetsmumu += events_Z_jets_mumu.at(0);
            Zjetsmumu_error += events_Z_jets_mumu.at(1);
            ZJETSMM = Zjetsmumu;
            ZJETSMM_er = Zjetsmumu_error;
            file->Close();
        }
        if (fileName == "/Ztt.root")
        {
            cout << "===" << fileName << "===" << endl;
            TFile *file = new TFile(filepath.c_str());
            std::vector<float> events_Ztt;
            events_Ztt = DoReco(file, solution);
            Ztt += events_Ztt.at(0);
            Ztt_error += events_Ztt.at(1);
            ZTT = Ztt;
            ZTT_er = Ztt_error;
            file->Close();
        }
    }

    // Some useful outputs
    Double_t totalBKG = fourl + fourljj + doublelqq + doublelvvjjWW + TOP + ttbar + tripleV + WJ + WT + ww + wz + WZJJ + ZJETSEE + ZJETSMM + ZTT;
    //Double_t totalBKG_er = sqrt(pow(fourl_er, 2) + pow(fourljj_er, 2) + pow(doublelqq_er, 2) + pow(doublelvvjjWW_er, 2) + pow(TOP_er, 2) + pow(ttbar_er, 2) + pow(tripleV_er, 2) + pow(WJ_er, 2) + pow(WT_er, 2) + pow(ww_er, 2) + pow(wz_er, 2) + pow(WZJJ_er, 2) + pow(ZJETSEE_er, 2) + pow(ZJETSMM_er, 2) + pow(ZTT_er, 2));
    cout << "B == total Bkg=" << totalBKG /*<< "+-" << totalBKG_er */<< endl;
    cout << "S/B=" << (doublelvv + doublelvvjj) / totalBKG << endl;
    Double_t S = doublelvv + doublelvvjj;
    //Double_t S_er = sqrt(pow(doublelvv_er,2) + pow(doublelvvjj_er,2));
    Double_t B = totalBKG;
    //Double_t B_er = totalBKG_er;
    // Calculate the significance if both B and S are greater than zero
    if (B > 0 && S > 0)
    {
        Double_t Z = sqrt(2 * ((S + B) * log(1 + (S / B)) - S));
        //Double_t Z_er = sqrt(2) * sqrt(pow(((log(1 + (S/B)) - (S/(S+B))) * S_er), 2) + pow((log(1 + (S/B)) + (S/(B * (S+ B)))) * B_er, 2));
        cout << " Significance = " << Z /*<< " +- " << Z_er */<< endl;
        return Z;
    }
    else
    {
        cout << "Significance calculation skipped: B or S is not greater than zero." << endl;
        return 0.0; // You can choose an appropriate value to represent "no significance."
    }
}

// Function to initialize a population of random solutions-cuts
std::vector<std::vector<double>> initializePopulation(int populationSize)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<int> dist_var1(90, 110);  // Integer range for the first variable
    std::uniform_int_distribution<int> dist_var2(222, 320); // Integer range for the second variable (2.15 - 3.20 * 100)
    //std::uniform_int_distribution<int> dist_var3(20, 185);  // Integer range for the third variable (0.20 - 1.85 * 100)
    std::uniform_int_distribution<int> dist_var4(50, 80);   // Integer range for the fourth variable (0.00 - 2.50 * 100)

    // Create initial population with random individuals
    std::vector<std::vector<double>> population;

    for (int i = 0; i < populationSize; ++i)
    {
        // Generate random values for each variable
        int random_var1 = dist_var1(gen); // Integer value
        double random_var2 = static_cast<double>(dist_var2(gen)) / 100.0;
        //double random_var3 = static_cast<double>(dist_var3(gen)) / 100.0;
        double random_var4 = static_cast<double>(dist_var4(gen)) / 100.0;

        // Create a vector with the generated random values
        std::vector<double> inputValues = {static_cast<double>(random_var1),
                                           random_var2,
                                           /*random_var3,*/
                                           random_var4};
        population.push_back(inputValues);
    }

    return population;
}

// Function to evaluate the fitness of each individual in the population
std::vector<double> evaluatePopulation(const std::vector<std::vector<double>> &population)
{
    int Counter = 1;
    std::vector<double> fitnessValues;
    for (const std::vector<double> &solution : population)
    {
        cout << "Iteration Number : " << Counter << endl;
        fitnessValues.push_back(Significance(solution));
        printf("Applied cuts were: met_tst > %f, dMetZPhi > %f, MetOHT > %f\n", solution[0], solution[1], solution[2]);
        Counter += 1;
    }
    return fitnessValues;
}

// Function to select parents based on fitness (simple roulette wheel selection)
std::pair<int, int> rouletteWheelSelection(const std::vector<std::vector<double>> &population, const std::vector<double> &fitnessValues)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> dis(fitnessValues.begin(), fitnessValues.end());

    int parentIndex1 = dis(gen);
    int parentIndex2 = dis(gen);

    return std::make_pair(parentIndex1, parentIndex2);
}

// Function to select parents using tournament selection
std::pair<int, int> selectParentsTournament(const std::vector<std::vector<double>> &population, const std::vector<double> &fitnessValues, int tournamentSize)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<> indexDis(0, population.size() - 1);

    int bestParentIndex1 = indexDis(gen);
    int bestParentIndex2 = indexDis(gen);

    for (int i = 1; i < tournamentSize; ++i)
    {
        int currentIndex = indexDis(gen);
        if (fitnessValues[currentIndex] > fitnessValues[bestParentIndex1])
        {
            bestParentIndex1 = currentIndex;
        }
    }

    for (int i = 1; i < tournamentSize; ++i)
    {
        int currentIndex = indexDis(gen);
        if (fitnessValues[currentIndex] > fitnessValues[bestParentIndex2])
        {
            bestParentIndex2 = currentIndex;
        }
    }

    return std::make_pair(bestParentIndex1, bestParentIndex2);
}

// Function to perform uniform crossover
std::pair<std::vector<double>, std::vector<double>> uniformCrossover(const std::vector<double> &parent1, const std::vector<double> &parent2)
{
    if (parent1.size() != parent2.size())
    {
        std::cerr << "Parent vectors must have the same size for uniform crossover." << std::endl;
        exit(1);
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);

    std::vector<double> offspring1(parent1.size());
    std::vector<double> offspring2(parent2.size());

    for (size_t i = 0; i < parent1.size(); ++i)
    {
        if (dis(gen) == 0)
        {
            offspring1[i] = parent1[i];
            offspring2[i] = parent2[i];
        }
        else
        {
            offspring1[i] = parent2[i];
            offspring2[i] = parent1[i];
        }
    }

    return std::make_pair(offspring1, offspring2);
}

// Function to apply mutation to a solution vector
void mutateSolution(std::vector<double> &solution, double mutationRate) // solution == child
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> dist_var1(-3, 3); // Mutation range for the first variable-cut
    std::uniform_real_distribution<double> dist_var234(-0.05, 0.05); // Mutation range for the second-third and forth variables-cuts

    for (double &value : solution)
    {
        if (std::rand() / double(RAND_MAX) < mutationRate)
        {
            if (&value == &solution[0])
            { // Apply mutation range for the first variable-cut
                value += dist_var1(gen);
            }
            else if (&value == &solution[1] || /*&value == &solution[2] ||*/ &value == &solution[2])
            { // Apply mutation range for the second-third-forth variables-cuts
                value += dist_var234(gen);
            }
        }
    }
}

int GA_0()
{
    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    int populationSize = 50;
    int numGenerations = 10;
    double mutationRate = 0.25;
    int tournamentSize = 4;

    // Create histograms
    //TH1F *histMaxSignificance = new TH1F("histMaxSignificance", "Max Significance vs Generation", numGenerations, 0.5, numGenerations + 0.5);
    //TH1F *histAvgSignificance = new TH1F("histAvgSignificance", "Average Significance vs Generation", numGenerations, 0.5, numGenerations + 0.5);

    std::vector<std::vector<double>> population = initializePopulation(populationSize);

    for (int generation = 0; generation < numGenerations; ++generation)
    {
        cout << "//////////////// This is Generation " << (generation + 1) << " /////////////////" << endl;

        double maxSignificance = 0.0, sumSignificane = 0.0;

        // Evaluate the fitness of the current population
        std::vector<double> fitnessValues = evaluatePopulation(population);

        std::vector<std::vector<double>> newPopulation;

        for (int i = 0; i < populationSize; ++i)
        {
            // Options for selection
            //std::pair<int, int> parentsIndices = rouletteWheelSelection(population, fitnessValues);

            std::pair<int, int> parentsIndices = selectParentsTournament(population, fitnessValues, tournamentSize);

            // Extract the selected parent vectors using the indices
            std::vector<double> parent1 = population[parentsIndices.first];
            std::vector<double> parent2 = population[parentsIndices.second];

            // Perform uniform crossover on the selected parents
            std::pair<std::vector<double>, std::vector<double>> offspring = uniformCrossover(parent1, parent2);

            // Choose one of the offspring randomly to add to the new population
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(0, 1);
            std::vector<double> child = dis(gen) == 0 ? offspring.first : offspring.second;

            // Apply mutation to the child
            mutateSolution(child, mutationRate);

            newPopulation.push_back(child);
        }
        // Calculate average and max Significance for each generation
        for (int i = 0; i < populationSize; i++)
        {
            sumSignificane += fitnessValues[i];
            if (fitnessValues[i] > maxSignificance)
            {
                maxSignificance = fitnessValues[i];
            }
        }
        cout << "Average Significance for generation: " << (generation + 1) << " is " << (sumSignificane / populationSize) << endl;
        cout << "Max Significance for generation: " << (generation + 1) << " is " << maxSignificance << endl;

        // Store the information in a text file in append mode
        std::ofstream outFile("./rouletteWheel_Selection.txt", std::ios_base::app);
        outFile << "Generation " << (generation + 1) << ": Max Significance = " << maxSignificance << ", Avg Significance = " << (sumSignificane / populationSize) << "\n";
        outFile.close();

        // Fill histograms
        // histMaxSignificance->SetBinContent(generation + 1, maxSignificance);
        // histAvgSignificance->SetBinContent(generation + 1, sumSignificane / populationSize);

        population = newPopulation;
    }

    // Create a canvas for max significance plot
    // TCanvas *canvasMaxSignificance = new TCanvas("canvasMaxSignificance", "Max Significance vs Generation", 800, 600);
    // Set marker size for histograms
    // histMaxSignificance->SetMarkerStyle(20); // Circle marker
    // histMaxSignificance->SetMarkerSize(1.5); // Adjust the size as needed
    // histMaxSignificance->GetXaxis()->SetTitle("Generation");
    // histMaxSignificance->GetYaxis()->SetTitle("Max Significance");

    // histMaxSignificance->Draw("P");

    // Create a canvas for average significance plot
    // TCanvas *canvasAvgSignificance = new TCanvas("canvasAvgSignificance", "Average Significance vs Generation", 800, 600);
    // Set marker size for histograms
    // histAvgSignificance->SetMarkerStyle(20); // Circle marker
    // histAvgSignificance->SetMarkerSize(1.5); // Adjust the size as needed
    // histAvgSignificance->GetXaxis()->SetTitle("Generation");
    // histAvgSignificance->GetYaxis()->SetTitle("Average Significance");

    // histAvgSignificance->Draw("P");

    // histAvgSignificance->Fit("pol1");
    // gStyle->SetOptFit(1011);

    // canvasAvgSignificance->Update();

    // Find the best solution in the final population
    std::vector<double> bestSolution = population[0]; // Initialize with the first solution
    cout << "/////////This is just to initiate a ''Best solution''- Best significance ////////// " << endl;
    double bestFitness = Significance(bestSolution);

    cout << "///////////////////// Extra " << population.size() << " Iterations in new population to derive the max significance (Like an extra generation)/////////////////////////" << endl;

    for (const std::vector<double> &solution : population)
    {
        double fitness = Significance(solution);
        printf("Applied cuts were: met_tst > %f, dMetZPhi > %f, MetOHT > %f\n", solution[0], solution[1], solution[2]);
        if (fitness > bestFitness)
        {
            bestSolution = solution;
            bestFitness = fitness;
        }
    }

    // Print the best solution and its fitness
    std::cout << "Best solution: ";
    for (double value : bestSolution)
    {
        std::cout << value << " ";
    }
    std::cout << std::endl;
    std::cout << "Best fitness - Max Significance : " << bestFitness << std::endl;

    auto end1 = std::chrono::high_resolution_clock::now();
    // Calculate the duration
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start).count();

    std::cout << "Time taken by script to run: " << duration1 / pow(10, 6) << " seconds"
              << " or " << duration1 / (60 * pow(10, 6)) << " minutes" << std::endl;

    return 0;
}