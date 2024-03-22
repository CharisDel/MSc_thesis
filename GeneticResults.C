#include <TROOT.h>
#include <iostream>
#include <TFile.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"

using namespace std;

// Function to read data from a text file
void readData(const string &filename, vector<int> &generations, vector<double> &maxSignificance, vector<double> &avgSignificance)
{
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    string line;
    while (getline(file, line))
    {
        int generation;
        double maxSig, avgSig;
        if (sscanf(line.c_str(), "Generation %d: Max Significance = %lf, Avg Significance = %lf", &generation, &maxSig, &avgSig) == 3)
        {
            generations.push_back(generation);
            maxSignificance.push_back(maxSig);
            avgSignificance.push_back(avgSig);
        }
    }

    file.close();
}

int GeneticResults()
{
    vector<int> tournamentGenerations, rouletteGenerations;
    vector<double> tournamentMaxSignificance, rouletteMaxSignificance, tournamentAvgSignificance, rouletteAvgSignificance;

    // Read data from text files
    readData("./Tournament_Selection.txt", tournamentGenerations, tournamentMaxSignificance, tournamentAvgSignificance);
    readData("./rouletteWheel_Selection.txt", rouletteGenerations, rouletteMaxSignificance, rouletteAvgSignificance);

    // Create canvas and graphs
    TCanvas *c1 = new TCanvas("c1", "Avg Significance vs Generation", 1000, 600);
    const int numPoints = tournamentGenerations.size();
    double *xValues = new double[numPoints];
    std::copy(tournamentGenerations.begin(), tournamentGenerations.end(), xValues);
    double *yValues = new double[numPoints];
    std::copy(tournamentAvgSignificance.begin(), tournamentAvgSignificance.end(), yValues);

    const int numPointsR = rouletteGenerations.size();
    double *xValuesR = new double[numPointsR];
    std::copy(rouletteGenerations.begin(), rouletteGenerations.end(), xValuesR);
    double *yValuesR = new double[numPointsR];
    std::copy(rouletteAvgSignificance.begin(), rouletteAvgSignificance.end(), yValuesR);

    TGraph *tournamentAvgGraph = new TGraph(numPoints, xValues, yValues);
    TGraph *rouletteAvgGraph = new TGraph(numPointsR, xValuesR, yValuesR);

    const char *additionalText[] = {
        "Genetic Algorithm Features :",
        "Population: 20",
        "Generations: 10",
        "Uniform Cross Over",
        "Mutation Rate: 0.25",
        "Roulette Wheel Selection",
        "Tournament Selection (Tournament Size: 4)"
    };

    double LinePosition = 0.75; // Initial position for the first line 

    // Customize graphs
    tournamentAvgGraph->SetTitle("Avg Significance vs Generation");
    tournamentAvgGraph->SetMarkerStyle(20);
    tournamentAvgGraph->SetMarkerColor(kBlue);
    tournamentAvgGraph->SetLineColor(kBlue);
    rouletteAvgGraph->SetMarkerStyle(20);
    rouletteAvgGraph->SetMarkerColor(kRed);
    rouletteAvgGraph->SetLineColor(kRed);
    tournamentAvgGraph->GetXaxis()->SetTitle("Generation");
    tournamentAvgGraph->GetYaxis()->SetTitle("Average Significance");

    // Draw graphs
    tournamentAvgGraph->Draw("APL");
    rouletteAvgGraph->Draw("PL SAME");
    for (const char *text1 : additionalText) {
        TLatex *textLabel1 = new TLatex(0.1, LinePosition, text1);
        textLabel1->SetNDC();
        textLabel1->SetTextFont(12);
        textLabel1->SetTextSize(0.04);
        textLabel1->Draw();

        LinePosition -= 0.07; // Adjust this value to control the spacing between lines
    }
    c1->Update();    

    // Create canvas and graphs for Max Significance
    TCanvas *c2 = new TCanvas("c2", "Max Significance vs Generation", 1000, 600);
    double *yValuesT = new double[numPoints];
    std::copy(tournamentMaxSignificance.begin(), tournamentMaxSignificance.end(), yValuesT);

    double *yValuesRR = new double[numPointsR];
    std::copy(rouletteMaxSignificance.begin(), rouletteMaxSignificance.end(), yValuesRR);

    TGraph *tournamentMaxGraph = new TGraph(numPoints, xValues, yValuesT);
    TGraph *rouletteMaxGraph = new TGraph(numPointsR, xValuesR, yValuesRR);

    // Customize graphs
    tournamentMaxGraph->SetTitle("Max Significance vs Generation");
    tournamentMaxGraph->SetMarkerStyle(20);
    tournamentMaxGraph->SetMarkerColor(kBlue);
    tournamentMaxGraph->SetLineColor(kBlue);
    rouletteMaxGraph->SetMarkerStyle(20);
    rouletteMaxGraph->SetMarkerColor(kRed);
    rouletteMaxGraph->SetLineColor(kRed);
    tournamentMaxGraph->GetXaxis()->SetTitle("Generation");
    tournamentMaxGraph->GetYaxis()->SetTitle("Max Significance");

     // Add lines of text using TLatex
    double linePosition = 0.75; // Initial position for the first line

    // Draw graphs
    tournamentMaxGraph->Draw("APL");
    rouletteMaxGraph->Draw("PL SAME");
    for (const char *text : additionalText) {
        TLatex *textLabel = new TLatex(0.1, linePosition, text);
        textLabel->SetNDC();
        textLabel->SetTextFont(12);
        textLabel->SetTextSize(0.04);
        textLabel->Draw();

        linePosition -= 0.07; // Adjust this value to control the spacing between lines
    }
    c2->Update();    

    return 0;
}
