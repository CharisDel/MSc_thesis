#include <TRandom3.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <vector>
#include <TLatex.h>
#include "TProfile.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <chrono> //Timer
using namespace std;
//
//
//
//
// MAIN
void zjets_splitter()
{
    vector<string> filenames = {"DATA", "WZ", "Z_jets_ee", "Z_jets_mumu",
                                "top", "ttbarV_ttbarVV", "Wt", "WW",
                                "llll", "llqq", "VVV", "W_jets", "Ztt",
                                "lllljj", "llvv", "llvvjj", "llvvjj_WW",
                                "WZ_jj"};

    string filepath = "./SAMPLES/Zjets/";

    for (string &filename : filenames)
    {
        cout << endl
             << "   Creating file " << filename << "..." << endl
             << endl;

        Double_t signal = 0.0;
        Double_t weight, n_jets;

        // Open the original TFile and TTree
        TFile *file = new TFile((string(filepath) + filename + ".root").c_str(), "READ");
        TTree *tree = (TTree *)file->Get("tree");

        // Create a new TFile and an empty output TTree
        TFile *file_output = new TFile(("./SAMPLES/Zjets3/" + filename + ".root").c_str(), "RECREATE");
        TTree *tree_output = tree->CloneTree(0); // Create an empty output tree with the same branches as the original tree

        // Retrieve the list of branches
        TObjArray *branches = tree->GetListOfBranches();
        tree->SetBranchAddress("n_jets", &n_jets);
        tree->SetBranchAddress("global_weight", &weight);

        Int_t eventCounter = 0;
        const Int_t batchSize = 1000; // Adjust the batch size as needed

        for (Int_t i = 0; i < tree->GetEntriesFast(); i++)
        {
            tree->GetEntry(i);
            if (n_jets > 1)
            {
                // Copy the desired event(s) from the original TTree to the new TTree
                tree->GetEntry(i);
                tree_output->Fill();
                signal = signal + weight;

                eventCounter++;

                // Save the output TTree periodically
                if (eventCounter % batchSize == 0)
                {
                    tree_output->AutoSave();
                    cout << "   Processed " << eventCounter << " events." << endl;
                }
            }
        }
        cout << "   TOTAL WEIGHT:   " << signal << endl;
        cout << "   TOTAL EVENTS:   " << tree_output->GetEntries() << endl;

        // Save the final output TTree
        tree_output->AutoSave();

        // Close the TFiles
        file_output->Close();
        file->Close();

        delete file_output;
        delete file;

        cout << "   -------------------------------------" << endl;
    }

    return;
}

