#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <vector>
#include "TTree.h"

void MergeRootFiles(const std::string& outputFileName, const std::vector<std::string>& inputFiles) {
    TChain chain("tree"); // Assuming the tree name is "tree"

    // Add input files to the TChain
    for (const std::string& inputFile : inputFiles) {
        chain.Add(inputFile.c_str());
    }

    // Create an output ROOT file
    TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");

    // Clone the TChain and write it to the output file
    TTree* outputTree = chain.CloneTree(-1, "fast");
    outputTree->Write();

    // Close the output file
    outputFile->Close();
    delete outputFile;
}

int Merging() {
    // List of ROOT files to merge
    std::vector<std::string> fileNames = {
        "DATA.root",
        "llll.root",
        "lllljj.root",
        "llqq.root",
        "llvv.root",
        "llvvjj.root",
        "llvvjj_WW.root",
        "top.root",
        "ttbarV_ttbarVV.root",
        "VVV.root",
        "W_jets.root",
        "Wt.root",
        "WW.root",
        "WZ.root",
        "WZ_jj.root",
        "Z_jets_ee.root",
        "Z_jets_mumu.root",
        "Ztt.root"
    };

    std::string outputDir = "./SAMPLES/emuNRCR";

    for (const std::string& fileName : fileNames) {
        std::string outputFileName = outputDir + "/" + fileName;

        std::vector<std::string> inputFiles = {
            "./SAMPLES/emCR_A/" + fileName,
            "./SAMPLES/emCR_B/" + fileName
        };

        MergeRootFiles(outputFileName, inputFiles);

        std::cout << "Merged " << fileName << " to " << outputFileName << std::endl;
    }

    return 0;
}

