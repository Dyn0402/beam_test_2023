#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "RtypesCore.h"

void convertTree(const char* inputFileName = "ftest.root",
                 const char* outputFileName = "ftest_vec.root",
                 int maxSamples = -1,
                 int maxChannels = -1) {

    TFile fin(inputFileName);
    TTree* inputTree = dynamic_cast<TTree*>(fin.Get("nt"));
    if (!inputTree) {
        std::cerr << "Error: Input tree not found." << std::endl;
        return;
    }

    ULong64_t timestamp = 0;
    ULong64_t delta_timestamp = 0;
    UShort_t  fine_timestamp = 0;
    ULong64_t eventID = 0;
    std::vector<UShort_t>* sample = nullptr;
    std::vector<UShort_t>* channel = nullptr;
    std::vector<UShort_t>* amplitude = nullptr;

    inputTree->SetBranchAddress("eventId", &eventID);
    inputTree->SetBranchAddress("timestamp", &timestamp);
    inputTree->SetBranchAddress("delta_timestamp", &delta_timestamp);
    inputTree->SetBranchAddress("ftst", &fine_timestamp);
    inputTree->SetBranchAddress("sample", &sample);
    inputTree->SetBranchAddress("channel", &channel);
    inputTree->SetBranchAddress("amplitude", &amplitude);

    // --- Determine maxSamples and maxChannels if needed ---
    if (maxSamples == -1 || maxChannels == -1) {
        Long64_t nEntries = inputTree->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            inputTree->GetEntry(i);
            for (size_t j = 0; j < sample->size(); ++j) {
                if ((*sample)[j] > maxSamples) maxSamples = (*sample)[j];
                if ((*channel)[j] > maxChannels) maxChannels = (*channel)[j];
            }
        }
    }

    std::cout << "maxSamples = " << maxSamples << ", maxChannels = " << maxChannels << std::endl;

    // --- Create output file and tree ---
    TFile fout(outputFileName, "recreate");
    TTree nt("nt", "nt");

    // ✅ Use ROOT typedefs for the array as well
    std::vector<std::vector<UShort_t>> ampVec(maxChannels + 1, std::vector<UShort_t>(maxSamples + 1, 0));
    UShort_t (*amp)[10000] = nullptr;  // dummy pointer for TBranch

    // Convert 2D vector to C-style array for branch definition
    std::vector<std::vector<UShort_t>> amp2D(maxChannels + 1, std::vector<UShort_t>(maxSamples + 1));

    // allocate a static array dynamically on heap (safe for large size)
    UShort_t** ampArray = new UShort_t*[maxChannels + 1];
    for (int ch = 0; ch <= maxChannels; ++ch)
        ampArray[ch] = new UShort_t[maxSamples + 1];

    nt.Branch("eventId", &eventID, "eventId/l");
    nt.Branch("timestamp", &timestamp, "timestamp/l");
    nt.Branch("delta_timestamp", &delta_timestamp, "delta_timestamp/l");
    nt.Branch("ftst", &fine_timestamp, "ftst/s");

    // --- You can define a single flat array if you prefer: ---
    std::vector<UShort_t> amp_flat((maxChannels + 1) * (maxSamples + 1), 0);
    nt.Branch("amp", amp_flat.data(),
              Form("amp[%d][%d]/s", maxChannels + 1, maxSamples + 1));

    // --- Loop over entries and fill ---
    Long64_t nEntries = inputTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        inputTree->GetEntry(i);

        // Reset amplitudes
        std::fill(amp_flat.begin(), amp_flat.end(), 0);

        for (size_t j = 0; j < sample->size(); ++j) {
            int ch = (*channel)[j];
            int s = (*sample)[j];
            if (ch <= maxChannels && s <= maxSamples)
                amp_flat[ch * (maxSamples + 1) + s] = (*amplitude)[j];
        }

        nt.Fill();
    }

    fout.Write();
    fout.Close();

    // cleanup
    for (int ch = 0; ch <= maxChannels; ++ch)
        delete[] ampArray[ch];
    delete[] ampArray;
}

int main(int argc, char* argv[]) {
    const char* inputFileName = "ftest.root";
    const char* outputFileName = "ftest_array.root";
    int maxSamples = -1;
    int maxChannels = -1;

    if (argc >= 2) inputFileName = argv[1];
    if (argc >= 3) outputFileName = argv[2];
    if (argc >= 4) maxSamples = std::atoi(argv[3]);
    if (argc >= 5) maxChannels = std::atoi(argv[4]);

    convertTree(inputFileName, outputFileName, maxSamples, maxChannels);
    return 0;
}