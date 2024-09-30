#include "TDCProcessor.h"

// Constructor
TDCProcessor::TDCProcessor() : fFileName(nullptr), fWindowSize(0.1) {}

// Destructor
TDCProcessor::~TDCProcessor() {}

// Set file name
void TDCProcessor::SetFileName(const char* filename) {
    fFileName = filename;
}

// Set window size
void TDCProcessor::SetWindowSize(double window_size) {
    fWindowSize = window_size;
}

// Process the data
void TDCProcessor::Process() {
    if (fFileName == nullptr) {
        std::cerr << "File name not set!" << std::endl;
        return;
    }
    ProcessTDC();
}

// Process the TDC data
void TDCProcessor::ProcessTDC() {
    // Open the ROOT file in update mode
    TFile *file = TFile::Open(fFileName, "UPDATE");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << fFileName << std::endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("mrdmonitor_tree;1");
    if (!tree) {
        std::cerr << "Error accessing tree!" << std::endl;
        return;
    }

    // Set up branch for 'tdc'
    std::vector<double> *tdc = nullptr;
    tree->SetBranchAddress("tdc", &tdc);

    // Collect all 'tdc' values
    std::vector<double> all_values;
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        all_values.insert(all_values.end(), tdc->begin(), tdc->end());
    }

    // Sort the collected values
    std::sort(all_values.begin(), all_values.end());

    // Group data points that are close enough in time
    std::vector<double> averaged_values;
    if (!all_values.empty()) {
        double current_sum = all_values[0];
        int count = 1;
        double current_start = all_values[0];

        for (size_t i = 1; i < all_values.size(); ++i) {
            if (all_values[i] - current_start <= fWindowSize) {
                // Add to the current group
                current_sum += all_values[i];
                ++count;
            } else {
                // Average the current group and push to averaged_values
                averaged_values.push_back(current_sum / count);
                // Start a new group
                current_sum = all_values[i];
                count = 1;
                current_start = all_values[i];
            }
        }
        // Average the last group
        averaged_values.push_back(current_sum / count);
    }

    // Define histogram parameters
    double min_value = *std::min_element(averaged_values.begin(), averaged_values.end());
    double max_value = *std::max_element(averaged_values.begin(), averaged_values.end());
    int num_bins = static_cast<int>((max_value - min_value) / fWindowSize) + 1;

    // Create histogram for averaged values
    TH1D *average_hist = new TH1D("average_hist", "Averaged TDC values", num_bins, min_value, max_value);
    for (double value : averaged_values) {
        average_hist->Fill(value);
    }

    // Write the average histogram to the file
    file->cd(); // Ensure the file is the current directory
    average_hist->Write(); // Write the histogram to the file

    // Plot the results
    TCanvas *canvas = new TCanvas("canvas", "Average TDC Values", 800, 600);
    average_hist->Draw();
    canvas->SaveAs("average_tdc_plot.png");

    // Clean up
    file->Close(); // Save and close the file
    delete file;
    delete average_hist;
    delete canvas;
}

// Global function implementation
extern "C" void process_tdc(const char* filename, double window_size) {
    TDCProcessor processor;
    processor.SetFileName(filename);
    processor.SetWindowSize(window_size);
    processor.Process();
}

