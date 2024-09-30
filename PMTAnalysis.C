/*///////////////////////////////////////////////////////////////
// Tool to be used for analysis of tank PMTs during laser runs //
// For the ANNIE Collaboration                                 //
// By: Eli Brunner-Huber                                       //
// Last revised: 07/18/24                                      //
///////////////////////////////////////////////////////////////*/



#define PMTAnalysis_cxx
#include "PMTAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TLegend.h>
#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TTree.h>

// Constructor
PMTAnalysis::PMTAnalysis() : fChain(nullptr), fFileName(""), fStartTime(0), fEndTime(0), fCutoffVoltage(0) {
    
    // Initialize histograms in the constructor

    allhitfreq = new TH1D("allhitfreq", "Combined PMT Hits", 132, 332, 464);
    allhitfreq->GetYaxis()->SetRangeUser(0, 130000);

    LUXhitfreq = new TH1D("LUXhitfreq", "LUX PMT Hit Frequencies", 132, 332, 464);
    LUXhitfreq->SetFillColor(kBlue);

    ETELhitfreq = new TH1D("ETELhitfreq", "ETEL PMT Hit Frequencies", 132, 332, 464);
    ETELhitfreq->SetFillColor(kCyan);

    HAMAMATSUhitfreq = new TH1D("HAMAMATSUhitfreq", "Hamamatsu PMT Hit Frequencies", 132, 332, 464);
    HAMAMATSUhitfreq->SetFillColor(kGreen);

    WATCHMANhitfreq = new TH1D("WATCHMANhitfreq", "Watchman PMT hit Frequencies", 132, 332, 464);
    WATCHMANhitfreq->SetFillColor(kOrange);

    WATCHBOYhitfreq = new TH1D("WATCHBOYhitfreq", "Watchboy PMT Hit Frequencies", 132, 332, 464);
    WATCHBOYhitfreq->SetFillColor(kYellow);

    RADIUShitfreq = new TH1D("RADIUShitfreq", "Hits as a function of radius", 100, 1, 2);
    THETAhitfreq = new TH1D("THETAhitfreq", "Hits as a function of polar angle", 180, 0, 180);
    PHIhitfreq = new TH1D("PHIhitfreq", "Hits as a function of azimuthal angle", 360, 0, 360);

    // Initialize PHD and timing plot histograms for all PMTs in the constructor
    for (int j = 332; j <= 464; ++j) {
        char name[50];
        sprintf(name, "PulseHeightDistribution_PMT%d", j);
        TH1D* hist = new TH1D(name, name, 10000, 0, 10);
        pulseHeightDistributions.push_back(hist);
    }
    
    for (int j = 332; j <= 464; ++j) {
        char name[50];
        sprintf(name, "TimingPlot_PMT%d", j);
        TH1D* histo = new TH1D(name, name, 10000, fStartTime, fEndTime);
        timingPlots.push_back(histo);
    }
}

// Destructor
PMTAnalysis::~PMTAnalysis() {
    // Delete histograms in the destructor
    if (LUXhitfreq) delete LUXhitfreq;
    if (ETELhitfreq) delete ETELhitfreq;
    if (HAMAMATSUhitfreq) delete HAMAMATSUhitfreq;
    if (WATCHMANhitfreq) delete WATCHMANhitfreq;
    if (WATCHBOYhitfreq) delete WATCHBOYhitfreq;
    if (RADIUShitfreq) delete RADIUShitfreq;
    if (THETAhitfreq) delete THETAhitfreq;
    if (PHIhitfreq) delete PHIhitfreq;
    for (auto hist : pulseHeightDistributions) {
        delete hist;
    }
    
    for (auto histo : timingPlots) {
        delete histo;
    }
            
    pulseHeightDistributions.clear();
    timingPlots.clear();
}

// Analyze function declaration and setup

void PMTAnalysis::Analyze(const char* fileName, double startTime, double endTime, double cutoffVoltage)
{
    // Store parameters in variables
    fFileName = fileName;
    fStartTime = startTime;
    fEndTime = endTime;
    fCutoffVoltage = cutoffVoltage;


    // Open the ROOT file containing your TTree or data structure
    TFile* file = TFile::Open(fFileName);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Failed to open input file: " << fFileName << std::endl;
        return;
    }

    fChain = dynamic_cast<TTree*>(file->Get("phaseIITriggerTree"));
    if (!fChain) {
        std::cerr << "Error: Failed to retrieve TTree 'phaseIITriggerTree' from file." << std::endl;
        file->Close();
        return;
    }

    // Variables to hold data from the TTree
    std::vector<double>* hitX = nullptr;
    std::vector<double>* hitY = nullptr;
    std::vector<double>* hitZ = nullptr;
    std::vector<int>* hitDetID = nullptr;
    std::vector<double>* hitQ = nullptr;
    std::vector<double>* hitT = nullptr;
    std::vector<double>* hitPE = nullptr;

    // Linking TBranches to variables
    fChain->SetBranchAddress("hitX", &hitX);
    fChain->SetBranchAddress("hitY", &hitY);
    fChain->SetBranchAddress("hitZ", &hitZ);
    fChain->SetBranchAddress("hitDetID", &hitDetID);
    fChain->SetBranchAddress("hitQ", &hitQ);
    fChain->SetBranchAddress("hitT", &hitT);
    fChain->SetBranchAddress("hitPE", &hitPE);

    Long64_t nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    
    // Start loop

    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = fChain->LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        for (int k = 0; k < hitX->size(); k++) {
            // Set conditions
            if ((*hitQ)[k] < fCutoffVoltage && (*hitT)[k] >= fStartTime && (*hitT)[k] <= fEndTime) {
                // Fill PHDs and timing plots
                int pmtID = (*hitDetID)[k];
                if (pmtID >= 332 && pmtID <= 464) {
                    pulseHeightDistributions[pmtID - 332]->Fill((*hitPE)[k]);
                    timingPlots[pmtID - 332]->Fill((*hitT)[k]);
                }

                // Fill other histograms based on hitDetID
                if (pmtID >= 332 && pmtID <= 351) {
                    LUXhitfreq->Fill(pmtID);
                } else if (pmtID >= 352 && pmtID <= 371) {
                    ETELhitfreq->Fill(pmtID);
                } else if (((pmtID >= 372 && pmtID <= 381) || (pmtID >= 383 && pmtID <= 392) ||
                            (pmtID >= 394 && pmtID <= 403) || (pmtID >= 406 && pmtID <= 415))) {
                    HAMAMATSUhitfreq->Fill(pmtID);
                } else if (pmtID >= 416 && pmtID <= 464) {
                    WATCHBOYhitfreq->Fill(pmtID);
                } else if (pmtID == 382 || pmtID == 393 || pmtID == 404 || pmtID == 405) {
                    WATCHMANhitfreq->Fill(pmtID);
                }

                // Fill spherical histograms based on hit coordinates
                std::vector<double> spherical = getSpherical((*hitX)[k], (*hitY)[k], (*hitZ)[k]);
                RADIUShitfreq->Fill(spherical[0]);
                THETAhitfreq->Fill(spherical[1]);
                PHIhitfreq->Fill(spherical[2]);
            }
        }
    }

    // Construct output file name with input file name included
    std::string outputFileName = "PMTAnalysis_" + std::string(fileName);
    TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Failed to create output file: " << outputFileName << std::endl;
        file->Close();
        return;
    }

    // Create a directory to store all histograms
    TDirectory* histDir = outputFile->mkdir("Histograms");
    histDir->cd();

    // Write all histograms
    LUXhitfreq->Write();
    ETELhitfreq->Write();
    HAMAMATSUhitfreq->Write();
    WATCHBOYhitfreq->Write();
    WATCHMANhitfreq->Write();
    RADIUShitfreq->Write();
    THETAhitfreq->Write();
    PHIhitfreq->Write();
    
    
    // Save all histograms to the same Canvas
    TCanvas* tc = new TCanvas("Combined Canvas", "Combined Hit Frequency Histograms of all PMTs");

    // Draw all histograms and remove stats box
    allhitfreq->SetStats(0);
    allhitfreq->Draw();
    LUXhitfreq->Draw("SAME");
    ETELhitfreq->Draw("SAME");
    HAMAMATSUhitfreq->Draw("SAME");
    WATCHBOYhitfreq->Draw("SAME");
    WATCHMANhitfreq->Draw("SAME");
    PHIhitfreq->Draw("SAME");
    THETAhitfreq->Draw("SAME");
    RADIUShitfreq->Draw("SAME");

    // Legend of combined histogram
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(LUXhitfreq, "LUX (top)", "f");
    legend->AddEntry(ETELhitfreq, "ETEL (bottom)", "f");
    legend->AddEntry(HAMAMATSUhitfreq, "Hamamatsu (bottom)", "f");
    legend->AddEntry(WATCHBOYhitfreq, "Watchboy (bottom)", "f");
    legend->AddEntry(WATCHMANhitfreq, "Watchman (tank)", "f");
    legend->Draw();
    tc->Write();

    // Create a directory to store pulse height distributions and timing plots
    TDirectory* phdDir = outputFile->mkdir("PulseHeightDistributions");
    TDirectory* timingDir = outputFile->mkdir("TimingPlots");

    // Write pulse height distributions and timing plots for each PMT
    for (int j = 0; j < pulseHeightDistributions.size(); ++j) {
        phdDir->cd();
        pulseHeightDistributions[j]->Write();
        timingDir->cd();
        timingPlots[j]->Write();
    }

    // Create a directory to store summary canvases
    TDirectory* summaryDir = outputFile->mkdir("Summary");

    // Create PDF file to save all summary canvases
    TCanvas* summaryPDFCanvas = new TCanvas("SummaryPDFCanvas", "Summary Canvases in PDF", 800, 600);
    summaryPDFCanvas->Print(Form("%s_SummaryCanvases.pdf[", outputFileName.c_str()));

    // Create a canvas for allhitfreq histogram and save it as the first page
    TCanvas* allHitsCanvas = new TCanvas("AllHitsCanvas", "Combined PMT Hits", 800, 600);
    allHitsCanvas->SetBatch(kTRUE); // Ensure canvas is not opened upon creation
    // Write all histograms
    LUXhitfreq->Write();
    ETELhitfreq->Write();
    HAMAMATSUhitfreq->Write();
    WATCHBOYhitfreq->Write();
    WATCHMANhitfreq->Write();
    RADIUShitfreq->Write();
    THETAhitfreq->Write();
    PHIhitfreq->Write();
    
    

    // Draw all histograms and remove stats box
    allhitfreq->SetStats(0);
    allhitfreq->Draw();
    LUXhitfreq->Draw("SAME");
    ETELhitfreq->Draw("SAME");
    HAMAMATSUhitfreq->Draw("SAME");
    WATCHBOYhitfreq->Draw("SAME");
    WATCHMANhitfreq->Draw("SAME");
    PHIhitfreq->Draw("SAME");
    THETAhitfreq->Draw("SAME");
    RADIUShitfreq->Draw("SAME");

    // Legend of combined histogram
    TLegend* legendp = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendp->AddEntry(LUXhitfreq, "LUX (top)", "f");
    legendp->AddEntry(ETELhitfreq, "ETEL (bottom)", "f");
    legendp->AddEntry(HAMAMATSUhitfreq, "Hamamatsu (bottom)", "f");
    legendp->AddEntry(WATCHBOYhitfreq, "Watchboy (bottom)", "f");
    legendp->AddEntry(WATCHMANhitfreq, "Watchman (tank)", "f");
    legendp->Draw();
    tc->Write();
    
    allHitsCanvas->Print(Form("%s_SummaryCanvases.pdf", outputFileName.c_str()));

    // Loop over all canvases in the "Summary" folder
    for (int j = 0; j < 132; j += 3) {
        TCanvas* summaryCanvas = new TCanvas(Form("SummaryCanvas_%d_to_%d", j, j + 2), "Summary Canvases", 1200, 800);
        summaryCanvas->Divide(3, 2);
        
        // Plot PHDs for 3 PMTs
        for (int k = 0; k < 3; ++k) {
            if (j + k < 132) {
                summaryCanvas->cd(k + 1);
                pulseHeightDistributions[j + k]->Draw();
            }
        }

        // Plot timing plots for the same 3 PMTs
        for (int k = 0; k < 3; ++k) {
            if (j + k < 132) {
                summaryCanvas->cd(k + 4);
                timingPlots[j + k]->Draw();
            }
        }

        // Save canvas to PDF
        summaryPDFCanvas->cd();
        summaryCanvas->Print(Form("%s_SummaryCanvases.pdf", outputFileName.c_str()));

        delete summaryCanvas;
    }

    // Finalize and close the PDF file
    summaryPDFCanvas->Print(Form("%s_SummaryCanvases.pdf]", outputFileName.c_str()));
    delete summaryPDFCanvas;
    delete allHitsCanvas;

    // Close output file and input file
    outputFile->Close();
    file->Close();

    // Delete dynamically allocated memory
    delete outputFile;
    delete file;
}

// Function to return Spherical co-ordinates to vectors
std::vector<double> PMTAnalysis::getSpherical(double x, double y, double z) {
    double r = sqrt(x*x + y*y + z*z);
    double theta = acos(z / r) * 180. / M_PI;
    double phi = atan2(y, x) * 180. / M_PI;
    std::vector<double> spherical = { r, theta, phi };
    return spherical;
}


