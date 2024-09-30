/*///////////////////////////////////////////////////////////////
// Tool to be used for analysis of tank PMTs during laser runs //
// For the ANNIE Collaboration                                 //
// By: Eli Brunner-Huber                                       //
// Last revised: 07/18/24                                      //
///////////////////////////////////////////////////////////////*/



#define pmtLAPPD_cxx
#include "pmtLAPPD.h"
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
pmtLAPPD::pmtLAPPD() : fChain(nullptr), fFileName(""), fStartTime(0), fEndTime(0), fCutOnPE(0), fCutoffPE(0) {
    
    // Initialize histograms in the constructor
    
    WATCHMANpmtsByLAPPD58hitfreq = new TH1D("WATCHMANpmtsByLAPPD58hitfreq", "Frequency of hits on WATCHMAN PMTs by LAPPD 58", 132, 332, 464);
    WATCHMANpmtsByLAPPD58hitfreq->SetFillColor(kOrange);
    WATCHBOYpmtsByLAPPD58hitfreq = new TH1D("WATCHBOYpmtsByLAPPD58hitfreq", "Frequency of hits on WATCHBOY PMTs by LAPPD 58", 132, 332, 464);
    WATCHBOYpmtsByLAPPD58hitfreq->SetFillColor(kYellow);
    HAMAMATSUpmtsByLAPPD58hitfreq = new TH1D("HAMAMATSUpmtsByLAPPD58hitfreq", "Frequency of hits on WATCHBOY PMTs by LAPPD 58", 132, 332, 464);
    HAMAMATSUpmtsByLAPPD58hitfreq->SetFillColor(kGreen);
    
    WATCHMANpmtsByLAPPD39hitfreq = new TH1D("WATCHMANpmtsByLAPPD39hitfreq", "Frequency of hits on WATCHMAN PMTs by LAPPD 39", 132, 332, 464);
    WATCHMANpmtsByLAPPD39hitfreq->SetFillColor(kOrange);
    WATCHBOYpmtsByLAPPD39hitfreq = new TH1D("WATCHBOYpmtsByLAPPD39hitfreq", "Frequency of hits on WATCHBOY PMTs by LAPPD 39", 132, 332, 464);
    WATCHBOYpmtsByLAPPD39hitfreq->SetFillColor(kYellow);
    HAMAMATSUpmtsByLAPPD39hitfreq = new TH1D("HAMAMATSUpmtsByLAPPD39hitfreq", "Frequency of hits on WATCHBOY PMTs by LAPPD 39", 132, 332, 464);
    HAMAMATSUpmtsByLAPPD39hitfreq->SetFillColor(kGreen);
    
    WATCHMANpmtsByLAPPD64hitfreq = new TH1D("WATCHMANpmtsByLAPPD64hitfreq", "Frequency of hits on WATCHMAN PMTs by LAPPD 64", 132, 332, 464);
    WATCHMANpmtsByLAPPD64hitfreq->SetFillColor(kOrange);
    WATCHBOYpmtsByLAPPD64hitfreq = new TH1D("WATCHBOYpmtsByLAPPD64hitfreq", "Frequency of hits on WATCHBOY PMTs by LAPPD 64", 132, 332, 464);
    WATCHBOYpmtsByLAPPD64hitfreq->SetFillColor(kYellow);
    HAMAMATSUpmtsByLAPPD64hitfreq = new TH1D("HAMAMATSUpmtsByLAPPD64hitfreq", "Frequency of hits on WATCHBOY PMTs by LAPPD 64", 132, 332, 464);
    HAMAMATSUpmtsByLAPPD64hitfreq->SetFillColor(kGreen);
    
    /*pmtsByLAPPD64probabilityOfAhit = new TH1D("prob64", "prob64", 132, 332, 464);
    pmtsByLAPPD39probabilityOfAhit = new TH1D("prob39", "prob39", 132, 332, 464);
    pmtsByLAPPD58probabilityOfAhit = new TH1D("prob58", "prob58", 132, 332, 464);*/
    
    
    // Initialize PHD and timing plot histograms for all PMTs in the constructor
    for (int j = 332; j <= 464; ++j) {
        if (j == 377 || j == 407 || j == 374 || j == 445 || j == 462 || j == 428 || j == 406 || j == 412 || j == 411 || j == 400 || j == 404 || j == 464) {
            char name[50];
            sprintf(name, "PulseHeightDistribution_PMT%d", j);
            TH1D* hist = new TH1D(name, name, 10000, 0, 10);
            pulseHeightDistributions.push_back(hist);
        }
    }
    
    for (int j = 332; j <= 464; ++j) {
        if (j == 377 || j == 407 || j == 374 || j == 445 || j == 462 || j == 428 || j == 406 || j == 412 || j == 411 || j == 400 || j == 404 || j == 464) {
            char name[50];
            sprintf(name, "TimingPlot_PMT%d", j);
            TH1D* histo = new TH1D(name, name, 10000, fStartTime, fEndTime);
            timingPlots.push_back(histo);
        }
    }
}

// Destructor
pmtLAPPD::~pmtLAPPD() {
    // Delete histograms in the destructor
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

void pmtLAPPD::Analyze(const char* fileName, double startTime, double endTime, double cutOnPE, double cutoffPE) {
    // Store parameters in variables
    fFileName = fileName;
    fStartTime = startTime;
    fEndTime = endTime;
    fCutOnPE = cutOnPE;
    fCutoffPE = cutoffPE;
    
    
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
    std::vector<int>* hitDetID = nullptr;
    std::vector<double>* hitQ = nullptr;
    std::vector<double>* hitT = nullptr;
    std::vector<double>* hitPE = nullptr;

    
    // Linking TBranches to variables
    fChain->SetBranchAddress("hitDetID", &hitDetID);
    fChain->SetBranchAddress("hitQ", &hitQ);
    fChain->SetBranchAddress("hitT", &hitT);
    fChain->SetBranchAddress("hitPE", &hitPE);
    
    Long64_t nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    
   double totalEntriesHitQ = 0.0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = fChain->LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        
        // Sum up all entries in hitQ
        for (int k = 0; k < hitQ->size(); k++) {
            totalEntriesHitQ += (*hitQ)[k];
        }
    }

    
    // Start loop
    
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = fChain->LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        
        for (int k = 0; k < hitQ->size(); k++) {
            // Set conditions
            if ((*hitPE)[k] < fCutoffPE && (*hitPE)[k] > fCutOnPE && (*hitT)[k] >= fStartTime && (*hitT)[k] <= fEndTime) {
                // Fill PHDs and timing plots
                int pmtID = (*hitDetID)[k];
                if (pmtID == 374) {
                    pulseHeightDistributions[0]->Fill((*hitPE)[k]);
                    timingPlots[0]->Fill((*hitT)[k]);
                }
                
                if (pmtID == 377) {
                    pulseHeightDistributions[1]->Fill((*hitPE)[k]);
                    timingPlots[1]->Fill((*hitT)[k]);
                    
                }
                
                if (pmtID == 400) {
                    pulseHeightDistributions[2]->Fill((*hitPE)[k]);
                    timingPlots[2]->Fill((*hitT)[k]);
                    
                }
                
                if (pmtID == 404) {
                    pulseHeightDistributions[3]->Fill((*hitPE)[k]);
                    timingPlots[3]->Fill((*hitT)[k]);
                    
                }
                
                if (pmtID == 406) {
                    pulseHeightDistributions[4]->Fill((*hitPE)[k]);
                    timingPlots[4]->Fill((*hitT)[k]);
                
                }
                
                if (pmtID == 407) {
                    pulseHeightDistributions[5]->Fill((*hitPE)[k]);
                    timingPlots[5]->Fill((*hitT)[k]);
                
                }
                
                if (pmtID == 406) {
                    pulseHeightDistributions[6]->Fill((*hitPE)[k]);
                    timingPlots[6]->Fill((*hitT)[k]);
                
                }
                
                if (pmtID == 411) {
                    pulseHeightDistributions[7]->Fill((*hitPE)[k]);
                    timingPlots[7]->Fill((*hitT)[k]);
                
                }
                
                if (pmtID == 412) {
                    pulseHeightDistributions[8]->Fill((*hitPE)[k]);
                    timingPlots[8]->Fill((*hitT)[k]);
                
                }
                
                if (pmtID == 428) {
                    pulseHeightDistributions[9]->Fill((*hitPE)[k]);
                    timingPlots[9]->Fill((*hitT)[k]);
                
                }
                
                if (pmtID == 445) {
                    pulseHeightDistributions[10]->Fill((*hitPE)[k]);
                    timingPlots[10]->Fill((*hitT)[k]);
                
                }
                
                if (pmtID == 462) {
                    pulseHeightDistributions[11]->Fill((*hitPE)[k]);
                    timingPlots[11]->Fill((*hitT)[k]);
                
                }
                
                if (pmtID == 464) {
                    pulseHeightDistributions[12]->Fill((*hitPE)[k]);
                    timingPlots[12]->Fill((*hitT)[k]);
                
                }
                
                if (pmtID == 377 || pmtID == 407 || pmtID == 374) {
                    HAMAMATSUpmtsByLAPPD64hitfreq->Fill(pmtID);
                }
                
                if (pmtID == 445) {
                    WATCHBOYpmtsByLAPPD64hitfreq->Fill(pmtID);
                }
                
                if (pmtID == 462 || pmtID == 428) {
                    WATCHBOYpmtsByLAPPD39hitfreq->Fill(pmtID);
                    
                }
                
                if (pmtID == 406 || pmtID == 412) {
                    HAMAMATSUpmtsByLAPPD39hitfreq->Fill(pmtID);
                }
                
                if (pmtID == 411 || pmtID == 400) {
                    HAMAMATSUpmtsByLAPPD58hitfreq->Fill(pmtID);
                    
                }
                
                if (pmtID == 404) {
                    WATCHMANpmtsByLAPPD58hitfreq->Fill(pmtID);
                }
                
               if (pmtID == 464) {
                    WATCHBOYpmtsByLAPPD58hitfreq->Fill(pmtID);
                }
                
                
            }
        }
    }
    
    // Construct output file name with input file name included
    std::string outputFileName = "pmtLAPPDAnalysis_" + std::string(fileName);
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
    /*pmtsByLAPPD64probabilityOfAhit->Write();
    pmtsByLAPPD39probabilityOfAhit->Write();
    pmtsByLAPPD58probabilityOfAhit->Write();*/
    
    
    // Normalize histograms
    // Assuming you have already filled the histograms, now normalize them
    double scaleFactor = 1120000;  // You may want to adjust this scale factor based on your needs
    if (scaleFactor > 0) {
        HAMAMATSUpmtsByLAPPD64hitfreq->Scale(1.0 / scaleFactor);
        WATCHBOYpmtsByLAPPD64hitfreq->Scale(1.0 / scaleFactor);
        WATCHMANpmtsByLAPPD64hitfreq->Scale(1.0 / scaleFactor);
        
        HAMAMATSUpmtsByLAPPD39hitfreq->Scale(1.0 / scaleFactor);
        WATCHBOYpmtsByLAPPD39hitfreq->Scale(1.0 / scaleFactor);
        WATCHMANpmtsByLAPPD39hitfreq->Scale(1.0 / scaleFactor);
        
        HAMAMATSUpmtsByLAPPD58hitfreq->Scale(1.0 / scaleFactor);
        WATCHBOYpmtsByLAPPD58hitfreq->Scale(1.0 / scaleFactor);
        WATCHMANpmtsByLAPPD58hitfreq->Scale(1.0 / scaleFactor);
        
      /*  pmtsByLAPPD64probabilityOfAhit->Scale(1.0 / scaleFactor);
        pmtsByLAPPD58probabilityOfAhit->Scale(1.0 / scaleFactor);
        pmtsByLAPPD39probabilityOfAhit->Scale(1.0 / scaleFactor);*/
        
    }

    
    TCanvas* tc64 = new TCanvas("PMTs by LAPPD 64 hit frequencies", "PMTs by LAPPD 64 hit frequencies");
    
    HAMAMATSUpmtsByLAPPD64hitfreq->SetStats(0);
    HAMAMATSUpmtsByLAPPD64hitfreq->SetTitle("Frequency of hits on PMTs by LAPPD 64");
    HAMAMATSUpmtsByLAPPD64hitfreq->Draw("HIST");
    WATCHBOYpmtsByLAPPD64hitfreq->Draw("HIST SAME");
    WATCHMANpmtsByLAPPD64hitfreq->Draw("HIST SAME");
    
    // Legend of combined histogram
    TLegend* legend64 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend64->AddEntry(HAMAMATSUpmtsByLAPPD64hitfreq, "Hamamatsu (tank)", "f");
    legend64->AddEntry(WATCHBOYpmtsByLAPPD64hitfreq, "Watchboy (tank)", "f");
    legend64->AddEntry(WATCHMANpmtsByLAPPD64hitfreq, "Watchman (tank)", "f");
    legend64->Draw();
    tc64->Write();
    
    TCanvas* tc39 = new TCanvas("PMTs by LAPPD 39 hit frequencies", "PMTs by LAPPD 39 hit frequencies");
    
    HAMAMATSUpmtsByLAPPD39hitfreq->SetStats(0);
    HAMAMATSUpmtsByLAPPD39hitfreq->SetTitle("Frequency of hits on PMTs by LAPPD 39");
    HAMAMATSUpmtsByLAPPD39hitfreq->Draw("HIST");
    WATCHBOYpmtsByLAPPD39hitfreq->Draw("HIST SAME");
    WATCHMANpmtsByLAPPD39hitfreq->Draw("HIST SAME");
    
    // Legend of combined histogram
    TLegend* legend39 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend39->AddEntry(HAMAMATSUpmtsByLAPPD39hitfreq, "Hamamatsu (tank)", "f");
    legend39->AddEntry(WATCHBOYpmtsByLAPPD39hitfreq, "Watchboy (tank)", "f");
    legend39->AddEntry(WATCHMANpmtsByLAPPD39hitfreq, "Watchman (tank)", "f");
    legend39->Draw();
    tc39->Write();
    
    TCanvas* tc58 = new TCanvas("PMTs by LAPPD 58 hit frequencies", "PMTs by LAPPD 58 hit frequencies");
    
    HAMAMATSUpmtsByLAPPD58hitfreq->SetStats(0);
    HAMAMATSUpmtsByLAPPD58hitfreq->SetTitle("Frequency of hits on PMTs by LAPPD 58");
    HAMAMATSUpmtsByLAPPD58hitfreq->Draw("HIST");
    WATCHBOYpmtsByLAPPD58hitfreq->Draw("HIST SAME");
    WATCHMANpmtsByLAPPD58hitfreq->Draw("HIST SAME");
    
    // Legend of combined histogram
    TLegend* legend58 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend58->AddEntry(HAMAMATSUpmtsByLAPPD58hitfreq, "Hamamatsu (tank)", "f");
    legend58->AddEntry(WATCHBOYpmtsByLAPPD58hitfreq, "Watchboy (tank)", "f");
    legend58->AddEntry(WATCHMANpmtsByLAPPD58hitfreq, "Watchman (tank)", "f");
    legend58->Draw();
    tc58->Write();
    

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

    
    // Close output file and input file
    outputFile->Close();
    file->Close();

    // Delete dynamically allocated memory
    delete outputFile;
    delete file;
}

