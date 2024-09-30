/*///////////////////////////////////////////////////////////////
// Tool to be used for analysis of tank PMTs during laser runs //
// For the ANNIE Collaboration                                 //
// By: Eli Brunner-Huber                                       //
// Last revised: 07/18/24                                      //
///////////////////////////////////////////////////////////////*/

#ifndef PMTAnalysis_H
#define PMTAnalysis_H

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <string>

class PMTAnalysis {
public:
    PMTAnalysis();  // Constructor
    virtual ~PMTAnalysis(); // Destructor

    void Analyze(const char* fileName, double startTime, double endTime, double cutoffVoltage);
    std::vector<double> getSpherical(double hitX, double hitY, double hitZ);
    
    void saveHistogramsToPDF(const char* rootFileName, const char* pdfFileName);
    
    void GeneratePulseHeightDistributions();
    
    void GenerateTimingPlots();

private:
    TTree* fChain;
    
    TH1D* allhitfreq;
    TH1D* LUXhitfreq;
    TH1D* ETELhitfreq;
    TH1D* HAMAMATSUhitfreq;
    TH1D* WATCHMANhitfreq;
    TH1D* WATCHBOYhitfreq;
    TH1D* RADIUShitfreq;
    TH1D* THETAhitfreq;
    TH1D* PHIhitfreq;
    std::vector<TH1D*> pulseHeightDistributions;
    std::vector<TH1D*> timingPlots;

    // Private member variables to store parameters
    const char* fFileName;
    double fStartTime;
    double fEndTime;
    double fCutoffVoltage;
};

void saveHistogramsToPDF(const char* rootFileName, const char* pdfFileName);

#endif // PMTAnalysis_H


