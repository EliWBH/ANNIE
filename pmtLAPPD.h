/*///////////////////////////////////////////////////////////////
// Tool to be used for analysis of tank PMTs during laser runs //
// For the ANNIE Collaboration                                 //
// By: Eli Brunner-Huber                                       //
// Last revised: 07/18/24                                      //
///////////////////////////////////////////////////////////////*/

#ifndef pmtLAPPD_H
#define pmtLAPPD_H

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <string>

class pmtLAPPD {
public:
    pmtLAPPD();  // Constructor
    virtual ~pmtLAPPD(); // Destructor

    void Analyze(const char* fileName, double startTime, double endTime, double cutOnPE, double cutoffPE);
    
    void GeneratePulseHeightDistributions();
    
    void GenerateTimingPlots();

private:
    TTree* fChain;
    
    TH1D* HAMAMATSUpmtsByLAPPD64hitfreq;
    TH1D* WATCHBOYpmtsByLAPPD64hitfreq;
    TH1D* WATCHMANpmtsByLAPPD64hitfreq;
    TH1D* HAMAMATSUpmtsByLAPPD58hitfreq;
    TH1D* WATCHBOYpmtsByLAPPD58hitfreq;
    TH1D* WATCHMANpmtsByLAPPD58hitfreq;
    TH1D* HAMAMATSUpmtsByLAPPD39hitfreq;
    TH1D* WATCHBOYpmtsByLAPPD39hitfreq;
    TH1D* WATCHMANpmtsByLAPPD39hitfreq;
    TH1D* pmtsByLAPPD64probabilityOfAhit;
    TH1D* pmtsByLAPPD39probabilityOfAhit;
    TH1D* pmtsByLAPPD58probabilityOfAhit;
    std::vector<TH1D*> pulseHeightDistributions;
    std::vector<TH1D*> timingPlots;

    // Private member variables to store parameters
    const char* fFileName;
    double fStartTime;
    double fEndTime;
    double fCutoffPE;
    double fCutOnPE;
};

void saveHistogramsToPDF(const char* rootFileName, const char* pdfFileName);

#endif // pmtLAPPD_H


