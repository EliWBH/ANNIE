#ifndef TDCPROCESSOR_H
#define TDCPROCESSOR_H

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <algorithm>

class TDCProcessor {
public:
    TDCProcessor();
    ~TDCProcessor();

    void SetFileName(const char* filename);
    void SetWindowSize(double window_size);
    void Process();

private:
    const char* fFileName; // File name
    double fWindowSize;    // Time window size

    void ProcessTDC();
};

// Global function to use in ROOT
extern "C" void process_tdc(const char* filename, double window_size);


#endif // TDCPROCESSOR_H
