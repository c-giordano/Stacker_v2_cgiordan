#ifndef PROCESS_H
#define PROCESS_H

#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>

#include <iostream>
#include <memory>

#include "../../Helpers/interface/StringTools.h"
#include "Histogram.h"

// Processes as linkedlist: first element of list is pointed to by ProcessList
// Should make it relatively easy to work in a specific order?
// 

class Process {
    private:
        Process* next = nullptr;
        Process* prev = nullptr;

        TString name;
        TString cleanedName;
        std::vector<std::vector<const char*>*> subdirectoriesPerFile;
        Color_t color;
        TFile* rootFile;
        std::vector<TFile*> inputfiles;
        TFile* outputFile;

        bool isSignal;
        bool isData;
    public:
        Process(TString& procName, int procColor, TFile* procInputfile, TFile* outputFile, bool signal, bool data, bool OldStuff);
        Process(TString& procName, int procColor, std::vector<TFile*>& inputfiles, TFile* outputFile, bool signal, bool data, bool OldStuff);
        
        virtual ~Process() {};

        void setNext(Process* newNext) {next = newNext;}
        void setPrev(Process* newPrev) {prev = newPrev;}
        Process* getNext() const {return next;}
        Process* getPrev() const {return prev;}
        TString const getName() {return name;}
        TString const getCleanedName()  {return cleanedName;}
        Color_t getColor() {return color;}
        bool isSignalProcess() {return isSignal;}

        TFile* GetOutputFile() {return outputFile;}

        virtual std::shared_ptr<TH1D> getHistogram(Histogram* histogram);
        virtual std::shared_ptr<TH1D> getHistogramUncertainty(std::string& uncName, std::string& upOrDown, Histogram* hist, std::string& outputFolder, bool envelope);

        virtual std::vector<std::shared_ptr<TH1D>> GetAllVariations(Histogram* histogram, int nVars, std::string& uncName);

        virtual TH2D* get2DHistogram(TString& histName, TLegend* legend);

};

#endif