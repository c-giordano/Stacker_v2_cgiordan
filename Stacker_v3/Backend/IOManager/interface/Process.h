#ifndef PROCESS_H
#define PROCESS_H

// C++ includes
#include <iostream>
#include <sstream>

// ROOT includes
#include <TString.h>
#include <TFile.h>
#include <TH1.h>

// Includes from code
#include "SubProcess.h"
#include "ProcessInfo.h"

#include "../Tools/interface/stringTools.h"

class Process {
    private:
        // Identification of process
        TString pName;

        // Subprocesses collection: responsible for loading elements. Can be only 1 item.
        std::vector<SubProcess*> pSubProcesses;

        // Link to the processinfo -> contains aesthatic information.
        ProcessInfo* pProcessInfo;

        // Private methods
        void mGenerateProcessInfo(std::istringstream& settingsline);
    public:
        Process(std::string& settingsline, std::vector<TFile*>& inputfiles);
        ~Process();

        // Getters
        ProcessInfo* const GetProcessInfo() {return pProcessInfo;};

        // Setters

        // TH1 extraction
        // Most straight up extraction: just sum them up
        std::shared_ptr<TH1D> GetHistogram(TString& id, TString type, TString subdirStack);
        std::shared_ptr<TH1D> GetNominalHistogram(TString& id);
        std::shared_ptr<TH1D> GetUncertaintyHistogram(TString& id, TString& uncID, TString& upOrDown);
        
};


#endif