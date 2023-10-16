#include "../interface/Process.h"

Process::Process(std::string& settingsline, std::vector<TFile*>& inputfiles) { 
    std::istringstream stream(settingsline);
    
    std::string processName;
    stream >> processName;
    pName = TString(processName);
    mGenerateProcessInfo(stream);

    std::vector<TString> components;

    // check for special stuff
    std::string part;
    if (stream >> part) {
        // check if it contains a plus. If so, the current process must be modified to take into account multiple subdirectories
        // maybe make a process object containing multiple processes
        std::vector<std::string> processNamesStrings = stringTools::split(part, "+");
        for (auto it : processNamesStrings) {
            components.push_back(TString(it));
        }
    } else {
        components.push_back(pName);
    }

    for (auto procName : components) {
        SubProcess* add = new SubProcess(procName, inputfiles);
        pSubProcesses.push_back(add);
    }

}

Process::~Process() {
    delete pProcessInfo;

    for (auto ptr : pSubProcesses) {
        delete ptr;
    }
}

void Process::mGenerateProcessInfo(std::istringstream& settingsline) {
    std::string color;
    std::string type;

    settingsline >> color >> type;

    pProcessInfo = new ProcessInfo(pName, color, type);
}


// Get Histogram function
// Loops all input files and extracts all relevant information.
// Adding them together or not:
// For nominal: never a problem
// For uncertainties: yes also fine for normal stuff. Envelopes: already per process
// Generic get function that allows to add subdirectories  etc
// return std::shared_ptr<TH1D>

std::shared_ptr<TH1D> Process::GetNominalHistogram(TString& id) {
    std::shared_ptr<TH1D> ret = GetHistogram(id, "Nominal", "");

    ret->SetLineColor(pProcessInfo->GetColor());
    ret->SetFillColor(pProcessInfo->GetColor());
    ret->SetMarkerColor(pProcessInfo->GetColor());
    
    return ret;
}

std::shared_ptr<TH1D> Process::GetUncertaintyHistogram(TString& id, TString& uncID, TString& upOrDown) {
    TString subdirs = uncID + "/" + upOrDown;
    std::shared_ptr<TH1D> ret = GetHistogram(id, "Uncertainties", subdirs);

    return ret;
}


std::shared_ptr<TH1D> Process::GetHistogram(TString& id, TString type, TString subdirStack) {
    std::vector<std::shared_ptr<TH1D>> subProcessOutput;
    std::shared_ptr<TH1D> ret = nullptr; 

    for (auto subProcess : pSubProcesses) {
        std::shared_ptr<TH1D> histogram = subProcess->ExtractHistogramFromFiles(id, type, subdirStack);

        if (ret == nullptr) {
            ret = std::make_shared<TH1D>(*histogram.get());
        } else {
            ret->Add(histogram.get());
        }
    }

    ret->SetName(id + pName + subdirStack);
    ret->SetTitle(id + pName + subdirStack);
    
    return ret;
}