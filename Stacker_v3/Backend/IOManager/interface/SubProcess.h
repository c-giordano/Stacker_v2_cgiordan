#ifndef SUBPROCESS_H
#define SUBPROCESS_H

// C++ includes
#include <iostream>
#include <sstream>

// ROOT includes
#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>

// Includes from code

class SubProcess {
    private:
        // Details relevant when loading THs
        TString pName;
        const std::vector<TFile*>& pInputfiles;
        std::vector<std::vector<const char*>*>* pSubdirectoriesPerFile;
        
        // Link to the processinfo -> used for drawing.

        // Private methods
    public:
        SubProcess(TString& name, std::vector<TFile*>& inputfiles);
        ~SubProcess();

        // Getters

        // Setters

        // TH1 extraction
        // Potentially use custom class for this that takes some responsibilities...
        std::shared_ptr<TH1D> ExtractHistogramFromFiles(TString& id, TString& type, TString& subdirStack); // can not work for envelope components
        // either add boolean for break after first loop OR sep fction for envelopes....
        
};


#endif