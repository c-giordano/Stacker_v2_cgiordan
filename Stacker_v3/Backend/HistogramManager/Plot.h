#ifndef PLOT_H
#define PLOT_H

// C++ Includes
#include <vector>
#include <memory>

// ROOT Includes
#include <TH1.h>
#include <THStack.h>

// Stacker includes
#include "UncertaintyTHs.h"

/*
Base class managing a single plot. For each plot that is made, this object is created containing the relevant info.
It contains all special settings for a plot, manages binning, ...
*/

class Plot
{
    private:
        // Object IDs
        TString pID;

        // Settings
        bool pLogscale;

        bool pRebin;
        int pRebinNb;
        double* pRebinEdges;
        std::vector<std::string> *pXBinLabels;

        // Plot memory: loaded objects
        std::vector<std::shared_ptr<TH1D>>* pProcessTHs;
        std::vector<UncertaintyTHs*>* pUncertaintyTHs;

    public:
        Plot(/* args */);
        ~Plot();

        // Getters
        TString const GetID() { return pID; };
        bool const GetLogscale() { return pLogscale; };
        bool const GetRebin() { return pRebin; };
        int const GetRebinNb() { return pRebinNb; };
        double* const GetRebinEdges() { return pRebinEdges; };
        std::vector<std::string> *const GetXBinLabels() { return pXBinLabels; };

        // Setters
        void SetID(TString id) { pID = id; };
        void SetLogscale(bool logscale) { pLogscale = logscale; };
        void SetRebin(bool rebin) { pRebin = rebin; };
        void SetRebinNb(int rebinNb) { pRebinNb = rebinNb; };
        void SetRebinEdges(double* rebinEdges) { pRebinEdges = rebinEdges; };
        void SetXBinLabels(std::vector<std::string> xBinLabels) { pXBinLabels = new std::vector<std::string>(xBinLabels); };

        // Able to combine information in an intelligent way: building stacks, sums, ... -> idea to extract any info desired by other classes from Plot class -> sort of an overhead of a lot of backend and no need to access process etc
        THStack* BuildStack();
        THStack* BuildStack();

        // TBD: End product production -> i.e. DrawFigure prints a png/pdf/whatever
        void DrawFigure();

    


};

#endif
