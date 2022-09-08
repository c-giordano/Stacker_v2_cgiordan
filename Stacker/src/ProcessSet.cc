#include "../interface/ProcessSet.h"


ProcessSet::ProcessSet(TString& name, std::vector<TString>& procNames, int procColor, TFile* procInputfile, TFile* outputFile, bool signal, bool data, bool OldStuff) :
    Process(name, procColor, procInputfile, outputFile, signal, data, OldStuff)
    {
    for (auto it : procNames) {
        Process* newProc = new Process(it, procColor, procInputfile, outputFile, signal, data, OldStuff);
        subProcesses.push_back(newProc);
    }
}

ProcessSet::ProcessSet(TString& name, std::vector<TString>& procNames, int procColor, std::vector<TFile*>& inputfiles, TFile* outputFile, bool signal, bool data, bool OldStuff) :
    Process(name, procColor, inputfiles, outputFile, signal, data, OldStuff)
    {
    for (auto it : procNames) {
        Process* newProc = new Process(it, procColor, inputfiles, outputFile, signal, data, OldStuff);
        subProcesses.push_back(newProc);
    }
}


std::shared_ptr<TH1D> ProcessSet::getHistogram(Histogram* histogram) {
    std::shared_ptr<TH1D> output = nullptr;
    TString histName = histogram->getID();

    for (auto it : subProcesses) {
        std::shared_ptr<TH1D> tmp = it->getHistogram(histogram);
        if (tmp == nullptr) {
            continue;
        }

        if (output == nullptr) {
            output = tmp;
        } else {
            output->Add(tmp.get());
        }
    }

    output->SetName(histName + getName());
    output->SetTitle(histName + getName());
    
    output->SetLineColor(getColor());
    output->SetFillColor(getColor());
    output->SetMarkerColor(getColor());

    return output;
}

std::shared_ptr<TH1D> ProcessSet::getHistogramUncertainty(std::string& uncName, std::string& upOrDown, Histogram* hist, std::string& outputFolder, bool envelope) {
    std::shared_ptr<TH1D> output = nullptr;
   // bool printToFile = hist->getPrintToFile();
   // hist->setPrintToFile(false);

    for (auto it : subProcesses) {
        std::shared_ptr<TH1D> tmp = it->getHistogramUncertainty(uncName, upOrDown, hist, outputFolder, envelope);
        
        if (output == nullptr) {
            output = tmp;
        } else {
            output->Add(tmp.get());
        }
    }

    output->SetName(hist->getID() + getName() + TString(uncName + upOrDown));
    output->SetTitle(hist->getID() + getName() + TString(uncName + upOrDown));

    //if (printToFile) {
    //    GetOutputFile()->cd();
    //    GetOutputFile()->cd((hist->getCleanName()).c_str());
//
    //    gDirectory->cd((outputFolder + upOrDown).c_str());
//
    //    output->Write(getName());
    //    hist->setPrintToFile(printToFile);
    //}

    return output;
}

TH2D* ProcessSet::get2DHistogram(TString& histName, TLegend* legend) {
    TH2D* output = nullptr;

    for (auto it : subProcesses) {
        TH2D* tmp = it->get2DHistogram(histName, legend);

        if (output == nullptr) {
            output = tmp;
        } else {
            output->Add(tmp);
        }
    }

    output->SetName(histName + getName());
    output->SetTitle(histName + getName());
    
    output->SetLineColor(getColor());
    output->SetFillColor(getColor());
    output->SetMarkerColor(getColor());

    return output;
}

