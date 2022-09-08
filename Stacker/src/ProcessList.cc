#include "../interface/ProcessList.h"
#include <iomanip>

void ProcessList::addProcess(TString& name, int color, TFile* inputfile, TFile* outputfile, bool signal, bool data, bool oldStuff) {
    // TODO: Create new process object
    allProcessNames.push_back(name);
    Process* brandNewObj = new Process(name, color, inputfile, outputfile, signal, data, oldStuff);

    if (tail) {
        brandNewObj->setPrev(tail);
        tail->setNext(brandNewObj); // check if tail already exists
    }
    tail = brandNewObj;

    if (! head) {
        head = brandNewObj;
    }
}

void ProcessList::addProcess(TString& name, int color, std::vector<TFile*>& inputfiles, TFile* outputfile, bool signal, bool data, bool oldStuff) {
    // TODO: Create new process object
    allProcessNames.push_back(name);
    Process* brandNewObj = new Process(name, color, inputfiles, outputfile, signal, data, oldStuff);

    if (tail) {
        brandNewObj->setPrev(tail);
        tail->setNext(brandNewObj); // check if tail already exists
    }
    tail = brandNewObj;

    if (! head) {
        head = brandNewObj;
    }
}

void ProcessList::addProcess(TString& name, std::vector<TString>& procNames, int color, TFile* inputfile, TFile* outputfile, bool signal, bool data, bool oldStuff) {
    // TODO: Create new process object
    allProcessNames.push_back(name);
    Process* brandNewObj = new ProcessSet(name, procNames, color, inputfile, outputfile, signal, data, oldStuff);

    if (tail) {
        brandNewObj->setPrev(tail);
        tail->setNext(brandNewObj); // check if tail already exists
    }
    tail = brandNewObj;

    if (! head) {
        head = brandNewObj;
    }
}

void ProcessList::addProcess(TString& name, std::vector<TString>& procNames, int color, std::vector<TFile*>& inputfiles, TFile* outputfile, bool signal, bool data, bool oldStuff) {
    // TODO: Create new process object
    allProcessNames.push_back(name);
    Process* brandNewObj = new ProcessSet(name, procNames, color, inputfiles, outputfile, signal, data, oldStuff);

    if (tail) {
        brandNewObj->setPrev(tail);
        tail->setNext(brandNewObj); // check if tail already exists
    }
    tail = brandNewObj;

    if (! head) {
        head = brandNewObj;
    }
}

Uncertainty* ProcessList::addUncertainty(std::string& name, bool flat, bool envelope, bool corrProcess, bool eraSpec, std::vector<TString>& processes, TFile* outputfile) {
    Uncertainty* brandNewObj = new Uncertainty(name, flat, envelope, corrProcess, eraSpec, processes, outputfile);

    if (tailUnc) tailUnc->setNext(brandNewObj); // check if tail already exists
    tailUnc = brandNewObj;

    if (! headUnc) {
        headUnc = brandNewObj;
    }

    return brandNewObj;
}


ProcessList::~ProcessList() {
    Process* toDel = head;
    while (toDel->getNext()) {
        Process* nextToDel = toDel->getNext();
        delete toDel;
        toDel = nextToDel;
    }
}

std::vector<std::shared_ptr<TH1D>> ProcessList::fillStack(THStack* stack, Histogram* hist, TLegend* legend, TFile* outfile, std::vector<std::shared_ptr<TH1D>>* signalHistograms, std::shared_ptr<TH1D>* sysUnc) {
    Process* current = head;
    std::vector<std::shared_ptr<TH1D>> histVec;

    TString histogramID = hist->getID();

    double signalYield = 0.;
    double bkgYield = 0.;
    if (hist->getPrintToFile()) outfile->mkdir(hist->getCleanName().c_str());

    if (verbose) std::cout << histogramID << std::endl;

    while (current) {
        std::shared_ptr<TH1D> histToAdd = current->getHistogram(hist);
        if (histToAdd == nullptr) {
            current = current->getNext();
        }
        legend->AddEntry(histToAdd.get(), current->getCleanedName());
        stack->Add(histToAdd.get());
        histVec.push_back(histToAdd);
        
        if (current->isSignalProcess()) {
            signalYield += histToAdd->Integral();
            std::shared_ptr<TH1D> signalHist = std::make_shared<TH1D>(TH1D(*histToAdd));
            signalHistograms->push_back(signalHist);
        } else {
            bkgYield += histToAdd->Integral();
        }

        if (verbose) {
            std::cout << current->getName();
            if (veryVerbose) {
                std::cout << " & ";
                for (int i=1; i < histToAdd->GetNbinsX() + 1; i++) {
                    std::cout << histToAdd->GetBinContent(i) << " & ";
                }
                std::cout << std::endl;
            } else {
                std::cout << ": " << histToAdd->Integral() << " events" << std::endl;
            }
        }

        for (int j=1; j < histToAdd->GetNbinsX() + 1; j++) {
            if (histToAdd->GetBinContent(j) <= 0.) {
                histToAdd->SetBinContent(j, 0.00001);
                histToAdd->SetBinError(j, 0.00001);
            }
        }

        if (hist->getPrintToFile()) {
            outfile->cd(hist->getCleanName().c_str());
            //for (int j=1; j < histToAdd->GetNbinsX() + 1; j++) {
            //    if (histToAdd->GetBinContent(j) <= 0.) {
            //        histToAdd->SetBinContent(j, 0.00001);
            //        histToAdd->SetBinError(j, 0.00001);
            //    }
            //}
            histToAdd->Write(current->getName(), TObject::kOverwrite);
        }
        current = current->getNext();
    }
    
    if (hist->getPrintToFile()) {
        std::shared_ptr<TH1D> allHistograms = sumVector(histVec);
        allHistograms->SetName("data_obs");
        allHistograms->SetTitle("data_obs");
        outfile->cd(hist->getCleanName().c_str());
        for (int j=1; j<allHistograms->GetNbinsX() + 1; j++) {
            allHistograms->SetBinError(j, sqrt(allHistograms->GetBinContent(j)));
        }
        allHistograms->Write("data_obs", TObject::kOverwrite);
    }
    
    // loop uncertainties as well if required
    Uncertainty* currUnc = headUnc;
    std::vector<std::shared_ptr<TH1D>> uncVec;
    while (currUnc && hist->getDrawUncertainties()) {
        // getShapeUncertainty or apply flat uncertainty
        std::shared_ptr<TH1D> newUncertainty = currUnc->getUncertainty(hist, head, histVec);

        uncVec.push_back(newUncertainty);

        if (*sysUnc == nullptr) {
            *sysUnc = std::make_shared<TH1D>( TH1D(*newUncertainty));
        } else {
            (*sysUnc)->Add(newUncertainty.get());
        }

        currUnc = currUnc->getNext();
    }
    
    if (verbose) std::cout << "S/B = " << signalYield << "/" << bkgYield << std::endl;

    if (veryVerbose) {
        std::cout << " & Signal & Bkg & S/B\\\\" << std::endl;
        for (int i=1; i < histVec[0]->GetNbinsX() + 1; i++) {
            double sig = 0.;
            double all = 0.;
            for (unsigned j=0; j < signalHistograms->size(); j++) sig += signalHistograms->at(j)->GetBinContent(i);
            for (unsigned j=0; j < histVec.size(); j++) all += histVec[j]->GetBinContent(i);
            
            std::cout << " & " << std::fixed << std::setprecision(2) << sig << " & " << all - sig << " & " << sig / (all - sig) << "\\\\" << std::endl;
        }
        //std::cout << std::endl;
    }

    return histVec;
}

std::map<TString, bool> ProcessList::printHistograms(Histogram* hist, TFile* outfile, bool isData, Process* dataProc) {
    Process* current = head;
    std::vector<std::shared_ptr<TH1D>> histVec;
    std::map<TString, bool> output;

    TString histogramID = hist->getID();

    if (hist->getPrintToFile()) outfile->mkdir(hist->getCleanName().c_str());

    if (verbose) std::cout << histogramID << std::endl;

    while (current) {

        std::shared_ptr<TH1D> histToAdd = current->getHistogram(hist);
        if (histToAdd == nullptr) {
            current = current->getNext();
        }
        histVec.push_back(histToAdd);

        output[current->getName()] = (histToAdd->Integral() > 0);

        if (verbose) {
            std::cout << current->getName();
            if (veryVerbose) {
                std::cout << " & ";
                for (int i=1; i < histToAdd->GetNbinsX() + 1; i++) {
                    std::cout << histToAdd->GetBinContent(i) << " & ";
                }
                std::cout << std::endl;
            } else {
                std::cout << ": " << histToAdd->Integral() << " events" << std::endl;
            }
        }

        if (hist->getPrintToFile()) {
            outfile->cd(hist->getCleanName().c_str());
            for (int j=1; j < histToAdd->GetNbinsX() + 1; j++) {
                if (histToAdd->GetBinContent(j) <= 0.) {
                    histToAdd->SetBinContent(j, 0.00001);
                    histToAdd->SetBinError(j, 0.00001);
                }
            }
            histToAdd->Write(current->getName(), TObject::kOverwrite);
        }
        current = current->getNext();
    }
    if (! isData && hist->getPrintToFile()) {
        std::shared_ptr<TH1D> allHistograms = sumVector(histVec);
        allHistograms->SetName("data_obs");
        allHistograms->SetTitle("data_obs");
        outfile->cd(hist->getCleanName().c_str());
        for (int j=1; j<allHistograms->GetNbinsX() + 1; j++) {
            allHistograms->SetBinError(j, sqrt(allHistograms->GetBinContent(j)));
        }
        allHistograms->Write("data_obs", TObject::kOverwrite);
    } else if (isData && hist->getPrintToFile()) {
        std::shared_ptr<TH1D> data = dataProc->getHistogram(hist);
        data->SetName("data_obs");
        data->SetTitle("data_obs");
        outfile->cd(hist->getCleanName().c_str());
        data->Write("data_obs", TObject::kOverwrite);
    }

    return output;
}


std::vector<TH2D*> ProcessList::fill2DStack(THStack* stack, TString& histogramID, TLegend* legend, TFile* outfile) {
    Process* current = head;
    std::vector<TH2D*> histVec;

    double signalYield = 0.;
    double bkgYield = 0.;

    //if (hist->getPrintToFile()) outfile->mkdir(histogramID);

    if (verbose) std::cout << histogramID << std::endl;

    while (current) {
        TH2D* histToAdd = current->get2DHistogram(histogramID, legend);
        stack->Add(histToAdd);
        histVec.push_back(histToAdd);
        
        if (current->isSignalProcess()) {
            signalYield += histToAdd->Integral();
        } else {
            bkgYield += histToAdd->Integral();
        }

        if (verbose) {
            std::cout << current->getName() << ": " << histToAdd->Integral() << " events" << std::endl;
        }

        //if (hist->getPrintToFile()) outfile->cd(histogramID);
        //if (hist->getPrintToFile()) histToAdd->Write(current->getName(), TObject::kOverwrite);

        current = current->getNext();
    }
    
    if (verbose) std::cout << "S/B = " << signalYield << "/" << bkgYield << std::endl;

    return histVec;
}

std::map<std::string, std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>>> ProcessList::UpAndDownHistograms(Histogram* hist, std::vector<std::shared_ptr<TH1D>>& nominalHists) {
    // loop uncertainties as well if required
    Uncertainty* currUnc = headUnc;
    std::map<std::string, std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>>> returnValue;
    while (currUnc && hist->getDrawUncertainties()) {
        if (currUnc->isFlat()) {
            currUnc = currUnc->getNext();
            continue;
        }
        // getShapeUncertainty or apply flat uncertainty
        std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> newUncertainty = currUnc->getUpAndDownShapeUncertainty(hist, head, nominalHists);

        returnValue[currUnc->getName()] = newUncertainty;

        currUnc = currUnc->getNext();
    }

    return returnValue;
}

std::vector<std::shared_ptr<TH1D>> ProcessList::CreateHistogramAllProcesses(Histogram* hist) {
    Process* current = getHead();
    std::vector<std::shared_ptr<TH1D>> histVec;

    TString histogramID = hist->getID();

    while (current) {
        std::shared_ptr<TH1D> histToAdd = current->getHistogram(hist);
        if (histToAdd == nullptr) {
            current = current->getNext();
        }
        histVec.push_back(histToAdd);

        current = current->getNext();
    }

    return histVec;
}
