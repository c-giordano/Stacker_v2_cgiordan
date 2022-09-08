#include "../interface/Uncertainty.h"

Uncertainty::Uncertainty(std::string& name, bool flat, bool envelope, bool corrProcess, bool eraSpec, std::vector<TString>& processes, TFile* outputfile) : 
    name(name), flat(flat), envelope(envelope), correlatedAmongProcesses(corrProcess), eraSpecific(eraSpec), relevantProcesses(processes), outfile(outputfile) {
    outputName = name;
    nameUp = name + "Up";
    nameDown = name + "Down";
}


std::shared_ptr<TH1D> Uncertainty::getUncertainty(Histogram* histogram, Process* head, std::vector<std::shared_ptr<TH1D>>& histVec) {
    if (flat) {
        return getFlatUncertainty(histogram, head, histVec);
    } else if (envelope && buildEnvelope) {
        return getEnvelope(histogram, head, histVec);
    } else {
        return getShapeUncertainty(histogram, head, histVec);
    }
}


std::shared_ptr<TH1D> Uncertainty::getShapeUncertainty(Histogram* histogram, Process* head, std::vector<std::shared_ptr<TH1D>>& histVec) {
    // Loop processes, ask to add stuff

    // create one final TH1D to return
    Process* current = head;
    std::shared_ptr<TH1D> ret = std::make_shared<TH1D>(TH1D(*histVec[0]));

    for (int bin = 1; bin < ret->GetNbinsX() + 1; bin++) {
        ret->SetBinContent(bin, 0.);
    }

    TString histogramID = histogram->getID();
    ret->SetName(histogramID + name);
    ret->SetTitle(name.c_str());

    //correlated case : linearly add up and down variations
    std::vector<double> varDown(histVec[0]->GetNbinsX(), 0.);
    std::vector<double> varUp(histVec[0]->GetNbinsX(), 0.);

    //uncorrelated case : quadratically add the maximum of the up and down variations
    std::vector<double> var(histVec[0]->GetNbinsX(), 0.);

    int histCount = 0;
    int procCount = 0;

    if (histogram->getPrintToFile()) {
        outfile->cd(histogram->getCleanName().c_str());
        gDirectory->mkdir(nameUp);
        gDirectory->mkdir(nameDown);
    }

    std::string up = "Up";
    std::string down = "Down";

    //std::cout << outputName << std::endl;
    while (current) {
        // TODO: check if uncertainty needs this process, otherwise continue and put stuff to next one;
        if (current->getName() != relevantProcesses[procCount]) {
            histCount++;
            current = current->getNext();

            continue;
        }
        //std::cout << current->getName();

        std::shared_ptr<TH1D> histNominal = histVec[histCount];
        std::shared_ptr<TH1D> histUp = current->getHistogramUncertainty(name, up, histogram, outputName, isEnvelope());
        std::shared_ptr<TH1D> histDown = current->getHistogramUncertainty(name, down, histogram, outputName, isEnvelope());

        if (histUp == nullptr && histDown == nullptr) {
            current = current->getNext();
            histCount++;
            procCount++;
            continue;
        }

        if (histogram->HasRebin()) {
            if (histogram->GetRebinVar() == nullptr) {
                histUp->Rebin(histogram->GetRebin());
                histDown->Rebin(histogram->GetRebin());
            } else {
                TString newNameUp = histUp->GetName();
                histUp->SetName(newNameUp + TString("OLD"));
                TString newNameDown = histDown->GetName();
                histDown->SetName(newNameDown + TString("OLD"));
                
                histUp = std::make_shared<TH1D> (*(TH1D*)histUp->Rebin(histogram->GetRebin()-1, newNameUp, histogram->GetRebinVar()));
                histDown = std::make_shared<TH1D> (*(TH1D*)histDown->Rebin(histogram->GetRebin()-1, newNameDown, histogram->GetRebinVar()));
            }
        }
        // do stuff
        // anyway
        for (int bin = 1; bin < histNominal->GetNbinsX() + 1; bin++) {

            double nominalContent = histNominal->GetBinContent( bin );
            double upVariedContent = histUp->GetBinContent( bin );
            double downVariedContent = histDown->GetBinContent( bin );
            double down = fabs( downVariedContent - nominalContent );
            double up = fabs( upVariedContent - nominalContent );

            //uncorrelated case : 
            if(! correlatedAmongProcesses ){
                double variation = std::max( down, up );
                var[bin - 1] += variation*variation;
            
            //correlated case :     
            } else {
                varDown[bin - 1] += down;
                varUp[bin - 1] += up;
            }

        }

        current = current->getNext();
        histCount++;
        procCount++;
    } 

    //std::cout << "bin content in Uncertainty:\t";

    for (int bin = 1; bin < histVec[0]->GetNbinsX() + 1; bin++) {
        //correlated case :
        if(correlatedAmongProcesses ){
            double varLocal = std::max( varDown[bin - 1], varUp[bin - 1] );
            var[bin - 1] = varLocal * varLocal;
        }

        ret->SetBinContent(bin, var[bin - 1]);
        //std::cout << ret->GetBinContent(bin) << "\t";
    }
    //std::cout << std::endl;
    // writeout uncertainty
    return ret;
    
    // set return value
}

std::shared_ptr<TH1D> Uncertainty::getEnvelope(Histogram* histogram, Process* head, std::vector<std::shared_ptr<TH1D>>& histVec) {
    // combine stuff of different variations here
    // use normal envelopes

    // create one final TH1D to return
    Process* current = head;
    std::shared_ptr<TH1D> ret = std::make_shared<TH1D>(TH1D(*histVec[0]));


    for (int bin = 1; bin < ret->GetNbinsX() + 1; bin++) {
        ret->SetBinContent(bin, 0.);
    }

    TString histogramID = histogram->getID();
    ret->SetName(histogramID + name);
    ret->SetTitle(name.c_str());

    //correlated case : linearly add up and down variations
    std::vector<double> varDown(histVec[0]->GetNbinsX(), 0.);
    std::vector<double> varUp(histVec[0]->GetNbinsX(), 0.);

    //uncorrelated case : quadratically add the maximum of the up and down variations
    std::vector<double> var(histVec[0]->GetNbinsX(), 0.);

    int histCount = 0;
    int procCount = 0;

    if (histogram->getPrintToFile()) {
        outfile->cd(histogram->getCleanName().c_str());
        gDirectory->mkdir(nameUp);
        gDirectory->mkdir(nameDown);
    }

    std::string up = "Up";
    std::string down = "Down";


    while (current) {
        // TODO: check if uncertainty needs this process, otherwise continue and put stuff to next one;
        if (current->getName() != relevantProcesses[procCount]) {
            histCount++;
            current = current->getNext();
        }

        std::shared_ptr<TH1D> histNominal = histVec[histCount];
        std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> histVars;

        if (name == "pdfShapeVar") {
            // alternative code is still available for building correct envelope here. Should in principle never be used anymore
            histVars = buildPDFFromSumSquaredCollections(histogram, current, histNominal);
        } else {
            histVars = buildEnvelopeForProcess(histogram, current, histNominal);
        }

        std::shared_ptr<TH1D> histUp = histVars.first;
        std::shared_ptr<TH1D> histDown = histVars.second;

        if (histUp == nullptr && histDown == nullptr) {
            current = current->getNext();
            histCount++;
            procCount++;
            continue;
        }


        if (histogram->getPrintToFile()) {
            for (int j=1; j < histUp->GetNbinsX() + 1; j++) {
                if (histUp->GetBinContent(j) <= 0.) histUp->SetBinContent(j, 0.00001);
                if (histDown->GetBinContent(j) <= 0.) histDown->SetBinContent(j, 0.00001);
            }

            TFile* outputFile = current->GetOutputFile();
            outputFile->cd();
            outputFile->cd((histogram->getCleanName()).c_str());

            gDirectory->cd((outputName + "Up").c_str());
            histUp->Write(current->getName());

            gDirectory->cd((outputName + "Down").c_str());
            histDown->Write(current->getName());
        }

        // do stuff
        // anyway
        for (int bin = 1; bin < histNominal->GetNbinsX() + 1; bin++) {

            double nominalContent = histNominal->GetBinContent( bin );
            double upVariedContent = histUp->GetBinContent( bin );
            double downVariedContent = histDown->GetBinContent( bin );
            double down = fabs( downVariedContent - nominalContent );
            double up = fabs( upVariedContent - nominalContent );

            //uncorrelated case : 
            if(! correlatedAmongProcesses ){
                double variation = std::max( down, up );
                var[bin - 1] += variation*variation;
            
            //correlated case :     
            } else {
                varDown[bin - 1] += down;
                varUp[bin - 1] += up;
            }

        }

        current = current->getNext();
        histCount++;
        procCount++;
    } 

    //std::cout << "bin content in Uncertainty:\t";

    for (int bin = 1; bin < histVec[0]->GetNbinsX() + 1; bin++) {
        //correlated case :
        if(correlatedAmongProcesses ){
            double varLocal = std::max( varDown[bin - 1], varUp[bin - 1] );
            var[bin - 1] = varLocal * varLocal;
        }

        ret->SetBinContent(bin, var[bin - 1]);
        //std::cout << ret->GetBinContent(bin) << "\t";
    }

    return ret;
    
    // set return value

}

std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> Uncertainty::buildEnvelopeForProcess(Histogram* histogram, Process* currentProcess, std::shared_ptr<TH1D> nominalHist) {
    // call this for each histogram (comes from earlier call. So only given histogram is relevant)
    // loop per process
    // get all n variations of this histogram for a given process, create each bincontent. 
    // postprocess for given process
    // should only be qcd uncertainty 
    std::shared_ptr<TH1D> upVar = std::make_shared<TH1D>(TH1D(*nominalHist));
    std::shared_ptr<TH1D> downVar = std::make_shared<TH1D>(TH1D(*nominalHist));

    upVar->SetName(histogram->getID() + currentProcess->getName() + TString(name + "Up"));
    upVar->SetTitle(histogram->getID() + currentProcess->getName() + TString(name + "Up"));

    downVar->SetName(histogram->getID() + currentProcess->getName() + TString(name + "Down"));
    downVar->SetTitle(histogram->getID() + currentProcess->getName() + TString(name + "Down"));

    unsigned variations = 6;

    std::vector<std::shared_ptr<TH1D>> allHistogramVariations = currentProcess->GetAllVariations(histogram, variations, name);

    for (auto currentVariation : allHistogramVariations) {
        // loop all possible up and down variations
        // get histograms for current variation
        // combine in correct way

        if (histogram->HasRebin()) {
            if (histogram->GetRebinVar() == nullptr) {
                currentVariation->Rebin(histogram->GetRebin());
            } else {
                std::string newName = currentVariation->GetName();
                currentVariation->SetName((newName + "OLD").c_str());
                currentVariation = std::make_shared<TH1D>(*(TH1D*) currentVariation->Rebin(histogram->GetRebin()-1, newName.c_str(), histogram->GetRebinVar()));
            }
        }

        for(int j=1; j < currentVariation->GetNbinsX()+1; j++){

            // for each up and down variation, we fix the content
            double bincontentCurrent = currentVariation->GetBinContent(j);

            if (bincontentCurrent > upVar->GetBinContent(j) ){
                upVar->SetBinContent(j, bincontentCurrent);
            } else if (bincontentCurrent < downVar->GetBinContent(j) ){
                downVar->SetBinContent(j, bincontentCurrent);
            }
        }
    }

    return {upVar, downVar};
}

std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> Uncertainty::buildPDFFromSumSquaredCollections(Histogram* histogram, Process* currentProcess, std::shared_ptr<TH1D> nominalHist) {
    // call this for each histogram (comes from earlier call. So only given histogram is relevant)
    // loop per process
    // get all n variations of this histogram for a given process, create each bincontent. 
    // postprocess for given process

    std::shared_ptr<TH1D> upVar = std::make_shared<TH1D>(TH1D(*nominalHist));
    upVar->SetName(histogram->getID() + currentProcess->getName() + TString(name + "Up"));
    upVar->SetTitle(histogram->getID() + currentProcess->getName() + TString(name + "Up"));

    std::shared_ptr<TH1D> downVar = std::make_shared<TH1D>(TH1D(*nominalHist));
    downVar->SetName(histogram->getID() + currentProcess->getName() + TString(name + "Down"));
    downVar->SetTitle(histogram->getID() + currentProcess->getName() + TString(name + "Down"));

    std::string up = "Up";
    std::shared_ptr<TH1D> totalVar = currentProcess->getHistogramUncertainty(name, up, histogram, outputName, isEnvelope());

    // get all histograms in up variation -> Should automatically sum in process

    for(int j=1; j < totalVar->GetNbinsX()+1; j++){
        // for each up and down variation, we fix the content
        totalVar->SetBinContent(j, sqrt(totalVar->GetBinContent(j)));
    }

    if (histogram->HasRebin()) {
        if (histogram->GetRebinVar() == nullptr) {
            totalVar->Rebin(histogram->GetRebin());
        } else {
            std::string newName = std::string(histogram->getID() + name + "Complete");
            totalVar = std::make_shared<TH1D> (*(TH1D*) totalVar->Rebin(histogram->GetRebin()-1, newName.c_str(), histogram->GetRebinVar()));
        }
    }
    
    upVar->Add(totalVar.get());
    downVar->Add(totalVar.get(), -1.);

    return {upVar, downVar};
}

std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> Uncertainty::buildSumSquaredEnvelopeForProcess(Histogram* histogram, Process* currentProcess, std::shared_ptr<TH1D> nominalHist) {
    // call this for each histogram (comes from earlier call. So only given histogram is relevant)
    // loop per process
    // get all n variations of this histogram for a given process, create each bincontent. 
    // postprocess for given process
    std::shared_ptr<TH1D> totalVar = std::make_shared<TH1D>(TH1D(*nominalHist));
    std::shared_ptr<TH1D> upVar = std::make_shared<TH1D>(TH1D(*nominalHist));
    std::shared_ptr<TH1D> downVar = std::make_shared<TH1D>(TH1D(*nominalHist));

    // make sure they are empty
    totalVar->Reset();

    unsigned variations = 100;

    std::vector<std::shared_ptr<TH1D>> allHistogramVariations = currentProcess->GetAllVariations(histogram, variations, name);

    for (auto currentVariation : allHistogramVariations) {
        // loop all possible up and down variations
        // get histograms for current variation
        // combine in correct way

        // for each contribution, take diff with nominal and square
        for(int j=1; j < currentVariation->GetNbinsX()+1; j++){
            // for each up and down variation, we fix the content
            double bincontentCurrent = currentVariation->GetBinContent(j) - nominalHist->GetBinContent(j);
            totalVar->SetBinContent(j, totalVar->GetBinContent(j) + bincontentCurrent * bincontentCurrent);
        }
    }

    for(int j=1; j < totalVar->GetNbinsX()+1; j++){
        // for each up and down variation, we fix the content
        totalVar->SetBinContent(j, sqrt(totalVar->GetBinContent(j)));
    }
    
    upVar->Add(totalVar.get());
    downVar->Add(totalVar.get(), -1.);

    return {upVar, downVar};
}

void Uncertainty::printOutShapeUncertainty(Histogram* histogram, Process* head) {
    // Loop processes, ask to add stuff
    Process* current = head;
    TString histogramID = histogram->getID();

    int histCount = 0;
    int procCount = 0;

    TString outputNameUp = outputName + "Up";
    TString outputNameDown = outputName + "Down";

    if (histogram->getPrintToFile()) {
        outfile->cd(histogram->getCleanName().c_str());
        gDirectory->mkdir(outputNameUp);
        gDirectory->mkdir(outputNameDown);
    }

    std::string up = "Up";
    std::string down = "Down";

    while (current) {
        // TODO: check if uncertainty needs this process, otherwise continue and put stuff to next one;
        if (relevantProcesses.size() == 0) return;
        if (current->getName() != relevantProcesses[procCount]) {
            histCount++;
            current = current->getNext();
        }

        //std::shared_ptr<TH1D> upVar = nullptr;
        //std::shared_ptr<TH1D> downVar = nullptr;
        //std::shared_ptr<TH1D> histNominal = nominalHists[histCount];
        ////current->getHistogramUncertainty(name, up, histogram, outputName, isEnvelope());
        ////current->getHistogramUncertainty(name, down, histogram, outputName, isEnvelope());
//
        //if (! envelope || (envelope && ! buildEnvelope)) {
        //    upVar = current->getHistogramUncertainty(name, up, histogram, outputName, isEnvelope());
        //    downVar = current->getHistogramUncertainty(name, down, histogram, outputName, isEnvelope());
//
        //    if (histogram->HasRebin()) {
        //        if (histogram->GetRebinVar() == nullptr) {
        //            upVar->Rebin(histogram->GetRebin());
        //            downVar->Rebin(histogram->GetRebin());
        //        } else {
        //            TString newNameUp = upVar->GetName();
        //            upVar->SetName(newNameUp + TString("OLD"));
        //            TString newNameDown = downVar->GetName();
        //            downVar->SetName(newNameDown + TString("OLD"));
        //            
        //            upVar = (std::shared_ptr<TH1D>) upVar->Rebin(histogram->GetRebin(), newNameUp, histogram->GetRebinVar());
        //            downVar = (std::shared_ptr<TH1D>) downVar->Rebin(histogram->GetRebin(), newNameDown, histogram->GetRebinVar());
        //        }
        //    }
        //} else if (name == "pdfShapeVar") {
        //    std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> histVars = buildPDFFromSumSquaredCollections(histogram, current, histNominal);
        //    upVar = histVars.first;
        //    downVar = histVars.second;
        //} else {
        //    std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> histVars = buildEnvelopeForProcess(histogram, current, histNominal);
        //    upVar = histVars.first;
        //    downVar = histVars.second;
        //}


        current = current->getNext();

        histCount++;
        procCount++;
    } 
}

std::shared_ptr<TH1D> Uncertainty::getFlatUncertainty(Histogram* histogram, Process* head, std::vector<std::shared_ptr<TH1D>>& histVec) {
    Process* current = head;
    //int histCount = 0;

    std::shared_ptr<TH1D> ret = std::make_shared<TH1D>(TH1D(*histVec[0]));
    TString histogramID = histogram->getID();
    ret->SetName(histogramID + name);
    ret->SetTitle(histogramID + name);

    std::vector<double> var(histVec[0]->GetNbinsX(), 0.);
    int procCount = 0;

    for (unsigned histCount = 0; histCount < histVec.size(); histCount++) {
        if (current->getName() != relevantProcesses[procCount]) {
            current = current->getNext();
            continue;
        }

        for (int bin = 1; bin < histVec[0]->GetNbinsX() + 1; bin++) {
            double variation = 0.;
            if (correlatedAmongProcesses) {
                double binContent = histVec[histCount]->GetBinContent( bin );
                variation = binContent * (flatUncertainty - 1.);
                var[bin - 1] += variation;
            } else {
                variation = histVec[histCount]->GetBinContent( bin )*( flatUncertainty - 1. );
                var[bin - 1] = variation * variation;
            }
        }

        current = current->getNext();
        procCount++;
    }

    for (int bin = 1; bin < histVec[0]->GetNbinsX() + 1; bin++) {
        //correlated case :
        
        if (correlatedAmongProcesses) {
            ret->SetBinContent(bin, var[bin - 1] * var[bin - 1]);
            //std::cout << "bin: " << bin << " error " << var[bin - 1] << std::endl;
        } else {
            ret->SetBinContent(bin, var[bin - 1]);
            //std::cout << "bin: " << bin << " error " << sqrt(var[bin - 1]) << std::endl;
        }
    }

    return ret;
}

std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> Uncertainty::getUpAndDownShapeUncertainty(Histogram* histogram, Process* head, std::vector<std::shared_ptr<TH1D>>& nominalHists) {
    // Loop processes, ask to add stuff
    Process* current = head;
    TString histogramID = histogram->getID();

    int histCount = 0;
    int procCount = 0;

    TString outputNameUp = outputName + "Up";
    TString outputNameDown = outputName + "Down";

    if (histogram->getPrintToFile()) {
        
        outfile->cd(histogram->getCleanName().c_str());
        gDirectory->mkdir(outputNameUp);
        gDirectory->mkdir(outputNameDown);
    }

    std::string up = "Up";
    std::string down = "Down";

    std::shared_ptr<TH1D> upVarReturn = nullptr;
    std::shared_ptr<TH1D> downVarReturn = nullptr;


    while (current) {
        // TODO: check if uncertainty needs this process, otherwise continue and put stuff to next one;
        if (current->getName() != relevantProcesses[procCount]) {
            histCount++;
            current = current->getNext();

            continue;
        }

        std::shared_ptr<TH1D> upVar = nullptr;
        std::shared_ptr<TH1D> downVar = nullptr;
        std::shared_ptr<TH1D> histNominal = nominalHists[histCount];

        if (! envelope || (envelope && ! buildEnvelope)) {
            upVar = current->getHistogramUncertainty(name, up, histogram, outputName, isEnvelope());
            downVar = current->getHistogramUncertainty(name, down, histogram, outputName, isEnvelope());

            if (histogram->HasRebin()) {
                if (histogram->GetRebinVar() == nullptr) {
                    upVar->Rebin(histogram->GetRebin());
                    downVar->Rebin(histogram->GetRebin());
                } else {
                    TString newNameUp = upVar->GetName();
                    upVar->SetName(newNameUp + TString("OLD"));
                    TString newNameDown = downVar->GetName();
                    downVar->SetName(newNameDown + TString("OLD"));
                    
                    upVar = std::make_shared<TH1D> (*(TH1D*) upVar->Rebin(histogram->GetRebin()-1, newNameUp, histogram->GetRebinVar()));
                    downVar = std::make_shared<TH1D> (*(TH1D*) downVar->Rebin(histogram->GetRebin()-1, newNameDown, histogram->GetRebinVar()));
                }
            }
        } else if (name == "pdfShapeVar") {
            std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> histVars = buildPDFFromSumSquaredCollections(histogram, current, histNominal);
            upVar = histVars.first;
            downVar = histVars.second;
        } else {
            std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> histVars = buildEnvelopeForProcess(histogram, current, histNominal);
            upVar = histVars.first;
            downVar = histVars.second;
        }

        if (histogram->getPrintToFile()) {
            outfile->cd();
            outfile->cd((histogram->getCleanName()).c_str());

            for (int j=1; j < upVar->GetNbinsX() + 1; j++) {
                if (upVar->GetBinContent(j) <= 0.) upVar->SetBinContent(j, 0.00001);
                if (downVar->GetBinContent(j) <= 0.) downVar->SetBinContent(j, 0.00001);
            }

            gDirectory->cd(outputNameUp);
            upVar->Write(current->getName());
            gDirectory->cd("..");

            gDirectory->cd(outputNameDown);
            downVar->Write(current->getName());
        }

        if (upVarReturn == nullptr) {
            upVarReturn = std::make_shared<TH1D>(TH1D(*upVar));
            downVarReturn = std::make_shared<TH1D>(TH1D(*downVar));

            if(! correlatedAmongProcesses ){
                upVarReturn->Multiply(upVarReturn.get());
                downVarReturn->Multiply(downVarReturn.get());
            }
        } else {
            if(! correlatedAmongProcesses ){
                upVar->Multiply(upVar.get());
                downVar->Multiply(downVar.get());

                upVarReturn->Add(upVar.get());
                downVarReturn->Add(downVar.get());
            } else {
                upVarReturn->Add(upVar.get());
                downVarReturn->Add(downVar.get());
            }
        }

        current = current->getNext();

        histCount++;
        procCount++;
    }

    return {upVarReturn, downVarReturn};
}