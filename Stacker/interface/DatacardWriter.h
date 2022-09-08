#ifndef DATACARDWRITER_H
#define DATACARDWRITER_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>

#include "Uncertainty.h"
#include "Histogram.h"
#include "ProcessList.h"

class DatacardWriter
{
private:
    std::ofstream datacard;
    std::string datacardName;
    std::string yearID;

    TFile* outfile;
    ProcessList* allProc;
    //std::vector<Process*> allProc;
    std::vector<Histogram*> allHistograms;
    std::vector<std::vector<std::shared_ptr<TH1D>>> allHistogramTH1Ds;

    Process* head;
    
    bool isData = false;
    Process* dataProc = nullptr;
public:
    DatacardWriter(std::string yearID, ProcessList* allProc, std::vector<Histogram*> histVec, TFile* outfile, Process* dataProcNew, std::string nameOverwrite);
    ~DatacardWriter();

    void addData(Process* dataProcNew) {
        if (dataProcNew != nullptr) {
            dataProc = dataProcNew;
            isData = true;
        }
    }

    void initDatacard();

    void writeEmptyLine(unsigned length);

    void writeInit();
    void writeShapeSource();
    void writeDataObs();
    void writeProcessHeader();
    void writeUncertainties(Uncertainty* uncertainty, bool eraSpecific = false);
    void write1718Uncertainty(Uncertainty* uncertainty);
    void writeMCStats();

    bool containsProcess(std::vector<TString>& vector, TString& process);

    void WriteDatacardVariations();
    void WriteDatacardVariationInit(Histogram* histogram);
    void WriteDatacardVariationUncertainties(std::ofstream& datacardVar, Uncertainty *uncertainty, Histogram *histogram, bool eraSpecific = false);
    void WriteDatacardVariation1718Uncertainties(std::ofstream& datacardVar, Uncertainty *uncertainty, Histogram *histogram);

};

#endif