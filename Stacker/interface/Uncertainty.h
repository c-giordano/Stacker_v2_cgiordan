#ifndef UNCERTAINTY_H
#define UNCERTAINTY_H

#include <string>
#include <vector>

#include <TString.h>
#include <TH1.h>

#include "Process.h"
#include "Histogram.h"

class Uncertainty {
    private:
        std::string name;
        std::string outputName;
        TString nameUp;
        TString nameDown;
        bool flat;
        bool envelope;
        bool buildEnvelope = true;
        bool correlatedAmongProcesses;
        bool eraSpecific;
        bool isBoth;

        double flatUncertainty; // value used in histograms
        double flatUncertaintyEra=0.0; // for datacards: this era of data
        double flatUncertaintyAll=0.0; // for datacards: all eras
        double flatUncertainty1718=0.0; // for datacards: all eras

        std::vector<TString> relevantProcesses;

        TFile* outfile;

        Uncertainty* next = nullptr;

    public:
        Uncertainty(std::string& name, bool flat, bool envelope, bool corrProcess, bool eraSpec, std::vector<TString>& processes, TFile* outputfile);
        ~Uncertainty();

        void setFlatRate(double rate) {flatUncertainty = rate;}
        
        void setFlatRateEra(double rate) {flatUncertaintyEra = rate;}
        void setFlatRateAll(double rate) {flatUncertaintyAll = rate;}
        void setFlatRate1718(double rate) {flatUncertainty1718 = rate;}
        void setEraSpec(bool newEraSpec) {eraSpecific = newEraSpec;}
        void setBoth(bool newBoth) {isBoth = newBoth;}

        void setCorrProcess(bool newCorr) {correlatedAmongProcesses = newCorr;}
        void setOutputName(std::string newName) {outputName = newName;}
        std::string getOutputName() const {return outputName;}


        bool isFlat() {return flat;}
        bool isEnvelope() {return envelope;}

        Uncertainty* getNext() {return next;}
        void setNext(Uncertainty* newNext) {next = newNext;}

        std::shared_ptr<TH1D> getShapeUncertainty(Histogram* histogram, Process* head, std::vector<std::shared_ptr<TH1D>>& histVec);
        std::shared_ptr<TH1D> getEnvelope(Histogram* histogram, Process* head, std::vector<std::shared_ptr<TH1D>>& histVec);
        std::shared_ptr<TH1D> getSumSquaredEnvelope(Histogram* histogram, Process* head, std::vector<std::shared_ptr<TH1D>>& histVec);

        std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> buildEnvelopeForProcess(Histogram* histogram, Process* currentProcess, std::shared_ptr<TH1D> nominalHist);
        std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> buildSumSquaredEnvelopeForProcess(Histogram* histogram, Process* currentProcess, std::shared_ptr<TH1D> nominalHist);
        std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> buildPDFFromSumSquaredCollections(Histogram* histogram, Process* currentProcess, std::shared_ptr<TH1D> nominalHist);

        std::pair<std::shared_ptr<TH1D>, std::shared_ptr<TH1D>> getUpAndDownShapeUncertainty(Histogram* histogram, Process* head, std::vector<std::shared_ptr<TH1D>>& nominalHists);
        std::shared_ptr<TH1D> getFlatUncertainty(Histogram* histogram, Process* head, std::vector<std::shared_ptr<TH1D>>& histVec);
        std::shared_ptr<TH1D> getUncertainty(Histogram* histogram, Process* head, std::vector<std::shared_ptr<TH1D>>& histVec);

        void printOutShapeUncertainty(Histogram* histogram, Process* head);

        std::string& getName() {return name;}
        std::vector<TString>& getRelevantProcesses() {return relevantProcesses;}
        
        bool isEraSpecific() {return eraSpecific;}
        bool isBothEraAndFull() {return isBoth;}
        bool is1718Specific() {return (flatUncertainty1718 != 0.);}
        double getFlatRateAll() const {return flatUncertaintyAll;}
        double getFlatRateEra() const {return flatUncertaintyEra;}
        double getFlatRate1718() const {return flatUncertainty1718;}
        std::string getOutputNam() const {return outputName;}

        bool getCorrelatedAmongProcesses() const {return correlatedAmongProcesses;}

};


#endif