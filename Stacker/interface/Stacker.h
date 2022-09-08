#ifndef STACKER_H
#define STACKER_H

#include <TROOT.h>
#include <TH1.h>
#include <THStack.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TGaxis.h>
#include <TColor.h>
#include <TString.h>
//#include <TRatioPlot.h>

#include <algorithm>

#include "ProcessList.h"
#include "Process.h"
#include "Histogram.h"
#include "Histogram2D.h"
#include "DatacardWriter.h"

#include "../../Helpers/interface/ParseTools.h"
#include "../../Helpers/interface/StringTools.h"
#include "../../Helpers/interface/thTools.h"
#include "../../Helpers/interface/Others.h"

#include "../../Styles/tdrStyle.h"
/*
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
*/
// Use global setting file which keeps track of colors, order, definitions of components... -> SettingParser
// Keep track of which files belong to which process with which color -> other class? or linked list with entry for each process and vector for applicable files
// 

class Stacker {
    private:

        DatacardWriter* dcwriter;
        ProcessList* processes;
        Process* dataProcess = nullptr;
        std::vector<Histogram*> histogramVec;
        std::vector<Histogram2D*> histogramVec2D;
        // settings
        int canvWidthX = 600; // Possibly changing directly using gStyle... or save in style object being made. 
        int canvWidthY = 600;

        // Build components of histogram based on settings
        TCanvas* getCanvas(TString& histID);
        TPad* getPad(TString& histID, int position=0);
        TLegend* getLegend();
        TLatex* getDatasetInfo(TPad* currentPad);

        // Main root file
        TFile* inputfile;
        std::vector<TFile*> inputfiles;
        // outputfile (structured for combine)
        TFile* outputfile = nullptr;

        // Histograms to consider (of all processes)

        // General settings
        std::string intLumi = "";
        std::string drawOpt = "HIST";
        std::string yAxisOverride = "";
        std::string yearID;

        bool noStack = false;

        bool verbose = false;
        bool veryVerbose = false;
        bool uncertainties = false;
        bool useData = false;
        bool fakeData = false;
        bool onlyDC = false;
        bool oldStuff = false;

        bool runT2B = true;
        std::string pathToOutput;
        std::string altOutput = "";

        bool isRatioPlot = false;
        bool isSignalLine = false;

    public:
        Stacker(const char* rootFile, std::string& settingFile);
        Stacker(std::vector<std::string>& cmdArgs);

        void ReadSettingFile(std::string& settingFile);
        ~Stacker();

        void setVerbosity(bool verbosity) {
            verbose = verbosity;
            processes->setVerbosity(verbosity);
        }
        void setVeryVerbosity(bool verbosity) {
            veryVerbose = verbosity;
            processes->setVeryVerbosity(veryVerbose);
        }
        //void useUncertainties(bool isUncertainties) {uncertainties = isUncertainties;}
        void readUncertaintyFile(std::string& filename);
        void useT2B(bool isT2B) {
            runT2B = isT2B;
            if (runT2B) pathToOutput = "/user/nivanden/public_html/Most_recent_plots/";
            else pathToOutput = "Output/";
        }
        void setData(bool dataExists) {useData = dataExists;}
        bool getData() {return useData;}
        void readData(std::vector<std::string>& dataFile, unsigned i);

        void setFakeData(bool isFD) {fakeData = isFD;}
        bool getFakeData() {return fakeData;}

        void setOnlyDC(bool set) {onlyDC = set;}
        bool getOnlyDC() {return onlyDC;}

        void printAllHistograms();
        void printHistogram(Histogram* histID);
        std::shared_ptr<TH1D> drawStack(Histogram* hist, THStack* histStack, std::vector<std::shared_ptr<TH1D>>& histVec, std::shared_ptr<TH1D>* sysUnc, std::shared_ptr<TH1D> data, TPad** mainPad);
        void drawSignalYield(TLegend* legend, std::vector<std::shared_ptr<TH1D>>& signalVec);
        std::shared_ptr<TH1D> drawRatioMC(Histogram* hist, std::vector<std::shared_ptr<TH1D>>& histoVec, std::vector<std::shared_ptr<TH1D>>& signalVec, TPad** smallPadPtr);
        std::shared_ptr<TH1D> drawRatioData(Histogram* hist, std::shared_ptr<TH1D> uncHist, std::shared_ptr<TH1D> data, TPad** smallPadPtr);

        void drawUncertaintyImpacts(Histogram* hist, std::vector<std::string>& uncToDraw);
        void drawAllUncertaintyImpacts();
        
        std::vector<std::shared_ptr<TH1D>> fillStack(THStack* stack, TString& histogramID, TLegend* legend, TFile* outfile);
        
        void printAll2DHistograms();
        void print2DHistogram(Histogram2D* histID);

        void setLumi(std::string& lumiSetting);
        void setDrawOpt(std::string& drawSetting);

        void stackSettingsPreDraw(THStack* stack, std::vector<std::shared_ptr<TH1D>>& histVec);
        void stackSettingsPostDraw(TPad* pad, THStack* stack, Histogram* hist, std::shared_ptr<TH1> first, std::shared_ptr<TH1> dataHist);

        void initDatacard() {
            dcwriter->initDatacard();
            dcwriter->writeUncertainties(processes->getUncHead());
        }
        DatacardWriter* GetDCWriter() {return dcwriter;}
        void SaveToVault();

        void plotDifference(std::vector<std::string>& argvStr);

        void GenerateSFs(std::string& SFFile);
        void GenerateSF(Histogram* histogram, TString& processName);
        void DrawSF(std::shared_ptr<TH1D> sfHistogram);

};

#endif