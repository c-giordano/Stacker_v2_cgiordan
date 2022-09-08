#include "../interface/Stacker.h"

void Stacker::printAllHistograms() {
    int tempCount =0 ;
    
    if (onlyDC) return;
    for (auto histogramID : histogramVec) {
        tempCount++;
        histogramID->setPrintToFile(false);
        //if (tempCount == 50) break;

        if (! onlyDC || (onlyDC && histogramID->getPrintToFile())) {
            printHistogram(histogramID);
        }
    }
}

void Stacker::printHistogram(Histogram* hist) {
    TString histID = hist->getID();

    THStack* histStack = new THStack(histID, histID);
    TLegend* legend = getLegend();
    std::vector<std::shared_ptr<TH1D>>* signalVector = new std::vector<std::shared_ptr<TH1D>>;
    std::shared_ptr<TH1D>* sysUnc = new std::shared_ptr<TH1D>();
    *sysUnc = nullptr; 
    std::vector<std::shared_ptr<TH1D>> histVec = processes->fillStack(histStack, hist, legend, outputfile, signalVector, sysUnc);
    
    if (onlyDC) return;

    TCanvas* canv = getCanvas(histID);
    canv->Draw();
    canv->cd();

    std::shared_ptr<TH1D> dataHistogram = nullptr;
    if (getData()) {
        // get data histogram. Either from specialised process(list) or fake it
        // pass to stackfiller and to ratiodrawer
        if (getFakeData()) {
            dataHistogram = sumVector(histVec);
            dataHistogram->SetTitle("Data (expected)");
            dataHistogram->SetName("Data (expected)");

            for (int bin = 1; bin < dataHistogram->GetNbinsX() + 1; ++bin) {
                dataHistogram->SetBinError(bin, sqrt(dataHistogram->GetBinContent(bin)));
            }
        } else {
            dataHistogram = dataProcess->getHistogram(hist);
            if (verbose) {
                std::cout << "Obs";
                if (veryVerbose) {
                    std::cout << " & ";
                    for (int i=1; i < dataHistogram->GetNbinsX() + 1; i++) {
                        std::cout << dataHistogram->GetBinContent(i) << " & ";
                    }
                    std::cout << std::endl;
                } else {
                    std::cout << ": " << dataHistogram->Integral() << " events" << std::endl;
                }
            }
            dataHistogram->SetTitle("Data");
            dataHistogram->SetName("Data");
        }

        legend->AddEntry(dataHistogram.get(), dataHistogram->GetName());
        // change dataHistogram settings
    }

    TPad** mainPad = new TPad*();
    std::shared_ptr<TH1D> totalUnc = drawStack(hist, histStack, histVec, sysUnc, dataHistogram, mainPad);

    drawSignalYield(legend, *signalVector);
    legend->Draw();
    
    canv->cd();

    TPad** smallPad = new TPad*(); 
    std::shared_ptr<TH1D> ratioPlot = nullptr;
    if (getData()) {
        if (! totalUnc) {
            std::shared_ptr<TH1D> allHistograms = sumVector(histVec);
            ratioPlot = drawRatioData(hist, allHistograms, dataHistogram, smallPad);

        } else {
            ratioPlot = drawRatioData(hist, totalUnc, dataHistogram, smallPad);
        }
    } else {
        ratioPlot = drawRatioMC(hist, histVec, *signalVector, smallPad);
    }

    // auto resize axis
    
    TH1* combiHist = (TH1*) histStack->GetStack()->Last();
    double xmin = combiHist->GetBinLowEdge(1);
    double xmax = combiHist->GetBinLowEdge(combiHist->GetNbinsX()) + combiHist->GetBinWidth(combiHist->GetNbinsX());
    bool change = false;
    
    double MinContent = 0.0002; // too small. Fix later because this does not take ratioplots into account
    int counter = 1;
    double currentBinContent = combiHist->GetBinContent(counter);
    if (dataHistogram) currentBinContent += dataHistogram->GetBinContent(counter);

    while (currentBinContent <= MinContent && counter <= combiHist->GetNbinsX()) {
        counter++;
        xmin = combiHist->GetBinLowEdge(counter);
        change = true;

        currentBinContent = combiHist->GetBinContent(counter);
        if (dataHistogram) currentBinContent += dataHistogram->GetBinContent(counter);
    }

    double MaxContent = 0.0002;
    counter = combiHist->GetNbinsX();

    currentBinContent = combiHist->GetBinContent(counter);
    if (dataHistogram) currentBinContent += dataHistogram->GetBinContent(counter);

    while (currentBinContent <= MaxContent && counter >= 1) {
        counter--;
        xmax = combiHist->GetBinLowEdge(counter) + combiHist->GetBinWidth(counter);
        change = true; 

        currentBinContent = combiHist->GetBinContent(counter);
        if (dataHistogram) currentBinContent += dataHistogram->GetBinContent(counter);
    }

    if (change) {
        // temp disabled -> dont really like the effect it has on some plots. Definitely a want for later versions but not now.
        //histStack->GetXaxis()->SetRangeUser(xmin, xmax);
        //if (ratioPlot) ratioPlot->GetXaxis()->SetRangeUser(xmin, xmax);
    }

    (*mainPad)->Update();
    (*mainPad)->Modified();

    if (*smallPad) {
        (*smallPad)->Update();
        (*smallPad)->Modified();
    }

    std::string fullPath = pathToOutput;
    if (runT2B) {
        std::string id = histID.Data();
        fullPath += getChannel(id);
        fullPath += "/";
    }

    canv->Print(fullPath + histID + ".png");
    canv->Print(fullPath + histID + ".pdf");

    delete mainPad;
    delete smallPad;
}


std::shared_ptr<TH1D> Stacker::drawStack(Histogram* hist, THStack* histStack, std::vector<std::shared_ptr<TH1D>>& histVec, std::shared_ptr<TH1D>* sysUnc, std::shared_ptr<TH1D> data, TPad** mainPad) {
    stackSettingsPreDraw(histStack, histVec);
    TString histID = hist->getID();
    TPad* pad = getPad(histID, 0);
    *mainPad = pad;

    pad->Draw();
    pad->cd();

    //stackSettingsPreDraw(histStack, histVec);

    histStack->Draw(drawOpt.c_str());

    stackSettingsPostDraw(pad, histStack, hist, histVec[0], data);

    std::shared_ptr<TH1D> allHistograms = sumVector(histVec);
    std::shared_ptr<TH1D> totalUnc = nullptr;
    if (*sysUnc) {
        totalUnc = std::make_shared<TH1D>(TH1D(*allHistograms));
        
        //std::cout << "printing uncertainties:\t";
        for(int bin = 1; bin < totalUnc->GetNbinsX() + 1; ++bin){
            double statError = allHistograms->GetBinError(bin);
            double systError = (*sysUnc)->GetBinContent(bin); // is already squared
            totalUnc->SetBinError(bin, sqrt( statError*statError + systError) );
            //totalUnc->SetBinError(bin, sqrt(systError) );

            //std::cout << totalUnc->GetBinError(bin) << "\t";

        }
        //std::cout << std::endl;
        totalUnc->SetFillStyle(3244); //3005  3244
        totalUnc->SetFillColor(kGray+2);
        totalUnc->SetMarkerStyle(0); //1
        totalUnc->Draw("E2 SAME");
    }

    if (data != nullptr) {
        data->SetMarkerColor(kBlack);
        data->SetFillColor(0);
        data->SetLineColor(kBlack);
        data->SetLineWidth(2);
        data->SetMarkerSize(1.1);

        data->Draw("E1 X0 SAME"); // E1 for more barsssss
    }


    TH1* combiHist = (TH1*) histStack->GetStack()->Last();
    double xmin = combiHist->GetBinLowEdge(1);
    double xmax = combiHist->GetBinLowEdge(combiHist->GetNbinsX()) + combiHist->GetBinWidth(combiHist->GetNbinsX());
    bool change = false;
    
    double MinContent = 0.00002; // too small. Fix later because this does not take ratioplots into account
    int counter = 1;

    while (combiHist->GetBinContent(counter) <= MinContent && counter <= combiHist->GetNbinsX()) {
        counter++;
        xmin = combiHist->GetBinLowEdge(counter);
        change = true;
    }

    double MaxContent = 0.00002;
    counter = combiHist->GetNbinsX();

    while (combiHist->GetBinContent(counter) <= MaxContent && counter >= 1) {
        counter--;
        xmax = combiHist->GetBinLowEdge(counter) + combiHist->GetBinWidth(counter);
        change = true; 
    }

    if (change) {
        histStack->GetXaxis()->SetRangeUser(xmin, xmax);
    }

    pad->Update();
    pad->Modified();

    TLatex* info = getDatasetInfo(pad);

    delete info;

    return totalUnc;
}

void Stacker::drawSignalYield(TLegend* legend, std::vector<std::shared_ptr<TH1D>>& signalVec) {
    if (! isSignalLine) return;

    std::shared_ptr<TH1D> signalTotal = sumVector(signalVec);
    signalTotal->SetTitle("Signal yield");

    legend->AddEntry(signalTotal.get(), "Signal x30");
    
    signalTotal->SetFillColor(0);
    signalTotal->SetLineColor(633);
    signalTotal->SetLineWidth(5);
    signalTotal->SetMarkerStyle(0);
    signalTotal->Scale(30.);
    signalTotal->Draw("SAME HIST");
}


std::shared_ptr<TH1D> Stacker::drawRatioMC(Histogram* hist, std::vector<std::shared_ptr<TH1D>>& histoVec, std::vector<std::shared_ptr<TH1D>>& signalVec, TPad** smallPadPtr) {
    /* 
    TODO:
        Make a function deciding if data or not, let it decide after what to do... ofzo
    */
    if (! isRatioPlot) return nullptr;

    TString histID = hist->getID() + "_ratio";
    TPad* smallPad = getPad(histID, 1);
    *smallPadPtr = smallPad;

    smallPad->Draw();
    smallPad->cd();

    std::shared_ptr<TH1D> signalTotal = sumVector(signalVec);
    std::shared_ptr<TH1D> allHistograms = sumVector(histoVec);
    
    signalTotal->Divide(allHistograms.get());

    signalTotal->SetTitleSize(0.192, "X");
    signalTotal->SetTitleSize(0.17, "Y");

    signalTotal->SetTitle("");
    signalTotal->SetYTitle("Signal fraction");

    //signalTotal->SetMaximum(1.1);
    signalTotal->SetFillColor(0);
    signalTotal->SetLineWidth(2);
    signalTotal->SetLineColor(1);
    signalTotal->SetLabelSize(0.16, "XY");
    signalTotal->SetTitleOffset(0.38, "Y");
    signalTotal->SetTickLength(0.096, "X");
    signalTotal->SetNdivisions(504, "Y");

    signalTotal->Draw("HIST");

    //signalTotal->GetXaxis()->SetTitle(first->GetXaxis()->GetTitle());

    if (hist->getXBinLabels()) {
        std::vector<std::string>* bins = hist->getXBinLabels();
        for (unsigned i = 1; i != bins->size() + 1; i++) {
            signalTotal->GetXaxis()->SetBinLabel(i, TString(bins->at(i - 1)));
        }
    }

    smallPad->Update();
    smallPad->Modified();
    //signalTotal->UseCurrentStyle();

    return signalTotal;
}

std::shared_ptr<TH1D> Stacker::drawRatioData(Histogram* hist, std::shared_ptr<TH1D> uncHist, std::shared_ptr<TH1D> data, TPad** smallPadPtr) {
    if (! isRatioPlot) return nullptr;
    
    TString histID = hist->getID() + "_ratio";
    TPad* smallPad = getPad(histID, 1);
    *smallPadPtr = smallPad;

    smallPad->Draw();
    smallPad->cd();

    std::shared_ptr<TH1D> dataTotal = std::make_shared<TH1D>(TH1D(*data));
    std::shared_ptr<TH1D> mcTotal = std::make_shared<TH1D>(TH1D(*uncHist));
    
    dataTotal->Divide(uncHist.get());

    dataTotal->SetTitleSize(0.192, "X");
    dataTotal->SetTitleSize(0.17, "Y");

    dataTotal->SetTitle("");
    dataTotal->SetYTitle("Data/pred.");

    //signalTotal->SetMaximum(1.1);
    dataTotal->SetFillColor(0);
    dataTotal->SetLineWidth(2);
    dataTotal->SetLineColor(1);
    dataTotal->SetLabelSize(0.16, "XY");
    dataTotal->SetTitleOffset(0.38, "Y");
    dataTotal->SetTickLength(0.096, "X");
    dataTotal->SetNdivisions(504, "Y");

    dataTotal->Draw("E1 X0"); // E1 for more barsssss

    //signalTotal->GetXaxis()->SetTitle(first->GetXaxis()->GetTitle());

    if (hist->getXBinLabels()) {
        std::vector<std::string>* bins = hist->getXBinLabels();
        for (unsigned i = 1; i != bins->size() + 1; i++) {
            dataTotal->GetXaxis()->SetBinLabel(i, TString(bins->at(i - 1)));
        }
    }


    int nrBins = mcTotal->GetNbinsX();
    for (int i = 1; i < nrBins + 1; i++) {
        dataTotal->SetBinError(i, sqrt(data->GetBinContent(i)) / uncHist->GetBinContent(i));

        if (mcTotal->GetBinContent(i) <= 0.0001) {
            mcTotal->SetBinError(i, 0.00001);
        } else {
            mcTotal->SetBinError(i, mcTotal->GetBinError(i) / mcTotal->GetBinContent(i));
        }
        mcTotal->SetBinContent(i, 1.);
    }
    mcTotal->Draw("SAME E2");

    TLine* line = new TLine(dataTotal->GetBinLowEdge(1), 1., dataTotal->GetXaxis()->GetBinUpEdge(dataTotal->GetNbinsX()), 1.);
    line->Draw("SAME");

    if (dataTotal->GetMaximum() > 2. || dataTotal->GetMaximum() + sqrt(dataTotal->GetMaximum()) > 2.) {
    }

    if (dataTotal->GetMinimum() < 0. || dataTotal->GetMinimum() - sqrt(dataTotal->GetMinimum()) < 0.) {
    }

    dataTotal->SetMaximum(2.45); 
    dataTotal->SetMinimum(0.5);

    smallPad->Update();
    smallPad->Modified();
    //signalTotal->UseCurrentStyle();

    return dataTotal;
}
