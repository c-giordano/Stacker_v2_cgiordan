#include "../interface/Stacker.h"
#include <functional>

void Stacker::GenerateSFs(std::string& SFFile) {
    // read input file ->parse to a vector which histograms we need
    // for each histogram, generate it
    std::ifstream infile(SFFile);

    std::string line;
    if (! infile.is_open()) {
        std::cout << "ERROR: SF inputfile " << SFFile << " not found" << std::endl;
        exit(1);
    }

    while (getline(infile, line)) {
        std::vector<Histogram*>::iterator it;

        // split line in process and histogram
        std::string processName;
        std::string histogram;
        std::stringstream stream(line);
        stream >> histogram >> processName;

        TString processNameTString = processName;
        it = std::find_if(histogramVec.begin(), histogramVec.end(), std::bind(Histogram::searchHist, std::placeholders::_1, histogram));
        GenerateSF(*it, processNameTString);
    }
}

void Stacker::GenerateSF(Histogram* histogram, TString& processName) {
    // generate stack + data info
    // also plot the inputdistribution and the output somewhere
    // generate outputroot file for each distribution. Renaming will be done manually anyway
    printHistogram(histogram); // improve on how outputdir is specified

    TString histID = histogram->getID();
    TLegend* legend = getLegend();
    THStack* histStack = new THStack(histID, histID);
    std::vector<std::shared_ptr<TH1D>>* signalVector = new std::vector<std::shared_ptr<TH1D>>;

    std::vector<std::shared_ptr<TH1D>> histVec = processes->fillStack(histStack, histogram, legend, outputfile, signalVector, nullptr);

    std::shared_ptr<TH1D> dataHistogram;
    if (getFakeData()) {
        dataHistogram = sumVector(histVec);
        dataHistogram->SetTitle("Data (expected)");
        dataHistogram->SetName("Data (expected)");

        for (int bin = 1; bin < dataHistogram->GetNbinsX() + 1; ++bin) {
            dataHistogram->SetBinError(bin, sqrt(dataHistogram->GetBinContent(bin)));
        }
    } else {
        dataHistogram = dataProcess->getHistogram(histogram);
        dataHistogram->SetTitle("Data");
        dataHistogram->SetName("Data");
    }

    std::shared_ptr<TH1D> sf = std::make_shared<TH1D>(TH1D(*dataHistogram));
    // take a deep copy of the process for which we extract a SF
    // histvector should be in the order of the processes
    // remove from histvec, sum histvec, extract from data
    // generate SF with contribution and backgrounds
    size_t index = 0;
    Process* currProc = processes->getHead();
    while (currProc && currProc->getName() != processName) {
        currProc = currProc->getNext();
        index++;
    }
    if (currProc == nullptr) {
        std::cerr << "SF generation failed. Process " << processName << " was not found." << std::endl;
        exit(1);
    }

    std::shared_ptr<TH1D> releventContribution = std::make_shared<TH1D>(TH1D(*histVec[index]));
    std::cout << releventContribution->GetName() << std::endl;
    histVec.erase(histVec.begin()+index);

    std::shared_ptr<TH1D> sum = sumVector(histVec);
    sf->Add(sum.get(), -1.);
    

    sf->SetName("SF_" + histogram->getID());
    sf->SetTitle("SF_" + histogram->getID());

    sf->Divide(releventContribution.get());
    
    DrawSF(sf);

    TFile* sfOutput = new TFile("ScaleFactors/Output/SF_" + histogram->getID() + ".root", "RECREATE");

    sfOutput->cd();

    sf->Write("SF");
    
    sfOutput->Close();
}

void Stacker::DrawSF(std::shared_ptr<TH1D> sfHistogram) {
    // canvas, pad, draw with usual setting but add text on it
    // lumi not specified generally speaking
    TString sfName = sfHistogram->GetName();
    isRatioPlot = false;
    TCanvas* canv = getCanvas(sfName);
    canv->Draw();
    canv->cd(); 
    TPad* pad = getPad(sfName, 0);
    pad->Draw();
    pad->cd();

    //gStyle->SetMarkerSize(2.);
    sfHistogram->SetYTitle("Scale factor");

    sfHistogram->Draw();
    for (int i=1; i<=sfHistogram->GetNbinsX(); i++) {
        auto t = new TText(sfHistogram->GetXaxis()->GetBinCenter(i), sfHistogram->GetBinContent(i)+0.5, Form("%4.2f",sfHistogram->GetBinContent(i)));
        t->SetTextAlign(22);
        t->Draw("SAME");
    }

    TString extraText = "Work in progress";

    const float l = pad->GetLeftMargin() + gStyle->GetTickLength()*1.2;
  	const float t = pad->GetTopMargin();
  	const float r = pad->GetRightMargin();
  	//const float b = pad->GetBottomMargin();

	float CMSTextSize = pad->GetTopMargin()*0.75;
	float lumiTextSize = pad->GetTopMargin()*0.6;

	//float CMSTextOffset = pad->GetTopMargin()*0.2;
	float lumiTextOffset = pad->GetTopMargin()*0.2;
	
	pad->cd();
	//Define latex text to draw on plot
	TLatex* latex = new TLatex(l,1+lumiTextOffset*t,"CMS");
	latex->SetNDC();
	latex->SetTextAngle(0);
	latex->SetTextColor(kBlack); 

	latex->SetTextFont(61);
	latex->SetTextAlign(11); 
	latex->SetTextSize(CMSTextSize);
	//float cmsX = latex->GetXsize();
	latex->DrawLatex(l,1  - (2 * t),"CMS");

	float extraTextSize = CMSTextSize*0.76;	 
	latex->SetTextFont(52);
	latex->SetTextSize(extraTextSize);
	latex->SetTextAlign(11);
	//std::cout << extraText << " " << l + 1.2 * cmsX << std::endl;
	// using cmsX gave strange results
	latex->DrawLatex(l, 1-(2 * t) - (extraTextSize*1.1), extraText);

    canv->Print("ScaleFactors/Output/" + sfName + ".pdf");
    canv->Print("ScaleFactors/Output/" + sfName + ".png");
}
