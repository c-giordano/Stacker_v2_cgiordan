#include "../interface/Stacker.h"

TCanvas* Stacker::getCanvas(TString& histID) {
	if (isRatioPlot) {
		gStyle->SetCanvasDefH(941);
	}
    return new TCanvas(histID + "_canvas");
}

TPad* Stacker::getPad(TString& histID, int position) {
	double margin = 0.013;
	TPad* outputPad;
	if (position == 0) {
		if (isRatioPlot) {
			outputPad = new TPad(histID + "_pad", histID + "_pad", 0., 0.3, 1., 1.);
			outputPad->SetBottomMargin(margin);		
		} else {
			outputPad = new TPad(histID + "_pad", histID + "_pad", 0., 0., 1., 1.);
		}
	}  else if (position == 1) {
		outputPad = new TPad(histID + "_pad", histID + "_pad", 0., 0., 1., 0.28);
		outputPad->SetTopMargin(margin);
		outputPad->SetBottomMargin(0.4);
	} else {
		outputPad = new TPad(histID + "_pad", histID + "_pad", 0., 0., 1., 1.);
	}

	return outputPad;
}

TLegend* Stacker::getLegend() {
    double x1 = 1. - gStyle->GetPadRightMargin() - 0.03 - 0.4; // 0.4;//0.65
    double y1 = 1. - gStyle->GetPadTopMargin() - 0.03 - 0.2;  // 0.69;
    double x2 = 1. - gStyle->GetPadRightMargin() - 0.03; // 1 - gStyle->GetPadLeftMargin() - gStyle->GetTickLength("X")*1.1 ;
    double y2 = 1. - gStyle->GetPadTopMargin() - 0.03; // 1 - gStyle->GetPadTopMargin() - gStyle->GetTickLength("Y")*1.1;//-0.06;
    int nColumns = 2; // 1 

    TLegend* legend = new TLegend(x1, y1, x2, y2, "", "NBNDC");
    legend->SetNColumns(nColumns);
    //legend->SetMargin(.4);
    //legend->SetColumnSeparation(0.1);
	legend->SetFillStyle(0);
	legend->SetTextSize(0.035);
	legend->SetTextFont(42);

    return legend;
}

TLatex* Stacker::getDatasetInfo(TPad* pad) {
    /*
    TLatex* datasetInfo = new TLatex();
    datasetInfo->SetTextSize(0.034);
    datasetInfo->SetNDC();
    datasetInfo->SetTextAlign(33);
    datasetInfo->SetTextFont(42);

    TString newString = intLumi + " fb^{-1} (13 TeV)";
    //newString += " fb^{-1} (13 TeV)";

    datasetInfo->DrawLatex(0.975, 0.99, newString);

    return datasetInfo;
    */

    TString lumiText = intLumi + " fb^{-1} (13 TeV)";
    TString extraText = "Work in progress";

    const float l = pad->GetLeftMargin() + 0.045 * (1 - pad->GetLeftMargin() - pad->GetRightMargin()); //+ gStyle->GetTickLength()*1.2; // 0.65;//
    const float t = pad->GetTopMargin();
    const float r = pad->GetRightMargin();
    //const float b = pad->GetBottomMargin();
    double pad_ratio = (pad->GetWh() * pad->GetAbsHNDC()) / (pad->GetWw() * pad->GetAbsWNDC());
    if (pad_ratio < 1.) {
        pad_ratio = 1.;
    }
    float CMSTextSize = pad->GetTopMargin() * 1. * pad_ratio;
    float lumiTextSize = pad->GetTopMargin() * 0.6 * pad_ratio;

	//float CMSTextOffset = pad->GetTopMargin()*0.2;
	float lumiTextOffset = pad->GetTopMargin()*0.2;
	
	pad->cd();
	//Define latex text to draw on plot
	TLatex* latex = new TLatex(l,1+lumiTextOffset*t,"CMS");
	latex->SetNDC();
	latex->SetTextAngle(0);
	latex->SetTextColor(kBlack); 

	latex->SetTextFont(61);
	latex->SetTextAlign(11);  // 13?
	latex->SetTextSize(CMSTextSize);
	//float cmsX = latex->GetXsize();
	double relPosY = 1. - t - 0.05 * (1. - t - pad->GetBottomMargin()); // 1 - (2*t)
	latex->DrawLatex(l, relPosY,"CMS"); 

	float extraTextSize = CMSTextSize*0.76;	 
	latex->SetTextFont(52);
	latex->SetTextSize(extraTextSize);
	latex->SetTextAlign(11);
	//std::cout << extraText << " " << l + 1.2 * cmsX << std::endl;
	// using cmsX gave strange results
	latex->DrawLatex(l, 1-(2 * t) - (extraTextSize*1.1), extraText);

	latex->SetTextFont(42);
	latex->SetTextAlign(31);
	latex->SetTextSize(lumiTextSize);  
	latex->DrawLatex(1-r,1-t+lumiTextOffset,lumiText);

	return latex;
}
