#include "../interface/thTools.h"

void normalizeHistograms(std::vector<std::shared_ptr<TH1D>>& histos) {
    for (auto curr : histos) {
        curr->Scale(1/curr->Integral());
    }
}

std::shared_ptr<TH1D> sumVector(std::vector<std::shared_ptr<TH1D>>& histoVec) {
    std::shared_ptr<TH1D> sum = nullptr;
    for (unsigned i = 0; i < histoVec.size(); i++) {
        std::shared_ptr<TH1D> currHist = histoVec[i];
        if (i == 0) {
            sum = std::make_shared<TH1D>(TH1D(*currHist));
        } else {
            sum->Add(currHist.get());
        }
    }
    return sum;
}