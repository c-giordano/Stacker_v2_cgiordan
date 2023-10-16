#ifndef UNCERTAINTY_H
#define UNCERTAINTY_H

class Uncertainty {
    private:
        // always need a name
        std::string name;

        // bunch of settings deciding a bunch of things. 
        // Should the shape uncertainties be built on top of this?
        // either way: some of these variables might be outdated.
        // Tactic comp to v2: use multple for correlated things of uncertainty
        bool flat;
        bool envelope;
        bool buildEnvelope;
        bool correlatedAmongProcesses;
        bool eraSpecific;

        double flatUncertainty;

        std::vector<TString> relevantProcesses;
    public:
        Uncertainty(/* args */);
        ~Uncertainty();


        // Build Uncertainty TH
        UncertaintyTHs* BuildUncertaintyTHs();
};

#endif
