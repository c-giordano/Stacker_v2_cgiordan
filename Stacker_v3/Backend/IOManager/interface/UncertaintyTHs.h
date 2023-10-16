#ifndef UNCERTAINTYTHS_H
#define UNCERTAINTYTHS_H

// C++ Includes
#include <vector>
#include <memory>

// ROOT Includes
#include <TH1.h>
#include <THStack.h>

// Framework
#include "Uncertainty.h"

class UncertaintyTHs {
    private:
        // Pointer to the object that made this object
        Uncertainty* pUncertainty;

        // THs related to the processes        
        std::vector<std::shared_ptr<TH1D>>* pProcessTHs;
        
    public:
        UncertaintyTHs(/* args */);
        ~UncertaintyTHs();

};

#endif