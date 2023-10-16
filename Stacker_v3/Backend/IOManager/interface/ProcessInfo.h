#ifndef PROCESSINFO_H
#define PROCESSINFO_H

// C++ includes
#include <iostream>

// ROOT includes
#include <TString.h>

class ProcessInfo {
    private:
        TString pCleanName;
        Color_t pColor;
        
        bool pIsSignal;
        bool pIsData;
    public:
        ProcessInfo(TString& name, std::string& color, std::string& type);
        ~ProcessInfo();

        // Getters
        Color_t const GetColor()             {return pColor;};
        TString const GetCleanName()         {return pCleanName;};
        bool const IsSignal()                {return pIsSignal;};
        bool const IsData()                  {return pIsData;};

        // Setters
        void SetColor(Color_t color)         {pColor = color;};
        void SetCleanName(TString cleanName) {pCleanName = cleanName;};
        void SetIsSignal(bool isSignal)      {pIsSignal = isSignal;};
        void SetIsData(bool isData)          {pIsData = isData;};

};



#endif