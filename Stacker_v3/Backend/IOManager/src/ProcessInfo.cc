#include "../interface/ProcessInfo.h"

ProcessInfo::ProcessInfo(TString& name, std::string& color, std::string& type) {
    pIsSignal = false;
    pIsData = false;

    if (type.find('S') != std::string::npos) {
        pIsSignal = true;
    } else if (type.find('D') != std::string::npos) {
        pIsData = true;
    } else if (type.find('B') == std::string::npos) {
        std::cerr << "Error: type identifier '" << type << "' unknown" << std::endl;
    }

    pColor = std::stoi(color);

    pCleanName = cleanTString(name);
}