#include "../interface/SubProcess.h"

SubProcess::SubProcess(TString& name, std::vector<TFile*>& inputfiles) : pName(name), pInputfiles(inputfiles) {
    pSubdirectoriesPerFile = new std::vector<std::vector<const char*>*>();
    
    for (auto it : inputfiles) {    
        it->cd("Nominal");

        if (! gDirectory->GetDirectory(pName)) {
            pSubdirectoriesPerFile->push_back(nullptr);
            continue;
        }
        
        std::vector<const char*>* subdirectories = new std::vector<const char*>;
        gDirectory->cd(pName);
        TList* keys = gDirectory->GetListOfKeys();
        for(const auto&& obj: *keys) {
            if (std::string(((TKey*) obj)->GetClassName()) != "TDirectoryFile" && std::string(((TKey*) obj)->GetClassName()) != "TDirectory") continue;
            subdirectories->push_back(obj->GetName());
        }
        if (subdirectories->size() == 0 && keys->GetEntries() != 0) {
            // identify cases where we should use "compatibilitymode" automatically
            subdirectories->push_back("../"+pName);
            // should do the trick
        }

        pSubdirectoriesPerFile->push_back(subdirectories);
    }
}

SubProcess::~SubProcess() {
    for (auto ptr : *pSubdirectoriesPerFile) {
        for (auto subPtr : *ptr) {
            delete subPtr; // not sure if this is correct...
        }
        delete ptr;
    }

    delete pSubdirectoriesPerFile;
}

std::shared_ptr<TH1D> SubProcess::ExtractHistogramFromFiles(TString& id, TString& type, TString& subdirStack) {
    std::shared_ptr<TH1D> ret = nullptr; 

    for (unsigned i = 0; i < pInputfiles.size(); i++) {
        TFile* currFile = pInputfiles[i];
        std::vector<const char*>* subdirectories = pSubdirectoriesPerFile->at(i);

        if (subdirectories == nullptr) {
            continue;
        }

        for(auto subdir : *subdirectories) {
            // no way to generalize this unless start knowing how deep we are compared to pName. Which is technically possible but not convinced of the use

            currFile->cd(type+"/"+pName+"/"+subdir+"/"+subdirStack);

            TH1D* inter;
            gDirectory->GetObject(id, inter);

            if (ret == nullptr) {
                ret = std::make_shared<TH1D>(*inter);
            } else {
                ret->Add(inter);
            }
        }
    }
    ret->SetName(id + pName + "Sub");
    ret->SetTitle(id + pName + "Sub");
    
    return ret;
}
