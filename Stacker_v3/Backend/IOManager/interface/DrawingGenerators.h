#include <memory>
#include <TCanvas.h>
#include <TPad.h>
#include <TString.h>

namespace DrawingGenerator
{
    // generic drawing
    std::shared_ptr<TCanvas> GenerateCanvas(TString& name);
    std::shared_ptr<TCanvas> GenerateCanvas(TString& name, double width, double height);

    std::shared_ptr<TPad> GenerateTPad(TString& name);
    std::shared_ptr<TPad> GenerateTPad(TString& name, double x1, double x2, double y1, double y2);
    
    // CMS Label
    void CMSLabel(TString& additionalText);

    // Lumi info
    void LumiInfo(double lumi);
    void LumiInfo(TString& lumi);

    // 
} // namespace DrawingGenerator
