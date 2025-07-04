import ROOT
import json
import argparse
import os, shutil

parser = argparse.ArgumentParser(description="Script to process ROOT files and plot histograms.")
parser.add_argument('--era', action='store', help="which era", default="2016")
parser.add_argument('--data', action='store_true', default=False, help='do you want to plot data?')
parser.add_argument('--version', action='store', default='1', help='which version is this?')
parser.add_argument('--fileDir', action='store', default='./output/bsmlimits/TopPhilicScalarOctet_data/600/', help="where the DC at")
parser.add_argument('--model', action='store', default='')
parser.add_argument('--mass', action='store', default='')
parser.add_argument('--BSM', action='store_true', default=False)

args = parser.parse_args()

def output_directory_setup(base_dir: str, version: str, era: str, index_file: str = "/user/nivanden/public_html/index.php"):
    if args.model!='' and args.mass!='':
        output_directory = os.path.join(os.path.expanduser(base_dir), "v" + version, args.model, args.mass, era)
    else:
        output_directory = os.path.join(os.path.expanduser(base_dir), "v" + version, era)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    os.system(f"cp {index_file} {output_directory}/index.php")

    return output_directory

output_dir = os.path.join("~/public_html/PlotsFromDC")
output_directory = output_directory_setup(output_dir,args.version, args.era)

region = {  "SR2L_ee_Sig": "H_{T}(GeV)",
            "SR2L_em_Sig": "H_{T}(GeV)",
            "SR2L_mm_Sig": "H_{T}(GeV)",
            "SR2L_ee_ttw": "BDT_TTX",
            "SR2L_em_ttw": "BDT_TTX",
            "SR2L_mm_ttw": "BDT_TTX",
            # "SR2L_Sig":     "H_{T}(GeV)",
            # "SR2L_ttw":     "BDT_TTX",
            "SR2L_NP":     "BDT_TT",
            "SR3L_Sig":    "H_{T}(GeV)",
            "SR3L_NP":     "BDT_TT",
            "SR3L_ttw":    "BDT_TTX",
            "SR4L_Sig":    "BDT_TTX",
            # "SR4L_ttw":    "BDT_TTX",
            "CR2LNP":      "BDT_TT",
            "CR2LTTW":     "BDT_TT",
            "CR3LNP":      "Sum Lep Charges",
            "CR4LZ":       "N_{jets}(GeV)",
            "CR3LZ":        "N_{jets}(GeV)"
             }
          
    
exclude_histograms = ["BSM_Quad", "BSM_Quartic"] # list of histograms to exclude
if not args.data:
    exclude_histograms += ['data_obs']

# colours from the json setting files (tttt looks beige but ok...)
process_colors = {
        "ttt": "#e9d98d", "Othert": "#92dadd", "Xg": "#3f90da", "ChargeMisID": "#832db6",
        "VVV": "#bd1f01", "WZ": "#bd1f01", "nonPromptElectron": "#e76300", "nonPromptMuon": "#a96b59",
        "ttH": "#b9ac70", "ttZ": "#717581", "ttW": "#94a4a2", "tttt": "#ffa90e"
    }


def plot_stacked_histograms(region_name, exclude_histograms, process_colors):
    file = ROOT.TFile(os.path.join(args.fileDir, "DC_"+args.era+".root "), "READ")
    if not file or file.IsZombie():
        print("Error opening file! Check your DC")
        return

    region_dir = file.Get(region_name)
    if not region_dir:
        print(f"Region {region_name} not found in the file!")
        return

    canvas = ROOT.TCanvas("canvas", region_name, 800, 800)
    canvas.cd()
    stack = ROOT.THStack("stack", region_name)

    legend = ROOT.TLegend(0.45, 0.45, 0.88, 0.88) 
    legend.SetTextSize(0.04)
    legend.SetLineColor(0)
    legend.SetFillStyle(0)

    histograms = []
    max_hist_value = 0
    margin = 0.3  # 10% margin above the highest value

    for key in region_dir.GetListOfKeys():
        hist = key.ReadObj()
        
        # skip histograms in the exclude list
        if hist.GetName() in exclude_histograms:
            continue

        if isinstance(hist, ROOT.TH1):
            # check the process name from the histogram name
            for process, color in process_colors.items():
                if process in hist.GetName():
                    
                    hist.SetLineColor(ROOT.TColor.GetColor(color))
                    hist.SetFillColor(ROOT.TColor.GetColor(color)) 
                    hist.SetFillStyle(1001)

                    histograms.append(hist)
                    stack.Add(hist) 
                    legend.AddEntry(hist, process, "f")

                    max_hist_value = max(max_hist_value, hist.GetMaximum())
                    break

    if histograms:
        stack.Draw("HIST")

    if args.data:
        data_hist = region_dir.Get("data_obs")
        if data_hist:
            data_hist.SetMarkerStyle(20)  
            data_hist.SetMarkerSize(1.2)  
            data_hist.SetMarkerColor(ROOT.kBlack)  
            data_hist.Draw("SAME E") 
            legend.AddEntry(data_hist, "Data", "p")

            max_hist_value = max(max_hist_value, data_hist.GetMaximum())
        
        stack.SetMaximum(max_hist_value * (1 + margin))  
        stack.SetMinimum(0)

    legend.Draw()  


    tex = ROOT.TLatex()
    tex.SetNDC(True)  
    tex.SetTextSize(0.04)

    x_axis_label = region.get(region_name, "X-axis Label")
    tex.DrawLatex(0.8, 0.03, x_axis_label )

    tex.SetTextAngle(90)
    tex.DrawLatex(0.05, 0.8, "Events")
    canvas.Update()
    
    canvas.SaveAs(f"{output_directory}/{region_name}_shape.png")
    canvas.SaveAs(f"{output_directory}/{region_name}_shape.pdf")

# def plot_bsm_quartic(region_name):
#     file = ROOT.TFile(os.path.join(args.fileDir, "DC_" + args.era + ".root"), "READ")
#     if not file or file.IsZombie():
#         print("Error opening file!")
#         return

#     region_dir = file.Get(region_name)
#     if not region_dir:
#         print(f"Region {region_name} not found!")
#         return

#     bsm_hist = region_dir.Get("BSM_Quartic")
#     if bsm_hist:
#         bsm_canvas = ROOT.TCanvas(region_name + "_BSM_variation", region_name + "_BSM_variation", 800, 800)
#         bsm_canvas.cd()
#         bsm_hist.SetLineColor(ROOT.kBlue)
#         bsm_hist.SetLineWidth(2)
#         bsm_hist.SetStats(0)
#         # bsm_hist.Draw("HIST")
#         bsm_hist.Draw("HIST")
#         bsm_hist.GetXaxis().SetTitle(region.get(region_name, "X-axis Label"))  # Use the same axis label
#         bsm_hist.GetYaxis().SetTitle("Events")
#         bsm_dir = os.path.join(output_directory, "BSM")
#         sys_file = "/user/nivanden/public_html/index.php"
#         if not os.path.exists(bsm_dir):
#             os.makedirs(os.path.join(output_directory, "BSM"))
#         shutil.copy(sys_file, os.path.join(bsm_dir, "index.php"))
#         bsm_canvas.SaveAs(f"{output_directory}/BSM/{region_name}_BSM_Quartic.png")
#         bsm_canvas.SaveAs(f"{output_directory}/BSM/{region_name}_BSM_Quartic.pdf")

def plot_bsm_quartic(region_name):
    file = ROOT.TFile(os.path.join(args.fileDir, "DC_" + args.era + ".root"), "READ")
    if not file or file.IsZombie():
        print("Error opening file!")
        return

    region_dir = file.Get(region_name)
    if not region_dir:
        print(f"Region {region_name} not found!")
        return

    bsm_quartic_hist = region_dir.Get("BSM_Quartic")
    bsm_quad_hist = region_dir.Get("BSM_Quad")
    
    if not bsm_quartic_hist or not bsm_quad_hist:
        print(f"Histograms not found in {region_name}!")
        return
    
    # Create canvas
    bsm_canvas = ROOT.TCanvas(region_name + "_BSM_variation", region_name + "_BSM_variation", 800, 800)
    bsm_canvas.cd()
    
    # Styling histograms
    bsm_quartic_hist.SetLineColor(ROOT.kBlue)
    bsm_quartic_hist.SetLineWidth(2)
    bsm_quartic_hist.SetStats(0)
    
    bsm_quad_hist.SetLineColor(ROOT.kRed)
    bsm_quad_hist.SetLineWidth(2)
    bsm_quad_hist.SetStats(0)
    
    # Draw histograms
    bsm_quartic_hist.Draw("HIST")
    bsm_quad_hist.Draw("HIST SAME")
    
    # Set axis labels
    bsm_quartic_hist.GetXaxis().SetTitle(region.get(region_name, "X-axis Label"))
    bsm_quartic_hist.GetYaxis().SetTitle("Events")
    
    # Create legend
    legend = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
    legend.AddEntry(bsm_quartic_hist, "BSM Quartic", "l")
    legend.AddEntry(bsm_quad_hist, "BSM Quad", "l")
    legend.Draw()
    
    # Save plots
    bsm_dir = os.path.join(output_directory, "BSM")
    sys_file = "/user/nivanden/public_html/index.php"
    if not os.path.exists(bsm_dir):
        os.makedirs(bsm_dir)
    shutil.copy(sys_file, os.path.join(bsm_dir, "index.php"))
    
    bsm_canvas.SaveAs(f"{output_directory}/BSM/{region_name}_BSM_Comparison.png")
    bsm_canvas.SaveAs(f"{output_directory}/BSM/{region_name}_BSM_Comparison.pdf")


for r in region.keys():
    plot_stacked_histograms(r, exclude_histograms, process_colors)

if args.BSM:
    for r in region.keys():
        plot_bsm_quartic(r)
