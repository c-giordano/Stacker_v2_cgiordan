import ROOT, os, argparse
from pseudoModelConfig import * # file dictionary for easy fetching 
general_directory = '/pnfs/iihe/cms/store/user/joknolle/topproduction' # base directory (Joscha)

parser = argparse.ArgumentParser(description="Script to process ROOT files and plot histograms.")
parser.add_argument('--scaled_ratio', action='store_true', help="Draw a pad with the scaled ratio.")
parser.add_argument('--version', action='store', default='3')
parser.add_argument('--points', action='store', default='old')
parser.add_argument('--which', action='store', default='TTTT')
args = parser.parse_args()

if args.points=="old":
    files_dictionary = old_MP
elif args.points=="new":
    files_dictionary = new_MP
else:
    raise ValueError("Invalid value for --points. Use 'old' or 'new'.")

if args.which:
    files_dictionary = {k: v for k, v in files_dictionary.items() if args.which in k}


def plot_genHt_comparison(scenario_name, files_dictionary, output_dir, genHtCuts=[280],):
    """
    Function to plot genHt histograms for PseudoScalar and Scalar models, including the ratio plot with statistical uncertainties & the normalization.

    Parameters:
    - scenario_name: The key in files_dictionary for the specific scenario
    - files_dictionary: A dictionary containing the pseudo and scalar root file paths
    - output_dir: Directory to save the resulting plot
    - genHtCuts: Min and max value for genHt (default 280 and then no upper cut)
    - x_range: Tuple specifying the range of the x-axis (see config)
    - n_bins: Number of bins for the histograms (see config)
    """

    pseudo_file_path = general_directory + '/' + files_dictionary[scenario_name]["PseudoScalar"]
    normal_file_path = general_directory + '/' + files_dictionary[scenario_name]["Scalar"]
    
    
    bin_edges_array = ROOT.std.vector('double')(len(files_dictionary[scenario_name]["bin_edges"]))
    for i, edge in enumerate(files_dictionary[scenario_name]["bin_edges"]):
        bin_edges_array[i] = edge

    pseudoModel = ROOT.TFile.Open(pseudo_file_path)
    normalModel = ROOT.TFile.Open(normal_file_path)

    pseudoTree = pseudoModel.Get("Events")
    normalTree = normalModel.Get("Events")

    h_pseudo = ROOT.TH1F('genHt_pseudo', '', len(files_dictionary[scenario_name]["bin_edges"]) - 1, bin_edges_array.data())
    h_normal = ROOT.TH1F('genHt_normal', '', len(files_dictionary[scenario_name]["bin_edges"]) - 1, bin_edges_array.data())

    h_pseudo.Sumw2()
    h_normal.Sumw2()

    pseudo_lhe_weight = None
    normal_lhe_weight = None

    event_limit = 5000 
    for tree, hist in [(pseudoTree, h_pseudo), (normalTree, h_normal)]:
        event_count = 0
        max_genHt = float('-inf')
        for event in tree:
            if event_count >= event_limit:
                break
            genHt = 0
            for i in range(len(event.GenJet_pt)):
                if event.GenJet_pt[i] > 25 and abs(event.GenJet_eta[i]) < 2.4:
                    genHt += event.GenJet_pt[i]
            
            if (genHt > genHtCuts[0]):
                hist.Fill(genHt)
            if genHt > max_genHt:
                max_genHt = genHt
            event_count += 1

    print(f"Maximum genHt for current tree: {max_genHt}")

    if args.scaled_ratio:
        canvas = ROOT.TCanvas('canvas', 'canvas', 600, 1100)
        canvas.SetLeftMargin(0.2)
        pad1 = ROOT.TPad("pad1", "Main Plot", 0.1, 0.4, 1, 1)      # top
        pad2 = ROOT.TPad("pad2", "Ratio Plot", 0.1, 0.2, 1, 0.4)   # middle
        pad3 = ROOT.TPad("pad3", "Scaled Ratio Plot", 0.1, 0, 1, 0.2)  # bottom
        # pad3.SetBottomMargin(0.3)
    else:
        canvas = ROOT.TCanvas('canvas', 'canvas', 600, 800)
        canvas.SetLeftMargin(0.2)
        pad1 = ROOT.TPad("pad1", "Main Plot", 0.1, 0.3, 1, 1)
        pad2 = ROOT.TPad("pad2", "Ratio Plot", 0.1, 0, 1, 0.3)
    
    pad1.SetTopMargin(0.015)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.15)
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.24)
    pad2.SetLeftMargin(0.15)
    if args.scaled_ratio:
        pad3.SetTopMargin(0.03)
        pad3.SetBottomMargin(0.24)
        pad3.SetLeftMargin(0.15)

    pad1.Draw()
    pad2.Draw()
    if args.scaled_ratio:
        pad3.Draw()

    pad1.cd()
    h_pseudo.GetXaxis().SetRangeUser(*files_dictionary[scenario_name]["x_range"])
    h_normal.GetXaxis().SetRangeUser(*files_dictionary[scenario_name]["x_range"])
    h_pseudo.GetXaxis().SetLabelSize(0)
    h_normal.GetXaxis().SetLabelSize(0)
    h_pseudo.GetYaxis().SetTitle("Entries")
    # h_pseudo.GetYaxis().SetTitleSize(0.1)
    h_pseudo.GetYaxis().SetTitleOffset(2.4)

    ROOT.gStyle.SetOptStat(0)

    h_pseudo.SetLineColor(ROOT.kRed)
    h_pseudo.SetLineWidth(2)
    h_pseudo.SetMarkerColor(ROOT.kRed)
    h_pseudo.Draw("E SAME")
    h_pseudo.Draw("HIST SAME")

    h_normal.SetLineColor(ROOT.kBlue)
    h_normal.SetLineWidth(2)
    h_normal.SetMarkerColor(ROOT.kBlue)
    h_normal.Draw("E SAME")
    h_normal.Draw("HIST SAME")

    # legend using pretty names
    pretty_names = files_dictionary[scenario_name]["pretty_names"]
    legend = ROOT.TLegend(0.5, 0.75, 0.85, 0.95)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetMargin(0.05)
    legend.SetTextSize(0.04)
    legend.SetEntrySeparation(0.001)
    legend.AddEntry(h_pseudo, pretty_names[0], "l")
    legend.AddEntry(h_normal, pretty_names[1], "l")
    legend.AddEntry(0, pretty_names[2], "") 

    legend.Draw()

    # ratio plot part
    pad2.cd()
    h_ratio = h_pseudo.Clone('genHt_ratio')
    h_ratio.SetTitle("")
    h_ratio.Divide(h_normal)
    h_ratio.SetMarkerStyle(20)
    h_ratio.SetMarkerSize(0.8)
    h_ratio.SetLineWidth(2)
    h_ratio.SetMinimum(0)
    h_ratio.SetMaximum(2)
    h_ratio.GetYaxis().SetNdivisions(2, 5, 0)
    h_ratio.GetYaxis().SetTickLength(0.03)
    h_ratio.GetYaxis().SetLabelSize(0.1)
    h_ratio.GetXaxis().SetRangeUser(*files_dictionary[scenario_name]["x_range"])
    h_ratio.SetLineColor(ROOT.kBlack)
    h_ratio.SetLineWidth(2)
    h_ratio.SetMarkerColor(ROOT.kBlack)
    h_ratio.Draw("E")
    h_ratio.Draw("HIST SAME")

    if args.scaled_ratio:
        # hides x-axis labels for the second pad if the third pad exists (contrite but works)
        h_ratio.GetXaxis().SetLabelSize(0)
        h_ratio.GetXaxis().SetTitleSize(0)
    else:
        # shows x-axis labels on the second pad if there is no third pad
        h_ratio.GetXaxis().SetTitle("Gen H_{T}")
        h_ratio.GetXaxis().SetTitleSize(0.1)
        h_ratio.GetXaxis().SetLabelSize(0.1)

    line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1, h_ratio.GetXaxis().GetXmax(), 1)
    line.SetLineColor(ROOT.kGray)
    line.SetLineStyle(2)
    line.Draw("SAME")

    h_ratio.GetYaxis().SetTitle("Ratio")
    h_ratio.GetYaxis().SetTitleSize(0.1)
    h_ratio.GetYaxis().SetTitleOffset(0.5)

    if args.scaled_ratio:
        pad3.cd()
        xSec_ratio = files_dictionary[scenario_name]["xSec_ratio"]
        h_scaled_ratio = h_ratio.Clone("h_scaled_ratio")
        h_scaled_ratio.Scale(xSec_ratio)  # scale by xSec_ratio

        h_scaled_ratio.SetLineColor(ROOT.kGreen+3)
        h_scaled_ratio.SetMarkerColor(ROOT.kGreen+3)
        h_scaled_ratio.SetTitle("")
        h_scaled_ratio.SetMarkerStyle(20)
        h_scaled_ratio.SetMarkerSize(0.8)
        h_scaled_ratio.SetLineWidth(2)
        h_scaled_ratio.SetMinimum(files_dictionary[scenario_name]["ratio_range"][0])
        h_scaled_ratio.SetMaximum(files_dictionary[scenario_name]["ratio_range"][1])
        h_scaled_ratio.GetXaxis().SetRangeUser(*files_dictionary[scenario_name]["x_range"])
        h_scaled_ratio.GetXaxis().SetTitle("Gen H_{T}")
        h_scaled_ratio.GetXaxis().SetTitleSize(0.1)
        h_scaled_ratio.GetXaxis().SetLabelSize(0.1)
        h_scaled_ratio.Draw("E")
        h_scaled_ratio.Draw("HIST SAME")

        line2 = ROOT.TLine(h_scaled_ratio.GetXaxis().GetXmin(), 1, h_scaled_ratio.GetXaxis().GetXmax(), 1)
        line2.SetLineColor(ROOT.kGray)
        line2.SetLineStyle(2)
        line2.Draw("SAME")

        h_scaled_ratio.GetXaxis().SetTitle("Gen H_{T}")
        h_scaled_ratio.GetXaxis().SetTitleSize(0.1)
        h_scaled_ratio.GetYaxis().SetTitle("Ratio x #frac{#sigma(PS)}{#sigma(S)}")
        h_scaled_ratio.GetYaxis().SetTitleSize(0.1)

    ratio_file = ROOT.TFile(f"{output_dir}/{scenario_name}_ratio.root", "RECREATE")
    h_ratio.Write()
    ratio_file.Close()
    if args.scaled_ratio:
        ratio_scaled = ROOT.TFile(f"{output_dir}/{scenario_name}_ratio_scaled.root", "RECREATE")
        h_scaled_ratio.Write()
        ratio_scaled.Close()

    canvas.SaveAs(f"{output_dir}/{scenario_name}.png")
    canvas.SaveAs(f"{output_dir}/{scenario_name}.pdf")
    canvas.SaveAs(f"{output_dir}/{scenario_name}.C")
    canvas.SaveAs(f"{output_dir}/{scenario_name}.root")


output_directory = os.path.expanduser('~/public_html/GenPseudo/finalPlots_v'+args.version)
if not os.path.exists(output_directory):
    os.makedirs(output_directory)


for scenario in files_dictionary:
    plot_genHt_comparison(scenario, files_dictionary, output_directory)


os.system(f"cp /user/nivanden/public_html/index.php {output_directory}/index.php")



