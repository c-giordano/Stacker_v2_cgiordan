#!/usr/bin/env python3
# this works in CMSSW_13_3_3
# OG directory is /ada_mnt/ada/user/cgiordan/CMSSW_13_3_3/bin/slc7_amd64_gcc12/plotLimits.py

from __future__ import absolute_import
from __future__ import print_function
import ROOT, os, datetime
import HiggsAnalysis.CombinedLimit.util.plotting as plot
import argparse
import CombineHarvester.CombineTools.maketable as maketable
from array import array
import json, numbers, sys

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs="+", help="""Input json files""")
parser.add_argument("--output", "-o", default="limit2018", help="""Name of the output plot without file extension""")
parser.add_argument("--show", default="exp,obs")
# parser.add_argument(
#     '--debug-output', '-d', help="""If specified, write the
#     TGraphs into this output ROOT file""")
parser.add_argument("--x-title", default="", help="""Title for the x-axis""")
parser.add_argument("--y-title", default=None, help="""Title for the y-axis""")
parser.add_argument("--limit-on", default="g", help="""Shortcut for setting the y-axis label""")
parser.add_argument("--cms-sub", default="Internal", help="""Text below the CMS logo""")
parser.add_argument(
    "--scenario-label",
    default="",
    help="""Scenario name to be drawn in top
    left of plot""",
)
parser.add_argument("--era", default="", help="Era of the analysis")
parser.add_argument("--title-right", default="", help="""Right header text above the frame""")
parser.add_argument("--title-left", default="", help="""Left header text above the frame""")
parser.add_argument("--logy", action="store_true", help="""Draw y-axis in log scale""")
parser.add_argument("--logx", action="store_true", help="""Draw x-axis in log scale""")
parser.add_argument("--ratio-to", default=None)
parser.add_argument("--pad-style", default=None, help="""Extra style options for the pad, e.g. Grid=(1,1)""")
parser.add_argument("--auto-style", nargs="?", const="", default=None, help="""Take line colors and styles from a pre-defined list""")
parser.add_argument("--table_vals", help="Amount of values to be written in a table for different masses", default=10)
parser.add_argument('--output-dir', '-d', default=os.path.expanduser('~/public_html/Limits/Setup_2025'), help="Base directory to save the output files")
parser.add_argument('--model', '-m', default='VectorOctet')
parser.add_argument('--suffix', default="")
parser.add_argument('--gamma-json', help='File JSON con il limite su Γ da plottare', default=None)
parser.add_argument('--gamma-style', default='LineColor=ROOT.kOrange+7,LineStyle=1,LineWidth=3',help='ROOT style (es. "LineColor=2,LineWidth=4")')
parser.add_argument('--gamma-key',default='gamma_10')
parser.add_argument('--onlyTTTT', action='store_true', help='Plot obs line from TTTT-only scenario')
# parser.add_argument('--fromJson', action='store_true', help='Import and plot from JSON (default behavior)')
# parser.add_argument('--fromRoot', help='ROOT file with TGraph to plot as gamma curve')
args = parser.parse_args()



if args.era == "2016":
    args.title_right = "36.3 fb^{-1} (13 TeV)"
elif args.era == "2017":
    args.title_right = "41.5 fb^{-1} (13 TeV)"
elif args.era == "2018":
    args.title_right = "59.8 fb^{-1} (13 TeV)"
elif args.era == "Run2":
    args.title_right = "138 fb^{-1} (13 TeV)"

base_gamma_dir = "/pnfs/iihe/cms/store/user/cgiordan/jsons/"
if args.gamma_json:
    gamma_file = os.path.join(base_gamma_dir, args.gamma_json)

model_titles = {
    "VectorOctet": ("m_{V_{8}} (GeV)", r"g_{8L/R}"),
    "VectorSinglet": ("m_{V_{1}} (GeV)", r"g_{1L/R}"),
    "ScalarOctet": ("m_{S_{8}} (GeV)", r"y_{8S}"),
    "ScalarSinglet": ("m_{S_{1}} (GeV)", r"y_{1S}"),
    "PseudoScalarOctet": ("m_{P_{8}} (GeV)", r"y_{8P}"),
    "PseudoScalarSinglet": ("m_{P_{1}} (GeV)", r"y_{1P}")
}

if args.model in model_titles:
    args.x_title, args.limit_on = model_titles[args.model]


def copy_index_html(folder):
    if os.getenv("CMSSW_VERSION"):
        outputfile = os.path.join(folder, "index.php")
        if os.path.exists(outputfile):
            return
        os.system(f"cp /user/nivanden/public_html/index.php {outputfile}")

output_dir = os.path.join(args.output_dir, args.model+"_"+args.era+args.suffix)


def graph_from_gamma_json(path, inner_key='gamma_10'):
    """
    Legge un JSON {mass: limit, ...} e restituisce un TGraph ordinato in massa.
    """
    with open(path) as f:
        data = json.load(f)

    # supporta chiavi "1000" o 1000
    masses = sorted(float(m) for m in data.keys())
    limits = []
    for m in masses:
        entry = data[str(m)]
        if not isinstance(entry, dict):
            raise TypeError(f"Voce non-dict per m={m}: {entry}")

        if inner_key not in entry or not isinstance(entry[inner_key], numbers.Real):
            raise ValueError(
                f"Per m={m} non trovo un valore numerico nella chiave '{inner_key}'. "
                f"Chiavi disponibili: {list(entry.keys())}"
            )
        limits.append(float(entry[inner_key]))

    # limits = [float(data[str(int(m))] if str(int(m)) in data else data[str(m)]) for m in masses]

    return ROOT.TGraph(len(masses), array('d', masses), array('d', limits))


# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

copy_index_html(output_dir)

def DrawAxisHists(pads, axis_hists, def_pad=None):
    for i, pad in enumerate(pads):
        pad.cd()
        axis_hists[i].Draw("AXIS")
        axis_hists[i].Draw("AXIGSAME")
    if def_pad is not None:
        def_pad.cd()


## Boilerplate
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
plot.ModTDRStyle()
ROOT.gStyle.SetNdivisions(510, "XYZ")  # probably looks better

canv = ROOT.TCanvas(args.output, args.output)


if args.ratio_to is not None:
    pads = plot.TwoPadSplit(0.30, 0.01, 0.01)
else:
    pads = plot.OnePad()

# Set the style options of the pads
for padx in pads:
    # Use tick marks on oppsite axis edges
    plot.Set(padx, Tickx=1, Ticky=1, Logx=args.logx)
    if args.pad_style is not None:
        settings = {x.split("=")[0]: eval(x.split("=")[1]) for x in args.pad_style.split(",")}
        print("Applying style options to the TPad(s):")
        print(settings)
        plot.Set(padx, **settings)

graphs = []
graph_sets = []

legend = plot.PositionedLegend(0.45, 0.10, 3, 0.015)
plot.Set(legend, NColumns=2)

axis = None

defcols = [
    ROOT.kGreen + 3,
    # '#607641',
    ROOT.kRed,
    ROOT.kBlue,
    ROOT.kBlack,
    # '#F5BB54',
    ROOT.kYellow + 2,
    ROOT.kOrange + 10,
    ROOT.kCyan + 3,
    ROOT.kMagenta + 2,
    ROOT.kViolet - 5,
    ROOT.kGray,
]

deflines = [1, 2, 3]

if args.auto_style is not None:
    icol = {x: 0 for x in args.auto_style.split(",")}
    icol["default"] = 0
    iline = {}
    iline["default"] = 1
    for i, x in enumerate(args.auto_style.split(",")):
        iline[x] = i + 1

# Process each input argument
for src in args.input:
    splitsrc = src.split(":")
    file = splitsrc[0]
    # limit.json => Draw as full obs + exp limit band
    if len(splitsrc) == 1:
        graph_sets.append(plot.StandardLimitsFromJSONFile(file, args.show.split(",")))
        print("Debug: Contents of graph_sets[-1]:", graph_sets[-1])
        if axis is None:
            axis = plot.CreateAxisHists(len(pads), list(graph_sets[-1].values())[0], True)
            DrawAxisHists(pads, axis, pads[0])
        plot.StyleLimitBand(graph_sets[-1])
        plot.DrawLimitBand(pads[0], graph_sets[-1], legend=legend)
        pads[0].RedrawAxis()
        pads[0].RedrawAxis("g")
        pads[0].GetFrame().Draw()

    # limit.json:X => Draw a single graph for entry X in the json file
    # 'limit.json:X:Title="Blah",LineColor=4,...' =>
    # as before but also apply style options to TGraph
    elif len(splitsrc) >= 2:
        settings = {}
        settings["Title"] = src
        if args.auto_style is not None:
            nm = "default"
            for x in icol.keys():
                if x in splitsrc[1]:
                    nm = x
            i = icol[nm]  # take the next default color...
            j = iline[nm]  # take the next default line style...
            settings["LineColor"] = defcols[i]
            settings["MarkerColor"] = defcols[i]
            settings["LineStyle"] = j
            icol[nm] = (i + 1) if (i + 1) < len(defcols) else 0
        graphs.append(plot.LimitTGraphFromJSONFile(file, splitsrc[1]))
        if len(splitsrc) >= 3:
            settings.update({x.split("=")[0]: eval(x.split("=")[1]) for x in splitsrc[2].split(",")})
        plot.Set(graphs[-1], **settings)
        if axis is None:
            axis = plot.CreateAxisHists(len(pads), graphs[-1], True)
            DrawAxisHists(pads, axis, pads[0])
        graphs[-1].Draw("PLSAME")
        legend.AddEntry(graphs[-1], "", "PL")


axis[0].GetYaxis().SetTitle("95%% CL limit on %s" % args.limit_on)
if args.y_title is not None:
    axis[0].GetYaxis().SetTitle(args.y_title)
axis[0].GetXaxis().SetTitle(args.x_title)
axis[0].GetXaxis().SetLabelOffset(axis[0].GetXaxis().GetLabelOffset() * 2)

if args.logy:
    axis[0].SetMinimum(0.1)  # we'll fix this later
    pads[0].SetLogy(True)
    # axis[0].GetYaxis().SetMoreLogLabels()
    # axis[0].SetNdivisions(50005, "X")

y_min, y_max = (plot.GetPadYMin(pads[0]), plot.GetPadYMax(pads[0]))
plot.FixBothRanges(pads[0], y_min if args.logy else 0, 0.05 if args.logy else 0, y_max, 0.25)

ratio_graph_sets = []
ratio_graphs = []



if args.ratio_to is not None:
    pads[1].cd()
    plot.SetupTwoPadSplitAsRatio(pads, axis[0], axis[1], "Ratio_{}", True, 0.1, 2.4)
    axis[1].SetNdivisions(506, "Y")
    splitsrc = args.ratio_to.split(":")
    ref = plot.LimitTGraphFromJSONFile(splitsrc[0], splitsrc[1])
    for gr_set in graph_sets:
        ratio_set = {}
        for key in gr_set:
            ratio_set[key] = plot.GraphDivide(gr_set[key], ref)
        ratio_graph_sets.append(ratio_set)
        plot.DrawLimitBand(pads[1], ratio_graph_sets[-1])
        pads[1].RedrawAxis()
        pads[1].RedrawAxis("g")
        pads[1].GetFrame().Draw()
    for gr in graphs:
        ratio_graphs.append(plot.GraphDivide(gr, ref))
        ratio_graphs[-1].Draw("LP")
    ry_min, ry_max = (plot.GetPadYMin(pads[1]), plot.GetPadYMax(pads[1]))
    plot.FixBothRanges(pads[1], ry_min, 0.1, ry_max, 0.1)


if args.gamma_json:
    

    gr_gamma = graph_from_gamma_json(gamma_file, inner_key=args.gamma_key)

    # default_style = "LineColor=ROOT.kOrange+7,LineStyle=2,LineWidth=3"
    # style_str = args.gamma_style if hasattr(args, 'gamma_style') else default_style
    # style_dict = {kv.split('=')[0]: eval(kv.split('=')[1]) for kv in style_str.split(',')}
    # plot.Set(gr_gamma, Title='', **style_dict)

    # # Disegna sopra il pad principale
    if args.model != "PseudoScalarSinglet":
        pad_main = pads[0]
        pad_main.cd()
        pads[0].cd()
        else_model_ymax = pad_main.GetUymax()
        npts = gr_gamma.GetN()
        xs = [gr_gamma.GetX()[i] for i in range(npts)]
        ys = [gr_gamma.GetY()[i] for i in range(npts)]

        # Costruisci contorno verso l’alto
        xs_shade = [xs[0]] + xs + xs[::-1] + [xs[0]]
        ys_shade = [else_model_ymax] + ys + [else_model_ymax] * len(xs) + [else_model_ymax]

        # TGraph shading
        gr_shade = ROOT.TGraph(len(xs_shade))
        for i in range(len(xs_shade)):
            gr_shade.SetPoint(i, xs_shade[i], ys_shade[i])

        gr_shade.SetFillColorAlpha(ROOT.kGray, 0.2)
        gr_shade.SetLineWidth(0)
        pad_main.cd()
        gr_shade.Draw("F SAME")

        # Curva sopra
        style_dict = {kv.split('=')[0]: eval(kv.split('=')[1]) for kv in args.gamma_style.split(',')}
        plot.Set(gr_gamma, Title='', **style_dict)
        gr_gamma.Draw("L SAME")

        # Legenda
        legend.AddEntry(gr_gamma, '#Gamma(t#bar{t}) > 10 GeV', 'L')


    else:
        pads[0].cd()
        axis[0].SetMaximum(1.4)
        ymax = 1.4

        npts = gr_gamma.GetN()
        xs = [gr_gamma.GetX()[i] for i in range(npts)]
        ys = [gr_gamma.GetY()[i] for i in range(npts)]


        xs_shade = [xs[0]] + xs + xs[::-1] + [xs[0]]
        ys_shade = [ymax] + ys + [ymax] * len(xs) + [ymax]


        gr_shade = ROOT.TGraph(len(xs_shade))
        for i in range(len(xs_shade)):
            gr_shade.SetPoint(i, xs_shade[i], ys_shade[i])

        gr_shade.SetFillColorAlpha(ROOT.kGray, 0.2)
        gr_shade.SetLineWidth(0)
        pads[0].cd()
        gr_shade.Draw("F SAME")

        plot.Set(gr_gamma, Title='', LineColor=ROOT.kGray, LineStyle=1, LineWidth=3)
        gr_gamma.Draw("L SAME")

        legend.AddEntry(gr_gamma, '#Gamma(t#bar{t}) > 10 GeV', 'L')



    # --- ratio-plot -----------------------------------------------------------
    if args.ratio_to is not None:
        pads[1].cd()
        ratio_gamma = plot.GraphDivide(gr_gamma, ref)
        plot.Set(ratio_gamma, **style_dict)
        ratio_gamma.Draw('L SAME')

if args.onlyTTTT:
    tttt_dir = f"../TopPhilic{args.model}_data_newUnc_onlyTTTT"
    tttt_file = os.path.join(tttt_dir, "total_default.json")
    if os.path.exists(tttt_file):
        print(f"Adding obs + exp0 + band from {tttt_file}")

        with open(tttt_file) as f:
            data = json.load(f)

        masses = sorted(float(m) for m in data.keys())
        obs = []
        exp0 = []

        for m in masses:
            entry = data[str(m)]
            obs.append(entry["obs"])
            exp0.append(entry["exp0"])

        # --- 1. Banda tra exp0 e obs ---
        xs_shade = masses + masses[::-1]
        ys_shade = obs + exp0[::-1]

        gr_band = ROOT.TGraph(len(xs_shade))
        for i in range(len(xs_shade)):
            gr_band.SetPoint(i, xs_shade[i], ys_shade[i])
        gr_band.SetFillStyle(3004)  # Tratteggio orizzontale/diagonale
        gr_band.SetFillColor(ROOT.kGreen + 1)
        gr_band.SetLineColor(0)
        gr_band.SetLineWidth(0)

        # --- 2. Curva obs ---
        gr_obs = ROOT.TGraph(len(masses), array('d', masses), array('d', obs))
        plot.Set(gr_obs,
                LineColor=ROOT.kCyan,
                LineStyle=1,        # linea continua
                LineWidth=3,
                MarkerStyle=20,     # punti visibili
                MarkerSize=1.0,
                MarkerColor=ROOT.kCyan)

        # --- 3. Punti exp0 ---
        gr_exp0 = ROOT.TGraph(len(masses), array('d', masses), array('d', exp0))
        plot.Set(gr_exp0,
                LineColor=ROOT.kOrange+7,
                LineStyle=2,        # tratteggiata
                LineWidth=2,
                MarkerStyle = 0,
                MarkerColor=ROOT.kOrange+7)
        gr_exp0.SetMarkerStyle(0) 
        gr_exp0.Draw("L SAME")

        # --- 4. Disegno nel giusto ordine ---
        pads[0].cd()
        gr_band.Draw("F SAME")     # prima la banda
        gr_obs.Draw("L SAME")      # poi la curva
        gr_exp0.Draw("L SAME")     # infine i punti
        legend.AddEntry(gr_obs, "Observed (only #it{t}#bar{#it{t}}#it{t}#bar{#it{t}})", "L")
        legend.AddEntry(gr_exp0, "Expected (only #it{t}#bar{#it{t}}#it{t}#bar{#it{t}})", "L")

    else:
        print(f"WARNING: File not found: {tttt_file}")



pads[0].cd()
if legend.GetNRows() == 1:
    legend.SetY1(legend.GetY2() - 0.5 * (legend.GetY2() - legend.GetY1()))
legend.Draw()

# line = ROOT.TLine()
# line.SetLineColor(ROOT.kBlue)
# line.SetLineWidth(2)
# plot.DrawHorizontalLine(pads[0], line, 1)

box = ROOT.TPave(pads[0].GetLeftMargin(), 0.81, 1 - pads[0].GetRightMargin(), 1 - pads[0].GetTopMargin(), 1, "NDC")
box.Draw()

legend.Draw()

plot.DrawCMSLogo(pads[0], "CMS", args.cms_sub, 11, 0.095, 0.035, 1.4, "", 0.8)
plot.DrawTitle(pads[0], args.title_right, 3)
plot.DrawTitle(pads[0], args.title_left, 1)

canv.Print(".pdf")
canv.Print(".png")

output_pdf = os.path.join(output_dir, f"{args.output}.pdf")
output_png = os.path.join(output_dir, f"{args.output}.png")
canv.Print(output_pdf)
canv.Print(output_png)
# canv.Print(output_pdf)
# canv.Print(output_png)
# maketable.TablefromJson(args.input[0], "TablefromJson.txt")
