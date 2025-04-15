#!/usr/bin/env python

import CombineHarvester.CombineTools.plotting as plot
import os, sys, json
import ROOT
import math
import argparse
from array import array

ROOT.gROOT.SetBatch(ROOT.kTRUE)
parser = argparse.ArgumentParser()
parser.add_argument('files', default='input.root', nargs="+", help='Input files')
parser.add_argument('--output', '-o', default='limit', help="""Name of the output plot without file extension""")
parser.add_argument('--sm-exp', default=[], nargs="+", help="""Input files for the SM expectation""")
parser.add_argument('--bg-exp', default=[], nargs="+", help="""Input files for the backgroud only expectation""")
parser.add_argument('--cms-sub', default='Internal', help="""Text below the CMS logo""")
parser.add_argument('--mass', default='', help="""Mass label on the plot""")
parser.add_argument('--title-right', default='', help="""Right header text above the frame""")
parser.add_argument('--title-left', default='', help="""Left header text above the frame""")
parser.add_argument('--x', default='at', help="""Variable on x axis""")
parser.add_argument('--y', default='bt', help="""Variable on y axis""")
parser.add_argument('--bestfitpoint', '-bf', default=None, help="""Best fit point""")
parser.add_argument('--debug_output', '-d', help="""If specified, write the contour TH2s and TGraphs into this output ROOT file""")
args = parser.parse_args()

def SetDeepSeaPalette():
    nRGBs = 5
    stops = array('d', [0.0000, 0.1000, 0.2000, 0.3000, 1.0000])
    red = array('d', [8./255.,  37./255., 72./255., 115./255., 171./255.])
    green = array('d', [52./255., 76./255.,  104./255.,  136./255., 180./255.])
    blue = array('d', [85./255., 106./255., 128./255., 152./255., 186./255.])
    ROOT.TColor.CreateGradientColorTable(nRGBs, stops, red, green, blue, 255, 1)        

def contourFromTH2(h2in, threshold, minPoints=10, frameValue=1000.):

    contoursList = [threshold]
    contours = array('d', contoursList)

#    h2 = frameTH2D(h2in, threshold, frameValue)
    h2 = h2in.Clone()
    h2.SetContour(1, contours)

    # Draw contours as filled regions, and Save points
    # backup = R.gPad # doesn't work in pyroot, backup behaves like a ref to gPad
    canv = ROOT.TCanvas('tmp', 'tmp')
    canv.cd()
    h2.Draw('CONT Z LIST')
    ROOT.gPad.Update()  # Needed to force the plotting and retrieve the contours in

    conts = ROOT.gROOT.GetListOfSpecials().FindObject('contours')
    contLevel = None

    if conts is None or conts.GetSize() == 0:
        print('*** No Contours Were Extracted!')
        return None
    ret = ROOT.TList()
    for i in range(conts.GetSize()):
        contLevel = conts.At(i)
        print('>> Contour %d has %d Graphs' % (i, contLevel.GetSize()))
        for j in range(contLevel.GetSize()):
            gr1 = contLevel.At(j)
            print('\t Graph %d has %d points' % (j, gr1.GetN()))
            if gr1.GetN() > minPoints:
                ret.Add(gr1.Clone())
            # // break;
    # backup.cd()
    canv.Close()
    return ret

xvar = args.x.split('_')[-1]
yvar = args.y.split('_')[-1]
vtitle = ['', '']
for iv, v in enumerate([xvar, yvar]):
    print(v)
    if v == 'TTTT': vvar = 't#bar{t}t#bar{t}'
    elif v == 'TTT': vvar = 'ttt'
    elif v == 'TTH': vvar = 't#bar{t}H'
    elif v == 'TTW': vvar = 't#bar{t}W'
    elif v == 'TTZ': vvar = 't#bar{t}Z'
    elif v == 'cpeven': vvar = '#kappa_{t}'
    elif v == 'at': vvar = '#kappa_{Htt}'
    elif v == 'cpodd': vvar = '#kappa_{t}'
    elif v == 'bt': vvar = '#kappa_{Att}'
    vtitle[iv] = '#mu_{'+vvar+'}'

plot.ModTDRStyle(width=1000, height=1000, l=0.12)
ROOT.gStyle.SetNdivisions(510, 'XYZ')
ROOT.gStyle.SetTitleYOffset(1.5)
SetDeepSeaPalette()
canv = ROOT.TCanvas(args.output, args.output)
pads = plot.OnePad()

if args.debug_output is not None:
    print("Create debug output file %s" % args.debug_output)
    debug = ROOT.TFile(args.debug_output, 'RECREATE')
else:
    debug = None

procx = args.x.replace('r_', '')
procy = args.y.replace('r_', '')
xsec = {'TTTT': 13.4, 'TTT': 2.0, 'TTW': 722, 'TTZ': 859, 'TTH': 504, 'yuk_cpeven': 1., 'yuk_cpodd': 0., 'at': 1., 'bt': 1.}

limit = plot.MakeTChain(args.files, 'limit')
graph = plot.TGraph2DFromTree(limit, args.x, args.y, '2*deltaNLL', 'quantileExpected > -0.5 && deltaNLL > 0 && deltaNLL < 1000')
# this is specific to changing signal strength, not relevant for us
# for ip in range(graph.GetN()):
#     graph.SetPoint(ip, graph.GetX()[ip]*xsec[procx], graph.GetY()[ip]*xsec[procy], graph.GetZ()[ip])
best = plot.TGraphFromTree(limit, args.x, args.y, 'deltaNLL == 0')
plot.RemoveGraphXDuplicates(best)
hists = plot.TH2FromTGraph2D(graph, method='BinCenterAligned')
plot.fastFillTH2(hists, graph, interpolateMissing=True)
if args.bg_exp:
    limit_bg = plot.MakeTChain(args.bg_exp, 'limit')
    best_bg = plot.TGraphFromTree(limit_bg, args.x, args.y, 'deltaNLL == 0')
    plot.RemoveGraphXDuplicates(best_bg)
if args.sm_exp:
    limit_sm = plot.MakeTChain(args.sm_exp, 'limit')
    best_sm = plot.TGraphFromTree(limit_sm, args.x, args.y, 'deltaNLL == 0')
    plot.RemoveGraphXDuplicates(best_sm)
hists.SetMaximum(20)
hists.SetMinimum(0)
hists.SetContour(255)

c2=ROOT.TCanvas()
c2.SetRightMargin(0.2)
c2.SetLeftMargin(0.17)
hists.Draw("COLZ")

cont68 = contourFromTH2(hists, ROOT.Math.chisquared_quantile_c(1 - 0.68, 2), 10, 20)
cont95 = contourFromTH2(hists, ROOT.Math.chisquared_quantile_c(1 - 0.95, 2), 10, 20)
c2.cd()
for i, p in enumerate(cont68):
    p.SetLineStyle(1)
    p.SetLineWidth(1)
    p.SetLineColor(ROOT.kWhite)
    p.Draw("C SAME")
for i, p in enumerate(cont95):
    p.SetLineStyle(2)
    p.SetLineWidth(1)
    p.SetLineColor(ROOT.kWhite)
    p.Draw("C SAME")
t = 'asi' if 'asi' in args.files[0] else 'obs'


# fcorrname = os.path.join(os.path.dirname(args.files[0]), 'corr.json')
# with open(fcorrname, 'r') as cfile:
#     corr = json.load(cfile)
# corrl = ROOT.TLatex()
# corrl.SetNDC()
# corrl.SetTextSize(0.033)
# corrl.SetTextColor(ROOT.kWhite)
# corrl.DrawLatex(0.65, 0.6, "\\rho = {:.2f}".format(corr['corr']))

legend = plot.PositionedLegend(0.3, 0.2, 3, 0.015)
legend.SetFillStyle(0)

smr = ROOT.TMarker(1., 0., 28)
smr.SetMarkerColor(ROOT.kWhite)
smr.SetMarkerSize(2.5)
smr.Draw()
legend.AddEntry(smr, "SM", "P")

if args.bestfitpoint is not None:
    fbs = ROOT.TFile(args.bestfitpoint, 'READ')
    tr = fbs.Get("limit")
    tr.GetEntry(0)
    rx = eval('tr.'+args.x)
    ry = eval('tr.'+args.y)

    c2.cd()
    print(rx)
    print(ry)
    print(xsec)
    bfr = ROOT.TMarker(rx*xsec[procx], ry*xsec[procy], 29)
    bfr.SetMarkerColor(ROOT.kWhite)
    bfr.SetMarkerSize(2.5)
    bfr.Draw()
    legend.AddEntry(bfr, "Best fit", "P")
"""

fname = 'higgsCombine.'+t+'_multidim.MultiDimFit.mH120.root'
if 'asi' in args.files[0]: fname = fname.replace('root', '123456.root')
fbs = ROOT.TFile(os.path.dirname(args.files[0])+'/'+fname, 'READ')


"""



legend.AddEntry(cont68[0], "68% CL", "L")
legend.AddEntry(cont95[0], "95% CL", "L")
legend.Draw()

pt = ['', '']
for ip, p in enumerate([procx, procy]):
    if p == 'TTTT': pt[ip] = 't#bar{t}t#bar{t}'
    elif p == 'TTT': pt[ip] = 't#bar{t}t'
    elif p == 'TTW': pt[ip] = 't#bar{t}W'
    elif p == 'TTZ': pt[ip] = 't#bar{t}Z'
    elif p == 'TTH': pt[ip] = 't#bar{t}H'
hists.GetXaxis().SetTitle('#kappa_{t} cos #alpha')
hists.GetYaxis().SetTitle('#kappa_{t} sin #alpha')
hists.GetZaxis().SetTitle('-2#Delta ln(L)')



if debug is not None:
    debug.WriteTObject(hists, 'hist')
    for i, cont in enumerate(cont68):
        debug.WriteTObject(cont, 'cont_1sigma_%i' % i)
    for i, cont in enumerate(cont95):
        debug.WriteTObject(cont, 'cont_2sigma_%i' % i)


c2.Print("pics/"+args.output+"_heatmap.pdf")
c2.Print("pics/"+args.output+"_heatmap.png")


if debug is not None:
    debug.Close()

sys.exit()

axis = ROOT.TH2D(hists.GetName(),hists.GetName(),hists.GetXaxis().GetNbins(),0,hists.GetXaxis().GetXmax(),hists.GetYaxis().GetNbins(),0,hists.GetYaxis().GetXmax())
axis.Reset()
axis.GetXaxis().SetTitle(vtitle[0])
axis.GetXaxis().SetLabelSize(0.025)
axis.GetYaxis().SetLabelSize(0.025)
axis.GetYaxis().SetTitle(vtitle[1])
axis.GetXaxis().SetLimits(-0.5, 2.5)
axis.GetYaxis().SetLimits(-0.5, 2.5)

cont_1sigma = plot.contourFromTH2(hists, ROOT.Math.chisquared_quantile_c(1 - 0.68, 2), 10, frameValue=20)
cont_2sigma = plot.contourFromTH2(hists, ROOT.Math.chisquared_quantile_c(1 - 0.95, 2), 10, frameValue=20)

if debug is not None:
    debug.WriteTObject(hists, 'hist')
    for i, cont in enumerate(cont_1sigma):
        debug.WriteTObject(cont, 'cont_1sigma_%i' % i)
    for i, cont in enumerate(cont_2sigma):
        debug.WriteTObject(cont, 'cont_2sigma_%i' % i)

if args.sm_exp or args.bg_exp:
    legend = plot.PositionedLegend(0.5, 0.25, 3, 0.015)
else:
    legend = plot.PositionedLegend(0.3, 0.2, 3, 0.015)

pads[0].cd()
axis.Draw()

for i, p in enumerate(cont_2sigma):
    p.SetLineStyle(1)
    p.SetLineWidth(2)
    p.SetLineColor(ROOT.kBlack)
#    p.SetFillColor(45)
#    p.SetFillStyle(1001)
#    p.Draw("F SAME")
    p.Draw("L SAME")
legend.AddEntry(cont_2sigma[0], "95% CL", "F")

for i, p in enumerate(cont68):
    p.SetLineStyle(1)
    p.SetLineWidth(2)
    p.SetLineColor(ROOT.kBlack)
#    p.SetFillColor(42)
#    p.SetFillStyle(1001)
#    p.Draw("F SAME")
    p.Draw("L SAME")
legend.AddEntry(cont68[0], "68% CL", "F")

best.SetMarkerStyle(34)
best.SetMarkerSize(3)
best.Draw("P SAME")
legend.AddEntry(best, "Best fit", "P")
if args.sm_exp:
    best_sm.SetMarkerStyle(33)
    best_sm.SetMarkerColor(1)
    best_sm.SetMarkerSize(3.0)
    best_sm.Draw("P SAME")
    legend.AddEntry(best_sm, "Expected for 125 GeV SM Higgs", "P")
if args.bg_exp:
    best_bg.SetMarkerStyle(33)
    best_bg.SetMarkerColor(46)
    best_bg.SetMarkerSize(3)
    best_bg.Draw("P SAME")
    legend.AddEntry(best_bg, "Expected for background only", "P")


if args.mass:
    legend.SetHeader("m_{#phi} = "+args.mass+" GeV")
legend.Draw("SAME")
if args.sm_exp:
    overlayLegend,overlayGraphs = plot.getOverlayMarkerAndLegend(legend, {legend.GetNRows()-1 : best_sm}, {legend.GetNRows()-1 : {"MarkerColor" : 2}}, markerStyle="P")

plot.DrawCMSLogo(pads[0], 'CMS', args.cms_sub, 11, 0.045, 0.035, 1.2, '', 1.0)
plot.DrawTitle(pads[0], args.title_right, 3)
plot.DrawTitle(pads[0], args.title_left, 1)
plot.FixOverlay()
if args.sm_exp:
    best_sm.Draw("P SAME")
    for overlayGraph in overlayGraphs:
        overlayGraph.Draw("P SAME")
    overlayLegend.Draw("SAME")
canv.Print('.pdf')
canv.Close()

if debug is not None:
    debug.Close()
