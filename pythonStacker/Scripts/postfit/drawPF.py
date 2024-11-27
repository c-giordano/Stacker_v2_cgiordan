import CombineHarvester.CombineTools.plotting as plot
import ROOT
import math
import argparse
import json
import sys
import os ,shutil
import fnmatch
from array import array
import ctypes

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)

def turnToGraph(h1):
    x = []
    y = []
    ex = []
    ey = []
    n = h1.GetNbinsX()
    for i in range(1, n+1):
        val = h1.GetBinContent(i)
        err = h1.GetBinError(i)
        if math.isnan(val):
            h1.SetBinContent(i, 0.00001)
            h1.SetBinError(i, 1.0)
    h1.SetBinErrorOption( ROOT.TH1.kPoisson )
    grr = ROOT.TGraphAsymmErrors( h1 )
    for ip in range(grr.GetN()):
        grr.SetPointEXhigh(ip, 0.0)
        grr.SetPointEXlow(ip, 0.0)
#    for i in range(1, n+1):
#        val = h1.GetBinContent(i)
#        err = h1.GetBinError(i)
#        if math.isnan(val):
#            val = 0.00001f
#            err = 1.0
#        x.append(h1.GetBinCenter(i))
#        y.append(val)
#        ex.append(0.0)
#        ey.append(err)
#    grr = ROOT.TGraphAsymmErrors(n, array('d',x), array('d',y), array('d',ex), array('d',ex), array('d',ey), array('d',ey))
    return grr    

def createAxisHists(n,src,xmin=0,xmax=499):
    result = []
    for i in range(0,n):
        res = src.Clone()
        res.Reset()
        res.SetTitle("")
        res.SetName("axis%(i)d"%vars())
        res.SetAxisRange(xmin,xmax)
        res.SetStats(0)
        result.append(res)
    return result

def getHistogram(fname, histname, dirname='', postfitmode='prefit', allowEmpty=False, logx=False):
    outname = fname.GetName()
    for key in fname.GetListOfKeys():
        histo = fname.Get(key.GetName())
        dircheck = False
        if dirname == '' : dircheck=True
        elif dirname in key.GetName(): dircheck=True
        if isinstance(histo,ROOT.TH1F) and key.GetName()==histname:
            if logx:
                bin_width = histo.GetBinWidth(1)
                xbins = []
                xbins.append(bin_width - 1)
                axis = histo.GetXaxis()
                for i in range(1,histo.GetNbinsX()+1):
                    xbins.append(axis.GetBinUpEdge(i))
                rethist = ROOT.TH1F(histname,histname,histo.GetNbinsX(),array('d',xbins))
                rethist.SetBinContent(1,histo.GetBinContent(1)*(histo.GetBinWidth(1)-(bin_width - 1))/(histo.GetBinWidth(1)))
                rethist.SetBinError(1,histo.GetBinError(1)*(histo.GetBinWidth(1)-(bin_width - 1))/(histo.GetBinWidth(1)))
                for i in range(2,histo.GetNbinsX()+1):
                    rethist.SetBinContent(i,histo.GetBinContent(i))
                    rethist.SetBinError(i,histo.GetBinError(i))
                histo = rethist
            return [histo,outname]
        elif isinstance(histo,ROOT.TDirectory) and postfitmode in key.GetName() and dircheck:
            return getHistogram(histo,histname, allowEmpty=allowEmpty, logx=logx)
    print('Failed to find %(postfitmode)s histogram with name %(histname)s in file %(fname)s '%vars())
    if allowEmpty: return [ROOT.TH1F('empty', '', 1, 0, 1), outname]
    else: return None
  
plot.ModTDRStyle(r=0.04, l=0.14)
ROOT.gStyle.SetHatchesLineWidth(1)
ROOT.gStyle.SetHatchesSpacing(0.2)
ROOT.gStyle.SetEndErrorSize(0.0)
ROOT.gStyle.SetCanvasBorderMode(0)
ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
ROOT.gStyle.SetCanvasDefH(941)
ROOT.gStyle.SetCanvasDefW(800)
ROOT.gStyle.SetCanvasDefX(0)
ROOT.gStyle.SetCanvasDefY(0)
###ROOT.gStyle.SetErrorX(0.0001)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', help='Input file if shape file has already been created')
parser.add_argument('--ratio', default=False, action='store_true', help='Draw ratio plot')
parser.add_argument('--custom_x_range', help='Fix x axis range', action='store_true', default=False)
parser.add_argument('--x_axis_min', help='Fix x axis minimum', default=0.0)
parser.add_argument('--x_axis_max', help='Fix x axis maximum', default=1000.0)
parser.add_argument('--custom_y_range', help='Fix y axis range', action='store_true', default=False)
parser.add_argument('--y_axis_min', help='Fix y axis minimum', default=0.5)
parser.add_argument('--y_axis_max', help='Fix y axis maximum', default=100000.0)
parser.add_argument('--log_y', action='store_true', help='Use log for y axis')
parser.add_argument('--log_x', action='store_true', help='Use log for x axis')
parser.add_argument('--extra_pad', help='Fraction of extra whitespace at top of plot', default=0.4)
parser.add_argument('--outname', default='sbordered', help='Output plot name')
parser.add_argument('--mode', default='postfit', help='Prefit or postfit')
parser.add_argument('--ratio_range', help='y-axis range for ratio plot in format MIN,MAX', default="0.5,2.0")
parser.add_argument('--x_title', default='log_{10}(S/B)', help='Title for the x-axis')
parser.add_argument('--y_title', default='Events / bin', help='Title for the y-axis')
parser.add_argument('--region', default='', help='Name of kinematic region')
parser.add_argument('--lumi', default='138 fb^{-1} (13 TeV)', help='Lumi label')
parser.add_argument('--split', default=False, action='store_true', help='Split into different processes')
parser.add_argument('--asimov', default=False, action='store_true', help='Do not plot data')
parser.add_argument('--plot_directory', default='/SS2L_3L_fit/post-fit/')

args = parser.parse_args()
extra_pad = float(args.extra_pad)
custom_x_range = args.custom_x_range
custom_y_range = args.custom_y_range
x_axis_min = float(args.x_axis_min)
x_axis_max = float(args.x_axis_max)
y_axis_min = float(args.y_axis_min)
y_axis_max = float(args.y_axis_max)
log_y = args.log_y
log_x = args.log_x
ratio_range = args.ratio_range
reg = args.region

#hashcol = ROOT.TColor.GetColor('#898d8d')
#hashcol = ROOT.TColor.GetColor('#a2acab')
hashcol = ROOT.TColor.GetColor('#bec6c4')
sigcol = ROOT.TColor.GetColor("#ffa90e")

add = bool('njets' in args.file or 'loose' in args.file or 'medium' in args.file or 'ht' in args.file)

if reg == "bdt_3L_ttbar": extra_pad = 0.40
elif reg == "bdt_3L_sig": 
    extra_pad = 0.40
    if args.mode == 'prefit':
        ratio_range = "0.0,2.5"
elif reg == "bdt_3L_ttw":
    extra_pad = 0.40
    ratio_range = "0.61,1.39"
    if args.mode == 'prefit':
        ratio_range = "0.5,1.5"
elif reg == "bdt_4L_sig": extra_pad = 0.45
elif reg == "bdt_4L_ttw": extra_pad = 0.45
elif reg == "bdt_DL_ttbar": 
    extra_pad = 0.45
    if args.mode == 'prefit':
        ratio_range = "0.0,2.5"
elif reg == "bdt_DLem_sig":
    extra_pad = 0.35
    ratio_range = "0.01,1.99"
    if args.mode == 'prefit':
        ratio_range = "0.0,2.5"
elif reg == "bdt_DLmm_sig":
    extra_pad = 0.40
    if args.mode == 'prefit':
        ratio_range = "0.0,2.5"
elif reg == "bdt_DLmm_ttw":
    extra_pad = 0.45
elif reg == "bdt_DLee_sig":
    extra_pad = 0.45
    if args.mode == 'prefit':
        ratio_range = "0.0,2.5"
elif reg == "bdt_DLee_ttw":
    extra_pad = 0.45
    ratio_range = "0.41,1.59"
    if args.mode == 'prefit':
        ratio_range = "0.21,1.79"
elif reg == "bdt_DLem_ttw":
    extra_pad = 0.45
    ratio_range = "0.45,1.55"
elif reg == "cro":
    if not add:
        extra_pad = 0.40
        ratio_range = "0.85,1.15"
        if args.mode == 'prefit':
            ratio_range = "0.50,1.50"
elif reg == "cro_crw":
    if not add:
        extra_pad = 0.40
        ratio_range = "0.85,1.15"
        if args.mode == 'prefit':
            ratio_range = "0.55,1.45"
elif reg == "cro_3L":
    if not add:
        extra_pad = 0.40
        ratio_range = "0.75,1.25"
        if args.mode == 'prefit':
            ratio_range = "0.50,1.50"
elif reg == "crw":
    if not add:
        extra_pad = 0.40
        ratio_range = "0.42,1.58"
        if args.mode == 'prefit':
            ratio_range = "0.25,1.75"
elif reg == "crz":
    if not add:
        extra_pad = 0.35
        ratio_range = "0.1,1.9"
        if args.mode == 'prefit':
            ratio_range = "0.0,2.2"
elif reg == "crz_4L":
    if not add:
        extra_pad = 0.35
        ratio_range = "0.01,2.8"
        if args.mode == 'prefit':
            ratio_range = "0.01,3.0"

        
if add:
    extra_pad = 0.40
    if 'nb_loose' in args.file:
        if reg == "bdt_DLem_sig":
            ratio_range = "0.61,1.39"
    if 'njets' in args.file:
        if reg == "bdt_3L_ttw":
            ratio_range = "0.0,3.0"
        elif reg == "bdt_3L_sig":
            ratio_range = "0.0,5.0"
        elif reg == "bdt_3L_ttw":
            ratio_range = "0.0,5.0"
        elif reg == "bdt_3L_ttbar":
            ratio_range = "0.0,3.0"
        elif reg == "bdt_DLee_sig":
            ratio_range = "0.0,3.0"
        elif reg == "bdt_DLee_ttw":
            ratio_range = "0.0,3.0"
        elif reg == "bdt_DLem_ttw":
            ratio_range = "0.0,2.0"
        elif reg == "bdt_DLmm_sig":
            ratio_range = "0.0,2.0"
        elif reg == "bdt_DLmm_ttw":
            ratio_range = "0.0,2.5"
        elif reg == "cro":
            ratio_range = "0.8,1.2"
            if args.mode == 'prefit':
                ratio_range = "0.5,1.5"
        elif reg == "crw":
            ratio_range = "0.9,1.1"
        elif reg == "crz_4L":
            ratio_range = "0.0,2.5"
            if args.mode == 'prefit':
                ratio_range = "0.0,3.5"
        if reg == "sig":
            if args.mode == 'prefit':
                ratio_range = "0.0,2.5"

histo_file = ROOT.TFile(args.file)

bkghist = getHistogram(histo_file, 'TotalBkg', '', args.mode, logx=log_x)[0]
splusbhist = getHistogram(histo_file, 'TotalProcs', '', args.mode, logx=log_x)[0]
total_datahist = getHistogram(histo_file, 'data_obs', '', args.mode, logx=log_x)[0]
sighist = getHistogram(histo_file, 'TotalSig', '', args.mode, logx=log_x)[0]
sighist_forratio = sighist.Clone()
sighist_forratio.SetName('sighist_forratio')

dirs = []
for d in histo_file.GetListOfKeys():
    dname = d.GetName()
    dir = histo_file.Get(dname)
    if not dname.endswith('_'+args.mode): continue
    dirs.append(dname)

procSel, procSelCol, procName = [], [], []
if '_sig' in args.outname:
    procs = ['TTT', 'other_t', 'Xgamma', 'ChargeMisID', 'WZ', 'VVV', 'nonPromptElectron', 'nonPromptMuon', 'TTH', 'TTZ', 'TTW', 'TotalSig']
    procnames = ['ttt', 'Other t', 'X#gamma', 'ChargeMisID', 'WZ', 'VV(V)', 'NP e', 'NP #mu', 't#bar{t}H', 't#bar{t}Z', 't#bar{t}W', 't#bar{t}t#bar{t}']
    proccol = ["#e9d98d", '#92dadd', "#3f90da", "#832db6", "#bd1f01", "#bd1f01", "#e76300", "#a96b59", "#b9ac70", "#717581", "#94a4a2", "#ffa90e"]
else:
    procs = ['TTT', 'other_t', 'Xgamma', 'ChargeMisID', 'WZ', 'VVV', 'nonPromptElectron', 'nonPromptMuon', 'TTH', 'TTZ', 'TTW', 'TotalSig']
    procnames = ['ttt', 'Other t', 'X#gamma', 'ChargeMisID', 'WZ', 'VV(V)', 'NP e', 'NP #mu', 't#bar{t}H', 't#bar{t}Z', 't#bar{t}W', 't#bar{t}t#bar{t}']
    proccol = ["#e9d98d", '#92dadd', "#3f90da", "#832db6", "#bd1f01", "#bd1f01", "#e76300", "#a96b59", "#b9ac70", "#717581", "#94a4a2", "#ffa90e"]

prochist = []
for ip, p in enumerate(procs):
    h = None
#    print(p)
    for i, d in enumerate(dirs):
#        print(d)
        hc = getHistogram(histo_file, p, '', d, logx=log_x)
        if hc == None: continue
        else: hc = hc[0]
        if h == None:
            h = hc
        else:
            h.Add(hc)
    if h != None:
        prochist.append(h)
        procSel.append(p)
        procSelCol.append(proccol[ip])
        procName.append(procnames[ip])

        print("yield for {}: {} & {}".format(p, h.Integral(), h.GetBinError(1)))

bkghist.SetFillColor(hashcol)
bkghist.SetLineColor(ROOT.kBlack)
bkghist.SetMarkerSize(0)

sighist.SetFillColor(sigcol)
sighist.SetLineWidth(0)
sighist.SetMarkerSize(0)

stack = ROOT.THStack("hs", "")
stack.Add(bkghist)
stack.Add(sighist)

data = {}
procstack = ROOT.THStack("hsproc", "")
for ip, p in enumerate(prochist):
    if isinstance(procSelCol[ip], str):
        ci = ROOT.TColor.GetColor(procSelCol[ip])
        p.SetLineColor(ci)
        p.SetFillColor(ci)
    else:
        p.SetLineColor(procSelCol[ip])
        p.SetFillColor(procSelCol[ip])
    p.SetMarkerSize(0)
    procstack.Add(p)
    error = ctypes.c_double(0.)
    integral = p.IntegralAndError(0, p.GetNbinsX()+1, error)
    data[procSel[ip]] = [integral, error.value]
error_data = ctypes.c_double(0.)
integral_data = total_datahist.IntegralAndError(0, total_datahist.GetNbinsX()+1, error_data)
data['data'] = [integral_data, error_data.value]
with open(args.outname+'.json', 'w') as f:
    json.dump(data, f)

c2 = ROOT.TCanvas()
c2.cd()
if args.ratio: pads=plot.TwoPadSplit(0.29,0.01,0.01)
else: pads=plot.OnePad()
pads[0].cd()
if(log_y): pads[0].SetLogy(1)
if(log_x): pads[0].SetLogx(1)
##if ('nb_loose' in args.file or 'nb_medium' in args.file) and 'OneMedB' not in args.file:
##    pads[0].SetLogy(1)
if custom_x_range: 
    if x_axis_max > bkghist.GetXaxis().GetXmax(): x_axis_max = bkghist.GetXaxis().GetXmax()
if args.ratio:
    if(log_x): pads[1].SetLogx(1)
    axish = createAxisHists(2,bkghist,bkghist.GetXaxis().GetXmin(),bkghist.GetXaxis().GetXmax()-0.01)
#    axish[1].GetXaxis().SetTitle(args.x_title)
    axish[1].GetXaxis().SetTitle(bkghist.GetXaxis().GetTitle())
    axish[1].GetYaxis().SetNdivisions(4)
    #axish[1].GetYaxis().SetTitle("Data / Pred.")
    axish[1].GetYaxis().SetTitle("")
    #axish[1].GetYaxis().SetTitleSize(0.04)
    #axish[1].GetYaxis().SetLabelSize(0.04)
    axish[1].GetYaxis().SetTitleOffset(1.45)
    #axish[0].GetYaxis().SetTitleSize(0.04)
    #axish[0].GetYaxis().SetLabelSize(0.04)
    axish[0].GetYaxis().SetTitleOffset(1.45)
    axish[0].GetXaxis().SetTitleSize(0)
    axish[0].GetXaxis().SetLabelSize(0)
    if custom_x_range:
        axish[0].GetXaxis().SetRangeUser(x_axis_min,x_axis_max-0.01)
        axish[1].GetXaxis().SetRangeUser(x_axis_min,x_axis_max-0.01)
    if custom_y_range:
        axish[0].GetYaxis().SetRangeUser(y_axis_min,y_axis_max)
        axish[1].GetYaxis().SetRangeUser(y_axis_min,y_axis_max)
    if 'njets_srFit_pf_variables_sig_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(4.0, 10.0)
        axish[1].GetXaxis().SetRangeUser(4.0, 10.0)
        ratio_range = '0.0,2.9'
    elif 'njets_srFit_pf_variables_ttw_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(3.0, 8.0)
        axish[1].GetXaxis().SetRangeUser(3.0, 8.0)
        ratio_range = '0.4,1.6'
    elif 'njets_srFit_pf_variables_ttbar_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(3.0, 8.0)
        axish[1].GetXaxis().SetRangeUser(3.0, 8.0)
        ratio_range = '0.4,1.6'
    elif 'njets_srFit_pf_variables_crz_4L_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(2.0, 6.0)
        axish[1].GetXaxis().SetRangeUser(2.0, 6.0)
        ratio_range = '0.0,2.9'
    elif 'nb_loose_srFit_pf_variables_ttw_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(2.0, 5.0)
        axish[1].GetXaxis().SetRangeUser(2.0, 5.0)
        ratio_range = '0.4,1.6'
    elif 'nb_loose_srFit_pf_variables_ttbar_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(2.0, 4.0)
        axish[1].GetXaxis().SetRangeUser(2.0, 4.0)
        ratio_range = '0.6,1.4'
    elif 'nb_loose' in args.outname and 'crz_4L' in args.outname:
        axish[0].GetXaxis().SetRangeUser(1.0, 3.0)
        axish[1].GetXaxis().SetRangeUser(1.0, 3.0)
    elif 'nb_medium' in args.outname and 'crz_4L' in args.outname:
        axish[0].GetXaxis().SetRangeUser(1.0, 2.0)
        axish[1].GetXaxis().SetRangeUser(1.0, 2.0)
    elif 'nb_medium' in args.outname and 'cro_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(0.0, 4.0)
        axish[1].GetXaxis().SetRangeUser(0.0, 4.0)
    elif 'nb_medium' in args.outname and 'sig' in args.outname:
        ratio_range = '0.0,2.0'
    elif 'nb_medium_srFit_pf_variables_sig_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(1.0, 4.0)
        axish[1].GetXaxis().SetRangeUser(1.0, 4.0)
        ratio_range = '0.5,1.5'
    elif 'nb_medium_srFit_pf_variables_ttw_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(0.0, 3.0)
        axish[1].GetXaxis().SetRangeUser(0.0, 3.0)
        ratio_range = '0.5,1.5'
    elif 'nb_medium_srFit_pf_variables_ttbar_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(0.0, 3.0)
        axish[1].GetXaxis().SetRangeUser(0.0, 3.0)
        ratio_range = '0.5,1.5'
    elif 'ht_srFit_pf_variables_sig_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(270, 1150)
        axish[1].GetXaxis().SetRangeUser(270, 1150)
        ratio_range = '0.0,2.5'
        if 'prefit' in args.outname:
            ratio_range = '0.0,3.5'
    elif 'ht_srFit_pf_variables_ttw_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(200, 1150)
        axish[1].GetXaxis().SetRangeUser(200, 1150)
        ratio_range = '0.3,1.7'
    elif 'ht_srFit_pf_variables_ttbar_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(270, 1100)
        axish[1].GetXaxis().SetRangeUser(270, 1100)
        ratio_range = '0.0,3.9'        
    elif 'met_srFit_pf_variables_sig_p' in args.outname:
        extra_pad = 0.40
        ratio_range = '0.0,3.5'
    elif 'met_srFit_pf_variables_ttw_p' in args.outname:
        extra_pad = 0.40
        ratio_range = '0.0,4.0'
    elif 'met_srFit_pf_variables_ttbar_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(0, 210)
        axish[1].GetXaxis().SetRangeUser(0, 210)
        extra_pad = 0.40
        ratio_range = '0.1,1.9'
    elif 'nb_medium_crz_OneMedB_srFit_pf_variables_crz_p' in args.outname:
        axish[0].GetXaxis().SetRangeUser(1, 3)
        axish[1].GetXaxis().SetRangeUser(1, 3)
    if reg == "bdt_DLem_sig":
        axish[0].GetXaxis().SetRangeUser(0.1,1.0)
        axish[1].GetXaxis().SetRangeUser(0.1,1.0)
        
###    if reg == "bdt_DLem_sig":
###        axish[0].GetXaxis().SetRangeUser(0.75,3.0)
###        axish[1].GetXaxis().SetRangeUser(0.75,3.0)
#    if ('nb_loose' in args.file and 'DLem' in args.file) or \
#    ('nb_medium' in args.file and 'DLem' in args.file):
#        axish[0].GetXaxis().SetRangeUser(1.5,3.5)
#        axish[1].GetXaxis().SetRangeUser(1.5,3.5)
else:
    axish = createAxisHists(1,bkghist,bkghist.GetXaxis().GetXmin(),bkghist.GetXaxis().GetXmax()-0.01)
    #  axish[0].GetYaxis().SetTitleOffset(1.4)
    if custom_x_range: axish[0].GetXaxis().SetRangeUser(x_axis_min,x_axis_max-0.01)
    if custom_y_range: axish[0].GetYaxis().SetRangeUser(y_axis_min,y_axis_max)

axish[0].GetYaxis().SetTitle(args.y_title)
axish[0].GetXaxis().SetTitle(args.x_title)
if not custom_y_range: axish[0].SetMaximum(extra_pad*bkghist.GetMaximum())
if not custom_y_range: 
###    if(log_y): axish[0].SetMinimum(0.0009)
    if(log_y): axish[0].SetMinimum(1.0)
    else:
        if not(('nb_loose' in args.file or 'nb_medium' in args.file) and 'OneMedB' not in args.file):
            axish[0].SetMinimum(0)
        else: axish[0].SetMinimum(1.0)
axish[0].Draw()

if not args.split: stack.Draw("histsame")
else: procstack.Draw("histsame")
#splusbhist.SetFillColor(plot.CreateTransparentColor(12,0.4))
splusbhist.SetFillColor(hashcol)
splusbhist.SetLineColor(0)
splusbhist.SetMarkerSize(0)
splusbhist.SetFillStyle(3944)
splusbhist.DrawCopy("e2same")
gr_data = turnToGraph(total_datahist)
#gr_data = ROOT.TGraph()
#for i in range(1, total_datahist.GetNbinsX()+1): 
#    gr_data.SetPoint(i-1, total_datahist.GetBinCenter(i), total_datahist.GetBinContent(i))
if not args.asimov:
    gr_data.SetMarkerStyle(20)
    gr_data.SetMarkerSize(1.4)
    gr_data.Draw("PS")
#total_datahist.SetMarkerStyle(20)
#total_datahist.Draw("PSAME")

if not args.split: legend = plot.PositionedLegend(0.40,0.10,3,0.03)
elif 'ttbar' in args.outname and 'njets' in args.outname: legend = plot.PositionedLegend(0.42,0.20,3,0.03, 0.05)
else: legend = plot.PositionedLegend(0.40,0.20,3,0.03, 0.05)
plot.Set(legend, NColumns=2)
legend.SetTextFont(42)
legend.SetTextSize(0.035)
legend.SetFillColor(0)
if not args.split:
    if not args.asimov:
        legend.AddEntry(gr_data,"Data", "PE")
    legend.AddEntry(bkghist, "Background", "f")
    legend.AddEntry(sighist, "t#bar{t}t#bar{t}", "f")
    legend.AddEntry(splusbhist, "Total unc.", "f")
else:
    if not args.asimov:
        legend.AddEntry(gr_data,"Data", "PE")
    legend.AddEntry(sighist, "t#bar{t}t#bar{t}", "f")
    for ip in range(len(procSel)):
        ipp = len(procSel)-ip-1
        p = procSel[ipp]
#    for ip, p in enumerate(procSel):
        if procSel[ipp] == 'TotalSig' or p not in procSel: continue
        if procSel[ipp] == 'TTT' and 'sig' not in args.outname: continue
        if procSel[ipp] == 'nonPromptMuon' and 'nonPromptElectron' in procSel: continue
        if procSel[ipp] == 'ttt' and 'Other t' in procSel and 'sig' not in args.outname: continue
        if procSel[ipp] == 'VVV' and 'WZ' in procSel: continue
        legend.AddEntry(prochist[ipp], procName[ipp], "f")
    legend.AddEntry(splusbhist, "Total unc.", "f")

legend.Draw("same")

plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), extra_pad if extra_pad>0 else 0.30)
#plot.DrawCMSLogo(pads[0], 'CMS', 'Preliminary', 11, 0.045, 0.05, 1.0, '', 1.0)
plot.DrawCMSLogo(pads[0], 'CMS', '', 11, 0.045, 0.05, 1.0, '', 1.0)

lumi =  args.lumi
if 'y1' in args.file:
    if 'y16' in args.file: lumi = '36.3 fb^{-1} (13 TeV)'
    elif 'y17' in args.file: lumi = '41.5 fb^{-1} (13 TeV)'
    elif 'y18' in args.file: lumi = '59.8 fb^{-1} (13 TeV)'
    plot.DrawTitle(pads[0], lumi, 3)
else:
    plot.DrawTitle(pads[0], args.lumi, 3)

if args.ratio:
    ratio_bkghist = plot.MakeRatioHist(bkghist,bkghist,True,False)
    ratio_bkghist.SetFillColor(hashcol)
    ratio_bkghist.SetFillStyle(3944)
    ratio_datahist = plot.MakeRatioHist(total_datahist,splusbhist,True,False)
    gr_ratio_data = turnToGraph(ratio_datahist)
    gr_ratio_data.SetMarkerStyle(20)
    gr_ratio_data.SetMarkerSize(1.4)
#    ratio_sighist = plot.MakeRatioHist(sighist_forratio,bkghist,True,False)
#    ratio_sighist.SetLineColor(ROOT.kRed)
#    ratio_sighist.SetMarkerSize(0)
#    for i in range(1,ratio_sighist.GetNbinsX()+1): ratio_sighist.SetBinContent(i,ratio_sighist.GetBinContent(i)+1)
    pads[1].cd()
    pads[1].SetGrid(0,1)
    axish[1].Draw("axis")
    axish[1].SetMinimum(float(ratio_range.split(',')[0]))
    axish[1].SetMaximum(float(ratio_range.split(',')[1]))
    ratio_bkghist.SetMarkerSize(0)
    ratio_bkghist.DrawCopy("e2same")
##    ratio_datahist.Draw("e0same")
#    ratio_sighist.Draw("histsame")
    gr_ratio_data.Draw("PS")
    pads[1].RedrawAxis("G")

# "CR-2l-23j1b", "CR-2l-45j2b", "CR-3l-2j1b", "CR-3l-Z", "CR-4l-Z"
xname = axish[1].GetXaxis().GetTitle()
yname = ""
regName = reg

if reg == "bdt_DLee_sig":
    regName = "\\splitline{\\text{SR-2}\\ell,\,\\text{ee}}{\\mathrm{t\\bar{t}t\\bar{t}\\;class}}"
elif reg == "bdt_DLmm_sig":
    regName = "\\splitline{\\text{SR-2}\\ell,\,\\mu\\mu}{\\mathrm{t\\bar{t}t\\bar{t}\\;class}}"
    if 'njets_DL' in args.file or 'nb_medium_DL' in args.file or 'nb_loose_DL' in args.file or 'ht_DL' in args.file:
        regName = "\\splitline{\\text{SR-2}\\ell}{\\mathrm{t\\bar{t}t\\bar{t}\\;class}}"
elif reg == "bdt_DLem_sig":
    regName = "\\splitline{\\text{SR-2}\\ell,\,\\text{e}\\mu}{\\mathrm{t\\bar{t}t\\bar{t}\\;class}}"
elif reg == "bdt_DLee_ttw":
    regName = "\\splitline{\\text{SR-2}\\ell,\,\\text{ee}}{\\mathrm{t\\bar{t}X\\;class}}"
elif reg == "bdt_DLmm_ttw":
    regName = "\\splitline{\\text{SR-2}\\ell,\,\\mu\\mu}{\\mathrm{t\\bar{t}X\\;class}}"
    if 'njets_DL' in args.file or 'nb_medium_DL' in args.file or 'nb_loose_DL' in args.file or 'ht_DL' in args.file:
        regName = "\\splitline{\\text{SR-2}\\ell}{\\mathrm{t\\bar{t}X\\;class}}"
elif reg == "bdt_DLem_ttw":
    regName = "\\splitline{\\text{SR-2}\\ell,\,\\text{e}\\mu}{\\mathrm{t\\bar{t}X\\;class}}"
elif reg == "bdt_DL_ttbar":
    regName = "\\splitline{\\text{SR-2}\\ell}{\\mathrm{t\\bar{t}\\;class}}"
elif reg == "bdt_3L_sig":
    regName = "\\splitline{\\text{SR-3}\\ell}{\\mathrm{t\\bar{t}t\\bar{t}\\;class}}"
elif reg == "bdt_3L_ttw":
    regName = "\\splitline{\\text{SR-3}\\ell}{\\mathrm{t\\bar{t}X\\;class}}"
elif reg == "bdt_3L_ttbar":
    regName = "\\splitline{\\text{SR-3}\\ell}{\\mathrm{t\\bar{t}\\;class}}"
elif reg == "bdt_4L_sig":
    regName = "\\splitline{\\text{SR-4}\\ell}{\\mathrm{t\\bar{t}t\\bar{t}\\;class}}"
elif reg == "bdt_4L_ttw":
    regName = "\\splitline{\\text{SR-4}\\ell}{\\mathrm{t\\bar{t}X\\;class}}"
elif reg == "cro":
    regName = "\\splitline{\\text{CR-2}\\ell\\text{-23j1b}}{}"
elif reg == "cro_crw":
    regName = "\\splitline{\\text{CR-2}\\ell\\text{-23j1b}}{\\text{CR-2}\\ell\\text{-45j2b}}"
elif reg == "cro_3L":
    regName = "\\splitline{\\text{CR-3}\\ell\\text{-2j1b}}{}"
elif reg == "crw":
    regName = "\\splitline{\\text{CR-2}\\ell\\text{-45j2b}}{}"
elif reg == "crz":
    regName = "\\splitline{\\text{CR-3}\\ell\\text{-Z}}{}"
elif reg == "crz_4L":
    regName = "\\splitline{\\text{CR-4}\\ell\\text{-Z}}{}"
elif reg == "sig":
    regName = "\\splitline{\\text{SR-2}\\ell,\\,\\text{SR-3}\\ell}{\\mathrm{t\\bar{t}t\\bar{t}\\;class}}"
elif reg == "ttw":
    regName = "\\splitline{\\text{SR-2}\\ell,\\,\\text{SR-3}\\ell}{\\mathrm{t\\bar{t}X\\;class}}"
elif reg == "ttbar":
    regName = "\\splitline{\\text{SR-2}\\ell,\\,\\text{SR-3}\\ell}{\\mathrm{t\\bar{t}\\;class}}"

axish[0].GetYaxis().SetLabelSize(0.05)
axish[1].GetYaxis().SetLabelSize(0.05)
axish[1].GetXaxis().SetLabelSize(0.05) # default is 0.04
    
if args.file.split('/')[-1].startswith('srFit_pf') or args.file.split('/')[-1].startswith('nb_loose_srFit_pf') or args.file.split('/')[-1].startswith('nb_medium_srFit_pf') or args.file.split('/')[-1].startswith('ht_srFit_pf'):
    if reg == "bdt_DLee_sig":
        xname = "BDT score t#bar{t}t#bar{t}"
        yname = "Events / 0.25 units"
    elif reg == "bdt_DLmm_sig":
        xname = "BDT score t#bar{t}t#bar{t}"
        yname = "Events / 0.25 units"
    elif reg == "bdt_DLem_sig":
        xname = "BDT score t#bar{t}t#bar{t}"
        yname = "Events / 0.11 units"
    elif reg == "bdt_DLee_ttw":
        xname = "BDT score t#bar{t}X"
        yname = "Events / bin"
    elif reg == "bdt_DLmm_ttw":
        xname = "BDT score t#bar{t}X"
        yname = "Events / 0.2 units"
    elif reg == "bdt_DLem_ttw":
        xname = "BDT score t#bar{t}X"
        yname = "Events / bin"
    elif reg == "bdt_DL_ttbar":
        xname = "BDT score t#bar{t}"
        yname = "Events / bin"
    elif reg == "bdt_3L_sig":
        xname = "BDT score t#bar{t}t#bar{t}"
        yname = "Events / 0.2 units"
    elif reg == "bdt_3L_ttw":
        xname = "BDT score t#bar{t}X"
        yname = "Events / 0.33 units"
    elif reg == "bdt_3L_ttbar":
        xname = "BDT score t#bar{t}"
        yname = "Events / bin"
    elif reg == "bdt_4L_sig":
        xname = "BDT score t#bar{t}t#bar{t}"
        yname = "Events / 0.5 units"
    elif reg == "bdt_4L_ttw":
        xname = "BDT score t#bar{t}X"
        yname = "Events / 0.5 units"
    elif reg == "cro":
        xname = "BDT score t#bar{t}"
        yname = "Events / 0.2 units"
    elif reg == "cro_crw":
        xname = "BDT score t#bar{t}"
        yname = "Events / 0.2 units"
    elif reg == "cro_3L":
        xname = "Sum of lepton charges"
        yname = "Events / 1 unit"
        axish[0].GetXaxis().SetNdivisions(-102)
        axish[1].GetXaxis().SetNdivisions(-102)
        axish[1].GetXaxis().SetBinLabel(1, '<0')
        axish[1].GetXaxis().SetBinLabel(2, '>0')
        axish[1].GetXaxis().SetLabelSize(0.065)
    elif reg == "crw":
        xname = "BDT score t#bar{t}"
        yname = "Events / 0.2 units"
    elif reg == "crz":
        xname = "Number of jets"
        yname = "Events / 1 unit"
    elif reg == "crz_4L":
        xname = "Number of jets"
        yname = "Events / bin"

if axish[1].GetXaxis().GetTitle() == 'N_{jets}' and add:
    xname = 'Number of jets'
    yname = 'Events / 1 unit'
        
if add:
    if 'nb' in args.file.split('/')[-1]:
        xname = 'Number of b jets'
        if 'nb_medium' in args.file.split('/')[-1]:
            xname = 'N_{b}^{medium}'
        yname = 'Events / 1 unit'
    elif 'njets' in args.file.split('/')[-1]:
        xname = 'Number of jets'
        yname = 'Events / 1 unit'
    elif 'ht' in args.file.split('/')[-1]:
        xname = 'H_{T} [GeV]'
        yname = 'Events / 50 GeV'

if 'met' in args.file.split('/')[-1]:
    yname = 'Events / 50 GeV'

if 'Number of' in xname:
    axish[0].GetXaxis().SetNdivisions(8)
    axish[1].GetXaxis().SetNdivisions(8)
    if 'ttw' in args.outname and 'nb_loose' in args.outname:
        axish[0].GetXaxis().SetNdivisions(4)
        axish[1].GetXaxis().SetNdivisions(4)
    if 'cro_3L' in args.outname and 'nb' in args.outname:
        axish[0].GetXaxis().SetNdivisions(2)
        axish[1].GetXaxis().SetNdivisions(2)
    if 'crz_OneMedB' in args.outname and 'nb_loose' in args.outname:
        axish[0].GetXaxis().SetNdivisions(4)
        axish[1].GetXaxis().SetNdivisions(4)
    if ('cro_p' in args.outname or 'crz_4L' in args.outname) and 'nb' in args.outname:
        axish[0].GetXaxis().SetNdivisions(4)
        axish[1].GetXaxis().SetNdivisions(4)
    if 'ttbar' in args.outname and 'nb' in args.outname:
        axish[0].GetXaxis().SetNdivisions(3)
        axish[1].GetXaxis().SetNdivisions(3)

if 'cro_3L' in args.outname and 'nb_medium' in args.outname:
    axish[0].GetXaxis().SetNdivisions(2)
    axish[1].GetXaxis().SetNdivisions(2)
if 'crz_4L' in args.outname and 'nb_medium' in args.outname:
    axish[0].GetXaxis().SetNdivisions(2)
    axish[1].GetXaxis().SetNdivisions(2)
        
label = ROOT.TMathText()
label.SetNDC()
label.SetTextAlign(13)
label.SetTextFont(42)
label.SetTextColor(ROOT.kBlack)
label.SetTextSize(0.067)
ylabel = 0.81
xlabel = 0.22
if '\\,' in regName:
    ylabel = 0.83
    xlabel = 0.20
label.DrawMathText(xlabel, ylabel, regName)
label.Draw()

#if add and 'OneMedB' in args.file:
#    labelb = ROOT.TMathText()
#    labelb.SetNDC()
#    labelb.SetTextAlign(13)
#    labelb.SetTextFont(42)
#    labelb.SetTextColor(ROOT.kBlack)
#    labelb.SetTextSize(0.055)
#    labelb.DrawMathText(0.65, 0.65, '\mathrm{N_{b}^{medium} \geq 1}')
#    labelb.Draw()
#
#pfl = ROOT.TMathText()
#pfl.SetNDC()
#pfl.SetTextAlign(13)
#pfl.SetTextFont(42)
#pfl.SetTextColor(ROOT.kBlack)
#pfl.SetTextSize(0.047)
#ylabel = 0.60
#xlabel = 0.70
#if 'nb_medium' in args.outname and 'sig' in args.outname:
#    xlabel = 0.20
#if 'nb_medium' in args.outname and 'ttw' in args.outname:
#    xlabel = 0.20
#if args.mode == 'prefit':
#    if 'bdt_4L' in args.file:
#        ylabel = 0.60
#        xlabel = 0.50
#    if 'cro_3L' in args.file:
#        ylabel = 0.70
#        xlabel = 0.70
#    if 'cro_crw' in args.file or ('cro' in args.file and not 'cro_3L' in args.file):
#        ylabel = 0.70
#        xlabel = 0.20
#    if 'nb_medium' in args.outname:
#        xlabel = 0.20
#    if 'nb_loose' in args.outname and 'cro_p' in args.outname:
#        xlabel = 0.70
#    pfl.DrawMathText(xlabel, ylabel, '\mathrm{Prefit}')
#else:
#    pfl.DrawMathText(xlabel, ylabel, '\mathrm{Postfit}')
#pfl.Draw()

xname = ""
yname = ""
axish[1].GetXaxis().SetTitle(xname)
axish[0].GetYaxis().SetTitle(yname)

pads[0].cd()
pads[0].GetFrame().Draw()
pads[0].RedrawAxis()

plot_directory_ = "/user/mshoosht/public_html/Interpretations/Plots/"+args.plot_directory+"/"
if not os.path.exists(plot_directory_):
	try:
		os.makedirs(plot_directory_)
	except OSError: # Resolve rare race condition
		pass

shutil.copyfile('/user/mshoosht/public_html/index.php', os.path.join( plot_directory_, 'index.php' ) )
print plot_directory_

c2.Print(plot_directory_+args.outname+".png")
c2.Print(plot_directory_+args.outname+".pdf")
c2.Print(plot_directory_+args.outname+".root")


