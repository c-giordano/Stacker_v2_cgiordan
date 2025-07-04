#!/bin/env python

import os, sys, ROOT, json

ROOT.gROOT.SetBatch(1)

f = ROOT.TFile(sys.argv[1])
hcorr = f.Get("correlation_matrix")
ibx, iby = -1, -1
for ib in range(1, hcorr.GetXaxis().GetNbins()+1):
    if sys.argv[2] == hcorr.GetXaxis().GetBinLabel(ib):
        ibx = ib
        break
for ib in range(1, hcorr.GetYaxis().GetNbins()+1):
    if sys.argv[3] == hcorr.GetYaxis().GetBinLabel(ib):
        iby = ib
        break
if ibx < 0 or iby < 0:
    print('Cannot find the requested parameters')
    sys.exit()
corr = {'corr': hcorr.GetBinContent(ibx, iby)}
with open('corr.json', 'w') as outf:
    json.dump(corr, outf)

