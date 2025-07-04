#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
cross_sec_reweight_legacy.py
Calcola xSec nominali + reweight per i modelli SS/SO/VS/VO.
Robusto contro file ROOT rotti: se un mass-point non contiene
dati validi, nel JSON compare comunque con campi null.
Compatibile con Python 2.7 e 3.x < 3.6 (niente type-hint, niente f-string).
"""

from __future__ import print_function, division
import glob
import json
import math
import sys
import argparse
import logging

import uproot
import awkward as ak
import numpy as np

# ----------------------------------------------------------------------
#  Dizionari dei pattern (li importi dal tuo repo)
# ----------------------------------------------------------------------
from xSec_setting import (
    dictionary_for_xSec_ss as dict_SS,
    dictionary_for_xSec_so as dict_SO,
    dictionary_for_xSec_vs as dict_VS,
    dictionary_for_xSec_vo as dict_VO,
)

MODEL_MAP = {
    "SS": (dict_SS, "results_SS_TTTJ.json"),
    "SO": (dict_SO, "results_SO_TTTJ.json"),
    "VS": (dict_VS, "results_VS_TTTJ.json"),
    "VO": (dict_VO, "results_VO_TTTJ.json"),
}

COUPLING_VALUES = [0.1, 0.2, 0.5, 2.0, 5.0]
COUPLING_KEYS   = ["c_0p1", "c_0p2", "c_0p5", "c_2", "c_5"]

# ----------------------------------------------------------------------
#  Helper
# ----------------------------------------------------------------------
def _none_if_nan(x):
    """Restituisce None (→ null in JSON) per NaN o None, altrimenti float(x)."""
    try:
        return None if x is None or math.isnan(float(x)) else float(x)
    except Exception:
        return None

def filter_valid_files(file_list, branches):
    """Ritorna i file che contengono tutti i rami con ≥ 1 entry."""
    out = []
    for f in file_list:
        try:
            tree = uproot.open(f)["Events"]
            arr  = tree.arrays(branches, entry_stop=1)
            if any(len(arr[b]) == 0 for b in branches):
                continue
            out.append(f)
        except Exception:
            logging.debug("Scarto %s (errore uproot)", f, exc_info=True)
    return out

def calc_nominal(files):
    if not files:
        return float("nan"), float("nan")
    arrays = []
    for f in files:
        arr = uproot.open(f)["Events"] \
              .arrays(["LHEWeight_originalXWGTUP"])["LHEWeight_originalXWGTUP"]
        arrays.append(arr)
    allw = ak.concatenate(arrays)
    n    = len(files)
    return ak.sum(allw) / n, math.sqrt(ak.sum(allw ** 2)) / n

def calc_reweighted(files):
    if not files:
        return [(float("nan"), float("nan"))] * len(COUPLING_VALUES)

    acc = [[] for _ in COUPLING_VALUES]
    for f in files:
        arr = uproot.open(f)["Events"].arrays(
            ["LHEWeight_originalXWGTUP", "LHEReweightingWeight"]
        )
        nom = arr["LHEWeight_originalXWGTUP"]
        rw  = arr["LHEReweightingWeight"]
        for i in range(len(COUPLING_VALUES)):
            acc[i].append(nom * rw[:, i])

    n = len(files)
    out = []
    for lst in acc:
        comb = ak.concatenate(lst)
        out.append((ak.sum(comb) / n, math.sqrt(ak.sum(comb ** 2)) / n))
    return out

# ----------------------------------------------------------------------
#  Core
# ----------------------------------------------------------------------
def run(model, verbose=False, max_files=0):
    logging.basicConfig(level=logging.INFO if verbose else logging.WARNING,
                        format="%(levelname)s - %(message)s", stream=sys.stdout)

    if model not in MODEL_MAP:
        logging.critical("Model '%s' non valido (%s)", model, ",".join(MODEL_MAP))
        sys.exit(1)

    patterns, out_json = MODEL_MAP[model]
    logging.info("Modello %s → %s", model, out_json)

    result = {}

    for mp, pattern in patterns.items():
        files = glob.glob(pattern, recursive=True)
        if max_files > 0:
            files = files[:max_files]
        logging.info("[%s] trovati %d file", mp, len(files))

        # nominali
        valid_nom = filter_valid_files(files, ["LHEWeight_originalXWGTUP"])
        x0, e0    = calc_nominal(valid_nom)

        # reweight
        valid_rw  = filter_valid_files(
            valid_nom,
            ["LHEWeight_originalXWGTUP", "LHEReweightingWeight"]
        )
        rews = calc_reweighted(valid_rw)

        entry = {"c_nominal": {"xsec": _none_if_nan(x0),
                               "error": _none_if_nan(e0)}}
        for ck, (x, e) in zip(COUPLING_KEYS, rews):
            entry[ck] = {"xsec": _none_if_nan(x),
                         "error": _none_if_nan(e)}
        result[mp] = entry

    with open(out_json, "w") as f:
        json.dump(result, f, indent=4, allow_nan=False)
    print("\nJSON scritto in {}".format(out_json))

# ----------------------------------------------------------------------
#  CLI
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Calcola xSec nominali + reweight per SS/SO/VS/VO"
    )
    ap.add_argument("--model", required=True, choices=list(MODEL_MAP.keys()))
    ap.add_argument("-v", "--verbose", action="store_true",
                    help="Log dettagliato")
    ap.add_argument("-n", "--max-files", type=int, default=0,
                    help="Numero massimo di file ROOT da usare per ciascun mass-point "
                         "(0 = tutti)")
    args = ap.parse_args()
    run(model=args.model, verbose=args.verbose, max_files=args.max_files)

if __name__ == "__main__":
    main()
