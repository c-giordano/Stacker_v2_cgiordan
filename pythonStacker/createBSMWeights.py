"""
For files with eft weights, prepare and store the weights per event class, since it takes a lot of computing time.
"""

import plugins.bsm as bsm

import awkward as ak
import argparse
import os
import uproot
import re

import src, fnmatch
import numpy as np
import json
import pandas as pd

from correctionlib.convert import from_uproot_THx
from correctionlib import CorrectionSet
from correctionlib.schemav2 import CorrectionSet as SchemaCorrectionSet

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-f", "--file", dest="inputfile", action="store", required=True,
                        help="Inputfile for which variations must be stored")
    parser.add_argument("-e", "--eventClass", dest="eventclass", action="store", default=11,
                        help="which event class to check")
    parser.add_argument("-p", "--process", dest="process", action="store", default="TTTT_BSM")
    parser.add_argument("--storage", dest="storage", type=str,
                        default="Intermediate", help="Path at which the histograms are stored")
    parser.add_argument("--pseudo", dest='pseudo', action="store_true", default=False,
                        help="activate pseudo reweighting option")

    opts, opts_unknown = parser.parse_known_args()
    return opts


def get_bsmvariations_filename(storage: str, filename: str, eventclass, suffix="") -> str:
    basefilename = filename.split("/")[-1].split(".")[0].split("_localSub")[0]
    if not os.path.exists(os.path.join(storage, "bsmWeights")):
        os.makedirs(os.path.join(storage, "bsmWeights"))
    return os.path.join(storage, "bsmWeights", f"{basefilename}_evClass_{eventclass}_wgts{suffix}.parquet")


def reweight_and_write(reweighter, eventclass, tree, storage, filename, process):
    # should become a record rather than a single array, using the naming scheme I devised in the init
    that_file = "/pnfs/iihe/cms/store/user/cgiordan/AnalysisOutput/ReducedTuples/*/Tree_TTTT_13TeV_LO_TopPhilicScalarSinglet_M0p4_C1p0e00_240611_092611_Run2SIM_UL2018NanoAOD_240612_213908_MCPrompt_2018_base.root"
    if fnmatch.fnmatch(filename, that_file):
        print("SS 0p4 2018, change the BSM variation reading")
        arrs = tree.arrays(["eftVariationsSel"], "eventClass==" + str(eventclass), aliases={"eftVariationsSel": "eftVariations[:, 1:3]"})
    else:
        arrs = tree.arrays(["eftVariationsSel"], "eventClass==" + str(eventclass), aliases={"eftVariationsSel": "eftVariations[:, 2:4]"})
    # var = np.array([reweighter.transform_weights(entry[1:]) for entry in arrs.eftVariationsNorm])
    var = np.transpose(np.array(reweighter.transform_weights(ak.to_numpy(arrs.eftVariationsSel))))
    if len(var) == 0:
        return

    # postprocess: get names, then make dict, then make record
    final_prerecord = dict()
    eft_names = bsm.getBSMVariationsGroomed()
    print(eft_names)
    for i, name in enumerate(eft_names):
        final_prerecord[name] = var[:, i]
    outputfile = get_bsmvariations_filename(storage, filename, eventclass)
    ak.to_parquet(ak.Record(final_prerecord), outputfile)
    return

def write_pseudo_nominal(eventclass,storage, filename, process):
    match = re.search(r"Tree_(TTT[TJW])_.*_TopPhilicScalar(Singlet|Octet)_M(0p\d+|1p[05])", filename)
    if not match:
        raise ValueError(f"Filename '{filename}' does not match the expected pattern.")
    base_name, type_name, mass_value = match.groups()
    root_file_name = f"{base_name}_{mass_value}_{type_name}_c1_ratio_scaled.root"
    print(f"Importing file {root_file_name} for creating parquet")
    root_file_path = f"/user/cgiordan/public_html/GenPseudo/finalPlots_v3/{root_file_name}"
    histogram_name = "h_scaled_ratio"
    try:
        correction = from_uproot_THx(f"{root_file_path}:{histogram_name}")
        correction.data.flow = "clamp" # for overflow/underflow thingy
    except FileNotFoundError:
        raise FileNotFoundError(f"ROOT file not found: {root_file_path}")
        
    correction_set = SchemaCorrectionSet(schema_version=2,corrections=[correction])
    with open(f"{root_file_name}_tmp.json", "w") as f:
        json.dump(correction_set.dict(), f)
    genHt_file = uproot.open(filename)
    matching_trees = [key.strip(";1") for key in genHt_file.keys() if re.match(r"^TTT.*", key)]
    if not matching_trees:
        raise ValueError("No matching trees found in the file.")
    tree_name = matching_trees[0]
    print(f"Found matching tree: {tree_name}")
    genHt_tree = genHt_file[tree_name]
    if "genJetHT" not in genHt_tree.keys():
        raise ValueError(f"The branch 'genJetHT' does not exist in the tree '{tree_name}'.")
    genHt_values = genHt_tree.arrays("genJetHT", "eventClass==" + str(eventclass)).genJetHT
    cset = CorrectionSet.from_file(f"{root_file_name}_tmp.json")
    correction_name = correction.name
    genHt_correction = cset[correction_name]
    weights = genHt_correction.evaluate(genHt_values)
    ak_array = ak.Array(weights)
    df = pd.DataFrame({"weights": weights})

    outputfile = get_bsmvariations_filename(storage, filename, eventclass, suffix="_Pseudo_nominal")
    print(outputfile)
    ak.to_parquet(ak_array, outputfile)#df.to_parquet(outputfile)
    os.remove(f"{root_file_name}_tmp.json")
    print("Pseudo-parquet created!")



if __name__ == "__main__":
    args = parse_arguments()

    reweighter = bsm.bsm_reweighter(2)
    # tree = uproot.open(args.inputfile)

    tree = src.get_tree_from_file(args.inputfile, args.process)

    if (int(args.eventclass) == -1):
        for i in range(15):
            reweight_and_write(reweighter, i, tree, args.storage, args.inputfile, args.process)
    else:
        print(f"running weights for class {args.eventclass}")
        reweight_and_write(reweighter, int(args.eventclass), tree, args.storage, args.inputfile, args.process)
    if args.pseudo:# and args.eventclass in ["3","4","5","6","8","9","10","11","12"]:
        write_pseudo_nominal(int(args.eventclass), args.storage, args.inputfile, args.process)
    print("Success!")
