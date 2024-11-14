"""
For files with eft weights, prepare and store the weights per event class, since it takes a lot of computing time.
"""

import plugins.bsm as bsm

import awkward as ak
import argparse
import os

import src
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-f", "--file", dest="inputfile", action="store", required=True,
                        help="Inputfile for which variations must be stored")
    parser.add_argument("-e", "--eventClass", dest="eventclass", action="store", default=11,
                        help="which event class to check")
    parser.add_argument("-p", "--process", dest="process", action="store", default="TTTT_BSM")
    parser.add_argument("--storage", dest="storage", type=str,
                        default="Intermediate", help="Path at which the histograms are stored")

    opts, opts_unknown = parser.parse_known_args()
    return opts


def get_bsmvariations_filename(storage: str, filename: str, eventclass) -> str:
    basefilename = filename.split("/")[-1].split(".")[0].split("_localSub")[0]
    if not os.path.exists(os.path.join(storage, "bsmWeights")):
        os.makedirs(os.path.join(storage, "bsmWeights"))
    return os.path.join(storage, "bsmWeights", f"{basefilename}_evClass_{eventclass}_wgts.parquet")


def reweight_and_write(reweighter, eventclass, tree, storage, filename, process):
    # should become a record rather than a single array, using the naming scheme I devised in the init
    if filename=="/pnfs/iihe/cms/store/user/cgiordan/AnalysisOutput/ReducedTuples/2024-09-09_09-40/Tree_TTTT_13TeV_LO_TopPhilicScalarSinglet_M0p4_C1p0e00_240611_092611_Run2SIM_UL2018NanoAOD_240612_213908_MCPrompt_2018_base.root":
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

    print("Success!")
