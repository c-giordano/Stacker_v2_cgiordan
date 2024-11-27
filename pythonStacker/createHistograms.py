import numpy as np
np.finfo(np.dtype("float32"))
np.finfo(np.dtype("float64"))
# import awkward as ak
import argparse
import sys
import os
import json
import glob

import uproot
import awkward as ak
import src
from src.variables.variableReader import VariableReader, Variable
from src.variables.weightManager import WeightManager
from src.histogramTools import HistogramManager
from src.configuration import load_uncertainties, Channel, Uncertainty
from createEFTWeights import get_eftvariations_filename
from createBSMWeights import get_bsmvariations_filename
import src.arguments as arguments

"""
Script that takes as input a file or set of files, applies cross sections and necessary normalizations (if still needed), and then creates a histogram.
The histogram is taken from settingfiles/Histogramming, so maybe do multiple histograms for the same inputfile? But take range, nbins etc from there as a start.

This is likely code that will be submitted.
Alternatively: create for a single histogram all inputs, so read all files etc.

Must be generic enough: ie. weight variations are easy to support.
More systematic variations is a different question but we can work something out.
"""


def parse_arguments() -> argparse.Namespace:
    """
    For now directly takne from create_histograms in interpretations.
    """
    parser = argparse.ArgumentParser(description='Process command line arguments.')

    # add file arguments
    arguments.add_settingfiles(parser)
    arguments.select_specifics(parser)
    arguments.add_toggles(parser)
    arguments.add_tmp_storage(parser)

    # Parse arguments
    args = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help(sys.stderr)
        exit(1)

    if len(args.years) > 1:
        print("Only accept a single year in CreateHistograms. Exiting...")
        exit(1)

    if args.process is None:
        print("Need to specify a process to produce. Exiting...")
        exit(1)

    if args.channel is None:
        print("Need to specify a channel to produce. Exiting...")
        exit(1)

    return args


def prepare_histogram(data, wgts, variable: Variable):
    hist_content, binning, hist_unc = src.histogram_w_unc_flow(ak.to_numpy(data), range=variable.range, wgts=ak.to_numpy(wgts), nbins=variable.nbins)
    return hist_content, binning, hist_unc


def get_histogram_data(variable: Variable, tree, channel: Channel):
    method = variable.get_method()
    data = method(tree, variable.branch_name, channel.selection)
    return data


def create_histogram(variable: Variable, data, weights):
    '''
    TODO
    '''
    # to histogram:
    hist_content, _, hist_unc = prepare_histogram(data, weights, variable)
    return hist_content, hist_unc


def clean_systematics(systematics: dict[str, Uncertainty], process):
    ret = {name: unc for name, unc in systematics.items() if unc.is_process_relevant(process)}
    return ret


def create_histograms_singledata(output_histograms: dict, args, files, channel: Channel, variables: VariableReader, systematics: dict, globalEFTToggle, globalBSMToggle):
    for filename in files:
        print(f"opening file {filename}")
        processes = ["nonPromptElectron", "nonPromptMuon"] if args.process == "nonPrompt" else [args.process]
        for process in processes:
            print(f"working on process {process}")
            current_tree: uproot.TTree = src.get_tree_from_file(filename, process)
    
            # generates masks for subchannels
            subchannelmasks, subchannelnames = channel.produce_masks(current_tree)
    
            # no structure yet to load systematics weights. Weightmanager? Or move to systematics loop?
            # Does imply more overhead in reading, can try here, see what memory effect it has, otherwise move to systematics.
            # weights = ak.to_numpy(current_tree.arrays(["weights"], cut=channel.selection, aliases={"weights": "nominalWeight"}).weights)
            print("Loading weights...")
            weights = WeightManager(current_tree, channel.selection, systematics)
            if globalEFTToggle:
                eventclass = channel.selection.split("==")[-1]
                weights.add_eftvariations(get_eftvariations_filename(args.storage, filename, eventclass))
            if globalBSMToggle:
                eventclass = channel.selection.split("==")[-1]
                weights.add_bsmvariations(get_bsmvariations_filename(args.storage, filename, eventclass))
            print("Done!")
            for _, variable in variables.get_variable_objects().items():
                if not variable.is_channel_relevant(args.channel):
                    continue
                # load data:
                data = get_histogram_data(variable, current_tree, channel)
    
                for name, syst in systematics.items():
                    if name == "stat_unc":
                        # don't need a dedicated run for this; should be filled before
                        continue
                    if name == "nominal":
                        hist_content, _, hist_unc = prepare_histogram(data, weights["nominal"], variable)
                        output_histograms[args.channel][variable.name]["nominal"] += hist_content
                        output_histograms[args.channel][variable.name]["stat_unc"] += hist_unc
                    else:
                        if not syst.is_process_relevant(process):
                            continue
                        keys = syst.get_weight_keys()
                        hist_content_up, _, _ = prepare_histogram(data, weights[keys[0]], variable)
                        output_histograms[args.channel][variable.name][name]["Up"] += hist_content_up
                        if keys[1] is not None:
                            hist_content_down, _, _ = prepare_histogram(data, weights[keys[1]], variable)
                            output_histograms[args.channel][variable.name][name]["Down"] += hist_content_down
    
                    for subchannel_name in subchannelnames:
                        if name == "stat_unc":
                            # don't need a dedicated run for this; should be filled before
                            continue
                        if name == "nominal":
                            hist_content, _, hist_unc = prepare_histogram(data[subchannelmasks[subchannel_name]], weights["nominal"][subchannelmasks[subchannel_name]], variable)
                            output_histograms[subchannel_name][variable.name]["nominal"] += hist_content
                            output_histograms[subchannel_name][variable.name]["stat_unc"] += hist_unc
                        else:
                            if not syst.is_process_relevant(process):
                                continue
                            keys = syst.get_weight_keys()
                            hist_content_up, _, _ = prepare_histogram(data[subchannelmasks[subchannel_name]], weights[keys[0]][subchannelmasks[subchannel_name]], variable)
                            output_histograms[subchannel_name][variable.name][name]["Up"] += hist_content_up
                            if keys[1] is not None:
                                # TODO: check if this works
                                hist_content_down, _, _ = prepare_histogram(data[subchannelmasks[subchannel_name]], weights[keys[1]][subchannelmasks[subchannel_name]], variable)
                                output_histograms[subchannel_name][variable.name][name]["Down"] += hist_content_down


def create_histogram_shapevar(output_histograms: dict, args, files, channel: Channel, variables: VariableReader, systematics: dict, globalEFTToggle, globalBSMToggle):
    for filename in files:
        name, syst = list(systematics.items())[0]
        for variation in ["Up", "Down"]:
            current_tree: uproot.TTree = src.get_tree_from_file(filename, "Unc_" + name + "_" + variation)
            weights = WeightManager(current_tree, channel.selection, systematics)
            subchannelmasks, subchannelnames = channel.produce_masks(current_tree)

            for _, variable in variables.get_variable_objects().items():
                data = get_histogram_data(variable, current_tree, channel)

                hist_content, _, hist_unc = prepare_histogram(data, weights["nominal"], variable)
                output_histograms[args.channel][variable.name][name][variation] += hist_content

                for subchannel_name in subchannelnames:
                    hist_content, _, _ = prepare_histogram(data[subchannelmasks[subchannel_name]], weights["nominal"][subchannelmasks[subchannel_name]], variable)
                    output_histograms[subchannel_name][variable.name][name][variation] += hist_content


if __name__ == "__main__":
    # parse arguments
    args = parse_arguments()

    # make sure we have correct storagepath
    storagepath = args.storage

    # initialize variable class:
    variables = VariableReader(args.variablefile, args.variable, channel=args.channel)

    # prob do same with systematics
    if args.systematic == "weight" or args.systematic == "shape":
        systematics: dict = load_uncertainties(args.systematicsfile, typefilter=args.systematic, allowflat=False)
    elif args.systematic is not None:
        systematics: dict = load_uncertainties(args.systematicsfile, namefilter=args.systematic, allowflat=False)
    else:
        systematics: dict = dict()

    systematics = clean_systematics(systematics, args.process)
    base_run = args.systematic == "weight" or args.systematic is None

    if args.systematic in systematics:
        base_run = base_run and systematics[args.systematic].type != "shape"

    if base_run:
        systematics["nominal"] = Uncertainty("nominal", {})
        systematics["stat_unc"] = Uncertainty("stat_unc", {})

    print(systematics.keys())
    if len(systematics) == 0:
        print(f"Nothing to do for process {args.process}, syst {args.systematic}, channel {args.channel}.")
        exit(1)

    # load process list:
    with open(args.processfile, 'r') as f:
        processfile = json.load(f)
        processinfo = processfile["Processes"][args.process]
        basedir = processfile["Basedir"]
        subbasedir = basedir.split("/")[-1]
        storagepath = os.path.join(storagepath, subbasedir)

    # prepare channels:
    with open(args.channelfile, 'r') as f:
        # A single job per channel, so we only need to load the current one.
        # This ensures memory usage is refreshed at the end.
        channelinfo = json.load(f)
        channel = Channel(channelinfo[args.channel], channelinfo)
        if channel.isSubchannel:
            print("Current channel is a subchannel. Nothing to be done. Exiting...")
            exit(0)

    # TODO: should this be a global variable?
    globalEFTToggle = args.UseEFT and processinfo.get("hasEFT", 0) > 0 and base_run
    if globalEFTToggle:
        print("USING EFT")
        import plugins.eft as eft
        # load names of EFT Variatioons
        # add to weightmanager -> later on
        # already add to systematics:
        eft_variations = eft.getEFTVariationsGroomed()
        # print(eft_variations)
        for eft_var in eft_variations:
            tmp_dict = {
                "processes": args.process,
                "isEFT": 1
            }
            systematics[eft_var] = Uncertainty(eft_var, tmp_dict)

    globalBSMToggle = args.UseBSM and processinfo.get("hasBSM", 0) > 0 and base_run
    if globalBSMToggle:
        print("USING BSM")
        import plugins.bsm as bsm
        bsm_variations = bsm.getBSMVariationsGroomed()
        for bsm_var in bsm_variations:
            tmp_dict = {
                "processes": args.process,
                "isBSM": 1
            }
            systematics[bsm_var] = Uncertainty(bsm_var, tmp_dict)


    # Now also do this for systematic variations? How though?
    # Maybe add an argument for the systematic variations to produce, can define these somewhere.
    # Should be weightvariations here? Maybe also systematic variations or something? idk let's see first for weight variations

    # update storagepath to include channel
    storagepath_main = os.path.join(storagepath, args.channel)
    output_histograms = dict()
    output_histograms[args.channel] = HistogramManager(storagepath_main, args.process, variables, list(systematics.keys()), year=args.years[0], channel=args.channel)
    for subchannel_name in channel.get_subchannels():
        channel_name = args.channel + subchannel_name
        storagepath_tmp = os.path.join(storagepath, channel_name)
        output_histograms[subchannel_name] = HistogramManager(storagepath_tmp, args.process, variables, list(systematics.keys()), year=args.years[0], channel=channel_name)

    # TODO: get files based on process names -> processmanager can return this, depending on the sys unc?
    filesuffix = "base"
    if not base_run:
        filesuffix = systematics[args.systematic].fileglob
    files = src.get_file_from_globs(basedir, processinfo["fileglobs"], args.years[0], filesuffix)

    if len(files) == 0:
        # make sure nothing is written in files without content to avoid unnecessary saves
        print("No files found for globs. Make sure this is expected!", file=sys.stderr)
        exit(1)
    # print(files)
    # this will not work for requiring a specifc systematic! Need to check.
    if base_run:
        create_histograms_singledata(output_histograms, args, files, channel, variables, systematics, globalEFTToggle, globalBSMToggle)
    else:
        create_histogram_shapevar(output_histograms, args, files, channel, variables, systematics, globalEFTToggle, globalBSMToggle)

    subchannelnames = channel.get_subchannels()
    output_histograms[args.channel].save_histograms()
    for subchannel_name in subchannelnames:
        output_histograms[subchannel_name].save_histograms()
    print("Finished")
    exit(0)
