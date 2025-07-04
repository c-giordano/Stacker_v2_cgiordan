import numpy as np
np.finfo(np.dtype("float32"))
np.finfo(np.dtype("float64"))
import argparse
import sys
import os
import json
import awkward as ak
import src.arguments as arguments
from src.configuration import load_uncertainties, load_channels_and_subchannels, Uncertainty
from src.variables.variableReader import VariableReader, Variable
from src.histogramTools import HistogramManager

"""
Script to reduce all systematics to a single band. Allows for quicker plotting
"""


def parse_arguments() -> argparse.Namespace:
    """
    For now directly takne from create_histograms in interpretations.
    """
    parser = argparse.ArgumentParser(description='Process command line arguments.')

    # add file arguments
    arguments.add_settingfiles(parser)
    arguments.select_specifics(parser)
    arguments.add_tmp_storage(parser)

    # Parse arguments
    args = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help(sys.stderr)
        exit(1)

    return args


def get_systematics_per_process(systematics: dict, processes):
    ret = dict()
    for process in processes:
        ret[process] = {name: unc for name, unc in systematics.items() if unc.is_process_relevant(process)}
    return ret


def get_uncertainty_variation_shape(variable: Variable, uncertainty: Uncertainty, histograms_proc: dict[str, dict[str, HistogramManager]]):
    up = np.zeros(variable.nbins)
    down = np.zeros(variable.nbins)
    if uncertainty.isFlat:
        return up, down

    for process, hists_per_year in histograms_proc.items():
        print(f"\t\t\t{process}")
        if not uncertainty.is_process_relevant(process):
            print("\t\t\tskip process")
            continue

        var_process_up = np.zeros(variable.nbins)
        var_process_down = np.zeros(variable.nbins)

        for year, hists in hists_per_year.items():
            print(f"\t\t\t{year}")
            # get the up and down variation:

            up_diff = np.array(ak.to_numpy(hists[variable.name][uncertainty.name]["Up"] - hists[variable.name]["nominal"]))

            if uncertainty.weight_key_down is None:
                down_diff = np.zeros(len(up_diff))
            else:
                down_diff = np.array(ak.to_numpy(hists[variable.name][uncertainty.name]["Down"] - hists[variable.name]["nominal"]))

            # keep only variations that are truly showing
            true_up_diff = np.where(up_diff > down_diff, up_diff, down_diff)
            true_up_diff[true_up_diff < 0] = 0.

            true_down_diff = np.where(down_diff < up_diff, down_diff, up_diff)
            true_down_diff[true_down_diff > 0] = 0.

            print(f"\t\t\t Process {process} for uncertainty {uncertainty.name} has relative updiff {true_up_diff / hists[variable.name]['nominal']} and downdiff {true_down_diff / hists[variable.name]['nominal']}")

            var_process_up += true_up_diff
            var_process_down -= true_down_diff

        if uncertainty.correlated_process:
            up += var_process_up
            down -= var_process_down
        else:
            up += (var_process_up * var_process_up)
            down += (var_process_down * var_process_down)

    if uncertainty.correlated_process:
        up_sq = np.power(up, 2)
        down_sq = np.power(down, 2)
    else:
        up_sq = up
        down_sq = down

    return up_sq, down_sq


def get_uncertainty_variation_flat(variable: Variable, uncertainty: Uncertainty, histograms_proc: dict[str, dict[str, HistogramManager]]):
    variation = np.zeros(variable.nbins)
    if not uncertainty.isFlat:
        return variation, variation

    for process, hists_per_year in histograms_proc.items():
        if not uncertainty.is_process_relevant(process):
            continue

        var_process = np.zeros(variable.nbins)
        for year, hists in hists_per_year.items():

            # get the up and down variation:
            rate = uncertainty.rate
            if type(rate) is str and "/" in str(rate):
                rate = rate.split("/")[-1]
            rate = float(rate)
            diff = ak.to_numpy(hists[variable.name]["nominal"]) * (1. - rate)
            print(f"\t\t\t {variable.name} has rate {rate} with variations {diff} over {hists[variable.name]['nominal']}")
            var_process += diff

        if uncertainty.correlated_process:
            variation += var_process
        else:
            variation += (var_process * var_process)

    if uncertainty.correlated_process:
        return np.power(variation, 2), np.power(variation, 2)
    else:
        return variation, variation

    return


def uncertaintyloop(variable: Variable, histograms_proc: dict[str, dict[str, HistogramManager]], uncertainties: dict[str, Uncertainty], channel: str):
    # NOTE: can automatically do asymm uncertainties in matplotlib!
    uncertainties_squared_up = np.zeros(variable.nbins)
    uncertainties_squared_down = np.zeros(variable.nbins)
    for name, uncertainty in uncertainties.items():
        print(f"\t\tuncertainty {name} Start -----------------------------------")
        if not uncertainty.is_channel_relevant(channel):
            print("skip channel")
            # skip irrelevant channels
            continue
        if name == "nominal" or name == "stat_unc":
            continue

        if uncertainty.isFlat:
            up_tmp, down_tmp = get_uncertainty_variation_flat(variable, uncertainty, histograms_proc)
        else:
            # For a single uncertainty, loop processes and add up differences to nominal
            up_tmp, down_tmp = get_uncertainty_variation_shape(variable, uncertainty, histograms_proc)

        print(f"\t\t\t Uncertainty band before new contributions: up {np.sqrt(uncertainties_squared_up)}, and down {np.sqrt(uncertainties_squared_down)}")
        uncertainties_squared_up += up_tmp
        uncertainties_squared_down += down_tmp
        print(f"\t\t\t Uncertainty band after new contributions: up {np.sqrt(uncertainties_squared_up)}, and down {np.sqrt(uncertainties_squared_down)}")
        print(f"\t\tuncertainty {name} End -----------------------------------")
    # unsquare
    unc_up = np.sqrt(uncertainties_squared_up)
    unc_down = np.sqrt(uncertainties_squared_down)
    return {"Up": unc_up, "Down": unc_down}


def variableloop(variables: VariableReader, histograms_proc: dict[str, dict[str, HistogramManager]], uncertainties: dict, channel: str):
    ret = dict()
    for var_name, variable in variables.get_variable_objects().items():
        print(f"\t{var_name} Start -----------------------------------")
        if not variable.is_channel_relevant(channel):
            continue
        ret[variable.name] = uncertaintyloop(variable, histograms_proc, uncertainties, channel)
        print(f"\t{var_name} End -----------------------------------")

    return ret


def channelloop(channels: list, variables: VariableReader, systematics_shape: dict, systematics_flat: dict, years: list, processes, storagepath: str, args):
    uncertainties = {**systematics_shape, **systematics_flat}

    unc_per_process = get_systematics_per_process(systematics_shape, list(processes.keys()))

    outputfilename = "total_systematic_"
    if len(years) == 1:
        outputfilename += years[0]
    elif len(years) < 4:
        outputfilename += "_".join(years)
    else:
        outputfilename += "all"
    outputfilename += ".parquet"

    for channelname in channels:
        if args.channel is not None and args.channel not in channelname:
            continue
        print(channelname + "\n--------------------------------------------------------------------")
        # update storagepath to include channel
        storagepath_channel = os.path.join(storagepath, channelname)

        histograms_proc: dict[str, dict[str, HistogramManager]] = dict()
        for process, _ in processes.items():
            histograms_proc[process] = dict()
            for year in years:
                histograms_proc[process][year] = HistogramManager(storagepath_channel, process, variables, list(unc_per_process[process].keys()), year=year, channel=channelname)
                histograms_proc[process][year].load_histograms()

        results = variableloop(variables, histograms_proc, uncertainties, channelname)

        # store results for each variable:
        for variable, var_object in variables.get_variable_objects().items():
            if not var_object.is_channel_relevant(channelname):
                continue
            tmp_path = os.path.join(storagepath_channel, variable, outputfilename)
            ak.to_parquet(ak.Record(results[variable]), tmp_path)

        print("End channel " + channelname + "\n--------------------------------------------------------------------")



if __name__ == "__main__":
    args = parse_arguments()

    # set up histogrammanager for all processes to be plotted. Load in data. Don't care about variables themselves. Whenever it's defined, it's loaded.

    # go over the systematics, add quadratically unless correlated, then linearly
    # for that part load systematics stuff
    systematics_shape: dict = load_uncertainties(args.systematicsfile, allowflat=False)
    systematics_flat: dict = load_uncertainties(args.systematicsfile, typefilter="flat")
    # should contain all. Define "nominal" and "stat".
    systematics_shape["nominal"] = Uncertainty("nominal", {})
    systematics_shape["stat_unc"] = Uncertainty("stat_unc", {})
    # initialize variable class:
    variables = VariableReader(args.variablefile, args.variable)

    # load process list:
    with open(args.processfile, 'r') as f:
        processfile = json.load(f)
        processes = processfile["Processes"]
        basedir = processfile["Basedir"]
        subbasedir = basedir.split("/")[-1]

    storagepath = os.path.join(args.storage, subbasedir)
    # prepare channels:
    channels = load_channels_and_subchannels(args.channelfile)
    channelloop(channels, variables, systematics_shape, systematics_flat, args.years, processes, storagepath, args)
    print("Finished!")
