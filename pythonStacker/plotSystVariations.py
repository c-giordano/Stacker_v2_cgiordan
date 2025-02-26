import argparse
import json
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import os

from src.variables.variableReader import VariableReader, Variable
from src.configuration import load_channels, load_channels_and_subchannels
from src.histogramTools import HistogramManager

import src.plotTools.figureCreator as fg
from src import generate_binning
from plotHistograms import modify_yrange_shape, generate_outputfolder, copy_index_html
from src.configuration import load_uncertainties
import src.arguments as arguments


def parse_arguments():
    parser = argparse.ArgumentParser(description='Process command line arguments.')

    arguments.add_settingfiles(parser)
    arguments.select_specifics(parser)
    arguments.add_tmp_storage(parser)
    arguments.add_plot_output(parser)

    args = parser.parse_args()
    return args


def filter_systematics(systematics, process):
    filtered_systematics = {}
    for key, value in systematics.items():
        if not value.is_process_relevant(process):
            continue
        if "PDF_" in key:
            continue
        filtered_systematics[key] = value
    return filtered_systematics


def batch_systematic_keys(systematics, batchsize=4):
    keys = list(systematics.keys())
    batched = []
    i = 0
    while i * batchsize < len(keys):
        batch = []
        if (i + 1) * batchsize > len(keys):
            batch.extend(keys[i * batchsize:])
        else:
            batch.extend(keys[i * batchsize:(i + 1) * batchsize])
        i += 1
        batched.append(batch)
    return batched


def plot_systematicsset(variable: Variable, plotdir: str, histograms, setnb, plotlabel: str, systematics: list):
    fig, ax_main = fg.create_singleplot()

    binning = generate_binning(variable.range, variable.nbins)
    nominal_content = np.array(ak.to_numpy(histograms[variable.name]["nominal"]))
    stat_unc_var = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name]["stat_unc"])) / nominal_content, nan=0., posinf=1., neginf=1.)
    nominal_weights = np.ones(len(nominal_content))
    ax_main.hist(binning[:-1], binning, weights=nominal_weights, histtype="step", color="k", label="SM")

    minim = 1.
    maxim = 1.
    colors = ["r", "g", "b", "y", "c", "m"]
    for i, syst in enumerate(systematics):
        # TODO: actually load the variation
        upvar = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][syst]["Up"])) / nominal_content, nan=-1., posinf=-1., neginf=-1.)
        downvar = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][syst]["Down"])) / nominal_content, nan=-1., posinf=-1., neginf=-1.)

        minim = min(minim, min(np.min(downvar), np.min(upvar)))
        maxim = max(maxim, max(np.max(downvar), np.max(upvar)))

        # plot
        ax_main.hist(binning[:-1], binning, weights=upvar, histtype="step", label=f"{syst}", color=colors[i])
        ax_main.hist(binning[:-1], binning, weights=downvar, histtype="step", color=colors[i], linestyle="dashed")

    ax_main.errorbar(x=binning[:-1] + 0.5 * np.diff(binning), y=np.ones(len(nominal_content)), yerr=np.abs(stat_unc_var), ecolor='k', label="stat unc.")
    ax_main.set_xlim((variable.range[0], variable.range[-1]))
    ax_main.set_ylabel("Unc / nom.")
    # modify_yrange_shape((minim, maxim), ax_main, minscale=0.95, maxscale=1.4)
    diff_minmax = maxim - minim
    ax_main.set_ylim((minim - 0.1 * diff_minmax, maxim + 0.8 * diff_minmax))
    # print((minim, maxim))
    # print((minim - 0.02 * diff_minmax, maxim + 0.4 * diff_minmax))

    ax_main.legend(ncol=1)
    ax_main.set_xlabel(variable.axis_label)
    ax_main.text(0.049, 0.77, plotlabel, transform=ax_main.transAxes)

    # fix output name
    print(f"Created figure {os.path.join(plotdir, f'{variable.name}_{setnb}.png')}")
    fig.savefig(os.path.join(plotdir, f"{variable.name}_{setnb}.png"))
    fig.savefig(os.path.join(plotdir, f"{variable.name}_{setnb}.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    args = parse_arguments()
    np.seterr(divide='ignore', invalid='ignore')

    # load process specifics
    # need a set of processes
    with open(args.processfile, 'r') as f:
        processfile = json.load(f)
        processinfo = processfile["Processes"][args.process]
        subbasedir = processfile["Basedir"].split("/")[-1]

    variables = VariableReader(args.variablefile, args.variable)
    channels = load_channels_and_subchannels(args.channelfile)
    storagepath = os.path.join(args.storage, subbasedir)

    outputfolder_base = generate_outputfolder(args.years, args.outputfolder, subbasedir, suffix="_Syst_Variations")

    # load systematics
    systematics = load_uncertainties(args.systematicsfile, allowflat=False)
    # batched systematics
    filtered_systematics = filter_systematics(systematics, args.process)
    batched_systematics = batch_systematic_keys(filtered_systematics)

    systematics_list = list(systematics.keys()) + ["nominal", "stat_unc"]
    for channel in channels:
        if args.channel is not None and not args.channel in channel:
            continue
        print(channel)
        storagepath_tmp = os.path.join(storagepath, channel)

        outputfolder = os.path.join(outputfolder_base, channel, args.process)
        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)
        copy_index_html(outputfolder)
        copy_index_html("/".join(outputfolder.split("/")[:-1]))

        histograms = HistogramManager(storagepath_tmp, args.process, variables, systematics_list, args.years[0])
        histograms.load_histograms()

        for _, variable in variables.get_variable_objects().items():
            if not variable.is_channel_relevant(channel):
                continue
            for i, syst_set in enumerate(batched_systematics):
                plot_systematicsset(variable, outputfolder, histograms, i, channel, syst_set)
        

        # for subchannel in channels[channel].subchannels.keys():
        #     print(subchannel)
        #     storagepath_tmp = os.path.join(storagepath, channel + subchannel)
        #     outputfolder = os.path.join(outputfolder_base, channel, args.process, subchannel)
        #     if not os.path.exists(outputfolder):
        #         os.makedirs(outputfolder)
        #     copy_index_html(outputfolder)
        # 
        #     # histograms = HistogramManager(storagepath_tmp, args.process, variables, systematics_list, args.years[0])
        #     # histograms.load_histograms()
        #     # for _, variable in variables.get_variable_objects().items():
        #     #     if not variable.is_channel_relevant(channel):
        #     #         continue
        #     #     for i, syst_set in enumerate(batched_systematics):
        #     #         plot_systematicsset(variable, outputfolder, histograms, i, channel + subchannel, syst_set)
    
    print("Finished!")
