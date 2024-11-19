import numpy as np
np.finfo(np.dtype("float32"))
np.finfo(np.dtype("float64"))
import awkward as ak
import json
import argparse
import matplotlib.pyplot as plt
import os

from src.variables.variableReader import VariableReader, Variable
from src.configuration import load_channels
from src.histogramTools import HistogramManager

import src.plotTools.figureCreator as fg
from src import generate_binning
from plotHistogramsRedo import modify_yrange_updown, generate_outputfolder, copy_index_html

import plugins.eft as eft

import src.arguments as arguments


def parse_arguments():
    parser = argparse.ArgumentParser(description='Process command line arguments.')

    arguments.add_settingfiles(parser)
    arguments.select_specifics(parser)
    arguments.add_tmp_storage(parser)
    arguments.add_plot_output(parser)

    args = parser.parse_args()
    return args




def main_plot_EFT(variable: Variable, plotdir: str, nominal_content, sm_content, process_info: dict, plotlabel: str, processname: str):
    fig, (ax_main, ax_ratio_one) = fg.create_multi_ratioplot(n_subplots=1)#fig, ax_main = fg.create_singleplot()

    binning = generate_binning(variable.range, variable.nbins)
    #print(binning)
    # first plot nominal, then start adding variations
    #nominal_content = np.array(ak.to_numpy(histograms[variable.name]["nominal"]))
    print(nominal_content)
    #stat_unc_var = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name]["stat_unc"])) / nominal_content, nan=0.)
    # pretty_name = generate_process_name("SM", info)
    nominal_weights = np.ones(len(nominal_content))
    ax_main.hist(binning[:-1], binning, weights=nominal_content, histtype="step", color="mediumorchid", label="private simulation (LO)")
    
    #sm_content = np.array(ak.to_numpy(sm_histograms[variable.name]["nominal"]))

    ax_main.hist(binning[:-1], binning, weights=sm_content, histtype="step", color="k", label="central simulation (LO QCD)")

    eft_variations = eft.getEFTVariationsLinear()
    minim = 1.
    maxim = 1.
    minim = min(np.min(sm_content), np.min(nominal_content))
    maxim = max(np.max(sm_content), np.max(nominal_content))
#        pretty_eft_name = eft_var + f" = {wc_factor}"
#        ax_main.hist(binning[:-1], binning, weights=current_variation, histtype="step", label=pretty_eft_name)


    #ax_main.errorbar(x=binning[:-1] + 0.5 * np.diff(binning), y=np.ones(len(nominal_content)), yerr=stat_unc_var, ecolor='k', label="stat unc.")
    ax_main.set_xlim(variable.range)
    ax_main.set_ylabel("Number of Events")
    modify_yrange_updown(ax_main, (minim, maxim), up_scale=1.4)
    ax_main.legend(ncol=1)
    #ax_main.set_xlabel(variable.axis_label)
    ax_main.text(0.049, 0.77, plotlabel, transform=ax_main.transAxes)

    ratio = nominal_content / sm_content
    ax_ratio_one.hist(binning[:-1], binning, weights=nominal_weights, histtype="step", color="black")
    ax_ratio_one.hist(binning[:-1], binning, weights=ratio, histtype="step", color="mediumorchid", label="LO/LO QCD")#deepskyblue lightcoral mediumorchid
    ax_ratio_one.set_xlim(variable.range)
    ax_ratio_one.set_ylabel("LO/LO QCD")
    ax_ratio_one.set_xlabel(variable.axis_label)
    modify_yrange_updown(ax_ratio_one, (0.1,1.2), up_scale=1)


    # fix output name
    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}.png"))
    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}.pdf"))
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
    channels = load_channels(args.channelfile)
    storagepath = os.path.join(args.storage, subbasedir)

    outputfolder_base = generate_outputfolder(args.years, args.outputfolder, subbasedir, suffix="_EFT_Validation")

    # first plot nominal, then start adding variations
    # load variables, want to do this for all processes

    # also load channels

    # contrary to plotHistograms: load uncertainties if no args.UseEFT
    # do need a selection somewhere defined for the uncertainties needed/desired
    for channel in channels:
        print(channel)
        if args.channel is not None and channel != args.channel:
            continue

        storagepath_tmp = os.path.join(storagepath, channel)
        systematics = ["nominal", "stat_unc"]

        systematics.extend(eft.getEFTVariationsGroomed())
        outputfolder = os.path.join(outputfolder_base, channel)
        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)
        copy_index_html(outputfolder)
        
        for _, variable in variables.get_variable_objects().items():
            if not variable.is_channel_relevant(channel):
                continue
            nominal_content = np.zeros((0,)) 
            sm_content = np.zeros((0,))

            first = True
            for year in args.years :
    
                histograms = HistogramManager(storagepath_tmp, args.process, variables, systematics, year)
                histograms.load_histograms()
                
                sm_histograms: dict[str, HistogramManager] = dict()
                sm_histograms = HistogramManager(storagepath_tmp, "TTTJ", variables,systematics, year)
                sm_histograms.load_histograms()
                
                if first : 
                    first = False
                    theShape = np.zeros((len(np.array(ak.to_numpy(histograms[variable.name]["nominal"]))),))
                    nominal_content = np.concatenate([nominal_content, theShape], axis=0)
                    sm_content = np.concatenate([sm_content, theShape], axis=0)
                nominal_content = np.add(nominal_content, np.array(ak.to_numpy(histograms[variable.name]["nominal"])))
                sm_content = np.add(sm_content, np.array(ak.to_numpy(sm_histograms[variable.name]["nominal"])))
    
            main_plot_EFT(variable, outputfolder, nominal_content, sm_content, processinfo, channel, args.process)

#        for subchannel in channels[channel].subchannels.keys():
#            if not variable.is_channel_relevant(channel + subchannel):
#                continue
#            storagepath_tmp = os.path.join(storagepath, channel + subchannel)
#            
#            outputfolder = os.path.join(outputfolder_base, channel, subchannel)
#            if not os.path.exists(outputfolder):
#                os.makedirs(outputfolder)
#            if not os.path.exists(outputfolder):
#                 os.makedirs(outputfolder)
#            copy_index_html(outputfolder)
#            
#            # for process, info in processinfo.items():
#            histograms = HistogramManager(storagepath_tmp, args.process, variables, systematics, args.years[0])
#            histograms.load_histograms()
#            for _, variable in variables.get_variable_objects().items():
#                #lin_quad_plot_EFT(variable, outputfolder, histograms, processinfo, channel + subchannel, args.process)
#                main_plot_EFT(variable, outputfolder, histograms, sm_histograms , processinfo, channel + subchannel, args.process)
#                #mix_plot_EFT(variable, outputfolder, histograms, processinfo, channel + subchannel, args.process)
#                #mix_SM_plot_EFT(variable, outputfolder, histograms, processinfo, channel + subchannel, args.process)
