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


def lin_quad_plot_EFT(variable: Variable, plotdir: str, histograms, process_info: dict, plotlabel: str, processname: str):
    # first plot nominal, then start adding variations
    print(variable.name)
    nominal_content = np.array(ak.to_numpy(histograms[variable.name]["nominal"]))
    nominal_weights = np.ones(len(nominal_content))

    eft_variations = eft.getEFTVariationsLinear()
    
    # Create lists to store wcRe, wcIm, and Everything_Integrated values
    wcRe_list = []
    wcIm_list = []
    Everything_Integrated_list = []

    for wcRe in range(-50,51):
        ctHRe_array = nominal_content + wcRe * np.array(ak.to_numpy(histograms[variable.name]["EFT_ctHRe"]["Up"])) + wcRe*wcRe* np.array(ak.to_numpy(histograms[variable.name]["EFT_ctHRe_ctHRe"]["Up"]))
        #ctHRe_Integrated = np.sum(ctHRe_array)
        #print ("Events integrated for ctHRe = ",wcRe, " are :" , ctHRe_Integrated)
        for wcIm in range(-50,51):
            ctHIm_array = nominal_content + wcIm * np.array(ak.to_numpy(histograms[variable.name]["EFT_ctHIm"]["Up"])) + wcIm*wcIm* np.array(ak.to_numpy(histograms[variable.name]["EFT_ctHIm_ctHIm"]["Up"]))
            Everything_array = ctHRe_array + ctHIm_array - nominal_content + (2*wcIm*wcRe*np.array(ak.to_numpy(histograms[variable.name]["EFT_ctHRe_ctHIm"]["Up"])))
            Everything_Integrated = np.sum(Everything_array)
            #print("Events integrated for ctHRe = ",wcRe,"and ctHIm = ",wcIm, " are :" , Everything_Integrated)
    
            # Store the values
            wcRe_list.append(wcRe)
            wcIm_list.append(wcIm)
            Everything_Integrated_list.append(Everything_Integrated)
    
    # Convert to numpy arrays
    wcRe_array = np.array(wcRe_list)
    wcIm_array = np.array(wcIm_list)
    Everything_Integrated_array = np.array(Everything_Integrated_list)
    
    # Create 2D plot using scatter
    plt.scatter(wcRe_array, wcIm_array, c=Everything_Integrated_array, cmap='viridis')
    plt.colorbar(label='Integrated xSec')
    plt.xlabel('ctHRe')
    plt.ylabel('ctHIm')
    plt.title('2D Plot of tttt xSec Integrated of ctHRe vs ctHIm')
    plt.grid(True)
    plt.savefig(outputfolder+'/plot_name.png')  # Save the plot as a PNG file
    plt.close()
 
    # Reshape Everything_Integrated_array into a 2D grid if necessary
    grid_size = int(np.sqrt(len(Everything_Integrated_array)))  # Assuming it's a square grid
    Everything_grid = Everything_Integrated_array.reshape((grid_size, grid_size))
    
    # Create a contour plot
    plt.contourf(wcRe_array.reshape((grid_size, grid_size)), wcIm_array.reshape((grid_size, grid_size)), Everything_grid, cmap='viridis')
    plt.colorbar(label='Integrated xSec')
    plt.xlabel('ctHRe')
    plt.ylabel('ctHIm')
    plt.title('2D Contour Plot of Integrated xSec')
    plt.grid(True)
    plt.savefig(outputfolder+'/grid_plot.png')  # Save the plot as a PNG file
    plt.close()
      
    
    
    
    
#    minim = 1.
#    maxim = 1.
#    for eft_var in eft_variations:
#        wc_factor = 1
#        if "ctHRe" in eft_var:
#            wc_factor = 20
#        if "ctHIm" in eft_var:
#            wc_factor = 20
#        lin_name = "EFT_" + eft_var
#        quad_name = lin_name + "_" + eft_var
#
#        current_variation = nominal_content
#        current_variation = current_variation + wc_factor * np.array(ak.to_numpy(histograms[variable.name][lin_name]["Up"]))
#        current_variation = current_variation + wc_factor * wc_factor * np.array(ak.to_numpy(histograms[variable.name][quad_name]["Up"]))
#
#        current_variation = np.nan_to_num(current_variation / nominal_content, nan=1.)
#
#        minim = min(minim, np.min(current_variation))
#        maxim = max(maxim, np.max(current_variation))
#        pretty_eft_name = eft_var + f" = {wc_factor}"
#        ax_main.hist(binning[:-1], binning, weights=current_variation, histtype="step", label=pretty_eft_name)
#
#        lin_ratio = np.nan_to_num(wc_factor * np.array(ak.to_numpy(histograms[variable.name][lin_name]["Up"])) / nominal_content, nan=0.)
#        quad_ratio = np.nan_to_num(wc_factor * wc_factor * np.array(ak.to_numpy(histograms[variable.name][quad_name]["Up"])) / nominal_content, nan=0.)
#        ax_ratio_one.hist(binning[:-1], binning, weights=lin_ratio, histtype="step")
#        ax_ratio_two.hist(binning[:-1], binning, weights=quad_ratio, histtype="step")
#    
#    mix_list = ["cQQ1_cQt1","cQQ1_cQt8","cQQ1_ctHIm","cQQ1_ctHRe","cQQ1_ctt",
#		    "cQQ8_cQQ1","cQQ8_cQt1",
#		    #"cQQ8_cQt8","cQQ8_ctHIm","cQQ8_ctHRe","cQQ8_ctt","cQt1_cQt8","cQt1_ctHIm","cQt1_ctHRe","cQt1_ctt","cQt8_ctHIm","cQt8_ctHRe","ctHRe_ctHIm","ctt_cQt8","ctt_ctHIm","ctt_ctHRe"
#		    ]
#    for eft_var in mix_list:
#        mix_name = "EFT_"+eft_var
#        mix_ratio = np.nan_to_num(np.array(ak.to_numpy(histograms[variable.name][mix_name]["Up"])) / nominal_content, nan=0.)
#        #ax_ratio_two.hist(binning[:-1], binning, weights=mix_ratio, histtype="step")
# 
#
#    ax_main.errorbar(x=binning[:-1] + 0.5 * np.diff(binning), y=np.ones(len(nominal_content)), ecolor='k', label="stat unc.")
#
#    ax_main.set_xlim(variable.range)
#    ax_main.set_ylabel("SM + EFT / SM")
#    modify_yrange_updown(ax_main, (minim, maxim), up_scale=1.2)
#    ax_main.legend(ncol=2)
#    ax_main.text(0.059, 0.74, plotlabel, transform=ax_main.transAxes)
#
#    ax_ratio_one.set_xlim(variable.range)
#    ax_ratio_two.set_xlim(variable.range)
#
#    ax_ratio_one.set_ylabel("Lin / SM")
#    ax_ratio_two.set_ylabel("quad / SM")
#    ax_ratio_two.set_xlabel(variable.axis_label)
#
#    # fix output name
#    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}_ratios.png"))
#    fig.savefig(os.path.join(plotdir, f"{processname}_{variable.name}_ratios.pdf"))
#    plt.close(fig)
#
#


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

    outputfolder_base = "/user/mshoosht/public_html/Interpretations/Plots/Debug/"#generate_outputfolder(args.years, args.outputfolder, subbasedir, suffix="_test")

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

        histograms = HistogramManager(storagepath_tmp, args.process, variables, systematics, args.years[0])
        histograms.load_histograms()

        for _, variable in variables.get_variable_objects().items():
            if not variable.is_channel_relevant(channel):
                continue
#            lin_quad_plot_EFT(variable, outputfolder, histograms, processinfo, channel, args.process)

        for subchannel in channels[channel].subchannels.keys():
            print(subchannel)
            if not variable.is_channel_relevant(channel + subchannel):
                continue
            storagepath_tmp = os.path.join(storagepath, channel + subchannel)
            
            outputfolder = os.path.join(outputfolder_base, channel, subchannel)
            if not os.path.exists(outputfolder):
                os.makedirs(outputfolder)
            if not os.path.exists(outputfolder):
                 os.makedirs(outputfolder)
            copy_index_html(outputfolder)
            
            # for process, info in processinfo.items():
            histograms = HistogramManager(storagepath_tmp, args.process, variables, systematics, args.years[0])
            histograms.load_histograms()
            for _, variable in variables.get_variable_objects().items():
                lin_quad_plot_EFT(variable, outputfolder, histograms, processinfo, channel + subchannel, args.process)
