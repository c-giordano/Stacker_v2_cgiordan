# Standard libraries
import numpy as np
np.finfo(np.dtype("float32"))
np.finfo(np.dtype("float64"))

import awkward as ak
import matplotlib.pyplot as plt

import os
import argparse
import json
import itertools

# Own code
from src import generate_binning, load_prepared_histograms
from src.dataprocessing import DataManager
from src.configuration import load_channels
from src.histogramTools import HistogramManager
from src.variables.variableReader import VariableReader, Variable
import src.arguments as arguments
import src.plotTools.figureCreator as fg

import plugins.eft as eft


def parse_arguments():
    parser = argparse.ArgumentParser(description='Script to plot histograms from the analysis.')

    parser.add_argument("--no_unc", dest="no_unc", action="store_true",
                        default=False, help="Switch to indicate whether the \
                        total uncertainty should be plotted")

    parser.add_argument("--EFTsignal", action="store", default="TTTT_EFT")
    parser.add_argument("--wc", action="store", default="ctt")
    parser.add_argument("--EFT_ratio", dest="EFT_ratio", action="store_true", default=False)
    parser.add_argument("--EFT_fullbkg", dest="EFT_fullbkg", action="store_true", default=False)
    parser.add_argument("--BSM_ratio", dest="BSM_ratio", action="store_true", default=False)
    parser.add_argument("--BSM_fullbkg", dest="BSM_fullbkg", action="store_true", default=False)
    parser.add_argument("--SBRatio", dest="SBRatio", action="store_true", default=False)
    parser.add_argument("--shapes", dest="shapes", action="store_true", default=False)
    parser.add_argument("--suffix", action="store", default="")

    arguments.add_settingfiles(parser)
    arguments.select_specifics(parser)
    arguments.add_tmp_storage(parser)
    arguments.add_plot_output(parser)
    arguments.add_toggles(parser)
    arguments.add_EFT_choice(parser)
    arguments.add_BSM_choices(parser)
    args = parser.parse_args()
    return args


def generate_process_name(name, processinfo):
    pretty_name = processinfo.get("pretty_name", name)
    if "\\" in pretty_name:
        pretty_name = r"${}$".format(pretty_name)
    elif "t" in pretty_name:
        pretty_name = r"${}$".format(pretty_name)
    return pretty_name


def modify_yrange(axis, content, up_scale=1.8):
    cont_max = up_scale * np.max(content)
    cont_min = 0.
    axis.set_ylim((cont_min, cont_max))


def modify_yrange_updown(axis, content, up_scale=0.1, down_scale=0.02):
    cont_max = np.max(content)
    cont_min = np.min(content)
    cont_diff = cont_max - cont_min
    cont_max = cont_max + up_scale * cont_diff
    cont_min = cont_min - down_scale * cont_diff
    axis.set_ylim((cont_min, cont_max))


def copy_index_html(folder):
    if os.getenv("CMSSW_VERSION"):
        outputfile = os.path.join(folder, "index.php")
        if os.path.exists(outputfile):
            return
        os.system(f"cp /user/nivanden/public_html/index.php {outputfile}")


def generate_outputfolder(years, outputfolder, subdir, suffix=""):
    # outputfolder with year:
    outputsubfolder = ""
    if len(years) == 1:
        outputsubfolder += years[0]
    elif len(years) < 4:
        outputsubfolder += "_".join(years)
    else:
        outputsubfolder += "Run2"

    outputsubfolder += suffix

    # outputsubfolder += args.bsm_model
    if args.UseBSM:
        # coupling_name = 
        outputfolder_base = os.path.join(outputfolder, subdir, outputsubfolder, args.bsm_model[0])
        print(outputfolder_base)
    else:
        outputfolder_base = os.path.join(outputfolder, subdir, outputsubfolder)
    if not os.path.exists(outputfolder_base):
        os.makedirs(outputfolder_base)
    copy_index_html(outputfolder_base)
    copy_index_html(os.path.join(outputfolder, subdir))
    return outputfolder_base


def plot_histograms_base(axis, histograms: dict, variable: Variable, processes: dict, years: list, shapes=False):
    signal = np.zeros(variable.nbins)
    bkg = np.zeros(variable.nbins)
    sum_of_content = np.zeros(variable.nbins)

    binning = generate_binning(variable.range, variable.nbins)
    for name, info in processes.items():
        # print(name, info)
        content = np.zeros(variable.nbins)
        for year in years:
            # fixes ordering based on ordering in json file. Honestly sufficient
            content_tmp = ak.to_numpy(histograms[name][year][variable.name]["nominal"])
            content_tmp = np.array(content_tmp)
            content += content_tmp

        # then add to figure, fix label and color, as well as legend
        pretty_name = generate_process_name(name, info)
        if shapes:
            axis.hist(binning[:-1], binning, weights=content, histtype="step", color=info["color"], label=pretty_name)
        else:
            axis.bar(x=binning[:-1], height=content, bottom=sum_of_content, width=np.diff(binning), align='edge', linewidth=0, color=info["color"], zorder=-1, label=pretty_name)

        sum_of_content += content
        if int(info.get("isSignal", 0)) > 0:
            signal += content
        else:
            bkg += content

    modify_yrange(axis, sum_of_content)
    if shapes:
        axis.set_ylabel("a.u.")
    else:
        axis.set_ylabel("Events")

    return {"signal": signal, "bkg": bkg, "sum": sum_of_content, "binning": binning}


def plot_data(axis, datamanager: DataManager, variable: Variable, years: list):
    data = np.zeros(variable.nbins)
    binning = generate_binning(variable.range, variable.nbins)

    # load data or make sure it is loaded?
    content = np.zeros(variable.nbins)
    content, stat_unc = datamanager.get_histogram_and_uncertainties(years, variable)
    stat_unc = np.abs(stat_unc - content)

    # then add to figure, fix label and color, as well as legend
    pretty_name = "Data"
    axis.errorbar(x=binning[:-1] + (np.diff(binning)/2), y=content, yerr=stat_unc, ecolor="#000000", linewidth=1,
                  fmt='o', markersize=5., c="#000000",# mec="#000000",
                  label=pretty_name)



def plot_systematics_band(axis, nominal_content, variable: Variable, storagepath: str, years: list):
    # get binning:
    binning = generate_binning(variable.range, variable.nbins)

    # load content:
    outputfilename = "total_systematic_"
    if len(years) == 1:
        outputfilename += years[0]
    elif len(years) < 4:
        outputfilename += "_".join(years)
    else:
        outputfilename += "all"
    outputfilename += ".parquet"

    path_to_unc = os.path.join(storagepath, variable.name, outputfilename)
    uncertainty_content = ak.from_parquet(path_to_unc)
    unc_up = np.array(ak.to_numpy(uncertainty_content["Up"]))
    unc_down = np.array(ak.to_numpy(uncertainty_content["Down"]))
    axis.bar(x=binning[:-1], height=unc_up + unc_down, bottom=nominal_content - unc_down, width=np.diff(binning),
             align='edge', linewidth=0, edgecolor='#BEC6C4', alpha=1., zorder=-1, hatch="xxxx", fill=False,
             label="Total unc.")
    return unc_up, unc_down


def plot_signal_bkg_ratio(axis, binning, signal, background):
    ratio_content = np.nan_to_num(np.divide(signal, background))
    axis.hist(binning[:-1], binning, weights=ratio_content, histtype="step", color="k")
    modify_yrange_updown(axis, ratio_content)
    return {"ratio": ratio_content}


def plot_EFT_line(axis, histograms, variable: Variable, years, operator: str, normalization_contribution=None):
    lin_name = "EFT_" + operator
    quad_name = lin_name + "_" + operator
    binning = generate_binning(variable.range, variable.nbins)
    nominal_content = np.zeros(variable.nbins)
    if normalization_contribution is None:
        normalization_contribution = np.ones(variable.nbins)

    for year in years:
        nominal_content += np.array(ak.to_numpy(histograms[year][variable.name]["nominal"]))

    all_variations = []
    wc_points = [1. , 2.]
    if (args.wc=="ctHRe" or args.wc=="ctHIm"): 
        wc_points = [20. , 30.]
    for wc_factor in wc_points:
        current_variation = nominal_content # nominal_content
        for year in years:
            current_variation = current_variation + wc_factor * np.array(ak.to_numpy(histograms[year][variable.name][lin_name]["Up"]))
            current_variation = current_variation + wc_factor * wc_factor * np.array(ak.to_numpy(histograms[year][variable.name][quad_name]["Up"]))
        pretty_eft_name = operator + f" = {wc_factor}"

        current_variation = np.nan_to_num(current_variation / normalization_contribution, nan=1., posinf=1., neginf=1.)

        axis.hist(binning[:-1], binning, weights=current_variation, histtype="step",
                  label=pretty_eft_name, linewidth=2.)

        all_variations.append(current_variation)
    return all_variations


def plot_BSM_line(axis, histograms, variable: Variable, years, models: list, masses: list, couplings: list, normalization_contribution=None):
    binning = generate_binning(variable.range, variable.nbins)
    if normalization_contribution is None:
        normalization_contribution = np.ones(variable.nbins)
    # nominal_content = np.zeros(variable.nbins)
    all_variations = []
    for model, mass, coupling in itertools.product(models, masses, couplings):
        name_string = f"{model}_{mass}"
        current_histograms = histograms[name_string]
        # print(histograms.keys())
        current_variation = np.zeros(variable.nbins)
        for year in years:
            current_variation += float(coupling)* float(coupling) * np.array(ak.to_numpy(current_histograms[year][variable.name]["BSM_Quad"]["Up"]))
            current_variation += (float(coupling) ** 4) * np.array(ak.to_numpy(current_histograms[year][variable.name]["BSM_Quartic"]["Up"]))
        
        current_variation = np.nan_to_num(current_variation / normalization_contribution, nan=1., posinf=1., neginf=1.)

        pretty_bsm_name = model.split("Philic")[-1]
        pretty_bsm_name = pretty_bsm_name.replace("Vector", "V")
        pretty_bsm_name = pretty_bsm_name.replace("Singlet", "S")
        pretty_bsm_name = pretty_bsm_name.replace("Octet", "O")
        pretty_bsm_name = pretty_bsm_name.replace("Scalar", "S")
        pretty_bsm_name = pretty_bsm_name.replace("Pseudo", "P")
        mass_name = float(mass/1000)
        coupling_name = "(g="+str(coupling)+")"
        pretty_bsm_name += " " + str(mass_name) + " TeV " + coupling_name
        axis.hist(binning[:-1], binning, weights=current_variation, histtype="step",
                  label=pretty_bsm_name, linewidth=2.)
        all_variations.append(current_variation)
    return all_variations


def finalize_plot(figure, axes, variable: Variable, plotdir: str, plotlabel=""):
    for ax in axes:
        # Not assuming a list with 2 entries ensures variable range doesn't clash with xlim setting 
        ax.set_xlim(variable.range[0], variable.range[-1])

    axes[-1].set_xlabel(variable.axis_label)
    axes[0].text(0.049, 0.76, plotlabel, transform=axes[0].transAxes)
    axes[0].legend(ncols=2)

    figure.savefig(os.path.join(plotdir, f"{variable.name}.png"))
    figure.savefig(os.path.join(plotdir, f"{variable.name}.pdf"))
    print(f"Created figure {os.path.join(plotdir, f'{variable.name}.[png/pdf]')}.")
    plt.close(figure)


def get_lumi(years):
    total_lumi = 0.
    lumi = {
        "2016PreVFP": 19.5,
        "2016PostVFP": 16.8,
        "2016": 36.3,
        "2017": 41.5,
        "2018": 59.8,
        # "Run2": 167.0,
    }
    for year in years:
        total_lumi += lumi[year]

    if total_lumi >= 100.:
        total_lumi = round(total_lumi)
    # total_lumi += lumi[years]
    return total_lumi


def plotting_sequence(args, histograms, variable, processinfo, plotdir, channel, storagepath, datahistograms=None):
    n_ratios = 0
    if args.SBRatio:
        n_ratios = 1
    if args.EFT_ratio or args.EFT_fullbkg:
        n_ratios = 1
    if args.EFT_fullbkg:
        n_ratios = 2
    if args.BSM_ratio:
        n_ratios = 1
    if args.BSM_fullbkg:
        n_ratios = 2
    
    lumi = get_lumi(args.years)
    if n_ratios == 0:
        fig, axes = fg.create_singleplot(lumi, True)
        axes = (axes,)
    elif n_ratios == 1:
        fig, axes = fg.create_ratioplot(lumi, True)
    else:
        fig, axes = fg.create_multi_ratioplot(lumi, True, n_subplots=n_ratios)

    main_plot_out = plot_histograms_base(axes[0], histograms, variable, processinfo, args.years, shapes=args.shapes)
    if datahistograms is not None:    
        data_out = plot_data(axes[0], datahistograms, variable, args.years)
    else:
        data_out = None

    if not args.no_unc:
        plot_systematics_band(axes[0], main_plot_out["sum"], variable, storagepath, args.years)

    if args.UseEFT:
        plot_EFT_line(axes[0], histograms[args.EFTsignal], variable, args.years, args.wc)
    if args.UseBSM:
        plot_BSM_line(axes[0], histograms, variable, args.years, args.bsm_model, args.bsm_mass, args.bsm_coupling)
    if args.BSM_ratio:
        bsm_content = plot_BSM_line(axes[1], histograms, variable, args.years, args.bsm_model, args.bsm_mass, args.bsm_coupling, None)
        axes[1].set_ylabel(r"Sig", fontsize="small")
        modify_yrange_updown(axes[1], bsm_content)

    if args.SBRatio:
        ratiocontent = plot_signal_bkg_ratio(axes[1], main_plot_out["binning"], main_plot_out["signal"], main_plot_out["bkg"])
    if args.EFT_ratio or args.EFT_fullbkg:
        eft_content = plot_EFT_line(axes[1], histograms[args.EFTsignal], variable, args.years, args.wc, main_plot_out["signal"])
        axes[1].set_ylabel(r"EFT / SM $t\bar{t}t\bar{t}$", fontsize="small")
        modify_yrange_updown(axes[1], eft_content)
    if args.EFT_fullbkg:
        eft_content = plot_EFT_line(axes[2], histograms[args.EFTsignal], variable, args.years, args.wc, main_plot_out["sum"])
        axes[2].set_ylabel("EFT / SM", fontsize="small")
        modify_yrange_updown(axes[2], eft_content)
    if args.BSM_fullbkg:
        bsm_content = plot_BSM_line(axes[2], histograms, variable, args.years, args.bsm_model, args.bsm_mass, args.bsm_coupling, main_plot_out["bkg"])
        axes[2].set_ylabel("Sig / BKG", fontsize="small")
        modify_yrange_updown(axes[2], bsm_content)

    finalize_plot(fig, axes, variable, plotdir, channel)


if __name__ == "__main__":
    args = parse_arguments()
    np.seterr(divide='ignore', invalid='ignore')

    # need a set of processes
    with open(args.processfile, 'r') as f:
        processfile = json.load(f)
        processinfo = processfile["Processes"]
        subbasedir = processfile["Basedir"].split("/")[-1]

    # need a set of variables to plot
    variables = VariableReader(args.variablefile, args.variable)

    # need a set of channels to plot
    # prepare channels:
    channels = load_channels(args.channelfile)

    # set up histogrammanager per process
    storagepath = os.path.join(args.storage, subbasedir)

    # TODO: outputfolder with year and suffix:
    outputfolder_base = generate_outputfolder(args.years, args.outputfolder, subbasedir, suffix=args.suffix)

    for channel in channels:
        if args.channel is not None and channel != args.channel:
            continue

        storagepath_tmp = os.path.join(storagepath, channel)
        systematics = ["nominal"]

        outputfolder = os.path.join(outputfolder_base, channel)

        if args.UseEFT : 
            outputfolder = os.path.join(outputfolder_base, args.EFTsignal, args.wc, channel)

        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)
        copy_index_html(outputfolder)

        histograms = load_prepared_histograms(processinfo, channel, variables, systematics, storagepath_tmp, args)
        data_histograms = None
        if args.UseData:
            data_histograms = DataManager(storagepath_tmp, variables, channel, args.years, eras_split=False)

        for variable_name, variable in variables.get_variable_objects().items():
            if not variable.is_channel_relevant(channel):
                continue
            if args.variable and args.variable != variable_name:
                continue
            plotting_sequence(args, histograms, variable, processinfo, outputfolder, channel, storagepath_tmp, data_histograms)
            # exit()
            # plot_variable_base(variable, outputfolder, processinfo, histograms, storagepath=storagepath_tmp, years=args.years, drawEFT=args.UseEFT, no_uncertainty=args.no_unc, plotlabel=channel)

        for subchannel in channels[channel].get_subchannels():
            storagepath_tmp = os.path.join(storagepath, channel + subchannel)
            histograms = dict()

            outputfolder = os.path.join(outputfolder_base, channel, subchannel)

            if args.UseEFT : 
                 outputfolder = os.path.join(outputfolder_base, args.EFTsignal, args.wc, channel, subchannel)

            if not os.path.exists(outputfolder):
                os.makedirs(outputfolder)
            copy_index_html(outputfolder)

            histograms = load_prepared_histograms(processinfo, channel, variables, systematics, storagepath_tmp, args)
            if args.UseData:
                data_histograms = DataManager(storagepath_tmp, variables, channel, args.years, eras_split=False)

            for _, variable in variables.get_variable_objects().items():
                if not variable.is_channel_relevant(channel + subchannel):
                    continue
                plotting_sequence(args, histograms, variable, processinfo, outputfolder, channel+"_"+subchannel, storagepath_tmp, data_histograms)
                # plot_variable_base(variable, outputfolder, processinfo, histograms, storagepath=storagepath_tmp, years=args.years, drawEFT=args.UseEFT, no_uncertainty=args.no_unc, plotlabel=channel+"_"+subchannel)

    print("Finished!")
