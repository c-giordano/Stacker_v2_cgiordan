import numpy as np
import uproot
import os
import glob
import itertools

from src.histogramTools import HistogramManager
import plugins.eft as eft
import plugins.bsm as bsm

"""
Some general tools to create histograms
"""


def generate_binning(range, nbins):
    raw_bins = np.linspace(range[0], range[1], nbins + 1)
    return raw_bins


def histogram_w_flow(data, range, wgts, nbins):
    """
    Creates a numpy histogram but takes into account the under- and overflow.
    """

    raw_bins = np.linspace(range[0], range[1], nbins + 1)
    use_bins = [np.array([-np.inf]), raw_bins, np.array([np.inf])]
    use_bins = np.concatenate(use_bins)

    binned_data, _ = np.histogram(data, bins=use_bins, weights=wgts)
    binned_data[1] += binned_data[0]   # add underflow to first bin
    binned_data[-2] += binned_data[-1]  # add overflow to last bin
    binned_data = binned_data[1:-1]     # chop off the under/overflow

    return binned_data, raw_bins


def histogram_w_unc_flow(data, range, wgts, nbins):
    """
    Uses:
        - the data to bin
        - range along which to bin
        - Weight for each data point
        - number of bins
    Returns:
        - The binned data
        - The generated binning
        - The absolute uncertainty on each bin

    Creates a numpy histogram. Takes into account over- and underflow. Calculates the MC uncertainties.
    """

    raw_bins = np.linspace(range[0], range[1], nbins + 1)
    use_bins = [np.array([-np.inf]), raw_bins, np.array([np.inf])]
    use_bins = np.concatenate(use_bins)

    binned_data, _ = np.histogram(data, bins=use_bins, weights=wgts)
    sq_binned_data, _ = np.histogram(data, bins=use_bins, weights=wgts * wgts)
    binned_data[1] += binned_data[0]   # add underflow to first bin
    binned_data[-2] += binned_data[-1]  # add overflow to last bin
    binned_data = binned_data[1:-1]     # chop off the under/overflow

    sq_binned_data[1] += sq_binned_data[0]   # add underflow to first bin
    sq_binned_data[-2] += sq_binned_data[-1]  # add overflow to last bin
    sq_binned_data = sq_binned_data[1:-1]     # chop off the under/overflow
    uncertainty = np.sqrt(sq_binned_data)

    return binned_data, raw_bins, uncertainty


"""
ROOTFile management tools:
"""


def get_file_from_globs(basedir: str, fileglobs: list[str], year: str, suffix: str = "base") -> list:
    files = []
    for filebase in fileglobs:
        fileglob = os.path.join(basedir, filebase)
        fileglob += f"*{year}"
        # if args.systematic == "weight" or args.systematic is None:
        fileglob += f"*{suffix}.root"
        print(fileglob)
        true_files = glob.glob(fileglob)
        files.extend(true_files)
    return files


def get_tree_from_file(filename, processname, fatalfail=False) -> uproot.TTree:
    current_rootfile = uproot.open(filename)

    try:
        current_tree: uproot.TTree = current_rootfile[processname]
    except KeyError:
        print(f"{processname} not found in the file {filename}. Trying other keys.")
        if fatalfail:
            exit(1)
        for key, classname in current_rootfile.classnames().items():
            if classname != "TTree":
                continue
            if "Namingscheme" in key:
                continue
            current_tree: uproot.TTree = current_rootfile[key]
            break

    return current_tree


"""
Tools to call histograms
"""


def load_prepared_histograms(processinfo, channel, variables, systematics, storagepath, args):
    histograms = dict()
    for process, info in processinfo.items():
        histograms[process] = dict()
        if args.UseEFT:
            histograms["TTTT_EFT"] = dict()

            histograms["TTTW_EFT"] = dict()
            histograms["TTTJ_EFT"] = dict()
            histograms["TTT_EFT"] = dict()
        if args.UseBSM:
            for model, mass in itertools.product(args.bsm_model, args.bsm_mass):
                histograms[f"{model}_{mass}"] = dict()
        for year in args.years:
            histograms[process][year] = HistogramManager(storagepath, process, variables, systematics, year, channel=channel)
            histograms[process][year].load_histograms()
            if args.UseEFT:
                systematics_tmp = list(systematics)
                systematics_tmp.extend(eft.getEFTVariationsGroomed())
                histograms["TTTT_EFT"][year] = HistogramManager(storagepath, "TTTT_EFT", variables, systematics_tmp, year, channel=channel)
                histograms["TTTT_EFT"][year].load_histograms()
                histograms["TTTW_EFT"][year] = HistogramManager(storagepath, "TTTW_EFT", variables, systematics_tmp, year, channel=channel)
                histograms["TTTW_EFT"][year].load_histograms()
                histograms["TTTJ_EFT"][year] = HistogramManager(storagepath, "TTTJ_EFT", variables, systematics_tmp, year, channel=channel)
                histograms["TTTJ_EFT"][year].load_histograms()
                histograms["TTT_EFT"][year] = HistogramManager(storagepath, "TTT_EFT", variables, systematics_tmp, year, channel=channel)
                histograms["TTT_EFT"][year].load_histograms()
            if args.UseBSM:
                systematics_tmp = list(systematics)
                systematics_tmp.extend(bsm.getBSMVariationsGroomed())
                for model, mass in itertools.product(args.bsm_model, args.bsm_mass):
                    histograms[f"{model}_{mass}"][year] = HistogramManager(storagepath, f"{model}_M{round(mass/1000)}p{int(mass/100)}", variables, systematics_tmp, year, channel=channel)
                    histograms[f"{model}_{mass}"][year].load_histograms()

    return histograms
