import numpy as np
np.finfo(np.dtype("float32"))
np.finfo(np.dtype("float64"))
import argparse
import sys
import os
import json
import uproot
import ROOT
import awkward as ak
import itertools

from src.histogramTools import HistogramManager
from src.variables.variableReader import VariableReader, Variable
from src import generate_binning
from src.configuration import load_channels_and_subchannels, load_uncertainties, Uncertainty
from src.configuration.ScaleVariation import make_envelope
import src.histogramTools.converters as cnvrt
from src.datacardTools import DatacardWriter

import src.arguments as arguments
"""
Script to generate a bunch of datacards.
If the inputhistograms don't exist, it will create these.
Need a writer class I guess, needs to 
"""


def parse_arguments():
    parser = argparse.ArgumentParser(description='Process some integers.')

    arguments.add_settingfiles(parser)
    arguments.add_tmp_storage(parser)
    arguments.add_toggles(parser)
    arguments.add_EFT_choice(parser)
    # args_select_specifics(parser)

    parser.add_argument("-dcf", "--datacardfile", dest="datacardfile",
                        type=str, help="JSON file with info to create the datacards")
    parser.add_argument("-op", "--outputpath", dest="outputpath",
                        type=str, help="Path to outputfolder for DC.",
                        default="output/datacards/")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    return args


def convert_and_write_histogram(input_histogram, variable: Variable, outputname: str, rootfile: uproot.WritableDirectory, statunc=None):
    raw_bins = generate_binning(variable.range, variable.nbins)
    if statunc is None:
        statunc = np.zeros(len(input_histogram))

    # safety functions:
    ret_th1: ROOT.TH1D = cnvrt.numpy_to_TH1D(input_histogram, raw_bins, err=statunc)
    for i in range(1, ret_th1.GetNbinsX() + 1):
        if ret_th1.GetBinContent(i) > 0.001:
            continue
        if ret_th1.GetBinContent(i) < -0.5:
            print(f"WARNING: Significant negative value in {outputname} bin {i}! Setting to 0.001")

        ret_th1.SetBinError(i, 0.001)
        ret_th1.SetBinContent(i, 0.001)
    rootfile[outputname] = ret_th1


def get_pretty_channelnames(dc_settings):
    pretty_names = [channel_DC_setting["prettyname"] for _, channel_DC_setting in dc_settings["channelcontent"].items()]
    return pretty_names


def patch_scalevar_correlations(systematics, processes):
    if not "ScaleVarEnvelope" in systematics:
        return

    envelope_relevant_mod = list(systematics["ScaleVarEnvelope"].processes)
    if "TTTT" in envelope_relevant_mod:
        for process in processes:
            if not ("sm" in process or "quad" in process):
                continue
            if process == "sm":
                continue
            envelope_relevant_mod.append(process)

    systematics["ScaleVarEnvelope"].processes = envelope_relevant_mod
    envelope_relevant = list(systematics["ScaleVarEnvelope"].processes)
    envelope_relevant += ["nonPromptElectron", "nonPromptMuon", "ChargeMisID"]
    for process in processes:
        if (("sm" in process) or ("quad" in process) or ("TTTT" in process)): envelope_relevant.append(process)

    relevant_processes = [process for process in processes if not process in envelope_relevant]
    
    #print ("\n","list of envelope relevant stuf:", envelope_relevant)
    #print ("\n","list of all relevant stuff:", relevant_processes)

    for i in range(6):
        key = f"ScaleVar_{i}"
        systematics[key].processes = relevant_processes



def nominal_datacard_creation(rootfile: uproot.WritableDirectory, datacard_settings: dict, channels: dict, processes: list, shape_systematics: dict, args: argparse.Namespace):
    """
    Code to create nominal datacard content (ie the really SM processes).
    For SM stuff, this is the only thing that should be used.
    """
    all_asimovdata = dict()
    for channelname, channel_DC_setting in datacard_settings["channelcontent"].items():
        # load histograms for this specific channel and the variable with HistogramManager
        histograms = dict()

        # asimov set
        # setup variable reader for the single variable
        var_name = channel_DC_setting["variable"]
        variables = VariableReader(args.variablefile, [var_name])
        storagepath = os.path.join(args.storage, channelname)
        # setup systematics for current channel

        asimov_data = np.zeros(variables.get_properties(var_name).nbins)

        for process in processes:
            if channels[channelname].is_process_excluded(process):
                continue
            histograms = HistogramManager(storagepath, process, variables, list(shape_systematics.keys()), args.years[0])
            histograms.load_histograms()

            # write nominal
            asimov_data += np.array(ak.to_numpy(histograms[var_name]["nominal"]))
            path_to_histogram = f"{channel_DC_setting['prettyname']}/{process}"
            convert_and_write_histogram(histograms[var_name]["nominal"], variables.get_properties(var_name), path_to_histogram, rootfile, statunc=histograms[var_name]["stat_unc"])

            # loop and write systematics
            for systname, syst in shape_systematics.items():
                if systname == "nominal" or systname == "stat_unc":
                    continue
                if not syst.is_process_relevant(process):
                    continue

                upvar = histograms[var_name][systname]["Up"]
                if syst.weight_key_down is None:
                    downvar = histograms[var_name]["nominal"]
                else:
                    downvar = histograms[var_name][systname]["Down"]
                if systname == "ScaleVarEnvelope":
                    upvar, downvar = make_envelope(histograms[var_name])

                rootpath_systname = syst.technical_name
                if not syst.correlated_process:
                    rootpath_systname += process

                path_to_histogram_systematic_up = f"{channel_DC_setting['prettyname']}/{rootpath_systname}Up/{process}"
                path_to_histogram_systematic_down = f"{channel_DC_setting['prettyname']}/{rootpath_systname}Down/{process}"
                convert_and_write_histogram(upvar, variables.get_properties(var_name), path_to_histogram_systematic_up, rootfile)
                convert_and_write_histogram(downvar, variables.get_properties(var_name), path_to_histogram_systematic_down, rootfile)

        # write data_obs:
        all_asimovdata[channel_DC_setting['prettyname']] = asimov_data
        # asimov_data_path = f"{channel_DC_setting['prettyname']}/data_obs"
        # convert_and_write_histogram(asimov_data, variables.get_properties(var_name), asimov_data_path, rootfile, statunc=np.sqrt(asimov_data))
    return all_asimovdata


def eft_datacard_creation(rootfile: uproot.WritableDirectory, datacard_settings: dict, eft_variations: list, shape_systematics: dict, args: argparse.Namespace):
    """
    return should be a list of the added "processes" in the datacard.
    """

    print (eft_variations)

    ret = [["sm", 0]]
    # do some stuff
    # first add sm to the list and load those histograms -> normal histogrammanager
    # load SM stuff:
    sm_histograms: dict[str, HistogramManager] = dict()
    all_asimovdata = dict()
    for channelname, channel_DC_setting in datacard_settings["channelcontent"].items():
        storagepath = os.path.join(args.storage, channelname)

        var_name = channel_DC_setting["variable"]

        variables = VariableReader(args.variablefile, [var_name])

        sm_histograms[channelname] = HistogramManager(storagepath, "TTTT", variables, list(shape_systematics.keys()), args.years[0])
        sm_histograms[channelname].load_histograms()
        all_asimovdata[channel_DC_setting['prettyname']] = sm_histograms[channelname][var_name]["nominal"]

        path_to_histogram = f"{channel_DC_setting['prettyname']}/sm"
        convert_and_write_histogram(sm_histograms[channelname][var_name]["nominal"], variables.get_properties(var_name), path_to_histogram, rootfile, statunc=sm_histograms[channelname][var_name]["stat_unc"])
        # loop and write systematics
        for systname, syst in shape_systematics.items():
            if systname == "nominal" or systname == "stat_unc":
                continue
            if not syst.is_process_relevant("TTTT"):
                continue

            upvar = sm_histograms[channelname][var_name][systname]["Up"]
            if syst.weight_key_down is None:
                downvar = sm_histograms[channelname][var_name]["nominal"]
            else:
                downvar = sm_histograms[channelname][var_name][systname]["Down"]
            if systname == "ScaleVarEnvelopeTTTT":
                upvar, downvar = make_envelope(sm_histograms[channelname][var_name])

            rootpath_systname = syst.technical_name
            if not syst.correlated_process:
                rootpath_systname += "TTTT"

            path_to_histogram_systematic_up = f"{channel_DC_setting['prettyname']}/{rootpath_systname}Up/sm"
            path_to_histogram_systematic_down = f"{channel_DC_setting['prettyname']}/{rootpath_systname}Down/sm"
            convert_and_write_histogram(upvar, variables.get_properties(var_name), path_to_histogram_systematic_up, rootfile)
            convert_and_write_histogram(downvar, variables.get_properties(var_name), path_to_histogram_systematic_down, rootfile)

    # next lin and quad terms per eft variation
    counter = 0
    for eft_var in eft_variations:
        ret.append([f"sm_lin_quad_{eft_var}", -(1+counter)])
        ret.append([f"quad_{eft_var}", -(2+counter)])
        counter +=2
        lin_name = "EFT_" + eft_var
        quad_name = "EFT_" + eft_var + "_" + eft_var
        for channelname, channel_DC_setting in datacard_settings["channelcontent"].items():
            # print(channelname)
            storagepath = os.path.join(args.storage, channelname)
            # write nominal part

            var_name = channel_DC_setting["variable"]
            variables = VariableReader(args.variablefile, [var_name])

            # load nominal and stat unc? stat unc from EFT sample itself, then nominal is the variation, so:
            content_to_load = ["nominal", "stat_unc", lin_name, quad_name]
            histograms_eft = HistogramManager(storagepath, "TTTT_EFT", variables, content_to_load, args.years[0])
            histograms_eft.load_histograms()
            # print(channelname)
            # print(var_name)
            
            rel = np.nan_to_num(sm_histograms[channelname][var_name]["nominal"] / histograms_eft[var_name]["nominal"], nan=1.)
            rel = np.where(np.abs(rel) > 1e10, 1., rel)

            content_sm_lin_quad_nominal_beforeRW = histograms_eft[var_name]["nominal"] + histograms_eft[var_name][lin_name]["Up"] + histograms_eft[var_name][quad_name]["Up"]
            content_sm_lin_quad_nominal = rel * content_sm_lin_quad_nominal_beforeRW
            statunc_sm_lin_quad_nominal = np.nan_to_num(ak.to_numpy(content_sm_lin_quad_nominal * histograms_eft[var_name]["stat_unc"] / histograms_eft[var_name]["nominal"]))
            path_to_sm_lin_quad = f"{channel_DC_setting['prettyname']}/sm_lin_quad_{eft_var}"
            
            content_quad_nominal_beforeRW = histograms_eft[var_name][quad_name]["Up"]
            content_quad_nominal = rel * content_quad_nominal_beforeRW
            statunc_quad_nominal = np.nan_to_num(ak.to_numpy(content_quad_nominal * histograms_eft[var_name]["stat_unc"] / histograms_eft[var_name]["nominal"]))
            path_to_quad = f"{channel_DC_setting['prettyname']}/quad_{eft_var}"
            
            
            convert_and_write_histogram(content_sm_lin_quad_nominal, variables.get_properties(var_name), path_to_sm_lin_quad, rootfile, statunc=statunc_sm_lin_quad_nominal)
            convert_and_write_histogram(content_quad_nominal, variables.get_properties(var_name), path_to_quad, rootfile, statunc=statunc_quad_nominal)


	    # print(content_sm_lin_quad_nominal)
            # loop and write systematics
            for systname, syst in shape_systematics.items():
                #print("\n",f"running systematic {systname}")
                if systname == "nominal" or systname == "stat_unc":
                    continue
                if not syst.is_process_relevant("TTTT"):
                    continue

                upvar_sm = sm_histograms[channelname][var_name][systname]["Up"]
                if syst.weight_key_down is None:
                    downvar_sm = sm_histograms[channelname][var_name]["nominal"]
                else:
                    downvar_sm = sm_histograms[channelname][var_name][systname]["Down"]

                if systname == "ScaleVarEnvelopeTTTT":
                    upvar_sm, downvar_sm = make_envelope(sm_histograms[channelname][var_name])

                rel_syst_up = np.nan_to_num(upvar_sm / histograms_eft[var_name]["nominal"], nan=1.)
                rel_syst_up = np.where(np.abs(rel_syst_up) > 1e10, 1., rel_syst_up)

                content_sm_lin_quad_syst_up = rel_syst_up * content_sm_lin_quad_nominal_beforeRW
                content_quad_syst_up = rel_syst_up * content_quad_nominal_beforeRW

#                if syst.weight_key_down is None:
#                    content_sm_lin_quad_syst_down = content_sm_lin_quad_nominal_beforeRW
#                    content_quad_syst_down = content_quad_nominal_beforeRW
#                else:
                rel_syst_down = np.nan_to_num(downvar_sm / histograms_eft[var_name]["nominal"], nan=1.)
                rel_syst_down = np.where(np.abs(rel_syst_down) > 1e10, 1., rel_syst_down)

                content_sm_lin_quad_syst_down = rel_syst_down * content_sm_lin_quad_nominal_beforeRW
                content_quad_syst_down = rel_syst_down * content_quad_nominal_beforeRW

                rootpath_smlinquad = syst.technical_name
                rootpath_quad = syst.technical_name
                if not syst.correlated_process:
                    rootpath_smlinquad += f"sm_lin_quad_{eft_var}"
                    rootpath_quad += f"quad_{eft_var}"
 
                path_to_sm_lin_quad_syst_up = f"{channel_DC_setting['prettyname']}/{rootpath_smlinquad}Up/sm_lin_quad_{eft_var}"
                path_to_sm_lin_quad_syst_down = f"{channel_DC_setting['prettyname']}/{rootpath_smlinquad}Down/sm_lin_quad_{eft_var}"
                #print("\n","I am here. I am doing this : ",path_to_sm_lin_quad_syst_up)
                convert_and_write_histogram(content_sm_lin_quad_syst_up, variables.get_properties(var_name), path_to_sm_lin_quad_syst_up, rootfile)
                convert_and_write_histogram(content_sm_lin_quad_syst_down, variables.get_properties(var_name), path_to_sm_lin_quad_syst_down, rootfile)
 
                path_to_quad_syst_up = f"{channel_DC_setting['prettyname']}/{rootpath_quad}Up/quad_{eft_var}"
                path_to_quad_syst_down = f"{channel_DC_setting['prettyname']}/{rootpath_quad}Down/quad_{eft_var}"
                convert_and_write_histogram(content_quad_syst_up, variables.get_properties(var_name), path_to_quad_syst_up, rootfile)
                convert_and_write_histogram(content_quad_syst_down, variables.get_properties(var_name), path_to_quad_syst_down, rootfile)


    if len(eft_variations)>=2 :
        mix_list = ["cQQ1_cQt1","cQQ1_cQt8","cQQ1_ctHIm","cQQ1_ctHRe","cQQ1_ctt","cQQ8_cQQ1","cQQ8_cQt1","cQQ8_cQt8","cQQ8_ctHIm","cQQ8_ctHRe","cQQ8_ctt","cQt1_cQt8","cQt1_ctHIm","cQt1_ctHRe","cQt1_ctt","cQt8_ctHIm","cQt8_ctHRe","ctHRe_ctHIm","ctt_cQt8","ctt_ctHIm","ctt_ctHRe"]
        combinations = list(itertools.combinations(eft_variations, 2))
        for combi in combinations:
            combistring = "_".join(combi)
            if not combistring in mix_list:
                combistring = "_".join(reversed(combi))
                if not combistring in mix_list:
                    print(f"ERROR: combination {combi} not found in mix_list")
                    exit(1)
            eft_var = combistring

            ret.append([f"sm_lin_quad_mixed_{eft_var}", -(1+counter)])
            counter +=1
            mixName = "EFT_"+eft_var
            for channelname, channel_DC_setting in datacard_settings["channelcontent"].items():
                storagepath = os.path.join(args.storage, channelname)
                var_name = channel_DC_setting["variable"]
                variables = VariableReader(args.variablefile, [var_name])
                
                # load nominal and stat unc? stat unc from EFT sample itself, then nominal is the variation, so:
                content_to_load = ["nominal", "stat_unc", mixName]
                histograms_eft = HistogramManager(storagepath, "TTTT_EFT", variables, content_to_load, args.years[0])
                histograms_eft.load_histograms()
                
                rel = np.nan_to_num(sm_histograms[channelname][var_name]["nominal"] / histograms_eft[var_name]["nominal"], nan=1.)
                rel = np.where(np.abs(rel) > 1e10, 1., rel)
               
                content_mix_nominal_beforeRW = histograms_eft[var_name][mixName]["Up"] 
                content_mix_nominal = rel * content_mix_nominal_beforeRW
                statunc_mix_nominal = np.nan_to_num(ak.to_numpy(content_mix_nominal * histograms_eft[var_name]["stat_unc"] / histograms_eft[var_name]["nominal"]))
                path_to_mix = f"{channel_DC_setting['prettyname']}/sm_lin_quad_mixed_{eft_var}"
            
                convert_and_write_histogram(content_mix_nominal, variables.get_properties(var_name), path_to_mix, rootfile, statunc=statunc_mix_nominal)
                
                # loop and write systematics
                for systname, syst in shape_systematics.items():
                    if systname == "nominal" or systname == "stat_unc":
                            continue
                    if not syst.is_process_relevant("TTTT"):
                            continue
                    
                    upvar_sm = sm_histograms[channelname][var_name][systname]["Up"]
                    if syst.weight_key_down is None:
                            downvar_sm = sm_histograms[channelname][var_name]["nominal"]
                    else:
                            downvar_sm = sm_histograms[channelname][var_name][systname]["Down"]
                    
                    if systname == "ScaleVarEnvelopeTTTT":
                            upvar_sm, downvar_sm = make_envelope(sm_histograms[channelname][var_name])
                    
                    rel_syst_up = np.nan_to_num(upvar_sm / histograms_eft[var_name]["nominal"], nan=1.)
                    rel_syst_up = np.where(np.abs(rel_syst_up) > 1e10, 1., rel_syst_up)
                    
                    content_mix_syst_up = rel_syst_up * content_mix_nominal_beforeRW
                    
#                    if syst.weight_key_down is None:
#                        content_mix_syst_down = content_mix_nominal
#                    else:
                    rel_syst_down = np.nan_to_num(downvar_sm / histograms_eft[var_name]["nominal"], nan=1.)
                    rel_syst_down = np.where(np.abs(rel_syst_down) > 1e10, 1., rel_syst_down)
                    
                    content_mix_syst_down = rel_syst_down * content_mix_nominal_beforeRW
                    
                    rootpath_mix = syst.technical_name
                    if not syst.correlated_process:
                            rootpath_mix += f"sm_lin_quad_mixed_{eft_var}"
                    
                    path_to_mix_syst_up = f"{channel_DC_setting['prettyname']}/{rootpath_mix}Up/sm_lin_quad_mixed_{eft_var}"
                    path_to_mix_syst_down = f"{channel_DC_setting['prettyname']}/{rootpath_mix}Down/sm_lin_quad_mixed_{eft_var}"
                    # print(path_to_mix_syst_up)
                    convert_and_write_histogram(content_mix_syst_up, variables.get_properties(var_name), path_to_mix_syst_up, rootfile)
                    convert_and_write_histogram(content_mix_syst_down, variables.get_properties(var_name), path_to_mix_syst_down, rootfile)



    return ret, all_asimovdata


def bsm_datacard_creation(rootfile: uproot.WritableDirectory, datacard_settings: dict, bsm_variations: list, shape_systematics: dict, bsm_process: dict, args: argparse.Namespace):
    ret = []
    bsm_process_name = list(bsm_process.keys())[0]
    bsm_process_info = bsm_process[bsm_process_name]

    for channel_num, (channelname, channel_DC_setting) in enumerate(datacard_settings["channelcontent"].items()):
        storagepath = os.path.join(args.storage, channelname)
        var_name = channel_DC_setting["variable"]
        variables = VariableReader(args.variablefile, [var_name])

        ### Now we need the histograms:
        # Load nominal, stat_unc, BSM variations and systematic variations:
        content_to_load = ["nominal", "stat_unc"]
        content_to_load += bsm_variations
        content_to_load += list(shape_systematics.keys())
        histograms_bsm = HistogramManager(storagepath, bsm_process_name, variables, content_to_load, args.years[0])
        histograms_bsm.load_histograms()

        # relative size of stat unc:
        stat_unc_rel = np.nan_to_num(histograms_bsm[var_name]["stat_unc"] / histograms_bsm[var_name]["nominal"])

        ### Next, loop the existing variations.
        for count, bsm_variation in enumerate(bsm_variations):
            if channel_num == 0:
                ret.append([bsm_variation, -count])

            # For each variation, first write the nominal component to file:
            path_to_histogram = f"{channel_DC_setting['prettyname']}/{bsm_variation}"
            # recalc stat unc:
            stat_unc = stat_unc_rel * histograms_bsm[var_name][bsm_variation]["Up"]
            convert_and_write_histogram(histograms_bsm[var_name][bsm_variation]["Up"], variables.get_properties(var_name), path_to_histogram, rootfile, statunc=stat_unc)

            # Loop systematics:
            for systname, syst in shape_systematics.items():
                if systname == "nominal" or systname == "stat_unc":
                    continue
                if not syst.is_process_relevant(bsm_process_name):
                    continue

                upvar = histograms_bsm[var_name][systname]["Up"]
                if syst.weight_key_down is None:
                    downvar = histograms_bsm[var_name]["nominal"]
                else:
                    downvar = histograms_bsm[var_name][systname]["Down"]
                if systname == "ScaleVarEnvelope":
                    upvar, downvar = make_envelope(histograms_bsm[var_name])
                
                # make them relative so that we can recast these to variations on the nominal contribution:
                rel_syst_up = np.nan_to_num(upvar / histograms_bsm[var_name]["nominal"], nan=1.)
                rel_syst_up = np.where(np.abs(rel_syst_up) > 1e10, 1., rel_syst_up)

                content_var_syst_up = rel_syst_up * histograms_bsm[var_name][bsm_variation]["Up"]

                if syst.weight_key_down is None:
                    content_var_syst_down = histograms_bsm[var_name][bsm_variation]["Up"]
                else:
                    rel_syst_down = np.nan_to_num(downvar / histograms_bsm[var_name]["nominal"], nan=1.)
                    rel_syst_down = np.where(np.abs(rel_syst_down) > 1e10, 1., rel_syst_down)

                    content_var_syst_down = rel_syst_down * histograms_bsm[var_name][bsm_variation]["Up"]

                rootpath_systname = syst.technical_name
                if not syst.correlated_process:
                    rootpath_systname += bsm_variation

                path_to_histogram_systematic_up = f"{channel_DC_setting['prettyname']}/{rootpath_systname}Up/{bsm_variation}"
                path_to_histogram_systematic_down = f"{channel_DC_setting['prettyname']}/{rootpath_systname}Down/{bsm_variation}"
                convert_and_write_histogram(upvar, variables.get_properties(var_name), path_to_histogram_systematic_up, rootfile)
                convert_and_write_histogram(downvar, variables.get_properties(var_name), path_to_histogram_systematic_down, rootfile)
    return ret

if __name__ == "__main__":
    np.seterr(divide='ignore', invalid='ignore')

    args = parse_arguments()

    # TODO: for systematics, add a year filter or something, so that we don't introduce 10 different config files.
    with open(args.datacardfile, 'r') as f:
        datacard_settings = json.load(f)

    # Load channels
    channels = load_channels_and_subchannels(args.channelfile)

    # Load all processes:
    with open(args.processfile, 'r') as f:
        processfile = json.load(f)
        processes = list(processfile["Processes"].keys())
        basedir = processfile["Basedir"]
        subbasedir = basedir.split("/")[-1]
        args.storage = os.path.join(args.storage, subbasedir)

    path_to_rootfile = os.path.join(args.outputpath, f"{datacard_settings['DC_name']}.root")
    rootfile = uproot.recreate(path_to_rootfile)

    shape_systematics = load_uncertainties(args.systematicsfile, allowflat=False)
    shape_systematics["nominal"] = Uncertainty("nominal", {})
    shape_systematics["stat_unc"] = Uncertainty("stat_unc", {})
    patch_scalevar_correlations(shape_systematics, processes)

    if args.UseEFT:
        eft_part, asimov_signal = eft_datacard_creation(rootfile, datacard_settings, args.eft_operator, shape_systematics, args)
        processes = [process for process in processes if process != "TTTT"]
        processes_write = [[process, i + 1] for i, process in enumerate(processes)]
    elif args.UseBSM:
        # find BSM process in processlist:
        for process in processes:
            if processfile["Processes"][process].get("isSignal", 0) > 0:
                bsm_process = {process: processfile["Processes"][process]}
                bsm_processname = process
                break
        bsm_part = bsm_datacard_creation(rootfile, datacard_settings, ["BSM_Lin", "BSM_Quad", "BSM_Cubic", "BSM_Quartic"], shape_systematics, bsm_process, args)
        processes = [process for process in processes if process != bsm_processname]
        processes_write = [[process, i + 1] for i, process in enumerate(processes)]
    else:
        processes_write = []
        sig_nb = -1
        bkg_nb = 1
        for process in processes:
            if processfile["Processes"][process].get("isSignal", 0) > 0:
                processes_write.append([process, sig_nb])
                sig_nb -= 1
            else:
                processes_write.append([process, bkg_nb])
                bkg_nb += 1
        # processes_write = [[process, i + 1] for i, process in enumerate(processes)]

    asimov_bkg = nominal_datacard_creation(rootfile, datacard_settings, channels, processes, shape_systematics, args)
    if args.UseBSM:
        processes_write.extend(bsm_part)

    if args.UseEFT:
        processes_write.extend(eft_part)
        for channelname, channel_DC_setting in datacard_settings["channelcontent"].items():
            asimov_data_path = f"{channel_DC_setting['prettyname']}/data_obs"
            var_name = channel_DC_setting["variable"]
            variables = VariableReader(args.variablefile, [var_name])
            asimov_data = asimov_bkg[channel_DC_setting['prettyname']] + asimov_signal[channel_DC_setting['prettyname']]
            convert_and_write_histogram(asimov_data, variables.get_properties(var_name), asimov_data_path, rootfile, statunc=np.sqrt(asimov_data))
    else:
        for channelname, channel_DC_setting in datacard_settings["channelcontent"].items():
            asimov_data_path = f"{channel_DC_setting['prettyname']}/data_obs"
            var_name = channel_DC_setting["variable"]
            variables = VariableReader(args.variablefile, [var_name])
            convert_and_write_histogram(asimov_bkg[channel_DC_setting['prettyname']], variables.get_properties(var_name), asimov_data_path, rootfile, statunc=np.sqrt(asimov_bkg[channel_DC_setting['prettyname']]))

    rootfile.close()

    # start txt writing
    path_to_txtfile = os.path.join(args.outputpath, f"{datacard_settings['DC_name']}.txt")
    dc_writer = DatacardWriter(path_to_txtfile)
    dc_writer.initialize_datacard(len(datacard_settings["channelcontent"]), f"{datacard_settings['DC_name']}.root")

    relevant_channels = {datacard_settings["channelcontent"][name]["prettyname"]: channel for name, channel  in channels.items() if name in list(datacard_settings["channelcontent"].keys())}
    dc_writer.add_channels(relevant_channels)
    dc_writer.add_processes(processes_write)

    systematics = load_uncertainties(args.systematicsfile)
    if args.UseEFT:
        processes.append("TTTT")
        for proc, nb in processes_write:
            if "quad" in proc:
                processes.append(proc)

    patch_scalevar_correlations(systematics, processes)
    if args.UseEFT:
        systematics["tttt_norm"] = Uncertainty("TTTTNorm", {"rate": "0.88/1.04", "processes": ["sm"]})

    for syst_name, syst_info in systematics.items():
        dc_writer.add_systematic(syst_info)

    dc_writer.write_card()

    exit()
