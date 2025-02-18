#!/usr/bin/env python3
import ROOT
import os
import argparse

parser = argparse.ArgumentParser(description="Script to process ROOT files and plot histograms.")
parser.add_argument('--data', action='store_true', default=False, help='do you want to plot data?')
parser.add_argument('--version', action='store', default='1', help='which version is this?')
parser.add_argument('--fileDir', action='store', default='./output/bsmlimits/TopPhilicScalarOctet_data/600/', help="where the DC at")

args = parser.parse_args()



# SR2L is split into subdirectories for each category in the new DC
new_sig_dirs = ["SR2L_ee_Sig", "SR2L_em_Sig", "SR2L_mm_Sig"]
new_ttw_dirs = ["SR2L_ee_ttw", "SR2L_em_ttw", "SR2L_mm_ttw"]
# SR2L not split in the old DC
old_sig_dir = "SR2L_Sig"
old_ttw_dir = "SR2L_ttw"

regions_to_compare = ["SR3L_Sig", "SR3L_ttw", "SR3L_NP", "CR2LNP", "CR2LTTW", "CR3LZ", "CR4LZ", "CR3LNP"]
skip_histograms = {}

years = ["2018", "2017", "2016"]

def accumulate_integrals(file, dir_list, category=""):
    """
    Loop over a list of directories in the given file, accumulate integrals for each histogram,
    and return a dictionary of the accumulated integrals.
    """
    integrals = {}  # {histogram_name: accumulated_integral}
    for dname in dir_list:
        print(f"\nProcessing directory: {dname} for category: {category}")
        dir_obj = file.Get(dname)
        if not dir_obj:
            print(f"  [!] Directory {dname} not found in the file!")
            continue

        # Loop over keys directly in the directory (ignoring subdirectories)
        for key in dir_obj.GetListOfKeys():
            hname = key.GetName()
            if hname in skip_histograms:
                print(f"  Skipping histogram: {hname} in {dname}")
                continue
            obj = key.ReadObj()
            if obj.InheritsFrom("TH1") and not obj.InheritsFrom("TDirectory"):
                h_integral = obj.Integral()
                integrals[hname] = integrals.get(hname, 0) + h_integral
    return integrals

def get_integrals_from_directory(file, directory, category=""):
    """
    Get the integrals of histograms from a single directory in the given file, returning a dictionary of integrals.
    """
    integrals = {}
    print(f"\nProcessing directory: {directory} for category: {category}")
    dir_obj = file.Get(directory)
    if not dir_obj:
        print(f"  [!] Directory {directory} not found in the file!")
        return integrals

    for key in dir_obj.GetListOfKeys():
        hname = key.GetName()
        if hname in skip_histograms:
            print(f"  Skipping histogram: {hname} in {directory}")
            continue
        obj = key.ReadObj()
        if obj.InheritsFrom("TH1") and not obj.InheritsFrom("TDirectory"):
            h_integral = obj.Integral()
            print(f"  Found histogram: {hname} in {directory} with integral: {h_integral}")
            integrals[hname] = h_integral
    return integrals

def print_single_region_table(region, ints):
    """
    Print a table for a single region from new files only.
    """
    print(f"\nAggregated integrals for {region} (New Files):")
    header = "{:<20} {:<20}".format("Histogram", "New Integral")
    print(header)
    print("-" * len(header))
    for hname, val in sorted(ints.items()):
        print("{:<20} {:<20.3f}".format(hname, val))

def process_new_SR4L_regions(file):
    """
    Process the SR4L regions in the new file.
    Returns a tuple of dictionaries: (SR4L_Sig integrals, SR4L_ttw integrals)
    """
    sig_ints = get_integrals_from_directory(file, "SR4L_Sig", category="SR4L_Sig")
    ttw_ints = get_integrals_from_directory(file, "SR4L_ttw", category="SR4L_ttw")
    return sig_ints, ttw_ints

def print_comparison_table(category, new_ints, old_ints):
    """
    Print a table comparing the integrals from the new file with those from the old file.
    (No Match column is printed.)
    """
    print(f"\nComparison Table for {category}:")
    header = "{:<20} {:<20} {:<20}".format("Histogram", "New Integral", "Old Integral")
    print(header)
    print("-" * len(header))
    all_hist_names = sorted(set(list(new_ints.keys()) + list(old_ints.keys())))
    for hname in all_hist_names:
        new_val = new_ints.get(hname, 0)
        old_val = old_ints.get(hname, 0)
        print("{:<20} {:<20.3f} {:<20.3f}".format(hname, new_val, old_val))

def add_dicts(target, source):
    """
    Helper function to add the values from source dictionary into target dictionary.
    """
    for key, val in source.items():
        target[key] = target.get(key, 0) + val


agg_new_sig = {}
agg_new_ttw = {}
agg_old_sig = {}
agg_old_ttw = {}

# for region comparisons, initialize a dictionary for each region
agg_new_regions = {region: {} for region in regions_to_compare}
agg_old_regions = {region: {} for region in regions_to_compare}

agg_new_SR4L_sig = {}
agg_new_SR4L_ttw = {}


# --- Main Processing Loop Over Years ---
for year in years:
    print(f"\n================== Processing Year: {year} ==================")
    
    # Change your paths
    filename_new = f"./output/bsmlimits/TopPhilicScalarOctet_data/600/DC_{year}.root"
    filename_old = f"../../newPlots/plots/pythonStacker/output/bsmlimits/TopPhilicScalarOctet_v2/600/DC_{year}.root"
    
    fnew = ROOT.TFile.Open(filename_new, "READ")
    if not fnew or fnew.IsZombie():
        print("Error opening new file:", filename_new)
        continue

    fold = ROOT.TFile.Open(filename_old, "READ")
    if not fold or fold.IsZombie():
        print("Error opening old file:", filename_old)
        fnew.Close()
        continue


    print("\n=== New File: Accumulating across subdirectories ===")
    new_sig_integrals = accumulate_integrals(fnew, new_sig_dirs, category="Sig")
    new_ttw_integrals = accumulate_integrals(fnew, new_ttw_dirs, category="ttw")

    print("\n=== Old File: Single directories ===")
    old_sig_integrals = get_integrals_from_directory(fold, old_sig_dir, category="Sig")
    old_ttw_integrals = get_integrals_from_directory(fold, old_ttw_dir, category="ttw")

    add_dicts(agg_new_sig, new_sig_integrals)
    add_dicts(agg_new_ttw, new_ttw_integrals)
    add_dicts(agg_old_sig, old_sig_integrals)
    add_dicts(agg_old_ttw, old_ttw_integrals)

    for region in regions_to_compare:
        new_reg = get_integrals_from_directory(fnew, region, category=region)
        old_reg = get_integrals_from_directory(fold, region, category=region)
        add_dicts(agg_new_regions[region], new_reg)
        add_dicts(agg_old_regions[region], old_reg)
    
    new_SR4L_sig, new_SR4L_ttw = process_new_SR4L_regions(fnew)
    add_dicts(agg_new_SR4L_sig, new_SR4L_sig)
    add_dicts(agg_new_SR4L_ttw, new_SR4L_ttw)

    fnew.Close()
    fold.Close()

# --- Final Aggregated Comparison Tables ---

print("\n================== FINAL AGGREGATED COMPARISONS (All Years) ==================")

print_comparison_table("SR2L Sig (Aggregated)", agg_new_sig, agg_old_sig)
print_comparison_table("SR2L ttw (Aggregated)", agg_new_ttw, agg_old_ttw)

for region in regions_to_compare:
    print_comparison_table(f"Region {region} (Aggregated)", agg_new_regions[region], agg_old_regions[region])

print_single_region_table("SR4L_Sig", agg_new_SR4L_sig)
print_single_region_table("SR4L_ttw", agg_new_SR4L_ttw)