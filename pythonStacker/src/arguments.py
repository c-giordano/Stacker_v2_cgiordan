import argparse
import os
import sys

def add_settingfiles(parser: argparse.ArgumentParser):
    # Add arguments here
    parser.add_argument('-vf', '--variablefile', dest="variablefile",
                        type=str, help='JSON file with variables.')
    parser.add_argument('-sf', '--systematicsfile', dest="systematicsfile",
                        type=str, help='JSON file with systematics.')
    parser.add_argument('-pf', '--processfile', dest='processfile',
                        type=str, help='JSON file with process definitions.')
    parser.add_argument('-cf', '--channelfile', dest='channelfile',
                        type=str, help='JSON file with channel definitions.')
    parser.add_argument('-y', '--years', dest='years', default=["2016", "2017", "2018"], nargs='+',
                        help='Specific years.')


def select_specifics(parser: argparse.ArgumentParser):
    # Add arguments here
    parser.add_argument('-v', '--variable', dest="variable", default=None,
                        type=str, help='Specific variable.')
    # Syst can be "weight", "shape", None, or a specific name
    parser.add_argument('-s', '--systematic', dest="systematic", default=None,
                        type=str, help='Specific systematic.')
    parser.add_argument('-p', '--process', dest='process', default=None,
                        type=str, help='Specific process.')
    parser.add_argument('-c', '--channel', dest='channel', default=None,
                        type=str, help='Specific channel.')


def add_toggles(parser: argparse.ArgumentParser):
    parser.add_argument('--EFT', '--eft', dest="UseEFT", default=False, action="store_true",
                        help="toggle to include EFT variations")
    parser.add_argument('--BSM', '--bsm', dest="UseBSM", default=False, action="store_true",
                        help="toggle to include BSM reweighted variations.")


def add_EFT_choice(parser: argparse.ArgumentParser):
    parser.add_argument('--EFTop', dest="eft_operator",nargs='+', default=["cQQ1"],
                        help="string to choose which WC to consider.")


def add_tmp_storage(parser: argparse.ArgumentParser):
    parser.add_argument("--storage", dest="storage", type=str,
                        default="Intermediate", help="Path at which the \
                        histograms are stored")


def add_plot_output(parser: argparse.ArgumentParser):
    default_output = f"/user/{os.getenv('USER')}/public_html/Interpretations/Plots/"
    if os.getenv("CMSSW_VERSION") is None:
        default_output = "output/"
    parser.add_argument("-o", "--output", dest="outputfolder", action="store", required=False,
                        default=default_output, help="outputfolder to use for plots")
