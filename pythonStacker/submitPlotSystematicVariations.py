import argparse
import json

import src.jobSubmission.condorTools as ct
import src.arguments as arguments
from src.variables.variableReader import VariableReader, Variable


def parse_arguments():
    parser = argparse.ArgumentParser(description='Script to submit plotting of histograms')
    arguments.add_settingfiles(parser)
    arguments.select_specifics(parser)
    arguments.add_toggles(parser)
    parser.add_argument("--wc", action="store", default="ctt")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()

    basecommand = "python3 plotSystVariations.py"
    basecommand += f" --variablefile {args.variablefile}"
    basecommand += f" --processfile {args.processfile}"
    basecommand += f" --systematicsfile {args.systematicsfile}"
    basecommand += f" --channelfile {args.channelfile}"

    with open(args.channelfile, 'r') as f:
        channels = json.load(f)
        channellist = list(channels.keys())
    with open(args.processfile, 'r') as f:
        processlist = list(json.load(f)["Processes"].keys())
    # Load the relevant variables:
    variables = VariableReader(args.variablefile, args.variable)

    cmds = []
    for year in args.years:
        for process in processlist:
            if args.process is not None and process != args.process:
                continue
            for channel in channels:
                if args.channel is not None and channel != args.channel:
                    continue
                if channels[channel].get("isSubchannel", 0) > 0:
                    continue
                for variable, var_obj in variables.get_variable_objects().items():
                    if not var_obj.is_channel_relevant(channel):
                        continue
                    cmd = basecommand + f" -y {year}"
                    cmd += f" --process {process}"
                    cmd += f" -c {channel}"
                    cmd += f" -v {variable}"
                    cmds.append([cmd])


    ct.submitCommandsetsAsCondorCluster("plotsystematics", cmds, scriptfolder="Scripts/condor/")
