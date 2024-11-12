import argparse
import json

import src.jobSubmission.condorTools as ct
import src.arguments as arguments


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

    basecommand = "python3 plotHistogramsRedo.py"
    basecommand += f" --variablefile {args.variablefile}"
    basecommand += f" --processfile {args.processfile}"
    basecommand += f" --systematicsfile {args.systematicsfile}"
    basecommand += f" --channelfile {args.channelfile}"
    basecommand += " --no_unc"

    with open(args.channelfile, 'r') as f:
        channels = json.load(f)
        channellist = list(channels.keys())

    if args.UseEFT:
        basecommand += " --EFT"
        basecommand += " --EFT_ratio --EFT_fullbkg"
        basecommand += f" --wc {args.wc}"
    cmds = []
    for year in args.years:
        for channel in channels:
            if args.channel is not None and channel != args.channel:
                continue
            if channels[channel].get("isSubchannel", 0) > 0:
                continue
            cmd = basecommand + f" -y {year}"
            cmd += f" -c {channel}"
            if args.UseData:
                cmd += " --data"
            cmds.append([cmd])

     #if "2016PreVFP" in args.years and "2016PostVFP" in args.years:
    #    cmd = basecommand + " -y 2016PreVFP 2016PostVFP"
    #    for channel in channels:
    #        if args.channel is not None and channel != args.channel:
    #            continue
    #        if channels[channel].get("isSubchannel", 0) > 0:
    #            continue
    #        cmd_tmp = cmd + f" -c {channel}"
    #        cmds.append([cmd_tmp])

    if len(args.years) >= 2:
        cmd = basecommand + " -y 2016 2017 2018"
        for channel in channels:
            if args.channel is not None and channel != args.channel:
                continue
            if channels[channel].get("isSubchannel", 0) > 0:
                continue
            cmd_tmp = cmd + f" -c {channel}"
            cmds.append([cmd_tmp])

    ct.submitCommandsetsAsCondorCluster("plothistograms", cmds, scriptfolder="Scripts/condor/")
