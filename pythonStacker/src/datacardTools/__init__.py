from src.configuration.Uncertainty import Uncertainty
from copy import deepcopy

class DatacardWriter():
    def __init__(self, filename) -> None:
        self.file = filename
        self.outputstring = ""
        self.channels = []
        self.processes = []

    def commentline(self, length=200):
        self.outputstring += "-" * length + "\n"

    def initialize_datacard(self, number_of_regions, rootfile):
        self.outputstring = f"imax {number_of_regions}\n"
        self.outputstring += "jmax *\n"
        self.outputstring += "kmax *\n"
        self.commentline()
        self.outputstring += f"shapes * * {rootfile} $CHANNEL/$PROCESS $CHANNEL/$SYSTEMATIC/$PROCESS\n"
        self.commentline()

    def add_channels(self, channels: dict):
        self.channels = channels
        self.channel_names = list(channels.keys())
        self.outputstring += "{:>25s}".format("bin")
        for channel in self.channel_names:
            self.outputstring += f"\t{channel:>15s}"
        self.outputstring += "\n"
        self.outputstring += "{:>25s}".format("observation")
        observationline = f"\t{-1:>15d}" * len(channels)
        self.outputstring += observationline
        self.outputstring += "\n"
        self.commentline()

    def add_processes(self, processes: list):
        if len(self.channels) == 0:
            raise Exception("No channels defined! Run add_channels before add_processes!")

        channelline = ""
        processline = ""
        processnumber = ""
        rate = ""
        self.processes = processes
        # process should be an object itself containing a lot of information
        # same for channel -> not super lightweight but ok
        # need to refix this
        for channelname, channel in self.channels.items():
            for processname, number in self.processes:
                if channel.is_process_excluded(processname):
                    continue

                channelline += "{:>20s}\t".format(channelname)
                processline += "{:>20s}\t".format(processname)
                processnumber += "{:>20d}\t".format(number)
                rate += "{:>20d}\t".format(-1)

        self.outputstring += "{:>30s}\t".format("bin") + channelline + "\n"
        self.outputstring += "{:>30s}\t".format("process") + processline + "\n"
        self.outputstring += "{:>30s}\t".format("process") + processnumber + "\n"
        self.outputstring += "{:>30s}\t".format("rate") + rate + "\n"
        self.commentline()

    def add_MCstats(self):
        self.commentline()
        self.outputstring += "* autoMCStats 0 1 1\n"

    def add_RateParamNormalization(self, processname, range):
        # TODO implement
        pass

    def add_systematic(self, systematic: Uncertainty):
        if systematic.correlated_process:
            self.add_systematic_correlated(systematic)
        else:
            self.add_systematic_uncorrelated(systematic)

    def add_systematic_uncorrelated(self, systematic: Uncertainty):
        for process, _ in self.processes:
            if not systematic.is_process_relevant(process):
                continue
            processname = process
            if "sm" == processname:
                processname = "TTTT"
            # make a copy
            systematic_mod = deepcopy(systematic)

            # change name, change rrelevant processes
            systematic_mod.pretty_name = systematic.pretty_name + processname
            systematic_mod.technical_name = systematic.technical_name + processname
            systematic_mod.name = systematic.name + processname
            systematic_mod.set_processes(["^"+process+"$"])
            # then use nominal addition
            self.add_systematic_correlated(systematic_mod)

        # relevant = False
        # print(systematic.technical_name)
        # # get important processes:
        # outputlines = {}
        # for process, _ in self.processes:
        #     if process == "sm" and systematic.is_process_relevant("TTTT"):
        #         outputlines["TTTT"] = ""
# 
        #     if not systematic.is_process_relevant(process):
        #         continue
        #     outputlines[process] = ""
# 
        # for _, channel in self.channels.items():
        #     if not systematic.is_channel_relevant(channel):
        #         continue
        #     for process, number in self.processes:
        #         print(process)
        #         if channel.is_process_excluded(process):
        #             continue
        #         is_relevant_eft_variation = ("sm" in process or "quad" in process) and systematic.is_process_relevant("TTTT")
        #         print(process)
        #         filled = False
        #         if systematic.is_process_relevant(process):
        #             if isinstance(systematic.rate, str):
        #                 outputlines[process] += "\t{:>20s}".format(systematic.rate)
        #             else:
        #                 outputlines[process] += "\t{:>20.2f}".format(systematic.rate)
        #             relevant = True
        #             filled = True
        #         elif is_relevant_eft_variation:
        #             outputlines["TTTT"] += "\t{:>20.2f}".format(systematic.rate)
        #             relevant = True
        #             filled = True
        #         for proc_key in outputlines.keys():
        #             if (process == proc_key and filled) or (is_relevant_eft_variation and proc_key == "TTTT"):
        #                 continue
        #             # outputlines[process] += "\t{:>20s}".format("-")
        #             outputlines[proc_key] += "\t{:>20s}".format("-")
# 
        # if relevant:
        #     for proc, outputline in outputlines.items():
        #         systematic_name = systematic.technical_name + proc
        #         self.outputstring += "{:<23s} ".format(systematic_name)
        #         if systematic.isFlat:
        #             self.outputstring += "{:>6s}".format("lnN")
        #         else:
        #             self.outputstring += "{:>6s}".format("shape")
# 
        #         self.outputstring += outputline + "\n"
# 

    def add_systematic_correlated(self, systematic: Uncertainty):
        # similar to add process stuff
        systematic_line = ""
        relevant = False
        for _, channel in self.channels.items():
            if not systematic.is_channel_relevant(channel):
                continue
            for process, number in self.processes:
                if channel.is_process_excluded(process):
                    continue
                if systematic.is_process_relevant(process):
                    if isinstance(systematic.rate, str):
                        systematic_line += "\t{:>20s}".format(systematic.rate)
                    else:
                        systematic_line += "\t{:>20.2f}".format(systematic.rate)
                    relevant = True
                else:
                    systematic_line += "\t{:>20s}".format("-")

        if relevant:
            self.outputstring += "{:<23s} ".format(systematic.technical_name)
            if systematic.isFlat:
                self.outputstring += "{:>6s}".format("lnN")
            else:
                self.outputstring += "{:>6s}".format("shape")

            self.outputstring += systematic_line + "\n"

    def write_card(self):
        # write autoMCStats line
        self.add_MCstats()
        with open(self.file, "w") as f:
            f.write(self.outputstring)


if __name__ == "__main__":
    def test_DatacardWrite():
        # Create a DatacardWrite object
        datacard = DatacardWriter('/home/njovdnbo/Documents/Stacker_v2/pythonStacker/output/test/test.txt')

        # Initialize the datacard with 3 regions
        datacard.initialize_datacard(3, "blub.root")

        # Add a comment line
        datacard.add_channels(["pretty"], ["BDT_EM_region"])
        datacard.add_processes(["blub"])

        # Write the datacard to the file
        datacard.write_card()

    # Run the test
    test_DatacardWrite()
