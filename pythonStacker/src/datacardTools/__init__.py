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
        self.outputstring += "* autoMCStats 0 0 1\n"
        self.commentline()
        self.outputstring += "PDF group = pdf_1 pdf_10 pdf_0 pdf_11 pdf_12 pdf_13 pdf_14 pdf_15 pdf_16 pdf_17 pdf_18 pdf_19 pdf_2 pdf_20 pdf_21 pdf_22 pdf_23 pdf_24 pdf_25 pdf_26 pdf_27 pdf_28 pdf_29 pdf_3 pdf_30 pdf_31 pdf_32 pdf_33 pdf_34 pdf_35 pdf_36 pdf_37 pdf_38 pdf_39 pdf_4 pdf_40 pdf_41 pdf_42 pdf_43 pdf_44 pdf_45 pdf_46 pdf_47 pdf_48 pdf_49 pdf_5 pdf_50 pdf_51 pdf_52 pdf_53 pdf_54 pdf_55 pdf_56 pdf_57 pdf_58 pdf_59 pdf_6 pdf_60 pdf_61 pdf_62 pdf_63 pdf_64 pdf_65 pdf_66 pdf_67 pdf_68 pdf_69 pdf_7 pdf_70 pdf_71 pdf_72 pdf_73 pdf_74 pdf_75 pdf_76 pdf_77 pdf_78 pdf_79 pdf_8 pdf_80 pdf_81 pdf_82 pdf_83 pdf_84 pdf_85 pdf_86 pdf_87 pdf_88 pdf_89 pdf_9 pdf_90 pdf_91 pdf_92 pdf_93 pdf_94 pdf_95 pdf_96 pdf_97 pdf_98 pdf_99\n"
        self.outputstring += "Lumi group = lumi_13TeV lumi_13TeV_1718 lumi_2016 lumi_2017 lumi_2018\nBtag group = CMS_btag_cferr1 CMS_btag_cferr2 CMS_btag_hf CMS_btag_hfstats1_2016 CMS_btag_hfstats1_2017 CMS_btag_hfstats1_2018 CMS_btag_hfstats2_2016 CMS_btag_hfstats2_2017 CMS_btag_hfstats2_2018 CMS_btag_lf CMS_btag_lfstats1_2016 CMS_btag_lfstats1_2017 CMS_btag_lfstats1_2018 CMS_btag_lfstats2_2016 CMS_btag_lfstats2_2017 CMS_btag_lfstats2_2018\nleptonID group = CMS_eff_e_id_2016 CMS_eff_e_id_2017 CMS_eff_e_id_2018 CMS_eff_e_id_correlated CMS_eff_e_reco CMS_eff_m_id_2016 CMS_eff_m_id_2017 CMS_eff_m_id_2018 CMS_eff_m_id_correlated\n"
        self.outputstring += "PS group = ps_isr_Othert ps_isr_Xg ps_isr_ttH ps_isr_ttW ps_isr_ttZ ps_isr_BSMsignal ps_fsr\nJEC group = CMS_scale_j CMS_res_j_1p93_2016 CMS_res_j_1p93_2017 CMS_res_j_1p93_2018 CMS_res_j_2p5_2016 CMS_res_j_2p5_2017 CMS_res_j_2p5_2018\nWZNjets group = CMS_TOP24008_WZ_j_2 CMS_TOP24008_WZ_j_3 CMS_TOP24008_WZ_j_4 CMS_TOP24008_WZ_j_5 CMS_TOP24008_WZ_j_6\n"
        self.outputstring += "Norm group = CMS_fake_m_2018 CMS_fake_m_2017 CMS_fake_m_2016 CMS_fake_m CMS_fake_e_2018 CMS_fake_e_2017 CMS_fake_e_2016 CMS_fake_e Norm_Xg Norm_ttH Norm_ttW Norm_ttZ Norm_Othert Norm_ChargeMisID Norm_BSMsignal\nrest group = CMS_HEM_2018 CMS_l1_ecal_prefiring CMS_scale_met CMS_pileup_13TeV\nscale group = QCDscale_fac_Othert QCDscale_fac_Xg QCDscale_fac_ttH QCDscale_fac_ttW QCDscale_fac_ttZ QCDscale_fac_BSMsignal QCDscale_ren_Othert QCDscale_ren_Xg QCDscale_ren_ttH QCDscale_ren_ttW QCDscale_ren_ttZ QCDscale_ren_BSMsignal\n"

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
            #systematic_mod.set_processes([process])
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
