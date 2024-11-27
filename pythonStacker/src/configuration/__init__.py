import json
import re

from src.configuration.Uncertainty import Uncertainty
import src.configuration.PDFVariation as PDFVariation
import src.configuration.ScaleVariation as ScaleVariation


def load_uncertainties(jsonfile, typefilter=None, namefilter=None, allowflat=True) -> dict:
    with open(jsonfile, 'r') as f:
        all_data = json.load(f)
    ret = dict()
    for key, val in all_data.items():
        syst_type = val.get("type", "flat")
        if typefilter and syst_type != typefilter:
            continue
        if namefilter and val.get("name", key) != namefilter:
            continue
        if not allowflat and syst_type == "flat":
            continue

        if key.lower() == "pdfuncertainty":
            pdfs = PDFVariation.generate_pdfvariations(val)
            for pdf in pdfs:
                ret[pdf.pretty_name] = pdf
            continue
        if key.lower() == "scalevar":
            scalevars = ScaleVariation.generate_scalevariations(val)
            for scalevar in scalevars:
                ret[scalevar.pretty_name] = scalevar
            continue

        ret[key] = Uncertainty(key, val)

    return ret


def load_channels(channelfile) -> dict:
    ret = dict()
    with open(channelfile, 'r') as f:
        data = json.load(f)
    for key, val in data.items():
        channel_tmp = Channel(val, data)
        if channel_tmp.isSubchannel:
            continue
        ret[key] = channel_tmp
    return ret


def load_channels_and_subchannels(channelfile) -> dict:
    channels = load_channels(channelfile)
    ret = dict(channels)
    for channelname, channelinfo in channels.items():
        subchannels = channelinfo.get_subchannels_dict()
        for subchannelname, subchannelinfo in subchannels.items():
            ret[channelname + subchannelname] = subchannelinfo

    return ret


# TODO: fix init ?
class Channel:
    def __init__(self, channelinfo, full_channelfile, ignore_proc=[]) -> None:
        self.channelinfo = channelinfo
        self._selection = channelinfo["selection"]
        self._isSubchannel = bool(channelinfo.get("isSubchannel", 0))

        # build process ignore list
        self.ign_processes = channelinfo.get("ignore_processes", [])
        if type(self.ign_processes) is not list:
            self.ign_processes = [self.ign_processes]
        self.ign_processes.extend(ignore_proc)
        self.reg_ign_processes = [re.compile(process) for process in self.ign_processes]


        # subchannel structure
        subchannels: list[str] = channelinfo.get("subchannels", [])
        self.subchannels: dict[str, Channel] = {subchannel: Channel(full_channelfile[subchannel], full_channelfile, self.ign_processes) for subchannel in subchannels}

        # load other config options

    def is_process_excluded(self, process: str):
        # return process in self.ign_processes
        return any([bool(ign_process.match(process)) for ign_process in self.reg_ign_processes])    


    def get_subchannels(self) -> list[str]:
        ret = []
        for name, subchannel in self.subchannels.items():
            ret.append(name)
            for subname in subchannel.get_subchannels():
                ret.append(name + "_" + subname)
        return ret
    
    def get_subchannels_dict(self) -> dict:
        ret = dict(self.subchannels)
        for name, subchannel in self.subchannels.items():
            for subname in subchannel.get_subchannels():
                ret[name + "_" + subname] = subchannel.subchannels[subname]
        return ret

    def produce_masks(self, tree):
        # use produce aliases to produce boolean masks for the subchannels
        aliases, keys = self.produce_aliases()
        masks = tree.arrays(keys, cut=self.selection, aliases=aliases)
        return masks, keys

    def produce_aliases(self) -> tuple[dict, list]:
        # for the chosen channel, load the info
        aliases: dict = {}
        for name, subchannel in self.subchannels.items():
            aliases[name] = subchannel.selection
            for subname in subchannel.get_subchannels():
                aliases[name + "_" + subname] = f"({subchannel.selection}) & ({subchannel.subchannels[subname].selection})"

        alias_names = list(aliases.keys())
        return aliases, alias_names

    @property
    def selection(self) -> str:
        return self._selection

    @selection.setter
    def selection(self, value):
        if not isinstance(value, str):
            raise ValueError("selection must be a str value")
        self._selection = value

    @property
    def isSubchannel(self) -> bool:
        return self._isSubchannel

    @isSubchannel.setter
    def isSubchannel(self, value):
        if not isinstance(value, bool):
            raise ValueError("isSubchannel must be a bool value")
        self._isSubchannel = value


# might replace with a class just keeping track of a few arrays.
# TODO: start from json file or at least dict
class Process:
    def __init__(self) -> None:
        self._isSignal = False
        self.name = ""

    @property
    def isSignal(self):
        return self._isSignal

    @isSignal.setter
    def isSignal(self, value):
        if not isinstance(value, bool):
            raise ValueError("isSignal must be a boolean value")
        self._isSignal = value
