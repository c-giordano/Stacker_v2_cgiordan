import re


class Uncertainty:
    def __init__(self, name, dict_entry):
        self.name = name

        # extract dictionary and add defaults for every property
        self.pretty_name = dict_entry.get("pretty_name", name)
        self.technical_name = dict_entry.get("technical_name", name)
        self.fileglob = dict_entry.get("fileglob", name)

        self.channels = dict_entry.get("channels", ["all"])
        if type(self.channels) is not list:
            self.channels = [self.channels]

        self.processes = dict_entry.get("processes", ["all"])
        if type(self.processes) is not list:
            self.processes = [self.processes]
        self.reg_processes = [re.compile(process) for process in self.processes]

        self._correlated_process = bool(dict_entry.get("corr_proc", True))

        self.correlated_years = True
        # Possible types: flat, weight, shape
        self.type = dict_entry.get("type", "flat")
        if self.type == "flat":
            self._isFlat = True
        else:
            self._isFlat = False

        self.rate = dict_entry.get("rate", 1.)

        self.weight_alias_up = dict_entry.get("aliasUp", "nominalWeight")
        if self.weight_alias_up == "nominalWeight":
            self.weight_key_up = "nominal"
        else:
            self.weight_key_up = self.name + "Up"

        self.weight_alias_down = dict_entry.get("aliasDown", "nominalWeight")
        if self.weight_alias_down == "nominalWeight":
            self.weight_key_down = "nominal"
        else:
            self.weight_key_down = self.name + "Down"

        if (dict_entry.get("isEFT", 0) > 0) or (dict_entry.get("isBSM", 0) > 0):
            self.weight_key_up = self.name
       
        self.exact = dict_entry.get("exact",False)

    def set_processes(self, processes):
        self.processes = processes
        self.reg_processes = [re.compile(process) for process in self.processes]

    @property
    def isFlat(self):
        return self._isFlat

    @isFlat.setter
    def isFlat(self, value):
        if isinstance(value, bool):
            self._isFlat = value
        else:
            raise ValueError("isFlat must be a boolean")

    @property
    def correlated_process(self):
        return self._correlated_process

    @correlated_process.setter
    def correlated_process(self, value):
        if isinstance(value, bool):
            self._correlated_process = value
        else:
            raise ValueError("correlated_process must be a boolean")

    def is_channel_relevant(self, channel):
        if self.channels[0] == "all":
            return True
        return channel in self.channels

    def is_process_relevant(self, process):
        if self.processes[0] == "all":
            return True
        elif self.processes[0].lower() in ["mconly", "mcall"]:
            return not ("nonPrompt" in process or "ChargeMisID" in process)
        return any([bool(relevant_process.match(process)) for relevant_process in self.reg_processes])    
#        if self.exact :
#            return any([relevant_process == process for relevant_process in self.processes])
#        else:
#            return any([relevant_process in process for relevant_process in self.processes])
#        #this does not work for me. Change for your case if needed
#        #return any([bool(relevant_process.match(process)) for relevant_process in self.reg_processes])    

    def get_weight_aliases(self):
        return (self.weight_alias_up, self.weight_alias_down)

    def get_weight_keys(self):
        return (self.weight_key_up, self.weight_key_down)
