# import numpy as np
import awkward as ak
import numpy as np
import os
from src.variables.variableReader import VariableReader


class HistogramManager:
    def __init__(self, storagepath: str, process: str, variables: VariableReader, systematics: list, year="2018", channel=""):
        self.store: str = storagepath
        self.process: str = process
        self.variables: list = []  # variables.get_variables()

        self.systematics = systematics

        filename = process.replace(" ", "_") + "_" + year + "_"  # + sys + ".parquet"
        self.cache_folder = dict()
        self.base_name = dict()
        self.histograms = dict()
        for var, prop in variables.variable_objects.items():
            if not prop.is_channel_relevant(channel):
                continue
            self.variables.append(var)
            self.cache_folder[var] = os.path.join(storagepath, var)
            if not os.path.exists(self.cache_folder[var]):
                try:
                    os.makedirs(self.cache_folder[var])
                except FileExistsError:
                    print(f"Folder {self.cache_folder[var]} exists. Makedirs failed.")

            self.base_name[var] = os.path.join(self.cache_folder[var], filename)

            if self.systematics is None:
                continue

            self.histograms[var] = dict()
            for sysname in self.systematics:
                if sysname == "nominal" or sysname == "stat_unc":
                    self.histograms[var][sysname] = np.zeros(prop.nbins)
                else:
                    self.histograms[var][sysname] = dict()
                    self.histograms[var][sysname]["Up"] = np.zeros(prop.nbins)
                    self.histograms[var][sysname]["Down"] = np.zeros(prop.nbins)

    def __getitem__(self, key):
        return self.histograms[key]

    def __setitem__(self, key, value):
        self.histograms[key] = value

    def get_name(self, var: str, sys: str = "nominal"):
        # print(self.base_name[var] + sys + ".parquet")
        return self.base_name[var] + sys + ".parquet"

    def load_histograms(self):
        self.histograms = self.load_all_histograms()

    def load_all_histograms(self):
        ret = dict()
        for var in self.variables:
            ret[var] = self.load_histograms_for_variable(var)
        return ret

    def load_histograms_for_variable(self, var):
        ret = dict()
        pdfs = None

        for sys in self.systematics:
            if "PDF" in sys and pdfs is None:
                if os.path.exists(self.get_name(var, "PDF")):
                    pdfs = self.load_histogram(var, "PDF")
                    ret[sys] = pdfs[sys]
                continue
            elif "PDF" in sys:
                ret[sys] = pdfs[sys]
                continue
            ret[sys] = self.load_histogram(var, sys)
        return ret

    def load_histogram(self, var: str, sys: str = "nominal"):
        if os.path.exists(self.get_name(var, sys)):
            content = ak.from_parquet(self.get_name(var, sys))
        else:
            content = dict()
            if ("sm" in self.base_name[var] or "quad" in self.base_name[var]):  
                print(f"Note: Variable {var} for systematic {sys} does not exist with base {self.base_name[var]}.")
            content["Up"] = np.zeros(1)
            content["Down"] = np.zeros(1)

        # if sys == "nominal" or sys == "stat_unc":
        #     content = ak.from_parquet(self.get_name(var, sys))
        # else:
        #     content = ak.from_parquet(self.get_name(var, sys))
        return content

    # reconsider saving and loading routine. Can reprocess to use ak.Record to combine either many variables or many systematics
    # Better to have 1 file per variable? But then confusing with shape and weight systematics. 1 file per systematic makes sense
    # until we have to reproduce everything.
    # Can add something to still make a single file if a specific variable is desired.
    def save_histograms(self):
        self.save_all_histograms(self.histograms)

    def save_all_histograms(self, content):
        for var in self.variables:
            self.save_histograms_for_variable(content[var], var)

    def save_histograms_for_variable(self, content_variations, var):
        pdf_uncertainties = dict()
        for sys in self.systematics:
            if "PDF_" in sys:
                pdf_uncertainties[sys] = ak.Record(content_variations[sys])
                continue
            self.save_histogram(content_variations[sys], var, sys)
        if len(pdf_uncertainties) != 0:
            self.save_histogram(pdf_uncertainties, var, "PDF")

    def save_histogram(self, content, var: str, sys: str = "nominal"):
        if sys == "nominal" or sys == "stat_unc":
            ak.to_parquet(content, self.get_name(var, sys))
        elif type(content) is dict:
            record = ak.Record(content)
            ak.to_parquet(record, self.get_name(var, sys))
