import numpy as np
from scipy.stats import chi2
import awkward as ak

from src.histogramTools import HistogramManager
from src.variables.variableReader import Variable


def data_stat_uncertainties(n_events, alpha=1 - 0.6827):
    lower_unc = np.where(n_events == 0, 0, chi2.ppf(q=alpha / 2, df=2 * n_events) / 2)
    upper_unc = chi2.ppf(q=1 - (alpha / 2), df=(2 * n_events) + 2) / 2
    return (lower_unc, upper_unc)


class DataManager:
    def __init__(self, data_path, variables, channel, years, eras_split=False):
        if type(years) is not list:
            years = [years]
        # object to load the actual data
        self.histograms = dict()

        systematics_dummy = ["nominal"]
        for year in years:
            self.histograms[year] = HistogramManager(data_path, "Data", variables, systematics_dummy, year, channel=channel)
            self.histograms[year].load_histograms()
        self.eras_split = eras_split
        
        self.histograms_per_era = dict()
        if eras_split:
            runs_per_era = {
                "2016PreVFP": ["B", "C", "D", "E", "F"],
                "2016PostVFP": ["F", "G", "H"],
                "2017": ["B", "C", "D", "E", "F"],
                "2018": ["A", "B", "C", "D"],
            }
            for year in years:
                runs = runs_per_era[year]
                for run in runs:
                    self.histograms_per_era[year+run] = HistogramManager(data_path, "Data"+year+run, variables, systematics_dummy, year, channel=channel)
                    self.histograms_per_era[year+run].load_histograms()

    def get_histogram_and_uncertainties(self, years, variable: Variable):
        if type(years) is not list:
            years = [years]
        content = np.zeros(variable.nbins)
        for year in years:
            content += ak.to_numpy(self.histograms[year][variable.name]["nominal"])

        uncertainty = data_stat_uncertainties(content)
        return content, uncertainty

