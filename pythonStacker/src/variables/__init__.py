import awkward as ak
import uproot
import numpy as np


def load_variable(tree: uproot.TTree, branchname: str, selection: str):
    arr = tree.arrays(["variable"], cut=selection, aliases={"variable": branchname})
    return ak.to_numpy(arr.variable)


def bdt_tanh_transfrom_nonflat(tree: uproot.TTree, branchname: str, selection: str):
    """
    Note: variable should be treatable as a regular-shaped numpy-like array (ideally)
    """
    variable = tree.arrays(["variable"], cut=selection, aliases={"variable": branchname}).variable
    # variable = ak.where(variable > 1., ak.full_like(variable, 0.9999), variable)
    old_shape = ak.num(variable, axis=-1)
    flattened = ak.to_numpy(ak.flatten(variable))
    flattened[flattened > 1.] = 0.9999
    transformed = np.arctanh((flattened + 1.) * 0.5)
    unflattened = ak.unflatten(transformed, old_shape)

    return unflattened


def bdt_tanh_transfrom_flat(tree: uproot.TTree, branchname: str, selection: str):
    """
    Note: variable should be treatable as a regular-shaped numpy-like array (ideally)
    """
    variable = tree.arrays(["variable"], cut=selection, aliases={"variable": branchname}).variable
    # variable = ak.where(variable > 1., ak.full_like(variable, 0.9999), variable)
    flattened = ak.to_numpy(variable)
    flattened[flattened > 1.] = 0.9999
    transformed = np.arctanh((flattened + 1.) * 0.5)

    return transformed


def sum_branch(tree: uproot.TTree, branchname: str, selection: str):
    variable = tree.arrays(["variable"], cut=selection, aliases={"variable": branchname}).variable
    return ak.to_numpy(ak.sum(variable, axis=1))


def count_btag(tree: uproot.TTree, branchname: str, selection: str, wp: int):
    variable = tree.arrays(["variable"], cut=selection, aliases={"variable": branchname}).variable
    return ak.to_numpy(ak.sum(variable >= wp, axis=1))


def count_btag_loose(tree: uproot.TTree, branchname: str, selection: str):
    return count_btag(tree, branchname, selection, 1)


def count_btag_med(tree: uproot.TTree, branchname: str, selection: str):
    return count_btag(tree, branchname, selection, 2)


def count_btag_tight(tree: uproot.TTree, branchname: str, selection: str):
    return count_btag(tree, branchname, selection, 3)


def sum_of_lepton_charges(tree: uproot.TTree, branchname: str, selection: str):
    variables = tree.arrays(["electronCharge", "muonCharge"], cut=selection)
    charges = ak.concatenate([variables.electronCharge, variables.muonCharge], axis=1)
    final_var = ak.sum(charges, axis=1)
    return ak.to_numpy(final_var)


def get_method_from_str(method: str):
    tmp = {
        "load": load_variable,
        "default": load_variable,
        "tanh_nonflat": bdt_tanh_transfrom_nonflat,
        "tanh_flat": bdt_tanh_transfrom_flat,
        "sum_branch": sum_branch,
        "btag_tight": count_btag_tight,
        "btag_med": count_btag_med,
        "btag_loose": count_btag_loose,
        "sum_of_lepton_charges": sum_of_lepton_charges
    }
    return tmp[method]
