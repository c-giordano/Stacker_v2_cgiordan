from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List
import os

import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator
from pathlib import Path  

# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------
hep.style.use(hep.style.CMS)
plt.rcParams.update({
    "axes.grid"      : True,
    "grid.linestyle" : ":",
    "grid.alpha"     : 0.4,
    "legend.frameon" : False,
})

# ---------------------------------------------------------------------------
# Files & process labels
# ---------------------------------------------------------------------------
INPUT_FILES: Dict[str, Dict[str, str]] = {
    # "PseudoScalarOctet"  : {"directory": "/pnfs/iihe/cms/store/user/cgiordan/jsons/results_PO.json", "tex": r"$P_8$"},
    "PseudoScalarOctet"  : {"directory": "./results_PO.json", "tex": r"$P_8$"},   
    # "PseudoScalarSinglet": {"directory": "/pnfs/iihe/cms/store/user/cgiordan/jsons/results_PS.json", "tex": r"$P_1$"},
    "PseudoScalarSinglet": {"directory": "./results_PS.json", "tex": r"$P_1$"},
    # "ScalarOctet"        : {"directory": "/pnfs/iihe/cms/store/user/cgiordan/jsons/results_SO.json", "tex": r"$S_8$"},
    "ScalarOctet"        : {"directory": "./results_SO.json", "tex": r"$S_8$"},    
    # "ScalarSinglet"      : {"directory": "/pnfs/iihe/cms/store/user/cgiordan/jsons/results_SS.json", "tex": r"$S_1$"},
    "ScalarSinglet"      : {"directory": "./results_SS.json", "tex": r"$S_1$"},
    "VectorOctet"        : {"directory": "/pnfs/iihe/cms/store/user/cgiordan/jsons/results_VO.json", "tex": r"$V_8$"},
    "VectorSinglet"      : {"directory": "/pnfs/iihe/cms/store/user/cgiordan/jsons/results_VS.json", "tex": r"$V_1$"},  
    
}

PROCESSES: List[str] = ["TTTT", "TTTW", "TTTJ"]
MASS_POINTS: List[str] = ["0p4", "0p6", "0p8", "1p0", "1p2", "1p4", "1p6"]

# Dynamic y‑axis label
_dict_process_label = {
    "TTTT": r"$\sigma(pp \to t\bar{t}t\bar{t}) - \sigma_{SM}\;[\mathrm{pb}]$",
    "TTTJ": r"$\sigma(pp \to t\bar{t}tq) - \sigma_{SM}\;[\mathrm{pb}]$",
    "TTTW": r"$\sigma(pp \to t\bar{t}tW) - \sigma_{SM}\;[\mathrm{pb}]$",
}

# Visual palette
# COLORS = [
#     "#1f77b4",  # blue
#     "#ff7f0e",  # orange
#     "#2ca02c",  # green
#     "#d62728",  # red
#     "#9467bd",  # purple
#     "#8c564b",  # brown
# ]

COLORS = [
    "#0072B2",  # blue
    "#D55E00",  # burnt orange
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
]
MARKER = "o"

# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _json_path(info: Dict[str, str]) -> Path:
    p = Path(info["directory"])
    if not p.exists():
        raise FileNotFoundError(f"Input JSON not found: {p}")
    return p


def _mass_to_float(label: str) -> float:
    """Convert '0p4' → 0.4 (TeV)."""
    return float(label.replace("p", "."))

# ---------------------------------------------------------------------------
# Core plotting routine
# ---------------------------------------------------------------------------

def make_plots(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    # Cache JSON content ----------------------------------------------------
    cache: Dict[str, Dict[str, dict]] = {
        model: json.load(_json_path(info).open())
        for model, info in INPUT_FILES.items()
    }
    # y_ref_fb = 13.37              # fb
    # y_ref_pb = y_ref_fb * 1e-3
    for proc in PROCESSES:
        fig, ax = plt.subplots(figsize=(7, 5.5))

        xtick_vals = sorted({_mass_to_float(mp) for mp in MASS_POINTS})
        xtick_lbls = [f"{x:.1f}" for x in xtick_vals]

        ymax, ymin = 0.0, np.inf
        legend_handles: List[Line2D] = []

        # if proc == "TTTT":
        #     ax.axhline(y_ref_pb,
        #                color="k", linestyle="--", linewidth=1.4,
        #                label=f"{y_ref_fb:g} fb")

        # Loop over simplified models --------------------------------------
        for (model, info), colour in zip(INPUT_FILES.items(), COLORS):
            data = cache[model]
            masses, xsecs, errs = [], [], []

            for mp in MASS_POINTS:
                key = next((k for k in data if k.startswith(proc) and k.endswith(mp)), None)
                if key is None:
                    continue
                entry = data[key].get("c_0p2")
                if not entry:
                    continue

                xs_raw  = entry.get("xsec")
                err_raw = entry.get("error")
                if xs_raw is None:  # skip pathological None values
                    continue

                masses.append(_mass_to_float(mp))
                xsecs.append(float(xs_raw))
                errs.append(float(err_raw) if err_raw is not None else 0.0)

            if not masses:
                continue

            # Sort by mass -------------------------------------------------
            order  = np.argsort(masses)
            masses = np.array(masses)[order]
            xsecs  = np.array(xsecs)[order]
            errs   = np.array(errs)[order]

            ymax = max(ymax, xsecs.max())
            ymin = min(ymin, xsecs.min())

            # Plot --------------------------------------------------------
            ax.errorbar(
                masses,
                xsecs,
                yerr=errs,
                color=colour,
                marker=MARKER,
                ms=6,
                linestyle="-",
                lw=1.1,
            )

            legend_handles.append(
                Line2D(
                    [], [],
                    color=colour,
                    linestyle="-",   # linea continua
                    linewidth=1.5,   # spessore (regolalo come preferisci)
                    marker="None",   # niente marker
                    label=info["tex"],
                )
            )

        # Cosmetics ---------------------------------------------------------
        ax.set_xlabel(r"$M_{X}$ [TeV]", fontsize="xx-small")
        ax.set_ylabel(_dict_process_label[proc], fontsize="xx-small")
        ax.set_xlim(0.35, 1.65)
        ax.set_xticks(xtick_vals)
        # ax.set_xticklabels(xtick_lbls, fontsize="xx-small")
        ax.tick_params(axis="y", which="both", labelsize="x-small")
        ax.tick_params(axis="x", which="both", labelsize="x-small")
        if ymin == np.inf:
            ymin = 1e-6  # fallback
        ax.set_yscale("log")
        ax.set_ylim(ymin * 0.6, ymax * 2.0)

        ax.grid(True, which="both", axis="y", linestyle=":", alpha=0.4)
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(1, 10) * 0.1))

        # Process label
        texts = [
            r"$y_{1P}=y_{8P}=0.2$",
            r"$y_{1S}=y_{8S}=0.2$",
            r"$g_{1L}=g_{1R}=g_{8L}=g_{8R}=0.2$",
        ]
        # ax.text(0.02, 0.96, r"$y_{1P}=g_{8P}=0.2$", transform=ax.transAxes,
        #         fontsize=15, va="left")
        # ax.text(0.02, 0.96, r"$y_{1S}=g_{8S}=0.2$", transform=ax.transAxes,
        #         fontsize=15, va="left")        
        # ax.text(0.02, 0.96, r"$g_{1L}=g_{1R}=g_{8L}=g_{8R}=1$", transform=ax.transAxes,
        #         fontsize=15, va="left")
        x = 0.96          # ~allineato con la legenda a destra
        y0 = 0.80         # un po’ sotto la legenda
        dy = 0.06         # passo verticale

        for i, t in enumerate(texts):
            ax.text(x, y0-i*dy, t,
                    transform=ax.transAxes,
                    ha="right", va="top",
                    fontsize=12)      # o "x-small"

        # CMS label
        hep.cms.label("", data=False, loc=0, com="13", fontsize=20)

        # Legend
        ax.legend(handles=legend_handles, ncols=3, fontsize="xx-small", loc="upper right",
                  handletextpad=0.2, columnspacing=0.2)

        fig.tight_layout()
        fig.subplots_adjust(bottom=0.12)
        for ext in ("pdf", "png"):
            fig.savefig(outdir / f"xsec_vs_mass_{proc}.{ext}")
        plt.close(fig)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot σ(m) at c_0p2 (CMS style)")
    parser.add_argument("-o", "--outdir", default="/user/cgiordan/public_html/theory_plots/cross_section", type=Path,
                        help="Directory for output figures")
    args = parser.parse_args()
    output_dir = Path(args.outdir) / "all_models"  # Path instead of str
    make_plots(output_dir) 
    print(f"Plots written to → {args.outdir.resolve()}")
