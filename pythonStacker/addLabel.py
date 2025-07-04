# add_label_png.py  –  timbra solo PNG/JPG
import numpy as np
from pathlib import Path
from PIL import Image
import matplotlib.pyplot as plt
from io import BytesIO

# ↳ mappa “modello → simbolo g” (modifica a piacere)
G_MAP = {
    "TopPhilicScalarOctet"  : r"y_{8S}",
    "TopPhilicVectorOctet"  : r"g_{8R}=g_{8L}",
    "TopPhilicScalarSinglet"  : r"y_{1S}",
    "TopPhilicVectorSinglet"  : r"g_{1R}=g_{1L}",
    "TopPhilicPseudoScalarSinglet":  r"y_{1P}",
    "TopPhilicPseudoScalarOctet":  r"y_{8P}"
}

def _pretty(name: str) -> str:
    m = name.split("Philic")[-1]
    for o, n in [("PseudoScalar","P"),("Scalar","S"),("Vector","V"),
                 ("Singlet",r"_{1}"),("Octet",r"_{8}")]:
        m = m.replace(o, n)
    return m

def _render_math(text, fontsize=22, dpi=300):
    fig = plt.figure(figsize=(0.01,0.01), dpi=dpi)
    fig.patch.set_alpha(0)
    ax = fig.add_axes([0,0,1,1]); ax.axis("off")
    t = ax.text(0,0,text, fontsize=fontsize, ha="left", va="bottom")
    fig.canvas.draw(); bbox = t.get_window_extent()
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, transparent=True,
                bbox_inches=bbox.transformed(fig.dpi_scale_trans.inverted()),
                pad_inches=0.0)
    plt.close(fig); buf.seek(0)
    return Image.open(buf)

def stamp_png(path,
              model: str,
              coupling,
              *,
              xfrac: float = 0.3,   # orizz. (0–1) – regola se serve
              margin_px: int = 5,    # distanza dal bordo grafico
              fontsize: int = 8,
              suffix: str = "_stamp",
              g_map: dict = G_MAP):
    """
    Crea <file>_stamp.png con la label latex in alto.
    """
    path = Path(path).expanduser().resolve()
    if path.suffix.lower() not in (".png",".jpg",".jpeg"):
        print("Ø  skip:", path.name); return

    latex = rf"${_pretty(model)}\;({g_map.get(model,'g')}={coupling})$"

    # --- carica immagine e trova top del grafico ---
    im = Image.open(path).convert("RGBA"); W,H = im.size
    arr = np.asarray(im)[:,:,:3]
    white = np.array([255,255,255], np.uint8)
    rows  = np.where(~np.all(arr == white, axis=2))[0]
    top   = int(rows[0]) if rows.size else 0

    # --- overlay latex ---
    ov = _render_math(latex, fontsize)
    ow, oh = ov.size
    x = int(W*xfrac - ow/2)
    y = top + margin_px

    im.alpha_composite(ov, (x, y))

    out_png = path.with_name(path.stem + suffix + path.suffix)
    im.save(out_png)
    print("✓  creato:", out_png)

# ---------------- esempio CLI rapido ----------------
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--model",    required=True)
    ap.add_argument("--coupling", required=True)
    ap.add_argument("images",     nargs="+")
    args = ap.parse_args()

    for img in args.images:
        stamp_png(img, args.model, args.coupling)
