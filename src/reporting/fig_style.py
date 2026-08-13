"""Publication figure standard for this repo.

Implements the `publication-figures` skill contract referenced by CLAUDE.md:
  - vector PDF + 600 DPI PNG for every figure
  - Arial, fonts embedded in the PDF (fonttype 42, editable text)
  - Okabe-Ito colorblind-safe palette, legible in grayscale
  - error bars, n per group, and significance markers on comparisons
  - no chartjunk: minimal gridlines, direct labels over legends where feasible

Usage in a notebook or script:

    import sys; sys.path.insert(0, str(PROJ / 'src/reporting'))
    from fig_style import apply_pub_style, save_figure, OKABE_ITO
    apply_pub_style()
    ...
    save_figure(fig, FIG / 'aim1_forest')   # writes .pdf (vector) + .png (600 dpi)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

# Okabe-Ito colorblind-safe palette, in the canonical order used by the skill.
OKABE_ITO = [
    "#000000",  # black
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
]

# Alias kept so callers can import either name.
PALETTE = OKABE_ITO

# Two-group comparisons (e.g. brain-met vs no brain-met). Blue/vermillion keeps
# a wide luminance gap so the pair survives grayscale conversion.
PALETTE_BINARY = ["#0072B2", "#D55E00"]

FONT_STACK = ["Arial", "Helvetica", "Helvetica Neue", "DejaVu Sans"]

EXPORT_DPI = 600


def apply_pub_style(base_size: int = 9) -> None:
    """Install the repo figure standard. Call before creating any figure."""
    mpl.rcParams.update({
        # --- typography: Arial, >= 8 pt everywhere ---
        "font.family": "sans-serif",
        "font.sans-serif": FONT_STACK,
        "font.size": base_size,
        "axes.titlesize": base_size + 1,
        "axes.labelsize": base_size,
        "xtick.labelsize": base_size - 1,
        "ytick.labelsize": base_size - 1,
        "legend.fontsize": base_size - 1,
        "figure.titlesize": base_size + 2,

        # --- no chartjunk ---
        "axes.grid": False,
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.spines.left": True,
        "axes.spines.bottom": True,
        "axes.linewidth": 0.8,
        "axes.edgecolor": "#000000",

        # --- ticks outward ---
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,

        # --- white background; 600 dpi on export ---
        "figure.facecolor": "white",
        "figure.edgecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.edgecolor": "white",
        "figure.dpi": 150,          # on-screen only
        "savefig.dpi": EXPORT_DPI,
        "savefig.pad_inches": 0.02,

        # --- constrained layout per the skill snippet ---
        "figure.constrained_layout.use": True,

        # --- legend ---
        "legend.frameon": False,
        "legend.handlelength": 1.6,
        "legend.columnspacing": 1.2,

        # --- lines / error bars ---
        "lines.linewidth": 1.4,
        "lines.markersize": 4.0,
        "errorbar.capsize": 2.5,

        # --- embed real TrueType fonts in vector output (editable, no Type-3) ---
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    mpl.rcParams["axes.prop_cycle"] = mpl.cycler(color=OKABE_ITO)

    # Keep seaborn from re-adding a grid or grey axes background.
    try:
        import seaborn as sns
        sns.set_style("ticks", rc={
            "axes.grid": False,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
        })
        sns.set_palette(OKABE_ITO)
    except ImportError:
        pass


def save_figure(fig, stem, dpi: int = EXPORT_DPI) -> list[Path]:
    """Export `fig` as vector PDF plus 600 DPI PNG.

    `stem` is a path without an extension; parent dirs are created. Returns the
    written paths. Both formats are always produced - PNG-only export violates
    the repo figure standard.
    """
    stem = Path(stem)
    if stem.suffix in {".png", ".pdf", ".svg", ".jpg", ".jpeg", ".tif", ".tiff"}:
        stem = stem.with_suffix("")
    stem.parent.mkdir(parents=True, exist_ok=True)

    written = []
    for ext, kw in (("pdf", {}), ("png", {"dpi": dpi})):
        out = stem.with_suffix(f".{ext}")
        fig.savefig(out, facecolor="white", edgecolor="white", **kw)
        written.append(out)
    return written


# Backwards-compatible alias for the earlier helper name.
savefig_pub = save_figure


def finalize(ax=None):
    """Strip any grid or grey facecolor a helper (seaborn/lifelines) re-added."""
    ax = ax or plt.gca()
    ax.grid(False)
    ax.set_facecolor("white")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return ax


def p_to_stars(p: float) -> str:
    """Significance marker for an annotated comparison. 'ns' when not significant."""
    if p is None:
        return ""
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def annotate_n(ax, labels_to_n: dict, y: float = -0.08) -> None:
    """Write 'n=<k>' under each categorical x tick.

    The skill requires the per-group n to appear on the figure itself rather
    than only in the caption.
    """
    ticks = ax.get_xticks()
    labels = [t.get_text() for t in ax.get_xticklabels()]
    for x, lab in zip(ticks, labels):
        n = labels_to_n.get(lab)
        if n is None:
            continue
        ax.annotate(f"n={n}", xy=(x, y), xycoords=("data", "axes fraction"),
                    ha="center", va="top", fontsize=mpl.rcParams["font.size"] - 2)


def bracket(ax, x1: float, x2: float, y: float, text: str,
            h: float | None = None, lw: float = 0.8) -> None:
    """Draw a significance bracket spanning x1..x2 at height y, labelled `text`."""
    if h is None:
        lo, hi = ax.get_ylim()
        h = (hi - lo) * 0.02
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=lw, c="black",
            solid_joinstyle="miter", clip_on=False)
    ax.text((x1 + x2) / 2, y + h, text, ha="center", va="bottom",
            fontsize=mpl.rcParams["font.size"] - 1)
