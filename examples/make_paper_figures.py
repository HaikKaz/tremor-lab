"""Draw Figures 1 to 3 of the paper from docs/tremor_lab_regeneration_data.json.

    python examples/regenerate_paper_figures.py      # once: writes the JSON
    python examples/make_paper_figures.py [OUTPUT_DIR]

Needs matplotlib, which is not a dependency of the package:

    pip install matplotlib

Writes Figure_1, Figure_2 and Figure_3 as vector PDF and as 300 dpi PNG into
OUTPUT_DIR (default: figures/ at the top of the repository). Every plotted value is
read from the JSON, so a figure cannot disagree with the table it illustrates.

Colours are the Okabe-Ito set, chosen to stay distinguishable under the common forms
of colour-vision deficiency, and each catalogue keeps its colour across figures.
Meaning is also carried by marker shape and by direct labels, so the figures read
correctly in greyscale.
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "docs" / "tremor_lab_regeneration_data.json"

BLUE = "#0072B2"
VERMILLION = "#D55E00"
GREEN = "#009E73"
INK = "#1c1c1c"
MUTED = "#595959"
GRID = "#e4e4e4"

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7.5,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "axes.linewidth": 0.8,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "text.color": INK,
        "legend.fontsize": 7.5,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.dpi": 300,
    }
)

# One entry per catalogue: panel letter, title, colour, and the x range shown.
CATALOGUES = {
    "kahramanmaras": ("A", "Kahramanmaraş 2023 (KOERI)", BLUE, (2.9, 5.0)),
    "hectormine": ("B", "Hector Mine 1999 (ComCat), control", VERMILLION, (1.4, 4.1)),
    "ridgecrest": ("C", "Ridgecrest 2019 (ComCat)", GREEN, (1.0, 3.7)),
}


def lift(ax, points):
    """A transform that draws a marker `points` above its data position."""
    return ax.transData + ScaledTranslation(0, points / 72, ax.figure.dpi_scale_trans)


def curve_point(rows, threshold):
    for row in rows:
        if abs(row["threshold"] - threshold) < 1e-9:
            return row
    raise KeyError(f"no point at threshold {threshold}")


def note(ax, text, xy, xytext, **kwargs):
    """Annotation in neutral ink with a thin leader to the marked point."""
    ax.annotate(
        text,
        xy=xy,
        xytext=xytext,
        fontsize=7,
        color=INK,
        linespacing=1.15,
        arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.6, "shrinkA": 1},
        **kwargs,
    )


def figure_1(data):
    """b against the completeness threshold, one panel per catalogue."""
    widths = [hi - lo for _, _, _, (lo, hi) in CATALOGUES.values()]
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(7.5, 3.3),
        sharey=True,
        gridspec_kw={"width_ratios": widths, "wspace": 0.06},
    )
    fig.subplots_adjust(left=0.075, right=0.995, top=0.90, bottom=0.25)
    for ax, (key, (letter, title, color, xlim)) in zip(
        axes, CATALOGUES.items(), strict=True
    ):
        rows = data[key]["b_sweep"]
        x = np.array([r["threshold"] for r in rows])
        b = np.array([r["b"] for r in rows])
        s = np.array([r["sigma"] for r in rows])
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color=GRID, lw=0.6)
        ax.errorbar(
            x, b, yerr=s, fmt="none", ecolor=color, elinewidth=0.7, alpha=0.55, zorder=2
        )
        ax.plot(x, b, color=color, lw=1.0, alpha=0.9, zorder=3)
        ax.plot(x, b, "o", color=color, ms=2.4, zorder=4)

        mc = data[key]["mc"]
        top = curve_point(rows, mc["maxcurvature"])
        ax.plot(
            top["threshold"],
            top["b"],
            marker="v",
            ms=7.5,
            mfc=color,
            mec="white",
            mew=0.6,
            ls="none",
            transform=lift(ax, 6.5),
            zorder=6,
            clip_on=False,
        )
        if mc["b_stability"] is not None:
            stab = curve_point(rows, mc["b_stability"])
            ax.plot(
                stab["threshold"],
                stab["b"],
                marker="o",
                ms=9,
                mfc="none",
                mec=color,
                mew=1.5,
                ls="none",
                zorder=5,
            )

        ax.set_xlim(*xlim)
        ax.set_ylim(0.65, 1.45)
        ax.set_xticks(np.arange(np.ceil(xlim[0] * 2) / 2, xlim[1], 0.5))
        ax.set_title(
            f"{letter}  {title}",
            loc="left",
            fontsize=8,
            fontweight="bold",
            pad=5,
        )
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        if ax is not axes[0]:
            ax.tick_params(axis="y", length=0)
            ax.spines["left"].set_visible(False)

        if key == "kahramanmaras":
            top_b = top["b"]
            note(
                ax,
                f"maximum curvature\nM {mc['maxcurvature']}",
                (top["threshold"], top_b),
                (3.45, 0.70),
                ha="left",
                va="center",
            )
            stab = curve_point(rows, mc["b_stability"])
            note(
                ax,
                f"b-stability\nM {mc['b_stability']}",
                (stab["threshold"], stab["b"] + 0.03),
                (3.2, 1.22),
                ha="left",
                va="center",
            )
        elif key == "hectormine":
            ax.axvspan(2.2, 3.0, color=color, alpha=0.08, lw=0, zorder=0)
            ax.text(
                2.6,
                1.415,
                "plateau, M 2.2 to 3.0",
                ha="center",
                va="top",
                fontsize=7,
                color=INK,
            )
            note(
                ax,
                f"maximum curvature\nM {mc['maxcurvature']}, b {top['b']:.3f}",
                (top["threshold"], top["b"]),
                (1.95, 0.71),
                ha="left",
                va="center",
            )
            stab = curve_point(rows, mc["b_stability"])
            note(
                ax,
                f"b-stability\nM {mc['b_stability']}, b {stab['b']:.3f}",
                (stab["threshold"] - 0.03, stab["b"] + 0.03),
                (1.45, 1.17),
                ha="left",
                va="center",
            )
        else:
            note(
                ax,
                f"maximum curvature\nM {mc['maxcurvature']}",
                (top["threshold"], top["b"]),
                (1.5, 0.69),
                ha="left",
                va="center",
            )
            ax.text(
                1.05,
                1.30,
                "b-stability returns no value:\nb does not stabilise",
                ha="left",
                va="center",
                fontsize=7,
                color=INK,
            )
        ax.set_xlabel("Completeness threshold, M")

    axes[0].set_ylabel("Gutenberg-Richter b-value")
    handles = [
        Line2D([], [], marker="v", ls="none", ms=7, mfc=INK, mec="white"),
        Line2D([], [], marker="o", ls="none", ms=8, mfc="none", mec=INK, mew=1.4),
        Line2D([], [], color=INK, lw=0.8, alpha=0.6),
    ]
    labels = [
        "maximum-curvature Mc",
        "b-stability Mc",
        "Shi and Bolt (1982) standard error",
    ]
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, 0.0),
        handletextpad=0.5,
        columnspacing=2.0,
    )
    return fig


def figure_2(data):
    """Completeness magnitude by time band after the Kahramanmaraş mainshock."""
    bands = data["kahramanmaras"]["time_bands"]
    early = [b for b in bands if b["band_hi"] <= 5.0]
    late = next(b for b in bands if b["band_lo"] >= 30.0)
    threshold = 3.5

    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    fig.subplots_adjust(left=0.085, right=0.985, top=0.93, bottom=0.17)
    ax.set_xscale("log")
    ax.set_xlim(0.01, 260)
    ax.set_ylim(3.0, 4.7)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, lw=0.6)

    ax.axhline(threshold, color=VERMILLION, lw=1.1, ls=(0, (5, 3)), zorder=2)
    ax.text(
        258,
        threshold + 0.015,
        "M 3.5, published threshold",
        ha="right",
        va="bottom",
        fontsize=7,
        color=INK,
    )

    def count_label(lo, hi, mc, n):
        """Event count beside a band, below it where the M 3.5 line is above."""
        below = abs(threshold - mc) < 0.15 and mc < threshold
        ax.text(
            np.sqrt(lo * hi),
            mc - 0.05 if below else mc + 0.05,
            str(n),
            ha="center",
            va="top" if below else "bottom",
            fontsize=6.5,
            color=MUTED,
        )

    previous = None
    for i, band in enumerate(early):
        # The first band starts at t = 0, which a logarithmic axis cannot show, so
        # it is drawn from the left edge of the axis.
        lo = 0.01 if i == 0 else band["band_lo"]
        hi = band["band_hi"]
        mc = band["mc"]
        ax.plot([lo, hi], [mc, mc], color=BLUE, lw=2.8, solid_capstyle="butt", zorder=4)
        if previous is not None:
            ax.plot([lo, lo], [previous, mc], color=BLUE, lw=0.9, alpha=0.7, zorder=3)
        previous = mc
        if abs(band["band_lo"] - 0.38) > 1e-9:
            count_label(lo, hi, mc, band["n"])
    ax.text(
        0.0115,
        early[0]["mc"] - 0.07,
        "begins at the mainshock, t = 0",
        ha="left",
        va="top",
        fontsize=6.5,
        color=MUTED,
    )

    ax.plot(
        [late["band_lo"], late["band_hi"]],
        [late["mc"], late["mc"]],
        color=BLUE,
        lw=2.8,
        solid_capstyle="butt",
        zorder=4,
    )
    count_label(late["band_lo"], late["band_hi"], late["mc"], late["n"])
    ax.text(
        np.sqrt(5.0 * 30.0),
        3.215,
        "no band computed\nbetween 5 and 30 days",
        ha="center",
        va="bottom",
        fontsize=6.5,
        color=MUTED,
        style="italic",
        linespacing=1.15,
    )

    elbistan = next(b for b in early if abs(b["band_lo"] - 0.38) < 1e-9)
    ax.plot(
        0.38, elbistan["mc"], "o", ms=7, mfc=VERMILLION, mec="white", mew=1.0, zorder=6
    )
    ax.annotate(
        f"M 7.6 Elbistan at t = 0.380 d\nleast complete band, Mc {elbistan['mc']}, "
        f"n = {elbistan['n']}",
        xy=(0.38, elbistan["mc"]),
        xytext=(0.62, 4.52),
        fontsize=7,
        color=INK,
        ha="left",
        va="center",
        linespacing=1.15,
        arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.6, "shrinkA": 1},
    )

    ticks = [0.2, 0.38, 0.6, 1, 2, 5, 30, 180]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_yticks(np.arange(3.0, 4.61, 0.2))
    ax.set_xlabel("Time after the M 7.8 mainshock (days, logarithmic axis)")
    ax.set_ylabel("Completeness Mc by maximum curvature")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return fig


def figure_3(data):
    """Omori offset c against the threshold used for the fit, log scale on c."""
    fig, ax = plt.subplots(figsize=(5.6, 3.5))
    fig.subplots_adjust(left=0.12, right=0.975, top=0.95, bottom=0.14)
    ax.set_yscale("log")
    ax.set_xlim(1.5, 4.5)
    ax.set_ylim(0.1, 6.5)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, lw=0.6)
    series = [
        ("hectormine", "Hector Mine 1999 (ComCat), control", VERMILLION),
        ("kahramanmaras", "Kahramanmaraş 2023 (KOERI)", BLUE),
    ]
    handles = []
    for key, label, color in series:
        rows = data[key]["omori_sweep"]
        x = np.array([r["threshold"] for r in rows])
        c = np.array([r["c"] for r in rows])
        ax.plot(x, c, color=color, lw=1.4, zorder=3)
        ax.plot(x, c, "o", color=color, ms=4, zorder=4)
        (line,) = ax.plot([], [], color=color, lw=1.4, marker="o", ms=4, label=label)
        handles.append(line)
        ax.annotate(
            f"{c[0]:.2f} d",
            (x[0], c[0]),
            xytext=(6, 3),
            textcoords="offset points",
            ha="left",
            va="bottom",
            fontsize=7,
            color=INK,
        )
        ax.annotate(
            f"{c[-1]:.2f} d",
            (x[-1], c[-1]),
            xytext=(0, -6),
            textcoords="offset points",
            ha="center",
            va="top",
            fontsize=7,
            color=INK,
        )
    ax.set_xticks(np.arange(1.5, 4.51, 0.5))
    ax.set_yticks([0.1, 0.2, 0.5, 1, 2, 5])
    ax.set_yticklabels(["0.1", "0.2", "0.5", "1", "2", "5"])
    ax.minorticks_off()
    ax.set_xlabel("Completeness threshold used for the fit, M")
    ax.set_ylabel("Omori-Utsu offset c (days)")
    ax.legend(handles=handles, loc="upper right", handlelength=2.2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "outdir",
        nargs="?",
        default=str(REPO / "figures"),
        help="where to write Figure_1 to Figure_3 (default: figures/)",
    )
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    data = json.loads(DATA.read_text(encoding="utf-8"))
    for number, draw in enumerate((figure_1, figure_2, figure_3), start=1):
        fig = draw(data)
        for suffix in ("pdf", "png"):
            path = outdir / f"Figure_{number}.{suffix}"
            fig.savefig(path)
            print(f"wrote {path}")
        plt.close(fig)


if __name__ == "__main__":
    main()
