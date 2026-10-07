"""Pure Matplotlib helpers matching the current formal Figure 4 artwork."""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea, HPacker, TextArea
from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter
import numpy as np

COLORS = ("#0072B2", "#D55E00", "#CC79A7")
MARKERS = ("o", "s", "^")


def style():
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8, "mathtext.fontset": "dejavusans", "axes.linewidth": .6,
        "axes.labelpad": 2, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white", "figure.facecolor": "white"})


def _ticks(ax):
    ax.tick_params(direction="in", top=True, right=True, length=2.8, width=.6, pad=2)
    ax.tick_params(which="minor", direction="in", top=True, right=True, length=1.6, width=.5)
    ax.grid(False)


def _energy_key(ax, energies, unit):
    parts = [TextArea(r"$E_F/" + unit + "$", textprops={"fontsize": 8})]
    for energy, color, marker in zip(energies, COLORS, MARKERS):
        drawing = DrawingArea(6, 8, 0, 0)
        drawing.add_artist(Line2D([0, 6], [4, 4], color=color, lw=.95))
        drawing.add_artist(Line2D([3], [4], color=color, ls="none", marker=marker,
                                   ms=3, mfc="white", mew=.75))
        text = f"{energy:.2f}".removeprefix("0") if unit == "t_2" else f"{energy:.1f}"
        parts.append(HPacker(children=[drawing, TextArea(text, textprops={"fontsize": 8})],
                             align="center", pad=0, sep=1.5))
    box = AnchoredOffsetbox(loc="center", child=HPacker(children=parts, align="center", pad=0, sep=3),
        pad=.05, borderpad=0, frameon=True, bbox_to_anchor=(.5, .865), bbox_transform=ax.transAxes)
    box.patch.set(facecolor="white", edgecolor="none", alpha=1)
    ax.add_artist(box)
    return box


def draw_coupling(ax, data, model, *, label=True):
    """Draw a or b from saved arrays only; return its auxiliary energy key."""
    unit = "t_2" if model == "rm" else "t"
    energies, hopping = np.asarray(data["energies"]), np.asarray(data["tau"])
    full, bare = np.asarray(data["equal_kappa_full"]), np.asarray(data["bare_kappa_predicted"])
    if energies.shape != (3,) or full.shape != (3, len(hopping)) or bare.shape != (3,):
        raise ValueError("Expected three Fermi energies and matching saved coupling curves.")
    if np.any(hopping <= 0) or not all(np.isfinite(x).all() for x in (energies, hopping, full, bare)):
        raise ValueError("Nonfinite or nonpositive coupling-panel data.")
    for index, (color, marker) in enumerate(zip(COLORS, MARKERS)):
        ax.plot(hopping, full[index], color=color, lw=1.05, marker=marker, markevery=5,
                ms=3, mfc="white", mew=.75, zorder=3)
        ax.axhline(bare[index], color=color, lw=.85, ls=(0, (4, 2.5)), zorder=2)
    ax.set_xscale("log")
    ax.set_xlim(.004, 1.12)
    ax.set_xticks([.01, .1, 1], labels=["0.01", "0.1", "1"])
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs=(2, 5)))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    _ticks(ax)
    ax.set_xlabel(r"$\tau/" + unit + "$", labelpad=2)
    if model == "rm":
        ax.set_ylim(-1.4, .95)
        ax.set_yticks([-1, -.5, 0, .5], labels=["-1", "-0.5", "0", "0.5"])
        ax.set_ylabel(r"$\widetilde\kappa_2$", labelpad=1)
    else:
        ax.set_ylim(.018, .073)
        ax.set_yticks([.02, .04, .06], labels=["0.02", "0.04", "0.06"])
    if label:
        heading = "(a)  Rice-Mele" if model == "rm" else "(b)  QSH"
        ax.text(0., 1.025, heading, transform=ax.transAxes, va="bottom", ha="left",
                fontsize=8.5, fontweight="bold", clip_on=False)
    return _energy_key(ax, energies, unit)


def draw_disorder(ax, data, model, *, label=True):
    """Draw c or d using saved mean and sample-SD arrays, without transport imports."""
    unit = "t_2" if model == "rm" else "t"
    strength = np.asarray(data["strengths"])
    mean, spread = np.asarray(data["four_normalized_mean"]), np.asarray(data["four_normalized_std"])
    if strength.shape != (4,) or mean.shape != (4,) or spread.shape != (4,):
        raise ValueError("Expected one model's four-strength disorder panel.")
    if not all(np.isfinite(x).all() for x in (strength, mean, spread)) or np.any(spread < 0):
        raise ValueError("Invalid saved disorder statistics.")
    color = COLORS[1]
    ax.fill_between(strength, mean-spread, mean+spread, color=color, alpha=.20, linewidth=0, zorder=1)
    ax.plot(strength, mean, color=color, lw=1.15, zorder=3)
    ax.axhline(1., color=".48", lw=.70, ls=(0, (3, 2.5)), zorder=2)
    ax.set(xlim=(0., float(strength[-1])), ylim=(-4., 5.))
    if model == "rm":
        ax.set_xticks([0., .05, .10], labels=["0", "0.05", "0.10"])
        ax.set_ylabel(r"$\kappa_2(W_{\rm dis})/\kappa_2(0)$", labelpad=1)
    else:
        ax.set_xticks([0., .25, .50], labels=["0", "0.25", "0.50"])
    ax.set_yticks([-4, -2, 0, 2, 4])
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    _ticks(ax)
    ax.set_xlabel(r"$W_{\rm dis}/" + unit + "$", labelpad=2)
    if label:
        ax.text(0., 1.025, "(c)" if model == "rm" else "(d)", transform=ax.transAxes,
                va="bottom", ha="left", fontsize=8.5, fontweight="bold", clip_on=False)
    text = (rf"$E_F/{unit}={float(data['energy']):.2f}$" + "\n"
            + rf"$\tau/{unit}={float(data['tau']):.3f}$")
    ax.text(.035, .055, text, transform=ax.transAxes, ha="left", va="bottom", fontsize=8, linespacing=1.2)
