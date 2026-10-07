"""Fig. 3(c): draw saved scalar-disorder means and SEM bars."""
from pathlib import Path
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from panels.common import plot_main
from panels.fig3._style import ENERGY_COLORS, MARKERS, style, axis_style, panel_label


def draw(ax, data, *, label=True):
    for energy, mean, sem, color, marker in zip(data["selected_energies"], data["kappa2_tilde_mean"],
            data["kappa2_tilde_sem"], ENERGY_COLORS, MARKERS):
        ax.errorbar(data["strengths"], mean, yerr=sem, color=color, lw=1,
            marker=marker, ms=3.2, mfc="white", mew=.8, elinewidth=.65,
            capsize=1.8, capthick=.65, label=f"{energy:g}")
    ax.set(xlim=(-.04,1.04), xticks=[0,.5,1], ylim=data["y_limits"], yticks=data["y_ticks"],
           xlabel=r"$W_{\rm dis}/t$", ylabel=r"$\langle\widetilde\kappa_2\rangle$")
    axis_style(ax)
    ax.yaxis.labelpad = 3
    if label:
        panel_label(ax, "c")
    return ax


if __name__ == "__main__":
    style()
    plot_main("fig3", "c", draw, size_cm=(5.6,5.), axes_box=(.23,.18,.72,.72))
