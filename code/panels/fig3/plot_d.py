"""Fig. 3(d): draw saved spin-preserving elastic-probe response."""
from pathlib import Path
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from panels.common import plot_main
from panels.fig3._style import ENERGY_COLORS, MARKERS, style, axis_style, panel_label


def draw(ax, data, *, label=True):
    for energy, values, color, marker in zip(data["selected_energies"], data["kappa2_tilde"], ENERGY_COLORS, MARKERS):
        ax.plot(data["gamma_phi"], values, color=color, marker=marker,
                ms=3.2, mfc="white", mew=.8, lw=1, label=f"{energy:g}")
    ax.set(xlim=(-.04,1.04), xticks=[0,.5,1], ylim=data["y_limits"], yticks=data["y_ticks"],
           xlabel=r"$\gamma_\phi/t$", ylabel=r"$\widetilde\kappa_2$")
    axis_style(ax)
    ax.yaxis.labelpad = -2
    if label:
        panel_label(ax, "d")
    return ax


if __name__ == "__main__":
    style()
    plot_main("fig3", "d", draw, size_cm=(5.,5.), axes_box=(.19,.18,.76,.72))
