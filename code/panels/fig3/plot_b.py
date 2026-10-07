"""Fig. 3(b): draw saved fixed-EF finite-temperature curves, without solving."""
from pathlib import Path
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from panels.common import plot_main
from panels.fig3._style import ENERGY_COLORS, MARKERS, style, axis_style, panel_label, inline_energy_legend


def draw(ax, data, *, label=True):
    handles = []
    for energy, values, color, marker in zip(data["selected_energies"], data["kappa2_tilde"], ENERGY_COLORS, MARKERS):
        handle, = ax.plot(data["temperatures"], values, color=color, marker=marker,
                         ms=3.2, mfc="white", mew=.8, lw=1., label=f"{energy:g}")
        handles.append(handle)
    ax.set(xlim=(-.0007,.0207), xticks=[0,.01,.02], ylim=(.049,.0605), yticks=[.05,.06],
           xlabel=r"$k_BT/t$", ylabel=r"$\widetilde\kappa_2$")
    axis_style(ax)
    ax.yaxis.labelpad = -1
    ax.xaxis.labelpad = 0
    ax.tick_params(axis="x", which="major", pad=.8)
    inline_energy_legend(ax, handles)
    if label:
        panel_label(ax, "b")
    return ax


if __name__ == "__main__":
    style()
    plot_main("fig3", "b", draw, size_cm=(5.,5.), axes_box=(.19,.18,.76,.72))
