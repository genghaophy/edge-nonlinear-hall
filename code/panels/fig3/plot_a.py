"""Fig. 3(a): draw only the saved, classified QSH ribbon curves."""
from pathlib import Path
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from matplotlib.lines import Line2D
from panels.common import plot_main
from panels.fig3._style import TOP, BOTTOM, ENERGY_COLORS, style, axis_style, panel_label


def draw(ax, data, *, label=True):
    line_styles = {0:dict(color="#D4D9DD", lw=.5, zorder=1),
                  1:dict(color=TOP, ls="-", lw=1.05, zorder=3),
                 -1:dict(color=BOTTOM, ls=(0,(3,1.8)), lw=1.05, zorder=3)}
    for curve, kind in zip(data["curves"], data["curve_kinds"]):
        ax.plot(data["kx"], curve, **line_styles[kind])
    for energy, color in zip(data["selected_energies"], ENERGY_COLORS):
        ax.axhline(energy, color=color, lw=.55, ls=(0,(2.5,2.5)), zorder=2)
    ax.set(xlim=(-1.25,1.25), ylim=(-1.45,1.45), xticks=[-1,0,1], yticks=[-1,0,1],
           xlabel=r"$k_xa$", ylabel=r"$E/t$")
    axis_style(ax)
    ax.yaxis.labelpad = 0
    ax.xaxis.labelpad = 0
    ax.tick_params(axis="x", which="major", pad=.8)
    handles = [Line2D([],[],color=TOP,lw=1.05,label="Top"),
               Line2D([],[],color=BOTTOM,lw=1.05,ls=(0,(3,1.8)),label="Bottom")]
    ax.legend(handles=handles, loc="upper center", ncol=2, frameon=True,
              facecolor="white", edgecolor="none", framealpha=1,
              handlelength=.75, handletextpad=.25, columnspacing=.4,
              borderpad=.05, borderaxespad=.2)
    if label:
        panel_label(ax, "a")
    return ax


if __name__ == "__main__":
    style()
    plot_main("fig3", "a", draw, size_cm=(5.,5.), axes_box=(.19,.18,.76,.72))
