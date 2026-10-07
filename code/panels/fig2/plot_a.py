"""Draw only Fig. 2(a), from a.npz; no numerical calculation is imported."""
from pathlib import Path
import sys
import numpy as np
from matplotlib.lines import Line2D
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from panels.common import plot_main, approved_style
from panels.fig2._style import TOP, BOTTOM, FONT, style, letter


def draw(ax, data, *, label=True):
    gap, ef = float(data["bulk_gap"]), float(data["EF"])
    ax.axhspan(-gap, gap, color="#F6F1E6", zorder=0)
    for i, edge in enumerate(data["edge_labels"]):
        ax.plot(data["kx"]/np.pi, data["strip_bands"][:, i],
                color=TOP if edge == 1 else BOTTOM if edge == -1 else "#D1D5D9",
                lw=1.1 if edge else .55, zorder=4 if edge else 1)
    ax.axhline(ef, lw=.6, ls="--", color="#444444")
    ax.axvline(float(data["kstar"])/np.pi, lw=.55, ls=":", color="#999999")
    for k, color in zip(data["fermi_momenta"], (TOP, BOTTOM)):
        ax.plot(np.array([-k, k])/np.pi, [ef, ef], "o", ms=3.2,
                mfc="white", mec=color, mew=.7)
    ax.set(xlim=(-1, 1), ylim=(-1.15, 1.15), xticks=[-1, 0, 1], yticks=[-1, 0, 1],
           xlabel=r"$k_xa_x/\pi$", ylabel=r"$E$ ($t_2$)")
    ax.text(.06, .05, r"$\Omega_{xy}=0$", transform=ax.transAxes, fontsize=FONT)
    handles = [Line2D([], [], color=TOP, lw=1.1), Line2D([], [], color=BOTTOM, lw=1.1)]
    ax.legend(handles, ["Top", "Bottom"], loc="upper center", ncol=2,
              frameon=True, facecolor="white", framealpha=1, edgecolor="none",
              fontsize=FONT, handlelength=.85, columnspacing=.8, handletextpad=.3, borderpad=.15)
    style(ax)
    if label:
        letter(ax, "a")
    return ax


if __name__ == "__main__":
    approved_style(FONT)
    plot_main("fig2", "a", draw, size_cm=(8.6, 7.5), axes_box=(.18, .16, .77, .75))
