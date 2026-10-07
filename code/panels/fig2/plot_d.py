"""Draw only Fig. 2(d): full-response Delta curves, with no scatter controls."""
from pathlib import Path
import sys
import numpy as np
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from panels.common import plot_main, approved_style
from panels.fig2._style import FONT, style, letter


def draw(ax, data, *, label=True):
    energies, deltas, values = (np.asarray(data[name]) for name in ("energies", "deltas", "kappa_full"))
    if values.shape != (len(deltas), len(energies)) or not np.isfinite(values).all():
        raise ValueError("Delta scan must contain finite dense energy curves")
    if not np.allclose(deltas, [-.2, 0, .1, .2]):
        raise ValueError("Unexpected production Delta scan values")
    if not np.allclose(values, -np.asarray(data["kH_full_code"])/4, rtol=0, atol=1e-14):
        raise ValueError("Full-bias convention normalization must be applied exactly once")
    colors = ["#8B6AB1", "#9099A2", "#D55E00", "#0072B2"]
    for delta, curve, color, linestyle in zip(deltas, values, colors, ("--", ":", "-", "-")):
        ax.plot(energies, curve, color=color, ls=linestyle, lw=1.1, marker=None,
                label=f"{delta:g}", zorder=3 if delta else 2)
    span = max(1., float(np.max(np.abs(values))))
    ax.set(xlim=(-.187, .187), xticks=[-.18, 0, .18], ylim=(-1.12*span, 1.65*span),
           yticks=[-1, 0, 1], xlabel=r"$E_F$ ($t_2$)", ylabel=r"$\widetilde\kappa_2$")
    ax.legend(title=r"$\Delta/t_2$", loc="upper center", ncol=2, fontsize=FONT,
              title_fontsize=FONT, handlelength=1.1, columnspacing=.7, handletextpad=.3,
              labelspacing=.12, borderpad=.15, frameon=True,
              facecolor="white", edgecolor="none", framealpha=1)
    style(ax)
    ax.yaxis.labelpad = .5
    if label:
        letter(ax, "d", horizontal_offset=-.04)
    return ax


if __name__ == "__main__":
    approved_style(FONT)
    plot_main("fig2", "d", draw, size_cm=(8.6, 7.5), axes_box=(.18, .16, .77, .75))
