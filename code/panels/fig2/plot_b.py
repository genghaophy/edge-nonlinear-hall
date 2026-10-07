"""Draw only Fig. 2(b), including its approved integrated-spectrum inset."""
from pathlib import Path
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from panels.common import plot_main, approved_style
from panels.fig2._style import TOP, BOTTOM, FONT, style, letter


def draw(ax, data, *, label=True):
    inset = ax.inset_axes([.22, .27, .53, .52])
    ax.plot(data["y_cell"], data["k_spectral_cell"], color="#A9B0B7", lw=.6, zorder=1)
    for i, color in enumerate((TOP, BOTTOM)):
        ax.plot(data["y_cell"], data["k_spectral_edge_cell"][:, i], "-o",
                color=color, ms=3.2, mfc="white", mew=.7, lw=1, zorder=2)
        inset.plot(data["y_cell"], data["edge_ldos_cell"][:, i], "-o",
                   color=color, ms=3.2, mfc="white", mew=.7, lw=1)
    for target in (ax, inset):
        # Padding below zero keeps the nonnegative near-zero markers visible.
        target.set(xlim=(-.2, 8), ylim=(-.05, .85), xticks=[0, 4, 8], yticks=[0, .8])
        style(target)
    ax.set_ylabel(r"$\rho$ ($t_2^{-1}$)")
    ax.set_xlabel(r"$y$ ($a_y$)")
    ax.text(.04, .55, r"$k_x^*$", transform=ax.transAxes, fontsize=FONT)
    inset.set_title("Integrated", fontsize=FONT, pad=1.5)
    if label:
        letter(ax, "b")
    return ax, inset


if __name__ == "__main__":
    approved_style(FONT)
    plot_main("fig2", "b", draw, size_cm=(8.6, 7.5), axes_box=(.18, .16, .77, .75))
