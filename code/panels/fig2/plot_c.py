"""Draw only Fig. 2(c): unmodified 30-cell left LPDOS with four contacts."""
from pathlib import Path
import sys
import numpy as np
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from panels.common import plot_main, approved_style
from panels.fig2._style import FONT, style, letter, add_contacts


def draw(ax, data, *, label=True):
    rho = np.asarray(data["rho1_cell"])
    if rho.shape != (30, 31) or np.any(rho < 0) or not np.isfinite(rho).all():
        raise ValueError("Panel c requires the real, finite 30 by 31 cell injectivity")
    cells, length = rho.shape
    image = ax.imshow(rho, origin="lower", extent=(-.5, length-.5, 0, cells),
                      aspect="auto", interpolation="nearest", cmap="Blues", vmin=0, vmax=.4)
    add_contacts(ax, length, cells)
    ax.set(xlim=(-1.4, length+.3), ylim=(-.4, cells+.4), xticks=[0, 15, 30],
           yticks=[0, 15, 30], xlabel=r"$x$ ($a_x$)", ylabel=r"$y$ ($a_y$)")
    x, y, width, height = ax.get_position().bounds
    cbax = ax.figure.add_axes([x+width+.015, y, .013, height])
    colorbar = ax.figure.colorbar(image, cax=cbax)
    colorbar.set_ticks([0, .4], labels=["0", "0.4"])
    colorbar.ax.tick_params(labelsize=FONT, length=2, width=.55, pad=1)
    colorbar.ax.set_title(r"$\rho_1$", fontsize=FONT, pad=4)
    style(ax, image=True)
    if label:
        letter(ax, "c")
    assert np.array_equal(np.asarray(image.get_array()), rho)
    return ax, cbax


if __name__ == "__main__":
    approved_style(FONT)
    plot_main("fig2", "c", draw, size_cm=(8.6, 7.5), axes_box=(.18, .16, .68, .75))
