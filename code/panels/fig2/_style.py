"""Pure drawing helpers preserving the approved Figure 2 panel styling."""
from matplotlib.patches import Rectangle
from matplotlib.ticker import AutoMinorLocator

TOP, BOTTOM = "#0072B2", "#D55E00"
FONT = 8.2


def style(ax, image=False):
    for spine in ax.spines.values():
        spine.set_linewidth(.65)
    ax.tick_params(direction="in", top=True, right=True, length=2.4,
                   width=.65, labelsize=FONT, pad=1.8)
    ax.xaxis.label.set_fontsize(8.6)
    ax.yaxis.label.set_fontsize(8.6)
    ax.xaxis.labelpad = ax.yaxis.labelpad = 2
    if not image:
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def letter(ax, text, horizontal_offset=0):
    x, y, w, h = ax.get_position().bounds
    ax.figure.text(x-.004+horizontal_offset, y+h+.008, f"({text})",
                   fontsize=9, weight="bold", ha="left", va="bottom")


def add_contacts(ax, length, cells):
    """Four reservoirs in cell coordinates, exactly as in the approved map."""
    for x, y, width, height in ((-1.2, 0, .7, cells), (length-.5, 0, .7, cells),
                                 (length//2-.5, cells, 1, .3), (length//2-.5, -.3, 1, .3)):
        ax.add_patch(Rectangle((x, y), width, height, facecolor="#C9D1D7",
                               edgecolor="#667581", lw=.55, clip_on=False))
