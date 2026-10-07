"""Pure Matplotlib style shared by Figure 3 drawing scripts and assembly."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea, HPacker, TextArea
from matplotlib.ticker import AutoMinorLocator

TOP, BOTTOM = "#167CB0", "#D8792B"
ENERGY_COLORS = ("#0072B2", "#D55E00", "#CC79A7")
MARKERS = ("o", "s", "^")


def style():
    plt.rcParams.update({"font.family":"sans-serif", "font.sans-serif":["Arial","DejaVu Sans"],
        "font.size":8, "axes.labelsize":8, "axes.titlesize":8,
        "xtick.labelsize":8, "ytick.labelsize":8, "legend.fontsize":8,
        "mathtext.fontset":"dejavusans", "axes.linewidth":.6,
        "axes.labelpad":2, "pdf.fonttype":42, "ps.fonttype":42,
        "svg.fonttype":"none", "svg.hashsalt":"qsh-panel-workflow",
        "savefig.facecolor":"white", "figure.facecolor":"white"})


def axis_style(ax):
    ax.tick_params(direction="in", top=True, right=True, length=2.6, width=.6, pad=1.5)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="minor", length=1.5, width=.5)
    for spine in ax.spines.values():
        spine.set_linewidth(.6)
    ax.grid(False)


def panel_label(ax, letter):
    box = ax.get_position()
    ax.figure.text(box.x0-.012, box.y1+.012, f"({letter})", fontsize=8.5,
                   fontweight="bold", ha="left", va="bottom", zorder=100)


def inline_energy_legend(ax, handles):
    """Approved single-row EF/t prefix and color/marker key inside panel (b)."""
    parts = [TextArea(r"$E_F/t$", textprops={"fontsize":8})]
    for handle in handles:
        drawing = DrawingArea(6, 8, 0, 0)
        drawing.add_artist(Line2D([0,6], [4,4], color=handle.get_color(), lw=1))
        drawing.add_artist(Line2D([3], [4], color=handle.get_color(), ls="none",
                                  marker=handle.get_marker(), ms=3.2, mfc="white", mew=.8))
        label = TextArea(handle.get_label(), textprops={"fontsize":8})
        parts.append(HPacker(children=[drawing,label], align="center", pad=0, sep=1))
    row = HPacker(children=parts, align="center", pad=0, sep=2)
    legend = AnchoredOffsetbox(loc="center", child=row, pad=.05, borderpad=0,
        frameon=True, bbox_to_anchor=(.5,.5), bbox_transform=ax.transAxes)
    legend.patch.set(facecolor="white", edgecolor="none", alpha=1, linewidth=0)
    ax.add_artist(legend)
    return legend
