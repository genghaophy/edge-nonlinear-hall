"""Draw the device schematic using only prepared coordinates and labels."""
from pathlib import Path
import sys
from matplotlib.patches import Circle, Rectangle
if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from panels.common import approved_style, plot_main


def draw(ax, data, *, label=True):
    ax.set(xlim=data['x_limits'], ylim=data['y_limits'], aspect='equal')
    ax.set_axis_off()
    x, y, w, h = data['body']
    ax.add_patch(Rectangle((x, y), w, h, facecolor='#F4F5F6', edgecolor='#505A62', lw=.8, zorder=1))
    for (x, y, w, h), color in zip(data['edges'], ('#99C5DF', '#267BB2')):
        ax.add_patch(Rectangle((x, y), w, h, facecolor=color, edgecolor='none', zorder=2))
    for terminal, (x, y, w, h) in enumerate(data['terminals'], start=1):
        ax.add_patch(Rectangle((x, y), w, h, facecolor='#CCD2D7', edgecolor='#667581', lw=.7, zorder=3))
        ax.text(x+w/2, y+h/2, str(terminal), ha='center', va='center',
                fontsize=9, fontweight='bold', color='#303840', zorder=4)
    for x, y in zip(data['wire_x'], data['wire_y']):
        ax.plot(x, y, color='#5A6269', lw=.8)
    x, y, radius = data['voltmeter']
    ax.add_patch(Circle((x, y), radius, facecolor='white', edgecolor='#5A6269', lw=.8))
    for (x, y), text, size in zip(data['annotation_xy'], data['annotation_text'], data['annotation_fontsize']):
        ax.text(x, y, str(text), ha='center', va='center', fontsize=float(size))
    return ax


if __name__ == '__main__':
    approved_style(9.)
    plot_main('fig1', 'a', draw, size_cm=(8.6, 7.5), axes_box=(.025, .025, .95, .95))
