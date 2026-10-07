"""Fig. 4(d): plot only saved QSH disorder statistics; no numerical backend."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from panels.common import plot_main
from panels.fig4 import _style


def draw(ax, data, *, label=True):
    return _style.draw_disorder(ax, data, "qsh", label=label)


if __name__ == "__main__":
    _style.style()
    plot_main("fig4", "d", draw, size_cm=(4.9, 4.7),
              axes_box=(.23, .22, 3.182/4.9, 3.216/4.7))
