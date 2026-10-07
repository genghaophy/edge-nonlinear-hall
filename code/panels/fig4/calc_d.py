"""Fig. 4(d): QSH disorder mean and SD; reuse validated cache unless --recompute."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from panels.common import calc_main
from panels.fig4._calc import disorder_data


def calculate(recompute=False):
    return disorder_data("qsh", recompute=recompute)


if __name__ == "__main__":
    calc_main("fig4", "d", calculate)
