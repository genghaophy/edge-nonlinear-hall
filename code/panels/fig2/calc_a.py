"""Calculate/extract only Fig. 2(a): Rice-Mele ribbon dispersion."""
from pathlib import Path
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from panels.common import calc_main, source_record
from panels.fig2._data import RIBBON_SOURCE, load, provenance, select
from panels.fig2._ribbon import calculate_ribbon

KEYS = ("kx", "strip_bands", "edge_labels", "bulk_gap", "EF", "kstar",
        "fermi_momenta", "fermi_velocities", "tx", "t1", "t2", "delta", "N")


def calculate(recompute=False):
    if recompute:
        raw, checks = calculate_ribbon()
        arrays = select(raw, KEYS)
    else:
        arrays, checks = load(RIBBON_SOURCE, KEYS), {}
    if arrays["strip_bands"].shape != (len(arrays["kx"]), len(arrays["edge_labels"])):
        raise ValueError("Ribbon dispersion axes do not match")
    return arrays, provenance("a", RIBBON_SOURCE, recompute,
                              observable="exact ribbon eigenenergies; top/bottom edge branches",
                              validation=checks, normalization="none",
                              calculation_helpers=[source_record(Path(__file__).with_name("_ribbon.py"))])


if __name__ == "__main__":
    calc_main("fig2", "a", calculate)
