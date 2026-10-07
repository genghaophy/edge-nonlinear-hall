"""Calculate/extract only Fig. 2(b): fixed-k and k-integrated edge spectra."""
from pathlib import Path
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from panels.common import calc_main, source_record
from panels.fig2._data import RIBBON_SOURCE, load, provenance, select
from panels.fig2._ribbon import calculate_ribbon

KEYS = ("y_cell", "k_spectral_cell", "k_spectral_edge_cell", "edge_ldos_cell",
        "ldos_cell", "kstar", "eta", "spectral_energy", "EF", "N",
        "tx", "t1", "t2", "delta")


def calculate(recompute=False):
    if recompute:
        raw, checks = calculate_ribbon()
        arrays = select(raw, KEYS)
    else:
        arrays, checks = load(RIBBON_SOURCE, KEYS), {}
    if arrays["k_spectral_edge_cell"].shape != (len(arrays["y_cell"]), 2):
        raise ValueError("Spectral weights must have top/bottom edge columns")
    return arrays, provenance("b", RIBBON_SOURCE, recompute,
                              observable="cell-resolved fixed-(kx, EF) spectral density and k-integrated edge LDOS",
                              validation=checks, normalization="A/B cell sum, no edgewise normalization",
                              broadening="eta applies only to fixed-k spectrum; integrated DOS is analytic",
                              density_units="t2^-1",
                              calculation_helpers=[source_record(Path(__file__).with_name("_ribbon.py"))])


if __name__ == "__main__":
    calc_main("fig2", "b", calculate)
