"""Calculate/extract only Fig. 2(c): left injectivity of the 30-cell strip."""
from pathlib import Path
import sys
import numpy as np
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from panels.common import calc_main, source_record
from panels.fig2._data import MAP_SOURCE, RM_ROOT, PARAMETERS, load, provenance

KEYS = ("rho1_cell", "rho1_atomic", "x", "y_cell", "cells", "length", "atomic_rows",
        "N", "L", "EF", "tx", "t1", "t2", "delta", "tc")


def _recompute():
    """Run just this equilibrium map with the established Device Hamiltonian."""
    import kwant
    from threadpoolctl import threadpool_limits
    sys.path.insert(0, str(RM_ROOT/"code"))
    from compute_initial_study import Device
    p = PARAMETERS.copy()
    p["N"] = 30
    dev = Device("four_x", cells=p["N"], length=p["L"], tx=p["tx"],
                 t1=p["t1"], t2=p["t2"], delta=p["delta"], tc=p["tc"])
    with threadpool_limits(limits=1):
        wave = kwant.wave_function(dev.system, p["EF"], params=dev.params())
        nu = np.array([sum((dev.density(state, params=dev.params()) for state in wave(lead)),
                           np.zeros(len(dev.positions)))/(2*np.pi) for lead in range(4)])
        rho = kwant.ldos(dev.system, p["EF"], params=dev.params())
    error = float(np.max(abs(nu.sum(axis=0)-rho)))
    if error > 1e-10*max(1., float(np.max(rho))):
        raise RuntimeError("Sum of contact injectivities disagrees with total Kwant LDOS")
    grid = np.zeros((2*p["N"], p["L"]))
    for i in np.flatnonzero(dev.sample):
        x, y = dev.positions[i]
        grid[y, x] = nu[0, i]
    arrays = {name: np.asarray(p[name]) for name in ("N", "L", "EF", "tx", "t1", "t2", "delta", "tc")}
    arrays.update(rho1_atomic=grid, rho1_cell=grid.reshape(p["N"], 2, p["L"]).sum(axis=1),
                  x=np.arange(p["L"]), y_cell=np.arange(p["N"])+.5,
                  cells=np.asarray(p["N"]), length=np.asarray(p["L"]), atomic_rows=np.asarray(2*p["N"]))
    return arrays, dict(injectivity_sum_vs_ldos_max_error=error)


def calculate(recompute=False):
    arrays, checks = _recompute() if recompute else (load(MAP_SOURCE, KEYS), {})
    if arrays["rho1_cell"].shape != (30, 31) or np.any(arrays["rho1_cell"] < 0):
        raise ValueError("Expected the unrescaled nonnegative 30 by 31 left injectivity")
    metadata = provenance("c", MAP_SOURCE, recompute,
                          observable="rho1=sum_incoming_left |psi|^2/(2*pi)",
                          normalization="unit-flux incoming states; A/B cell sum; no interpolation or edgewise rescaling",
                          density_units="t2^-1", validation=checks,
                          shared_model_core="code/compute_initial_study.py",
                          calculation_helpers=[source_record(RM_ROOT/"code/compute_initial_study.py")],
                          contact_geometry="identical scalar metallic probes, tc=.15 at x=L//2")
    metadata["parameters"]["N"] = 30
    return arrays, metadata


if __name__ == "__main__":
    calc_main("fig2", "c", calculate)
