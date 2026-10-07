"""Calculate/extract only Fig. 2(d): dense quadratic Hall-voltage Delta curves."""
from pathlib import Path
import sys
import numpy as np
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from panels.common import calc_main, source_record
from panels.fig2._data import DELTA_SOURCE, RM_ROOT, PARAMETERS, load, provenance

DELTAS = np.array([-.2, 0., .1, .2])
ENERGIES = np.unique(np.r_[np.linspace(-.18, .18, 161), .1])
STEP = 1e-5
KEYS = ("energies", "deltas", "kappa_full", "kH_full_code", "tx", "t1", "t2", "tc3", "tc4", "N", "L")


def _recompute():
    """Compute only four energy curves, reusing the verified transport kernel."""
    from threadpoolctl import threadpool_limits
    sys.path.insert(0, str(RM_ROOT/"code"))
    from compute_initial_study import Device
    from compute_floating_study import floating_response
    p = PARAMETERS
    raw = np.zeros((len(DELTAS), len(ENERGIES)))
    common_gap = min(np.sqrt((p["t2"]-p["t1"])**2+d**2)-2*abs(p["tx"]) for d in DELTAS)
    if max(abs(ENERGIES)) >= common_gap:
        raise ValueError("The Delta scan must remain inside the common global bulk gap")
    with threadpool_limits(limits=1):
        for j, delta in enumerate(DELTAS):
            dev = Device("four_x", cells=p["N"], length=p["L"], delta=float(delta),
                         tc=p["tc"], tx=p["tx"], t1=p["t1"], t2=p["t2"])
            for i, energy in enumerate(ENERGIES):
                raw[j, i] = floating_response(dev, float(energy), step=STEP)[0]["kH"]
                if i % 40 == 0 or i == len(ENERGIES)-1:
                    print(f"Fig2(d): Delta={delta:g}, {i+1}/{len(ENERGIES)}", flush=True)
    arrays = dict(energies=ENERGIES.copy(), deltas=DELTAS.copy(), kappa_full=-raw/4,
                  kH_full_code=raw, tc3=np.asarray(p["tc"]), tc4=np.asarray(p["tc"]))
    arrays.update({name: np.asarray(p[name]) for name in ("tx", "t1", "t2", "N", "L")})
    checks = dict(mirror_sign_error=float(np.max(abs(arrays["kappa_full"][0]+arrays["kappa_full"][-1]))),
                  symmetric_zero_error=float(np.max(abs(arrays["kappa_full"][1]))),
                  common_bulk_gap_halfwidth=float(common_gap))
    if max(checks["mirror_sign_error"], checks["symmetric_zero_error"]) > 2e-5:
        raise RuntimeError("Delta mirror/zero controls failed")
    return arrays, checks


def calculate(recompute=False):
    arrays, checks = _recompute() if recompute else (load(DELTA_SOURCE, KEYS), {})
    if arrays["kappa_full"].shape != (len(arrays["deltas"]), len(arrays["energies"])):
        raise ValueError("Delta curves have inconsistent axes")
    if not np.allclose(arrays["kappa_full"], -arrays["kH_full_code"]/4, rtol=0, atol=1e-14):
        raise ValueError("Apply the historical full-bias convention conversion exactly once")
    return arrays, provenance("d", DELTA_SOURCE, recompute,
                              observable="tilde kappa2=t2*kappa2/e, full floating four-terminal voltage response",
                              normalization="tilde kappa2=-kH_full_code/4 exactly once",
                              protocol="V1=+V/2,V2=-V/2,I3=I4=0,Vperp=V3-V4",
                              electrostatics="first-order local neutrality with internal and electrode feedback",
                              deltas_t2=arrays["deltas"].tolist(), energy_points=len(arrays["energies"]),
                              derivative_step=STEP, validation=checks,
                              calculation_helpers=[source_record(RM_ROOT/"code/compute_floating_study.py"),
                                                   source_record(RM_ROOT/"code/compute_initial_study.py")],
                              shared_solver_core="code/compute_floating_study.py")


if __name__ == "__main__":
    calc_main("fig2", "d", calculate)
