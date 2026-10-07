"""Fig. 3(a): compute only the Delta=.2 QSH ribbon spectrum and edge labels.

Default: extract approved ribbon data. --recompute: diagonalize the same
width-64 ribbon, with no transport scan or plotting side effects.
"""
from __future__ import annotations
from pathlib import Path
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
from panels.common import calc_main
from panels.fig3._calc import PROJECT, ENERGIES, load, provenance

SOURCE = PROJECT / "data/raw/minimal_model_bands/ribbon_bands.npz"
WIDTH, NK, N_EDGE, EDGE_CUTOFF = 64, 481, 4, .55
DELTA = .2


def diagonalize():
    from scipy.linalg import eigh
    from threadpoolctl import threadpool_limits
    kx = np.linspace(-np.pi, np.pi, NK)
    spins = np.array([1, -1])
    tau_x = np.array([[0., 1.], [1., 0.]])
    tau_z = np.diag([1., -1.])
    hop = .5 * np.array([[1., -1.], [1., -1.]])
    position = np.repeat(np.arange(WIDTH, dtype=float), 2)
    shape = (len(spins), NK, 2 * WIDTH)
    energies, top, bottom = (np.zeros(shape) for _ in range(3))
    residual = 0.
    with threadpool_limits(limits=1):
        for si, spin in enumerate(spins):
            for ki, momentum in enumerate(kx):
                onsite = (spin*np.sin(momentum)+DELTA)*tau_x+(-1.+np.cos(momentum))*tau_z
                h = np.kron(np.eye(WIDTH), onsite)
                for y in range(WIDTH-1):
                    h[2*y:2*y+2, 2*y+2:2*y+4] = hop
                    h[2*y+2:2*y+4, 2*y:2*y+2] = hop.T
                e, vectors = eigh(h, check_finite=False, driver="evr")
                boundaries = np.r_[0, np.flatnonzero(np.diff(e)>1e-11)+1, len(e)]
                for first, last in zip(boundaries[:-1], boundaries[1:]):
                    if last-first > 1:
                        block = vectors[:, first:last]
                        _, rotation = eigh(block.T@(position[:, None]*block), check_finite=False)
                        vectors[:, first:last] = block@rotation
                energies[si, ki] = e
                rho = (vectors*vectors).reshape(WIDTH, 2, 2*WIDTH).sum(axis=1)
                top[si, ki] = rho[-N_EDGE:].sum(axis=0)
                bottom[si, ki] = rho[:N_EDGE].sum(axis=0)
                if ki in (NK//2, NK//2+32):
                    residual = max(residual, float(np.max(abs(h@vectors-vectors*e[None, :]))))
            print(f"Fig3(a): completed spin {spin}", flush=True)
    localized = top+bottom >= EDGE_CUTOFF
    for si, spin in enumerate(spins):
        projected_min = np.sqrt((spin*np.sin(kx)+DELTA)**2+(np.abs(-1.+np.cos(kx))-1.)**2)
        localized[si] &= np.abs(energies[si]) < projected_min[:, None]-1e-8
    labels = np.zeros(shape, dtype=np.int8)
    labels[localized & (top>bottom)] = 1
    labels[localized & (bottom>top)] = -1
    trs_error = float(np.max(abs(energies[0]-energies[1, ::-1])))
    if residual > 1e-9 or trs_error > 1e-10:
        raise RuntimeError("Ribbon eigenpair/TRS check failed")
    return dict(kx=kx, spins=spins, energies=energies, labels=labels,
                top_weight=top, bottom_weight=bottom), dict(
                eigenpair_residual_max=residual, TRS_spectrum_max_error=trs_error)


def calculate(recompute=False):
    if recompute:
        arrays, checks = diagonalize()
        sources = [Path(__file__)]
    else:
        archive = load(SOURCE)
        indices = np.flatnonzero(np.isclose(archive["deltas"], DELTA))
        if len(indices) != 1:
            raise ValueError("Exactly one Delta=.2 ribbon must be present")
        i = int(indices[0])
        arrays = {key: archive[key][i] for key in ("energies", "labels", "top_weight", "bottom_weight")}
        arrays.update(kx=archive["kx"], spins=archive["spins"])
        checks = {}
        sources = [SOURCE, SOURCE.with_name("validation.json")]
    # Plot-ready masks are prepared here; drawing never classifies eigenstates.
    for name, edge in (("bulk_curves", 0), ("top_curves", 1), ("bottom_curves", -1)):
        arrays[name] = np.where(arrays["labels"] == edge, arrays["energies"], np.nan).transpose(0, 2, 1).reshape(-1, arrays["kx"].size)
        if edge != 0:
            arrays[name] = arrays[name][np.any(np.isfinite(arrays[name]), axis=1)]
    # Preserve the production artist order too: intersections of equal-zorder
    # top/bottom branches then have precisely the same overlaid line colors.
    curves, kinds = [], []
    for spin in range(arrays["energies"].shape[0]):
        for state in range(arrays["energies"].shape[2]):
            energy, lab = arrays["energies"][spin, :, state], arrays["labels"][spin, :, state]
            for kind in (0, 1, -1):
                curve = np.where(lab == kind, energy, np.nan)
                if kind == 0 or np.any(np.isfinite(curve)):
                    curves.append(curve)
                    kinds.append(kind)
    arrays["curves"] = np.asarray(curves)
    arrays["curve_kinds"] = np.asarray(kinds, dtype=np.int8)
    arrays["selected_energies"] = ENERGIES.copy()
    arrays["delta"] = np.array(DELTA)
    metadata = provenance("a", recompute, sources, energy_unit="t", momentum_unit="1/a",
        width=WIDTH, Nk=NK, edge_rows=N_EDGE, edge_weight_cutoff=EDGE_CUTOFF,
        boundary="periodic x, open y; top y=W-1, bottom y=0", checks=checks,
        classification="outer-four-row weight and exclusion of spin-resolved projected bulk continuum")
    return arrays, metadata


if __name__ == "__main__":
    calc_main("fig3", "a", calculate)
