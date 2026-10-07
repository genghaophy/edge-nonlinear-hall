"""Exact finite-y Rice-Mele ribbon calculation; no plotting or transport."""
import numpy as np
from ._data import PARAMETERS, parameter_arrays


def calculate_ribbon():
    """Reproduce the production eigenproblem and displayed spectral weights.

    The k-integrated edge DOS is analytic, with no Lorentzian broadening.
    eta broadens only the fixed-(kx, EF) spectral function of panel b.
    """
    p = PARAMETERS
    cells, rows = p["N"], 2*p["N"]
    hy = np.diag(p["delta"] * (-1.)**np.arange(rows))
    for y in range(rows-1):
        hy[y, y+1] = hy[y+1, y] = p["t1"] if y % 2 == 0 else p["t2"]
    eig, psi = np.linalg.eigh(hy)
    prob = abs(psi)**2
    boundary = prob[-2:].sum(axis=0)-prob[:2].sum(axis=0)
    edge_ids = np.argsort(abs(eig))[:2]
    top = int(edge_ids[np.argmax(boundary[edge_ids])])
    bottom = int(edge_ids[np.argmin(boundary[edge_ids])])
    labels = np.zeros(rows, int)
    labels[top], labels[bottom] = 1, -1
    kx = np.linspace(-np.pi, np.pi, p["nk"])
    bands = eig[None, :]-2*p["tx"]*np.cos(kx[:, None])
    shifted = p["EF"]-eig
    propagating = abs(shifted) < 2*p["tx"]
    dos = np.zeros(rows)
    dos[propagating] = 1/(np.pi*np.sqrt((2*p["tx"])**2-shifted[propagating]**2))
    edge_ldos = prob[:, [top, bottom]]*dos[[top, bottom]][None, :]
    ek = eig-2*p["tx"]*np.cos(p["kstar"])
    spectral = p["eta"]/(np.pi*((p["EF"]-ek)**2+p["eta"]**2))
    cell_sum = lambda z: z.reshape(cells, 2, *z.shape[1:]).sum(axis=1)
    kf = np.arccos((eig[[top, bottom]]-p["EF"])/(2*p["tx"]))
    arrays = parameter_arrays()
    arrays.update(kx=kx, strip_bands=bands, edge_labels=labels,
                  bulk_gap=np.asarray(np.sqrt((p["t2"]-p["t1"])**2+p["delta"]**2)-2*p["tx"]),
                  y_cell=np.arange(cells)+.5, ribbon_eigenvalues=eig,
                  edge_profile_atomic=prob[:, [top, bottom]],
                  edge_ldos_cell=cell_sum(edge_ldos), ldos_cell=cell_sum(prob@dos),
                  k_spectral_cell=cell_sum(prob@spectral),
                  k_spectral_edge_cell=cell_sum(prob[:, [top, bottom]]*spectral[[top, bottom]][None, :]),
                  kstar=np.asarray(p["kstar"]), eta=np.asarray(p["eta"]),
                  spectral_energy=np.asarray(p["EF"]), fermi_momenta=kf,
                  fermi_velocities=2*p["tx"]*np.sin(kf),
                  top_bottom_edge_ids=np.array([top, bottom]))
    checks = dict(eigenvector_sum_error=float(np.max(abs(prob.sum(axis=0)-1))),
                  mirror_profile_error=float(np.max(abs(prob[:, top]-prob[::-1, bottom]))))
    if checks["eigenvector_sum_error"] > 1e-12 or checks["mirror_profile_error"] > 1e-8:
        raise RuntimeError("Ribbon eigenvector validation failed")
    return arrays, checks
