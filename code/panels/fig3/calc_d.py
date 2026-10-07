"""Fig. 3(d): independent spin-preserving elastic-probe scan.

--recompute calculates only gamma_phi dependence at zero T and disorder.
The paired virtual probes enforce separate up/down energy-resolved zero
currents. gamma_phi is a phenomenological linewidth, not a fitted lifetime.
"""
from __future__ import annotations
from pathlib import Path
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
from panels.common import calc_main
from panels.fig3._calc import (ENERGIES, SOURCE, EXTENSION, GAMMA, STEP,
    environment_archives, provenance, backend, backend_sources, check_response)


def calculate(recompute=False):
    if not recompute:
        original, extension, validation = environment_archives()
        raw = np.concatenate((original["kH_dephasing"], extension["kH_dephasing"]), axis=0)
        gxx = np.concatenate((original["Gxx_dephasing"], extension["Gxx_dephasing"]), axis=0)
        sources = [SOURCE, EXTENSION, validation]
        checks = dict(approved_rows_preserved=True)
    else:
        Device, _, _, thread_limits = backend()
        raw, gxx = (np.zeros((len(ENERGIES), len(GAMMA))) for _ in range(2))
        records = []
        with thread_limits(limits=1):
            for gi, gamma in enumerate(GAMMA):
                device = Device(gamma_phi=float(gamma), probe_stride=5)
                for ei, energy in enumerate(ENERGIES):
                    report, _ = device.response_zero(float(energy), step=STEP)
                    check_response(report)
                    raw[ei, gi], gxx[ei, gi] = report["kH"], report["Gxx"]
                    records.append(dict(gamma_phi=float(gamma), probe_positions=device.probe_positions,
                                        virtual_spin_channels=len(device.virtual_spins), **report))
                print(f"Fig3(d): gamma={gamma:g} completed", flush=True)
        sources = [Path(__file__), *backend_sources()]
        checks = dict(raw_response_records=records)
    paper = -raw/4.
    low, high = float(paper.min()), float(paper.max())
    span = max(high-low, 1e-6)
    ylim = np.array([low-.12*span, high+.12*span])
    arrays = dict(gamma_phi=GAMMA.copy(), selected_energies=ENERGIES.copy(),
        kH_raw=raw, kappa2_tilde=paper, Gxx=gxx, y_limits=ylim,
        y_ticks=np.arange(.05, ylim[1]+.000001, .01))
    metadata = provenance("d", recompute, sources, checks=checks, geometry=[31, 12],
        contact_coupling=.15, temperature=0., disorder=0.,
        virtual_probe_model="paired spin-preserving WBA elastic reservoirs; separate I_up(E)=I_down(E)=0",
        virtual_probe_positions=[[x,y] for x in (5,10,20,25) for y in (0,1,10,11)],
        gamma_interpretation="energy linewidth, not calibrated hbar/tau_phi",
        axes_order=["Fermi energy", "elastic-probe linewidth"])
    return arrays, metadata


if __name__ == "__main__":
    calc_main("fig3", "d", calculate)
