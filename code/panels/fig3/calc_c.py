"""Fig. 3(c): independent scalar-disorder scan, with eight-sample SEM.

--recompute calculates only zero-temperature nonmagnetic disorder response
at the three EF values. No thermal or elastic-probe scan is dispatched.
"""
from __future__ import annotations
from pathlib import Path
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
from panels.common import calc_main
from panels.fig3._calc import (ENERGIES, SOURCE, EXTENSION, STRENGTHS, SEEDS, STEP,
    environment_archives, provenance, backend, backend_sources, check_response)


def calculate(recompute=False):
    if not recompute:
        original, extension, validation = environment_archives()
        mean, sem, sd, gxx = [np.concatenate((original[key], extension[key]), axis=0)
            for key in ("kH_disorder_mean", "kH_disorder_sem", "kH_disorder_sd", "Gxx_disorder_mean")]
        counts = original["disorder_samples"].copy()
        sources = [SOURCE, EXTENSION, validation]
        checks = dict(approved_rows_preserved=True)
    else:
        Device, _, _, thread_limits = backend()
        mean, sem, sd, gxx = (np.zeros((len(ENERGIES), len(STRENGTHS))) for _ in range(4))
        counts = np.zeros(len(STRENGTHS), dtype=int)
        records = []
        with thread_limits(limits=1):
            for wi, strength in enumerate(STRENGTHS):
                seeds = [0] if strength == 0 else SEEDS
                values, gxvalues = [], []
                for seed in seeds:
                    device = Device(disorder=float(strength), seed=int(seed))
                    sample, gxsample = [], []
                    for energy in ENERGIES:
                        report, _ = device.response_zero(float(energy), step=STEP)
                        check_response(report)
                        sample.append(report["kH"])
                        gxsample.append(report["Gxx"])
                        records.append(dict(strength=float(strength), seed=int(seed), **report))
                    values.append(sample)
                    gxvalues.append(gxsample)
                values = np.asarray(values)
                counts[wi] = len(seeds)
                mean[:, wi] = values.mean(axis=0)
                gxx[:, wi] = np.asarray(gxvalues).mean(axis=0)
                sd[:, wi] = values.std(axis=0, ddof=1) if len(seeds)>1 else 0.
                sem[:, wi] = sd[:, wi]/np.sqrt(len(seeds))
                print(f"Fig3(c): W={strength:g}, {len(seeds)} realizations completed", flush=True)
        sources = [Path(__file__), *backend_sources()]
        checks = dict(raw_response_records=records)
    paper_mean, paper_sem, paper_sd = -mean/4., sem/4., sd/4.
    lower, upper = paper_mean-paper_sem, paper_mean+paper_sem
    c_low = min(.039, np.floor((lower.min()-.001)/.005)*.005)
    c_high = max(.063, np.ceil((upper.max()+.001)/.005)*.005)
    arrays = dict(strengths=STRENGTHS.copy(), selected_energies=ENERGIES.copy(),
        kH_raw_mean=mean, kH_raw_sem=sem, kH_raw_sd=sd,
        kappa2_tilde_mean=paper_mean, kappa2_tilde_sem=paper_sem, kappa2_tilde_sd=paper_sd,
        samples=counts, seeds=SEEDS.copy(), Gxx_mean=gxx,
        y_limits=np.array([c_low, c_high]), y_ticks=np.arange(.04, c_high+.000001, .02))
    metadata = provenance("c", recompute, sources, checks=checks, geometry=[31, 12],
        contact_coupling=.15, temperature=0., disorder="sample-only IID scalar onsite W*uniform[-1/2,1/2]",
        disorder_seeds=SEEDS.tolist(), disorder_samples=counts.tolist(),
        error_statistic="SEM=sample SD(ddof=1)/sqrt(n); zero for one clean point",
        screening="first-order local neutrality; no Poisson solver", axes_order=["Fermi energy", "disorder strength"])
    return arrays, metadata


if __name__ == "__main__":
    calc_main("fig3", "c", calculate)
