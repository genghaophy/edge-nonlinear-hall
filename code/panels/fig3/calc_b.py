"""Fig. 3(b): fixed-EF finite-temperature floating Hall voltage.

Default extracts the four actual computed temperatures at EF/t=.3,.4,.5.
--recompute builds only the clean thermal spectrum/reduction needed here;
it does not calculate disorder or elastic-probe scans.
"""
from __future__ import annotations
from pathlib import Path
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
from panels.common import calc_main
from panels.fig3._calc import (ENERGIES, SOURCE, STEP, load, provenance,
                              backend, backend_sources, check_response)

TEMPERATURES = np.array([0., .005, .01, .02])


def calculate(recompute=False):
    if not recompute:
        archive = load(SOURCE)
        indices = []
        for energy in ENERGIES:
            matches = np.flatnonzero(np.isclose(archive["energies"], energy, rtol=0, atol=1e-14))
            if len(matches) != 1:
                raise ValueError("A target EF is not an exact archive grid column")
            indices.append(int(matches[0]))
        if not np.array_equal(archive["temperatures"], TEMPERATURES):
            raise ValueError("The four approved temperature points changed")
        raw = archive["kH_temperature"][:, indices].T.copy()
        gxx = archive["Gxx_temperature"][:, indices].T.copy()
        sources = [SOURCE, SOURCE.with_name("thermal_convergence.json")]
        checks = dict(archive_energy_indices=indices, interpolation=False)
    else:
        Device, quadrature, thermal_response, thread_limits = backend()
        history = []
        with thread_limits(limits=1):
            device = Device()
            zero_reports = [device.response_zero(float(e), step=STEP)[0] for e in ENERGIES]
            for report in zero_reports:
                check_response(report)
            previous = None
            accepted = False
            for order in (4, 8, 16):
                nodes, weights = quadrature(-.15, .79, .025, order)
                records = []
                for index, energy in enumerate(nodes):
                    point = device.spectral_record(float(energy), step=STEP)
                    if max(point["checks"].values()) > 5e-8:
                        raise RuntimeError("Thermal spectral gauge/reciprocity checks failed")
                    records.append(point)
                    if (index+1) % 40 == 0:
                        print(f"Fig3(b): order {order}, {index+1}/{len(nodes)} spectral nodes", flush=True)
                spectrum = {key: np.asarray([r[key] for r in records])
                            for key in ("C", "rho", "D", "qlead", "qwindow")}
                rows = [[thermal_response(device, nodes, weights, spectrum, float(mu), float(temp))[0]
                         for mu in ENERGIES] for temp in TEMPERATURES[1:]]
                for row in rows:
                    for report in row:
                        # Thermal reduction reports do not contain an on-shell
                        # reciprocity dictionary; retain floating/current checks.
                        if max(report["probe_linear_error"], report["probe_quadratic_error"]) > 1e-8:
                            raise RuntimeError("Finite-T floating currents are nonzero")
                values = np.asarray([[r["kH"] for r in row] for row in rows])
                mass_error = max(abs(r["thermal_mass"]-1.) for row in rows for r in row)
                error = None if previous is None else float(np.max(abs(values-previous)))
                scale = max(float(np.max(abs(values))), 1e-3)
                # The approved scan accepted order 16. Keep its final numerical
                # quadrature even when just these three EF values would satisfy
                # the coarse/fine criterion at a lower order.
                accepted = order == 16 and error is not None and error < 2e-4+.003*scale and mass_error < 2e-5
                history.append(dict(order=order, nodes=len(nodes), max_absolute_difference=error,
                                    missing_thermal_mass=mass_error, accepted=accepted))
                if accepted:
                    break
                previous = values
            if not accepted:
                raise RuntimeError("Declared thermal coarse/fine acceptance failed")
            all_rows = [zero_reports, *rows]
            raw = np.asarray([[r["kH"] for r in row] for row in all_rows]).T
            gxx = np.asarray([[r["Gxx"] for r in row] for row in all_rows]).T
        sources = [Path(__file__), *backend_sources()]
        checks = dict(thermal_convergence=history, interpolation=False)
    arrays = dict(temperatures=TEMPERATURES.copy(), selected_energies=ENERGIES.copy(),
                  kH_raw=raw, kappa2_tilde=-raw/4., Gxx=gxx)
    metadata = provenance("b", recompute, sources, checks=checks,
        axes_order=["Fermi energy", "temperature"], geometry=[31, 12], contact_coupling=.15,
        screening="thermally averaged injectivities before the local-neutrality ratio",
        temperature_unit="kBT/t", energy_unit="t", calculated_points_only=True)
    return arrays, metadata


if __name__ == "__main__":
    calc_main("fig3", "b", calculate)
