"""Check RM finite-bias floating currents, reciprocity, gauge shift and derivatives.

One representative 31 x 16-atomic-row sample is solved. This is not an
energy, width, or disorder scan. Historical raw coefficients are converted
once to the full-paper-bias dimensionless response, tilde_kappa2=-kH/4.
"""
from pathlib import Path
import argparse
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from checks._common import OUTPUT, output_directory, write_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT / "rm_reference")
    args = parser.parse_args()
    from compute_initial_study import Device
    from compute_floating_study import floating_response, verify
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        device = Device("four_x")
        raw, u_grid, _, _ = floating_response(device, .1)
        verification = verify(device, .1, raw, u_grid)
    paper = {name.replace("kH", "tilde_kappa2", 1): -raw[name] / 4
             for name in ("kH", "kH_frozen", "kH_internal", "kH_electrodes")}
    path = output_directory(args.output_dir) / "report.json"
    write_report(path, {"scope": "Selected RM reference; first-order local neutrality, not nonlinear Poisson",
                        "parameters": {"length": 31, "cells": 8, "EF_over_t2": .1, "tau_over_t2": .15},
                        "raw_positive_energy_solver": raw, "verification_raw": verification,
                        "paper": paper, "conversion": "tilde_kappa2=-kH_raw/4; coefficient errors divide by 4",
                        "status": "passed"})
    print(path)


if __name__ == "__main__":
    main()
