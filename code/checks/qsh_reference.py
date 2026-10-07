"""Independent QSH clean, disorder, elastic-probe and finite-temperature checks.

Calls the original verified selected-point algorithms. It never refreshes
or overwrites production caches. Finite-temperature validation uses a 9 x 6
device and does not establish size convergence of the 31 x 12 device.
"""
from pathlib import Path
import argparse
import json
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from checks._common import DATA, OUTPUT, output_directory, protect_inputs, write_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT / "qsh_reference")
    parser.add_argument("--skip-thermal", action="store_true", help="Skip the independent 9 x 6 finite-temperature current integration")
    args = parser.parse_args()
    import numpy as np
    from compute_qsh_figure3 import validate_device, validate_thermal
    from qsh_environment_response import EnvironmentDevice
    from threadpoolctl import threadpool_limits
    reference_path = DATA / "supplementary" / "floating_probe_reference.json"
    with protect_inputs([reference_path]) as fingerprints, threadpool_limits(limits=1):
        clean = EnvironmentDevice()
        checks = {"clean": validate_device(clean)}
        equilibrium = clean.effective(.4)
        checks["independent_kwant_matrix_error"] = float(np.max(abs(equilibrium - clean.base.matrix(.4))))
        reference = json.loads(reference_path.read_text(encoding="utf-8"))["report"]
        checks["historical_reference_kH_error_raw"] = abs(checks["clean"]["kH"] - reference["kH"])
        checks["gamma_zero_matrix_error"] = float(np.max(abs(EnvironmentDevice(gamma_phi=0.).effective(.4) - equilibrium)))
        if checks["independent_kwant_matrix_error"] > 1e-9 or checks["gamma_zero_matrix_error"] > 1e-9:
            raise RuntimeError("Independent QSH matrix or zero-coupling check failed")
        if checks["historical_reference_kH_error_raw"] > 1e-7:
            raise RuntimeError("Selected-point QSH result differs from the historical reference")
        checks["disorder"] = validate_device(EnvironmentDevice(disorder=.5, seed=1700))
        checks["elastic_probes"] = validate_device(EnvironmentDevice(gamma_phi=.1))
        if not args.skip_thermal:
            checks["finite_temperature_9x6"] = validate_thermal()
    path = output_directory(args.output_dir) / "report.json"
    write_report(path, {"status": "passed", "input_sha256": fingerprints,
                        "conversion": "All kH fields remain historical raw amplitudes; paper tilde_kappa2=-kH/4 and errors divide by4",
                        "scope": "Selected samples within first-order local neutrality; finite-T check is not device-size convergence",
                        "checks": checks})
    print(path)


if __name__ == "__main__":
    main()
