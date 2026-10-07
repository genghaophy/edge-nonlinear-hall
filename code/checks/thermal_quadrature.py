"""Compare orders 4, 8, 16 on the article's fixed thermal integration window.

Default reduces exact shipped spectra; --recompute regenerates equilibrium
spectra using the existing verified kernel. Fresh arrays go to outputs/checks.
This is quadrature convergence on [-0.15,0.79]t, not window enlargement or
device-size convergence. It does not replace the finite-T direct-current test.
"""
from pathlib import Path
import argparse
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from checks._common import DATA, OUTPUT, output_directory, protect_inputs, write_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DATA / "supplementary" / "thermal_spectra")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT / "thermal_quadrature")
    parser.add_argument("--recompute", action="store_true", help="Regenerate equilibrium spectra for orders 4/8/16 instead of reducing saved spectra")
    args = parser.parse_args()
    import numpy as np
    from qsh_environment_response import EnvironmentDevice, quadrature, thermal_response
    from threadpoolctl import threadpool_limits
    output = output_directory(args.output_dir)
    inputs = [] if args.recompute else [args.data_dir / f"clean_o{order}.npz" for order in (4, 8, 16)]
    history, previous = [], None
    with protect_inputs(inputs) as fingerprints, threadpool_limits(limits=1):
        device = EnvironmentDevice()
        for order in (4, 8, 16):
            if args.recompute:
                nodes, weights = quadrature(-.15, .79, .025, order)
                points = [device.spectral_record(float(energy)) for energy in nodes]
                spectrum = {key: np.asarray([point[key] for point in points])
                            for key in ("C", "rho", "D", "qlead", "qwindow")}
                spectrum.update(energies=nodes, weights=weights)
                np.savez_compressed(output / f"clean_o{order}.npz", **spectrum)
            else:
                with np.load(args.data_dir / f"clean_o{order}.npz", allow_pickle=False) as source:
                    spectrum = {name: source[name].copy() for name in source.files}
            rows = [[thermal_response(device, spectrum["energies"], spectrum["weights"], spectrum,
                                      float(mu), float(temperature))[0]
                     for mu in np.linspace(.25, .55, 13)] for temperature in (.005, .01, .02)]
            values = np.array([[row["kH"] for row in group] for group in rows])
            mass_error = max(abs(row["thermal_mass"]-1) for group in rows for row in group)
            difference = None if previous is None else float(np.max(abs(values-previous)))
            history.append({"order": order, "maximum_missing_thermal_mass": mass_error,
                            "maximum_response_difference_raw": difference,
                            "maximum_response_difference_paper": None if difference is None else difference/4,
                            "response_raw_kH": values, "response_paper_tilde_kappa2": -values/4})
            previous = values
    write_report(output / "report.json", {"status": "computed", "recomputed_spectra": args.recompute,
                 "input_sha256": fingerprints, "window": [-.15, .79], "panel_width": .025,
                 "energies": np.linspace(.25, .55, 13), "temperatures": [.005, .01, .02],
                 "conversion": "paper tilde_kappa2=-raw kH/4; absolute differences divide by4",
                 "scope": "Only quadrature order convergence on a fixed window; missing thermal mass is not a bound on coefficient error",
                 "history": history})
    print(output / "report.json")


if __name__ == "__main__":
    main()
