"""Fixed-window thermal boundary diagnostic with explicit raw and paper units.

The saved production spectra are never modified. The integration algorithm
is the original diagnostic; the full source-drain voltage conversion is
applied once: paper difference = -raw kH difference / 4.
"""
from pathlib import Path
import numpy as np
from checks._common import protect_inputs, write_report


def audit(source, output):
    source = Path(source)
    with protect_inputs([source]) as fingerprints:
        with np.load(source, allow_pickle=False) as data:
            e, weights, matrices = data["energies"], data["weights"], data["C"]
            qwindow = data["qwindow"]
        alpha = np.array([1., -1., 0., 0.])
        cases = []
        for mu in np.linspace(.25, .55, 13):
            for temperature in (.005, .01, .02):
                arg = (e-mu)/(2*temperature)
                w = 1/(4*temperature*np.cosh(arg)**2)
                wp = -w*np.tanh(arg)/temperature
                cbar = np.einsum("e,eab->ab", weights*w, matrices)
                q_occupation = -.5*np.einsum("e,eab,b->a", weights*wp,
                                             matrices, alpha**2)
                q_derivative = np.einsum("e,ea->a", weights*w, qwindow)
                delta_beta_code = -np.linalg.solve(
                    cbar[2:, 2:], (q_occupation-q_derivative)[2:])
                raw_difference = float(delta_beta_code[0]-delta_beta_code[1])
                cases.append({"EF_over_t": float(mu), "kBT_over_t": temperature,
                              "source_form_difference_max_raw_code": float(np.max(abs(q_occupation-q_derivative))),
                              "occupation_minus_energy_derivative_kH_raw": raw_difference,
                              "occupation_minus_energy_derivative_tilde_kappa2_paper": -raw_difference/4})
        report = {"input_sha256": fingerprints, "finite_window_over_t": [-.15, .79],
                  "quadrature_order": 16, "cases": cases,
                  "conversion": "tilde_kappa2_paper=-kH_raw/4 exactly once; raw source amplitude is half the full paper bias with opposite onsite sign",
                  "meaning": "Integration-by-parts boundary diagnostic on a fixed window; not all omitted-energy error or window-enlargement convergence",
                  "maximum_abs_paper_tilde_kappa2_difference": max(abs(item["occupation_minus_energy_derivative_tilde_kappa2_paper"]) for item in cases),
                  "status": "computed"}
    write_report(output, report)
    return report
