# Supplemental numerical diagnostics

These seven entry points expose selected supplemental checks independently of the main-figure workflow. Run them from the repository root. Their defaults write new reports under `outputs/checks/` and retain the supplied scientific inputs. They do not refresh the main-figure datasets in `data/panels/`.

Unlike the main-figure default cache extraction, most diagnostics perform fresh numerical work when invoked. All except `thermal_window.py` require the computation dependencies, including Kwant; the cached-spectrum reduction in `thermal_quadrature.py` also imports the QSH model backend. `thermal_window.py` needs only NumPy and the supplied spectrum.

## Entry points and default outputs

| Entry point | Scope | Default output directory | Principal files |
|---|---|---|---|
| `rm_reference.py` | One clean RM reference at `E_F/t2=0.10`, `tau/t2=0.15`; floating currents, reciprocity, gauge shift, derivative and finite-bias checks | `outputs/checks/rm_reference/` | `report.json` |
| `qsh_reference.py` | Clean, one disorder realization, elastic probes, and an optional independent small-device thermal check | `outputs/checks/qsh_reference/` | `report.json` |
| `rm_spatial.py` | One-energy RM injectivity/LDOS reconstruction and local-neutrality potential diagnosis | `outputs/checks/rm_spatial/` | `spatial_diagnostic.npz`, `report.json` |
| `rm_weak_probe.py` | RM seven-coupling audit at the saved reference energy; full/frozen responses and feedback decomposition | `outputs/checks/rm_weak_probe/` | `rm_tc_scan.npz`, `rm_tc_scan.json` |
| `rm_bare_injectivity.py` | Bare two-terminal RM scattering-state injection derivative and comparison with the selected coupling report | `outputs/checks/rm_bare_injectivity/` | `bare_injectivity_response.json` |
| `thermal_quadrature.py` | QSH thermal reduction at orders 4, 8, and 16 on the fixed production energy window | `outputs/checks/thermal_quadrature/` | `report.json`; fresh `clean_o4.npz`, `clean_o8.npz`, `clean_o16.npz` only with `--recompute` |
| `thermal_window.py` | Integration-by-parts boundary diagnostic using the saved order-16 thermal spectrum | `outputs/checks/thermal_window/` | `report.json` |

Each entry supports `--output-dir`. Use a directory under `outputs/` for new results. Output directories inside `data/raw/` are rejected. Repeating a command in the same output directory replaces that diagnostic's previous files.

## RM checks

```sh
python code/checks/rm_reference.py
python code/checks/rm_spatial.py
python code/checks/rm_weak_probe.py
python code/checks/rm_bare_injectivity.py
```

`rm_spatial.py` and `rm_weak_probe.py` accept `--input`; their default is `data/raw/rm_figure2/figure2_data.npz`. The spatial diagnostic uses the saved eight-cell RM geometry, rather than the wider 30-cell Fig. 2(c) map. It separates local spectral density, source-injectivity imbalance, and the local-neutrality potential. It does not solve Poisson/Hartree electrostatics.

The weak-probe audit computes seven couplings, `tau/t2 = 0.20, 0.15, 0.10, 0.05, 0.02, 0.01, 0.005`, at the reference energy. It retains source-electrode and internal-potential feedback, checks a frozen-transmission prediction, and refines derivatives at selected couplings. This is not the full energy scan.

The bare-injectivity check defaults to the historical report `data/raw/supplementary/weak_probe_audit/rm_tc_scan.json`. To compare with a freshly generated seven-coupling audit:

```sh
python code/checks/rm_weak_probe.py
python code/checks/rm_bare_injectivity.py --coupling-report outputs/checks/rm_weak_probe/rm_tc_scan.json
```

The bare check has no side probes. It computes equilibrium injectivities, the physical-voltage derivative of their left/right imbalance, and the weak-probe prediction at two derivative steps. It records differences against each finite coupling rather than assuming that weak coupling removes electrostatic or electrode feedback.

## QSH selected-point check

```sh
python code/checks/qsh_reference.py
```

The zero-temperature checks use `E_F/t=0.40` on the established 31 × 12 device: a clean sample, one scalar-disorder sample with `W_dis/t=0.50` and seed 1700, and a spin-conserving elastic-probe case with `gamma_phi/t=0.10`. Checks include kernel differentiation, derivative refinement, direct finite-bias floating currents, and independent equilibrium-matrix comparisons. The clean coefficient is compared with `data/raw/supplementary/floating_probe_reference.json`.

The default also runs an independent finite-temperature direct-current test on a **9 × 6** device at `E_F/t=0.40` and `kBT/t=0.02`, using its own fixed integration window. This tests the selected finite-temperature reduction; it does not establish size convergence of the 31 × 12 device. To omit this more expensive part:

```sh
python code/checks/qsh_reference.py --skip-thermal
```

`--skip-thermal` leaves the clean, disorder, and elastic-probe checks active.

## Thermal quadrature and boundary checks

```sh
python code/checks/thermal_quadrature.py
python code/checks/thermal_window.py
```

By default, `thermal_quadrature.py` reduces the saved `clean_o4.npz`, `clean_o8.npz`, and `clean_o16.npz` spectra from `data/raw/supplementary/thermal_spectra/`. It evaluates 13 Fermi energies from `E_F/t=0.25` to `0.55` at `kBT/t=0.005, 0.01, 0.02`, records changes between quadrature orders, and reports missing thermal mass. Its `--data-dir` selects an alternate spectrum directory.

To regenerate the equilibrium spectra before reduction:

```sh
python code/checks/thermal_quadrature.py --recompute --output-dir outputs/checks/thermal_quadrature_recomputed
```

This performs numerous spectral transport solves at all three orders and stores new spectra in the selected output directory. It uses the fixed window `[-0.15, 0.79]t` and panel width `0.025t`. It tests quadrature order on that window, rather than window enlargement or device-size convergence, and does not replace the independent finite-temperature current test. The report records differences; it does not declare a universal convergence bound.

`thermal_window.py` accepts `--input`, defaulting to the saved order-16 spectrum. It compares the occupation and energy-derivative source forms over the same 13 energies and three temperatures. The reported boundary difference is not an estimate of all omitted energies, and missing thermal mass alone is not a bound on the quadratic-coefficient error.

## Units and interpretation

The physical convention is `e>0`, `mu=E_F-eV`, `H=H0-eU`, source voltages `+V/2,-V/2`, and `V_perp=V3-V4`, with `V` the full source–drain voltage. Historical raw `kH` fields use a single-terminal source-energy amplitude. The plotted dimensionless coefficient is obtained once:

```text
tilde_kappa2 = -kH_raw/4
RM:  tilde_kappa2 = t2*kappa2/e
QSH: tilde_kappa2 = t*kappa2/e
```

Absolute raw coefficient errors divide by four; a signed raw coefficient difference acquires the minus sign as well. `rm_reference.py` includes separate raw and paper fields. `qsh_reference.py` retains historical raw `kH` fields and states their conversion. The weak-probe and thermal reports explicitly distinguish raw and paper quantities. Do not convert fields that already carry paper units a second time.

For the RM spatial diagnostic, `paper_u_V = reconstructed_ud/2` is the physical characteristic potential `U/V` for the full bias; the historical raw drive field is also retained. Its edge localization should not be identified with an induced charge profile.

These are checks of selected coherent models with first-order local charge neutrality, not full nonlinear electrostatics or complete disorder/size convergence studies. Thermal direct-current checks and fresh equilibrium-spectrum generation can be costly; no wall-time estimate applies to every machine. The numerical wrappers limit native thread pools to one thread where specified in their algorithms.

The source relocation and wrapper scope are recorded in [the supplemental source map](../../provenance/supplementary_source_map.json). See [conventions](../../docs/conventions.md), [provenance](../../docs/provenance.md), and [package verification](../../docs/verification.md) for the distinction between saved-data integrity and fresh numerical validation.
