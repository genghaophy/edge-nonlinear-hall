# Edge nonlinear Hall transport

[中文说明](README.zh-CN.md)

Python code and saved numerical data for the four main figures of **Nonlinear Hall Effect from Inequivalent Edge Bands in Bulk Insulators**. The package contains 13 independent panel workflows, shared Rice–Mele (RM) and quantum spin Hall (QSH) transport backends, and the reference figure PDFs. All runtime inputs are inside this repository.

## Draw the figures from the supplied data

Run these commands from the repository root:

```sh
python -m pip install -r requirements-plot.txt
python code/run_panel_workflow.py --stage plot
```

This draws all 13 panels and assembles the four complete figures. It reads `data/panels/`, writes PDF/PNG/SVG files and metadata to `outputs/`, and does not load a transport solver. NumPy and Matplotlib are sufficient. The reference PDFs in `figures/reference/` are retained separately.

For one panel or one complete figure:

```sh
python code/run_panel_workflow.py --figure 3 --panel c --stage plot
python code/run_panel_workflow.py --figure 4 --stage assemble
```

Selecting a single panel does not assemble a complete figure. Missing or inconsistent panel data raise an error; plotting never starts a calculation automatically.

## Workflow stages

| Stage | Action |
|---|---|
| `plot` | Draw supplied panel data; also assemble each selected complete figure when `--panel` is omitted |
| `assemble` | Draw complete figures directly from all their panel datasets |
| `data` | Prepare panel datasets from the saved input caches; no plotting |
| `all` (default) | Prepare data, draw panels, and assemble complete figures |

```sh
python code/run_panel_workflow.py
```

The default prepares data from existing caches. It writes the selected NPZ/JSON pairs in `data/panels/`; it does not rerun transport scans. Only an explicit `--recompute` requests a new calculation. The recent workflow record is `outputs/last_run.json`.

## Recompute without replacing the supplied data

Install the computation dependencies using [the environment guide](docs/installation.md). A direct calculation command can save results in a separate location:

```sh
python code/panels/fig2/calc_a.py --recompute --output outputs/recomputed/fig2/a.npz
python code/panels/fig2/plot_a.py --data outputs/recomputed/fig2/a.npz --output outputs/recomputed/fig2/a
```

The first command writes both `a.npz` and `a.json`. The second reads that pair and writes the plot. The shared workflow's `--recompute` option uses the default `data/panels/` destination and replaces the selected supplied panel data, so use the direct commands above for comparisons. Transport, temperature, and disorder scans can take substantial time; a selected panel runs only its own scan. Fig. 1 regenerates a qualitative schematic rather than solving a transport model.

## Contents

```text
code/run_panel_workflow.py    Select figures, panels, and stages
code/panels/                 Independent calculation and plotting scripts
code/*.py                   Shared numerical models and readout backends
code/checks/                Supplemental numerical diagnostics
data/panels/                Display-ready NPZ arrays and JSON definitions
data/raw/                   Saved numerical inputs used to prepare panels
figures/reference/          Four unchanged reference PDFs
provenance/source_map.json  Source relocation map and SHA-256 records
scripts/                    Package verification tools
docs/                       Installation, workflow, conventions, and provenance
outputs/                    Generated results; ignored by Git
```

| Figure | Panels |
|---|---|
| 1 | Device schematic |
| 2 | RM bands; edge spectral density and LDOS; left-injected LPDOS; quadratic Hall response versus energy |
| 3 | QSH bands; temperature; nonmagnetic disorder; spin-conserving elastic probes |
| 4 | RM/QSH probe-coupling scans and weak-probe disorder response |

The manuscript sources and submission letter are not included. Panel commands and dataset details are in [the workflow guide](docs/workflow.md) and [the panel directory guide](code/panels/README.md).

## Conventions and provenance

The physical drive is `V1 = +V/2`, `V2 = -V/2`, with the full source–drain voltage `V`, and the measured transverse voltage is `V3 - V4`. The dimensionless quadratic coefficient is `t2*kappa2/e` for RM and `t*kappa2/e` for QSH. Historical raw `kH` values are converted by `-kH/4` once, on the calculation side; display-ready data are not converted again. Fig. 3(c) shows SEM over eight disorder realizations, while Fig. 4(c,d) show sample SD over sixteen realizations. See [the conventions guide](docs/conventions.md).

The reference PDFs and saved numerical arrays are preserved from the source workspace. Paths and wrappers were adapted for this independent package, with original and release hashes recorded in `provenance/source_map.json`. Some historical generator versions differ from the included source versions; those differences remain in the JSON records. Verified cache extraction reproduces saved inputs and is not a fresh validation of the current solver. See [the provenance guide](docs/provenance.md) before comparing a new calculation with an archived result.

The tested source runtime is **Python 3.13.5 (Anaconda)** with NumPy 2.1.3, SciPy 1.15.3, Matplotlib 3.10.0, Kwant 1.5.0, tinyarray 1.2.5, and threadpoolctl 3.5.0. Dependency files record these package versions. Fonts and rendering libraries can affect the appearance and bytes of newly exported graphics.

## Verify the package

```sh
python scripts/verify_release.py
python scripts/verify_release.py --relocate --plot
```

The first checks packaged data and release records. The second also copies the repository to an independent temporary directory, prepares all 13 datasets, draws all 13 panels, and assembles four figures with solver imports prohibited. Reports are written under `outputs/validation/`. An optional `--compute` check evaluates one clean weak-probe point per model and requires the computation dependencies. See [verification scope](docs/verification.md); these checks do not rerun full scientific scans.

## License and citation

The original code is available under the [MIT License](LICENSE), copyright 2026 H. Geng and contributors. Dependencies are installed separately and retain their own licenses; none are vendored here.

Use [CITATION.cff](CITATION.cff) for software attribution. It records the four authors using their supplied initials and the [public repository](https://github.com/genghaophy/edge-nonlinear-hall). No software DOI or article publication identifier has been assigned.
