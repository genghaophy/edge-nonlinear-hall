# Calculation and plotting workflow

Commands below assume the repository root as the working directory. The common driver selects one of four figures and its panels:

```sh
python code/run_panel_workflow.py --figure all --stage plot
```

## Driver options

| Option | Values | Meaning |
|---|---|---|
| `--figure` | `all`, `1`, `2`, `3`, `4` | Default: all figures |
| `--panel` | `a`, `b`, `c`, `d` | One panel; requires a single figure; Fig. 1 only has a |
| `--stage` | `data`, `plot`, `assemble`, `all` | Default: all |
| `--recompute` | flag | Use a new calculation in the data/all stages |

The driver writes the most recent operation log to `outputs/last_run.json`, including the selected stages, child commands, and captured output. A single-panel selection suppresses whole-figure assembly. `--stage assemble` requires a whole figure rather than a panel selection.

## Prepare data from the saved inputs

```sh
python code/run_panel_workflow.py --figure 3 --panel c --stage data
```

Without `--recompute`, the calculation script extracts the required arrays from `data/raw/`, checks the panel's expected parameters, performs any required historical convention conversion, and writes the panel NPZ/JSON pair. Fig. 1 generates saved schematic coordinates. This step writes `data/panels/`; it does not launch the transport scans.

The default all-stage command does the same preparation for each selected panel, draws each one, then assembles complete figures:

```sh
python code/run_panel_workflow.py --stage all
```

## Plot and assemble

```sh
python code/run_panel_workflow.py --figure 2 --panel b --stage plot
python code/run_panel_workflow.py --figure 2 --stage assemble
```

Plotting loads saved data only. `load_panel` checks the figure/panel identifiers, NPZ SHA-256, array names, shapes, and dtypes against the JSON. Missing data or mismatched definitions cause an error. The arrays are loaded with `allow_pickle=False`.

| Output | Default location |
|---|---|
| Single panel | `outputs/panels/figN/<panel>.pdf`, `.png`, `.svg` |
| Single-panel metadata | `outputs/panels/figN/<panel>_metadata.json` |
| Fig. 1 assembly | `outputs/assembled/fig1_device.*` |
| Fig. 2 assembly | `outputs/assembled/fig2_revised.*` |
| Fig. 3 assembly | `outputs/assembled/fig3_revised.*` |
| Fig. 4 assembly | `outputs/assembled/fig_probe_rm_qsh.*` |

Assemblies read NPZ/JSON pairs rather than embedding the individual rendered panels. They retain the reference panel layout. The exporter checks visible text against the canvas and checks cross-row text overlap in assemblies. The four reference PDFs are kept in `figures/reference/`.

## Keep new calculations separate

The shared driver does not accept an alternate data output root. Its `--recompute` flag replaces selected default panel datasets. For a comparison, run the individual calculation and plot scripts with explicit paths:

```sh
python code/panels/fig2/calc_a.py --recompute --output outputs/recomputed/fig2/a.npz
python code/panels/fig2/plot_a.py --data outputs/recomputed/fig2/a.npz --output outputs/recomputed/fig2/a
```

Each calculation writes metadata beside its NPZ, including the parameters, source hashes, mode (`reuse` or `recompute`), and array contract. Each plot records input and output hashes. Transport scans may be expensive; an individual panel does not start the other panels' scans.

To assemble a complete figure from alternate datasets:

```sh
python code/panels/assemble.py --figure 2 --data-root outputs/recomputed --output-dir outputs/recomputed/assembled
```

Provide all four Fig. 2 pairs under `outputs/recomputed/fig2/`. The alternate root is used for every panel of the selected figure; there is no fallback to the bundled data.

## Parameter changes

The panel scripts and shared model backends contain the established parameters and explicit expected grids. The driver selects an operation; it is not a general parameter-sweep interface. Changing a scientific parameter may require updating the calculation, data validation, metadata, and plotting range together. Save exploratory outputs outside `data/panels/` and record the change before comparing with the reference figures.

The supplemental diagnostics in `code/checks/` are separate from the main-figure driver. They do not run during plotting, assembly, or default cache extraction. Their scope is described in [that directory's README](../code/checks/README.md).
