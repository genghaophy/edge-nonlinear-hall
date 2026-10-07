# Independent main-figure panels

Each panel has `calc_<panel>.py` to prepare numeric data and `plot_<panel>.py` to draw it. Calculation scripts do not draw figures. Plot scripts load a matching NPZ/JSON pair and do not import the transport solvers. `assemble.py` draws the complete figures directly from saved panel data at the reference layout.

Run commands from the repository root. The common entry point is:

```sh
python code/run_panel_workflow.py --stage plot
```

| Directory | Panel | Content |
|---|---|---|
| `fig1` | a | Qualitative device schematic |
| `fig2` | a | RM ribbon bands and edge localization |
| `fig2` | b | Fixed-momentum spectral density; integrated edge LDOS |
| `fig2` | c | Left-injected LPDOS on 30 × 31 RM cells |
| `fig2` | d | Quadratic Hall coefficient for four values of Delta |
| `fig3` | a | QSH ribbon bands and selected Fermi energies |
| `fig3` | b | Fixed-Fermi-energy temperature scan |
| `fig3` | c | Nonmagnetic disorder scan, mean and SEM |
| `fig3` | d | Spin-conserving elastic-probe scan |
| `fig4` | a | RM symmetric probe-coupling scan and two-terminal prediction |
| `fig4` | b | QSH symmetric probe-coupling scan and two-terminal prediction |
| `fig4` | c | RM weak-probe disorder response, mean and sample SD |
| `fig4` | d | QSH weak-probe disorder response, mean and sample SD |

## Panel interface

Default calculation reads the included saved inputs and writes `data/panels/figN/<panel>.npz` plus the matching JSON. `--recompute` is an explicit new calculation. Use `--output` to preserve the supplied data:

```sh
python code/panels/fig2/calc_a.py --recompute --output outputs/recomputed/fig2/a.npz
python code/panels/fig2/plot_a.py --data outputs/recomputed/fig2/a.npz --output outputs/recomputed/fig2/a
```

The plot `--output` is a stem without an extension; PDF, PNG, SVG, and `_metadata.json` are written. The calculation `--output` must end in `.npz`. Keep its sibling JSON beside it.

For an alternate complete-figure dataset, use the same `figN/<panel>.npz` structure under another root and provide every panel of that figure:

```sh
python code/panels/assemble.py --figure 2 --data-root outputs/recomputed --output-dir outputs/recomputed/assembled
```

The example above needs all four Fig. 2 NPZ/JSON pairs. It does not reuse missing panels from `data/panels/`.

## Shared helpers

- `common.py`: dataset contracts, SHA-256, command-line options, plotting style, and exports.
- `assemble.py`: fixed figure layouts and drawing from saved data.
- `fig2/_data.py`, `fig2/_ribbon.py`: RM inputs and ribbon calculations.
- `fig3/_calc.py`, `fig4/_calc.py`: calculation-side helpers with lazy transport imports.
- `_style.py`: drawing styles; these do not run calculations.

All numerical backends are inside `code/`, and raw inputs are inside `data/raw/`. Historical source aliases in JSON are provenance records rather than external runtime dependencies. See [workflow details](../../docs/workflow.md), [conventions](../../docs/conventions.md), and [source-version limits](../../docs/provenance.md).
