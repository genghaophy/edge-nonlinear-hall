# Data provenance and verification scope

## Three artifact layers

| Layer | Location | Purpose |
|---|---|---|
| Saved numerical inputs | `data/raw/` | Historical numerical arrays and accompanying records used by default calculation commands |
| Prepared panel data | `data/panels/figN/` | Display-ready arrays with matching definitions, parameters, units, and hashes |
| Reference artwork | `figures/reference/fig1.pdf` through `fig4.pdf` | Preserved figure PDFs from the source workspace |

New plots, calculations, and checks write under `outputs/` by default. They are not reference artifacts.

The prepared panel array values are preserved from the source results. All 13 default cache-extraction commands were run during packaging, so prepared NPZ files and their JSON records were written again; NPZ bytes can change even when the stored array values are identical. In the raw `figure3_ef05_extension/environment_data.npz`, only private path strings inside the serialized `metadata_json` field were sanitized, and the accompanying output checksum was updated. No physical numerical arrays were changed. The four reference PDFs were copied byte-for-byte.

Repository-relative paths replace private-workspace locations. Shared RM files and data formerly stored in a sibling project are now internal. No runtime access to the original workspace is required.

`provenance/source_map.json` records each original workspace-relative name, its release destination, byte count, original SHA-256, and release SHA-256. A hash difference for a Python, JSON, or NPZ file can reflect path adaptation, metadata sanitization, or cache-file regeneration. It must not be confused with a numerical recomputation. The source map is a packaging record, not a claim that every included backend freshly regenerates every cached result. Supplemental diagnostic sources and saved inputs have a separate record in `provenance/supplementary_source_map.json`.

## Reference PDFs

| File | SHA-256 |
|---|---|
| `figures/reference/fig1.pdf` | `e4df77e2462eda4e8dad9b1c35cde4277d2cba1833b208d751ecd65fe3b29ef5` |
| `figures/reference/fig2.pdf` | `797c9f7eef6b4d3cd6fe76bf66663be4938a4e430cc5c5df6517b1b5941b1bef` |
| `figures/reference/fig3.pdf` | `62ae1137eda4271e139e9b5c5ca6e4386092dbfeedce3f677af90a9a771a71e6` |
| `figures/reference/fig4.pdf` | `ba9812b7270309de49735e6c17c0f1792a21a75e1439ae4cbc74185ed0e260fc` |

The Fig. 2 reference includes the axis padding added to panel (b) and its inset to reveal hollow near-zero markers. The saved data, ticks, and panel geometry were retained in that revision. These hashes identify the supplied reference bytes; newly exported figures need not be byte-identical because metadata, fonts, or rendering libraries can differ.

## Panel contracts

Each panel has an NPZ plus a sibling JSON. The JSON records the figure and panel identifiers, array names/shapes/dtypes, NPZ SHA-256, calculation mode, parameters, observables, units, and source records. Plot loading checks the NPZ contract and rejects a mismatch. It does not run a transport solver and does not establish the numerical correctness of an arbitrary replacement dataset.

Source hashes in a panel JSON describe provenance. They are not all revalidated by the plotting loader. The complete packaging map and the numerical cache records serve different checks and should be examined separately.

## Historical source-version differences

Historical cache JSON may name an old generator version or a `historical_sources/` alias that is not included as a runtime file. Some generator hashes differ from current files, including files adapted for this release. The Fig. 4 extraction helpers verify the saved dataset checksum and retain these differences in `legacy_source_version_differences`; they do not silently relabel historical data as a fresh current-code run.

Therefore:

- Successful plotting demonstrates that the prepared data contract can be read and rendered.
- Successful default data extraction demonstrates that the included saved inputs can supply the expected panel arrays.
- A checksum match demonstrates artifact integrity, not a new numerical solution.
- A successful `--recompute` validates only the selected calculation and its returned checks. Agreement with a historical curve requires an explicit array comparison at stated tolerances.

The independent packaging process did not rerun all scientific scans. Historical diagnostic records and supplemental checks remain useful within their recorded parameter and convergence scope; they should not be described as an exhaustive fresh validation of this release.

## Comparing a new result

Save recomputed NPZ/JSON pairs under `outputs/recomputed/` as shown in [the workflow guide](workflow.md). Compare the observable definitions, units, parameter grids, seed lists, calculation mode, and source versions before comparing numeric arrays. For floating-point results, record maximum absolute/relative differences and an appropriate tolerance. Do not use matching plot appearance as a substitute for those checks.

When deliberately replacing a supplied dataset, retain its original record and update the panel JSON and release provenance consistently. A parameter change, a code correction, and a rendering change should each have their own documented reason.
