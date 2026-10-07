# Package verification

The verifier has separate artifact, relocation, and representative-computation checks. Run it from the repository root. Generated reports and logs belong in `outputs/validation/`.

## Artifact integrity

```sh
python scripts/verify_release.py
```

This validates all 13 prepared NPZ/JSON contracts, checks numerical arrays for nonfinite values with an explicit allowance for QSH band-plot NaN masks, and checks embedded NPZ strings for private absolute paths. When `provenance/manifest.json` is present, it also checks each listed release-file SHA-256. It does not evaluate transport responses or draw figures.

## Independent-directory operation

```sh
python scripts/verify_release.py --relocate --plot
```

This copies the package to a temporary directory, excluding Git history, previous outputs, and bytecode. It runs all 13 default data preparations and checks that their arrays match the supplied panel data exactly, including intentional NaN masks. It then runs all 13 panel plots and four complete-figure assemblies.

A subprocess import guard rejects Kwant, tinyarray, and selected transport-backend imports during these stages. This tests that cache preparation and drawing work independently of the original workspace and without the solver runtime. The temporary copy is cleaned up; logs and the generated figures are retained in `outputs/validation/`. To select a writable temporary parent:

```sh
python scripts/verify_release.py --relocate --plot --temp-root outputs/temporary
```

This verifies operation and data integrity. It does not prove that a fresh solver reproduces every historical scan or that figure exports are byte-identical to the reference PDFs.

## Representative numerical computation

```sh
python scripts/verify_release.py --compute
```

This requires the computation dependencies and evaluates one clean point per model:

| Model | Fermi energy | Equal probe coupling |
|---|---|---|
| RM | `E_F/t2 = 0.10` | `tau/t2 = 0.003` |
| QSH | `E_F/t = 0.40` | `tau/t = 0.003` |

The check compares the full four-terminal coefficient with the packaged clean disorder-reference coefficient at an absolute tolerance of `1e-7`, checks returned transport diagnostics, and records an independently computed bare two-terminal weak-probe prediction. The finite-coupling coefficient and limiting prediction are distinct observables; the report does not require them to be exactly equal. No full coupling, thermal, energy, or disorder scan is launched.

Options can be combined. The most recent report is `outputs/validation/last_validation.json` and is overwritten by the next verification run.

## Supplemental diagnostics

`code/checks/` contains the separate supplemental calculation workflows. Consult [its README](../code/checks/README.md) for the scope, runtime, and output of each diagnostic. They are not started by the main-figure driver or the default package verifier.

Historical validation records, fresh representative checks, and cache-integrity checks address different questions. Record the command, runtime, parameters, and result when reporting a new validation. See [provenance](provenance.md) for the source-version limits of cached data.
