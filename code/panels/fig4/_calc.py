"""Figure 4 panel data, with read-only cache reuse and lazy transport imports.

The legacy datasets already use the paper's full source--drain voltage.
Only an explicit recompute loads Kwant through the existing model APIs.
No legacy calculator main() or legacy output writer is called here.
"""
from importlib import import_module
from pathlib import Path
import hashlib
import json

import numpy as np

from panels.common import source_record

PROJECT = Path(__file__).resolve().parents[3]
PROTOCOL = "mu=EF-eV,H=H0-eU,V1=+VSD/2,V2=-VSD/2,Vperp=V3-V4"
PARAMETERS = {
    "rm": dict(tx=.25, t1=.3, t2=1., delta=.2, cells=8, length=31),
    "qsh": dict(m=-1., delta=.2, t=1., A=1., length=31, width=12),
}
ENERGIES = {"rm": np.array([.05, .10, .15]),
            "qsh": np.array([.3, .4, .5])}
STRENGTHS = {"rm": np.array([0., .025, .05, .10]),
             "qsh": np.array([0., .10, .25, .50])}
TAU = np.geomspace(.005, 1., 31)
SEEDS = np.arange(1700, 1716)
STEP = 1e-5


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _validated_cache(path, sidecar):
    """Verify the historical artifact and report source-version drift.

    An immutable, checksum-verified formal-figure dataset remains usable for
    reproduction when its old generator has subsequently changed. This is
    not evidence that the current backend reproduces that historical run.
    """
    metadata = json.loads(sidecar.read_text(encoding="utf-8"))
    if metadata.get("status", "complete") not in ("complete", "completed"):
        raise ValueError(f"Calculation is not complete: {sidecar}")
    expected = metadata.get("arrays_sha256", metadata.get("output_sha256"))
    if expected is None or _digest(path) != expected:
        raise ValueError(f"Validated dataset checksum differs: {path}")
    sources = metadata.get("source_files", metadata.get("source_sha256", {}))
    differences = []
    for filename, checksum in sources.items():
        current_path = PROJECT / filename
        current = _digest(current_path) if current_path.is_file() else None
        if current != checksum:
            differences.append(dict(path=filename, recorded_sha256=checksum,
                                    current_sha256=current))
    metadata["_current_source_differences"] = differences
    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: archive[name].copy() for name in archive.files}
    return arrays, metadata, [source_record(path), source_record(sidecar)]


def _numeric_checks(record, maxima):
    """Check returned diagnostics only; this does not launch extra calculations."""
    for key, value in record.get("checks", {}).items():
        if key.endswith("error") and np.ndim(value) == 0:
            number = float(value)
            if not np.isfinite(number) or abs(number) > 2e-7:
                raise ValueError(f"Transport diagnostic failed: {key}={number}")
            maxima[key] = max(maxima.get(key, 0.), abs(number))
    raw = record.get("raw", {})
    for key in ("probe_linear_error", "probe_quadratic_error",
                "current_sum_error", "characteristic_sum_error"):
        if key in raw:
            number = abs(float(raw[key]))
            if not np.isfinite(number) or number > 2e-7:
                raise ValueError(f"Floating-current diagnostic failed: {key}={number}")
            maxima[key] = max(maxima.get(key, 0.), number)
    if "contactformula_error" in record:
        error = float(record["contactformula_error"])
        if not np.isfinite(error) or error > 2e-6:
            raise ValueError(f"Contact readout formula failed: {error}")
        maxima["contactformula_error"] = max(maxima.get("contactformula_error", 0.), error)
    if "kH" in raw and not np.isclose(record["kappa_full"], -raw["kH"]/4,
                                     rtol=0, atol=1e-12):
        raise ValueError("Historical full-bias conversion differs from -kH/4.")


def _metadata(model, recompute, sources, checks, *, disorder=False):
    unit = "t2" if model == "rm" else "t"
    return dict(
        sources=sources, calculation_helper=source_record(__file__),
        parameters=dict(PARAMETERS[model]), model=model,
        observable=("kappa2(W_dis)/kappa2(0)" if disorder else "dimensionless kappa2"),
        units={"energy": unit, "probe_hopping": unit,
               "ordinate": "dimensionless", **({"disorder": unit} if disorder else {})},
        conversion=("API converts historical kH to paper kappa=-kH/4 once"
                    if recompute else "Already in paper units; copied without a sign or bias conversion"),
        physical_convention=PROTOCOL,
        electrostatics="First-order local charge neutrality; coherent, T=0, B=0",
        probe_constraint="Identical floating probes, t3=t4=tau, I3=I4=0",
        calculation="recomputed_this_panel" if recompute else "validated_cache_reuse",
        validation=checks,
        statistics=({"spread": "sample standard deviation, ddof=1; not SEM",
                     "nonzero_strength_samples": 16, "clean_independent_samples": 1,
                     "normalization": "signed clean four-terminal kappa at the same tau"}
                    if disorder else {"spread": "none", "bare_prediction": "independent, not fitted"}),
    )


def _validate_coupling(data, model):
    if not np.array_equal(data["energies"], ENERGIES[model]):
        raise ValueError("Unexpected Fermi-energy series.")
    if not np.array_equal(data["tau"], TAU):
        raise ValueError("Expected the validated 31-point probe-hopping scan.")
    if data["equal_kappa_full"].shape != (3, 31) or data["bare_kappa_predicted"].shape != (3,):
        raise ValueError("Coupling-scan array shapes differ.")
    if not all(np.isfinite(value).all() for value in data.values()):
        raise ValueError("Nonfinite coupling-panel data.")


def coupling_data(model, recompute=False):
    """Calculate just the RM or QSH coupling panel, or copy its validated cache."""
    names = ("energies", "tau", "equal_kappa_full", "bare_kappa_predicted")
    if not recompute:
        folder = "weak_probe_comparison" if model == "rm" else "qsh_weak_probe_comparison"
        path = PROJECT / "data/raw" / folder / "coupling_scan.npz"
        cached, original, sources = _validated_cache(path, path.with_suffix(".json"))
        if original["parameters"] != PARAMETERS[model]:
            raise ValueError("Cached model parameters differ from the formal panel.")
        data = {name: cached[name] for name in names}
        checks = {"dataset_checksum_verified": True,
                  "legacy_source_version_differences": original["_current_source_differences"],
                  "calculation_validation": original.get("validation", {})}
    else:
        # These imports must stay in this explicit branch. Plot modules never import them.
        from threadpoolctl import threadpool_limits
        api = import_module("edge_probe_readout" if model == "rm" else "qsh_edge_probe_readout")
        parameters = (api.RMParameters() if model == "rm" else api.QSHParameters())
        bare, rows, maxima = [], [], {}
        with threadpool_limits(limits=1):
            for energy in ENERGIES[model]:
                record = api.bare_edge_quantities(parameters, float(energy),
                                                  energy_step=STEP, voltage_step=STEP)
                _numeric_checks(record, maxima)
                bare.append(float(record["kappa_predicted"]))
                row = []
                for hopping in TAU:
                    if model == "rm":
                        record = api.four_terminal_readout(parameters, float(energy),
                                                          float(hopping), float(hopping), step=STEP)
                    else:
                        record = api.four_terminal_readout(parameters, float(energy),
                                                          float(hopping), step=STEP)
                    _numeric_checks(record, maxima)
                    row.append(float(record["kappa_full"]))
                rows.append(row)
                print(f"Figure 4 {model} coupling: completed E_F={energy:g}", flush=True)
        data = dict(energies=ENERGIES[model].copy(), tau=TAU.copy(),
                    equal_kappa_full=np.asarray(rows), bare_kappa_predicted=np.asarray(bare))
        sources = [source_record(api.__file__), source_record(__file__),
                   source_record(PROJECT / "code/compute_floating_study.py"),
                   source_record(PROJECT / "code/compute_initial_study.py")]
        if model == "qsh":
            sources.append(source_record(PROJECT / "code/edge_probe_readout.py"))
        checks = {"returned_diagnostic_maxima": maxima,
                  "additional_legacy_whole_figure_checks": "not rerun by this panel calculator"}
    _validate_coupling(data, model)
    return data, _metadata(model, recompute, sources, checks)


def _validate_disorder(data, model):
    if not np.array_equal(data["strengths"], STRENGTHS[model]):
        raise ValueError("Disorder strengths differ from this model's formal panel.")
    if data["four_kappa"].shape != (4, 16):
        raise ValueError("Expected four disorder strengths and sixteen stored seed slots.")
    if not np.array_equal(data["sample_counts"], [1, 16, 16, 16]):
        raise ValueError("Clean/nonzero independent sample counts differ.")
    if not np.array_equal(data["seeds"], SEEDS) or not np.isclose(data["tau"], .003, rtol=0, atol=1e-14):
        raise ValueError("Seed ordering or weak probe hopping differs.")
    if not np.isclose(data["energy"], .10 if model == "rm" else .40, rtol=0, atol=1e-14):
        raise ValueError("Disorder panel Fermi energy differs.")
    if not all(np.isfinite(value).all() for value in data.values()):
        raise ValueError("Nonfinite disorder-panel data.")
    if float(data["clean_kappa"]) == 0:
        raise ValueError("Cannot normalize by a zero clean coefficient.")
    normalized = data["four_kappa"] / data["clean_kappa"]
    if not np.allclose(normalized, data["four_normalized"], rtol=1e-12, atol=1e-12):
        raise ValueError("Signed clean normalization differs.")
    if not np.allclose(normalized.mean(axis=-1), data["four_normalized_mean"], rtol=1e-12, atol=1e-12):
        raise ValueError("Stored ensemble mean differs.")
    if not np.allclose(normalized.std(axis=-1, ddof=1), data["four_normalized_std"], rtol=1e-12, atol=1e-12):
        raise ValueError("Stored sample SD differs.")
    if not np.allclose(normalized[0], 1., rtol=0, atol=1e-12) or data["four_normalized_std"][0] != 0:
        raise ValueError("Clean reference must be one with no ensemble spread.")


def disorder_data(model, recompute=False):
    """Calculate only this model's four-terminal disorder observable."""
    if not recompute:
        path = PROJECT / "data/raw/fig4_disorder_asymmetric_ranges/ensemble.npz"
        cached, original, sources = _validated_cache(path, path.with_name("validation.json"))
        if original["parameters"][0 if model == "rm" else 1] != PARAMETERS[model]:
            raise ValueError("Cached disorder model parameters differ.")
        matches = np.flatnonzero(cached["models"] == model)
        if len(matches) != 1:
            raise ValueError("Expected exactly one selected model in the saved ensemble.")
        index = int(matches[0])
        data = {name: cached[name][index].copy() for name in
                ("strengths", "sample_counts", "four_kappa", "four_normalized",
                 "four_normalized_mean", "four_normalized_std")}
        data.update(energy=cached["energies"][index].copy(),
                    clean_kappa=cached["clean_kappa"][index].copy(),
                    tau=cached["tau"].copy(), seeds=cached["seeds"].copy())
        checks = {"dataset_checksum_verified": True,
                  "legacy_source_version_differences": original["_current_source_differences"],
                  "model_index": index, "original_status": original["status"]}
    else:
        from threadpoolctl import threadpool_limits
        api = import_module("disorder_probe_readout")
        parameters = api.RMParameters() if model == "rm" else api.QSHParameters(width=12)
        energy = .10 if model == "rm" else .40
        values, maxima = [], {}
        with threadpool_limits(limits=1):
            for strength in STRENGTHS[model]:
                row = []
                for seed in ([0] if strength == 0 else SEEDS):
                    grid = api.make_disorder(model, parameters, float(strength), int(seed))
                    record = api.four_quantities(model, parameters, energy, grid, .003, STEP)
                    _numeric_checks(record, maxima)
                    row.append(float(record["kappa_full"]))
                if strength == 0:
                    row *= len(SEEDS)  # one clean calculation; replicated slots carry no extra samples
                values.append(row)
                print(f"Figure 4 {model} disorder: completed W={strength:g}", flush=True)
        full = np.asarray(values)
        clean = np.asarray(full[0, 0])
        if float(clean) == 0:
            raise ValueError("Cannot normalize by a zero clean coefficient.")
        normalized = full / clean
        spread = normalized.std(axis=-1, ddof=1)
        spread[0] = 0.
        data = dict(strengths=STRENGTHS[model].copy(), energy=np.array(energy),
                    tau=np.array(.003), seeds=SEEDS.copy(), sample_counts=np.array([1, 16, 16, 16]),
                    four_kappa=full, clean_kappa=clean, four_normalized=normalized,
                    four_normalized_mean=normalized.mean(axis=-1), four_normalized_std=spread)
        sources = [source_record(filename) for filename in api.source_contract()]
        sources.extend([source_record(__file__),
                        source_record(PROJECT / "code/compute_initial_study.py")])
        checks = {"returned_diagnostic_maxima": maxima,
                  "additional_bare_reference_calculations": "not needed for this plotted observable"}
    _validate_disorder(data, model)
    metadata = _metadata(model, recompute, sources, checks, disorder=True)
    metadata["parameters"].update(energy=float(data["energy"]), tau=float(data["tau"]),
                                  strengths=data["strengths"].tolist(), seeds=SEEDS.tolist())
    return data, metadata
