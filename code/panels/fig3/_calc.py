"""Shared calculation-side provenance and the existing QSH numerical backend.

No Matplotlib imports. Historical archives are read-only. The numerical
Hamiltonian, contact and screening implementation remains in the validated
qsh_environment_response backend; individual panel scripts select their own
scan rather than dispatching the legacy whole-figure driver.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from panels.common import source_record

PROJECT = Path(__file__).resolve().parents[3]
CODE = PROJECT / "code"
SOURCE = PROJECT / "data/raw/figure3_environment/environment_data.npz"
EXTENSION = PROJECT / "data/raw/figure3_ef05_extension/environment_data.npz"
ENERGIES = np.array([.3, .4, .5])
STRENGTHS = np.array([0., .2, .5, 1.])
GAMMA = np.array([0., .01, .03, .1, .3, 1.])
SEEDS = np.arange(1700, 1708)
STEP = 1e-5


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def provenance(panel, recompute, sources, **extra):
    return dict(figure="fig3", panel=panel,
        mode="independent panel calculation" if recompute else "verified archive extraction",
        sources=[source_record(path) for path in sources],
        source_sha256={source_record(path)['path']: digest(path) for path in sources},
        response_convention="dimensionless tilde kappa2=t*kappa2/e=-historical kH/4",
        bias="V1=+V/2, V2=-V/2; identical floating probes 3(top),4(bottom)",
        time_reversal=True, magnetic_field=0., **extra)


def environment_archives():
    original, extension = load(SOURCE), load(EXTENSION)
    validation_path = EXTENSION.with_name("validation.json")
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if (validation["source_data_sha256"] != digest(SOURCE)
            or validation["output_sha256"] != digest(EXTENSION)
            or not validation["original_production_data_and_modules_unchanged"]):
        raise ValueError("EF=.5 extension must match the checked production archive")
    if not np.array_equal(original["selected_energies"], [.3, .4]):
        raise ValueError("Expected original EF=.3,.4 rows")
    if not np.array_equal(extension["selected_energies"], [.5]):
        raise ValueError("Expected extension EF=.5 row")
    for key in ("disorder_strengths", "gamma_phi", "disorder_samples"):
        if not np.array_equal(original[key], extension[key]):
            raise ValueError(f"Original and extension inputs differ: {key}")
    return original, extension, validation_path


def backend():
    # Lazily import Kwant only for --recompute. Archive extraction and plotting
    # work with ordinary NumPy/Matplotlib installations.
    if str(CODE) not in sys.path:
        sys.path.insert(0, str(CODE))
    from qsh_environment_response import EnvironmentDevice, quadrature, thermal_response
    from threadpoolctl import threadpool_limits
    metadata_path = SOURCE.with_name("metadata.json")
    parameters = json.loads(metadata_path.read_text(encoding="utf-8"))["parameters"]
    expected = dict(length=31, width=12, A=1., t=1., m=-1., delta=.2,
                    tc=.15, magnetic_field=0.)
    if parameters != expected:
        raise ValueError("The production QSH model/contact parameters changed")
    return EnvironmentDevice, quadrature, thermal_response, threadpool_limits


def check_response(report):
    if report["cpp_condition"] > 1e8:
        raise RuntimeError("Floating-probe conductance block is ill-conditioned")
    if max(report["probe_linear_error"], report["probe_quadratic_error"]) > 1e-8:
        raise RuntimeError("Floating probes have nonzero current")
    if report["current_sum_error"] > 1e-6 or max(report["checks"].values()) > 5e-8:
        raise RuntimeError("Gauge/conservation/elastic-probe reciprocity check failed")


def backend_sources():
    return [CODE / "qsh_environment_response.py", CODE / "compute_letter_preview.py",
            SOURCE.with_name("metadata.json")]
