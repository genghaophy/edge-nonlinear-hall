"""Figure 2 source paths and small read-only data utilities."""
from pathlib import Path
import hashlib
import numpy as np
from panels.common import source_record

PROJECT = Path(__file__).resolve().parents[3]
RM_ROOT = PROJECT
RIBBON_SOURCE = PROJECT / "data/raw/rm_figure2/figure2_data.npz"
MAP_SOURCE = PROJECT / "data/raw/figure2_width30/lpdos_data.npz"
DELTA_SOURCE = PROJECT / "data/raw/figure2_delta_scan/scan.npz"

PARAMETERS = dict(tx=.25, t1=.3, t2=1., delta=.2, tc=.15,
                  EF=.1, N=8, L=31, kstar=2., eta=.02, nk=801)


def load(path, keys=None):
    with np.load(path, allow_pickle=False) as archive:
        names = archive.files if keys is None else keys
        arrays = {name: archive[name].copy() for name in names}
    for name, array in arrays.items():
        if np.issubdtype(array.dtype, np.number) and not np.isfinite(array).all():
            raise ValueError(f"Nonfinite source array {path}:{name}")
    return arrays


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def provenance(panel, source, recompute, **extra):
    record = dict(figure="fig2", panel=panel,
                  data_mode="recomputed" if recompute else "verified_cache_extraction",
                  model="stacked Rice-Mele", magnetic_field=0.,
                  parameters=PARAMETERS.copy(), units="energies in t2, lengths in lattice cells",
                  source_path=source_record(source)['path'] if source.exists() else source.name, source_sha256=digest(source) if source.exists() else None,
                  sources=[source_record(source)] if source.exists() else [],
                  helpers=[source_record(Path(__file__))],
                  input_sources_modified=False)
    record.update(extra)
    return record


def select(arrays, keys):
    return {name: np.asarray(arrays[name]).copy() for name in keys}


def parameter_arrays(parameters=None):
    p = PARAMETERS if parameters is None else parameters
    return {name: np.asarray(p[name]) for name in ("tx", "t1", "t2", "delta", "tc", "EF", "N", "L")}
