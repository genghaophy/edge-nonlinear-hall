"""Check packaged data, independent relocation, and selected transport points.

No full production sweep is launched. Generated records belong in outputs/.
"""
from __future__ import annotations
import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
from uuid import uuid4

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'code'))
from panels.common import load_panel

PANELS = {'fig1': 'a', 'fig2': 'abcd', 'fig3': 'abcd', 'fig4': 'abcd'}
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def packaged_data():
    checked = []
    for figure, panels in PANELS.items():
        for panel in panels:
            data = load_panel(figure, panel)
            for key, value in data.items():
                if np.issubdtype(value.dtype, np.number) and not np.isfinite(value).all():
                    # Band-plot masks intentionally use NaN to separate bulk
                    # and edge pieces. They are not missing transport results.
                    masked_band = figure == 'fig3' and panel == 'a' and key.endswith('curves')
                    if not masked_band or np.isinf(value).any():
                        raise ValueError(f'Nonfinite panel data: {figure}{panel}:{key}')
            checked.append(figure + panel)
    embedded_checked = 0
    pattern = re.compile(r'[A-Za-z]:[\\/]|/lustre/|/home/')
    for path in (ROOT / 'data').rglob('*.npz'):
        with np.load(path, allow_pickle=False) as data:
            for key in data.files:
                if data[key].dtype.kind in 'US':
                    for value in data[key].ravel():
                        embedded_checked += 1
                        if pattern.search(str(value)):
                            raise ValueError(f'Absolute historical path in {path.relative_to(ROOT)}:{key}')
    manifest_path = ROOT / 'provenance/manifest.json'
    manifest_checked = 0
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        for item in manifest['files']:
            if digest(ROOT / item['path']) != item['sha256']:
                raise ValueError(f'Release file changed since manifest: {item["path"]}')
            manifest_checked += 1
    return dict(panel_contracts_checked=checked, embedded_string_values_checked=embedded_checked,
                manifest_files_checked=manifest_checked)


def run_workflow(repository, stage, cwd, env, logfile):
    command = [sys.executable, '-X', 'utf8', str(repository / 'code/run_panel_workflow.py'), '--stage', stage]
    result = subprocess.run(command, cwd=cwd, env=env, capture_output=True,
                            encoding='utf-8', errors='replace')
    logfile.write_text(result.stdout + result.stderr, encoding='utf-8')
    if result.returncode:
        raise RuntimeError(f'{stage} failed; see {logfile}')


@contextmanager
def temporary_workspace(parent):
    """Use ordinary directory permissions, including managed Windows hosts."""
    parent = Path(parent).resolve()
    path = parent / ('edge-nhe-release-' + uuid4().hex)
    path.mkdir()
    if not path.resolve().is_relative_to(parent):
        raise RuntimeError('Unexpected temporary directory location.')
    try:
        yield path
    finally:
        # Re-check the absolute target before recursive cleanup.
        if path.exists() and path.resolve().parent == parent:
            shutil.rmtree(path)


def relocated(plot=False, temp_root=None):
    output = ROOT / 'outputs/validation'
    output.mkdir(parents=True, exist_ok=True)
    temp_root = Path(temp_root or tempfile.gettempdir()).resolve()
    temp_root.mkdir(parents=True, exist_ok=True)
    with temporary_workspace(temp_root) as temporary:
        temporary = Path(temporary)
        if not temporary.resolve().is_relative_to(temp_root):
            raise RuntimeError('Unexpected temporary directory location.')
        repository = temporary / 'repository'
        shutil.copytree(ROOT, repository,
                        ignore=shutil.ignore_patterns('.git', 'outputs', '__pycache__', '*.pyc', '*.zip'))
        guard = temporary / 'solver_import_guard'
        guard.mkdir()
        (guard / 'sitecustomize.py').write_text(
            'import sys, importlib.abc\n'
            'class NoSolver(importlib.abc.MetaPathFinder):\n'
            '    def find_spec(self, fullname, path=None, target=None):\n'
            '        if fullname.split(".")[0] in {"kwant", "tinyarray", "edge_probe_readout", '
            '"qsh_environment_response", "compute_floating_study", "compute_initial_study"}:\n'
            '            raise ModuleNotFoundError("Solver import prohibited in cache/plot stages: " + fullname)\n'
            'sys.meta_path.insert(0, NoSolver())\n', encoding='utf-8')
        env = dict(os.environ, PYTHONPATH=str(guard), PYTHONUTF8='1', MPLBACKEND='Agg')
        print('Checking data preparation after relocation, with solver imports disabled.', flush=True)
        run_workflow(repository, 'data', temporary, env, output / 'relocated_data.log')
        compared = 0
        for figure, panels in PANELS.items():
            for panel in panels:
                left = ROOT / f'data/panels/{figure}/{panel}.npz'
                right = repository / f'data/panels/{figure}/{panel}.npz'
                with np.load(left, allow_pickle=False) as a, np.load(right, allow_pickle=False) as b:
                    if set(a.files) != set(b.files):
                        raise ValueError('Relocated panel array keys changed.')
                    for key in a.files:
                        equal = (np.array_equal(a[key], b[key], equal_nan=True)
                                 if np.issubdtype(a[key].dtype, np.number) else np.array_equal(a[key], b[key]))
                        if not equal:
                            raise ValueError(f'Relocated panel changed: {figure}{panel}:{key}')
                compared += 1
        result = dict(independent_directory=True, solver_imports_disabled=True,
                      prepared_panel_arrays_identical=compared)
        if plot:
            print('Checking all thirteen panel plots and four figure assemblies.', flush=True)
            run_workflow(repository, 'plot', temporary, env, output / 'relocated_plot.log')
            figures = repository / 'outputs'
            paths = list((figures / 'panels').rglob('*.pdf')) + list((figures / 'assembled').glob('*.pdf'))
            if len(paths) != 17:
                raise ValueError('Expected thirteen panel PDFs and four assembled figures.')
            destination = output / 'relocated_figures'
            shutil.copytree(figures, destination, dirs_exist_ok=True)
            result['panel_pdfs'] = 13
            result['assembled_pdfs'] = 4
        return result


def compute_points():
    from threadpoolctl import threadpool_limits
    import disorder_probe_readout as api
    records = {}
    with threadpool_limits(limits=1):
        for model, energy, panel in [('rm', .10, 'c'), ('qsh', .40, 'd')]:
            parameters = api.RMParameters() if model == 'rm' else api.QSHParameters(width=12)
            disorder = api.make_disorder(model, parameters, 0., 0)
            report = api.four_quantities(model, parameters, energy, disorder, .003, 1e-5)
            expected = float(load_panel('fig4', panel)['clean_kappa'])
            difference = abs(float(report['kappa_full']) - expected)
            if difference > 1e-7:
                raise ValueError(f'{model} clean coefficient differs from packaged production data: {difference}')
            for key, value in report.get('checks', {}).items():
                if key.endswith('error') and np.ndim(value) == 0 and abs(float(value)) > 2e-7:
                    raise ValueError(f'{model} transport diagnostic failed: {key}')
            bare = api.bare_quantities(model, parameters, energy, disorder, 1e-5)
            records[model] = {'energy':energy, 'tau':.003, 'clean_kappa':float(report['kappa_full']),
                              'production_difference':difference,
                              'weak_probe_prediction':float(bare['kappa_general'])}
            print(f'Checked {model.upper()} clean weak-probe point.', flush=True)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--relocate', action='store_true', help='Check cache preparation in a temporary independent copy.')
    parser.add_argument('--plot', action='store_true', help='Also check every plot and figure assembly in that copy.')
    parser.add_argument('--compute', action='store_true', help='Check one clean weak-probe point per model; requires computation dependencies.')
    parser.add_argument('--temp-root', type=Path, help='Optional writable parent for the independent temporary copy.')
    args = parser.parse_args()
    started = time.perf_counter()
    report = {'data_integrity':packaged_data()}
    if args.relocate or args.plot:
        report['relocation'] = relocated(args.plot, args.temp_root)
    if args.compute:
        report['representative_computation'] = compute_points()
    report['elapsed_seconds'] = time.perf_counter() - started
    output = ROOT / 'outputs/validation'
    output.mkdir(parents=True, exist_ok=True)
    (output / 'last_validation.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
