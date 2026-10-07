"""Panel I/O and CLIs. Matplotlib is imported only by plotting functions.

NPZ files contain display-ready numeric arrays; JSON records definitions and
sources. Loading a panel never imports a transport solver or computes data.
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import os
from pathlib import Path
import tempfile

import numpy as np

PROJECT = Path(__file__).resolve().parents[2]
CODE = PROJECT / 'code'
DATA_ROOT = PROJECT / 'data/panels'
OUTPUT_ROOT = PROJECT / 'outputs'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_record(path):
    path = Path(path).resolve()
    return {'path': os.path.relpath(path, PROJECT).replace('\\', '/'),
            'sha256': digest(path)}


def load_npz(path):
    with np.load(Path(path), allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def default_data(figure, panel):
    return DATA_ROOT / figure / f'{panel}.npz'


def default_plot(figure, panel):
    return OUTPUT_ROOT / 'panels' / figure / panel


def _json_value(value):
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, Path):
        return source_record(value)['path'] if value.is_file() else str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def save_panel(figure, panel, arrays, metadata, path=None):
    path = Path(path or default_data(figure, panel)).resolve()
    if path.suffix != '.npz':
        raise ValueError('Panel data output must end in .npz')
    prepared = {key: np.asarray(value) for key, value in arrays.items()}
    if any(value.dtype.hasobject for value in prepared.values()):
        raise ValueError('Panel data cannot contain pickled objects')
    for source in metadata.get('sources', []):
        if (PROJECT / source['path']).resolve() == path:
            raise ValueError('Do not overwrite an original source data file')
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix='.npz', dir=path.parent, delete=False) as temporary:
        temporary_path = Path(temporary.name)
    try:
        np.savez_compressed(temporary_path, **prepared)
        os.replace(temporary_path, path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()
    record = {'schema_version': 1, 'figure': figure, 'panel': panel,
              'arrays': {key: {'shape': list(value.shape), 'dtype': str(value.dtype)}
                         for key, value in prepared.items()},
              **_json_value(metadata), 'data_sha256': digest(path)}
    path.with_suffix('.json').write_text(json.dumps(record, indent=2, ensure_ascii=False), encoding='utf-8')
    return path


def load_panel(figure, panel, path=None):
    path = Path(path or default_data(figure, panel))
    if not path.is_file():
        raise FileNotFoundError(f'Missing {path}. Run calc_{panel}.py first.')
    metadata_path = path.with_suffix('.json')
    if not metadata_path.is_file():
        raise FileNotFoundError(f'Missing panel definitions: {metadata_path}')
    metadata = json.loads(metadata_path.read_text(encoding='utf-8'))
    if metadata['figure'] != figure or metadata['panel'] != panel:
        raise ValueError('Panel identifiers do not match the requested plot')
    if metadata['data_sha256'] != digest(path):
        raise ValueError('Panel NPZ hash differs from its metadata')
    data = load_npz(path)
    expected = metadata['arrays']
    if set(data) != set(expected):
        raise ValueError('Panel array names differ from metadata')
    for key, value in data.items():
        if list(value.shape) != expected[key]['shape'] or str(value.dtype) != expected[key]['dtype']:
            raise ValueError(f'Panel array contract changed: {key}')
    return data


def calc_main(figure, panel, calculate):
    parser = argparse.ArgumentParser(description=f'Prepare independent {figure}({panel}) data; never draws figures.')
    parser.add_argument('--recompute', action='store_true',
                        help='Recalculate this panel with the established solver; default reuses validated original results.')
    parser.add_argument('--output', type=Path, default=default_data(figure, panel),
                        help='Output NPZ; matching JSON definitions are written alongside it.')
    args = parser.parse_args()
    arrays, metadata = calculate(recompute=args.recompute)
    script = Path(inspect.getsourcefile(calculate))
    metadata = {**metadata, 'calculation_mode': 'recompute' if args.recompute else 'reuse',
                'calculation_script': source_record(script)}
    result = save_panel(figure, panel, arrays, metadata, args.output)
    print(json.dumps({'panel': f'{figure}{panel}', 'data': str(result),
                      'metadata': str(result.with_suffix('.json')), 'recomputed': args.recompute}, ensure_ascii=False))


def approved_style(font_size=9.):
    """Local copy of the approved Fig. 1/2 style; no legacy plot import."""
    import matplotlib as mpl
    size = float(font_size)
    mpl.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.size': size, 'mathtext.fontset': 'dejavusans',
        'axes.labelsize': size, 'axes.titlesize': size, 'axes.linewidth': .7,
        'axes.labelpad': 3., 'axes.unicode_minus': True,
        'axes.spines.top': True, 'axes.spines.right': True,
        'xtick.labelsize': size, 'ytick.labelsize': size,
        'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.major.width': .7, 'ytick.major.width': .7,
        'xtick.minor.width': .6, 'ytick.minor.width': .6,
        'xtick.major.size': 3.2, 'ytick.major.size': 3.2,
        'xtick.minor.size': 1.8, 'ytick.minor.size': 1.8,
        'xtick.top': True, 'ytick.right': True,
        'lines.linewidth': 1.1, 'lines.markersize': 3.,
        'legend.fontsize': size, 'legend.frameon': False,
        'legend.handlelength': 1.6, 'legend.borderaxespad': .3,
        'legend.labelspacing': .3, 'figure.facecolor': 'white',
        'axes.facecolor': 'white', 'savefig.facecolor': 'white',
        'savefig.edgecolor': 'white', 'savefig.dpi': 600,
        'savefig.bbox': None, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'svg.fonttype': 'none', 'svg.hashsalt': 'pnhe-prb', 'text.usetex': False,
    })


def export_figure(fig, stem, formats=('pdf', 'png', 'svg'), dpi=450):
    import matplotlib as mpl
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    result = []
    with mpl.rc_context({'savefig.bbox': None}):
        for extension in formats:
            path = stem.parent / f'{stem.name}.{extension}'
            fig.savefig(path, format=extension, dpi=dpi, bbox_inches=None)
            result.append({'path': str(path), 'sha256': digest(path)})
    return result


def check_canvas(fig):
    from matplotlib.text import Text
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    outside = []
    for text in fig.findobj(match=Text):
        if not text.get_visible() or not text.get_text():
            continue
        bounds = text.get_window_extent(renderer)
        if bounds.width <= 0 or bounds.height <= 0:
            continue
        if not fig.bbox.contains(*bounds.get_points()[0]) or not fig.bbox.contains(*bounds.get_points()[1]):
            outside.append(text.get_text())
    if outside:
        raise ValueError(f'Text outside output canvas: {outside}')
    return {'all_visible_text_inside_canvas': True}


def plot_main(figure, panel, draw, *, size_cm=(5.6, 5.0), axes_box=(.20, .19, .68, .70)):
    parser = argparse.ArgumentParser(description=f'Plot {figure}({panel}) from saved data; never runs numerical solvers.')
    parser.add_argument('--data', type=Path, default=default_data(figure, panel))
    parser.add_argument('--output', type=Path, default=default_plot(figure, panel), help='Output stem without extension.')
    args = parser.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    data = load_panel(figure, panel, args.data)
    fig = plt.figure(figsize=np.asarray(size_cm) / 2.54)
    ax = fig.add_axes(axes_box)
    draw(ax, data, label=True)
    qa = check_canvas(fig)
    outputs = export_figure(fig, args.output)
    metadata = {'figure': figure, 'panel': panel, 'size_cm': size_cm,
                'data': source_record(args.data), 'outputs': outputs, 'layout_checks': qa,
                'plot_script': source_record(inspect.getsourcefile(draw))}
    Path(str(args.output) + '_metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    plt.close(fig)
    print(json.dumps({'panel': f'{figure}{panel}', 'outputs': outputs}, ensure_ascii=False))
