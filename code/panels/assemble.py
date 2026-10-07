"""Assemble saved panels at the approved layout; no numerical solver imports."""
from __future__ import annotations
import argparse
import importlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
import numpy as np

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from panels.common import (DATA_ROOT, OUTPUT_ROOT, approved_style, check_canvas,
                           export_figure, load_panel, source_record)

LAYOUTS = {
    'fig1': {'size_cm': (8.6, 7.5), 'stem': 'fig1_device',
             'boxes': {'a': (.025, .025, .95, .95)}},
    'fig2': {'size_cm': (10.8, 10.6), 'publication_width_cm': 8.6, 'stem': 'fig2_revised',
             'boxes': {'a': (.105, .61, .34, .33), 'b': (.655, .61, .31, .33),
                       'c': (.105, .14, .36, .36), 'd': (.655, .14, .31, .36)}},
    'fig3': {'size_cm': (8.6, 8.9), 'stem': 'fig3_revised',
             'boxes': {'a': (.14, .60, .34, .34), 'b': (.63, .60, .34, .34),
                       'c': (.14, .10, .34, .34), 'd': (.63, .10, .34, .34)}},
    'fig4': {'size_cm': (8.6, 9.6), 'stem': 'fig_probe_rm_qsh',
             'boxes': {'a': (.125, .575, .370, .365), 'b': (.590, .575, .370, .365),
                       'c': (.125, .105, .370, .335), 'd': (.590, .105, .370, .335)}},
}


def build(figure, data_root=DATA_ROOT):
    """Return figure/axes drawn only from that figure's prepared NPZ files."""
    if figure not in LAYOUTS:
        raise ValueError(f'Unknown figure: {figure}')
    plt.rcdefaults()
    layout = LAYOUTS[figure]
    if figure in ('fig1', 'fig2'):
        approved_style(9. if figure == 'fig1' else 8.2)
    else:
        importlib.import_module(f'panels.{figure}._style').style()
    fig = plt.figure(figsize=np.asarray(layout['size_cm']) / 2.54)
    axes, inputs = {}, []
    try:
        for panel, box in layout['boxes'].items():
            path = Path(data_root) / figure / f'{panel}.npz'
            data = load_panel(figure, panel, path)
            ax = fig.add_axes(box)
            module = importlib.import_module(f'panels.{figure}.plot_{panel}')
            module.draw(ax, data, label=figure in ('fig1', 'fig4'))
            axes[panel] = ax
            inputs.append({'panel': panel, **source_record(path),
                           'definitions': source_record(path.with_suffix('.json')),
                           'plot_script': source_record(module.__file__)})
        if figure == 'fig2':
            from panels.fig2._style import letter
            for panel, ax in axes.items():
                letter(ax, panel, horizontal_offset=-.036 if panel == 'd' else 0)
            axes['d'].yaxis.labelpad = .5
            scale = layout['size_cm'][0] / layout['publication_width_cm']
            for text in fig.findobj(match=Text):
                text.set_fontsize(text.get_fontsize() * scale)
            for line in fig.findobj(match=Line2D):
                line.set_linewidth(line.get_linewidth() * scale)
                if line.get_marker() not in (None, 'None', '', ' '):
                    line.set_markersize(line.get_markersize() * scale)
                    line.set_markeredgewidth(line.get_markeredgewidth() * scale)
            for ax in fig.findobj(match=matplotlib.axes.Axes):
                for spine in ax.spines.values():
                    spine.set_linewidth(spine.get_linewidth() * scale)
                ax.xaxis.labelpad *= scale
                ax.yaxis.labelpad *= scale
                ax.tick_params(which='both', width=.65 * scale, pad=1.8 * scale)
        elif figure == 'fig3':
            from panels.fig3._style import panel_label
            for panel, ax in axes.items():
                panel_label(ax, panel)
        return fig, axes, inputs
    except Exception:
        plt.close(fig)
        raise


def check_layout(fig, axes):
    checks = check_canvas(fig)
    renderer = fig.canvas.get_renderer()
    row_gaps = {}
    def decorations(panel):
        texts = axes[panel].findobj(match=Text)
        texts += [t for t in fig.texts if t.get_text() == f'({panel})']
        return [t for t in texts if t.get_visible() and t.get_text()]
    for upper, lower in (('a', 'c'), ('b', 'd')):
        if lower not in axes:
            continue
        lower_top = axes[lower].get_tightbbox(renderer).y1
        for text in fig.texts:
            if text.get_text() == f'({lower})':
                lower_top = max(lower_top, text.get_window_extent(renderer).y1)
        gap = (axes[upper].get_tightbbox(renderer).y0 - lower_top) / fig.dpi * 72
        # Tight axis bounding boxes include empty horizontal space. A negative
        # vertical clearance alone does not imply overlap of the actual labels.
        for upper_text in decorations(upper):
            upper_bounds = upper_text.get_window_extent(renderer)
            for lower_text in decorations(lower):
                if upper_bounds.overlaps(lower_text.get_window_extent(renderer)):
                    raise ValueError(f'Overlapping labels: {upper_text.get_text()} / {lower_text.get_text()}')
        row_gaps[f'{upper}-{lower}'] = float(gap)
    checks['row_clearance_pt'] = row_gaps
    checks['no_text_overlap_between_rows'] = True
    return checks


def assemble(figure, data_root=DATA_ROOT, output_dir=OUTPUT_ROOT / 'assembled'):
    fig, axes, inputs = build(figure, data_root)
    try:
        checks = check_layout(fig, axes)
        layout = LAYOUTS[figure]
        stem = Path(output_dir) / layout['stem']
        outputs = export_figure(fig, stem, dpi=600 if figure in ('fig1', 'fig2', 'fig3') else 450)
        metadata = {'figure': figure, 'layout': layout, 'inputs': inputs, 'outputs': outputs,
                    'layout_checks': checks, 'assembly_script': source_record(__file__),
                    'policy': 'Only prepared panel data are read. Approved production files are not replaced.'}
        Path(str(stem) + '_metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
        return metadata
    finally:
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figure', choices=('all', '1', '2', '3', '4'), default='all')
    parser.add_argument('--data-root', type=Path, default=DATA_ROOT)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT_ROOT / 'assembled')
    args = parser.parse_args()
    figures = list(LAYOUTS) if args.figure == 'all' else [f'fig{args.figure}']
    for figure in figures:
        record = assemble(figure, args.data_root, args.output_dir)
        print(json.dumps({'figure': figure, 'outputs': record['outputs']}, ensure_ascii=False))


if __name__ == '__main__':
    main()
