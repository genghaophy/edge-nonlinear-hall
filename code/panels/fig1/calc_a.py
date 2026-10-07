"""Prepare the approved schematic geometry; this is not a DOS calculation."""
from pathlib import Path
import sys
import numpy as np
if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from panels.common import CODE, calc_main, source_record


def calculate(recompute=False):
    annotations = [r'$\mathcal{V}_\perp$', r'$I_3=0$ (floating)',
                   r'$I_4=0$ (floating)', r'$+\mathcal{V}/2$', r'$-\mathcal{V}/2$',
                   'Bulk insulator', r'$\rho_T(E_F)<\rho_B(E_F)$',
                   r'$B=0$; TRS; $\mathcal{M}_y$ broken', '+', '\N{MINUS SIGN}',
                   'Darker edge: higher DOS']
    arrays = {
        'body': [2., 2.7, 5., 3.],
        'edges': [[2., 5.52, 5., .18], [2., 2.7, 5., .18]],
        'terminals': [[1.37, 2.7, .63, 3.], [7., 2.7, .63, 3.],
                      [4.2, 5.7, .6, .7], [4.2, 2., .6, .7]],
        'wire_x': [[4.5, 4.5, 9.2, 9.2], [4.5, 4.5, 9.2, 9.2]],
        'wire_y': [[6.4, 7.05, 7.05, 4.72], [2., 1.35, 1.35, 3.68]],
        'voltmeter': [9.2, 4.2, .52],
        'annotation_xy': [[9.2, 4.2], [4.5, 7.59], [4.5, .82], [1.685, 2.14],
                          [7.315, 2.14], [4.5, 4.72], [4.5, 4.13], [4.5, 3.52],
                          [9.2, 5.04], [9.2, 3.36], [4.5, .16]],
        'annotation_text': np.asarray(annotations),
        # The production writer raises the bottom note from 8 to 8.5 pt.
        'annotation_fontsize': [9.] * 10 + [8.5],
        'x_limits': [0., 10.2], 'y_limits': [0., 8.15],
    }
    metadata = {'kind': 'qualitative device schematic',
                'sources': [source_record(Path(__file__))],
                'parameters': {'terminal_order': ['L', 'R', 'T', 'B'],
                               'bias': 'V1=+V/2; V2=-V/2; I3=I4=0'},
                'interpretation': 'Equal edge widths; shade indicates qualitative equilibrium DOS inequality, not calculated DOS.',
                'recompute_meaning': 'Regenerate schematic coordinates; no transport solver is required.'}
    return arrays, metadata


if __name__ == '__main__':
    calc_main('fig1', 'a', calculate)
