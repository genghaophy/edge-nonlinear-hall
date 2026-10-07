"""Select calculation, independent panel plotting, or approved whole-figure assembly."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
import sys

CODE = Path(__file__).resolve().parent
PANELS = {'1': 'a', '2': 'abcd', '3': 'abcd', '4': 'abcd'}


def execute(script, arguments):
    result = subprocess.run([sys.executable, '-X', 'utf8', str(script), *arguments],
                            cwd=CODE, capture_output=True, encoding='utf-8', errors='replace')
    if result.returncode:
        sys.stderr.write(result.stdout + result.stderr)
        raise SystemExit(result.returncode)
    # Keep numerical progress/provenance output available in the run record.
    return {'script': str(script), 'arguments': arguments,
            'stdout': result.stdout, 'stderr': result.stderr}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figure', choices=('all', *PANELS), default='all')
    parser.add_argument('--panel', choices=tuple('abcd'))
    parser.add_argument('--stage', choices=('data', 'plot', 'assemble', 'all'), default='all')
    parser.add_argument('--recompute', action='store_true',
                        help='Recalculate the selected panel(s); may run long transport scans. Default prepares saved production data.')
    args = parser.parse_args()
    if args.panel and args.figure == 'all':
        parser.error('--panel requires one --figure')
    if args.panel and args.panel not in PANELS[args.figure]:
        parser.error(f'Figure {args.figure} has no panel {args.panel}')
    if args.panel and args.stage == 'assemble':
        parser.error('Whole-figure assembly needs all panels; omit --panel')
    if args.recompute and args.stage not in ('data', 'all'):
        parser.error('--recompute applies only to the data/all stages')
    figures = list(PANELS) if args.figure == 'all' else [args.figure]
    operations = []
    for number in figures:
        selected = args.panel or PANELS[number]
        for panel in selected:
            directory = CODE / f'panels/fig{number}'
            if args.stage in ('data', 'all'):
                print(f"Fig.{number}({panel}): {'recompute' if args.recompute else 'prepare saved data'}", flush=True)
                operations.append(execute(directory / f'calc_{panel}.py', ['--recompute'] if args.recompute else []))
            if args.stage in ('plot', 'all'):
                print(f'Fig.{number}({panel}): plot saved panel data', flush=True)
                operations.append(execute(directory / f'plot_{panel}.py', []))
        if args.stage in ('assemble', 'plot', 'all') and args.panel is None:
            print(f'Fig.{number}: assemble approved layout', flush=True)
            operations.append(execute(CODE / 'panels/assemble.py', ['--figure', number]))
    from panels.common import OUTPUT_ROOT
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    report = OUTPUT_ROOT / 'last_run.json'
    report.write_text(json.dumps({'selection': vars(args), 'operations': operations}, indent=2), encoding='utf-8')
    print(f'Completed {len(operations)} operations; run record: {report}', flush=True)


if __name__ == '__main__':
    main()
