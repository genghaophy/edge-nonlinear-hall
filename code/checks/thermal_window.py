"""Audit the fixed-window integration-by-parts boundary term from saved spectra.

Reports both the raw-code amplitude difference and the paper full-bias
dimensionless difference. This is not an estimate of all omitted energies.
No transport solver or Kwant installation is needed.
"""
from pathlib import Path
import argparse
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from checks._common import DATA, OUTPUT, output_directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DATA / "supplementary" / "thermal_spectra" / "clean_o16.npz")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT / "thermal_window")
    args = parser.parse_args()
    from audit_thermal_window import audit
    target = output_directory(args.output_dir) / "report.json"
    audit(args.input, target)
    print(target)


if __name__ == "__main__":
    main()
