"""Reproduce the RM weak-probe audit using the verified original seven-coupling scan.

The full response retains source-electrode and internal-potential feedback.
This selected-energy diagnostic does not run the article's full energy scan.
"""
from pathlib import Path
import argparse
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from checks._common import DATA, OUTPUT, output_directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DATA / "rm_figure2" / "figure2_data.npz")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT / "rm_weak_probe")
    args = parser.parse_args()
    import compute_rm_weak_probe_audit as audit
    audit.SOURCE = args.input.resolve()
    audit.OUT = output_directory(args.output_dir)
    audit.main()


if __name__ == "__main__":
    main()
