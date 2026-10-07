"""Check the weak-probe formula from bare two-terminal RM scattering states.

By default, compare to the shipped historical seven-coupling audit. Use
--coupling-report to compare to a fresh checks/rm_weak_probe.py result.
"""
from pathlib import Path
import argparse
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from checks._common import DATA, OUTPUT, output_directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coupling-report", type=Path, default=DATA / "supplementary" / "weak_probe_audit" / "rm_tc_scan.json")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT / "rm_bare_injectivity")
    args = parser.parse_args()
    import compute_bare_injectivity_weak_probe_audit as audit
    audit.COUPLING_REPORT = args.coupling_report.resolve()
    audit.OUT = output_directory(args.output_dir)
    audit.main()


if __name__ == "__main__":
    main()
