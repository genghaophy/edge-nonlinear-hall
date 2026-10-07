"""Recompute RM injectivities/LDOS and distinguish density localization from potential."""
from pathlib import Path
import argparse
import sys
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from checks._common import DATA, OUTPUT, digest, output_directory, protect_inputs, write_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DATA / "rm_figure2" / "figure2_data.npz")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT / "rm_spatial")
    args = parser.parse_args()
    import numpy as np
    from compute_rm_spatial_diagnostic import diagnose
    from threadpoolctl import threadpool_limits
    output = output_directory(args.output_dir)
    with protect_inputs([args.input]) as fingerprints, threadpool_limits(limits=1):
        report, arrays = diagnose(args.input)
    # Original reconstructed_ud is the raw source-amplitude drive field.
    # Store the physical characteristic potential separately without modifying it.
    arrays["paper_u_V"] = arrays["reconstructed_ud"] / 2
    arrays["paper_u_V_top"] = arrays["ud_top"] / 2
    arrays["paper_u_V_bottom"] = arrays["ud_bottom"] / 2
    target = output / "spatial_diagnostic.npz"
    np.savez_compressed(target, **arrays)
    report.update(input_sha256=fingerprints, diagnostic_npz_sha256=digest(target),
                  paper_conversion="paper_u_V=raw reconstructed_ud/2; raw arrays retained unchanged")
    write_report(output / "report.json", report)
    print(target)


if __name__ == "__main__":
    main()
