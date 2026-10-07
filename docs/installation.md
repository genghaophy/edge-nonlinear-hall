# Installation

## Saved-data plotting

The plot and assembly stages need only NumPy and Matplotlib. From the repository root:

```sh
python -m pip install -r requirements-plot.txt
python code/run_panel_workflow.py --stage plot
```

Default cache extraction also uses the light dependency set. Its helpers read and hash included numerical sources without importing Kwant. An external LaTeX installation is not required for these plots.

## Numerical recomputation

The complete dependency set is recorded in `environment.yml` for Conda:

```sh
conda env create -f environment.yml
conda activate edge-nonlinear-hall
```

Alternatively, in an environment with the compiled Kwant dependencies available:

```sh
python -m pip install -r requirements-compute.txt
```

Kwant and tinyarray contain compiled components; installation requirements depend on the platform. The environment file requests conda-forge packages and exact versions, but is not a lockfile for every platform. Its Python version is left for the dependency solver to select. A successful install with another Python version is not the recorded tested runtime.

## Recorded source runtime

| Component | Tested version |
|---|---|
| Python | 3.13.5, Anaconda |
| NumPy | 2.1.3 |
| SciPy | 1.15.3 |
| Matplotlib | 3.10.0 |
| Kwant | 1.5.0 |
| tinyarray | 1.2.5 |
| threadpoolctl | 3.5.0 |

The requirements files pin these package versions. New installations and other platforms may select different BLAS/LAPACK implementations. Figure exports also depend on available fonts: drawing styles prefer Arial with a DejaVu Sans fallback. A numerical comparison should use arrays and stated tolerances rather than requiring identical exported file bytes.

Dependencies are not bundled in this repository and retain their own licenses. The MIT license covers the original project code.
