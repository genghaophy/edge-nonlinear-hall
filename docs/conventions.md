# Physical and statistical conventions

## Voltages and quadratic response

The positive symbol `e` is the magnitude of the electron charge. Physical electrical potentials obey:

```text
mu_a = E_F - e V_a
H = H_0 - e U
V_1 = +V/2, V_2 = -V/2
V_perp = V_3 - V_4
V_perp = kappa_1 V + kappa_2 V^2 + higher orders
```

Here `V` is the full source–drain voltage, not a single-terminal amplitude. Terminals are ordered 1 = left, 2 = right, 3 = top, 4 = bottom. Measurement probes float: their currents vanish. The physical `kappa_2` has units of inverse voltage.

The plotted dimensionless coefficient is:

| Model | Energy unit | Plotted coefficient |
|---|---|---|
| RM | `t2` | `tilde_kappa2 = t2*kappa2/e` |
| QSH | `t` | `tilde_kappa2 = t*kappa2/e` |

Some historic arrays use dimensionless electron-energy variables, `v_a = -e*V_a/t` and `u_code = -e*U/t`, with the corresponding RM scale `t2` when appropriate. Their source amplitude is a single-terminal amplitude, `v_code = -e*V/(2*t)`. Consequently the historical quadratic coefficient `kH` maps to the paper convention as:

```text
tilde_kappa2 = -kH/4
```

The minus sign comes from the electron-energy/electrical-voltage convention, while the factor of four comes from the full-bias square. This conversion occurs once in the calculation/readout layer. Fig. 4 caches already contain the plotted coefficient and are copied without another conversion. Plot scripts do not alter the sign or voltage normalization. Array names such as `kappa_full` may refer to the dimensionless plotted coefficient; consult their JSON units rather than inferring dimensions from a name.

## Density observables

Fig. 2(b) distinguishes fixed-momentum broadened spectral density from momentum-integrated edge LDOS. Fig. 2(c) is the **left-injected local partial density of states (LPDOS)** at equilibrium: the sum over unit-flux incoming left modes, divided by `2*pi`. Its unit is `t2^-1`. A and B atomic contributions are summed within each cell, using a shared absolute color scale without edge-by-edge rescaling or interpolation.

The Fig. 2(c) geometry is 30 transverse cells × 31 longitudinal cells, corresponding to 60 atomic rows. The other RM main-figure calculations use the established width of eight cells. This LPDOS is not the total equilibrium LDOS, a net current, or the finite-bias electron density. The calculation checks that the sum of all contact injectivities agrees with the total LDOS.

Fig. 1's shading illustrates unequal equilibrium edge DOS. It is a qualitative schematic, not a calculated density map.

## Disorder and error bars

The on-site nonmagnetic disorder strength uses the model's own energy scale. The disorder draws are uniform in `[-W_dis/2, W_dis/2]`.

| Panel | Nonzero-disorder samples | Displayed spread | Strengths in model units |
|---|---|---|---|
| Fig. 3(c), QSH | 8, seeds 1700–1707 | Standard error of the mean, `SD/sqrt(8)` | `W_dis/t = 0, 0.2, 0.5, 1.0` |
| Fig. 4(c), RM | 16, seeds 1700–1715 | Sample standard deviation, `ddof=1` | `W_dis/t2 = 0, 0.025, 0.05, 0.10` |
| Fig. 4(d), QSH | 16, seeds 1700–1715 | Sample standard deviation, `ddof=1` | `W_dis/t = 0, 0.10, 0.25, 0.50` |

The clean reference is calculated once, with zero statistical spread. Repeated clean array slots do not count as independent realizations. Fig. 4 disorder responses are normalized by the signed clean coefficient at the same Fermi energy and probe coupling `tau/t_model = 0.003`.

## Model and interpretation limits

The transport calculations use the stated coherent elastic models at zero magnetic field. Electrostatic feedback uses first-order local charge neutrality, a strong-screening closure rather than a nonlinear Poisson/Hartree solution. Finite temperature and spin-conserving elastic virtual probes are specific model extensions; the latter are not a general model of all inelastic or magnetic scattering.

The weak-probe voltage depends on both local spectral information and the bias dependence of injection. A pure spectral reduction requires additional ideal-contact conditions, including reflected amplitudes that vanish through first order in bias; perfect equilibrium transmission alone is insufficient. The ordinary RM example does not require QSH topology. The QSH calculations test improved stability for the specified nonmagnetic disorder, contacts, widths, and sample counts; the quadratic coefficient is neither quantized nor established to be topologically protected.

For the strongest Fig. 3(c) disorder, persistence of a finite-device signal does not itself prove that every realization retains a bulk gap and an edge-only transport sector. Temperature quadrature tests concern the stated finite energy window and do not establish unrestricted behavior near all thresholds.
