"""Independent, one-energy diagnosis of the RM Figure 2 spatial panels.

Read the production arrays without changing them. Recompute equilibrium
injectivities and the floating linear voltages, independently call kwant.ldos,
and distinguish spectral localization from the local-neutrality potential.
No Poisson solver, energy scan, or new production data is involved.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time

import kwant
import numpy as np
from threadpoolctl import threadpool_limits
from compute_initial_study import Device

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs' / 'checks' / 'rm_spatial'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def comparison(actual, reference):
    difference = np.asarray(actual) - np.asarray(reference)
    scale = np.maximum(np.abs(reference), np.finfo(float).tiny)
    return dict(max_absolute_error=float(np.max(np.abs(difference))),
                max_pointwise_relative_error=float(np.max(np.abs(difference)/scale)),
                max_error_relative_to_peak=float(np.max(np.abs(difference)) /
                                                 np.max(np.abs(reference))))


def diagnose(input_path):
    started = time.time()
    with np.load(input_path, allow_pickle=False) as source:
        saved = {key: source[key].copy() for key in source.files}
    p = dict(tx=float(saved['tx']), t1=float(saved['t1']), t2=float(saved['t2']),
             delta=float(saved['delta']), tc=float(saved['tc']),
             L=int(saved['L']), N=int(saved['N']), EF=float(saved['EF']))
    dev = Device('four_x', cells=p['N'], delta=p['delta'], tc=p['tc'],
                 tx=p['tx'], t1=p['t1'], t2=p['t2'], length=p['L'])
    with threadpool_limits(limits=1):
        c = dev.matrix(p['EF'])
        alpha = np.array([1., -1., 0., 0.])
        alpha[2:] = np.linalg.solve(c[2:, 2:], -c[2:, :2] @ alpha[:2])
        wave = kwant.wave_function(dev.system, p['EF'], params=dev.params())
        nu_sites = np.asarray([sum((dev.density(state, params=dev.params())
                                   for state in wave(a)),
                                  np.zeros(len(dev.positions)))
                               for a in range(4)]) / (2*np.pi)
        # This uses Kwant's independent LDOS entry point, not stored u*rho.
        builtin_sites = kwant.ldos(dev.system, p['EF'], params=dev.params())
    nu = np.zeros((4, dev.width, dev.length))
    builtin_ldos = dev.zero.copy()
    for i in np.flatnonzero(dev.sample):
        x, y = dev.positions[i]
        nu[:, y, x] = nu_sites[:, i]
        builtin_ldos[y, x] = builtin_sites[i]
    rho = nu.sum(axis=0)
    if np.any(rho <= 0):
        raise RuntimeError('The local-neutrality ratio is undefined at zero LDOS.')
    u = nu / rho[None, :, :]
    numerator = np.einsum('a,ayx->yx', alpha, nu)
    source_numerator = nu[0] - nu[1]
    ud = numerator / rho
    probe_nu = nu[2] + nu[3]
    checks = dict(
        builtin_ldos_vs_fresh_injectivity_sum=comparison(builtin_ldos, rho),
        builtin_ldos_vs_saved_map=comparison(builtin_ldos, saved['ldos_map']),
        fresh_injectivities_vs_saved_u_times_rho=comparison(
            nu, saved['u_grid'] * saved['ldos_map'][None, :, :]),
        reconstructed_ud_vs_saved_field=comparison(ud, saved['field']),
        characteristic_sum_max_error=float(np.max(abs(u.sum(axis=0)-1))),
        floating_linear_current_max_error=float(np.max(abs((c @ alpha)[2:]))),
        matrix_reciprocity_max_error=float(np.max(abs(c-c.T))),
        Mx_LDOS_max_error=float(np.max(abs(rho-rho[:, ::-1]))),
        Mx_ud_odd_max_error=float(np.max(abs(ud+ud[:, ::-1]))))
    for name in ('builtin_ldos_vs_fresh_injectivity_sum', 'builtin_ldos_vs_saved_map',
                 'fresh_injectivities_vs_saved_u_times_rho', 'reconstructed_ud_vs_saved_field'):
        if checks[name]['max_absolute_error'] > 1e-10:
            raise RuntimeError(f'Independent spatial check failed: {name}: {checks[name]}')
    if checks['floating_linear_current_max_error'] > 1e-10:
        raise RuntimeError('Floating voltage linear-current check failed.')

    cells, length = p['N'], p['L']
    rows = 2*cells
    cell_rho = rho.reshape(cells, 2, length).sum(axis=1)
    rho_total = float(rho.sum())
    injection_abs_total = float(abs(source_numerator).sum())

    def region(row_ids):
        row_ids = np.asarray(row_ids, int)
        selection = rho[row_ids]
        return dict(atomic_rows=row_ids.tolist(), sites=int(selection.size),
            ldos_sum=float(selection.sum()), ldos_mean=float(selection.mean()),
            ldos_min=float(selection.min()), ldos_max=float(selection.max()),
            ldos_fraction_of_sample=float(selection.sum()/rho_total),
            source_injectivity_difference_abs_sum=float(abs(source_numerator[row_ids]).sum()),
            source_injectivity_difference_abs_fraction=float(
                abs(source_numerator[row_ids]).sum()/injection_abs_total),
            probe_origin_spectral_fraction=float(probe_nu[row_ids].sum()/selection.sum()),
            ud_abs_max=float(abs(ud[row_ids]).max()),
            ud_signed_mean=float(ud[row_ids].mean()))

    central_cells = np.arange(cells//2-1, cells//2+1)
    central_rows = np.ravel(np.column_stack((2*central_cells, 2*central_cells+1)))
    regions = dict(
        all_sample=region(np.arange(rows)),
        outermost_one_cell_each_edge=region(np.r_[0:2, rows-2:rows]),
        outermost_two_cells_each_edge=region(np.r_[0:4, rows-4:rows]),
        central_two_cells=region(central_rows))
    regions['central_mean_to_outermost_one_cell_mean_ratio'] = (
        regions['central_two_cells']['ldos_mean'] /
        regions['outermost_one_cell_each_edge']['ldos_mean'])

    # Exact finite-y eigenvectors and x-dispersions provide a contact-free control.
    hy = np.diag(p['delta'] * (-1.)**np.arange(rows))
    for y in range(rows-1):
        hy[y, y+1] = hy[y+1, y] = p['t1'] if y % 2 == 0 else p['t2']
    eigenvalues, psi = np.linalg.eigh(hy)
    probabilities = abs(psi)**2
    boundary_polarization = probabilities[-2:].sum(axis=0)-probabilities[:2].sum(axis=0)
    edge_ids = np.argsort(abs(eigenvalues))[:2]
    top = int(edge_ids[np.argmax(boundary_polarization[edge_ids])])
    bottom = int(edge_ids[np.argmin(boundary_polarization[edge_ids])])
    propagating = abs(p['EF']-eigenvalues) < 2*abs(p['tx'])
    dos_1d = np.zeros(rows)
    dos_1d[propagating] = 1/(np.pi*np.sqrt((2*p['tx'])**2 -
                                              (p['EF']-eigenvalues[propagating])**2))
    intrinsic_rho = probabilities @ dos_1d
    intrinsic_cell_rho = intrinsic_rho.reshape(cells, 2).sum(axis=1)
    bulk_mask = np.ones(rows, bool)
    bulk_mask[[top, bottom]] = False
    bulk_half_gap = np.sqrt((abs(p['t2'])-abs(p['t1']))**2+p['delta']**2) - 2*abs(p['tx'])
    gap = dict(global_bulk_half_gap=float(bulk_half_gap),
        EF_gap_margin=float(bulk_half_gap-abs(p['EF'])),
        EF_inside_global_bulk_gap=bool(abs(p['EF']) < bulk_half_gap),
        transverse_eigenvalues=eigenvalues.tolist(),
        top_bottom_edge_ids=[top, bottom],
        propagating_transverse_ids=np.flatnonzero(propagating).tolist(),
        propagating_bulk_subbands=int(np.count_nonzero(propagating & bulk_mask)),
        propagating_edge_subbands=int(np.count_nonzero(propagating & ~bulk_mask)),
        propagation_criterion='|EF-lambda_n| < 2|tx| for E_n(kx)=lambda_n-2tx cos(kx)',
        probability_tail_factor_per_y_cell=float((p['t1']/p['t2'])**2),
        top_edge_energy_shift_from_minus_delta=float(eigenvalues[top]+p['delta']),
        bottom_edge_energy_shift_from_plus_delta=float(eigenvalues[bottom]-p['delta']))
    if not gap['EF_inside_global_bulk_gap'] or gap['propagating_bulk_subbands']:
        raise RuntimeError('Representative sample has a propagating bulk subband.')
    checks['recomputed_intrinsic_cell_LDOS_vs_saved'] = comparison(intrinsic_cell_rho,
                                                                  saved['ldos_cell'])
    xc = length//2
    probe_points = []
    # Inactive sublattice rows: contact-origin evanescent spectral tails.
    inactive = np.r_[np.arange(1, rows//2, 2), np.arange(rows//2, rows-1, 2)]
    for y in inactive:
        probe_points.append(dict(x=xc, y=int(y), cell=int(y//2),
            ldos=float(rho[y, xc]), intrinsic_ribbon_ldos=float(intrinsic_rho[y]),
            device_to_intrinsic_ldos_ratio=float(rho[y, xc]/intrinsic_rho[y]),
            probe_origin_spectral_fraction=float(probe_nu[y, xc]/rho[y, xc]),
            interpretation='Reservoir-origin local spectral weight; not a propagating bulk band'))

    far_x = np.r_[np.arange(5, 11), np.arange(length-11, length-5)]
    far_x = np.unique(far_x[(far_x >= 0) & (far_x < length)])
    far_rho = rho[:, far_x].mean(axis=1)
    max_y, max_x = np.unravel_index(np.argmax(abs(ud)), ud.shape)
    tail_points = []
    cut_x = max(0, xc-5)
    # The potential can stay finite after numerator and denominator both decay.
    for y in range(rows):
        tail_points.append(dict(x=int(cut_x), y=int(y), cell=int(y//2),
            ldos=float(rho[y, cut_x]),
            source_injectivity_difference=float(source_numerator[y, cut_x]),
            source_injectivity_difference_over_ldos=float(source_numerator[y, cut_x]/rho[y, cut_x]),
            ud=float(ud[y, cut_x]), injectivity_fractions=u[:, y, cut_x].tolist()))
    maximum = dict(x=int(max_x), y=int(max_y), ud=float(ud[max_y, max_x]),
                   ldos=float(rho[max_y, max_x]),
                   weighted_injectivity_numerator=float(numerator[max_y, max_x]),
                   injectivity_fractions=u[:, max_y, max_x].tolist())
    old_cutoff = 1e-8*float(rho.max())
    ratio_diagnostics = dict(
        formula='ud(x,y)=sum_a alpha_a nu_a(x,y) / rho(x,y), rho=sum_a nu_a',
        alpha=alpha.tolist(), denominator_min=float(rho.min()),
        denominator_max=float(rho.max()), denominator_min_to_max=float(rho.min()/rho.max()),
        previous_relative_display_cutoff=1e-8,
        previous_absolute_display_cutoff=old_cutoff,
        sites_below_previous_cutoff=int(np.count_nonzero(rho < old_cutoff)),
        maximum_abs_ud_site=maximum,
        center_Mx_plane_ud_abs_max=float(abs(ud[:, xc]).max()),
        source_only_ratio_vs_full_ud_abs_max_error=float(
            np.max(abs(source_numerator/rho-ud))),
        tail_cut_x=int(cut_x), tail_cut=tail_points,
        tail_cancellation_explanation=('Local neutrality divides the voltage-weighted injectivity '
            'by the local equilibrium DOS. A common exponentially small boundary tail can '
            'cancel in that ratio; ud is a potential coefficient, not a probability density.'))
    profiles = dict(atomic_ldos_mean_over_x=rho.mean(axis=1).tolist(),
        cell_ldos_mean_over_x=cell_rho.mean(axis=1).tolist(),
        intrinsic_atomic_ldos=intrinsic_rho.tolist(),
        intrinsic_cell_ldos=intrinsic_cell_rho.tolist(),
        intrinsic_central_two_cells_ldos_fraction=float(
            intrinsic_rho[central_rows].sum()/intrinsic_rho.sum()),
        far_from_probe_x=far_x.tolist(), atomic_ldos_far_from_probe=far_rho.tolist(),
        far_central_two_cells_ldos_fraction=float(far_rho[central_rows].sum()/far_rho.sum()),
        cell_probe_origin_fraction=(probe_nu.reshape(cells, 2, length).sum(axis=(1,2)) /
                                    rho.reshape(cells, 2, length).sum(axis=(1,2))).tolist(),
        contact_inactive_sublattice_points=probe_points,
        outer_atomic_cuts=dict(x=list(range(length)), bottom_y=0, top_y=rows-1,
                              bottom_ud=ud[0].tolist(), top_ud=ud[-1].tolist(),
                              bottom_ldos=rho[0].tolist(), top_ldos=rho[-1].tolist(),
                              bottom_source_injectivity_difference=source_numerator[0].tolist(),
                              top_source_injectivity_difference=source_numerator[-1].tolist()))
    report = dict(schema='rm_spatial_diagnostic_v1', parameters=p,
        protocol='VL=+V, VR=-V, IT=IB=0; linear floating voltages solved from C_PP',
        units='LDOS and injectivity per t2 per atomic site; U is positive onsite energy; v=eV/t2',
        spatial_definition='2N atomic rows; cell n comprises rows 2n(A),2n+1(B); bottom row0, top row2N-1',
        approximation='Zero-temperature coherent transport, first-order local neutrality; not full Poisson/Hartree',
        checks=checks, regions=regions, gap_and_subbands=gap,
        ratio_diagnostics=ratio_diagnostics, profiles=profiles,
        interpretation=[
            'Actual equilibrium LDOS and the absolute source-injectivity difference are localized at the boundaries.',
            'Color in a low-DOS region of ud is not evidence of a conducting bulk channel.',
            'The central LDOS is finite and independently reproduced by kwant.ldos; it is not simply machine noise.',
            'Local spectral tails from point contacts are evanescent in the bulk gap; reservoir origin does not establish bulk conduction.',
            'Numerical agreement of ud does not establish electrostatic validity. Low DOS weakens the local-neutrality strong-screening approximation.',
            'A true potential can extend through an insulating bulk. Its profile requires a specified dielectric/capacitance geometry or nonlocal Poisson response.',
            'The source-injectivity difference is an unscreened density response; the local-neutrality closure cancels the net first-order induced charge.'
        ],
        source_files={name:sha256(ROOT/'code'/name) for name in
                      [Path(__file__).name, 'compute_initial_study.py',
                       'compute_floating_study.py']},
        input_file=str(input_path.resolve()), input_sha256=sha256(input_path),
        kwant_version=kwant.__version__, numpy_version=np.__version__,
        runtime_seconds=time.time()-started)
    arrays = dict(x=np.arange(length), y_atomic=np.arange(rows), nu_grid=nu,
        independent_kwant_ldos=builtin_ldos, wavefunction_ldos=rho,
        reconstructed_ud=ud, alpha=alpha, characteristic_potentials=u,
        source_injectivity_difference=source_numerator, weighted_injectivity_numerator=numerator,
        probe_origin_fraction=probe_nu/rho, intrinsic_atomic_ldos=intrinsic_rho,
        ud_top=ud[-1], ud_bottom=ud[0], ldos_top=rho[-1], ldos_bottom=rho[0],
        source_injectivity_difference_top=source_numerator[-1],
        source_injectivity_difference_bottom=source_numerator[0])
    return report, arrays


# Use code/checks/rm_spatial.py for portable inputs and report output.
