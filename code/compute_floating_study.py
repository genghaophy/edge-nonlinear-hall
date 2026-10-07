"""Four-terminal voltage probes: source voltages +v,-v, I_T=I_B=0.

v=eV/t2 is the electrode amplitude (source-drain difference is 2v).
All terminal currents point from reservoirs into the sample. The readout is
v_H=v_T-v_B, with coefficient k_H in v_H=alpha_H*v+k_H*v**2+O(v**3).
Screening uses equilibrium injectivities and the local-neutrality closure;
finite-bias checks retain its linear voltage dependence, not full Poisson.
"""
from pathlib import Path
import importlib.util
import json
import time
import kwant
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import root
from threadpoolctl import threadpool_limits
from compute_initial_study import Device, ROOT, OUT, EF, STEP


def floating_response(dev, energy, step=STEP):
    if dev.nlead != 4:
        raise ValueError('Lead order must be L,R,T,B.')
    positions = dev.positions if hasattr(dev, 'positions') else dev.pos
    c = dev.matrix(energy)
    cpp = c[2:, 2:]
    if np.linalg.cond(cpp) > 1e12:
        raise RuntimeError('Floating voltages are not uniquely fixed by connected source channels.')
    alpha = np.array([1., -1., 0., 0.])
    alpha[2:] = np.linalg.solve(cpp, -c[2:, :2] @ alpha[:2])
    wave = kwant.wave_function(dev.system, energy, params=dev.params())
    nu = np.asarray([sum((dev.density(p) for p in wave(a)),
                        np.zeros(len(positions))) for a in range(4)]) / (2*np.pi)
    total = nu.sum(axis=0)
    if np.any(total <= 0):
        raise RuntimeError('Characteristic potentials undefined at zero LDOS.')
    u = nu / total
    for i in np.flatnonzero(~dev.sample):
        a = 2 + int(positions[i, 1] < 0)
        u[:, i] = 0.
        u[a, i] = 1.
    u_grid = np.zeros((4, dev.width, dev.length))
    ldos = dev.zero.copy()
    for i in np.flatnonzero(dev.sample):
        x, y = positions[i]
        u_grid[:, y, x] = u[:, i]
        ldos[y, x] = total[i]
    field = np.einsum('a,ayx->yx', alpha, u_grid)
    ce = (dev.matrix(energy+step)-dev.matrix(energy-step))/(2*step)
    cv = (dev.matrix(energy, step*field, step*alpha)
          - dev.matrix(energy, -step*field, -step*alpha))/(2*step)
    frozen = .5 * ce @ (alpha**2)
    q = cv @ alpha + frozen
    beta = np.zeros(4)
    beta[2:] = np.linalg.solve(cpp, -q[2:])
    beta_frozen = np.zeros(4)
    beta_frozen[2:] = np.linalg.solve(cpp, -frozen[2:])
    cs = (dev.matrix(energy, step*field)-dev.matrix(energy, -step*field))/(2*step)
    cl = (dev.matrix(energy, voltages=step*alpha)
          - dev.matrix(energy, voltages=-step*alpha))/(2*step)
    beta_internal = np.linalg.solve(cpp, -(cs@alpha)[2:])
    beta_leads = np.linalg.solve(cpp, -(cl@alpha)[2:])
    current_linear = c @ alpha
    current_quadratic = c @ beta + q
    report = dict(energy=float(energy), delta=float(dev.delta),
        width=int(dev.width), length=int(dev.length),
        alpha=alpha.tolist(), beta=beta.tolist(), beta_frozen=beta_frozen.tolist(),
        alpha_H=float(alpha[2]-alpha[3]), kH=float(beta[2]-beta[3]),
        kH_frozen=float(beta_frozen[2]-beta_frozen[3]),
        kH_internal=float(beta_internal[0]-beta_internal[1]),
        kH_electrodes=float(beta_leads[0]-beta_leads[1]),
        current_linear=current_linear.tolist(), current_quadratic=current_quadratic.tolist(),
        probe_linear_error=float(np.max(abs(current_linear[2:]))),
        probe_quadratic_error=float(np.max(abs(current_quadratic[2:]))),
        current_sum_error=float(max(abs(current_linear.sum()), abs(current_quadratic.sum()))),
        characteristic_sum_error=float(np.max(abs(u.sum(axis=0)-1))),
        reciprocity_error=float(np.max(abs(c-c.T))),
        gauge_shift_error=float(np.max(abs(dev.matrix(energy+.017, dev.zero+.017,
                                                     dev.vzero+.017)-c))),
        cpp_condition=float(np.linalg.cond(cpp)), cpp=cpp.tolist())
    return report, u_grid, field, ldos


def terminal_current(dev, energy, voltages, u_grid, frozen=False, order=12):
    nodes, weights = leggauss(order)
    potential = dev.zero if frozen else np.einsum('a,ayx->yx', voltages, u_grid)
    leadshifts = dev.vzero if frozen else voltages
    currents = np.zeros(4)
    for b, vb in enumerate(voltages):
        half = vb/2
        if half == 0:
            continue
        for node, weight in zip(nodes, weights):
            currents += half * weight * dev.matrix(energy+half*(1+node),
                                                   potential, leadshifts)[:, b]
    return currents


def solve_floating(dev, energy, amplitude, u_grid, report, frozen=False):
    va = np.array([amplitude, -amplitude, 0., 0.])
    beta = np.asarray(report['beta_frozen' if frozen else 'beta'])
    guess = np.asarray(report['alpha'])[2:]*amplitude + beta[2:]*amplitude**2
    def residual(probes):
        va[2:] = probes
        return terminal_current(dev, energy, va, u_grid, frozen)[2:]
    sol = root(residual, guess, method='hybr', options={'xtol': 1e-9})
    va[2:] = sol.x
    currents = terminal_current(dev, energy, va, u_grid, frozen)
    if np.max(abs(currents[2:])) > 1e-12:
        raise RuntimeError(f'Floating solve failed: {sol.message}; currents={currents}')
    return va.copy(), currents


def verify(dev, energy, report, u_grid, step=STEP):
    refined = floating_response(dev, energy, step/2)[0]
    v = 5e-5
    vals, currents = [], []
    checks = {}
    for frozen in (False, True):
        vp, ip = solve_floating(dev, energy, v, u_grid, report, frozen)
        vm, im = solve_floating(dev, energy, -v, u_grid, report, frozen)
        k_direct = ((vp[2]-vp[3])+(vm[2]-vm[3]))/(2*v*v)
        key = 'kH_frozen' if frozen else 'kH'
        checks[key+'_direct'] = float(k_direct)
        checks[key+'_finite_bias_error'] = float(abs(k_direct-report[key]))
        vals.extend((vp.tolist(), vm.tolist()))
        currents.extend((ip, im))
        assert abs(k_direct-report[key]) < max(2e-6, .005*abs(report[key])), checks
    checks.update(derivative_step_error=abs(refined['kH']-report['kH']),
        probe_current_error=float(np.max(abs(np.asarray(currents)[:, 2:]))),
        current_conservation_error=float(np.max(abs(np.asarray(currents).sum(axis=1)))),
        terminal_voltages_at_pm_v=vals, check_bias=v)
    assert checks['derivative_step_error'] < max(2e-6, .005*abs(report['kH'])), checks
    assert report['gauge_shift_error'] < 1e-9, report
    assert abs(report['alpha_H']) < 1e-8, report
    return checks


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.time()
    with threadpool_limits(limits=1):
        dev = Device('four_x')
        rep, u_grid, field, ldos = floating_response(dev, EF)
        print('RM floating: '+json.dumps(rep), flush=True)
        validation = verify(dev, EF, rep, u_grid)
        print('RM validation: '+json.dumps(validation), flush=True)
        energies = np.linspace(-.18, .18, 41)
        energy_reports = [floating_response(dev, float(e))[0] for e in energies]
        print('Energy scan completed.', flush=True)
        widths = np.array([2, 4, 6, 8, 10, 12])
        width_reports = [floating_response(Device('four_x', cells=int(n)), EF)[0]
                         for n in widths]
        deltas = np.array([-.2, -.1, 0., .1, .2])
        # A wider derivative step suppresses numerical near-degenerate lead-mode noise.
        delta_reports = [floating_response(Device('four_x', delta=float(d)), EF, 5e-5)[0]
                         for d in deltas]
        assert abs(delta_reports[0]['kH']+delta_reports[-1]['kH']) < 2e-5
        assert abs(delta_reports[2]['kH']) < 2e-5
        contacts = np.array([.05, .1, .15, .3, .6])
        contact_reports = [floating_response(Device('four_x', tc=float(t)), EF)[0]
                           for t in contacts]
        biases = np.array([2.5e-5, 5e-5, 1e-4, 2e-4])
        plus, minus = [], []
        max_probe_error = 0.
        for v in biases:
            vp, ip = solve_floating(dev, EF, v, u_grid, rep)
            vm, im = solve_floating(dev, EF, -v, u_grid, rep)
            plus.append(vp[2]-vp[3]); minus.append(vm[2]-vm[3])
            max_probe_error = max(max_probe_error, float(np.max(abs(np.r_[ip[2:], im[2:]]))))
        np.savez_compressed(OUT/'floating_study.npz',
            energies=energies, kH_full=[r['kH'] for r in energy_reports],
            kH_frozen=[r['kH_frozen'] for r in energy_reports],
            widths=widths, kH_width=[r['kH'] for r in width_reports],
            deltas=deltas, kH_delta=[r['kH'] for r in delta_reports],
            contact_tc=contacts, contact_kH=[r['kH'] for r in contact_reports],
            biases=biases, vH_plus=plus, vH_minus=minus,
            representative_kH=rep['kH'], field=field, ldos=ldos, u_grid=u_grid)
        # Verify the same corrected protocol once in the original QSH model.
        qsh_path = ROOT.parent/'QSH_Nonlinear_Transport'/'code'/'compute_letter_preview.py'
        spec = importlib.util.spec_from_file_location('qsh_preview', qsh_path)
        qsh_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(qsh_module)
        qdev = qsh_module.Device('four_x', .2)
        qrep, qu, qfield, qldos = floating_response(qdev, .4, 1e-5)
        qcheck = verify(qdev, .4, qrep, qu, 1e-5)
        print('QSH floating: '+json.dumps(qrep), flush=True)
        qout = ROOT.parent/'QSH_Nonlinear_Transport'/'results'/'floating_probe'
        qout.mkdir(parents=True, exist_ok=True)
        qmetadata = dict(protocol='VL=+V,VR=-V,IT=IB=0; VH=VT-VB; V_SD=2V',
            units='v=eV/t, i=hI/(e*t), eVH/t=kH*v**2',
            parameters=dict(A=1., t=1., m=-1., delta=.2, EF=.4, tc=.15,
                            length=31, width=12),
            electrostatics='Equilibrium injectivities, local neutrality, first-order potential',
            temperature=0., report=qrep, validation=qcheck)
        (qout/'representative.json').write_text(json.dumps(qmetadata, indent=2), encoding='utf-8')
        np.savez_compressed(qout/'representative.npz', field=qfield, ldos=qldos, u_grid=qu)
        report = dict(protocol='VL=+V,VR=-V,IT=IB=0; VH=VT-VB; V_SD=2V',
            units='v=eV/t2, i=hI/(e*t2), eVH/t2=kH*v**2',
            electrostatics='Equilibrium injectivities, local neutrality, first-order potential',
            temperature=0., representative=rep,
            validation=validation, width_scan=width_reports, delta_scan=delta_reports,
            contact_scan=contact_reports, bias_scan_probe_error=max_probe_error,
            qsh_representative=qrep, qsh_validation=qcheck,
            elapsed_seconds=time.time()-started)
        (OUT/'floating_validation.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        print(f'Completed in {time.time()-started:.1f} s.', flush=True)


if __name__ == '__main__':
    raise SystemExit("Support module. Use code/run_panel_workflow.py --recompute for a selected panel.")
