"""TRS QSH response with physical floating contacts and elastic spin probes.

Basis per site: (orbital 1 up, orbital 1 down, orbital 2 up, orbital 2 down).
Physical terminals are L,R,T,B. H=H0+U, mu_a=EF+v_a, v=eV/t.
Virtual WBA reservoirs are paired up/down, have separate energy-resolved
zero-current occupations, and do not relax the conserved spin population.
Physical T/B are Fermi reservoirs with integrated zero current. No field,
phonons, magnetic disorder, imaginary eta, or nonlinear Poisson solver.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'

from functools import lru_cache
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.special import expit
from scipy.optimize import root
from numpy.polynomial.legendre import leggauss
import kwant
from compute_letter_preview import Device, I4, S0, SY


def current_matrix(transmission):
    """T[out,in]; currents positive from the reservoirs into the sample."""
    return np.diag(transmission.sum(axis=0))-transmission


def quadrature(lo, hi, panel_width=.02, order=8):
    nodes, weights = leggauss(order)
    edges = np.linspace(lo, hi, max(1, int(np.ceil((hi-lo)/panel_width)))+1)
    energies, result_weights = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        energies.extend((a+b)/2+(b-a)*nodes/2)
        result_weights.extend((b-a)*weights/2)
    return np.array(energies), np.array(result_weights)


def fermi_weight(energies, mu, temperature):
    if temperature <= 0:
        raise ValueError('T=0 uses an on-shell calculation.')
    f = expit((mu-np.asarray(energies))/temperature)
    w = f*(1-f)/temperature
    return w, -w*(1-2*f)/temperature


class EnvironmentDevice:
    def __init__(self, length=31, width=12, delta=.2, disorder=0., seed=0,
                 gamma_phi=0., probe_stride=5, mirror_disorder=False):
        self.base = Device('four_x', delta, length, width)
        self.length, self.width, self.delta = length, width, float(delta)
        self.gamma_phi = float(gamma_phi)
        self.positions = self.base.pos.copy()
        self.sample = self.base.sample.copy()
        self.nsites = len(self.positions)
        self.n = 4*self.nsites
        self.h = self.base.system.hamiltonian_submatrix(
            params=self.base.params(), sparse=True).tocsc()
        draw = np.random.default_rng(seed).uniform(-.5, .5, (width, length))
        if mirror_disorder:
            for x in range((length+1)//2):
                draw[:, length-1-x] = draw[:, x]
        self.disorder_grid = float(disorder)*draw
        onsite = np.zeros(self.nsites)
        for i in np.flatnonzero(self.sample):
            x, y = self.positions[i]
            onsite[i] = self.disorder_grid[y, x]
        self.h += sp.diags(np.repeat(onsite, 4), format='csc')
        self.groups = [(4*np.asarray(ids)[:, None]+np.arange(4)).ravel()
                       for ids in self.base.system.lead_interfaces]
        self.probe_positions = []
        self.virtual_spins = []
        if self.gamma_phi > 0:
            xc = length//2
            xs = sorted({x for a in range(probe_stride, xc, probe_stride)
                         for x in (a, length-1-a)})
            ys = sorted({0, 1, width-2, width-1})
            sites = {tuple(p): i for i, p in enumerate(self.positions) if self.sample[i]}
            for x in xs:
                for y in ys:
                    i = sites[(x, y)]
                    self.probe_positions.append((x, y))
                    for spin in (0, 1):
                        self.groups.append(4*i+np.array([spin, 2+spin]))
                        self.virtual_spins.append(spin)
        self.nterm = len(self.groups)
        self.terminal_swap = np.r_[np.arange(4),
            np.arange(4, self.nterm).reshape(-1, 2)[:, ::-1].ravel()]
        tr = sp.kron(sp.eye(self.nsites, format='csc'),
                     np.kron(S0, 1j*SY), format='csc')
        residual = self.h-tr@self.h.conjugate()@tr.getH()
        self.trs_h_error = float(np.max(abs(residual.data))) if residual.nnz else 0.
        if self.trs_h_error > 1e-12:
            raise RuntimeError('Hamiltonian is not spinful TRS.')

    @lru_cache(maxsize=8192)
    def physical_sigma(self, lead, energy):
        # Shifting lead onsite by v is equivalent to Sigma(E-v).
        return self.base.system.leads[lead].selfenergy(
            float(energy), params=self.base.params())

    def _solve(self, energy, potential=None, voltages=None, density=False):
        va = np.zeros(4) if voltages is None else np.asarray(voltages)
        sitepotential = np.zeros(self.nsites) if potential is None else np.asarray(potential)
        if sitepotential.shape != (self.nsites,):
            raise ValueError('Potential must contain one scalar per finalized site.')
        matrix = (sp.diags(np.full(self.n, energy, dtype=complex)
                          -np.repeat(sitepotential, 4), format='csc')-self.h).tolil()
        roots, channel_terminal = [], []
        for a, ids in enumerate(self.groups):
            sigma = (self.physical_sigma(a, float(energy-va[a])) if a < 4
                     else -.5j*self.gamma_phi*np.eye(2))
            matrix[np.ix_(ids, ids)] -= sigma
            linewidth = 1j*(sigma-sigma.conj().T)
            ev, vec = np.linalg.eigh(linewidth)
            if ev.min() < -2e-8:
                raise RuntimeError('Noncausal retarded self energy.')
            keep = ev > 1e-10
            roots.append((ids, vec[:, keep]*np.sqrt(ev[keep])))
            channel_terminal.extend([a]*np.count_nonzero(keep))
        channels = np.asarray(channel_terminal, int)
        b = np.zeros((self.n, len(channels)), complex)
        offset = 0
        for ids, factor in roots:
            b[np.ix_(ids, np.arange(offset, offset+factor.shape[1]))] = factor
            offset += factor.shape[1]
        lu = splu(matrix.tocsc())
        x = lu.solve(b)
        amplitudes = b.conj().T@x
        members = np.zeros((len(channels), self.nterm))
        members[np.arange(len(channels)), channels] = 1.
        transmission = members.T@(abs(amplitudes)**2)@members
        np.fill_diagonal(transmission, 0.)
        c = current_matrix(transmission)
        w = np.eye(self.nterm, 4)
        if self.nterm > 4:
            w[4:] = -np.linalg.solve(c[4:, 4:], c[4:, :4])
        effective = c[:4]@w
        checks = dict(current_conservation=float(np.max(abs(c.sum(axis=0)))),
                      gauge=float(np.max(abs(c.sum(axis=1)))),
                      effective_conservation=float(np.max(abs(effective.sum(axis=0)))),
                      effective_gauge=float(np.max(abs(effective.sum(axis=1)))),
                      occupation_partition=float(np.max(abs(w.sum(axis=1)-1))),
                      virtual_probe_residual=float(np.max(abs((c@w)[4:]))) if self.nterm > 4 else 0.,
                      zero_bias_current=float(np.max(abs(effective@np.ones(4)))))
        if max(checks.values()) > 5e-8:
            raise RuntimeError(checks)
        result = dict(C=effective, W=w, C_all=c, checks=checks)
        if density:
            rho_all = ((abs(x)**2)@members).reshape(self.nsites, 4, self.nterm).sum(axis=1)/(2*np.pi)
            rho = rho_all@w
            result.update(rho=rho, x=x, y=lu.solve(b.conj(), trans='T').T,
                          amplitudes=amplitudes, channels=channels)
            result['checks']['injectivity_partition'] = float(np.max(
                abs(rho.sum(axis=1)-rho_all.sum(axis=1))))
            result['checks']['reciprocity'] = float(np.max(abs(
                c-c[np.ix_(self.terminal_swap, self.terminal_swap)].T)))
        return result

    def effective(self, energy, potential=None, voltages=None):
        return self._solve(energy, potential, voltages)['C']

    def sensitivity(self, data, alpha):
        """d[C_eff alpha]/d onsite U_r, re-eliminating virtual occupations."""
        c, w = data['C_all'], data['W']
        drive = w@alpha
        ell = np.eye(self.nterm, 4)
        if self.nterm > 4:
            ell[4:] = -np.linalg.solve(c[4:, 4:].T, c[:4, 4:].T)
        channel = data['channels']
        x, y, amplitude = data['x'], data['y'], data['amplitudes']
        kernels = np.empty((4, self.nsites))
        for a in range(4):
            z = (ell[channel, a, None]*(drive[channel, None]-drive[channel][None, :])
                 *amplitude.conj())
            orbital = 2*np.einsum('ir,ir->r', y, z@x.T, optimize=True).real
            kernels[a] = orbital.reshape(self.nsites, 4).sum(axis=1)
        return kernels

    def characteristic(self, rho):
        total = rho.sum(axis=1)
        if np.any(total <= 0):
            raise RuntimeError('Characteristic potential undefined at zero LDOS.')
        u = rho/total[:, None]
        for i in np.flatnonzero(~self.sample):
            u[i] = 0.
            u[i, 2+int(self.positions[i, 1] < 0)] = 1.
        return u

    def point(self, energy, alpha=None, step=1e-5):
        data = self._solve(energy, density=True)
        c = data['C']
        if alpha is None:
            alpha = np.array([1., -1., 0., 0.])
            alpha[2:] = np.linalg.solve(c[2:, 2:], -c[2:, :2]@alpha[:2])
        kernel = self.sensitivity(data, alpha)
        leadcv = (self.effective(energy, voltages=step*alpha)
                  -self.effective(energy, voltages=-step*alpha))/(2*step)
        ce = (self.effective(energy+step)-self.effective(energy-step))/(2*step)
        data.update(alpha=alpha, D=kernel, qlead=leadcv@alpha,
                    CE=ce, qwindow=.5*ce@(alpha**2))
        return data

    def response_zero(self, energy, step=1e-5):
        data = self.point(energy, step=step)
        u = self.characteristic(data['rho'])
        field = u@data['alpha']
        qint = data['D']@field
        report = reduce_floating(data['C'], data['alpha'], qint,
                                 data['qlead'], data['qwindow'])
        report.update(energy=float(energy), temperature=0., checks=data['checks'],
                      trs_h_error=self.trs_h_error)
        return report, dict(C=data['C'], rho=data['rho'], D=data['D'], CE=data['CE'],
                            qlead=data['qlead'], u=u, field=field,
                            disorder=self.disorder_grid,
                            probe_positions=np.asarray(self.probe_positions, int).reshape(-1, 2))

    def spectral_record(self, energy, step=1e-5):
        # Clean Mx geometry: alpha_P is zero at every E and after thermal averaging.
        p = self.point(energy, alpha=np.array([1., -1., 0., 0.]), step=step)
        return {k:p[k] for k in ('C', 'rho', 'D', 'qlead', 'qwindow', 'checks')}

    def terminal_current(self, energy, voltages, u, temperature=0., nodes=None, weights=None):
        va = np.asarray(voltages)
        field = u@va
        result = np.zeros(4)
        if temperature == 0:
            z, ws = leggauss(12)
            for b, vb in enumerate(va):
                if vb == 0:
                    continue
                for node, weight in zip(z, ws):
                    en = energy+vb*(1+node)/2
                    result += vb*weight/2*self.effective(en, field, va)[:, b]
        else:
            equilibrium = expit((energy-nodes)/temperature)
            for en, weight, f0 in zip(nodes, weights, equilibrium):
                df = expit((energy+va-en)/temperature)-f0
                result += weight*(self.effective(en, field, va)@df)
        return result

    def solve_floating(self, energy, amplitude, report, u, temperature=0., nodes=None, weights=None):
        va = np.array([amplitude, -amplitude, 0., 0.])
        guess = np.asarray(report['alpha'])[2:]*amplitude+np.asarray(report['beta'])[2:]*amplitude**2
        def residual(probes):
            va[2:] = probes
            return self.terminal_current(energy, va, u, temperature, nodes, weights)[2:]
        sol = root(residual, guess, method='hybr', options={'xtol':1e-9})
        va[2:] = sol.x
        currents = self.terminal_current(energy, va, u, temperature, nodes, weights)
        if np.max(abs(currents[2:])) > 2e-11:
            raise RuntimeError(f'Floating solve failed: {sol.message}; I={currents}')
        return va.copy(), currents


def reduce_floating(c, alpha, qint, qlead, qwindow):
    sources = np.array([qint, qlead, qwindow])
    q = sources.sum(axis=0)
    beta = np.zeros(4)
    beta[2:] = np.linalg.solve(c[2:, 2:], -q[2:])
    parts = -np.linalg.solve(c[2:, 2:], sources[:, 2:].T).T
    i1, i2 = c@alpha, c@beta+q
    return dict(kH=float(beta[2]-beta[3]), alpha_H=float(alpha[2]-alpha[3]),
        alpha=alpha.tolist(), beta=beta.tolist(), kH_internal=float(parts[0, 0]-parts[0, 1]),
        kH_electrodes=float(parts[1, 0]-parts[1, 1]),
        kH_frozen=float(parts[2, 0]-parts[2, 1]),
        Gxx=float(i1[0]/2), Gxx_bias_amplitude=float(i1[0]),
        cpp_condition=float(np.linalg.cond(c[2:, 2:])),
        probe_linear_error=float(np.max(abs(i1[2:]))),
        probe_quadratic_error=float(np.max(abs(i2[2:]))),
        current_sum_error=float(max(abs(i1.sum()), abs(i2.sum()))),
        current_linear=i1.tolist(), current_quadratic=i2.tolist())


def thermal_response(device, energies, weights, spectrum, mu, temperature):
    w, wp = fermi_weight(energies, mu, temperature)
    avg = lambda a: np.tensordot(weights*w, a, axes=(0, 0))
    c, rho, d, qlead = (avg(spectrum[k]) for k in ('C', 'rho', 'D', 'qlead'))
    alpha = np.array([1., -1., 0., 0.])
    alpha[2:] = np.linalg.solve(c[2:, 2:], -c[2:, :2]@alpha[:2])
    if max(abs(alpha[2:])) > 2e-8:
        raise ValueError('This cached thermal kernel is for clean Mx geometries only.')
    u = device.characteristic(rho)
    qwindow = -.5*np.einsum('e,eab,b->a', weights*wp, spectrum['C'], alpha**2)
    qint = d@(u@alpha)
    report = reduce_floating(c, alpha, qint, qlead, qwindow)
    report.update(energy=float(mu), temperature=float(temperature),
                  thermal_mass=float(np.sum(weights*w)),
                  window_ibp_error=float(np.max(abs(qwindow-avg(spectrum['qwindow'])))))
    return report, u
