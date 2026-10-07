"""Exploratory, gauge-consistent transport data for the three-figure Letter.

Canonical +sin(ky) model. True T/B two-terminal device and grounded T/B
four-terminal Hall device are finalized separately. Scalar local-neutrality
screening, zero temperature, no magnetic field. Reduced v=eV/t, i=hI/(e t).
Only U^(1) is retained in finite-bias checks; these verify the quadratic term.
"""
from pathlib import Path
import argparse
import json
import time
import kwant
import numpy as np
from numpy.polynomial.legendre import leggauss
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results'/'letter_preview'
L, W, TC = 31, 12, .15
M, STEP = -1., 2e-6
S0 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1., -1.]).astype(complex)
I4 = np.eye(4)


class Device:
    def __init__(self, geometry, delta, length=L, width=W):
        L, W = int(length), int(width)
        if L < 5 or L % 2 != 1 or W < 4:
            raise ValueError('Use odd length >=5 and width >=4.')
        self.length, self.width = L, W
        self.geometry, self.delta = geometry, delta
        self.weights = np.array([.5, -.5] if geometry == 'two_y' else [.5, -.5, 0, 0])
        self.project = np.array([1., 0.] if geometry == 'two_y' else [0, 0, 1., -1.])
        lat = kwant.lattice.square(norbs=4, name='sample')
        metal = kwant.lattice.square(norbs=4, name='metal')
        h0 = M*np.kron(SZ, S0)+delta*np.kron(SX, S0)
        hops = {(1, 0): (np.kron(SZ, S0)+1j*np.kron(SX, SZ))/2,
                (0, 1): (np.kron(SZ, S0)+1j*np.kron(SY, S0))/2}

        def sample_onsite(site, potential):
            x, y = site.tag
            return h0+potential[y, x]*I4

        def lead_onsite(index, base):
            def onsite(site, voltages):
                return base+voltages[index]*I4
            return onsite

        syst = kwant.Builder()
        for x in range(L):
            for y in range(W):
                syst[lat(x, y)] = sample_onsite
        for d, hop in hops.items():
            syst[kwant.builder.HoppingKind(d, lat)] = hop
        if geometry == 'four_x':
            for index, direction in enumerate((-1, 1)):
                lead = kwant.Builder(kwant.TranslationalSymmetry((direction, 0)))
                for y in range(W):
                    lead[lat(0, y)] = lead_onsite(index, h0)
                for d, hop in hops.items():
                    lead[kwant.builder.HoppingKind(d, lat)] = hop
                syst.attach_lead(lead)
        xc = L//2
        for index, side in enumerate((1, -1), start=0 if geometry == 'two_y' else 2):
            ybuf, ys = (W, W-1) if side == 1 else (-1, 0)
            syst[metal(xc, ybuf)] = lead_onsite(index, np.zeros((4, 4)))
            syst[metal(xc, ybuf), lat(xc, ys)] = TC*I4
            lead = kwant.Builder(kwant.TranslationalSymmetry((0, side)))
            lead[metal(xc, 0)] = lead_onsite(index, np.zeros((4, 4)))
            lead[kwant.builder.HoppingKind((0, 1), metal)] = -I4
            syst.attach_lead(lead)
        self.system = syst.finalized()
        self.pos = np.array([s.pos for s in self.system.sites], dtype=int)
        self.sample = np.array([s.family.name == 'sample' for s in self.system.sites])
        self.nlead = len(self.system.leads)
        self.zero = np.zeros((W, L))
        self.vzero = np.zeros(self.nlead)
        self.density = kwant.operator.Density(self.system)

    def params(self, potential=None, voltages=None):
        return dict(potential=self.zero if potential is None else potential,
                    voltages=self.vzero if voltages is None else voltages)

    def matrix(self, energy, potential=None, voltages=None):
        sm = kwant.smatrix(self.system, energy, params=self.params(potential, voltages))
        trans = np.zeros((self.nlead, self.nlead))
        for a in range(self.nlead):
            for b in range(self.nlead):
                if a != b:
                    trans[a, b] = sm.transmission(a, b)
        return np.diag(trans.sum(axis=0))-trans

    def response(self, energy, step=STEP, pieces=False):
        c = self.matrix(energy)
        wf = kwant.wave_function(self.system, energy, params=self.params())
        nu = np.array([sum((self.density(psi) for psi in wf(a)),
                          np.zeros(len(self.pos))) for a in range(self.nlead)])/(2*np.pi)
        total = nu.sum(axis=0)
        if np.any(total <= 0):
            raise RuntimeError('Zero LDOS; characteristic potential undefined.')
        u = nu/total
        for i in np.flatnonzero(~self.sample):
            a = (0 if self.geometry == 'two_y' else 2)+int(self.pos[i, 1] < 0)
            u[:, i] = 0
            u[a, i] = 1
        field = np.zeros((self.width, self.length))
        ldos = np.zeros((self.width, self.length))
        for i in np.flatnonzero(self.sample):
            x, y = self.pos[i]
            field[y, x] = self.weights@u[:, i]
            ldos[y, x] = total[i]
        ce = (self.matrix(energy+step)-self.matrix(energy-step))/(2*step)
        cv = (self.matrix(energy, step*field, step*self.weights)
              -self.matrix(energy, -step*field, -step*self.weights))/(2*step)
        frozen = .5*ce@(self.weights**2)
        full = cv@self.weights+frozen
        g1 = c@self.weights
        report = dict(energy=float(energy), delta=float(self.delta), geometry=self.geometry,
            g1=float(self.project@g1), q=float(self.project@full),
            q_frozen=float(self.project@frozen),
            currents_linear=g1.tolist(), currents_quadratic=full.tolist(),
            gauge_sum_error=float(np.max(np.abs(u.sum(axis=0)-1))),
            reciprocity_error=float(np.max(np.abs(c-c.T))),
            row_sum_error=float(np.max(np.abs(c.sum(axis=1)))),
            nonlinear_current_sum_error=float(abs(full.sum())),
            ldos_min=float(ldos.min()), ldos_max=float(ldos.max()))
        if pieces:
            cs = (self.matrix(energy, step*field)-self.matrix(energy, -step*field))/(2*step)
            cl = (self.matrix(energy, voltages=step*self.weights)
                  -self.matrix(energy, voltages=-step*self.weights))/(2*step)
            report.update(q_internal=float(self.project@cs@self.weights),
                          q_electrodes=float(self.project@cl@self.weights),
                          split_error=float(np.max(np.abs(cv-cs-cl))))
            offset = .017
            shifted = self.matrix(energy+offset, self.zero+offset, self.vzero+offset)
            report['global_gauge_shift_error'] = float(np.max(np.abs(c-shifted)))
        return dict(report=report, field=field, u=u, injectivity=nu, ldos=ldos)

    def current(self, energy, voltage, field, frozen=False, order=8):
        nodes, weights = leggauss(order)
        current = np.zeros(self.nlead)
        potential = self.zero if frozen else voltage*field
        voltages = self.vzero if frozen else voltage*self.weights
        for b, w in enumerate(self.weights):
            if w == 0:
                continue
            half = voltage*w/2
            for node, weight in zip(nodes, weights):
                e = energy+half*(1+node)
                current += half*weight*self.matrix(e, potential, voltages)[:, b]
        return current


def main(test=False):
    OUT.mkdir(parents=True, exist_ok=True)
    start = time.time()
    devices = {}
    cache = {}
    def get(geometry, delta, energy, pieces=False):
        dk = (geometry, round(float(delta), 8))
        if dk not in devices:
            devices[dk] = Device(*dk)
        key = (*dk, round(float(energy), 8), pieces)
        if key not in cache:
            cache[key] = devices[dk].response(energy, pieces=pieces)
        return cache[key]
    with threadpool_limits(limits=1):
        reps = [get(g, .2, .4, pieces=True) for g in ('two_y', 'four_x')]
        for rep in reps:
            print(json.dumps(rep['report']), flush=True)
        if test:
            print(f'Timing for two protocols: {time.time()-start:.2f}s', flush=True)
            return
        energies = np.linspace(.10, .70, 41)
        scans = []
        for g in ('two_y', 'four_x'):
            rows = []
            for i, energy in enumerate(energies):
                rows.append(get(g, .2, energy)['report'])
                if i % 10 == 0:
                    print(f'{g}: {i+1}/{len(energies)} energies', flush=True)
            scans.append(rows)
        deltas = np.linspace(-.3, .3, 13)
        reversals = [[get(g, d, .4)['report'] for d in deltas] for g in ('two_y', 'four_x')]
        validations = []
        voltage_points = np.array([.00005, .0001, .0002, .0004, .0008])
        even, odd, frozen_even = [], [], []
        for g, rep in zip(('two_y', 'four_x'), reps):
            dev = devices[(g, .2)]
            refined = dev.response(.4, step=STEP/2)['report']
            positive = np.array([dev.current(.4, v, rep['field']) for v in voltage_points])
            negative = np.array([dev.current(.4, -v, rep['field']) for v in voltage_points])
            fp = np.array([dev.current(.4, v, rep['field'], frozen=True) for v in voltage_points])
            fm = np.array([dev.current(.4, -v, rep['field'], frozen=True) for v in voltage_points])
            ev = (positive+negative)@dev.project/2
            od = (positive-negative)@dev.project/2
            fe = (fp+fm)@dev.project/2
            even.append(ev); odd.append(od); frozen_even.append(fe)
            ctrl = get(g, 0., .4)['report']
            mirror = get(g, -.2, .4)['report']
            q = rep['report']['q']
            scale = max(abs(q), 1e-8)
            check = dict(geometry=g, q=q,
                derivative_step_relative_error=abs(refined['q']-q)/scale,
                finite_bias_relative_error_smallest=abs(ev[0]/voltage_points[0]**2-q)/scale,
                delta_reversal_relative_error=abs(mirror['q']+q)/scale,
                symmetric_delta_zero_q=ctrl['q'],
                frozen_even_smallest=float(fe[0]),
                finite_bias_current_sum_error=float(max(np.max(abs(positive.sum(axis=1))),
                                                        np.max(abs(negative.sum(axis=1))))))
            assert check['derivative_step_relative_error'] < .005, check
            assert check['finite_bias_relative_error_smallest'] < .02, check
            assert check['delta_reversal_relative_error'] < 1e-4, check
            assert abs(ctrl['q']) < 1e-5, check
            assert rep['report']['global_gauge_shift_error'] < 1e-9
            validations.append(check)
            print(json.dumps(check), flush=True)
        records = [v['report'] for v in cache.values()]
        assert max(r['gauge_sum_error'] for r in records) < 1e-12
        assert max(r['reciprocity_error'] for r in records) < 1e-9
        assert max(r['nonlinear_current_sum_error'] for r in records) < 1e-6
        np.savez_compressed(OUT/'preview_data.npz', energies=energies, deltas=deltas,
            g1=np.array([[r['g1'] for r in rows] for rows in scans]),
            q=np.array([[r['q'] for r in rows] for rows in scans]),
            q_frozen=np.array([[r['q_frozen'] for r in rows] for rows in scans]),
            q_delta=np.array([[r['q'] for r in rows] for rows in reversals]),
            g1_delta=np.array([[r['g1'] for r in rows] for rows in reversals]),
            representative_fields=np.array([r['field'] for r in reps]),
            representative_ldos=np.array([r['ldos'] for r in reps]),
            voltage=voltage_points, even=even, odd=odd, frozen_even=frozen_even,
            q_at_representative=np.array([r['report']['q'] for r in reps]))
        report = dict(parameters=dict(L=L, W=W, m=M, A=1., t=1., tc=TC, magnetic_field=0,
                                      derivative_step=STEP, reference_energy=.4),
            geometry_order=['two_y', 'four_x'],
            approximation='zero-temperature local charge neutrality; finite bias uses U^(1) only',
            scan_status='Exploratory 41-point energy grid; no resonance-resolution or size convergence claim',
            representative=[r['report'] for r in reps], checks=validations,
            scans=scans, reversals=reversals, runtime_seconds=time.time()-start)
        (OUT/'validation.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        print(f'Saved preview data in {time.time()-start:.1f}s', flush=True)


if __name__ == '__main__':
    raise SystemExit("Support module. Use code/run_panel_workflow.py --recompute for a selected panel.")
