"""Stacked Rice-Mele boundary metal: an initial zero-field transport study.

Each physical site has one spinless orbital. Even/odd rows are A/B;
N cells give 2N rows. Bottom terminates on A (+Delta), top on B (-Delta).
Currents are positive from reservoirs into sample, mu_a=EF+v_a.
Reduced energy uses t2=1, i=hI/(e t2), v=eV/t2; q multiplies v**2.
Spin degeneracy, when applicable, multiplies all terminal currents by two.
Local neutrality is an explicitly approximate electrostatic closure.
"""
from pathlib import Path
import json
import time
import kwant
import numpy as np
from numpy.polynomial.legendre import leggauss
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
TX, T1, T2, DELTA, EF, TC = .25, .3, 1., .2, .1, .15
LENGTH, NREF, STEP = 31, 8, 1e-5


class Device:
    def __init__(self, geometry, cells=NREF, delta=DELTA, tc=TC,
                 tx=TX, t1=T1, t2=T2, length=LENGTH):
        self.geometry = geometry
        self.width, self.length = 2 * int(cells), int(length)
        self.delta, self.tc = float(delta), float(tc)
        self.weights = np.array([.5, -.5] if geometry == 'two_y'
                                else [.5, -.5, 0., 0.])
        self.project = np.array([1., 0.] if geometry == 'two_y'
                                else [0., 0., 1., -1.])
        lat = kwant.lattice.square(norbs=1, name='rm_sample')
        metal = kwant.lattice.square(norbs=1, name='probe_metal')

        def onsite(site, potential):
            x, y = site.tag
            return delta * (-1)**y + potential[y, x]

        def vertical(site_to, site_from):
            lower = min(site_to.tag[1], site_from.tag[1])
            return t1 if lower % 2 == 0 else t2

        def shifted(index, sample):
            def lead_onsite(site, voltages):
                base = delta * (-1)**site.tag[1] if sample else 0.
                return base + voltages[index]
            return lead_onsite

        builder = kwant.Builder()
        for x in range(self.length):
            for y in range(self.width):
                builder[lat(x, y)] = onsite
        builder[kwant.builder.HoppingKind((1, 0), lat)] = -tx
        builder[kwant.builder.HoppingKind((0, 1), lat)] = vertical
        if geometry == 'four_x':
            for index, direction in enumerate((-1, 1)):
                lead = kwant.Builder(kwant.TranslationalSymmetry((direction, 0)))
                for y in range(self.width):
                    lead[lat(0, y)] = shifted(index, True)
                lead[kwant.builder.HoppingKind((1, 0), lat)] = -tx
                lead[kwant.builder.HoppingKind((0, 1), lat)] = vertical
                builder.attach_lead(lead)
        xc = self.length // 2
        first = 0 if geometry == 'two_y' else 2
        for index, side in enumerate((1, -1), start=first):
            ybuf, edge = (self.width, self.width-1) if side == 1 else (-1, 0)
            builder[metal(xc, ybuf)] = shifted(index, False)
            builder[metal(xc, ybuf), lat(xc, edge)] = tc
            lead = kwant.Builder(kwant.TranslationalSymmetry((0, side)))
            lead[metal(xc, 0)] = shifted(index, False)
            lead[kwant.builder.HoppingKind((0, 1), metal)] = -1.
            builder.attach_lead(lead)
        self.system = builder.finalized()
        self.positions = np.asarray([s.pos for s in self.system.sites], int)
        self.sample = np.asarray([s.family == lat for s in self.system.sites])
        self.nlead = len(self.system.leads)
        self.zero = np.zeros((self.width, self.length))
        self.vzero = np.zeros(self.nlead)
        self.density = kwant.operator.Density(self.system)

    def params(self, potential=None, voltages=None):
        return dict(potential=self.zero if potential is None else potential,
                    voltages=self.vzero if voltages is None else voltages)

    def matrix(self, energy, potential=None, voltages=None):
        sm = kwant.smatrix(self.system, energy,
                          params=self.params(potential, voltages))
        trans = np.zeros((self.nlead, self.nlead))
        for a in range(self.nlead):
            for b in range(self.nlead):
                if a != b:
                    trans[a, b] = sm.transmission(a, b)
        return np.diag(trans.sum(axis=0))-trans

    def response(self, energy, step=STEP, pieces=False):
        c = self.matrix(energy)
        wave = kwant.wave_function(self.system, energy, params=self.params())
        nu = np.asarray([sum((self.density(p) for p in wave(a)),
                            np.zeros(len(self.positions)))
                         for a in range(self.nlead)]) / (2*np.pi)
        total = nu.sum(axis=0)
        if np.any(total <= 0):
            raise RuntimeError('Characteristic potential undefined at zero LDOS.')
        u = nu/total
        for i in np.flatnonzero(~self.sample):
            index = (0 if self.geometry == 'two_y' else 2)
            index += int(self.positions[i, 1] < 0)
            u[:, i] = 0.
            u[index, i] = 1.
        field, ldos = self.zero.copy(), self.zero.copy()
        for i in np.flatnonzero(self.sample):
            x, y = self.positions[i]
            field[y, x] = self.weights @ u[:, i]
            ldos[y, x] = total[i]
        ce = (self.matrix(energy+step)-self.matrix(energy-step))/(2*step)
        cv = (self.matrix(energy, step*field, step*self.weights)
              -self.matrix(energy, -step*field, -step*self.weights))/(2*step)
        frozen = .5 * ce @ (self.weights**2)
        quadratic = cv @ self.weights + frozen
        linear = c @ self.weights
        report = dict(geometry=self.geometry, cells=self.width//2,
                      energy=float(energy), delta=self.delta, tc=self.tc,
                      g1=float(self.project @ linear),
                      q=float(self.project @ quadratic),
                      q_frozen=float(self.project @ frozen),
                      currents_linear=linear.tolist(),
                      currents_quadratic=quadratic.tolist(),
                      reciprocity_error=float(np.max(abs(c-c.T))),
                      linear_sum_error=float(abs(linear.sum())),
                      quadratic_sum_error=float(abs(quadratic.sum())),
                      characteristic_sum_error=float(np.max(abs(u.sum(axis=0)-1))),
                      min_sample_ldos=float(ldos.min()))
        if pieces:
            cs = (self.matrix(energy, step*field)-self.matrix(energy, -step*field))/(2*step)
            cl = (self.matrix(energy, voltages=step*self.weights)
                  -self.matrix(energy, voltages=-step*self.weights))/(2*step)
            report.update(q_internal=float(self.project @ cs @ self.weights),
                          q_electrodes=float(self.project @ cl @ self.weights),
                          split_error=float(np.max(abs(cv-cs-cl))),
                          gauge_shift_error=float(np.max(abs(
                              self.matrix(energy+.017, self.zero+.017, self.vzero+.017)-c))))
        return report, field, ldos

    def current(self, energy, voltage, field, frozen=False, order=12):
        nodes, weights = leggauss(order)
        result = np.zeros(self.nlead)
        potential = self.zero if frozen else voltage*field
        voltages = self.vzero if frozen else voltage*self.weights
        for b, w in enumerate(self.weights):
            half = voltage*w/2
            if half == 0:
                continue
            for node, weight in zip(nodes, weights):
                result += half*weight*self.matrix(energy+half*(1+node),
                                                  potential, voltages)[:, b]
        return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.time()
    records, checks = [], []
    with threadpool_limits(limits=1):
        devices = {g: Device(g) for g in ('two_y', 'four_x')}
        representative = {}
        for g, dev in devices.items():
            rep, field, ldos = dev.response(EF, pieces=True)
            representative[g] = rep
            records.append(rep)
            print(json.dumps(rep), flush=True)
            refined = dev.response(EF, step=STEP/2)[0]
            v = 2e-5
            plus = dev.current(EF, v, field)
            minus = dev.current(EF, -v, field)
            q_direct = float(dev.project @ (plus+minus)/(2*v*v))
            fp = dev.current(EF, v, field, frozen=True)
            fm = dev.current(EF, -v, field, frozen=True)
            qf_direct = float(dev.project @ (fp+fm)/(2*v*v))
            check = dict(geometry=g, q=rep['q'], q_direct=q_direct,
                         derivative_step_error=abs(rep['q']-refined['q']),
                         finite_bias_error=abs(rep['q']-q_direct),
                         frozen_finite_bias_error=abs(rep['q_frozen']-qf_direct),
                         current_sum_error=float(max(abs(plus.sum()), abs(minus.sum()))))
            assert check['finite_bias_error'] < max(2e-7, .005*abs(rep['q'])), check
            assert check['derivative_step_error'] < max(1e-7, .002*abs(rep['q'])), check
            assert rep['gauge_shift_error'] < 1e-9, rep
            checks.append(check)
            np.savez_compressed(OUT/f'{g}_representative.npz', field=field, ldos=ldos)
        widths = np.array([2, 4, 6, 8, 10, 12])
        yw, hw = [], []
        for n in widths:
            y = Device('two_y', cells=n).response(EF)[0]
            h = Device('four_x', cells=n).response(EF)[0]
            yw.append(y); hw.append(h); records.extend((y, h))
            print(f'width cells={n}: Gy={y["g1"]:.8g}, qy={y["q"]:.8g}, qH={h["q"]:.8g}', flush=True)
        energies = np.linspace(-.18, .18, 41)
        hs = []
        for energy in energies:
            r = devices['four_x'].response(float(energy))[0]
            hs.append(r); records.append(r)
        deltas = np.array([-.2, -.1, 0., .1, .2])
        dh, dy = [], []
        for delta in deltas:
            h = Device('four_x', delta=delta).response(EF)[0]
            y = Device('two_y', delta=delta).response(EF)[0]
            dh.append(h); dy.append(y); records.extend((h,y))
        assert abs(dh[0]['q']+dh[-1]['q']) < 1e-7
        assert abs(dh[2]['q']) < 1e-7
        assert abs(dy[0]['q']+dy[-1]['q']) < 1e-7
        assert abs(dy[2]['q']) < 1e-7
        # Fixed bulk spectrum with the opposite dimer termination: no propagating
        # L/R strip channels in this gap. This is a boundary-source control only.
        trivial = Device('four_x', t1=T2, t2=T1).response(EF)[0]
        coupling = []
        for tc in (.05, .1, .15, .3, .6):
            r = Device('four_x', tc=tc).response(EF)[0]
            coupling.append(r); records.append(r)
        # Infinite-width edge projection: exact boundary weight, a side-coupled
        # one-channel metal self-energy, matched infinite x chain.
        def projected(energy, center):
            velocity = np.sqrt(4*TX**2-(energy-center)**2)
            sigma = (1-(T1/T2)**2)*TC**2*(energy-1j*np.sqrt(4-energy**2))/2
            gamma = -2*sigma.imag
            return gamma/velocity / abs(1+1j*sigma/velocity)**2
        def projected_q(energy):
            def difference(e):
                return projected(e, -DELTA)-projected(e, DELTA)
            return -(difference(energy+STEP)-difference(energy-STEP))/(8*STEP)
        edge_model_q = np.asarray([projected_q(e) for e in energies])
        # Strip bands are a tensor-sum exactly; identify the two y-edge modes.
        strip_cells = 24
        rows = 2*strip_cells
        hy = np.diag([DELTA*(-1)**y for y in range(rows)])
        for y in range(rows-1):
            hy[y,y+1] = hy[y+1,y] = T1 if y%2 == 0 else T2
        values, vectors = np.linalg.eigh(hy)
        itop, ibottom = int(np.argmin(abs(values+DELTA))), int(np.argmin(abs(values-DELTA)))
        labels = np.zeros(rows, int); labels[itop] = -1; labels[ibottom] = 1
        kx = np.linspace(-np.pi, np.pi, 401)
        np.savez_compressed(OUT/'initial_study.npz', kx=kx,
            strip_bands=values[None,:]-2*TX*np.cos(kx[:,None]), edge_labels=labels,
            analytic_top=-2*TX*np.cos(kx)-DELTA, analytic_bottom=-2*TX*np.cos(kx)+DELTA,
            edge_profile_top=abs(vectors[:,itop])**2,
            edge_profile_bottom=abs(vectors[:,ibottom])**2,
            widths=widths, g1_y=[r['g1'] for r in yw], q_y=[r['q'] for r in yw],
            qfrozen_y=[r['q_frozen'] for r in yw], qH_width=[r['q'] for r in hw],
            energies=energies, qH_full=[r['q'] for r in hs],
            qH_frozen=[r['q_frozen'] for r in hs], qH_projected=edge_model_q,
            deltas=deltas, qH_delta=[r['q'] for r in dh], q_y_delta=[r['q'] for r in dy],
            contact_tc=[r['tc'] for r in coupling], contact_qH=[r['q'] for r in coupling])
        rmin = float(np.sqrt((T2-T1)**2+DELTA**2))
        validation = dict(parameters=dict(tx=TX,t1=T1,t2=T2,delta=DELTA,EF=EF,
                            tc=TC,L=LENGTH,Nref=NREF,spinless=True,magnetic_field=0),
            approximation='zero temperature; local neutrality; first-order potential; no full Poisson',
            scope='Initial fixed-contact study, not a general size/contact/temperature convergence claim',
            bulk_gap_halfwidth=rmin-2*TX,
            density_decay_per_cell=(T1/T2)**2,
            eigenedge_energy_error=float(max(abs(values[itop]+DELTA),abs(values[ibottom]-DELTA))),
            representative=representative, checks=checks, width_two_y=yw, width_four_x=hw,
            delta_four_x=dh, delta_two_y=dy, termination_control=trivial, contact_scan=coupling,
            edge_projection_error_at_reference=float(abs(projected_q(EF)-representative['four_x']['q_frozen'])),
            max_reciprocity_error=max(r['reciprocity_error'] for r in records),
            max_quadratic_sum_error=max(r['quadratic_sum_error'] for r in records),
            elapsed_seconds=time.time()-started)
        assert validation['max_reciprocity_error'] < 1e-9
        assert validation['max_quadratic_sum_error'] < 1e-6
        (OUT/'validation.json').write_text(json.dumps(validation, indent=2), encoding='utf-8')
        print(json.dumps(dict(summary=validation['representative'],
                             edge_projection_error=validation['edge_projection_error_at_reference'],
                             elapsed=validation['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    raise SystemExit("Support module. Use code/run_panel_workflow.py --recompute for a selected panel.")
