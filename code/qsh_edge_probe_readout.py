"""QSH bare-edge quantities and symmetric floating-probe readout.

The Hamiltonian and source/probe geometry match ``compute_letter_preview``.
Bare means that both transverse probes and their buffer sites are absent.
Four-terminal probes are identical: tc3=tc4=tc.  Their scalar orbital
linewidth is Gamma=tc**2*sqrt(4-E**2), multiplying the 4x4 identity.

The returned kappa2 is t*kappa2_physical/e for mu=EF-eV and source
voltages +/-V_SD/2.  The historical floating solver uses onsite-positive
amplitudes +/-v; the conversion kappa2=-kH_code/4 occurs exactly once.
Contact densities sum all four spin/orbital components of incoming,
unit-flux states, divided by 2*pi.  Zero field, disorder and dephasing.
"""
from dataclasses import asdict, dataclass

import kwant
import numpy as np

from edge_probe_readout import (
    floating_response, predict_finite_probe_readout,
    predict_weak_probe_readout, sample_grid, scattering_checks,
    site_injectivities,
)
from compute_floating_study import solve_floating

S0 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.diag([1., -1.]).astype(complex)
I4 = np.eye(4, dtype=complex)


@dataclass(frozen=True)
class QSHParameters:
    """Orbital tensor spin basis, canonical +sin(ky) QSH model."""
    m: float = -1.
    delta: float = .2
    t: float = 1.
    A: float = 1.
    length: int = 31
    width: int = 12

    def __post_init__(self):
        if self.length < 5 or int(self.length) != self.length or self.length % 2 != 1:
            raise ValueError('Use odd length >=5 for contacts on the longitudinal reflection plane')
        if self.width < 4 or int(self.width) != self.width:
            raise ValueError('Use an integer width >=4')
        if not all(np.isfinite(value) for value in asdict(self).values()):
            raise ValueError('Parameters must be finite')
        if self.t != 1.:
            raise ValueError('Normalize all hoppings and energies to t, then set t=1')
        if self.A <= 0:
            raise ValueError('A must be positive')


class QSHReadoutDevice:
    """Matched QSH source leads with optional identical metallic probes."""
    def __init__(self, parameters, tc=None):
        if tc is not None and (not np.isfinite(tc) or tc <= 0):
            raise ValueError('Probe hopping amplitude must be finite and positive')
        self.parameters = parameters
        self.tc = tc
        self.tc3 = self.tc4 = tc
        self.length, self.width = parameters.length, parameters.width
        self.delta = parameters.delta
        self.geometry = 'bare_x' if tc is None else 'four_x'
        self.weights = np.array([.5, -.5] if tc is None else [.5, -.5, 0., 0.])
        self.project = np.array([1., -1.] if tc is None else [0., 0., 1., -1.])
        lat = kwant.lattice.square(norbs=4, name='qsh_readout_sample')
        metal = kwant.lattice.square(norbs=4, name='qsh_readout_probe')
        h0 = parameters.m*np.kron(SZ, S0)+parameters.delta*np.kron(SX, S0)
        hops = {
            (1, 0): (parameters.t*np.kron(SZ, S0)
                     +1j*parameters.A*np.kron(SX, SZ))/2,
            (0, 1): (parameters.t*np.kron(SZ, S0)
                     +1j*parameters.A*np.kron(SY, S0))/2,
        }

        def sample_onsite(site, potential):
            x, y = site.tag
            return h0+potential[y, x]*I4

        def shifted(index, base):
            def onsite(site, voltages):
                return base+voltages[index]*I4
            return onsite

        builder = kwant.Builder()
        for x in range(self.length):
            for y in range(self.width):
                builder[lat(x, y)] = sample_onsite
        for d, hop in hops.items():
            builder[kwant.builder.HoppingKind(d, lat)] = hop
        for index, direction in enumerate((-1, 1)):
            lead = kwant.Builder(kwant.TranslationalSymmetry((direction, 0)))
            for y in range(self.width):
                lead[lat(0, y)] = shifted(index, h0)
            for d, hop in hops.items():
                lead[kwant.builder.HoppingKind(d, lat)] = hop
            builder.attach_lead(lead)
        xc = self.length//2
        if tc is not None:
            for index, side in ((2, 1), (3, -1)):
                ybuf, edge = (self.width, self.width-1) if side == 1 else (-1, 0)
                builder[metal(xc, ybuf)] = shifted(index, np.zeros((4, 4)))
                builder[metal(xc, ybuf), lat(xc, edge)] = tc*I4
                lead = kwant.Builder(kwant.TranslationalSymmetry((0, side)))
                lead[metal(xc, 0)] = shifted(index, np.zeros((4, 4)))
                lead[kwant.builder.HoppingKind((0, 1), metal)] = -I4
                builder.attach_lead(lead)
        self.system = builder.finalized()
        self.positions = np.array([site.pos for site in self.system.sites], dtype=int)
        self.pos = self.positions
        self.sample = np.array([site.family == lat for site in self.system.sites])
        self.nlead = len(self.system.leads)
        self.zero = np.zeros((self.width, self.length))
        self.vzero = np.zeros(self.nlead)
        self.density = kwant.operator.Density(self.system)
        self.contact_positions = np.array([[xc, self.width-1], [xc, 0]], dtype=int)
        self.contact_indices = np.array([
            next(i for i, pos in enumerate(self.positions)
                 if self.sample[i] and np.array_equal(pos, point))
            for point in self.contact_positions
        ])

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


def build_qsh_device(parameters=QSHParameters(), tc=None):
    return QSHReadoutDevice(parameters, tc=tc)


def bare_edge_quantities(parameters, energy, energy_step=1e-5, voltage_step=1e-5):
    """Extract contact spectra and injection response without voltage probes.

    The first-order local-neutrality U/V=(nuL-nuR)/(2*rho) comes from the
    bare sample.  Physical positive source voltage uses Hamiltonian source
    shifts [-V/2,+V/2] and sample shift -U.  No four-port data enter this
    independently calculated weak-probe predictor.
    """
    if energy_step <= 0 or voltage_step <= 0:
        raise ValueError('Derivative steps must be positive')
    device = build_qsh_device(parameters)
    contacts = device.contact_indices
    nu = site_injectivities(device, energy)
    rho = nu.sum(axis=0)
    if np.any(rho <= 0):
        raise RuntimeError('Local neutrality undefined at zero DOS')
    eta = (nu[0]-nu[1])/rho
    u = eta/2
    u_grid = sample_grid(device, u)
    rho_plus = site_injectivities(device, energy+energy_step).sum(axis=0)[contacts]
    rho_minus = site_injectivities(device, energy-energy_step).sum(axis=0)[contacts]
    dE_ln_rho = (np.log(rho_plus)-np.log(rho_minus))/(2*energy_step)
    eta_values = []
    for voltage in (voltage_step, -voltage_step):
        nu_bias = site_injectivities(
            device, energy, potential=-voltage*u_grid,
            voltages=np.array([-voltage/2, voltage/2]))[:, contacts]
        eta_values.append((nu_bias[0]-nu_bias[1])/nu_bias.sum(axis=0))
    dV_eta = (eta_values[0]-eta_values[1])/(2*voltage_step)
    predicted = predict_weak_probe_readout(dE_ln_rho, dV_eta)
    # kwant.ldos returns orbital-resolved values; density sums all orbitals.
    ldos = kwant.ldos(device.system, energy, params=device.params()).reshape(-1, 4).sum(axis=1)
    checks = scattering_checks(device, energy)
    checks.update(
        ldos_vs_sum_injectivities_error=float(np.max(abs(ldos-rho))),
        equilibrium_eta_contact_abs_max=float(np.max(abs(eta[contacts]))),
        equilibrium_U_over_V_max=float(np.max(abs(u))),
        equilibrium_rho_Mx_error=float(np.max(abs(
            sample_grid(device, rho)-sample_grid(device, rho)[:, ::-1]))),
    )
    return dict(
        energy=float(energy), contact_order=['top', 'bottom'],
        contact_positions=device.contact_positions.copy(),
        nu_L=nu[0, contacts].copy(), nu_R=nu[1, contacts].copy(),
        rho=rho[contacts].copy(), eta=eta[contacts].copy(),
        dE_ln_rho=dE_ln_rho, dV_eta=dV_eta,
        eta_plus=eta_values[0], eta_minus=eta_values[1],
        nu_atomic=sample_grid(device, nu), rho_atomic=sample_grid(device, rho),
        U_over_V_atomic=u_grid, energy_step=float(energy_step),
        voltage_step=float(voltage_step), checks=checks, **predicted,
    )


def four_terminal_readout(parameters, energy, tc, step=1e-5):
    """Identical finite-coupling probes, with the full floating readout."""
    if abs(energy) >= 2:
        raise ValueError('The hopping1 metallic probes require |E|<2')
    device = build_qsh_device(parameters, tc=tc)
    raw, _, field, _ = floating_response(device, energy, step=step)
    checks = scattering_checks(device, energy)
    trans = checks['transmission']
    contacts = device.contact_indices
    source_nu = site_injectivities(device, energy)[:2, contacts]
    source_rho = source_nu.sum(axis=0)
    source_eta = (source_nu[0]-source_nu[1])/source_rho
    alpha_raw = np.asarray(raw['alpha'])
    eta_bias = []
    for voltage in (step, -step):
        nu_bias = site_injectivities(
            device, energy, potential=-voltage*field/2,
            voltages=-voltage*alpha_raw/2)[:2, contacts]
        eta_bias.append((nu_bias[0]-nu_bias[1])/nu_bias.sum(axis=0))
    dV_eta = (eta_bias[0]-eta_bias[1])/(2*step)
    rho_energy = [
        site_injectivities(device, probe_energy)[:2, contacts].sum(axis=0)
        for probe_energy in (energy+step, energy-step)
    ]
    contact_dE_ln_rho = (np.log(rho_energy[0])-np.log(rho_energy[1]))/(2*step)
    probe_gamma = np.repeat(tc**2*np.sqrt(4-energy**2), 2)
    F = trans[2:, :2].sum(axis=1)
    rho_from_F = F/(2*np.pi*probe_gamma)
    checks.update(
        source_contact_rho_vs_F_over_gamma_error=float(np.max(abs(source_rho-rho_from_F))),
        equilibrium_source_eta_contact_abs_max=float(np.max(abs(source_eta))),
    )
    contact_prediction = predict_finite_probe_readout(
        contact_dE_ln_rho, dV_eta, F, float(trans[2, 3]))
    contact_prediction['contactformula_error'] = abs(
        contact_prediction['predicted_kappa_full_from_contacts']+raw['kH']/4)
    return dict(
        energy=float(energy), tc=float(tc), tc3=float(tc), tc4=float(tc),
        kappa_full=-raw['kH']/4, kappa_frozen=-raw['kH_frozen']/4,
        kappa_internal=-raw['kH_internal']/4,
        kappa_electrodes=-raw['kH_electrodes']/4,
        F=F, interprobe_T=float(trans[2, 3]),
        source_to_probe_transmission=trans[2:, :2].copy(),
        source_contact_nu_L=source_nu[0], source_contact_nu_R=source_nu[1],
        source_contact_rho=source_rho, source_contact_eta=source_eta,
        contact_dV_eta=dV_eta, contact_dE_ln_rho=contact_dE_ln_rho,
        probe_gamma=probe_gamma, **contact_prediction,
        paper_alpha_contacts=alpha_raw[2:]/2,
        paper_beta_contacts=-np.asarray(raw['beta'])[2:]/4,
        field_code=field, raw=raw, checks=checks,
        convention='kappa=t*kappa_physical/e=-kH_code/4; V_perp=V3-V4, full source-drain voltage',
    )


def finite_bias_check(parameters, energy, tc, voltage=1e-4, step=1e-5):
    """Check the quadratic derivative with directly solved floating currents.

    ``voltage`` is the full, positive physical source-drain difference.
    The internal potential uses the equilibrium first-order neutrality
    closure, exactly as in the derivative calculation.  Both signs of the
    voltage are solved independently before taking the even Hall signal.
    """
    if not np.isfinite(voltage) or voltage <= 0:
        raise ValueError('Use a finite positive full source-drain voltage')
    device = build_qsh_device(parameters, tc=tc)
    raw, u_grid, _, _ = floating_response(device, energy, step=step)
    records = [
        solve_floating(device, energy, -physical_voltage/2, u_grid, raw)
        for physical_voltage in (voltage, -voltage)
    ]
    # Physical probe voltage is minus the onsite-positive solver voltage.
    hall_voltages = np.array([-(v[2]-v[3]) for v, _ in records])
    currents = np.array([current for _, current in records])
    direct = float(hall_voltages.sum()/(2*voltage**2))
    derivative = -float(raw['kH'])/4
    error = abs(direct-derivative)
    return dict(
        kappa_direct=direct, kappa_derivative=derivative,
        absolute_error=error, relative_error=error/max(abs(derivative), 1e-12),
        probe_current_error=float(np.max(abs(currents[:, 2:]))),
        current_sum_error=float(np.max(abs(currents.sum(axis=1)))),
        voltage=float(voltage), physical_probe_voltages=np.array([-v[2:] for v, _ in records]),
        physical_hall_voltages=hall_voltages,
        convention='Full source-drain voltage; physical voltages=-onsite-positive solver shifts',
    )
