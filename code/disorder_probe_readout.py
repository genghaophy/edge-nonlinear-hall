"""Scalar-disorder audit of independent two-port and floating four-port readout.

Disorder is present only in the sample. Source leads retain exactly the
clean sample Hamiltonian and hoppings. Energies and all returned quadratic
coefficients use t=1 (QSH) or t2=1 (RM), physical mu=EF-eV and the full
source-drain voltage. The historical four-port coefficient is converted
once, kappa2=-kH_code/4. Identical probes always satisfy t3=t4=tau.
"""
from dataclasses import asdict
from pathlib import Path
import hashlib
import json

import kwant
import numpy as np

from edge_probe_readout import (
    RMParameters, RMReadoutDevice, floating_response,
    sample_grid, scattering_checks, site_injectivities,
)
from qsh_edge_probe_readout import QSHParameters, QSHReadoutDevice
from compute_floating_study import solve_floating


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    return value


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(jsonable(value), indent=2, allow_nan=False) + '\n',
                         encoding='utf-8')
    temporary.replace(path)


class _SampleDisorder:
    def set_disorder(self, disorder):
        self.disorder = np.asarray(disorder, float).copy()
        if self.disorder.shape != self.zero.shape or not np.all(np.isfinite(self.disorder)):
            raise ValueError('Use a finite real sample disorder array of shape (width,length).')

    def params(self, potential=None, voltages=None):
        # The clean lead onsite functions receive only voltages, never this
        # array. The potential argument represents a bias-induced shift.
        shift = self.zero if potential is None else np.asarray(potential)
        return super().params(potential=self.disorder + shift, voltages=voltages)


class DisorderedQSH(_SampleDisorder, QSHReadoutDevice):
    def __init__(self, parameters, disorder, tau=None):
        QSHReadoutDevice.__init__(self, parameters, tc=tau)
        self.set_disorder(disorder)


class DisorderedRM(_SampleDisorder, RMReadoutDevice):
    def __init__(self, parameters, disorder, tau=None):
        RMReadoutDevice.__init__(self, parameters, tc3=tau, tc4=tau)
        self.set_disorder(disorder)


def make_disorder(model, parameters, strength, seed, mirror=False):
    width = parameters.width if model == 'qsh' else 2 * parameters.cells
    values = np.random.default_rng(int(seed)).uniform(-.5, .5,
                                                     (width, parameters.length))
    if mirror:
        # Copy a complete half-realization to keep the uniform single-site
        # distribution; do not average two draws and narrow the distribution.
        for x in range(parameters.length // 2):
            values[:, -1-x] = values[:, x]
    return float(strength) * values


def build(model, parameters, disorder, tau=None):
    return (DisorderedQSH if model == 'qsh' else DisorderedRM)(parameters, disorder, tau)


def source_contract():
    project = Path(__file__).resolve().parents[1]
    files = [Path(__file__), project/'code/edge_probe_readout.py',
             project/'code/qsh_edge_probe_readout.py',
             project/'code/compute_floating_study.py']
    return {str(p): digest(p) for p in files}


def bare_quantities(model, parameters, energy, disorder, step=1e-5):
    """General weak-tunneling predictor, allowing nonzero linear probe voltage.

    beta_p=-(1-eta_p^2)*dE ln(Gamma_p*nu_p)/8+dV eta_p/2.
    The identical chain probes have Gamma=tau^2*sqrt(4-E^2).
    In generic disorder eta need not vanish and their common linewidth
    energy derivative need not cancel between the two contacts.
    """
    dev = build(model, parameters, disorder)
    contacts = dev.contact_indices
    nu = site_injectivities(dev, energy)
    rho = nu.sum(axis=0)
    if np.any(rho <= 0):
        raise RuntimeError('Local neutrality is undefined at zero spectral density.')
    eta_all = (nu[0]-nu[1]) / rho
    eta = eta_all[contacts]
    u = sample_grid(dev, eta_all/2)
    rho_energy = [site_injectivities(dev, energy+s*step).sum(axis=0)[contacts]
                  for s in (1, -1)]
    slopes = (np.log(rho_energy[0])-np.log(rho_energy[1]))/(2*step)
    eta_driven = []
    for voltage in (step, -step):
        vnu = site_injectivities(dev, energy, potential=-voltage*u,
                                voltages=np.array([-voltage/2, voltage/2]))[:, contacts]
        eta_driven.append((vnu[0]-vnu[1])/vnu.sum(axis=0))
    eta_dot = (eta_driven[0]-eta_driven[1])/(2*step)
    gamma_slope = -float(energy)/(4-float(energy)**2)
    beta_spectral = -(1-eta**2)*slopes/8
    beta_gamma = -(1-eta**2)*gamma_slope/8
    beta_injection = eta_dot/2
    difference = lambda a: float(a[0]-a[1])
    beta = beta_spectral + beta_gamma + beta_injection
    checks = scattering_checks(dev, energy)
    norbs = 4 if model == 'qsh' else 1
    ldos = kwant.ldos(dev.system, energy, params=dev.params()).reshape(-1, norbs).sum(axis=1)
    checks['ldos_vs_injectivity_error'] = float(np.max(abs(ldos-rho)))
    sm = kwant.smatrix(dev.system, energy, params=dev.params())
    return dict(
        energy=float(energy), contact_positions=dev.contact_positions,
        rho=rho[contacts], nu_L=nu[0, contacts], nu_R=nu[1, contacts],
        eta=eta, dE_ln_rho=slopes, dV_eta=eta_dot,
        dE_ln_gamma=gamma_slope, alpha=eta/2, beta=beta,
        kappa1=difference(eta)/2,
        kappa_general=difference(beta),
        kappa_spectral=-difference(slopes)/8,
        kappa_mirror=-difference(slopes)/8 + difference(eta_dot)/2,
        kappa_spectral_general=difference(beta_spectral),
        kappa_gamma=difference(beta_gamma),
        kappa_injection=difference(beta_injection),
        source_reflection=float(sm.transmission(0, 0)),
        source_transmission=float(sm.transmission(1, 0)),
        source_modes=int(sm.num_propagating(0)),
        U_over_V_grid=u, step=float(step), checks=checks,
    )


def four_quantities(model, parameters, energy, disorder, tau, step=1e-5):
    dev = build(model, parameters, disorder, tau)
    raw, _, _, _ = floating_response(dev, energy, step=step)
    checks = scattering_checks(dev, energy)
    return dict(energy=float(energy), tau=float(tau),
                kappa_full=-float(raw['kH'])/4,
                kappa_frozen=-float(raw['kH_frozen'])/4,
                kappa_internal=-float(raw['kH_internal'])/4,
                kappa_electrodes=-float(raw['kH_electrodes'])/4,
                kappa1=float(raw['alpha_H'])/2,
                alpha=np.asarray(raw['alpha'])[2:]/2,
                beta=-np.asarray(raw['beta'])[2:]/4,
                raw=raw, checks=checks, step=float(step))


def direct_check(model, parameters, energy, disorder, tau,
                 voltage=1e-4, step=1e-5):
    dev = build(model, parameters, disorder, tau)
    raw, u_grid, _, _ = floating_response(dev, energy, step=step)
    records = [solve_floating(dev, energy, -v/2, u_grid, raw)
               for v in (voltage, -voltage)]
    hall = np.array([-(va[2]-va[3]) for va, _ in records])
    currents = np.array([ip for _, ip in records])
    direct = float(hall.sum()/(2*voltage**2))
    return dict(voltage=float(voltage), kappa_direct=direct,
                kappa_derivative=-float(raw['kH'])/4,
                absolute_error=abs(direct+float(raw['kH'])/4),
                linear_hall=float((hall[0]-hall[1])/(2*voltage)),
                derivative_linear_hall=float(raw['alpha_H'])/2,
                probe_current_error=float(np.max(abs(currents[:, 2:]))),
                current_sum_error=float(np.max(abs(currents.sum(axis=1)))))
