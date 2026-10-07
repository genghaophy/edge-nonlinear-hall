"""Reusable bare-edge to floating-probe readout interfaces.

The numerical Hamiltonian uses positive onsite energy and source energy
shifts. The paper convention is mu=EF-eV, H=H0-eU, source voltages
+/-V_SD/2, and V_perp=V_top-V_bottom. Energy is measured in t2 and the
returned dimensionless kappa2 is t2*kappa2_physical/e. Conversion of the
historical floating solver is exactly kappa2=-kH_code/4.

The pure ``predict_weak_probe_readout`` interface can also consume edge
quantities from other two-terminal models. Its contact ordering is always
[top, bottom]. The RM construction and extraction functions here provide
one concrete implementation. Weak-probe predictions and finite-coupling
four-terminal responses are separate outputs, never fitted to one another.
"""
from dataclasses import asdict, dataclass
from pathlib import Path
import sys

import kwant
import numpy as np

RM_CODE=Path(__file__).resolve().parent
if str(RM_CODE) not in sys.path:
    sys.path.insert(0,str(RM_CODE))
from compute_floating_study import floating_response


@dataclass(frozen=True)
class RMParameters:
    """One spinless orbital per atomic row; two rows make a y cell."""
    tx: float=.25
    t1: float=.3
    t2: float=1.
    delta: float=.2
    cells: int=8
    length: int=31

    def __post_init__(self):
        if self.cells < 1 or int(self.cells) != self.cells:
            raise ValueError('cells must be a positive integer')
        if self.length < 3 or int(self.length) != self.length:
            raise ValueError('length must be an integer >=3')
        if self.length % 2 == 0:
            raise ValueError('Use odd length so the contact lies on the longitudinal reflection plane')
        if self.tx <= 0 or self.t1 <= 0 or self.t2 <= 0:
            raise ValueError('This RM implementation requires positive hoppings')
        if not all(np.isfinite(value) for value in asdict(self).values()):
            raise ValueError('Parameters must be finite')
        if self.t2 != 1.:
            raise ValueError('Normalize all hoppings and energies to t2 first, then set t2=1')


class RMReadoutDevice:
    """Device-compatible RM strip with optional independently coupled probes.

    tc3 and tc4 are hopping amplitudes, not linewidths. Both omitted gives
    the completely bare two-terminal sample; both positive gives lead
    ordering L,R,top,bottom. Zero coupling is excluded because floating
    voltages would be singular. Fixed unequal tc3/tc4 is allowed.
    """
    def __init__(self,parameters,tc3=None,tc4=None):
        if (tc3 is None) != (tc4 is None):
            raise ValueError('Supply both probe couplings, or omit both')
        if tc3 is not None and (tc3 <= 0 or tc4 <= 0
                               or not np.isfinite(tc3+tc4)):
            raise ValueError('Floating probe hopping amplitudes must be finite and positive')
        self.parameters=parameters
        self.tc3,self.tc4=tc3,tc4
        self.width,self.length=2*parameters.cells,parameters.length
        self.delta=parameters.delta
        self.geometry='bare_x' if tc3 is None else 'four_x'
        self.weights=np.array([.5,-.5] if tc3 is None else [.5,-.5,0.,0.])
        self.project=np.array([1.,-1.] if tc3 is None else [0.,0.,1.,-1.])
        lat=kwant.lattice.square(norbs=1,name='rm_readout_sample')
        metal=kwant.lattice.square(norbs=1,name='rm_readout_probe')
        def onsite(site,potential):
            x,y=site.tag
            return parameters.delta*(-1)**y+potential[y,x]
        def vertical(site_to,site_from):
            y=min(site_to.tag[1],site_from.tag[1])
            return parameters.t1 if y % 2 == 0 else parameters.t2
        def shifted(index,sample):
            def lead_onsite(site,voltages):
                base=parameters.delta*(-1)**site.tag[1] if sample else 0.
                return base+voltages[index]
            return lead_onsite
        builder=kwant.Builder()
        for x in range(self.length):
            for y in range(self.width):
                builder[lat(x,y)]=onsite
        builder[kwant.builder.HoppingKind((1,0),lat)]=-parameters.tx
        builder[kwant.builder.HoppingKind((0,1),lat)]=vertical
        for index,direction in enumerate((-1,1)):
            lead=kwant.Builder(kwant.TranslationalSymmetry((direction,0)))
            for y in range(self.width):
                lead[lat(0,y)]=shifted(index,True)
            lead[kwant.builder.HoppingKind((1,0),lat)]=-parameters.tx
            lead[kwant.builder.HoppingKind((0,1),lat)]=vertical
            builder.attach_lead(lead)
        xc=self.length//2
        if tc3 is not None:
            for index,side,coupling in ((2,1,tc3),(3,-1,tc4)):
                ybuf,edge=(self.width,self.width-1) if side == 1 else (-1,0)
                builder[metal(xc,ybuf)]=shifted(index,False)
                builder[metal(xc,ybuf),lat(xc,edge)]=coupling
                lead=kwant.Builder(kwant.TranslationalSymmetry((0,side)))
                lead[metal(xc,0)]=shifted(index,False)
                lead[kwant.builder.HoppingKind((0,1),metal)]=-1.
                builder.attach_lead(lead)
        self.system=builder.finalized()
        self.positions=np.asarray([s.pos for s in self.system.sites],int)
        self.sample=np.array([s.family == lat for s in self.system.sites])
        self.nlead=len(self.system.leads)
        self.zero=np.zeros((self.width,self.length))
        self.vzero=np.zeros(self.nlead)
        self.density=kwant.operator.Density(self.system)
        self.contact_positions=np.array([[xc,self.width-1],[xc,0]],int)
        self.contact_indices=np.array([next(i for i,pos in enumerate(self.positions)
            if self.sample[i] and np.array_equal(pos,point))
            for point in self.contact_positions])

    def params(self,potential=None,voltages=None):
        return dict(potential=self.zero if potential is None else potential,
                    voltages=self.vzero if voltages is None else voltages)

    def matrix(self,energy,potential=None,voltages=None):
        sm=kwant.smatrix(self.system,energy,params=self.params(potential,voltages))
        trans=np.zeros((self.nlead,self.nlead))
        for a in range(self.nlead):
            for b in range(self.nlead):
                if a != b:
                    trans[a,b]=sm.transmission(a,b)
        return np.diag(trans.sum(axis=0))-trans


def build_rm_device(parameters=RMParameters(),tc3=None,tc4=None):
    return RMReadoutDevice(parameters,tc3=tc3,tc4=tc4)


def site_injectivities(device,energy,potential=None,voltages=None):
    """Sum every unit-flux incoming state from each lead; nu=|psi|²/(2pi)."""
    params=device.params(potential,voltages)
    wave=kwant.wave_function(device.system,energy,params=params)
    return np.array([sum((device.density(state,params=params) for state in wave(lead)),
                        np.zeros(len(device.positions)))/(2*np.pi)
                     for lead in range(device.nlead)])


def sample_grid(device,site_array):
    values=np.asarray(site_array)
    grid=np.zeros((*values.shape[:-1],device.width,device.length))
    for index in np.flatnonzero(device.sample):
        x,y=device.positions[index]
        grid[...,y,x]=values[...,index]
    return grid


def predict_weak_probe_readout(dE_ln_rho,dV_eta):
    """Map bare contact data [top,bottom] to a weak floating Hall voltage.

    Inputs use E/t2 and e*V_SD/t2 derivatives. Assume longitudinal
    reflection, equal equilibrium L/R injectivity at each contact, and
    negligible interprobe transmission. Identical probe reservoir spectral
    functions suffice; constant unequal probe hoppings cancel from the
    energy-logarithmic ratio as both hoppings approach zero.
    """
    slope=np.asarray(dE_ln_rho,float)
    eta_prime=np.asarray(dV_eta,float)
    if slope.shape != (2,) or eta_prime.shape != (2,):
        raise ValueError('Both inputs must have [top,bottom] shape (2,)')
    if not np.all(np.isfinite(np.r_[slope,eta_prime])):
        raise ValueError('Contact quantities must be finite')
    dos=-float(slope[0]-slope[1])/8
    injection=float(eta_prime[0]-eta_prime[1])/2
    return dict(kappa_DOS=dos,kappa_injectivity=injection,
                kappa_predicted=dos+injection)


def scattering_checks(device,energy):
    sm=kwant.smatrix(device.system,energy,params=device.params())
    trans=np.array([[sm.transmission(a,b) for b in range(device.nlead)]
                    for a in range(device.nlead)])
    incoming=np.array([sm.num_propagating(b) for b in range(device.nlead)])
    return dict(transmission=trans,incoming_modes=incoming,
        unitarity_error=float(np.max(abs(trans.sum(axis=0)-incoming))),
        reciprocity_error=float(np.max(abs(trans-trans.T))))


def predict_finite_probe_readout(dE_ln_rho,dV_eta,F,interprobe_T):
    """Map probe-dressed contact quantities to the finite-probe readout.

    This is distinct from extracting a bare two-port predictor. It assumes
    Mx/TRS, zero linear floating probe voltages, common probe spectral shape,
    and scalar point contacts. F=[Ftop,Fbottom] sums source-to-probe
    transmissions, while interprobe_T is their reciprocal mixing g.
    """
    F=np.asarray(F,float)
    if F.shape != (2,) or not np.all(np.isfinite(F)) or np.any(F <= 0):
        raise ValueError('F must contain two finite positive source absorption sums')
    if not np.isfinite(interprobe_T) or interprobe_T < 0:
        raise ValueError('Interprobe transmission must be finite and nonnegative')
    components=predict_weak_probe_readout(dE_ln_rho,dV_eta)
    mixing=1+float(interprobe_T)*float(np.sum(1/F))
    return dict(contact_kappa_spectral=components['kappa_DOS'],
        contact_kappa_injection=components['kappa_injectivity'],
        interprobe_mixing=mixing,
        predicted_kappa_full_from_contacts=components['kappa_predicted']/mixing)


def bare_edge_quantities(parameters,energy,energy_step=1e-5,voltage_step=1e-5):
    """Extract bare two-terminal contact DOS slopes and driven injectivities.

    The first-order local-neutrality U/V_SD=(nu_L-nu_R)/(2*rho) is computed
    from the bare equilibrium sample. For physical positive V_SD, source
    Hamiltonian shifts are [-V_SD/2,+V_SD/2] and sample shift is -U.
    This function does not add a voltage probe or fit a four-terminal result.
    """
    if energy_step <= 0 or voltage_step <= 0:
        raise ValueError('Derivative steps must be positive')
    device=build_rm_device(parameters)
    contacts=device.contact_indices
    nu=site_injectivities(device,energy)
    rho=nu.sum(axis=0)
    if np.any(rho <= 0):
        raise RuntimeError('Local neutrality is undefined at zero DOS')
    eta=(nu[0]-nu[1])/rho
    u=eta/2
    u_grid=sample_grid(device,u)
    rho_plus=site_injectivities(device,energy+energy_step).sum(axis=0)[contacts]
    rho_minus=site_injectivities(device,energy-energy_step).sum(axis=0)[contacts]
    dE_ln_rho=(np.log(rho_plus)-np.log(rho_minus))/(2*energy_step)
    eta_values=[]
    for voltage in (voltage_step,-voltage_step):
        nu_bias=site_injectivities(device,energy,potential=-voltage*u_grid,
                                  voltages=np.array([-voltage/2,voltage/2]))[:,contacts]
        eta_values.append((nu_bias[0]-nu_bias[1])/nu_bias.sum(axis=0))
    dV_eta=(eta_values[0]-eta_values[1])/(2*voltage_step)
    predicted=predict_weak_probe_readout(dE_ln_rho,dV_eta)
    ldos=kwant.ldos(device.system,energy,params=device.params())
    checks=scattering_checks(device,energy)
    checks['ldos_vs_sum_injectivities_error']=float(np.max(abs(ldos-rho)))
    checks['equilibrium_eta_contact_abs_max']=float(np.max(abs(eta[contacts])))
    checks['equilibrium_U_over_V_max']=float(np.max(abs(u)))
    checks['equilibrium_rho_Mx_error']=float(np.max(abs(
        sample_grid(device,rho)-sample_grid(device,rho)[:,::-1])))
    return dict(energy=float(energy),contact_order=['top','bottom'],
        contact_positions=device.contact_positions.copy(),
        nu_L=nu[0,contacts].copy(),nu_R=nu[1,contacts].copy(),
        rho=rho[contacts].copy(),eta=eta[contacts].copy(),
        dE_ln_rho=dE_ln_rho,dV_eta=dV_eta,
        eta_plus=eta_values[0],eta_minus=eta_values[1],
        nu_atomic=sample_grid(device,nu),rho_atomic=sample_grid(device,rho),
        U_over_V_atomic=u_grid,energy_step=float(energy_step),
        voltage_step=float(voltage_step),checks=checks,**predicted)


def four_terminal_readout(parameters,energy,tc3,tc4,step=1e-5):
    """Actual finite-coupling four-terminal result, converted exactly once."""
    if abs(energy) >= 2:
        raise ValueError('The hopping1 probe reservoirs require |E|<2')
    device=build_rm_device(parameters,tc3=tc3,tc4=tc4)
    raw,_,field,_=floating_response(device,energy,step=step)
    checks=scattering_checks(device,energy)
    trans=checks['transmission']
    contacts=device.contact_indices
    source_nu=site_injectivities(device,energy)[:2,contacts]
    source_rho=source_nu.sum(axis=0)
    source_eta=(source_nu[0]-source_nu[1])/source_rho
    alpha_raw=np.asarray(raw['alpha'])
    eta_bias=[]
    for voltage in (step,-step):
        # Physical source voltage +/-V/2 corresponds to code amplitude -V/2.
        # Keep the solved linear floating probe voltages (zero by Mx to
        # numerical precision), rather than silently grounding the probes.
        nu_bias=site_injectivities(device,energy,potential=-voltage*field/2,
                                  voltages=-voltage*alpha_raw/2)[:2,contacts]
        eta_bias.append((nu_bias[0]-nu_bias[1])/nu_bias.sum(axis=0))
    dV_eta=(eta_bias[0]-eta_bias[1])/(2*step)
    rho_energy=[]
    for probe_energy in (energy+step,energy-step):
        rho_energy.append(site_injectivities(device,probe_energy)[:2,contacts].sum(axis=0))
    contact_dE_ln_rho=(np.log(rho_energy[0])-np.log(rho_energy[1]))/(2*step)
    probe_gamma=np.array([tc3**2,tc4**2])*np.sqrt(4-energy**2)
    rho_from_F=trans[2:,:2].sum(axis=1)/(2*np.pi*probe_gamma)
    checks['source_contact_rho_vs_F_over_gamma_error']=float(np.max(abs(source_rho-rho_from_F)))
    checks['equilibrium_source_eta_contact_abs_max']=float(np.max(abs(source_eta)))
    contact_prediction=predict_finite_probe_readout(contact_dE_ln_rho,dV_eta,
        trans[2:,:2].sum(axis=1),float(trans[2,3]))
    contact_prediction['contactformula_error']=abs(
        contact_prediction['predicted_kappa_full_from_contacts']+raw['kH']/4)
    return dict(energy=float(energy),tc3=float(tc3),tc4=float(tc4),
        kappa_full=-raw['kH']/4,kappa_frozen=-raw['kH_frozen']/4,
        kappa_internal=-raw['kH_internal']/4,
        kappa_electrodes=-raw['kH_electrodes']/4,
        F=trans[2:,:2].sum(axis=1),interprobe_T=float(trans[2,3]),
        source_to_probe_transmission=trans[2:,:2].copy(),
        source_contact_nu_L=source_nu[0],source_contact_nu_R=source_nu[1],
        source_contact_rho=source_rho,source_contact_eta=source_eta,
        contact_dV_eta=dV_eta,contact_dE_ln_rho=contact_dE_ln_rho,
        probe_gamma=probe_gamma,**contact_prediction,
        paper_alpha_contacts=alpha_raw[2:]/2,
        paper_beta_contacts=-np.asarray(raw['beta'])[2:]/4,
        field_code=field,raw=raw,checks=checks,
        convention='kappa=t2*kappa_physical/e=-kH_code/4; V_perp=V3-V4, full source-drain voltage')


def analytic_rm_contact_quantities(parameters,energy):
    """Independent exact matched-strip DOS check, not the production prediction."""
    p=parameters
    rows=2*p.cells
    hy=np.diag(p.delta*(-1.)**np.arange(rows))
    for y in range(rows-1):
        hy[y,y+1]=hy[y+1,y]=p.t1 if y % 2 == 0 else p.t2
    eig,psi=np.linalg.eigh(hy)
    shifted=energy-eig
    active=abs(shifted)<2*p.tx
    denom=(2*p.tx)**2-shifted[active]**2
    dos=np.zeros(rows)
    derivative=np.zeros(rows)
    dos[active]=1/(np.pi*np.sqrt(denom))
    derivative[active]=shifted[active]/(np.pi*denom**1.5)
    weight=abs(psi[[-1,0]])**2
    rho=weight@dos
    return dict(rho=rho,dE_ln_rho=(weight@derivative)/rho,
                propagating_subbands=int(active.sum()))
