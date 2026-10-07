"""Independent bare-edge scattering-state check of weak-probe readout.

Compute the voltage derivative of eta=(nu_L-nu_R)/(nu_L+nu_R) directly
from scattering states of a completely bare two-terminal RM sample.
Physical source voltages are +/-V_SD/2; the numerical Hamiltonian onsite
shifts are -V_SD/2,+V_SD/2, and the sample shift is -U. Units e=t2=1.
The sample's first-order U follows the equilibrium bare local-neutrality
closure, retained explicitly even though it vanishes for the matched strip.
Only a new audit JSON is saved; production data are not modified.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

import kwant
import numpy as np
from threadpoolctl import threadpool_limits

PROJECT=Path(__file__).resolve().parents[1]
OUT=PROJECT/'outputs/checks/rm_bare_injectivity'
SOURCE=PROJECT/'data/raw/rm_figure2/figure2_data.npz'
COUPLING_REPORT=PROJECT/'data/raw/supplementary/weak_probe_audit/rm_tc_scan.json'
RM_ROOT=PROJECT
sys.path.insert(0,str(PROJECT/'code'))
from compute_rm_weak_probe_audit import analytic_contact_dos
from checks._common import write_report, relative


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def make_bare(p):
    lat=kwant.lattice.square(norbs=1,name='bare_biased_rm')
    b=kwant.Builder()
    def onsite(site,potential):
        x,y=site.tag
        return p['delta']*(-1)**y+potential[y,x]
    def hy(s1,s2):
        y=min(s1.tag[1],s2.tag[1])
        return p['t1'] if y % 2 == 0 else p['t2']
    for x in range(p['length']):
        for y in range(2*p['cells']):
            b[lat(x,y)]=onsite
    b[kwant.builder.HoppingKind((1,0),lat)]=-p['tx']
    b[kwant.builder.HoppingKind((0,1),lat)]=hy
    def source_onsite(index):
        def lead_onsite(site,voltages):
            return p['delta']*(-1)**site.tag[1]+voltages[index]
        return lead_onsite
    for index,direction in enumerate((-1,1)):
        lead=kwant.Builder(kwant.TranslationalSymmetry((direction,0)))
        for y in range(2*p['cells']):
            lead[lat(0,y)]=source_onsite(index)
        lead[kwant.builder.HoppingKind((1,0),lat)]=-p['tx']
        lead[kwant.builder.HoppingKind((0,1),lat)]=hy
        b.attach_lead(lead)
    return b.finalized()


def site_nu(system,energy,params):
    wave=kwant.wave_function(system,energy,params=params)
    density=kwant.operator.Density(system)
    return np.array([sum((density(state,params=params) for state in wave(lead)),
                        np.zeros(len(system.sites)))/(2*np.pi) for lead in (0,1)])


def main():
    start=time.perf_counter()
    source=SOURCE
    tracked=[source,COUPLING_REPORT]
    before={relative(path):digest(path) for path in tracked}
    tc_metadata=json.loads(COUPLING_REPORT.read_text(encoding='utf-8'))
    p=tc_metadata['parameters']
    zero=np.zeros((2*p['cells'],p['length']))
    with threadpool_limits(limits=1):
        system=make_bare(p)
        positions=np.array([s.tag for s in system.sites],int)
        contact_ids=[next(i for i,pos in enumerate(positions) if tuple(pos)==(p['length']//2,y))
                     for y in (2*p['cells']-1,0)]
        eqparams=dict(potential=zero,voltages=np.zeros(2))
        nu0=site_nu(system,p['EF'],eqparams)
        rho0=nu0.sum(axis=0)
        # U/V_SD=(nu_L-nu_R)/(2 rho), from V_L=+V_SD/2,V_R=-V_SD/2.
        u_sites=(nu0[0]-nu0[1])/(2*rho0)
        u_grid=zero.copy()
        for i,(x,y) in enumerate(positions):
            u_grid[y,x]=u_sites[i]
        eta0=(nu0[0]-nu0[1])/rho0
        results=[]
        for step in (1e-5,5e-6):
            nus=[]
            for sign in (1,-1):
                V=sign*step
                nus.append(site_nu(system,p['EF'],dict(potential=-V*u_grid,
                    voltages=np.array([-V/2,V/2]))))
            eta=[(nu[0,contact_ids]-nu[1,contact_ids])/nu[:,contact_ids].sum(axis=0)
                 for nu in nus]
            eta_prime=(eta[0]-eta[1])/(2*step)
            _,slope,_=analytic_contact_dos(p,p['EF'])
            kappa_dos=-slope/8
            correction=(eta_prime[0]-eta_prime[1])/2
            results.append(dict(physical_voltage_step=step,
                eta_plus=eta[0].tolist(),eta_minus=eta[1].tolist(),
                eta_voltage_derivative_top_bottom=eta_prime.tolist(),
                kappa_DOS=kappa_dos, kappa_injectivity_feedback=float(correction),
                kappa_predicted_full=float(kappa_dos+correction)))
        ldos_check=float(np.max(abs(kwant.ldos(system,p['EF'],params=eqparams)-rho0)))
    refined=results[-1]
    refinement=abs(results[0]['kappa_predicted_full']-refined['kappa_predicted_full'])
    comparisons=[dict(tc=rec['tc'],full=rec['kappa2_full'],
                     predicted_full=refined['kappa_predicted_full'],
                     difference_full_minus_prediction=rec['kappa2_full']-refined['kappa_predicted_full'],
                     difference_over_tc_squared=(rec['kappa2_full']-refined['kappa_predicted_full'])/rec['tc']**2)
                 for rec in tc_metadata['records']]
    after={relative(path):digest(path) for path in tracked}
    assert before==after
    if ldos_check>1e-10 or refinement>1e-7:
        raise RuntimeError('Independent scattering-state audit validation failed.')
    report=dict(schema='bare_injectivity_weak_probe_audit_v1',parameters=p,
        geometry='Only matched L/R sample leads; no side probes or probe buffer sites',
        physical_convention='mu=EF-eV,H=H0-eU; V_L=+V_SD/2,V_R=-V_SD/2,V_perp=V3-V4; e=t2=1',
        formula='kappa2=-d_E ln(rho_top/rho_bottom)/8 +(partial_V eta_top-partial_V eta_bottom)/2',
        eta_definition='(nu_L-nu_R)/(nu_L+nu_R), bare sample scattering-state injectivities',
        contact_site_indices=contact_ids,contact_positions=positions[contact_ids].tolist(),
        equilibrium_contact_injectivities_LR_by_top_bottom=nu0[:,contact_ids].tolist(),
        equilibrium_eta_top_bottom=eta0[contact_ids].tolist(),
        equilibrium_U_over_V_max=float(np.max(abs(u_sites))),
        equilibrium_U_over_V_at_contact=u_sites[contact_ids].tolist(),
        builtin_ldos_vs_injectivity_sum_max_error=ldos_check,
        steps=results,derivative_refinement_error=float(refinement),
        comparison_to_finite_tc_full=comparisons,
        production_unchanged=True,protected_sha256=after,
        tc_audit_sha256=digest(COUPLING_REPORT),source_sha256=digest(Path(__file__)),
        kwant_version=kwant.__version__,numpy_version=np.__version__,
        runtime_seconds=time.perf_counter()-start)
    path=OUT/'bare_injectivity_response.json'
    write_report(path, report)
    print(json.dumps(dict(output=str(path),equilibrium_U_over_V_max=report['equilibrium_U_over_V_max'],
        refined=refined,refinement=refinement,ldos_error=ldos_check,
        tc_comparisons=comparisons[-3:]),indent=2),flush=True)


# Use code/checks/rm_bare_injectivity.py for configurable audit inputs.
