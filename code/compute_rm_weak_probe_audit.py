"""Independent weak-probe DOS-readout audit; leaves production data intact.

All published dimensionless voltage coefficients use V_SD=V1-V2 and
V_perp=V3-V4, with mu=EF-e*V and H=H0-e*U. Thus kappa2=-kH_code/4
relative to the historical positive-energy solver with source energy
shifts +v,-v and VH=V3-V4. The scalar DOS prediction in these units is
-1/8*d_E log(rho_top/rho_bottom), with E in t2 units. Full and frozen
responses are compared separately; weak coupling is not assumed to remove
electrode or electrostatic feedback.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

import kwant
import numpy as np
from threadpoolctl import threadpool_limits

PROJECT = Path(__file__).resolve().parents[1]
RM_ROOT = PROJECT
sys.path.insert(0, str(RM_ROOT/'code'))
from compute_initial_study import Device
from compute_floating_study import floating_response
from checks._common import write_report, relative

OUT = PROJECT / 'outputs/checks/rm_weak_probe'
SOURCE = PROJECT / 'data/raw/rm_figure2/figure2_data.npz'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bare_system(p):
    """Matched L/R strip with no metallic side probes at all."""
    lat = kwant.lattice.square(norbs=1, name='bare_rm')
    b = kwant.Builder()
    for x in range(p['length']):
        for y in range(2*p['cells']):
            b[lat(x,y)] = p['delta'] * (-1)**y
    def hy(s1,s2):
        y = min(s1.tag[1],s2.tag[1])
        return p['t1'] if y % 2 == 0 else p['t2']
    b[kwant.builder.HoppingKind((1,0),lat)] = -p['tx']
    b[kwant.builder.HoppingKind((0,1),lat)] = hy
    for direction in (-1,1):
        lead = kwant.Builder(kwant.TranslationalSymmetry((direction,0)))
        for y in range(2*p['cells']):
            lead[lat(0,y)] = p['delta']*(-1)**y
        lead[kwant.builder.HoppingKind((1,0),lat)] = -p['tx']
        lead[kwant.builder.HoppingKind((0,1),lat)] = hy
        b.attach_lead(lead)
    return b.finalized()


def analytic_contact_dos(p, energy):
    rows = 2*p['cells']
    hy = np.diag(p['delta']*(-1.)**np.arange(rows))
    for y in range(rows-1):
        hy[y,y+1] = hy[y+1,y] = p['t1'] if y % 2 == 0 else p['t2']
    eig, psi = np.linalg.eigh(hy)
    shifted = energy-eig
    active = abs(shifted)<2*p['tx']
    denom = (2*p['tx'])**2-shifted[active]**2
    dos = np.zeros(rows)
    derivative = np.zeros(rows)
    dos[active] = 1/(np.pi*np.sqrt(denom))
    derivative[active] = shifted[active]/(np.pi*denom**1.5)
    contact_weight = abs(psi[[-1,0]])**2
    rho = contact_weight@dos
    drho = contact_weight@derivative
    slope = drho[0]/rho[0]-drho[1]/rho[1]
    return rho, float(slope), int(active.sum())


def main():
    start = time.perf_counter()
    protected = [SOURCE]
    before = {relative(path):sha(path) for path in protected}
    with np.load(SOURCE, allow_pickle=False) as data:
        p = dict(tx=float(data['tx']), t1=float(data['t1']), t2=float(data['t2']),
            delta=float(data['delta']), EF=float(data['EF']), length=int(data['L']),
            cells=int(data['N']))
    couplings=np.array([.2,.15,.1,.05,.02,.01,.005])
    records=[]
    step_checks=[]
    weak_formula_reference=None
    step=1e-5
    with threadpool_limits(limits=1):
        bare=bare_system(p)
        xpos=p['length']//2
        contact_ids=[next(i for i,s in enumerate(bare.sites) if tuple(s.tag)==(xpos,y))
                     for y in (2*p['cells']-1,0)]
        bare_rho=kwant.ldos(bare,p['EF'],params={})[contact_ids]
        rho_plus=kwant.ldos(bare,p['EF']+step,params={})[contact_ids]
        rho_minus=kwant.ldos(bare,p['EF']-step,params={})[contact_ids]
        bare_slope=(np.log(rho_plus[0]/rho_plus[1])-np.log(rho_minus[0]/rho_minus[1]))/(2*step)
        rho_analytic,slope_analytic,nactive=analytic_contact_dos(p,p['EF'])
        bare_check=dict(kwant_ldos=bare_rho.tolist(),analytic_ldos=rho_analytic.tolist(),
            contact_ldos_max_error=float(np.max(abs(bare_rho-rho_analytic))),
            kwant_finite_difference_log_slope=float(bare_slope),
            analytic_log_slope=slope_analytic,log_slope_error=float(abs(bare_slope-slope_analytic)),
            propagating_transverse_subbands=nactive)
        if bare_check['contact_ldos_max_error']>1e-10 or bare_check['log_slope_error']>1e-7:
            raise RuntimeError('Bare local DOS validation failed.')
        dos_prediction=-slope_analytic/8
        for tc in couplings:
            dev=Device('four_x',cells=p['cells'],delta=p['delta'],tc=float(tc),
                tx=p['tx'],t1=p['t1'],t2=p['t2'],length=p['length'])
            record,_,field,_=floating_response(dev,p['EF'],step=step)
            # F_i sums source-to-probe transmissions. This checks the
            # frozen transmission formula independently of the floating solve.
            cplus=dev.matrix(p['EF']+step)
            cminus=dev.matrix(p['EF']-step)
            fplus=-cplus[2:,:2].sum(axis=1)
            fminus=-cminus[2:,:2].sum(axis=1)
            fslope=(np.log(fplus[0]/fplus[1])-np.log(fminus[0]/fminus[1]))/(2*step)
            c=dev.matrix(p['EF'])
            f=-c[2:,:2].sum(axis=1)
            report=dict(tc=float(tc), kH_raw=record['kH'],
                kH_frozen_raw=record['kH_frozen'], kH_internal_raw=record['kH_internal'],
                kH_electrodes_raw=record['kH_electrodes'],
                kappa2_full=-record['kH']/4,kappa2_frozen=-record['kH_frozen']/4,
                kappa2_internal=-record['kH_internal']/4,
                kappa2_electrodes=-record['kH_electrodes']/4,
                electrode_frozen_firstorder_prediction=-(record['kH_frozen']+record['kH_internal'])/4,
                frozen_F_prediction=float(-fslope/8),bare_DOS_prediction=dos_prediction,
                frozen_F_error=float(abs(-record['kH_frozen']/4+fslope/8)),
                frozen_DOS_error=float(abs(-record['kH_frozen']/4-dos_prediction)),
                split_sum_error=float(abs(record['kH']-record['kH_frozen']
                                          -record['kH_internal']-record['kH_electrodes'])/4),
                probe_F=f.tolist(), probe_F_over_tc_squared=(f/tc**2).tolist(),
                interprobe_transmission=float(-c[2,3]),
                alpha_H=record['alpha_H'], probe_condition=record['cpp_condition'],
                current_conservation_error=record['current_sum_error'],
                probe_quadratic_current_error=record['probe_quadratic_error'],
                zero_field_reciprocity_error=record['reciprocity_error'],
                gauge_shift_error=record['gauge_shift_error'])
            records.append(report)
            print(json.dumps(report),flush=True)
            if tc == couplings[-1]:
                alpha=np.asarray(record['alpha'])
                cv=(dev.matrix(p['EF'],step*field,step*alpha)
                    -dev.matrix(p['EF'],-step*field,-step*alpha))/(2*step)
                cs=(dev.matrix(p['EF'],step*field)
                    -dev.matrix(p['EF'],-step*field))/(2*step)
                cl=(dev.matrix(p['EF'],voltages=step*alpha)
                    -dev.matrix(p['EF'],voltages=-step*alpha))/(2*step)
                # t'_L-t'_R is differentiated with respect to the code
                # source-energy amplitude v, including its induced field.
                D=-cv[2:,0]+cv[2:,1]
                Ds=-cs[2:,0]+cs[2:,1]
                Dl=-cl[2:,0]+cl[2:,1]
                Fprime=(fplus-fminus)/(2*step)
                beta_formula=.5*Fprime/f+D/f
                beta_block=np.asarray(record['beta'])[2:]
                kappa_formula=-(beta_formula[0]-beta_formula[1])/4
                weak_formula_reference=dict(tc=float(tc), F=f.tolist(), Fprime=Fprime.tolist(),
                    tprime_L_minus_tprime_R_code=D.tolist(),
                    tprime_internal_code=Ds.tolist(), tprime_electrodes_code=Dl.tolist(),
                    bias_derivative_split_max_error=float(np.max(abs(D-Ds-Dl))),
                    beta_formula_code=beta_formula.tolist(), beta_full_block_code=beta_block.tolist(),
                    beta_formula_vs_block_max_error=float(np.max(abs(beta_formula-beta_block))),
                    kappa_formula=float(kappa_formula),kappa_block=report['kappa2_full'],
                    kappa_formula_vs_block_error=float(abs(kappa_formula-report['kappa2_full'])),
                    assumptions='Mx implies alpha_P=0 and TL=TR; neglect T34. T34 was measured, not set to zero.',
                    raw_formula='beta_p = Fprime_p/(2 F_p) + (partial_v T_pL - partial_v T_pR)/F_p',
                    published_formula='kappa2=-[beta_top-beta_bottom]/4 in e/t2 units; first term is -d_E ln(Ftop/Fbottom)/8')
                if weak_formula_reference['beta_formula_vs_block_max_error']>1e-7:
                    raise RuntimeError('General weak-probe formula disagrees with the block solve.')
            if tc in (.1,.02):
                refined=floating_response(dev,p['EF'],step=step/2)[0]
                errors={name:abs(record[name]-refined[name])/4 for name in
                    ['kH','kH_frozen','kH_internal','kH_electrodes']}
                step_checks.append(dict(tc=float(tc), absolute_kappa2_errors=errors))
                if max(errors.values())>2e-6:
                    raise RuntimeError('Derivative refinement failed.')
    arrays=dict(tc=couplings,EF=np.array(p['EF']),bare_contact_ldos=bare_rho,
        bare_contact_log_ratio_derivative_analytic=np.array(slope_analytic),
        bare_DOS_prediction=np.array(dos_prediction))
    for name in ['kH_raw','kH_frozen_raw','kH_internal_raw','kH_electrodes_raw',
        'kappa2_full','kappa2_frozen','kappa2_internal','kappa2_electrodes',
        'electrode_frozen_firstorder_prediction',
        'frozen_F_prediction','probe_F','probe_F_over_tc_squared','interprobe_transmission']:
        arrays[name]=np.array([record[name] for record in records])
    OUT.mkdir(parents=True,exist_ok=True)
    npz_path=OUT/'rm_tc_scan.npz'
    np.savez_compressed(npz_path,**arrays)
    after={relative(path):sha(path) for path in protected}
    assert before==after
    metadata=dict(schema='rm_weak_probe_dos_audit_v1',parameters=p,
        probe_couplings=couplings.tolist(),derivative_step=step,
        protocol='Historical positive-energy solver source shifts +v,-v; probe energy shifts float. Paper mu=EF-eV, V1=+V_SD/2, V2=-V_SD/2, V_perp=V3-V4.',
        conversion='kappa2_published=-kH_raw/4 exactly once',
        DOS_formula='kappa2_DOS=-1/8*d_E ln(rho_bare_top/rho_bare_bottom)',
        units='E in t2; published kappa2 in e/t2; DOS per t2 per contacted atomic site',
        frozen_definition='Freeze all sample potentials and all reservoir Hamiltonian onsite shifts; vary reservoir distribution functions only.',
        full_definition='Existing first-order local-neutrality sample potential plus reservoir Hamiltonian onsite shifts.',
        bare_definition='Same finite sample and matched L/R leads, with both metallic probes entirely removed.',
        bare_check=bare_check,records=records,step_refinement=step_checks,
        weak_formula_reference=weak_formula_reference,
        electrode_frozen_control='This uses the existing linear-response split: frozen occupation term plus internal-potential correction, excluding electrode Hamiltonian feedback; not a new nonlinear finite-bias calculation.',
        source_data_unchanged=True,protected_sha256=after,
        output_sha256=sha(npz_path),source_files={relative(path):sha(path) for path in
            [Path(__file__),RM_ROOT/'code/compute_initial_study.py',
             RM_ROOT/'code/compute_floating_study.py']},
        kwant_version=kwant.__version__,numpy_version=np.__version__,
        runtime_seconds=time.perf_counter()-start)
    write_report(OUT/'rm_tc_scan.json', metadata)
    print(json.dumps(dict(output=str(npz_path),bare_check=bare_check,
        bare_DOS_prediction=dos_prediction,step_checks=step_checks,
        weak_formula_reference=weak_formula_reference,
        runtime_seconds=metadata['runtime_seconds']),indent=2),flush=True)


# Use code/checks/rm_weak_probe.py; no production data are overwritten.
