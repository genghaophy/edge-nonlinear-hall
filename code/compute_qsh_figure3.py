"""Verified selected-point QSH validators, extracted from the production driver.

Functions preserve the source algorithms. Use code/checks/qsh_reference.py;
this module never launches a production scan or writes production caches.
"""
import numpy as np
from qsh_environment_response import EnvironmentDevice, quadrature, thermal_response

def validate_device(device, energy=.4, finite_bias=True):
    report, arrays = device.response_zero(energy)
    refined = device.response_zero(energy, step=5e-6)[0]
    rng = np.random.default_rng(321)
    direction = rng.normal(size=device.nsites)
    step = 2e-6
    derivative = ((device.effective(energy, step*direction)
                   -device.effective(energy, -step*direction))@np.asarray(report['alpha']))/(2*step)
    predicted = arrays['D']@direction
    check = dict(kH=report['kH'], derivative_step_error=abs(refined['kH']-report['kH']),
        kernel_absolute_error=float(np.max(abs(derivative-predicted))),
        trs_h_error=device.trs_h_error, equilibrium=report['checks'],
        virtual_spin_channels=len(device.virtual_spins),
        zero_voltage_current=float(np.max(abs(device.terminal_current(
            energy, np.zeros(4), arrays['u'])))))
    if finite_bias:
        v=5e-5
        vp, ip=device.solve_floating(energy, v, report, arrays['u'])
        vm, im=device.solve_floating(energy, -v, report, arrays['u'])
        direct=((vp[2]-vp[3])+(vm[2]-vm[3]))/(2*v*v)
        check.update(finite_bias_kH=float(direct), finite_bias_error=abs(direct-report['kH']),
            probe_current_error=float(max(np.max(abs(ip[2:])),np.max(abs(im[2:])))),
            finite_bias_current_sum_error=float(max(abs(ip.sum()),abs(im.sum()))),
            bias_amplitude=v, terminal_voltages=[vp.tolist(),vm.tolist()])
        if check['finite_bias_error']>max(3e-5,.003*abs(report['kH'])):
            raise RuntimeError(check)
    if check['kernel_absolute_error']>3e-6 or check['derivative_step_error']>2e-5:
        raise RuntimeError(check)
    return check


def validate_thermal():
    """Independent finite-Fermi-window zero-current solve on a small device."""
    device=EnvironmentDevice(length=9,width=6)
    energy,temperature=.4,.02
    nodes,weights=quadrature(.04,.76,.04,8)
    points=[device.spectral_record(float(e)) for e in nodes]
    spectrum={k:np.asarray([p[k] for p in points])
              for k in ('C','rho','D','qlead','qwindow')}
    report,u=thermal_response(device,nodes,weights,spectrum,energy,temperature)
    errors=[];direct=[];residual=0.
    for v in (.001,.0005):
        vp,ip=device.solve_floating(energy,v,report,u,temperature,nodes,weights)
        vm,im=device.solve_floating(energy,-v,report,u,temperature,nodes,weights)
        value=((vp[2]-vp[3])+(vm[2]-vm[3]))/(2*v*v)
        errors.append(abs(value-report['kH']));direct.append(float(value))
        residual=max(residual,float(np.max(abs(ip[2:]))),float(np.max(abs(im[2:]))))
    check=dict(length=9,width=6,energy=energy,temperature=temperature,
               kH=report['kH'],finite_bias_kH=direct,finite_bias_errors=errors,
               biases=[.001,.0005],probe_current_error=residual,
               thermal_mass=report['thermal_mass'],window_ibp_error=report['window_ibp_error'],
               nodes=len(nodes))
    if errors[-1]>2e-5 or residual>2e-11:
        raise RuntimeError(check)
    return check

