lanczosx_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/turbo_lanczos.x'],
    'lr_control': {
        'itermax': 500,          # Number of Lanczos iterations
        'ipol': 4,               # Polarization: 1-3 for specific component, 4 for full tensor
        'n_ipol': 3,             # Number of zeta coefficients to calculate
        'ltammd': False,         # Tamm-Dancoff approximation
        'no_hxc': False,         # Independent electron approximation
        'lrpa': False,           # Random Phase Approximation (no XC)
        'scissor': 0.0,          # Scissor shift in Rydberg units
        'charge_response': 0,    # Set to 1 to compute charge density response
        'pseudo_hermitian': True,  # Use pseudo-Hermitian Lanczos algorithm
        'd0psi_rs': False,       # Dipole computed in real space
        'lshift_d0psi': True     # Shift position operator for periodicity
    },
    'lr_post': {
        'omeg': 0.0,              # Transition energy in Rydberg units
        'epsil': 0.02,            # Broadening/damping term in Rydberg units
        'beta_gamma_z_prefix': 'pwscf',  # Prefix for beta gamma zeta coefficients
        'w_T_npol': 3,            # Number of polarization directions in previous calc
        'plot_type': 1            # Format for charge density response output
    }
}
