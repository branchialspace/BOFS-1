eelsx_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/turbo_eels.x'],
    'q_vector': [0.0, 0.0, 0.0],    # q-vector for EELS calculations in units of 2pi/a
    'lr_control': {
        'approximation': 'TDDFT',   # Level of theory: 'TDDFT', 'IPA', or 'RPA_with_CLFE'
        'calculator': 'lanczos',    # Algorithm: 'lanczos' or 'sternheimer'
        'itermax': 500,             # Number of iterations (Lanczos or Sternheimer)
        'pseudo_hermitian': True,   # Use pseudo-Hermitian Lanczos algorithm
        'ethr_nscf': 1.0e-11        # Threshold for convergence of eigenvalues
        # The following parameters are used when calculator = 'sternheimer'
        # 'alpha_mix': 0.7,         # Mixing parameter for SCF potential
        # 'epsil': 0.02,            # Broadening/damping term in Rydberg units
        # 'units': 1,               # Units: 0=Rydbergs, 1=Electron volts
        # 'start': 0.0,             # Starting frequency
        # 'end': 2.5,               # Ending frequency
        # 'increment': 0.001        # Frequency step
    }
}
