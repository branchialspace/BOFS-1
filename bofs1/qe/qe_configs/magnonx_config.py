magnonx_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/turbo_magnon.x'],
    'q_vector': [0.0, 0.0, 0.0],    # q-vector for magnon calculations in units of 2pi/a
    'lr_control': {
        'itermax': 500,             # Number of Lanczos iterations
        'pseudo_hermitian': True,   # Use pseudo-Hermitian Lanczos algorithm
        'approximation': 'TDDFT',   # Level of theory: 'TDDFT' or 'IPA'
        'ipol': 4,                  # Polarization: 1-3 for specific component, 4 for full tensor
    }
}
