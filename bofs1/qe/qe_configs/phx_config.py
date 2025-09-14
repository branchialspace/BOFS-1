phx_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/ph.x'],
    'xq': [0.0, 0.0, 0.0],         # q-point for non-ldisp calculations
    'qpts_q_spacing': 0.25,        # q-point spacing for automatic grid generation (coarser than k-points)
    'inputph': {
        'tr2_ph': 1.0e-14,         # Convergence threshold for phonons
        'ldisp': True,             # Run phonons on a grid of q-points
        'epsil': False,            # Calculate dielectric constant
        'trans': True,             # Calculate phonons
        'elop': True,              # Calculate electro-optic tensor
        'electron_phonon': '',     # electron-phonon coefficient method
        'fildyn': 'dynmat',        # Prefix for dynamical matrices
        'fildrho': 'drho',         # File for charge density response
        'fildvscf': 'dvscf',       # File for potential variation
        'max_seconds': 86400,      # Maximum allowed run time in seconds
        'asr': True,               # Apply Acoustic Sum Rule
        'search_sym': True,        # Enable mode symmetry analysis
    }
}
