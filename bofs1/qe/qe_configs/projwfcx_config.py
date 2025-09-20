projwfcx_config = {
    'command': ['./bofs1_env/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', './bin/projwfc.x'],
    'projwfc': {
        'ngauss': 0,                # Type of gaussian broadening (0 = Gaussian, -99 = Fermi-Dirac)
        'degauss': 0.01,              # Gaussian broadening in Ry
        'DeltaE': 0.05,               # Energy grid step in eV
        'lsym': False,                # Symmetrize projections
        'diag_basis': False,           # Project onto global or local XYZ frame
        'pawproj': False,             # Use PAW projectors (for PAW pseudos only, not implemented in noncolinear case)
        'lwrite_overlaps': False,     # Orbital overlap matrix (for parallel, paste , '-nd', '1' to end of command)
        'kresolveddos': True,         # Sum over all k-points (not k-resolved)
    }
}
