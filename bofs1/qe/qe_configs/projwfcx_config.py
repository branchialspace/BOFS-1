projwfcx_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/projwfc.x'],
    'projwfc': {
        'ngauss': -99,                # Type of gaussian broadening (-99 = Fermi-Dirac)
        'degauss': 0.01,              # Gaussian broadening in Ry
        'DeltaE': 0.05,               # Energy grid step in eV
        'lsym': False,                # Symmetrize projections
        'diag_basis': False,           # Project onto global or local XYZ frame
        'pawproj': False,             # Use PAW projectors (for PAW pseudos only)
        'lwrite_overlaps': False,     # Orbital overlap matrix (for parallel, paste , '-nd', '1' to end of command)
        'kresolveddos': True,         # Sum over all k-points (not k-resolved)
    }
}
