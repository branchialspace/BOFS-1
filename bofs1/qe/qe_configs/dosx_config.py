dosx_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/dos.x'],
    'dos': {
        'bz_sum': 'smearing',        # Method for Brillouin zone summation
        'ngauss': -99,                 # Type of gaussian broadening (-99 = Fermi-Dirac)
        'degauss': 0.01,             # Gaussian broadening in Ry
        'DeltaE': 0.05,              # Energy grid step in eV
    }
}
