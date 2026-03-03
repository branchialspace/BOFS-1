projwfcx_config = {
    'command': ['./bofs1_env/bin/mpirun',
                '--allow-run-as-root',
                '--use-hwthread-cpus',
                '-x', 'OMP_NUM_THREADS=1',
                '-x', 'BLIS_NUM_THREADS=1',
                '-x', 'FLAME_NUM_THREADS=1',
                '-np', '32',
                'projwfc.x'],
    'projwfc': {
        'ngauss': 0,
        'degauss': 0.01,
        'DeltaE': 0.05,
        'lsym': False,
        'diag_basis': False,
        'pawproj': False,
        'lwrite_overlaps': False,
        'kresolveddos': False,
    }
}
