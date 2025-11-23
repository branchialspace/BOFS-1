hpx_config = {
    'command': ['./bofs1_env/bin/mpirun',
                '--allow-run-as-root',
                '--use-hwthread-cpus',
                '-x', 'OMP_NUM_THREADS=1',
                '-x', 'BLIS_NUM_THREADS=1',
                '-x', 'FLAME_NUM_THREADS=1',
                '-np', '32',
                './qe-7.5/bin/hp.x'],
    'inputhp': {
        'iverbosity': 2,           # Verbosity level
        'find_atpert': 1,          # Method for determining which atoms to perturb (4=Perturb all Hubbard atoms)
        'skip_equivalence_q': False, # If True, use full q-point grid; if False, use symmetry to reduce q-points
        'num_neigh': 30,            # Number of nearest neighbors for Hubbard V parameters
        'lmin': 0,                 # Minimum orbital quantum number for Hubbard V parameters
    }
}
