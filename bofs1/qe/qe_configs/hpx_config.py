hpx_config = {
    'command': ['./bofs1_env/bin/mpirun',
                '--allow-run-as-root',
                '--use-hwthread-cpus',
                '-x', 'OMP_NUM_THREADS=1',
                '-x', 'BLIS_NUM_THREADS=1',
                '-x', 'FLAME_NUM_THREADS=1',
                '-np', '32', './bin/hp.x'],
    'qpts_q_spacing': 0.065,        # q-point spacing
    'inputhp': {
        'iverbosity': 2,           # Verbosity level
        'find_atpert': 4,          # Method for determining which atoms to perturb (4=Perturb all Hubbard atoms)
        'skip_equivalence_q': True, # If True, use full q-point grid; if False, use symmetry to reduce q-points
        'num_neigh': 16,            # Number of nearest neighbors for Hubbard V parameters
        'lmin': 0,                 # Minimum orbital quantum number for Hubbard V parameters
        'rmax': 20,             # Maximum distance to search for neighbors
    }
}
