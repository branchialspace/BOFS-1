hpx_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/hp.x'],
    'qpts_q_spacing': 0.065,        # q-point spacing
    'inputhp': {
        'iverbosity': 2,           # Verbosity level
        'find_atpert': 4,          # Method for determining which atoms to perturb (4=Perturb all Hubbard atoms)
        'docc_thr': 5.0e-5,        # Threshold for comparing unperturbed occupations
        'skip_equivalence_q': True, # If True, use full q-point grid; if False, use symmetry to reduce q-points
        'conv_thr_chi': 1.0e-5,    # Convergence threshold for response function chi
        'thresh_init': 1.0e-14,    # Initial threshold for linear system solution
        'ethr_nscf': 1.0e-11,      # Threshold for NSCF eigenvalue convergence
        'niter_max': 100,          # Maximum number of iterations
        'alpha_mix': 0.3,          # Mixing parameter for updating SCF potential
        'nmix': 4,                 # Number of iterations used in Broyden potential mixing
        'compute_hp': False,       # If True, post-process previously calculated results
        'determine_num_pert_only': False, # If True, determine only the number of perturbations and exit
        'determine_q_mesh_only': False,   # If True, determine only the q mesh and exit
        'sum_pertq': False,        # If True, collect response matrices for all q points
        'max_seconds': 86400,      # Maximum allowed run time in seconds
        'num_neigh': 20,            # Number of nearest neighbors for Hubbard V parameters
        'lmin': 1,                 # Minimum orbital quantum number for Hubbard V parameters
        'rmax': 230.0,             # Maximum distance to search for neighbors
        'dist_thr': 2.0e-3,        # Threshold for comparing inter-atomic distances
    }
}
