pwx_nscf_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/pw.x'],
    'wfn_scalar': 1.15,
    'rho_scalar': 1.15,
    'kpts_k_spacing': 0.09, # scf: 0.13    nscf: 0.09
    'kpts_shift': (1,1,1),
    'initial_u_value': 0.1,
    'nbnd_scalar': 2.5,
    'n_manifolds': 1,
    'control': {
        'calculation': 'nscf', # scf     nscf     bands
        'restart_mode': 'from_scratch',
        'pseudo_dir': '/content/pslibrary/rel-pbe/PSEUDOPOTENTIALS',
        'disk_io': 'medium',
        'verbosity': 'high',
        'wf_collect': True,
        'tprnfor': True,
        'tstress': True
    },
    'system': {
        'ibrav': 0,
        'occupations': 'smearing',
        'smearing': 'gaussian', # gaussian     marzari-vanderbilt     fermi-dirac
        'degauss': 0.01,
        'noncolin': True,
        'lspinorb': True
    },
    'electrons': {
        'conv_thr': 1.0e-6,
        'mixing_beta': 0.3,
        'electron_maxstep': 300,
    }
}
