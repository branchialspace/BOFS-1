pwx_scf_config = {
    'command': ['./bofs1_env/bin/mpirun',
                '--allow-run-as-root',
                '-x', 'OMP_NUM_THREADS=1',
                '-x', 'BLIS_NUM_THREADS=1',
                '-x', 'FLAME_NUM_THREADS=1',
                '-np', '1',
                './bin/pw.x'],
    'wfn_scalar': 1.1,
    'rho_scalar': 1.1,
    'kpts_k_spacing': 0.13, # scf: 0.13    nscf: 0.09
    'kpts_shift': (1,1,1),
    'nbnd_scalar': 1.8,
    'initial_u_value': 0.1,
    "magnetization": {
        "d_3d": 0.6,
        "d_4d5d": 0.3,
        "f": 0.8
    },
    'control': {
        'calculation': 'scf', # scf     nscf     bands
        'restart_mode': 'from_scratch',
        'pseudo_dir': './pslibrary/rel-pbe/PSEUDOPOTENTIALS', # /content/ONCVPseudoPack/Abinit_v0.4/UPF/PBEsol   /content/pslibrary/rel-pbe/PSEUDOPOTENTIALS
        'disk_io': 'high',
        'tprnfor': True,
        'tstress': True
    },
    'system': {
        'input_dft': 'pbe',
        'ibrav': 0,
        'occupations': 'smearing',
        'smearing': 'marzari-vanderbilt', # gaussian     marzari-vanderbilt     fermi-dirac
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
