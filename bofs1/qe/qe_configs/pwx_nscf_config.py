pwx_nscf_config = {
    'command': ['./bofs1_env/bin/mpirun',
                '--allow-run-as-root',
                '--use-hwthread-cpus',
                '-x', 'OMP_NUM_THREADS=1',
                '-x', 'BLIS_NUM_THREADS=1',
                '-x', 'FLAME_NUM_THREADS=1',
                '-np', '32',
                './qe-7.5/bin/pw.x'],
    'wfn_scalar': 1.15,
    'rho_scalar': 1.15,
    'kpts_k_minimum': 6,
    'kpts_k_spacing': 0.05,
    'kpts_shift': (1,1,1),
    'nbnd_scalar': 2.3,
    'initial_u_value': "off",
    "magnetization": {
        "d_3d": 0.6,
        "d_4d5d": 0.3,
        "f": 0.8
    },
    'control': {
        'calculation': 'nscf',
        'restart_mode': 'from_scratch',
        'pseudo_dir': './pslibrary/rel-pbe/PSEUDOPOTENTIALS',    # ./ONCVPseudoPack/PseudoDojo/FR_v0.4/PBE_stringent   ./pslibrary/rel-pbe/PSEUDOPOTENTIALS 
        'disk_io': 'low',
        'verbosity': 'high',
        'tprnfor': True,
        'tstress': False
    },
    'system': {
        'input_dft': 'pbe',
        'vdw_corr': 'mbd',
        'nosym': 'true',
        'noinv': 'true',
        'ibrav': 0,
        'occupations': 'smearing',
        'smearing': 'fermi-dirac',
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
