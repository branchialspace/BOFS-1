pw2w90x_config = {
    'command': ['./bofs1_env/bin/mpirun',
                '--allow-run-as-root',
                '--use-hwthread-cpus',
                '-x', 'OMP_NUM_THREADS=1',
                '-x', 'BLIS_NUM_THREADS=1',
                '-x', 'FLAME_NUM_THREADS=1',
                '-np', '32',
                "pw2wannier90.x"],
    "inputpp": {
        "seedname": "wannier",
        "write_mmn": True,
        "write_amn": True,
        "write_unk": False,
        "scdm_proj": True, # necessary
        "scdm_entanglement": "erfc",
        "scdm_sigma": 0.2,
        # 'scdm_mu' is .pwo fermi energy
    }
}
