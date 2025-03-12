# Run QuantumESPRESSO ph.x

import os
import subprocess
from subprocess import CalledProcessError


def qe_phx(
    structure_path,
    config
):
    """
    Run QE ph.x calculation for phonon properties.
    structure_path : string
        Path of the structure file that was used in the previous pw.x calculation.
        Used to derive the structure name, which must match the prefix used in pw.x.
    config : dict
        Configuration dictionary containing required settings.
    """
    def write_phx_input(config, input_filename):
        """
        Write the QE ph.x input file from config settings.
        """
        with open(input_filename, 'w') as f:
            # INPUTPH namelist
            f.write('&INPUTPH\n')
            for key, value in config['inputph'].items():
                if isinstance(value, bool):
                    val = '.true.' if value else '.false.'
                elif isinstance(value, str):
                    val = f"'{value}'"
                else:
                    val = value
                f.write(f"  {key} = {val}\n")
            f.write('/\n')
            # Q-point specifications if not using ldisp or qplot
            if not config['inputph'].get('ldisp', False) and not config['inputph'].get('qplot', False):
                if 'xq' in config:
                    f.write(f"\n{config['xq'][0]:.8f} {config['xq'][1]:.8f} {config['xq'][2]:.8f}\n")
            # Q-points for qplot
            elif config['inputph'].get('qplot', False):
                if 'qpoints' in config:
                    f.write(f"\n{len(config['qpoints'])}\n")
                    for q_point in config['qpoints']:
                        f.write(f"{q_point[0]:.8f} {q_point[1]:.8f} {q_point[2]:.8f} {q_point[3]}\n")
            # Atom selection if nat_todo is specified
            if 'nat_todo' in config['inputph'] and config['inputph']['nat_todo'] > 0:
                if 'atoms' in config:
                    f.write('\n')
                    f.write(' '.join(str(atom) for atom in config['atoms']))
                    f.write('\n')
    
    # Args
    structure_name = os.path.splitext(os.path.basename(structure_path))[0]
    calculation = "ph"
    run_name = f"{structure_name}_{calculation}"
    command = config['command']
    config['inputph']['prefix'] = structure_name
    config['inputph']['outdir'] = structure_name
    os.makedirs(structure_name, exist_ok=True)
    # Write QE ph.x input file
    write_phx_input(config, f"{run_name}.phi")
    # Subprocess run
    try:
        with open(f"{run_name}.pho", 'w') as f_out:
            command_list = command + ['-in', f"{run_name}.phi"]
            subprocess.run(
                command_list,
                stdout=f_out,
                stderr=subprocess.STDOUT,
                check=True)
        print("Phonon calculation completed successfully.")
    except CalledProcessError as cpe:
        print(f"Error running phonon calculation: {cpe}")
        try:
            with open(f"{run_name}.pho", 'r') as f_out:
                print("\nPhonon Output:")
                print(f_out.read())
        except Exception as e:
            print(f"Could not read output file: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

phx_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/ph.x'],
    'xq': [0.0, 0.0, 0.0],         # q-point for non-ldisp calculations
    'inputph': {
        'tr2_ph': 1.0e-14,         # Convergence threshold for phonons
        'ldisp': True,             # Run phonons on a grid of q-points
        'epsil': True,             # Calculate dielectric constant
        'trans': True,             # Calculate phonons
        'electron_phonon': '',     # electron-phonon coefficient method
        'fildyn': 'dynmat',        # Prefix for dynamical matrices
        'fildrho': 'drho',         # File for charge density response
        'fildvscf': 'dvscf',       # File for potential variation
        'max_seconds': 86400,      # Maximum allowed run time in seconds
        'asr': True,               # Apply Acoustic Sum Rule
        'search_sym': True,        # Enable mode symmetry analysis
    }
}

# ASE structure
mof = "/content/mofs/SIWZOO_full_n2.cif"
# Run ph.x
qe_phx(mof, phx_config)
