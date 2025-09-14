# Run QuantumESPRESSO projwfc.x

import os
import subprocess
from subprocess import CalledProcessError


def projwfcx(
    structure_path,
    config
):
    """
    Run QE projwfc.x calculation to compute projected wavefunctions,
    Lowdin charges, spilling parameter, and projected DOS.
    structure_path : string
        Path of the structure file that was used in the previous pw.x calculation.
        Used to derive the structure name, which must match the prefix used in pw.x.
    config : dict
        Configuration dictionary containing required settings.
    """
    def write_projwfcx_input(config, input_filename):
        """
        Write the QE projwfc.x input file from config settings.
        """
        with open(input_filename, 'w') as f:
            # PROJWFC namelist
            f.write('&PROJWFC\n')
            for key, value in config['projwfc'].items():
                if isinstance(value, bool):
                    val = '.true.' if value else '.false.'
                elif isinstance(value, str):
                    val = f"'{value}'"
                else:
                    val = value
                f.write(f"  {key} = {val}\n")
            f.write('/\n')
    
    # Args
    structure_name = os.path.splitext(os.path.basename(structure_path))[0]
    calculation = "projwfc"
    run_name = f"{structure_name}_{calculation}"
    command = config['command']
    config['projwfc']['prefix'] = structure_name
    config['projwfc']['outdir'] = structure_name
    config['projwfc']['filpdos'] = f"{structure_name}"
    # Write QE projwfc.x input file
    write_projwfcx_input(config, f"{run_name}.wfci")
    # Subprocess run
    try:
        with open(f"{run_name}.wfco", 'w') as f_out:
            command_list = command + ['-in', f"{run_name}.wfci"]
            subprocess.run(
                command_list,
                stdout=f_out,
                stderr=subprocess.STDOUT,
                check=True)
        print("PROJWFC calculation completed successfully.")
    except CalledProcessError as cpe:
        print(f"Error running PROJWFC calculation: {cpe}")
        try:
            with open(f"{run_name}.wfco", 'r') as f_out:
                print("\nPROJWFC Output:")
                print(f_out.read())
        except Exception as e:
            print(f"Could not read output file: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")


wfcx_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/projwfc.x'],
    'projwfc': {
        'ngauss': -99,                # Type of gaussian broadening (-99 = Fermi-Dirac)
        'degauss': 0.01,              # Gaussian broadening in Ry
        'DeltaE': 0.05,               # Energy grid step in eV
        'lsym': False,                # Symmetrize projections
        'diag_basis': False,           # Project onto global or local XYZ frame
        'pawproj': False,             # Use PAW projectors (for PAW pseudos only)
        'lwrite_overlaps': False,     # Orbital overlap matrix (for parallel, paste , '-nd', '1' to end of command)
        'kresolveddos': True,         # Sum over all k-points (not k-resolved)
    }
}

# ASE structure
mof = "/content/mofs/SIWZOO_full_n2.cif"
# Run projwfc.x
qe_projwfcx(mof, wfcx_config)
