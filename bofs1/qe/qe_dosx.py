# Run QuantumESPRESSO dos.x

import os
import subprocess
from subprocess import CalledProcessError


def qe_dosx(
    structure_path,
    config
):
    """
    Run QE dos.x calculation to compute Density of States.
    structure_path : string
        Path of the structure file that was used in the previous pw.x calculation.
        Used to derive the structure name, which must match the prefix used in pw.x.
    config : dict
        Configuration dictionary containing required settings.
    """
    def write_dosx_input(config, input_filename):
        """
        Write the QE dos.x input file from config settings.
        """
        with open(input_filename, 'w') as f:
            # DOS namelist
            f.write('&DOS\n')
            for key, value in config['dos'].items():
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
    calculation = "dos"
    run_name = f"{structure_name}_{calculation}"
    command = config['command']
    config['dos']['prefix'] = structure_name
    config['dos']['outdir'] = structure_name
    config['dos']['fildos'] = f"{structure_name}.dos"
    # Write QE dos.x input file
    write_dosx_input(config, f"{run_name}.dosi")
    # Subprocess run
    try:
        with open(f"{run_name}.doso", 'w') as f_out:
            command_list = command + ['-in', f"{run_name}.dosi"]
            subprocess.run(
                command_list,
                stdout=f_out,
                stderr=subprocess.STDOUT,
                check=True)
        print("DOS calculation completed successfully.")
    except CalledProcessError as cpe:
        print(f"Error running DOS calculation: {cpe}")
        try:
            with open(f"{run_name}.doso", 'r') as f_out:
                print("\nDOS Output:")
                print(f_out.read())
        except Exception as e:
            print(f"Could not read output file: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")


dosx_config = {
    'command': ['/usr/bin/mpirun', '--allow-run-as-root', '-x', 'OMP_NUM_THREADS=2', '-np', '4', '/content/bin/dos.x'],
    'dos': {
        'bz_sum': 'smearing',        # Method for Brillouin zone summation
        'ngauss': -99,                 # Type of gaussian broadening (-99 = Fermi-Dirac)
        'degauss': 0.01,             # Gaussian broadening in Ry
        'DeltaE': 0.05,              # Energy grid step in eV
    }
}

# ASE structure
mof = "/content/mofs/SIWZOO_full_n2.cif"
# Run dosx
qe_dosx(mof, dosx_config)
