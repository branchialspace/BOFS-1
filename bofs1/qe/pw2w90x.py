# Run QuantumESPRESSO pw2wannier90.x

import re
import os
import subprocess
from subprocess import CalledProcessError


def pw2w90x(
        structure_path,
        config
):
    """
    Run QE pw2wannier90.x interface calculation.
    structure_path : string
        Path of the structure file that was used in the previous pw.x calculation.
        Used to derive the structure name, which must match the prefix used in pw.x.
    config : dict
        Configuration dictionary containing required settings.
    """

    def get_fermi_energy(pwo_path):
        """
        Extract Fermi energy from the .pwo file for SCDM mu.
        """
        with open(pwo_path, 'r') as f:
            lines = f.readlines()
            content = "".join(lines)
        fermi_match = re.search(r'the Fermi energy is\s+([-\d.]+)\s+ev', content, re.IGNORECASE)
        return float(fermi_match.group(1))

    def write_pw2w90_input(config, input_filename):
        """
        Write the pw2wannier90.x input file.
        """
        with open(input_filename, 'w') as f:
            f.write('&INPUTPP\n')
            outdir = config['inputpp'].get('outdir', f"./{structure_name}")
            prefix = config['inputpp'].get('prefix', structure_name)
            f.write(f"  outdir = '{outdir}'\n")
            f.write(f"  prefix = '{prefix}'\n")
            for key, value in config['inputpp'].items():
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
    config['inputpp']['seedname'] = structure_name
    pwo_path = config.get('nscf_output', f"{structure_name}_nscf.pwo")
    e_fermi = get_fermi_energy(pwo_path)
    print(f"Detected Fermi energy from {pwo_path}: {e_fermi} eV. Setting scdm_mu.")
    config['inputpp']['scdm_mu'] = e_fermi
    # Write input file
    write_pw2w90_input(config, f"{structure_name}.pw2win")
    # Subprocess run
    try:
        with open(f"{structure_name}.pw2wout", 'w') as f_out:
            command_list = config['command'] + ['-in', f"{structure_name}.pw2win"]
            subprocess.run(
                command_list,
                stdout=f_out,
                stderr=subprocess.STDOUT,
                check=True)
        print("pw2wannier90 calculation completed successfully.")
    except CalledProcessError as cpe:
        print(f"Error running pw2wannier90: {cpe}")
        try:
            with open(f"{structure_name}.pw2wout", 'r') as f_out:
                print("\npw2wannier90 Output:")
                print(f_out.read())
        except Exception as e:
            print(f"Could not read output file: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")
