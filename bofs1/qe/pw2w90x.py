# Run QuantumESPRESSO pw2wannier90.x

import re
import subprocess
from subprocess import CalledProcessError


def pw2w90x(
        structure_name,
        run_name,
        config
):
    """
    Run QE pw2wannier90.x interface calculation.
    structure_name : string
        Name of the structure (used for prefix/outdir matching).
    run_name : string
        The filename of the pw.x calculation.
        used to parse Fermi energy for SCDM and locate the outdir.
    config : dict
        Configuration dictionary containing required settings.
    """

    def get_fermi_energy(pwo_path):
        """
        Extract Fermi energy from the .pwo file for SCDM mu.
        """
        try:
            with open(pwo_path, 'r') as f:
                lines = f.readlines()
                content = "".join(lines)
            fermi_match = re.search(r'the Fermi energy is\s+([-\d.]+)\s+ev', content, re.IGNORECASE)
            return float(fermi_match.group(1))

    def write_pw2w90_input(
            config,
            input_filename
    ):
        """
        Write the pw2wannier90.x input file.
        """
        with open(input_filename, 'w') as f:
            f.write('&INPUTPP\n')
            outdir = config['control'].get('outdir', f"./{structure_name}")
            prefix = config['control'].get('prefix', structure_name)
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
    run_name = f"{structure_name}_pw2wan"
    pwo_path = f"{run_name}.pwo"
    e_fermi = get_fermi_energy(pwo_path)
    print(f"Detected Fermi energy from {pwo_path}: {e_fermi} eV. Setting scdm_mu.")
    config['inputpp']['scdm_mu'] = e_fermi
    # Write input file
    write_pw2w90_input(config, f"{run_name}.win")
    # Subprocess run
    try:
        with open(f"{run_name}.wout", 'w') as f_out:
            command_list = config['command'] + ['-in', f"{run_name}.win"]
            subprocess.run(
                command_list,
                stdout=f_out,
                stderr=subprocess.STDOUT,
                check=True)
        print("pw2wannier90 calculation completed successfully.")
    except CalledProcessError as cpe:
        print(f"Error running pw2wannier90: {cpe}")
        try:
            with open(f"{run_name}.wout", 'r') as f_out:
                print("\npw2wannier90 Output:")
                print(f_out.read())
        except Exception as e:
            print(f"Could not read output file: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")
