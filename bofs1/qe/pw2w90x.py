# Run QuantumESPRESSO pw2wannier90.x

import re
import os
import subprocess
from subprocess import CalledProcessError
from .scdm_fit import scdm_fit


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
        Extract Fermi energy from the .pwo file for SCDM mu fallback.
        """
        with open(pwo_path, 'r') as f:
            content = f.read()
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
    pwo_path = config.get('nscf_output', f"{structure_name}_nscf.pwo")
    config['inputpp'].setdefault('seedname', structure_name)

    # Derive SCDM mu/sigma from projectability fitting (Vitale et al.)
    atomic_proj_xml = f"{structure_name}/{structure_name}.save/atomic_proj.xml"
    try:
        scdm_params = scdm_fit(structure_name, atomic_proj_xml)
        config['inputpp']['scdm_mu'] = scdm_params['scdm_mu']
        config['inputpp']['scdm_sigma'] = scdm_params['scdm_sigma']
        print(f"Using fitted SCDM parameters: mu={scdm_params['scdm_mu']:.4f}, sigma={scdm_params['scdm_sigma']:.4f}")
    except Exception as e:
        print(f"Projectability fitting unavailable ({e}), falling back to Fermi energy for scdm_mu.")
        e_fermi = get_fermi_energy(pwo_path)
        config['inputpp']['scdm_mu'] = e_fermi
        print(f"Detected Fermi energy from {pwo_path}: {e_fermi} eV. Setting scdm_mu.")
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
