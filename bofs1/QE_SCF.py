# QuantumESPRESSO Single-Point SCF with SOC

import os
from pathlib import Path
import subprocess
from subprocess import CalledProcessError
from ase.data import atomic_masses, atomic_numbers
from ase.build import bulk


# ASE structure
bi = bulk('Bi', 'rhombohedral', a=4.75, c=12.36, orthorhombic=False)

def run_qe(
    structure,
    config,
    input_data,
    kpts=(4, 4, 4),
    input_filename='espresso.pwi',
    output_filename='espresso.pwo'
):
    """
    Run a Quantum ESPRESSO calculation for a given ASE structure, configuration parameters,
    and QE input data dictionary.

    :param structure: ASE Atoms object defining the system.
    :param config: Dictionary containing at least:
                   {
                     'command': path to pw.x or mpirun + pw.x,
                     'pseudo_dir': directory containing pseudopotential files,
                     'outdir': output directory for QE,
                   }
    :param input_data: Dictionary containing 'control', 'system', and 'electrons' keys
                       needed to write the QE input file.
    :param kpts: Tuple of (kx, ky, kz) for automatic K_POINTS.
    :param input_filename: Name of the input file to be generated.
    :param output_filename: Name of the output file to be generated by QE.
    """

    # Unpack config
    command = config['command']
    pseudo_dir = config['pseudo_dir']
    outdir = config['outdir']
    os.makedirs(outdir, exist_ok=True)

    # Pseudopotentials
    def find_pseudopotentials(ase_structure, pseudo_directory):
        """
        Determine the appropriate pseudopotential file for each atomic species
        from the given directory. Tries to pick up 'fr' (fully relativistic) files
        if they are present; otherwise, picks the first found.
        """
        species = set(atom.symbol for atom in ase_structure)
        pseudopot_dict = {}
        pseudo_path = Path(pseudo_directory)
        for symbol in species:
            candidates = []
            for pp_file in pseudo_path.glob('*.[uU][pP][fF]'):
                pp_stem = pp_file.stem.split('_')[0].split('.')[0]
                if pp_stem.lower() == symbol.lower():
                    candidates.append(pp_file)
            if not candidates:
                raise FileNotFoundError(f"No pseudopotential found for {symbol}")
            fr_files = [c for c in candidates if 'fr' in c.stem.lower()]
            selected = max(fr_files, key=lambda x: len(x.stem)) if fr_files else candidates[0]
            pseudopot_dict[symbol] = selected.name

        return pseudopot_dict

    pseudopotentials = find_pseudopotentials(structure, pseudo_dir)

    # Write QE input file
    def write_espresso_input(
        ase_structure,
        qe_input_data,
        pps,
        kpoints,
        filename
    ):
        """
        Write the Quantum ESPRESSO input file from an ASE structure, user input data,
        and selected pseudopotentials.
        """
        with open(filename, 'w') as f:
            # Control section
            f.write('&control\n')
            for key, value in qe_input_data['control'].items():
                if isinstance(value, bool):
                    val = '.true.' if value else '.false.'
                elif isinstance(value, str):
                    val = f"'{value}'"
                else:
                    val = value
                f.write(f"  {key} = {val}\n")
            f.write('/\n')

            # System section
            f.write('&system\n')
            for key, value in qe_input_data['system'].items():
                if isinstance(value, bool):
                    val = '.true.' if value else '.false.'
                elif isinstance(value, str):
                    val = f"'{value}'"
                else:
                    val = value
                f.write(f"  {key} = {val}\n")
            f.write('/\n')

            # Electrons section
            f.write('&electrons\n')
            for key, value in qe_input_data['electrons'].items():
                f.write(f"  {key} = {value}\n")
            f.write('/\n')

            # Atomic species
            f.write('\nATOMIC_SPECIES\n')
            unique_symbols = set(ase_structure.get_chemical_symbols())
            for symbol in unique_symbols:
                mass = atomic_masses[atomic_numbers[symbol]]
                f.write(f"  {symbol.title()} {mass:.4f} {pps[symbol]}\n")

            # Atomic positions
            f.write('\nATOMIC_POSITIONS angstrom\n')
            abs_pos = ase_structure.get_positions()
            for atom, pos in zip(ase_structure, abs_pos):
                f.write(f"  {atom.symbol.title()} {pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f}\n")

            # K-points grid
            f.write('\nK_POINTS automatic\n')
            f.write(f"  {kpoints[0]} {kpoints[1]} {kpoints[2]} 0 0 0\n")

            # Cell parameters in angstrom
            f.write('\nCELL_PARAMETERS angstrom\n')
            for vec in ase_structure.cell:
                f.write(f"  {vec[0]:.10f} {vec[1]:.10f} {vec[2]:.10f}\n")

    # Create QE input file
    write_espresso_input(structure, input_data, pseudopotentials, kpts, input_filename)

    # Run QE as subprocess
    try:
        with open(output_filename, 'w') as f_out:
            subprocess.run(
                [command, '-in', input_filename],
                stdout=f_out,
                stderr=subprocess.STDOUT,
                check=True
            )
        print("QE calculation completed successfully.")

    except CalledProcessError as cpe:
        print(f"Error running QE: {cpe}")
        try:
            with open(output_filename, 'r') as f_out:
                print("\nQE Output:")
                print(f_out.read())
        except Exception as e:
            print(f"Could not read output file: {e}")

    except Exception as e:
        print(f"Unexpected error: {e}")

config = {
    'command': '/content/bin/pw.x',
    'pseudo_dir': '/content/ONCVPSP/abinit/',
    'outdir': '/content/tmp'
}

input_data = {
    'control': {
        'calculation': 'scf',
        'restart_mode': 'from_scratch',
        'pseudo_dir': config['pseudo_dir'],
        'outdir': config['outdir'],
        'prefix': 'Bi',
        'disk_io': 'medium',
        'wf_collect': True,
        'tprnfor': True,
        'tstress': True,
    },
    'system': {
        'ibrav': 0, # fix
        'nat': 2,
        'ntyp': 1,
        'ecutwfc': 50,
        'ecutrho': 506,
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.01,
        'noncolin': True,
        'lspinorb': True,
        'nbnd': 32,
    },
    'electrons': {
        'conv_thr': 1.0e-6,
        'mixing_beta': 0.7,
        'electron_maxstep': 100,
    }
}

run_qe(bi, config, input_data, kpts=(4,4,4))
