# QuantumESPRESSO Single-Point SCF with SOC

from ase.build import bulk
from ase.data import atomic_masses, atomic_numbers
import os
from pathlib import Path
import subprocess
from subprocess import CalledProcessError


# Generate the structure using ASE
bi = bulk('Bi', 'rhombohedral', a=4.75, c=12.36, orthorhombic=False)

# Configuration parameters
command = '/content/bin/pw.x'
pseudo_dir = '/content/ONCVPSP/abinit/'
outdir = '/content/tmp'
os.makedirs(outdir, exist_ok=True)

def pseudopotentials(atomic_structure, pseudo_dir):
    species = set(atom.symbol for atom in atomic_structure)
    pseudopotentials = {}
    pseudo_path = Path(pseudo_dir)
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
        pseudopotentials[symbol] = selected.name
    return pseudopotentials

pseudopotentials = pseudopotentials(bi, pseudo_dir)

# Input parameters
input_data = {
    'control': {
        'calculation': 'scf',
        'restart_mode': 'from_scratch',
        'pseudo_dir': pseudo_dir,
        'outdir': outdir,
        'prefix': 'Bi',
        'disk_io': 'medium',
        'wf_collect': True,
        'tprnfor': True,
        'tstress': True,
    },
    'system': {
        'ibrav': 0,
        'nat': 2,
        'ntyp': 1,
        'ecutwfc': 60,
        'ecutrho': 240,
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

def write_espresso_input(structure, input_data, pseudopotentials, kpts, filename='espresso.pwi'):
    with open(filename, 'w') as f:
        # Control section
        f.write('&control\n')
        for key, value in input_data['control'].items():
            val = '.true.' if isinstance(value, bool) and value else \
                  '.false.' if isinstance(value, bool) else \
                  f"'{value}'" if isinstance(value, str) else value
            f.write(f"  {key} = {val}\n")
        f.write('/\n')

        # System section
        system_params = input_data['system'].copy()
        f.write('&system\n')
        for key, value in system_params.items():
            val = '.true.' if isinstance(value, bool) and value else \
                  '.false.' if isinstance(value, bool) else \
                  f"'{value}'" if isinstance(value, str) else value
            f.write(f"  {key} = {val}\n")
        f.write('/\n')

        # Electrons section
        f.write('&electrons\n')
        for key, value in input_data['electrons'].items():
            f.write(f"  {key} = {value}\n")
        f.write('/\n')

        # Atomic species
        f.write('\nATOMIC_SPECIES\n')
        symbols = set(structure.get_chemical_symbols())
        for symbol in symbols:
            mass = atomic_masses[atomic_numbers[symbol]]
            f.write(f"  {symbol.title()} {mass:.4f} {pseudopotentials[symbol]}\n")

        # Atomic positions
        f.write('\nATOMIC_POSITIONS angstrom\n')
        abs_pos = structure.get_positions()
        for atom, pos in zip(structure, abs_pos):
            f.write(f"  {atom.symbol.title()} {pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f}\n")

        # K-points grid
        f.write('\nK_POINTS automatic\n')
        f.write(f"  {kpts[0]} {kpts[1]} {kpts[2]} 0 0 0\n")

        # Cell parameters in angstroms
        f.write('\nCELL_PARAMETERS angstrom\n')
        for vec in structure.cell:
            f.write(f"  {vec[0]:.10f} {vec[1]:.10f} {vec[2]:.10f}\n")

# Generate QE input file
write_espresso_input(bi, input_data, pseudopotentials, kpts=(4,4,4), filename='espresso.pwi')

# Run QuantumESPRESSO calculation
input_file = 'espresso.pwi'
output_file = 'espresso.pwo'

try:
    with open(output_file, 'w') as f_out:
        subprocess.run([command, '-in', input_file], stdout=f_out, stderr=subprocess.STDOUT, check=True)
    print("QE calculation completed successfully.")

except CalledProcessError as e:
    print(f"Error running QE: {e}")
    try:
        with open(output_file, 'r') as f:
            print("\nQE Output:")
            print(f.read())
    except Exception as e:
        print(f"Could not read output file: {e}")

except Exception as e:
    print(f"Unexpected error: {e}")
