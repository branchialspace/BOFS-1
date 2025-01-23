# Single-Point SCF with SOC + magnetization in QuantumESPRESSO

from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile
import os
from pathlib import Path
import subprocess
from subprocess import PIPE, CalledProcessError


bi = bulk('Bi', 'rhombohedral', a=4.75, c=12.36, orthorhombic=False)

command = '/content/bin/pw.x'
pseudo_dir = '/content/rel_pbe'
outdir = '/content/tmp'; os.makedirs(outdir, exist_ok=True)

def pseudopotentials(atomic_structure, pseudo_dir):
    species = set(atom.symbol for atom in atomic_structure)
    pseudopotentials = {}
    pseudo_path = Path(pseudo_dir)    
    for symbol in species:
        matching_files = []
        for extension in ['.upf', '.UPF']:
            for file in pseudo_path.glob(f"*{extension}"):
                # Extract species name from filename (before first '.' or '_')
                filename = file.name
                species_name = filename.split('.')[0].split('_')[0]
                if species_name.lower() == symbol.lower():
                    matching_files.append(file)
        if matching_files:
            # Use the PAW pseudopotential with longest filename if available
            paw_matches = [f for f in matching_files if "paw" in f.name.lower()]
            if paw_matches:
                # Sort by filename length in descending order and take the first one
                longest_paw = max(paw_matches, key=lambda x: len(x.name))
                pseudopotentials[symbol] = longest_paw.name
            else:
                pseudopotentials[symbol] = matching_files[0].name
        else:
            raise FileNotFoundError(f"No pseudopotential file found for {symbol}")
    
    return pseudopotentials

pseudopotentials = pseudopotentials(bi, pseudo_dir)

# Input parameters: SCF calculation, Noncollinear magnetization , SOC, Orbital magnetization
input_data = {
    'control': {
        'calculation': 'scf',
        'restart_mode': 'from_scratch',
        'pseudo_dir': pseudo_dir,
        'outdir': outdir,
        'prefix': 'Bi',
        'disk_io': 'low',
        'wf_collect': True,
        'tprnfor': True,   # Print forces
        'tstress': True,   # Print stress
    },
    'system': {
        # Basic plane-wave cutoffs
        'ecutwfc': 50,
        'ecutrho': 400,

        # Smearing parameters (adjust as appropriate)
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.01,

        'noncolin': True, # spin-orbit coupling
        'lspinorb': True, # noncollinear magnetization
        # 'lorbm': True, # orbital magnetization is needed for magneto-optical properties
    },
    'electrons': {
        'conv_thr': 1.0e-6,
        'mixing_beta': 0.7,
        'electron_maxstep': 100,
    }
}

# Calculator profile
profile = EspressoProfile(
    command=command,
    pseudo_dir=pseudo_dir
)

# Attach calculator to our bismuth structure, specify k-point grid and pseudopotentials
calc = Espresso(
    profile=profile,
    pseudopotentials=pseudopotentials,
    input_data=input_data,
    kpts=(4, 4, 4),  # Adjust k-point mesh if needed
)

bi.calc = calc

# compute SCF
try:
    energy = bi.get_potential_energy()
    print(f"Single-point (SCF) energy: {energy:.6f} eV")

    # Forces (eV/Å)
    forces = bi.get_forces()
    print("\nForces (eV/Å):")
    print(forces)

    # Stress (eV/Å^3)
    stress = bi.get_stress()
    print("\nStress (eV/Å³):")
    print(stress)

except CalledProcessError as e:
    print(f"Error occurred: {str(e)}")
    
    # Try to read the QE output file directly
    try:
        with open('espresso.pwo', 'r') as f:
            qe_output = f.read()
        print("\nQE Output:")
        print(qe_output)
    except FileNotFoundError:
        print("Could not find QE output file")
