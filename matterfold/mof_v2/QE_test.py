from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.optimize import LBFGS
import os
from pathlib import Path

# Create the rocksalt structure
rocksalt = bulk('NaCl', crystalstructure='rocksalt', a=6.0)

# Get unique species from structure
species = set(atom.symbol for atom in rocksalt)

# Find pseudopotentials for each species
pseudo_dir = '/content/pseudo_sssp'
pseudopotentials = {}
for symbol in species:
    pseudo_path = Path(pseudo_dir)
    matching_files = []
    for extension in ['.upf', '.UPF']:
        for file in pseudo_path.glob(f"*{extension}"):
            # Extract the species name from filename (before first . or _)
            filename = file.name
            species_name = filename.split('.')[0].split('_')[0]
            if species_name.lower() == symbol.lower():
                matching_files.append(file)    
    if matching_files:
        # If multiple files found, use the first one
        pseudopotentials[symbol] = matching_files[0].name
    else:
        raise FileNotFoundError(f"No pseudopotential file found for {symbol}")

# Define input parameters explicitly
input_data = {
    'control': {
        'calculation': 'relax',  # Changed to relax to enable force calculations
        'restart_mode': 'from_scratch',
        'pseudo_dir': '/content/pseudo_sssp',
        'outdir': './tmp',
        'prefix': 'NaCl',
        'disk_io': 'low',
        'wf_collect': True,
        'tprnfor': True,     # Enable force calculation
        'tstress': True,     # Enable stress calculation
    },
    'system': {
        'ecutwfc': 50,
        'ecutrho': 400,
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.01,
    },
    'electrons': {
        'conv_thr': 1.0e-6,
        'mixing_beta': 0.7,
        'electron_maxstep': 100,
    },
    'ions': {
        'ion_dynamics': 'bfgs',
        'upscale': 100,      # Scale factor for BFGS
    },
}
# Create profile with explicit paths
profile = EspressoProfile(
    command='/content/bin/pw.x',
    pseudo_dir='/content/pseudo_sssp'
)

os.makedirs('./tmp', exist_ok=True)
# Initialize calculator with more detailed parameters
calc = Espresso(
    profile=profile,
    pseudopotentials=pseudopotentials,
    input_data=input_data,
    kpts=(4, 4, 4),
)
rocksalt.calc = calc
try:
    # Run single point calculation first
    energy = rocksalt.get_potential_energy()
    print(f"Single point energy: {energy} eV")

    # Get forces to verify they're being calculated
    forces = rocksalt.get_forces()
    print("\nInitial forces:")
    print(forces)

    # Setup and run geometry optimization
    opt = LBFGS(rocksalt, trajectory='NaCl_opt.traj')
    opt.run(fmax=0.005)

    # Print optimized lattice constant
    final_lattice = (8 * rocksalt.get_volume() / len(rocksalt)) ** (1.0 / 3.0)
    print(f"\nOptimized lattice constant: {final_lattice:.3f} Å")
except Exception as e:
    print(f"Error occurred: {str(e)}")
    # Print additional debugging information
    print("\nCalculator parameters:")
    print(calc.parameters)
