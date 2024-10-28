# NBO bonding site analysis in ORCA

import subprocess


def electron_bonding_sites(ligand, charge=0, mult=1, method='B3LYP', basis_set='def2-SVP'):
    # Define the ORCA executable path
    orca_path = '/root/orca_6_0_0/orca'

    # Define filenames
    input_filename = 'orca_input.inp'
    output_filename = 'orca.out'

    # Write the ORCA input file
    write_orca_input(ligand, input_filename, charge, mult, method, basis_set)

    # Run ORCA
    try:
        with open(output_filename, 'w') as f_out:
            subprocess.run([orca_path, input_filename], stdout=f_out, stderr=subprocess.STDOUT, check=True)
        print("ORCA calculation completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during ORCA calculation: {e}")
        return {}

def write_orca_input(ligand, filename, charge, mult, method, basis_set):
    with open(filename, 'w') as f:
        # Write method and basis set
        f.write(f"! {method} {basis_set} NBO\n")
        f.write("%output\n")
        f.write("    Print[P_NBO] 1\n")
        f.write("    Print[P_NPA] 1\n")
        f.write("end\n\n")
        # Write charge and multiplicity
        f.write(f"* xyz {charge} {mult}\n")
        # Write atomic coordinates
        for atom in ligand:
            symbol = atom.symbol
            x, y, z = atom.position
            f.write(f"  {symbol} {x:.6f} {y:.6f} {z:.6f}\n")
        f.write("*\n")
