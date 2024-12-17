# Ligtand Metal Docking / Generate Secondary Building Unit

import subprocess


def ligand_metal_docking(host, guest_xyz_file, charge=0, mult=1, docking_method='!XTB', docker_options=None):
    """
    Performs an ORCA docking calculation using the given HOST (as an ASE Atoms object)
    and a specified GUEST .xyz file.

    Parameters:
    - host: ASE Atoms object representing the host molecule.
    - guest_xyz_file: Path to the guest .xyz file (string).
    - charge: Overall charge of the host complex (usually 0).
    - mult: Spin multiplicity of the host complex (usually 1).
    - docking_method: PES method to use for docking (e.g., 'XTB').
    - docker_options: Dictionary of additional options to customize the DOCKER input.

    Returns:
    - None. The function runs ORCA and prints the docking results.
      The best docked structure and associated files are written to the current directory.
    """

    # Define the ORCA executable path (adjust as needed)
    orca_path = '/root/orca_6_0_1/orca'

    # Default DOCKER options if none provided
    if docker_options is None:
        docker_options = {}

    # Write the ORCA DOCKER input file
    input_filename = 'orca_docker_input.inp'
    output_filename = 'orca_docker_output.out'
    write_orca_docker_input(host, guest_xyz_file, input_filename, charge, mult, docking_method, docker_options)

    # Run ORCA docking calculation
    try:
        with open(output_filename, 'w') as f_out:
            subprocess.run([orca_path, input_filename], stdout=f_out, stderr=subprocess.STDOUT, check=True)
        print("ORCA DOCKER calculation completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during ORCA DOCKER calculation: {e}")

def write_orca_docker_input(host, guest_xyz_file, filename, charge, mult, docking_method, docker_options):

    # Prepare DOCKER block
    # The minimal DOCKER block requires a GUEST line:
    # %DOCKER GUEST "guest.xyz" END
    docker_lines = ["%DOCKER"]
    docker_lines.append(f'    GUEST "{guest_xyz_file}"')
    for key, value in docker_options.items():
        docker_lines.append(f'    {key} {value}')
    docker_lines.append("END")

    with open(filename, 'w') as f:
        # Write the docking method
        f.write(docking_method + "\n")

        # Write the DOCKER block
        f.write("\n".join(docker_lines))
        f.write("\n")

        # Write host coordinates
        f.write(f"* xyz {charge} {mult}\n")
        for atom in host:
            symbol = atom.symbol
            x, y, z = atom.position
            f.write(f"  {symbol} {x:.6f} {y:.6f} {z:.6f}\n")
        f.write("*\n")

if __name__ == "__main__":
  
    perform_docking(host=ligand, guest_xyz_file="/content/Bi_center.xyz", charge=-4,
                    docker_options={'GUESTCHARGE': 2, 'FIXHOST': 'FALSE', 'NREPEATGUEST': 2, 'DOCKLEVEL': 'COMPLETE'})
