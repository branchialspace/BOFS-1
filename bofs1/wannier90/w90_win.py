# Write Wannier90 .win input file

import re
import sys
from pathlib import Path
import numpy as np
from ase.io import read
import xml.etree.ElementTree as ET
import seekpath
from seekpath.util import atoms_num_dict


def w90_win(
    pwo_path,
    pwi_path,
    config
):
    """
    Generate a Wannier90 .win input file by parsing QE input/output, pseudopotential files, using SeeK-Path for k-point path.
    pwo_path : str
        Path to the Quantum ESPRESSO output file (for geometry, Fermi energy).
    pwi_path : str
        Path to the Quantum ESPRESSO input file (for pseudos, k-points, atom positions, lattice vectors).
    config : dict
        Configuration dictionary containing W90 job control parameters.
    """
  
    def parse_pwi_data(pwi_path):
        """
        Extract pseudo_dir, pseudo filenames, k points, lattice vectors, and atomic positions from .pwi.
        """
        with open(pwi_path, 'r') as f:
            content = f.read()
        # Pseudo Directory
        pseudo_dir_match = re.search(r'pseudo_dir\s*=\s*[\'"]([^\'"]+)[\'"]', content, re.IGNORECASE)
        pseudo_dir = pseudo_dir_match.group(1) if pseudo_dir_match else "./"
        # Atomic Species Mapping
        pseudo_dict = {}
        species_block = re.search(r'ATOMIC_SPECIES\s*\n(.*?)(?:\n\s*\n)', content, re.DOTALL | re.IGNORECASE)
        if species_block:
            for line in species_block.group(1).strip().split('\n'):
                parts = line.split()
                if len(parts) >= 3:
                    # Format: Symbol Mass Filename
                    pseudo_dict[parts[0]] = parts[2]
        # K-Points (MP Grid)
        mp_grid = [1, 1, 1]
        k_block = re.search(r'K_POINTS\s+([a-zA-Z]+)\s*\n\s*([\d\s]+)', content, re.IGNORECASE)
        if k_block and 'automatic' in k_block.group(1).lower():
            k_data = k_block.group(2).split()
            if len(k_data) >= 3:
                mp_grid = [int(k_data[0]), int(k_data[1]), int(k_data[2])]
        # Lattice Vectors
        structure = read(pwi_path, format='espresso-in')
        lattice = np.array(structure.get_cell())
        # Atomic Positions
        scaled_positions = structure.get_scaled_positions()
        symbols = structure.get_chemical_symbols()
        atoms = []
        for sym, pos in zip(symbols, scaled_positions):
            atoms.append((sym, pos))

        return pseudo_dir, pseudo_dict, mp_grid, lattice, atoms

    def parse_pwo_data(pwo_path):
        """
        Extract Fermi energy and number of bands from .pwo.
        """
        with open(pwo_path, 'r') as f:
            lines = f.readlines()
            content = "".join(lines)
        # Fermi Energy & Bands
        fermi_match = re.search(r'the Fermi energy is\s+([-\d.]+)\s+ev', content, re.IGNORECASE)
        e_fermi = float(fermi_match.group(1)) if fermi_match else 0.0
        nbnd_match = re.search(r'number of Kohn-Sham states=\s+(\d+)', content)
        num_bands = int(nbnd_match.group(1)) if nbnd_match else 0
        
        return e_fermi, num_bands
    
    def total_wannier(atoms_list, pseudo_dict, pseudo_dir, config):
        """
        Calculates total_wann for SCDM by parsing UPF files.
        Prioritizes PP_RELWFC (jchi) for fully relativistic PPs.
        Falls back to PP_PSWFC (l) for scalar PPs.
        Handles SOC doubling when necessary.
        """
        spinors = str(config.get("spinors", "false")).lower() == "true"
        unique_symbols = sorted(list(set(a[0] for a in atoms_list)))
        orb_counts = {}
        for symbol in unique_symbols:
            pp_filename = pseudo_dict.get(symbol)
            if not pp_filename:
                continue
            pp_path = Path(pseudo_dir) / pp_filename.strip("'").strip('"')
            if not pp_path.exists():
                print(f"Warning: UPF {pp_path} not found. Skipping.")
                orb_counts[symbol] = 0
                continue
            try:
                tree = ET.parse(pp_path)
                root = tree.getroot()
                atom_wfc_count = 0
                is_fully_relativistic = False
                # Check for Fully Relativistic Block <PP_SPIN_ORB>
                spin_orb_section = root.find('PP_SPIN_ORB')
                if spin_orb_section is not None:
                    is_fully_relativistic = True
                    # Find PP_RELWFC (Wavefunctions)
                    for child in spin_orb_section:
                        if "PP_RELWFC" in child.tag:
                            try:
                                # Extract jchi
                                jchi = float(child.get('jchi'))
                                # Degeneracy = 2j + 1
                                deg = int(2 * jchi + 1)
                                atom_wfc_count += deg
                            except (ValueError, TypeError):
                                continue                                
                # Fallback if Scalar <PP_PSWFC>
                else:
                    pswfc_section = root.find('PP_PSWFC')
                    if pswfc_section is None:
                        pswfc_section = root.find('PP_WAVEFUNCTIONS') # Support for older UPF formats
                    if pswfc_section is not None:
                        for wfc in pswfc_section:
                            l = int(wfc.get('l'))
                            # Degeneracy = 2l + 1 (Spatial only)
                            atom_wfc_count += (2 * l + 1)
                # If the UPF is scalar but we are running an SOC calculation, multiply by 2 to account for spin up/down.
                # If the UPF is fr, spin is included in jchi.
                if spinors and not is_fully_relativistic:
                    atom_wfc_count *= 2
                orb_counts[symbol] = atom_wfc_count                
            except Exception as e:
                print(f"Error parsing UPF for {symbol}: {e}")
                orb_counts[symbol] = 0
        total_wann = 0
        for sym, _ in atoms_list:
            total_wann += orb_counts.get(sym, 0)
    
        return total_wann

    def get_kpoint_path(lattice, atoms):
        """
        Use SeeK-path to determine the k-point path for band structure plotting.
        lattice : list
            3x3 list of lattice vectors in Angstrom.
        atoms : list
            List of tuples (symbol, [x, y, z]) with fractional coordinates.    
        Returns
        kpoint_path : list
            List of strings formatted for Wannier90 kpoint_path block.
        path_info : dict
            Dictionary containing additional path information (labels, coordinates, etc.)
        """
        # Convert atoms to seekpath format
        positions = [atom[1] for atom in atoms]
        numbers = [atoms_num_dict.get(atom[0], 0) for atom in atoms]
        structure = (lattice, positions, numbers)
        # Get k-path for the original cell (no standardization)
        path_result = seekpath.get_path_orig_cell(
            structure,
            with_time_reversal=True,
            recipe='hpkot',
            threshold=1.0e-7,
            symprec=1e-5,
            angle_tolerance=-1.0)
        # Extract path and point coordinates
        path = path_result['path']
        point_coords = path_result['point_coords']
        # Wannier90 format: "LABEL1 k1x k1y k1z  LABEL2 k2x k2y k2z"
        kpoint_path = []
        for start_label, end_label in path:
            start_coords = point_coords[start_label]
            end_coords = point_coords[end_label]
            # Format label (replace GAMMA with G for brevity)
            start_label_fmt = "G" if start_label == "GAMMA" else start_label
            end_label_fmt = "G" if end_label == "GAMMA" else end_label
            line = (
                f"{start_label_fmt} {start_coords[0]:.6f} {start_coords[1]:.6f} {start_coords[2]:.6f}  "
                f"{end_label_fmt} {end_coords[0]:.6f} {end_coords[1]:.6f} {end_coords[2]:.6f}")
            kpoint_path.append(line)
        # Collect path info for reference
        path_info = {
            'bravais_lattice': path_result.get('bravais_lattice', 'unknown'),
            'spacegroup_number': path_result.get('spacegroup_number', None),
            'spacegroup_international': path_result.get('spacegroup_international', None),
            'point_coords': point_coords,
            'path': path,
            'has_inversion_symmetry': path_result.get('has_inversion_symmetry', None),
            'is_supercell': path_result.get('is_supercell', False)}
        
        return kpoint_path, path_info

    def w90_kmeshpl(nk1, nk2, nk3):
        """
        kmesh.pl k-point grid method for consistency between nscf and wannier90.
        Returns
        content_nscf: Block for QE nscf input (includes weights)
        content_win:  Block for wannier90.win (coordinates only)
        """
        kpoints = []
        # Wannier90 kmesh.pl logic:
        # for ($x=0; $x<$ARGV[0]; $x++) {      <-- Outer
        #   for ($y=0; $y<$ARGV[1]; $y++) {    <-- Middle
        #     for ($z=0; $z<$ARGV[2]; $z++) {  <-- Inner (Fastest)
        for i in range(nk1):          # x direction
            for j in range(nk2):      # y direction
                for k in range(nk3):  # z direction
                    x_frac = i / nk1
                    y_frac = j / nk2
                    z_frac = k / nk3
                    kpoints.append((x_frac, y_frac, z_frac))
        # W90 kmesh.pl weights: 1 / totpts
        total_pts = nk1 * nk2 * nk3
        weight = 1.0 / total_pts
        # NSCF k-point grid
        nscf_lines = []
        nscf_lines.append("K_POINTS crystal")
        nscf_lines.append(f"{len(kpoints)}")
        for kp in kpoints:
            # W90 kmesh.pl format: %12.8f%12.8f%12.8f%14.6e
            line = f"  {kp[0]:12.8f}{kp[1]:12.8f}{kp[2]:12.8f}  {weight:14.6e}"
            nscf_lines.append(line)
        content_nscf = "\n".join(nscf_lines)
        # Wannier90 k-point grid
        win_lines = []
        win_lines.append("begin kpoints")
        for kp in kpoints:
            # W90 kmesh.pl format: %12.8f%12.8f%12.8f
            line = f"  {kp[0]:12.8f}{kp[1]:12.8f}{kp[2]:12.8f}"
            win_lines.append(line)
        win_lines.append("end kpoints")
        content_win = "\n".join(win_lines)
    
        return content_nscf, content_win

    def write_win_file(filename, config, mp_grid, e_fermi, num_bands, num_wann, lattice, atoms, kpoint_path, kpoint_block, path_info):
        """
        Write the final .win file.
        """
        with open(filename, 'w') as f:
            f.write("! Generated by w90_win toolkit\n\n")
            # System
            f.write(f"num_bands = {num_bands}\n")
            f.write(f"num_wann = {num_wann}\n")
            f.write(f"mp_grid = {mp_grid[0]} {mp_grid[1]} {mp_grid[2]}\n")
            # Config params
            for key, val in config.items():
                f.write(f"{key} = {val}\n")
            # Geometry
            f.write("\nbegin unit_cell_cart\nAng\n")
            for vec in lattice:
                f.write(f"{vec[0]:.6f}  {vec[1]:.6f}  {vec[2]:.6f}\n")
            f.write("end unit_cell_cart\n")
            f.write("\nbegin atoms_frac\n")
            for sym, pos in atoms:
                f.write(f"{sym:<4}  {pos[0]:>10.6f}  {pos[1]:>10.6f}  {pos[2]:>10.6f}\n")
            f.write("end atoms_frac\n")
            # K-path (from SeeK-path)
            if config.get("bands_plot") == "true" and kpoint_path:
                f.write("\n! K-path determined by SeeK-path (HPKOT convention)\n")
                if path_info:
                    f.write(f"! Bravais lattice: {path_info.get('bravais_lattice', 'unknown')}\n")
                    f.write(f"! Space group: {path_info.get('spacegroup_international', 'unknown')} "
                            f"(#{path_info.get('spacegroup_number', '?')})\n")
                    if path_info.get('is_supercell'):
                        f.write("! WARNING: Cell detected as supercell - k-path is for primitive cell\n")
                f.write("begin kpoint_path\n")
                for segment in kpoint_path:
                    f.write(f"{segment}\n")
                f.write("end kpoint_path\n")
            # K-points
            f.write(f"\n{kpoint_block}\n")
                
    # Args
    output_filename = Path(pwi_path).stem + ".win"
    pseudo_dir, pseudo_dict, mp_grid, lattice, atoms = parse_pwi_data(pwi_path)
    e_fermi, num_bands = parse_pwo_data(pwo_path)
    num_wann = total_wannier(atoms, pseudo_dict, pseudo_dir, config)
    kpoint_path, path_info = get_kpoint_path(lattice, atoms)
    _, kpoint_block = w90_kmeshpl(mp_grid[0], mp_grid[1], mp_grid[2])
    # Write .win
    write_win_file(output_filename, config, mp_grid, e_fermi, num_bands, num_wann, lattice, atoms, kpoint_path, kpoint_block, path_info)
    print(f"Successfully wrote {output_filename} with {num_wann} Wannier functions.")

if __name__ == "__main__":
    pwo_arg = sys.argv[1]
    pwi_arg = sys.argv[2]
    config_arg = sys.argv[3]
    config_scope = {}
    with open(config_arg, 'r') as f:
        exec(f.read(), {}, config_scope)
    config_dict = list(config_scope.values())[0]

    w90_win(pwo_arg, pwi_arg, config_dict)
