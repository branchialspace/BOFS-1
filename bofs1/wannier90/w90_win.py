# Write Wannier90 .win input file

import os
import re
from pathlib import Path


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
        # Lattice Vectors (CELL_PARAMETERS)
        lattice = []
        cell_block = re.search(r'CELL_PARAMETERS\s+[a-zA-Z]+\s*\n(.*?)(?=\n\s*[A-Z_]+|\Z)', content, re.DOTALL | re.IGNORECASE)
        if cell_block:
            for line in cell_block.group(1).strip().split('\n'):
                parts = line.split()
                if len(parts) >= 3:
                    lattice.append([float(x) for x in parts[:3]])
        # Atomic Positions
        atoms = []
        pos_block = re.search(r'ATOMIC_POSITIONS\s+[a-zA-Z]+\s*\n(.*?)(?=\n\s*[A-Z_]+|\Z)', content, re.DOTALL | re.IGNORECASE)
        if pos_block:
            for line in pos_block.group(1).strip().split('\n'):
                parts = line.split()
                if len(parts) >= 4:
                    # Format: Symbol x y z
                    sym = parts[0]
                    coords = [float(x) for x in parts[1:4]]
                    atoms.append((sym, coords))

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

    def get_wannier_projections(atoms_list, pseudo_dict, pseudo_dir):
            """
            Identify valence orbitals based on available UPF channels.
            """
            unique_symbols = sorted(list(set(a[0] for a in atoms_list)))
            projections = []
            orb_counts = {}
            for symbol in unique_symbols:
                # Locate Pseudo
                pp_filename = pseudo_dict.get(symbol)
                if not pp_filename:
                    continue
                pp_path = Path(pseudo_dir) / pp_filename.strip("'").strip('"')
                with open(pp_path, 'r') as f:
                    content = f.read()
                # Extract all available channels (n, l)
                pattern = (r'<PP_CHI\.(\d+).*?label\s*=\s*"([^"]+)"')
                matches = re.findall(pattern, content, re.DOTALL)
                # Select all orbitals found in the pseudo
                found_orbitals = set()
                for _, label in matches:
                    l_char = next((c.lower() for c in label if c.lower() in 'spdf'), None)
                    if l_char:
                        found_orbitals.add(l_char)
                selected_types = []
                for orb in ['s', 'p', 'd', 'f']:
                    if orb in found_orbitals:
                        selected_types.append(orb)
                # format
                if selected_types:
                    order = {'s':0, 'p':1, 'd':2, 'f':3}
                    selected_types.sort(key=lambda x: order[x])
                    projections.append(f"{symbol}: {'; '.join(selected_types)}")
                # Count orbitals for num_wann
                count = 0
                for t in selected_types:
                    count += {'s':1, 'p':3, 'd':5, 'f':7}.get(t, 0)
                orb_counts[symbol] = count
            # Calculate Total num_wann
            spinors = config.get("spinors", "false").lower() == "true"
            spin_factor = 2 if spinors else 1
            total_wann = 0
            for sym, _ in atoms_list:
                total_wann += (orb_counts.get(sym, 0) * spin_factor)

        return projections, total_wann

    def write_win_file(filename, config, mp_grid, e_fermi, num_bands, num_wann, lattice, atoms, projections):
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
            # Disentanglement Windows (Offsets)
            if "dis_win_min_offset" in config:
                f.write(f"\n! Fermi Energy: {e_fermi:.4f} eV\n")
                f.write(f"dis_win_min = {e_fermi + config['dis_win_min_offset']:.4f}\n")
                f.write(f"dis_win_max = {e_fermi + config['dis_win_max_offset']:.4f}\n")
                f.write(f"dis_froz_min = {e_fermi + config['dis_froz_min_offset']:.4f}\n")
                f.write(f"dis_froz_max = {e_fermi + config['dis_froz_max_offset']:.4f}\n")
            # Geometry
            f.write("\nbegin unit_cell_cart\nAng\n")
            for vec in lattice:
                f.write(f"{vec[0]:.6f}  {vec[1]:.6f}  {vec[2]:.6f}\n")
            f.write("end unit_cell_cart\n")
            f.write("\nbegin atoms_frac\n")
            for sym, pos in atoms:
                f.write(f"{sym:<4}  {pos[0]:>10.6f}  {pos[1]:>10.6f}  {pos[2]:>10.6f}\n")
            f.write("end atoms_frac\n")
            # Projections
            if projections:
                f.write("\nbegin projections\n")
                for line in projections:
                    f.write(f"{line}\n")
                f.write("end projections\n")
            # K-path
            if config.get("bands_plot") == "true" and "kpoint_path" in config:
                f.write("\nbegin kpoint_path\n")
                for segment in config["kpoint_path"]:
                    f.write(f"{segment}\n")
                f.write("end kpoint_path\n")
                
    # Args
    output_filename = Path(pwi_path).stem + ".win"
    pseudo_dir, pseudo_dict, mp_grid, lattice, atoms = parse_pwi_data(pwi_path)
    e_fermi, num_bands = parse_pwo_data(pwo_path)
    projections, num_wann = get_wannier_projections(atoms, pseudo_dict, pseudo_dir)
    write_win_file(output_filename, config, mp_grid, e_fermi, num_bands, num_wann, lattice, atoms, projections)
    print(f"Successfully wrote {output_filename} with {num_wann} Wannier functions.")
