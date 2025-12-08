import os
import re
from pathlib import Path
from mendeleev import element

def w90_win(
    pwo_path,
    pwi_path,
    config
):
    """
    Generate a Wannier90 .win input file by parsing QE input/output files.
    pwo_path : str
        Path to the Quantum ESPRESSO output file (for geometry, Fermi energy).
    pwi_path : str
        Path to the Quantum ESPRESSO input file (for pseudos, k-points).
    config : dict
        Configuration dictionary containing W90 job control parameters.
    """
  
    def parse_pwi_data(pwi_path):
        """
        Extract pseudo_dir, pseudo filenames, and k points for mp_grid from .pwi.
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

        return pseudo_dir, pseudo_dict, mp_grid

    def parse_pwo_data(pwo_path):
        """
        Extract Fermi energy, number of bands, lattice vectors, and atomic positions from .pwo.
        """
        with open(pwo_path, 'r') as f:
            lines = f.readlines()
            content = "".join(lines)
        # Fermi Energy & Bands
        fermi_match = re.search(r'the Fermi energy is\s+([-\d.]+)\s+ev', content, re.IGNORECASE)
        e_fermi = float(fermi_match.group(1)) if fermi_match else 0.0
        nbnd_match = re.search(r'number of Kohn-Sham states=\s+(\d+)', content)
        num_bands = int(nbnd_match.group(1)) if nbnd_match else 0
        # Lattice Parameters (alat conversion)
        alat_match = re.search(r'lattice parameter \(alat\)\s+=\s+([-\d.]+)\s+a.u.', content)
        alat = float(alat_match.group(1)) if alat_match else 1.0
        bohr_to_ang = 0.52917721067
        alat_ang = alat * bohr_to_ang
        # Lattice Vectors
        lattice = []
        for i, line in enumerate(lines):
            if "crystal axes:" in line:
                for j in range(1, 4):
                    vec_nums = re.findall(r'([-\d.]+)', lines[i+j])
                    # Format: a(1) = ( x y z ) -> take last 3
                    if len(vec_nums) >= 3:
                        lattice.append([float(x) * alat_ang for x in vec_nums[-3:]])
                break
        # 4. Atoms (Crystallographic/Fractional)
        atoms = [] # List of (Symbol, [x, y, z])
        atom_start = -1
        for i, line in enumerate(lines):
            if "positions (cryst. coord.)" in line:
                atom_start = i + 1
                break
        if atom_start != -1:
            for line in lines[atom_start:]:
                if not line.strip() or "k points" in line or "End" in line:
                    break
                parts = line.split()
                # Format: index Symbol tau(i) = ( x y z )
                if len(parts) > 5 and "tau" in line:
                    sym = parts[1]
                    coords_raw = re.findall(r'\((.*?)\)', line)[-1]
                    coords = [float(x) for x in coords_raw.split()]
                    atoms.append((sym, coords))

        return e_fermi, num_bands, lattice, atoms

    def get_wannier_projections(atoms_list, pseudo_dict, pseudo_dir):
        """
        Identify valence orbitals based on chemical block and available UPF channels.
        - s-block: s
        - p-block: s, p
        - d-block: s, d (and p if highly excited)
        - f-block: s, d, f
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
            if not pp_path.exists():
                print(f"Warning: UPF {pp_path} not found. Skipping projections for {symbol}.")
                continue
            with open(pp_path, 'r') as f:
                content = f.read()
            # Extract all available channels (n, l)
            pattern = (r'<PP_CHI\.(\d+).*?label\s*=\s*"([^"]+)"')
            matches = re.findall(pattern, content, re.DOTALL)
            # Map l-character to list of n values found in pseudo
            # e.g. {'s': [5, 6], 'p': [6], 'd': [5]}
            available_channels = {'s': [], 'p': [], 'd': [], 'f': []}
            for _, label in matches:
                n_digits = ''.join(filter(str.isdigit, label))
                n = int(n_digits) if n_digits else 0
                l_char = next((c.lower() for c in label if c.lower() in 'spdf'), None)
                if l_char:
                    available_channels[l_char].append(n)
            # Select Orbitals based on Blocks
            block = element(symbol).block
            selected_types = []
            # S-Block
            if block == 's':
                if available_channels['s']: selected_types.append('s')
            # P-Block
            elif block == 'p':
                # Always take s and p if available
                if available_channels['s']: selected_types.append('s')
                if available_channels['p']: selected_types.append('p')
                # Explicitly ignore d/f here (semicore for p-block)
            # D-Block
            elif block == 'd':
                if available_channels['s']: selected_types.append('s')
                if available_channels['d']: selected_types.append('d')
                # Optional: Include p if it's the outermost shell (n_p == n_s)
                # This helps with localization for some transition metals
                if available_channels['p'] and available_channels['s']:
                    if max(available_channels['p']) == max(available_channels['s']):
                        selected_types.append('p')
            # F-Block
            elif block == 'f':
                if available_channels['s']: selected_types.append('s')
                if available_channels['d']: selected_types.append('d')
                if available_channels['f']: selected_types.append('f')
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

    def write_win_file(filename, config, mp_grid, e_fermi, num_bands, 
                       num_wann, lattice, atoms, projections):
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
    pseudo_dir, pseudo_dict, mp_grid = parse_pwi_data(pwi_path)
    e_fermi, num_bands, lattice, atoms = parse_pwo_data(pwo_path)
    projections, num_wann = get_wannier_projections(atoms, pseudo_dict, pseudo_dir)
    write_win_file(output_filename, config, mp_grid, e_fermi, num_bands, num_wann, lattice, atoms, projections)
    print(f"Successfully wrote {output_filename} with {num_wann} Wannier functions.")
