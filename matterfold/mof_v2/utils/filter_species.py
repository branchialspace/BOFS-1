# Filter dataset directories for specified atomic species

!gdown  # QMOF
!gdown  # ARCMOF 1
!gdown  # ARCMOF 2
!gdown  # CSD MOF
!gdown  # CO2 MOF
!unzip QMOF.zip
!unzip qmof_database.zip
!unzip qmof_database/relaxed_structures.zip
!unzip CSD_MOF_Collection.zip
!tar -xzf ARCMOF_1.tar.gz
!tar -xzf ARCMOF_2.tar.gz
!mkdir /content/co2
!tar -xzf co2_MOF_database.tar.gz -C /content/co2

import os
from pathlib import Path
from tqdm import tqdm
import shutil

def find_species_in_cifs(directory_paths, target_species, output_dir):
    count = 0
    directory_counts = {}

    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # First, collect all .cif files from all directories
    all_cif_files = []
    print("Collecting CIF files from directories...")
    for directory in tqdm(directory_paths, desc="Scanning directories"):
        for root, _, files in os.walk(directory):
            cif_files = [Path(root) / f for f in files if f.endswith('.cif')]
            all_cif_files.extend(cif_files)

    # Process all CIF files
    print(f"\nAnalyzing {len(all_cif_files)} CIF files for {target_species}...")
    for cif_path in tqdm(all_cif_files, desc="Processing CIFs"):
        try:
            with open(cif_path, 'r', encoding='utf-8') as f:
                found_species = False
                for line in f:
                    # Skip empty lines and comment lines
                    if not line.strip() or line.strip().startswith('#'):
                        continue
                    # Split the line and check first or second column for species
                    parts = line.strip().split()
                    if parts and (target_species in parts[0] or (len(parts) > 1 and target_species in parts[1])):
                        found_species = True
                        break

                if found_species:
                    # Copy file to output directory
                    shutil.copy2(cif_path, output_path / cif_path.name)
                    count += 1
                    directory_counts[str(cif_path.parent)] = directory_counts.get(str(cif_path.parent), 0) + 1

        except Exception as e:
            print(f"\nError processing {cif_path}: {str(e)}")

    print(f"\nFound and copied {count} mofs containing {target_species} to {output_dir}")
    for directory, dir_count in directory_counts.items():
        print(f"  {directory}: Found {dir_count} mofs containing {target_species}")
    return count

directories = ['/content/all_structures_1', '/content/all_structures_2', '/content/relaxed_structures', '/content/CSD_MOF_Collection', '/content/co2']
output_directory = '/content/bimofs'
total_files = find_species_in_cifs(directories, 'Bi', output_directory)
