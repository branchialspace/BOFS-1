# OPTIMADE, MPContribs databases

!pip install optimade[http_client]
!pip install mpcontribs-client

from optimade.client import OptimadeClient

client = OptimadeClient(
    include_providers={"mp"},
)
# client.get('elements HAS "Bi"')
client.list_properties("structures")

from mpcontribs.client import Client

client = Client()
client.available_query_params()  # print list of available query parameters

query = {"formula__contains": "Bi", "project": "qmof"}
fields = ["id", "identifier", "formula", "structures"]
data = client.query_contributions(
    query=query, fields=fields, sort="", paginate=False
)

# QMOF, ARC-MOF databases

!gdown 1OcIWju1Sj0Q7ruojNlHXJRqabXJ18mIw
!gdown 1e32YGv96-ceEQAmuxy4pY4sgJN7B0tCL
!gdown 16ZjIqFg-MrXqB4O74H3dTe0toSiZTbfg
!unzip QMOF.zip
!unzip qmof_database.zip
!unzip qmof_database/relaxed_structures.zip
!tar -xzf ARCMOF_1.tar.gz
!tar -xzf ARCMOF_2.tar.gz

import os
from pathlib import Path
from tqdm import tqdm
import shutil


def find_species_in_cifs(directory_paths, target_species, output_dir):
    count = 0
    
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
    
    # Process all CIF files with a progress bar
    print(f"\nAnalyzing {len(all_cif_files)} CIF files for {target_species}...")
    for cif_path in tqdm(all_cif_files, desc="Processing CIFs"):
        try:
            with open(cif_path, 'r', encoding='utf-8') as f:
                found_species = False
                for line in f:
                    # Skip empty lines and comment lines
                    if not line.strip() or line.strip().startswith('#'):
                        continue
                    # Split the line and check first column for species
                    parts = line.strip().split()
                    if parts and target_species in parts[0]:
                        found_species = True
                        break
                
                if found_species:
                    # Copy file to output directory
                    shutil.copy2(cif_path, output_path / cif_path.name)
                    count += 1
                    
        except Exception as e:
            print(f"\nError processing {cif_path}: {str(e)}")
    
    print(f"\nFound and copied {count} files containing {target_species} to {output_dir}")
    return count

directories = ['/content/all_structures_1', '/content/all_structures_2', '/content/relaxed_structures']
output_directory = '/content/bimofs'
total_files = find_species_in_cifs(directories, 'Bi', output_directory)
