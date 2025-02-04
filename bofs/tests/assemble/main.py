import bofs1


if __name__ == "__main__":

    metal_center = bofs1.generate_metal_polyhedron(species= 'Bi', num_atoms = 6)
    ligand, mol = bofs1.generate_ligand('C(=C/C(=O)[O-])\C(=O)[O-]') #       [Te]CC#CC#CC[Te]         C1=C(S(=O)(=O)O)C(C2=CC(=CC(=C2)C(=S)S)C(=S)S)=CC(S(=O)(=O)O)=C1C3=CC(=CC(=C3)C(=S)S)C(=S)S
    ligand_electron_analysis = bofs1.ligand_electron_analysis(ligand)
    bonding_sites = bofs1.ligand_bonding_sites(ligand=ligand, ligand_electron_analysis=ligand_electron_analysis)
    combined_structure = bofs1.ligand_metal_docking(ligand=ligand, metal_center=metal_center, bonding_sites=bonding_sites)
    periodic_structure = bofs1.mof_cell(combined_structure, metal_center, ligand, bonding_sites)
    extended_lattice = bofs1.mof_nanoparticle(unit_cell=periodic_structure,
        combined_structure=combined_structure,
        metal_center=metal_center,
        ligand=ligand,
        bonding_sites=bonding_sites,
        target_size=125,
        target_shape='sc'
    )
