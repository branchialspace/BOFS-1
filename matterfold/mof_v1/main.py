import matterfold as mf


if __name__ == "__main__":

    metal_center = mf.generate_metal_polyhedron(species= 'Bi', num_atoms = 6)
    ligand, mol = mf.generate_ligand('C(=C/C(=O)[O-])\C(=O)[O-]') #       [Te]CC#CC#CC[Te]         C1=C(S(=O)(=O)O)C(C2=CC(=CC(=C2)C(=S)S)C(=S)S)=CC(S(=O)(=O)O)=C1C3=CC(=CC(=C3)C(=S)S)C(=S)S
    ligand_electron_analysis = mf.ligand_electron_analysis(ligand)
    bonding_sites = mf.ligand_bonding_sites(ligand=ligand, ligand_electron_analysis=ligand_electron_analysis)
    combined_structure = mf.ligand_metal_docking(ligand=ligand, metal_center=metal_center, bonding_sites=bonding_sites)
    periodic_structure = mf.mof_cell(combined_structure, metal_center, ligand, bonding_sites)
    extended_lattice = mf.mof_nanoparticle(unit_cell=periodic_structure,
        combined_structure=combined_structure,
        metal_center=metal_center,
        ligand=ligand,
        bonding_sites=bonding_sites,
        target_size=125,
        target_shape='sc'
    )
