
namespace coot {
   
   std::vector<fle_residues_helper_t>
   get_flev_residue_centres(CResidue *reference_residue,
			    CMMDBManager *mol_containing_residue_ligand, 
			    std::vector<CResidue *> residues,
			    CMMDBManager *flat_mol);

}


