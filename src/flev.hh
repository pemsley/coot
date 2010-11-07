
namespace coot {
   
   std::vector<fle_residues_helper_t>
   get_flev_residue_centres(CResidue *reference_residue,
			    CMMDBManager *mol_containing_residue_ligand, 
			    std::vector<CResidue *> residues,
			    CMMDBManager *flat_mol);

   // return a vector of the same size as filtered_residues.
   // 
   std::vector<int> make_add_reps_for_near_residues(std::vector<CResidue *> filtered_residues,
						    int imol);

}


