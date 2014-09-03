
namespace coot {
   
   std::vector<fle_residues_helper_t>
   get_flev_residue_centres(mmdb::Residue *reference_residue,
			    CMMDBManager *mol_containing_residue_ligand, 
			    std::vector<mmdb::Residue *> residues,
			    CMMDBManager *flat_mol);

   // return a vector of the same size as filtered_residues.
   // 
   std::vector<int> make_add_reps_for_near_residues(std::vector<mmdb::Residue *> filtered_residues,
						    int imol);

   void add_animated_ligand_interactions(int imol,
					 const std::vector<fle_ligand_bond_t> &ligand_bonds); 

}


