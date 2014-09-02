

namespace coot {
   
   bool
   compare_residue_torsions(CMMDBManager *mol1, CResidue *res_1,
			    CMMDBManager *mol2, CResidue *res_2,
			    double tollerance, // in degrees
			    coot::protein_geometry *geom_p);

   // return true if the atom names of the torsions match
   // 
   // if there are no torsionable bonds, return true
   bool compare_residue_torsion_atom_names(const std::vector<torsion_atom_quad> &tqv_1,
					   const std::vector<torsion_atom_quad> &tqv_2);
   
   
}

