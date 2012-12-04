

#include "protein-geometry.hh"
#include "coot-coord-utils.hh"
#include "simple-restraint.hh"
#include "coot-compare-residues.hh"

#include "coot-sysdep.h"

bool
coot::compare_residue_torsions(CMMDBManager *mol1, CResidue *res_1,
			       CMMDBManager *mol2, CResidue *res_2,
			       double tolerance,
			       coot::protein_geometry *geom_p) {

   bool similar_status = false;
   std::string resname_1 = res_1->GetResName();
   std::string resname_2 = res_2->GetResName();

   std::pair<bool, coot::dictionary_residue_restraints_t> restraints = 
      geom_p->get_monomer_restraints(resname_1);

   if (restraints.first) {

      PPCAtom residue_atoms_1 = 0;
      PPCAtom residue_atoms_2 = 0;
      int n_residue_atoms_1;
      int n_residue_atoms_2;
      res_1->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
      res_2->GetAtomTable(residue_atoms_2, n_residue_atoms_2);

      if (n_residue_atoms_1 && n_residue_atoms_2) { 
      
	 std::vector<torsion_atom_quad> tqv_1 = 
	    torsionable_quads(mol1, residue_atoms_1, n_residue_atoms_1, geom_p);
	 std::vector<torsion_atom_quad> tqv_2 = 
	    torsionable_quads(mol2, residue_atoms_2, n_residue_atoms_2, geom_p);

	 bool all_atom_names_match = compare_residue_torsion_atom_names(tqv_1, tqv_2);

	 if (all_atom_names_match) {
	    similar_status = true;
	    for (unsigned int i=0; i<tqv_1.size(); i++) { 
	       double t_1 = tqv_1[i].torsion();
	       double t_2 = tqv_2[i].torsion();
	       double delta = t_2 - t_1;
	       if (delta < -360.0)
		  delta += 360.0;
	       if (delta > 360.0)
		  delta -= 360.0;
	       if (delta > tolerance)
		  return false;
	       // std::cout << "   " << i << " " << t_1 << " " << t_2 << "    "
	       // << delta << std::endl;
	    }
	 }
      }
   }
   return similar_status;
}

bool
coot::compare_residue_torsion_atom_names(const std::vector<torsion_atom_quad> &tqv_1,
					 const std::vector<torsion_atom_quad> &tqv_2) {

   bool status = false;
   unsigned int n1 = tqv_1.size();
   unsigned int n2 = tqv_2.size();
   if (n1 != n2) {
      return false;
   } else {
      for (unsigned int i=0; i<tqv_1.size(); i++) {
	 if (! tqv_1[i].filled_p()) return false;
	 if (! tqv_2[i].filled_p()) return false;
	 std::string n11=tqv_1[i].atom_1->GetAtomName();
	 std::string n12=tqv_1[i].atom_2->GetAtomName();
	 std::string n13=tqv_1[i].atom_3->GetAtomName();
	 std::string n14=tqv_1[i].atom_4->GetAtomName();
	 std::string n21=tqv_2[i].atom_1->GetAtomName();
	 std::string n22=tqv_2[i].atom_2->GetAtomName();
	 std::string n23=tqv_2[i].atom_3->GetAtomName();
	 std::string n24=tqv_2[i].atom_4->GetAtomName();
	 if (! ((n11 == n21) && (n12 == n22) && (n13 == n23) && (n14 == n24))) {
	    return false;
	 } 
      }
      status = true;
   }
   return status;
}
   
