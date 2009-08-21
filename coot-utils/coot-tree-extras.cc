
#include "coot-coord-extras.hh"


// constructor can throw an exception
// 
// the constructor should not throw an exception if there are no
// tree in the restraints.  It should instead try the bonds in
// the restraints.  If there are no bonds then it throws an
// exception.
// 
coot::atom_tree_t::atom_tree_t(const coot::dictionary_residue_restraints_t &rest,
			       const coot::minimol::residue &res_in,
			       const std::string &altconf) {

   made_from_minimol_residue_flag = 1;

   CResidue *residue_p = new CResidue;

   std::string residue_name = res_in.name; 
   int seqnum = res_in.seqnum; 
   std::string ins_code = res_in.ins_code;
   residue->SetResID(residue_name.c_str(),  seqnum, ins_code.c_str());

   for (unsigned int i=0; i<res_in.atoms.size(); i++) {
      coot::minimol::atom mat = res_in.atoms[i];
      CAtom *at = new CAtom;
      at->SetAtomName(mat.name.c_str());
      at->SetElementName(mat.element.c_str());
      at->SetCoordinates(mat.pos.x(), mat.pos.y(), mat.pos.z(),
			 mat.occupancy, mat.temperature_factor);
      int new_length = mat.altLoc.length() +1;
      char *new_alt_loc = new char [new_length];
      strncpy(at->altLoc, mat.altLoc.c_str(), new_length);
      residue->AddAtom(at);
   }

   construct_internal(rest, residue_p, altconf);
}


// this can throw an exception
// 
// return the set of angles - should be the same that they were
// set to (for validation).
std::vector<double>
coot::atom_tree_t::set_dihedral_multi(const std::vector<tree_dihedral_info_t> &di) {

   std::vector<double> v(di.size());
   for (unsigned int id=0; id<di.size(); id++) {
      v[id] = set_dihedral(di[id].quad.atom1, di[id].quad.atom2,
			   di[id].quad.atom3, di[id].quad.atom4, di[id].dihedral_angle);
   }
   return v;
}



// For use with above function, where the class constructor is made
// with a minimol::residue.
// 
coot::minimol::residue
coot::atom_tree_t::GetResidue() const {

   return coot::minimol::residue(residue);
}
