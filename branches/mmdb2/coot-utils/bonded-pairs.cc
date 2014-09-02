

#include "bonded-pairs.hh"

bool
coot::bonded_pair_container_t::linked_already_p(CResidue *r1, CResidue *r2) const {

   bool r = 0;
   for (unsigned int i=0; i<bonded_residues.size(); i++) {
      if (((bonded_residues[i].res_1 == r1) &&
	   (bonded_residues[i].res_2 == r2)) ||
	  ((bonded_residues[i].res_1 == r2) &&
	   (bonded_residues[i].res_2 == r1))) {
	 r = 1;
	 break;
      }
   }
   return r;
}

bool
coot::bonded_pair_container_t::try_add(const coot::bonded_pair_t &bp) {

   bool found = 0;
   for (unsigned int i=0; i<bonded_residues.size(); i++) {
      if ( (bonded_residues[i].res_1 == bp.res_1 &&
	    bonded_residues[i].res_2 == bp.res_2) ||
	   (bonded_residues[i].res_1 == bp.res_2 &&
	    bonded_residues[i].res_2 == bp.res_1) ) {
	 found = 1;
	 break;
      }
   }
   
   if (! found) {
      bonded_residues.push_back(bp);
   }
   return found;
}

std::ostream&
coot::operator<<(std::ostream &s, coot::bonded_pair_container_t bpc) {

   s << "Bonded Pair Container contains " << bpc.bonded_residues.size() << " bonded residues"
     << "\n";

   for (unsigned int i=0; i<bpc.bonded_residues.size(); i++)
      s << "   " << i << "  [:"
	<< bpc[i].link_type << ": "
	<< bpc[i].res_1->GetChainID() << " " << bpc[i].res_1->GetSeqNum() << " "
	<< bpc[i].res_1->GetInsCode() << " to " << bpc[i].res_2->GetChainID() << " "
	<< bpc[i].res_2->GetSeqNum() << " " << bpc[i].res_2->GetInsCode() << "]"
	<< "\n";

   return s; 
}

std::ostream&
coot::operator<<(std::ostream &s, coot::bonded_pair_t bp) {
   s << "[:";
   if (bp.res_1)
      s << bp.res_1->GetChainID() <<  " " << bp.res_1->GetSeqNum() << " " << bp.res_1->GetInsCode();
   s << " ";
   if (bp.res_2)
      s << bp.res_2->GetChainID() <<  " " << bp.res_2->GetSeqNum() << " " << bp.res_2->GetInsCode();
   s << "]";
   return s;
}


void
coot::bonded_pair_t::apply_chem_mods(const coot::protein_geometry &geom) {

   if (res_2 && res_2) { 
      try { 
	 // apply the mods given the link type

	 // get the chem mods for each residue (can throw a runtime
	 // error if there is one - (not an error).
	 // 
	 std::pair<coot::protein_geometry::chem_mod, coot::protein_geometry::chem_mod>
	    mods = geom.get_chem_mods_for_link(link_type);
	 std::string res_1_name = res_1->GetResName();
	 std::string res_2_name = res_2->GetResName();
	 for (unsigned int i=0; i<mods.first.atom_mods.size(); i++) {
	    if (0)
	       std::cout << "  ====== applying chem_mod " << i << " of "
			 << mods.first.atom_mods.size() << std::endl;
	    if (mods.first.atom_mods[i].function == CHEM_MOD_FUNCTION_DELETE) {
	       std::string atom_name = mods.first.atom_mods[i].atom_id;
	       std::string at_name = geom.atom_id_expand(atom_name, res_1_name);
	       delete_atom(res_1, at_name);
	    }
	 }
	 for (unsigned int i=0; i<mods.second.atom_mods.size(); i++) {
	    if (0)
	       std::cout << "  ====== applying chem_mod " << i << " of "
			 << mods.second.atom_mods.size() << std::endl;
	    if (mods.second.atom_mods[i].function == CHEM_MOD_FUNCTION_DELETE) {
	       std::string atom_name = mods.second.atom_mods[i].atom_id;
	       std::string at_name = geom.atom_id_expand(atom_name, res_2_name);
	       delete_atom(res_2, at_name);
	    }
	 }
      }
      catch (std::runtime_error rte) {
	 // it's OK if we don't find a chem mod for this link
      }
   }
}

void
coot::bonded_pair_container_t::apply_chem_mods(const protein_geometry &geom) {

   std::vector<coot::bonded_pair_t>::iterator it;
   for (it=bonded_residues.begin(); it != bonded_residues.end(); it++) {
      it->apply_chem_mods(geom);
   }
} 


void
coot::bonded_pair_t::delete_atom(CResidue *res, const std::string &atom_name) {
   
   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   bool deleted = false;
   res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
      CAtom *at = residue_atoms[iat];
      if (at) {  // unneeded precaution?
	 std::string at_name(at->name);
	 if (at_name == atom_name) {
	    delete at;
	    at = NULL;
	    deleted = true;
	 }
      }
   }

   if (deleted)
      res->TrimAtomTable();

}
