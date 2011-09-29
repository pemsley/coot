/* coot-utils/coot-coord-extras.cc
 * 
 * Copyright 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2009 by The University of Oxford
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#include <stdexcept>
#include <iomanip>

#include "string.h"

#include "coot-utils.hh"
#include "coot-coord-utils.hh"
#include "coot-coord-extras.hh"

#include "coot-sysdep.h"


// Return 0 if any of the residues don't have a dictionary entry
// geom_p gets updated to include the residue restraints if necessary
// 
std::pair<int, std::vector<std::string> >
coot::util::check_dictionary_for_residues(PCResidue *SelResidues, int nSelResidues,
					  coot::protein_geometry *geom_p,
					  int read_number) {

   std::pair<int, std::vector<std::string> > r;

   int status;
   int fail = 0; // not fail initially.

   for (int ires=0; ires<nSelResidues; ires++) { 
      std::string resname(SelResidues[ires]->name);
      status = geom_p->have_dictionary_for_residue_type(resname, read_number);
      // This bit is redundant now that try_dynamic_add has been added
      // to have_dictionary_for_residue_type():
      if (status == 0) { 
	 status = geom_p->try_dynamic_add(resname, read_number);
	 if (status == 0) {
	    fail = 1; // we failed to find it then.
	    r.second.push_back(resname);
	 }
      }
   }

   if (fail)
      r.first = 0;
   return r;
}


// For use with wiggly ligands, constructed from a minimol residue,
// the get_contact_indices_from_restraints() needs a CResidue *.
// Caller disposes.
CResidue *
coot::GetResidue(const minimol::residue &res_in) {

   CResidue *res = new CResidue;

   std::string residue_name = res_in.name; 
   int seqnum = res_in.seqnum; 
   std::string ins_code = res_in.ins_code;
   res->SetResID(residue_name.c_str(),  seqnum, ins_code.c_str());

   for (unsigned int i=0; i<res_in.atoms.size(); i++) {
      coot::minimol::atom mat = res_in.atoms[i];
      CAtom *at = new CAtom;
      at->SetAtomName(mat.name.c_str());
      at->SetElementName(mat.element.c_str());
      at->SetCoordinates(mat.pos.x(), mat.pos.y(), mat.pos.z(),
			 mat.occupancy, mat.temperature_factor);
      int new_length = mat.altLoc.length() +1;
      char *new_alt_loc = new char [new_length];
      // reset new_alt_loc
      for (unsigned int ic=0; ic<new_length; ic++)
	 new_alt_loc[ic] = 0;
      strncpy(at->altLoc, mat.altLoc.c_str(), new_length);
      res->AddAtom(at);
   }
   
   return res;
} 
      



// 200900905 These days, the following may not be needed.
// 
// We also now pass regular_residue_flag so that the indexing of the
// contacts is inverted in the case of not regular residue.  I don't
// know why this is necessary, but I have stared at it for hours, this
// is a quick (ugly hack) fix that works.  I suspect that there is
// some atom order dependency in mgtree that I don't understand.
// Please fix (remove the necessity of depending on
// regular_residue_flag) if you know how.
// 
std::vector<std::vector<int> >
coot::util::get_contact_indices_from_restraints(CResidue *residue,
						coot::protein_geometry *geom_p,
						bool regular_residue_flag,
						bool add_reverse_contacts) {

   int nResidueAtoms = residue->GetNumberOfAtoms(); 
   std::vector<std::vector<int> > contact_indices(nResidueAtoms);
   std::string restype(residue->name);
   CAtom *atom_p;

   int n_restr = geom_p->size();

   for (int icomp=0; icomp<n_restr; icomp++) {
      if ((*geom_p)[icomp].residue_info.comp_id == restype) {
// 	 std::cout << "There are " << (*geom_p)[icomp].bond_restraint.size()
// 		   << " bond restraints " << "for " << restype << std::endl;
	 for (unsigned int ibr=0; ibr< (*geom_p)[icomp].bond_restraint.size(); ibr++) {
	    for (int iat=0; iat<nResidueAtoms; iat++) {
	       atom_p = residue->GetAtom(iat);
	       std::string at_name(atom_p->GetAtomName());
	       if ( (*geom_p)[icomp].bond_restraint[ibr].atom_id_1_4c() == at_name ) {
// 		  std::cout << "found a bond match "
// 			    << (*geom_p)[icomp].bond_restraint[ibr].atom_id_1_4c()
// 			    << " to "
// 			    << (*geom_p)[icomp].bond_restraint[ibr].atom_id_2_4c()
// 			    << std::endl;
		  int ibond_to = -1;  // initially unassigned.
		  std::string at_name_2;
		  for (int iat2=0; iat2<nResidueAtoms; iat2++) {
		     atom_p = residue->GetAtom(iat2);
		     at_name_2 = atom_p->GetAtomName();
		     if ( (*geom_p)[icomp].bond_restraint[ibr].atom_id_2_4c() == at_name_2 ) {
			ibond_to = iat2;
			break;
		     }
		  }
		  if (ibond_to > -1 ) {
		     if (add_reverse_contacts == 0) { 
			if (regular_residue_flag) {
			   contact_indices[iat].push_back(ibond_to);  // for ALA etc
			} else {
			   contact_indices[ibond_to].push_back(iat);  // ligands
			   // contact_indices[iat].push_back(ibond_to);  // ALA etc
			}
		     } else {
			// add reverse contacts.
			contact_indices[ibond_to].push_back(iat);
			contact_indices[iat].push_back(ibond_to);
		     }
		  } 
//                 This spits out the names of Hydrogens, often.
//		  else
//  		     std::cout << "failed to find bonded atom "
//  			       << (*geom_p)[icomp].bond_restraint[ibr].atom_id_2_4c()
//  			       << std::endl;
	       }
	    }
	 }
      }
   }
   return contact_indices;
}

// The atoms of residue_atoms are in the "right" order for not making
// a tree along the main chain.
// 
std::vector<std::vector<int> >
coot::util::get_contact_indices_for_PRO_residue(PPCAtom residue_atoms,
						int nResidueAtoms, 
						coot::protein_geometry *geom_p) { 

   std::vector<std::vector<int> > contact_indices(nResidueAtoms);
   CAtom *atom_p;
   int n_restr = geom_p->size();
   for (int icomp=0; icomp<n_restr; icomp++) {
      if ((*geom_p)[icomp].residue_info.comp_id == "PRO") {
	 for (unsigned int ibr=0; ibr< (*geom_p)[icomp].bond_restraint.size(); ibr++) {
	    for (int iat=0; iat<nResidueAtoms; iat++) {
	       atom_p = residue_atoms[iat];
	       std::string at_name(atom_p->GetAtomName());
	       if ( (*geom_p)[icomp].bond_restraint[ibr].atom_id_1_4c() == at_name ) {
		  int ibond_to = -1;  // initially unassigned.
		  std::string at_name_2;
		  for (int iat2=0; iat2<nResidueAtoms; iat2++) {
		     atom_p = residue_atoms[iat2];
		     at_name_2 = atom_p->GetAtomName();
		     if ( (*geom_p)[icomp].bond_restraint[ibr].atom_id_2_4c() == at_name_2 ) {
			ibond_to = iat2;
			break;
		     }
		  }
		  if (ibond_to != -1)
		     contact_indices[iat].push_back(ibond_to);
	       }
	    }
	 }
      }
   }
   return contact_indices;
}


coot::util::dict_residue_atom_info_t::dict_residue_atom_info_t(const std::string &residue_name_in,
							       coot::protein_geometry *geom_p) {

   residue_name = residue_name_in;

   std::pair<short int, dictionary_residue_restraints_t> p = 
      geom_p->get_monomer_restraints(residue_name);

   if (p.first) {
      for (unsigned int iat=0; iat<p.second.atom_info.size(); iat++) {
	 std::string atom_name = p.second.atom_info[iat].atom_id_4c;
	 short int isHydrogen = 0;
	 if (p.second.atom_info[iat].type_symbol == "H" ||
	     p.second.atom_info[iat].type_symbol == "D") {
	    isHydrogen = 1;
	 }
	 atom_info.push_back(coot::util::dict_atom_info_t(atom_name, isHydrogen));
      }
   }

}

// This one we can do a dynamic add.
// 
bool
coot::util::is_nucleotide_by_dict_dynamic_add(CResidue *residue_p, coot::protein_geometry *geom_p) {

   bool is_nuc = 0;
   bool ifound = 0;
   std::string residue_name = residue_p->GetResName();

   int n_restr = geom_p->size();
   for (int icomp=0; icomp<n_restr; icomp++) {
      if ((*geom_p)[icomp].residue_info.comp_id == residue_name) {
	 ifound = 1;
	 if ((*geom_p)[icomp].residue_info.group == "RNA" ||
	     (*geom_p)[icomp].residue_info.group == "DNA" ) {
	    is_nuc = 1;
	 }
	 break;
      }
   }

   int read_number = 40;
   if (ifound == 0) {
      int status = geom_p->try_dynamic_add(residue_name, read_number);
      if (status != 0) {
	 // we successfully added it, let's try to run this function
	 // again.  Or we could just test the last entry in
	 // geom_p->dict_res_restraints(), but it is not public, so
	 // it's messy.
	 // 
	 is_nuc = is_nucleotide_by_dict_dynamic_add(residue_p, geom_p);
      } 
   }
   return is_nuc;
}


// This one we can NOT do a dynamic add.
//
bool
coot::util::is_nucleotide_by_dict(CResidue *residue_p, const coot::protein_geometry &geom) {

   bool is_nuc = 0;
   std::string residue_name = residue_p->GetResName();

   int n_restr = geom.size();
   for (int icomp=0; icomp<n_restr; icomp++) {
      if (geom[icomp].residue_info.comp_id == residue_name) {
	 if (geom[icomp].residue_info.group == "RNA" ||
	     geom[icomp].residue_info.group == "DNA" ) {
	    is_nuc = 1;
	 }
	 break;
      }
   }

   return is_nuc;
}



// Move the atoms in res_moving.
// 
// Return the number of rotated torsions.
// 
coot::match_torsions::match_torsions(CResidue *res_moving_in, CResidue *res_ref_in,
				     const coot::dictionary_residue_restraints_t &rest) {

   res_moving = res_moving_in;
   res_ref = res_ref_in;
   moving_residue_restraints = rest;

}
   
// Move the atoms in res_moving.
// 
// Return the number of rotated torsions.
//
// tr_moving is not used - we use an atom name match to find the
// torsions in the moving residue.
// 
int
coot::match_torsions::match(const std::vector <coot::dict_torsion_restraint_t>  &tr_moving,
			    const std::vector <coot::dict_torsion_restraint_t>  &tr_ref) {

   int n_matched = 0; // return value.
   coot::graph_match_info_t match_info = coot::graph_match(res_moving, res_ref, 0, 0);

   if (! match_info.success) {
      std::cout << "WARNING:: Failed to match graphs " << std::endl;
   } else { 
      std::string alt_conf = ""; // kludge it in

      // std::map<std::pair<std::string, std::string>, std::pair<std::string, std::string> > atom_name_map;

      // match_info.matching_atom_names have moving molecule first and
      // reference atoms second.  We want to map from reference atom
      // names to moving atom names.
      // 
      std::map<std::string, std::string> atom_name_map;
      for (unsigned int i=0; i<match_info.matching_atom_names.size(); i++) { 
	 atom_name_map[match_info.matching_atom_names[i].second.first] =
	    match_info.matching_atom_names[i].first.first;
	 if (0) 
	    std::cout << "      name map construction  :"
		      << match_info.matching_atom_names[i].second.first
		      << ": -> :"
		      << match_info.matching_atom_names[i].first.first
		      << ":\n";
      }

      // for debugging
      std::vector<std::pair<coot::atom_name_quad, double> > check_quads;

      for (unsigned int itr=0; itr<tr_ref.size(); itr++) {
	 coot::atom_name_quad quad_ref(tr_ref[itr].atom_id_1_4c(),
				       tr_ref[itr].atom_id_2_4c(),
				       tr_ref[itr].atom_id_3_4c(),
				       tr_ref[itr].atom_id_4_4c());
					  
	 coot::atom_name_quad quad_moving(atom_name_map[tr_ref[itr].atom_id_1_4c()],
					  atom_name_map[tr_ref[itr].atom_id_2_4c()],
					  atom_name_map[tr_ref[itr].atom_id_3_4c()],
					  atom_name_map[tr_ref[itr].atom_id_4_4c()]);

	 if (quad_ref.all_non_blank()) {
	    if (quad_moving.all_non_blank()) { 
	       
	       std::cout << "  Reference torsion: "
			 << ":" << tr_ref[itr].format() << " maps to "
			 << quad_moving << std::endl;
	       std::pair<bool, double> result = apply_torsion(quad_moving, quad_ref, alt_conf);
	       if (result.first) { 
		  n_matched++;
		  // result.second is in radians
		  std::pair<coot::atom_name_quad, double> cq (quad_moving, result.second);
		  check_quads.push_back(cq);
	       }
	    }
	 }
      }

      // after matching, check the torsions:
      for (unsigned int iquad=0; iquad<check_quads.size(); iquad++) {
	 std::pair<bool, double> mtr = get_torsion(coot::match_torsions::MOVING_TORSION,
						   check_quads[iquad].first);
	 if (mtr.first) { 
	    std::cout << "   torsion check: " << check_quads[iquad].first
		      << " should be " << std::fixed << std::setw(7) << std::setprecision(2)
		      << check_quads[iquad].second * 180/M_PI
		      << " and is "  << std::fixed << std::setw(7) << std::setprecision(2)
		      << mtr.second * 180/M_PI;
	    if (fabs(check_quads[iquad].second - mtr.second) > M_PI/180)
	       std::cout << "  ----- WRONG!!!! ";
	    std::cout << "\n";
	 }
      }
   }
   return n_matched;
}

// return in radians
std::pair<bool, double>
coot::match_torsions::get_torsion(int torsion_type,
				  const coot::atom_name_quad &quad) const {

   switch (torsion_type) {

   case coot::match_torsions::REFERENCE_TORSION:
      return get_torsion(res_ref, quad);
      
   case coot::match_torsions::MOVING_TORSION:
      return get_torsion(res_moving, quad);

   default:
      return std::pair<bool, double> (0,0);
   }
}

std::pair<bool, double>
coot::match_torsions::get_torsion(CResidue *res, const coot::atom_name_quad &quad) const {

   bool status = 0;
   double tors = 0;
   std::vector<CAtom *> atoms(4, static_cast<CAtom *> (NULL));
   atoms[0] = res->GetAtom(quad.atom1.c_str());
   atoms[1] = res->GetAtom(quad.atom2.c_str());
   atoms[2] = res->GetAtom(quad.atom3.c_str());
   atoms[3] = res->GetAtom(quad.atom4.c_str());

   if (atoms[0] && atoms[1] && atoms[2] && atoms[3]) {
      clipper::Coord_orth pts[4];
      for (unsigned int i=0; i<4; i++)
	 pts[i] = clipper::Coord_orth(atoms[i]->x, atoms[i]->y, atoms[i]->z);
      tors = clipper::Coord_orth::torsion(pts[0], pts[1], pts[2], pts[3]); // radians
      status = 1;
   }
   return std::pair<bool, double> (status, tors);
}

// Move the atoms of res_moving to match the torsion of res_ref - and
// the torsion of res_ref is determined from the
// torsion_restraint_reference atom names.
//
// The alt conf is the alt conf of the moving residue.
// 
// return the torsion which we applied (in radians).
// 
std::pair<bool, double>
coot::match_torsions::apply_torsion(const coot::atom_name_quad &moving_quad,
				    const coot::atom_name_quad &reference_quad,
				    const std::string &alt_conf) {

   bool status = 0;
   double new_angle = 0;
   std::pair<bool, double> tors = get_torsion(res_ref, reference_quad);
   if (tors.first) { 
      try {
	 coot::atom_tree_t tree(moving_residue_restraints, res_moving, alt_conf);
	 
	 new_angle = tree.set_dihedral(moving_quad.atom1, moving_quad.atom2,
				       moving_quad.atom3, moving_quad.atom4,
				       tors.second * 180/M_PI);
	 status = 1; // may not happen if set_dihedral() throws an exception
      }
      catch (std::runtime_error rte) {
	 std::cout << "WARNING setting dihedral failed, " << rte.what() << std::endl;
      } 
   }
   return std::pair<bool, double> (status, new_angle * M_PI/180.0);
}
