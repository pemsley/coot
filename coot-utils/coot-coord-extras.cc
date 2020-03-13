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
#include <algorithm> // for std::find

#include "string.h"

#include "utils/coot-utils.hh"
#include "geometry/mol-utils.hh"
#include "coot-coord-utils.hh"
#include "coot-coord-extras.hh"

#include "compat/coot-sysdep.h"


// Return 0 if any of the residues don't have a dictionary entry
// geom_p gets updated to include the residue restraints if necessary
//
std::pair<int, std::vector<std::string> >
coot::util::check_dictionary_for_residues(mmdb::PResidue *SelResidues, int nSelResidues,
					  coot::protein_geometry *geom_p,
					  int read_number) {

   int imol_enc = protein_geometry::IMOL_ENC_ANY; // pass this?

   std::pair<int, std::vector<std::string> > r;

   int status;
   int fail = 0; // not fail initially.

   for (int ires=0; ires<nSelResidues; ires++) {
      std::string resname(SelResidues[ires]->name);
      status = geom_p->have_dictionary_for_residue_type(resname, imol_enc, read_number);
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
// the get_contact_indices_from_restraints() needs a mmdb::Residue *.
// Caller disposes.
mmdb::Residue *
coot::GetResidue(const minimol::residue &res_in) {

   mmdb::Residue *res = new mmdb::Residue;

   std::string residue_name = res_in.name;
   int seqnum = res_in.seqnum;
   std::string ins_code = res_in.ins_code;
   res->SetResID(residue_name.c_str(),  seqnum, ins_code.c_str());

   for (unsigned int i=0; i<res_in.atoms.size(); i++) {
      coot::minimol::atom mat = res_in.atoms[i];
      mmdb::Atom *at = new mmdb::Atom;
      at->SetAtomName(mat.name.c_str());
      at->SetElementName(mat.element.c_str());
      at->SetCoordinates(mat.pos.x(), mat.pos.y(), mat.pos.z(),
			 mat.occupancy, mat.temperature_factor);
      unsigned int new_length = mat.altLoc.length() +1;
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
coot::util::get_contact_indices_from_restraints(mmdb::Residue *residue,
						coot::protein_geometry *geom_p,
						bool regular_residue_flag,
						bool add_reverse_contacts) {

   int nResidueAtoms = residue->GetNumberOfAtoms();
   std::vector<std::vector<int> > contact_indices(nResidueAtoms);
   std::string restype(residue->name);
   mmdb::Atom *atom_p;

   int n_restr = geom_p->size();

   // this is a horrible and unconventional method of getting the restraints

   for (int icomp=0; icomp<n_restr; icomp++) {
      if ((*geom_p)[icomp].second.residue_info.comp_id == restype) {
// 	 std::cout << "There are " << (*geom_p)[icomp].bond_restraint.size()
// 		   << " bond restraints " << "for " << restype << std::endl;
	 for (unsigned int ibr=0; ibr< (*geom_p)[icomp].second.bond_restraint.size(); ibr++) {
	    for (int iat=0; iat<nResidueAtoms; iat++) {
	       atom_p = residue->GetAtom(iat);
	       std::string at_name(atom_p->GetAtomName());
	       if ( (*geom_p)[icomp].second.bond_restraint[ibr].atom_id_1_4c() == at_name ) {
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
		     if ( (*geom_p)[icomp].second.bond_restraint[ibr].atom_id_2_4c() == at_name_2 ) {
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

std::vector<std::vector<int> >
coot::util::get_contact_indices_from_restraints(mmdb::Residue *residue,
						const coot::dictionary_residue_restraints_t &restraints,
						bool regular_residue_flag,
						bool add_reverse_contacts) {

   int nResidueAtoms = residue->GetNumberOfAtoms();
   std::vector<std::vector<int> > contact_indices(nResidueAtoms);
   mmdb::Atom *atom_p;

   for (unsigned int ibr=0; ibr< restraints.bond_restraint.size(); ibr++) {
      for (int iat=0; iat<nResidueAtoms; iat++) {
	 atom_p = residue->GetAtom(iat);
	 std::string at_name(atom_p->GetAtomName());
	 if (restraints.bond_restraint[ibr].atom_id_1_4c() == at_name ) {
	    int ibond_to = -1;  // initially unassigned.
	    std::string at_name_2;
	    for (int iat2=0; iat2<nResidueAtoms; iat2++) {
	       atom_p = residue->GetAtom(iat2);
	       at_name_2 = atom_p->GetAtomName();
	       if (restraints.bond_restraint[ibr].atom_id_2_4c() == at_name_2 ) {
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
	 }
      }
   }
   return contact_indices;
}


// The atoms of residue_atoms are in the "right" order for not making
// a tree along the main chain.
//
std::vector<std::vector<int> >
coot::util::get_contact_indices_for_PRO_residue(mmdb::PPAtom residue_atoms,
						int nResidueAtoms,
						coot::protein_geometry *geom_p) {

   std::vector<std::vector<int> > contact_indices(nResidueAtoms);
   mmdb::Atom *atom_p;

   std::pair<bool, coot::dictionary_residue_restraints_t> rest =
      geom_p->get_monomer_restraints("PRO", protein_geometry::IMOL_ENC_ANY);

   if (rest.first) {
      for (unsigned int ibr=0; ibr< rest.second.bond_restraint.size(); ibr++) {
	 for (int iat=0; iat<nResidueAtoms; iat++) {
	    atom_p = residue_atoms[iat];
	    std::string at_name(atom_p->GetAtomName());
	    if (rest.second.bond_restraint[ibr].atom_id_1_4c() == at_name ) {
	       int ibond_to = -1;  // initially unassigned.
	       std::string at_name_2;
	       for (int iat2=0; iat2<nResidueAtoms; iat2++) {
		  atom_p = residue_atoms[iat2];
		  at_name_2 = atom_p->GetAtomName();
		  if (rest.second.bond_restraint[ibr].atom_id_2_4c() == at_name_2 ) {
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
   return contact_indices;
}


coot::util::dict_residue_atom_info_t::dict_residue_atom_info_t(const std::string &residue_name_in,
							       coot::protein_geometry *geom_p) {

   residue_name = residue_name_in;

   std::pair<short int, dictionary_residue_restraints_t> p =
      geom_p->get_monomer_restraints(residue_name, protein_geometry::IMOL_ENC_ANY);

   if (p.first) {
      for (unsigned int iat=0; iat<p.second.atom_info.size(); iat++) {
	 std::string atom_name = p.second.atom_info[iat].atom_id_4c;
	 short int isHydrogen = 0;
	 if (p.second.atom_info[iat].type_symbol == "H" ||
	     p.second.atom_info[iat].type_symbol == "D") {
	    isHydrogen = 1;
	 }
	 atom_info.push_back(util::dict_atom_info_t(atom_name, isHydrogen));
      }
   }

}

// This one we can do a dynamic add.
//
bool
coot::util::is_nucleotide_by_dict_dynamic_add(mmdb::Residue *residue_p, coot::protein_geometry *geom_p) {

   bool is_nuc = 0;
   std::string residue_name = residue_p->GetResName();

   std::pair<short int, dictionary_residue_restraints_t> p =
      geom_p->get_monomer_restraints(residue_name, protein_geometry::IMOL_ENC_ANY);

   if (p.first) {
      if (p.second.residue_info.group == "RNA" ||
	  p.second.residue_info.group == "DNA" ) {
	 is_nuc = 1;
      }
   } else {
      int read_number = 40;
      int status = geom_p->try_dynamic_add(residue_name, read_number);
      if (status != 0) {
	 // we successfully added it, let's try to run this function
	 // again.  Or we could just test the last entry in
	 // geom_p->dict_res_restraints(), but it is not public, so
	 // it's messy.
	 //
	 p = geom_p->get_monomer_restraints(residue_name, protein_geometry::IMOL_ENC_ANY);
	 if (p.first) {
	    if (p.second.residue_info.group == "RNA" ||
		p.second.residue_info.group == "DNA" ) {
	       is_nuc = 1;
	    }
	 }
      }
   }
   return is_nuc;
}


// This one we can NOT do a dynamic add.
//
bool
coot::util::is_nucleotide_by_dict(mmdb::Residue *residue_p, const coot::protein_geometry &geom) {

   bool is_nuc = 0;
   std::string residue_name = residue_p->GetResName();

   std::pair<short int, dictionary_residue_restraints_t> p =
      geom.get_monomer_restraints(residue_name, protein_geometry::IMOL_ENC_ANY);
   if (p.second.residue_info.group == "RNA" ||
       p.second.residue_info.group == "DNA" ) {
      is_nuc = 1;
   }
   return is_nuc;
}



// Move the atoms in res_moving.
//
// Return the number of rotated torsions.
//
coot::match_torsions::match_torsions(mmdb::Residue *res_moving_in, mmdb::Residue *res_ref_in,
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
      std::vector<double> starting_quad_torsions;

      for (unsigned int itr=0; itr<tr_ref.size(); itr++) {

	 coot::atom_name_quad quad_ref(tr_ref[itr].atom_id_1_4c(),
				       tr_ref[itr].atom_id_2_4c(),
				       tr_ref[itr].atom_id_3_4c(),
				       tr_ref[itr].atom_id_4_4c());

	 coot::atom_name_quad quad_moving(atom_name_map[tr_ref[itr].atom_id_1_4c()],
					  atom_name_map[tr_ref[itr].atom_id_2_4c()],
					  atom_name_map[tr_ref[itr].atom_id_3_4c()],
					  atom_name_map[tr_ref[itr].atom_id_4_4c()]);

	 // test if quad_moving/tr_moving (whatever that is) is a ring torsion
	 //
	 // if (moving_residue_restraints.is_ring_torsion(quad_moving))


	 if (moving_residue_restraints.is_ring_torsion(quad_moving)) {
	    // std::cout << "    ignore this ring torsion " << quad_moving << std::endl;
	 } else {
	    // std::cout << "    OK moving torsion " << quad_moving << " is not a ring torsion"
	    // << std::endl;


	    if (quad_ref.all_non_blank()) {
	       if (quad_moving.all_non_blank()) {

		  if (0)
		     std::cout << "  Reference torsion: "
			       << ":" << tr_ref[itr].format() << " maps to "
			       << quad_moving << std::endl;

		  double starting_quad_tor = quad_moving.torsion(res_moving); // debugging
		  std::pair<bool, double> result = apply_torsion(quad_moving, quad_ref, alt_conf);
		  if (! result.first) {
		     // no tree in restraints? Try without
		     result = apply_torsion_by_contacts(quad_moving, quad_ref, alt_conf);
		  }

		  if (result.first) {
		     n_matched++;
		     // result.second is in radians
		     std::pair<coot::atom_name_quad, double> cq (quad_moving, result.second);
		     check_quads.push_back(cq);
		     starting_quad_torsions.push_back(starting_quad_tor);
		  }
	       } else {
		  std::cout << "WARNING:: in torsion match() quad moving not all non-blank" << std::endl;
	       }
	    } else {
	       std::cout << "WARNING:: in torsion match() quad ref not all non-blank" << std::endl;
	    }
	 }
      }

      std::cout << "------ after matching, check the torsions " << std::endl;
      // after matching, check the torsions:
      for (unsigned int iquad=0; iquad<check_quads.size(); iquad++) {
	 std::pair<bool, double> mtr = get_torsion(coot::match_torsions::MOVING_TORSION,
						   check_quads[iquad].first);
	 if (mtr.first) {
	    std::cout << "   torsion check:  "
		      << check_quads[iquad].first << "  was "
		      << std::fixed << std::setw(7) << std::setprecision(2)
		      << starting_quad_torsions[iquad] << " "
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
coot::match_torsions::get_torsion(mmdb::Residue *res, const coot::atom_name_quad &quad) const {

   bool status = 0;
   double tors = 0;
   std::vector<mmdb::Atom *> atoms(4, static_cast<mmdb::Atom *> (NULL));
   atoms[0] = res->GetAtom(quad.atom_name(0).c_str());
   atoms[1] = res->GetAtom(quad.atom_name(1).c_str());
   atoms[2] = res->GetAtom(quad.atom_name(2).c_str());
   atoms[3] = res->GetAtom(quad.atom_name(3).c_str());

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

	 new_angle = tree.set_dihedral(moving_quad.atom_name(0), moving_quad.atom_name(1),
				       moving_quad.atom_name(2), moving_quad.atom_name(3),
				       tors.second * 180/M_PI);
	 status = 1; // may not happen if set_dihedral() throws an exception
      }
      catch (const std::runtime_error &rte) {
	 // std::cout << "WARNING tree-based setting dihedral failed, " << rte.what() << std::endl;
      }
   }
   return std::pair<bool, double> (status, new_angle * M_PI/180.0);
}


std::pair<bool, double>
coot::match_torsions::apply_torsion_by_contacts(const coot::atom_name_quad &moving_quad,
						const coot::atom_name_quad &reference_quad,
						const std::string &alt_conf) {

   bool status = false;
   double new_angle = 0.0;

   try {
      bool add_reverse_contacts = true;
      std::vector<std::vector<int> > contact_indices =
	 coot::util::get_contact_indices_from_restraints(res_moving, moving_residue_restraints, 1, add_reverse_contacts);
      std::pair<bool, double> tors = get_torsion(res_ref, reference_quad);

      int base_atom_index = 0; // hopefully this will work
      coot::minimol::residue ligand_residue(res_moving);
      coot::atom_tree_t tree(moving_residue_restraints, contact_indices, base_atom_index, ligand_residue, alt_conf);
      new_angle = tree.set_dihedral(moving_quad.atom_name(0), moving_quad.atom_name(1),
				    moving_quad.atom_name(2), moving_quad.atom_name(3),
				    tors.second * 180/M_PI);

      coot::minimol::residue wiggled_ligand_residue = tree.GetResidue();

      // transfer from the ligand_residue to res_moving
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      int n_transfered = 0;
      res_moving->GetAtomTable(residue_atoms, n_residue_atoms);
      if (int(wiggled_ligand_residue.atoms.size()) <= n_residue_atoms) {
	 for (unsigned int iat=0; iat<wiggled_ligand_residue.atoms.size(); iat++) {
	    mmdb::Atom *at = res_moving->GetAtom(wiggled_ligand_residue.atoms[iat].name.c_str(), NULL, alt_conf.c_str());
	    if (at) {
	       if (0)
		  std::cout << "transfering coords was "
			    << at->z << " " << at->y << " " << at->z << " to "
			    << ligand_residue.atoms[iat] << std::endl;
	       at->x = wiggled_ligand_residue.atoms[iat].pos.x();
	       at->y = wiggled_ligand_residue.atoms[iat].pos.y();
	       at->z = wiggled_ligand_residue.atoms[iat].pos.z();
	       n_transfered++;
	    }
	 }
      }
      if (0) {
	 std::cout << "-------------------------------- n_transfered " << n_transfered << "------------------- "
		   << std::endl;
	 std::cout << "in apply_torsion_by_contacts() new_angle is " << new_angle << std::endl;
      }
      status = 1;
   }
   catch (const std::runtime_error &rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   }
   return std::pair<bool, double> (status, new_angle * M_PI/180.0);
}


// Don't return any hydrogen torsions - perhaps we should make that a
// passed parameter.
//
std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> >
coot::torsionable_bonds_monomer_internal(mmdb::Residue *residue_p,
					 mmdb::PPAtom atom_selection, int n_selected_atoms,
					 bool include_pyranose_ring_torsions_flag,
					 coot::protein_geometry *geom_p) {

   std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > v;

   bool hydrogen_torsions = false;
   std::string rn = residue_p->GetResName();
   std::vector <dict_torsion_restraint_t> tors_restraints =
      geom_p->get_monomer_torsions_from_geometry(rn, hydrogen_torsions);
   bool is_pyranose = false; // reset maybe
   std::string group = geom_p->get_group(residue_p);
   // CCD dictionaries are marked at D-SACCHARIDE or SACCHARIDE (FUC)
   // SIA (sialic acid) is marked as NON-POLYMER.  Hmm..
   if (group == "pyranose" || group == "D-pyranose" || group == "L-pyranose" ||
       group == "D-SACCHARIDE" || group == "SACCHARIDE")
      is_pyranose = true;

   if (tors_restraints.size()) {
      for (unsigned int itor=0; itor<tors_restraints.size(); itor++) {

	 if (! tors_restraints[itor].is_const()) {
	    std::string tr_atom_name_2 = tors_restraints[itor].atom_id_2_4c();
	    std::string tr_atom_name_3 = tors_restraints[itor].atom_id_3_4c();

	    for (int iat1=0; iat1<n_selected_atoms; iat1++) {
	       mmdb::Residue *res_1 = atom_selection[iat1]->residue;
	       std::string atom_name_1 = atom_selection[iat1]->name;
	       for (int iat2=0; iat2<n_selected_atoms; iat2++) {
		  if (iat1 != iat2) {
		     mmdb::Residue *res_2 = atom_selection[iat2]->residue;
		     if (res_1 == res_2) {
			std::string atom_name_2 = atom_selection[iat2]->name;
			if (atom_name_1 == tr_atom_name_2) {
			   if (atom_name_2 == tr_atom_name_3) {

			      if ((include_pyranose_ring_torsions_flag == 1) ||
				  (is_pyranose && !tors_restraints[itor].is_pyranose_ring_torsion(rn)) ||
				  (! is_pyranose)) {

				 std::pair<mmdb::Atom *, mmdb::Atom *> p(atom_selection[iat1],
							       atom_selection[iat2]);
				 v.push_back(p);
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return v;
}

// The quad version of this (for actually setting torsions)
//
// Don't return any hydrogen torsions - perhaps we should make that a
// passed parameter.
//
std::vector<coot::torsion_atom_quad>
coot::torsionable_bonds_monomer_internal_quads(mmdb::Residue *residue_p,
					       mmdb::PPAtom atom_selection, int n_selected_atoms,
					       bool include_pyranose_ring_torsions_flag,
					       coot::protein_geometry *geom_p) {
   std::vector<coot::torsion_atom_quad> quads;
   bool hydrogen_torsions = false;
   std::string rn = residue_p->GetResName();
   std::vector <dict_torsion_restraint_t> tors_restraints =
      geom_p->get_monomer_torsions_from_geometry(rn, hydrogen_torsions);
   bool is_pyranose = false;
   std::string group = geom_p->get_group(residue_p);
   if (group == "pyranose" || group == "D-pyranose" || group == "L-pyranose" ||
       group == "D-SACCHARIDE" || group == "SACCHARIDE")
      is_pyranose = true;
   std::vector<std::string> residue_alt_confs = coot::util::get_residue_alt_confs(residue_p);
   for (unsigned int itor=0; itor<tors_restraints.size(); itor++) {
      if (! tors_restraints[itor].is_const()) {
	 std::string tor_atom_name[5];
	 std::vector<mmdb::Atom *> ats(5, static_cast<mmdb::Atom *>(NULL));
	 tor_atom_name[1] = tors_restraints[itor].atom_id_1_4c();
	 tor_atom_name[2] = tors_restraints[itor].atom_id_2_4c();
	 tor_atom_name[3] = tors_restraints[itor].atom_id_3_4c();
	 tor_atom_name[4] = tors_restraints[itor].atom_id_4_4c();
	 if ((include_pyranose_ring_torsions_flag == 1) ||
	     (is_pyranose && !tors_restraints[itor].is_pyranose_ring_torsion(rn)) ||
	     (! is_pyranose)) {
	    for (unsigned int ialt=0; ialt<residue_alt_confs.size(); ialt++) {
	       for (int iat=0; iat<n_selected_atoms; iat++) {
		  std::string atom_name = atom_selection[iat]->name;
		  std::string alt_conf  = atom_selection[iat]->altLoc;
		  if (alt_conf == residue_alt_confs[ialt]) {
		     for (unsigned int jtor=1; jtor<5; jtor++) {
			if (atom_name == tor_atom_name[jtor])
			   ats[jtor] = atom_selection[iat];
		     }
		  }
	       }
	       // yes we got for atoms (of matching alt confs)
	       if (ats[1] && ats[2] && ats[3] && ats[4]) {
		  coot::torsion_atom_quad q(ats[1],ats[2],ats[3],ats[4],
					    tors_restraints[itor].angle(),
					    tors_restraints[itor].esd(),
					    tors_restraints[itor].periodicity());
		  q.name = tors_restraints[itor].id();
		  q.residue_name = rn;
		  quads.push_back(q);

	       }
	    }
	 }
      }
   }
   return quads;
}


coot::bonded_pair_container_t
coot::linkrs_in_atom_selection(mmdb::Manager *mol, mmdb::PPAtom atom_selection, int n_selected_atoms,
			       protein_geometry *geom_p) {
   coot::bonded_pair_container_t bpc;
#ifdef MMDB_WITHOUT_LINKR
#else
   // normal case
   std::vector<mmdb::Residue *> residues;
   for (int i=0; i<n_selected_atoms; i++) {
      mmdb::Residue *r = atom_selection[i]->residue;
      if (std::find(residues.begin(), residues.end(), r) == residues.end())
	 residues.push_back(r);
   }

   bool found = false;
   mmdb::Model *model_p = mol->GetModel(1);
   int n_linkrs = model_p->GetNumberOfLinkRs();
   std::cout << "model has " << n_linkrs << " LINKR records"
	     << " and " << model_p->GetNumberOfLinks() << " LINK records"
	     << std::endl;
   for (int ilink=1; ilink<=n_linkrs; ilink++) {
      mmdb::PLinkR linkr = model_p->GetLinkR(ilink);
      coot::residue_spec_t link_spec_1(linkr->chainID1,
				       linkr->seqNum1,
				       linkr->insCode1);
      coot::residue_spec_t link_spec_2(linkr->chainID2,
				       linkr->seqNum2,
				       linkr->insCode2);
      for (unsigned int i=0; i<residues.size(); i++) {
	 coot::residue_spec_t spec_1(residues[i]);
	 if (spec_1 == link_spec_1) {
	    for (unsigned int j=0; j<residues.size(); j++) {
	       if (i != j) {
		  coot::residue_spec_t spec_2(residues[j]);
		  if (spec_2 == link_spec_2) {
		     found = true;
		     coot::bonded_pair_t pair(residues[i], residues[j], 0, 0, linkr->linkRID);
		     break;
		  }
	       }
	    }
	 }
	 if (found)
	    break;
      }
      if (found)
	 break;
   }

#endif
   return bpc;
}


// use residues-near-residue to find linked residues
std::vector<mmdb::Residue *>
coot::simple_residue_tree(mmdb::Residue *residue_centre, mmdb::Manager *mol, float close_dist_max) {

   double dist_crit = close_dist_max;
   std::vector<mmdb::Residue *> v;
   std::set<mmdb::Residue *> s;

   std::queue<mmdb::Residue *> q; // what is dequeue? (double-ended)

   q.push(residue_centre);
   s.insert(residue_centre);

   while (q.size()) {
      mmdb::Residue *test_residue = q.front();
      s.insert(test_residue);
      q.pop();

      // OK, what new ones shall we add?
      // Don't add residues that are already in the set. Everything that is in the queue
      // is in the set also.
      std::vector<mmdb::Residue *> residues = residues_near_residue(test_residue, mol, dist_crit);
      for (unsigned int ires=0; ires<residues.size(); ires++) {
	 mmdb::Residue *rnr = residues[ires];
	 std::set<mmdb::Residue *>::const_iterator it = s.find(rnr);
	 if (it == s.end()) {
	    q.push(rnr);
	    s.insert(rnr);
	 }
      }
   }

   std::set<mmdb::Residue *>::const_iterator its;
   for (its=s.begin(); its!=s.end(); its++)
      v.push_back(*its);

   return v;
}

// Is this the wrong file for this util function?
std::vector<std::pair<mmdb::Residue *, double> >
coot::util::CO_orientations(mmdb::Manager *mol) {

   std::vector<std::pair<mmdb::Residue *, double> > scores;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         if (nres > 2) {
            int res_idx_min = 1;
            int res_idx_max = nres-2;
            for (int ires=res_idx_min; ires<res_idx_max; ires++) {
               mmdb::Residue *prev_p = chain_p->GetResidue(ires-1);
               mmdb::Residue *this_p = chain_p->GetResidue(ires);
               mmdb::Residue *next_p = chain_p->GetResidue(ires+1);
               int delta_1 = this_p->GetSeqNum() - prev_p->GetSeqNum();
               int delta_2 = next_p->GetSeqNum() - this_p->GetSeqNum();
               if (delta_1 != 1) continue;
               if (delta_2 != 1) continue;
               if (prev_p && this_p && next_p) {
                  int n_atoms_prev = prev_p->GetNumberOfAtoms();
                  int n_atoms_this = this_p->GetNumberOfAtoms();
                  int n_atoms_next = next_p->GetNumberOfAtoms();
                  mmdb::Atom *prev_O = 0;
                  mmdb::Atom *prev_C = 0;
                  mmdb::Atom *this_O = 0;
                  mmdb::Atom *this_C = 0;
                  mmdb::Atom *next_O = 0;
                  mmdb::Atom *next_C = 0;
                  for (int iat=0; iat<n_atoms_prev; iat++) {
                     mmdb::Atom *at = prev_p->GetAtom(iat);
                     std::string atom_name(at->GetAtomName());
                     std::string alt_conf(at->altLoc);
                     if (alt_conf == "") {
                        if (atom_name == " C  ") prev_C = at;
                        if (atom_name == " O  ") prev_O = at;
                     }
                  }
                  if (! prev_C) continue;
                  if (! prev_O) continue;
                  for (int iat=0; iat<n_atoms_this; iat++) {
                     mmdb::Atom *at = this_p->GetAtom(iat);
                     std::string atom_name(at->GetAtomName());
                     std::string alt_conf(at->altLoc);
                     if (alt_conf == "") {
                        if (atom_name == " C  ") this_C = at;
                        if (atom_name == " O  ") this_O = at;
                     }
                  }
                  if (! this_C) continue;
                  if (! this_O) continue;
                  for (int iat=0; iat<n_atoms_next; iat++) {
                     mmdb::Atom *at = next_p->GetAtom(iat);
                     std::string atom_name(at->GetAtomName());
                     std::string alt_conf(at->altLoc);
                     if (alt_conf == "") {
                        if (atom_name == " C  ") next_C = at;
                        if (atom_name == " O  ") next_O = at;
                     }
                  }
                  if (! next_C) continue;
                  if (! next_O) continue;

                  // OK!
                  clipper::Coord_orth v1(co(prev_O) - co(prev_C));
                  clipper::Coord_orth v2(co(this_O) - co(this_C));
                  clipper::Coord_orth v3(co(next_O) - co(next_C));

                  clipper::Coord_orth v1n(v1.unit());
                  clipper::Coord_orth v2n(v2.unit());
                  clipper::Coord_orth v3n(v3.unit());

                  double dp_1 = clipper::Coord_orth::dot(v1n, v2n);
                  double dp_2 = clipper::Coord_orth::dot(v2n, v3n);
                  double sum = dp_1; // + dp_2;
                  // std::cout << "dp_1 " << dp_1 << " dp_2 " << dp_2 << "\n";
                  std::pair<mmdb::Residue *, double> s(this_p, sum);
                  scores.push_back(s);
               }
            }
         }
      }
   }
   return scores;
}

#include "atom-selection-container.hh"
void
coot::util::parse_prosmart_log_and_gen_CO_plot(const std::string &prosmart_log_file_helix,
                                               const std::string &prosmart_log_file_strand,
                                               const std::string &data_points_file_name,
                                               const std::string &pdb_file_name,
                                               const std::string &chain_id) {

   atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, false);
   if (asc.read_success) {
      std::vector<std::pair<mmdb::Residue *, double> > co_scores = CO_orientations(asc.mol);
      std::map<residue_spec_t, double> scores_map;
      std::vector<std::pair<mmdb::Residue *, double> >::const_iterator it;
      for (it=co_scores.begin(); it!=co_scores.end(); it++) {
         residue_spec_t spec(it->first);
         scores_map[spec] = it->second;
      }
      if (co_scores.size() > 0) {
         std::ifstream f_helix(prosmart_log_file_helix.c_str());
         std::ifstream f_strand(prosmart_log_file_strand.c_str());
         std::ofstream fo(data_points_file_name.c_str());
         std::string pdb_fn = file_name_non_directory(pdb_file_name);
         std::map<residue_spec_t, double> helix_scores;
         std::map<residue_spec_t, double> strand_scores;
         if (f_helix) {
            std::string line;
            while (std::getline(f_helix, line)) {
               std::vector<std::string> bits = util::split_string_on_whitespace_no_blanks(line);
               // std::cout << "read " << line << " from " << prosmart_log_file << " "
               //          << bits.size() << std::endl;
               if (bits.size() == 12) {
                  if (bits[3] == "ALA") { // does this work for strand test also?
                     try {
                        int res_no = string_to_int(bits[0]);
                        float deviation_flexi      = string_to_float(bits[6]);
                        float deviation_procrustes = string_to_float(bits[7]);
                        residue_spec_t res_spec(chain_id, res_no, "");
                        helix_scores[res_spec] = deviation_procrustes;
                     }
                     catch (const std::runtime_error &rte) {
                        // residue number was not a number - oh well
                        std::cout << "something bad parsing " << line  << " " << rte.what()
                                  << std::endl;
                     }
                  }
               }
            }
         }

         if (f_strand) {
            std::string line;
            while (std::getline(f_strand, line)) {
               std::vector<std::string> bits = util::split_string_on_whitespace_no_blanks(line);
               // std::cout << "read " << line << " from " << prosmart_log_file << " "
               //          << bits.size() << std::endl;
               if (bits.size() == 12) {
                  if (bits[3] == "ALA") { // does this work for strand test also?
                     try {
                        int res_no = string_to_int(bits[0]);
                        float deviation_flexi      = string_to_float(bits[6]);
                        float deviation_procrustes = string_to_float(bits[7]);
                        residue_spec_t res_spec(chain_id, res_no, "");
                        strand_scores[res_spec] = deviation_procrustes;
                     }
                     catch (const std::runtime_error &rte) {
                        // residue number was not a number - oh well
                        std::cout << "something bad parsing " << line  << " " << rte.what()
                                  << " " << prosmart_log_file_strand << std::endl;
                     }
                  }
               }
            }
         }

         for (it=co_scores.begin(); it!=co_scores.end(); it++) {
            residue_spec_t res_spec(it->first);
            const double &co_dp(it->second);
            std::map<residue_spec_t, double>::const_iterator it_helix;
            std::map<residue_spec_t, double>::const_iterator it_strand;
            it_strand = strand_scores.find(res_spec);
            it_helix  = helix_scores.find(res_spec);
            if (it_helix != helix_scores.end()) {
               if (it_strand != strand_scores.end()) {
                  const double &helix_score(it_helix->second);
                  const double &strand_score(it_strand->second);
                  fo << pdb_fn << " " << chain_id << " " << res_spec.res_no
                     << " CO-dp: " << co_dp
                     << " helix: " << helix_score
                     << " strand: " << strand_score
                     << "\n";
               }
            } else {
               if (false) {
                  std::cout << "debug:: failed to find residue " << res_spec
                            << " in " << " helix map of size " << helix_scores.size() << std::endl;
                  std::map<residue_spec_t, double>::const_iterator it;
                  for (it=helix_scores.begin(); it!=helix_scores.end(); it++) {
                     std::cout << "   " << it->first << " " << it->second << std::endl;
                  }
               }
            }
         }
      }
   }
}

void
coot::util::multi_parse_prosmart_log_and_gen_CO_plot() {
   std::string dir = "pdb";
   std::vector<std::string> sub_dirs = glob_files(dir, "*");
   for (unsigned int i=0; i<sub_dirs.size(); i++) {
      const std::string &sub_dir = sub_dirs[i];
      std::vector<std::string> sub_dir_files = glob_files(sub_dir, "*.pdb");
      for (unsigned int j=0; j<sub_dir_files.size(); j++) {
         const std::string &pdb_file = sub_dir_files[j];
         // std::cout << "pdb_file: " << pdb_file << std::endl;
         std::string code_pdb = util::file_name_non_directory(pdb_file);
         if (code_pdb.length() > 4) code_pdb = code_pdb.substr(0,4);
         std::string code = util::name_sans_extension(code_pdb);
         std::string code_star = code + "_*";
         // std::cout << "pdb file name: " << pdb_file << " code: " << code << "\n";

         std::string prosmart_1 = "../src/ProSMART_Output/Output_Files/Residue_Alignment_Scores";
         std::vector<std::string> chain_files = glob_files(prosmart_1, code_star);

         for (unsigned int k=0; k<chain_files.size(); k++) {
            std::string chain_file = chain_files[k];
            // ha! now I need to strip the directory
            std::string chain_file_file = util::file_name_non_directory(chain_file);
            // std::cout << "chain file file:  " << chain_file_file << std::endl;
            int cff_len = chain_file_file.length();
            if (cff_len == 6 || cff_len == 14) {
               std::string chain_id(chain_file_file.substr(5,1));
               if (cff_len == 14) chain_id = chain_file_file.substr(13,1);
               // std::cout << "chain file " << chain_file << " chain-id: " << chain_id << std::endl;
               std::string data_file_name = chain_file_file + ".data";
               std::string fn_helix  = chain_file_file + "_helix_A.txt";
               std::string fn_strand = chain_file_file + "_strand_A.txt";
               std::string log_file_name_helix  = append_dir_file(chain_file, fn_helix);
               std::string log_file_name_strand = append_dir_file(chain_file, fn_strand);

               parse_prosmart_log_and_gen_CO_plot(log_file_name_helix, log_file_name_strand,
                                                  data_file_name, pdb_file, chain_id);
            }
         }
      }
   }
}
