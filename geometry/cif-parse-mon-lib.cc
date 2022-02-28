/* geometry/dict-utils.cc
 * 
 * Copyright 2014, 2015 by Medical Research Council
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

#include <map>
#include <algorithm>
#include <iomanip> // setw()
#include "utils/coot-utils.hh"
#include "protein-geometry.hh"
#include "dict-mismatches.hh"
#include "dict-utils.hh"

void
coot::protein_geometry::mon_lib_add_chem_comp(const std::string &comp_id,
					      int imol_enc,
					      const std::string &three_letter_code,
					      const std::string &name,
					      const std::string &group,
					      int number_atoms_all, int number_atoms_nh,
					      const std::string &description_level) {


   if (false)
      std::cout << "DEBUG:: in mon_lib_add_chem_comp :"
		<< comp_id << ": imol_enc " << imol_enc << " :"
		<< three_letter_code << ": :"
		<< name << ": :" << group << ": :"
		<< description_level << ": :" << number_atoms_all << ": :"
		<< number_atoms_nh << std::endl;

   coot::dict_chem_comp_t ri(comp_id, three_letter_code, name, group,
			     number_atoms_all, number_atoms_nh,
			     description_level);
   bool ifound = false;

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
	 if (dict_res_restraints[i].first == imol_enc) {

	    if (false)
	       std::cout << "DEBUG:: matched imol_enc " << imol_enc
			 << " comparing read numbers: "
			 << dict_res_restraints[i].second.read_number << " "
			 << read_number << std::endl;
	 
	    if (dict_res_restraints[i].second.read_number == read_number) { 
	       ifound = true;
	       dict_res_restraints[i].second.residue_info = ri;
	       break;
	    } else {
	       // trash the old one then
	       // Message for Kevin.
	       std::cout << "INFO:: clearing old restraints for \"" << comp_id << "\"" << std::endl;
	       dict_res_restraints[i].second.clear_dictionary_residue();
	    }
	 }
      }
   }

   if (! ifound) {
      // std::cout << "DEBUG:: residue not found in mon_lib_add_chem_comp" << std::endl;
      dictionary_residue_restraints_t rest(comp_id, read_number);
      std::pair<int, dictionary_residue_restraints_t> p(imol_enc, rest);
      dict_res_restraints.push_back(p);
      dict_res_restraints[dict_res_restraints.size()-1].second.residue_info = ri;
   }
}


void
coot::protein_geometry::mon_lib_add_atom(const std::string &comp_id,
					 int imol_enc,
					 const std::string &atom_id,
					 const std::string &atom_id_4c,
					 const std::string &type_symbol,
					 const std::string &type_energy,
					 const std::pair<bool, mmdb::realtype> &partial_charge,
					 const std::pair<bool, int> &formal_charge,
					 dict_atom::aromaticity_t arom_in,
					 const std::pair<bool, clipper::Coord_orth> &model_pos,
					 const std::pair<bool, clipper::Coord_orth> &model_pos_ideal) { 

   // Are you sure that this is the version of mon_lib_add_atom() that you want?

   // debugging
   bool debug = false;
   
   if (debug) {
      std::cout << "   mon_lib_add_atom  " << comp_id << " atom-id:" << atom_id << ": :"
		<< atom_id_4c << ": " << type_symbol << " " << type_energy << " ( "
		<< partial_charge.first << "," << partial_charge.second << ")";
      std::cout << " model-pos: " << model_pos.first << " ";
      if (model_pos.first)
	 std::cout << "( "
		   << model_pos.second.x() << " "
		   << model_pos.second.y() << " "
		   << model_pos.second.z() << " ) ";
      std::cout << "model-pos-ideal: " << model_pos_ideal.first << " ";
      if (model_pos_ideal.first)
	 std::cout << "( "
		   << model_pos_ideal.second.x() << " "
		   << model_pos_ideal.second.y() << " "
		   << model_pos_ideal.second.z() << " ) ";
      std::cout << std::endl;
   }

   coot::dict_atom at_info(atom_id, atom_id_4c, type_symbol, type_energy, partial_charge);
   at_info.aromaticity = arom_in;
   at_info.formal_charge = formal_charge;
   
   if (debug) {
      std::cout << "   mon_lib_add_atom model_pos       " << model_pos.first << " "
		<< model_pos.second.format() << std::endl;
      std::cout << "   mon_lib_add_atom model_pos_ideal " << model_pos_ideal.first << " "
		<< model_pos_ideal.second.format() << std::endl;
   }
   
   if (model_pos.first)
      at_info.add_pos(coot::dict_atom::REAL_MODEL_POS, model_pos);
   if (model_pos_ideal.first)
      at_info.add_pos(coot::dict_atom::IDEAL_MODEL_POS, model_pos_ideal);

   bool ifound = 0;
   int this_index = -1; // unset

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
//       std::cout << "comparing comp_ids: :" << dict_res_restraints[i].comp_id
// 		<< ":  :" << comp_id << ":" << std::endl;
      
      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
	 if (dict_res_restraints[i].first == imol_enc) {

	    // 	 std::cout << "comparing read numbers: "
	    // 		   << dict_res_restraints[i].read_number << " and "
	    // 		   << read_number << std::endl;
	    if (dict_res_restraints[i].second.read_number == read_number) { 
	       ifound = true;
	       this_index = i;
	       dict_res_restraints[i].second.atom_info.push_back(at_info);
	       break;
	    } else {
	       // trash the old one then
	       std::cout << "######## trash the old one " << comp_id << std::endl;
	       dict_res_restraints[i].second.clear_dictionary_residue();
	    }
	 }
      }
   }

   if (! ifound) {
      // std::cout << "residue not found in mon_lib_add_atom" << std::endl;
      dictionary_residue_restraints_t rest(comp_id, read_number);
      std::pair<int, dictionary_residue_restraints_t> p(imol_enc, rest);
      dict_res_restraints.push_back(p);
      this_index = dict_res_restraints.size()-1;
      dict_res_restraints[this_index].second.atom_info.push_back(at_info);
   }

   if (debug) {
      std::cout << "   dictionary for " << dict_res_restraints[this_index].second.residue_info.comp_id
		<< " now contains " << dict_res_restraints[this_index].second.atom_info.size()
		<< " atoms" << std::endl;
      for (unsigned int i=0; i<dict_res_restraints[this_index].second.atom_info.size(); i++) { 
	 // 	  std::cout << "  " << i << "  " << dict_res_restraints[this_index].atom_info[i]
	 // 		    << std::endl;
      }
   }
}


void
coot::protein_geometry::mon_lib_add_atom(const std::string &comp_id,
					 int imol_enc,
					 const coot::dict_atom &atom_info) {

   bool debug  = false;
   bool ifound = false;
   int this_index = -1; // not found (not very useful)

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (debug)
	 std::cout << "comparing restraints [" << i << "] \""
		   << dict_res_restraints[i].second.residue_info.comp_id
		   << "\" with \"" << comp_id << "\"" << std::endl;
      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
	 if (dict_res_restraints[i].first == imol_enc) {
	    if (dict_res_restraints[i].second.read_number == read_number) {
	       ifound = true;
	       this_index = i;
	       dict_res_restraints[i].second.atom_info.push_back(atom_info);
	       break;
	    } else {
	       std::cout << "INFO:: delete old entry for " << comp_id << std::endl;
	       // trash the old one then
	       dict_res_restraints[i].second.clear_dictionary_residue();
	    }
	 }
      }
   }

   if (! ifound) {
      dictionary_residue_restraints_t rest(comp_id, read_number);
      rest.atom_info.push_back(atom_info);
      std::pair<int, dictionary_residue_restraints_t> p(imol_enc, rest);
      dict_res_restraints.push_back(p);
   }
}


void
coot::protein_geometry::mon_lib_add_tree(std::string comp_id,
					 int imol_enc,
					 std::string atom_id,
					 std::string atom_back,
					 std::string atom_forward,
					 std::string connect_type) {

   coot::dict_chem_comp_tree_t ac(atom_id, atom_back, atom_forward, connect_type);
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
	 if (dict_res_restraints[i].first == imol_enc) {
	    dict_res_restraints[i].second.tree.push_back(ac);
	    break;
	 }
      }
   }
}

void
coot::protein_geometry::mon_lib_add_bond(std::string comp_id,
					 int imol_enc,
					 std::string atom_id_1,
					 std::string atom_id_2,
					 std::string type,
					 mmdb::realtype value_dist,
					 mmdb::realtype value_dist_esd,
					 mmdb::realtype value_dist_nuclear,
					 mmdb::realtype value_dist_nuclear_esd,
					 dict_bond_restraint_t::aromaticity_t arom_in,
                                         dict_bond_restraint_t::bond_length_type_t type_in) {

   if (false)
      std::cout << "adding bond for " << comp_id << " " << atom_id_1
		<< " " << atom_id_2 << " " << type << " " << value_dist
		<< " " << value_dist_esd << std::endl;

   // add a bond restraint to the list for comp_id.
   // The list container for comp_id is a dictionary_residue_restraints_t

   bool value_dist_nuclear_was_set = false;
   if (value_dist_nuclear_esd > 0.0)
      value_dist_nuclear_was_set = true;

   add_restraint(comp_id, imol_enc, dict_bond_restraint_t(atom_id_1,
							  atom_id_2,
							  type,
							  value_dist,
							  value_dist_esd,
							  value_dist_nuclear,
							  value_dist_nuclear_esd,
                                                          value_dist_nuclear_was_set,
							  arom_in,
                                                          type_in));
}

void
coot::protein_geometry::mon_lib_add_bond_no_target_geom(std::string comp_id,
							int imol_enc,
							std::string atom_id_1,
							std::string atom_id_2,
							std::string type,
							dict_bond_restraint_t::aromaticity_t arom_in) { 

   add_restraint(comp_id, imol_enc, dict_bond_restraint_t(atom_id_1, atom_id_2, type, arom_in));
}


void
coot::protein_geometry::mon_lib_add_angle(std::string comp_id,
					  int imol_enc,
					  std::string atom_id_1,
					  std::string atom_id_2,
					  std::string atom_id_3,
					  mmdb::realtype value_angle,
					  mmdb::realtype value_angle_esd) {

//    std::cout << "adding angle " << comp_id <<  " " << atom_id_1
// 	     << " " << atom_id_2 << " " << atom_id_3 << " "
// 	     << value_angle << std::endl;

   add_restraint(comp_id, imol_enc, dict_angle_restraint_t(atom_id_1,
						atom_id_2,
						atom_id_3,
						value_angle,
						value_angle_esd));
}

void
coot::protein_geometry::mon_lib_add_torsion(std::string comp_id,
					    int imol_enc,
					    std::string torsion_id,
					    std::string atom_id_1,
					    std::string atom_id_2,
					    std::string atom_id_3,
					    std::string atom_id_4,
					    mmdb::realtype value_angle,
					    mmdb::realtype value_angle_esd,
					    int period) {

   if (0)
      std::cout << "adding torsion " << comp_id <<  " " << atom_id_1
		<< " " << atom_id_2 << " " << atom_id_3 << " "
		<< atom_id_4 << " "
		<< "value_angle: " << value_angle
		<< ", value_angle_esd: " << value_angle_esd
		<< ", period: " << period << std::endl;
   
   add_restraint(comp_id, imol_enc, dict_torsion_restraint_t(torsion_id,
						   atom_id_1,
						   atom_id_2,
						   atom_id_3,
						   atom_id_4,
						   value_angle,
						   value_angle_esd,
						   period));

}



void
coot::protein_geometry::mon_lib_add_chiral(std::string comp_id,
					   int imol_enc,
					   std::string id,
					   std::string atom_id_centre,
					   std::string atom_id_1,
					   std::string atom_id_2,
					   std::string atom_id_3,
					   std::string volume_sign) {

    
    int volume_sign_int = 0;

    volume_sign_int = coot::protein_geometry::chiral_volume_string_to_chiral_sign(volume_sign);

    // std::cout << "DEBUG:: " << comp_id << " " << atom_id_centre << " " << volume_sign
    // << " " << volume_sign_int << std::endl;

//     std::cout << "adding chiral " << comp_id <<  " " << id
// 	      << " " << atom_id_centre << " " << volume_sign 
// 	      << " " << volume_sign.substr(0,3) << " " << volume_sign_int << endl; 
    
    // We only want to know about dictionary restraints that are
    // certainly one thing or the other, not "both".  The CB of VAL is
    // labelled as "both".  This is because it is the formal
    // definition of a chiral centre - the xGn atoms in a VAL are
    // equivalent so the Cb is not a chiral centre. Only if the atoms
    // are not equivalent can the atom be a chiral centre.  However,
    // in the equivalent atom case, we can use the "chiral" volume to
    // find residues with nomenclature errors (the CG1 and CG2 atom
    // names are swapped).
    // 
    if (volume_sign_int != 0)
       if (volume_sign_int != coot::dict_chiral_restraint_t::CHIRAL_VOLUME_RESTRAINT_VOLUME_SIGN_UNASSIGNED)
	  add_restraint(comp_id, imol_enc, dict_chiral_restraint_t(id, atom_id_centre,
								   atom_id_1,
								   atom_id_2,
								   atom_id_3,
								   volume_sign_int));
    

}

// Add a plane restraint atom
//
// Add a plane atom, we look through restraints trying to find a
// comp_id, and then a plane_id, if we find it, simply push back
// the atom name, if not, we create a restraint.
// 
void
coot::protein_geometry::mon_lib_add_plane(const std::string &comp_id,
					  int imol_enc,
					  const std::string &plane_id,
					  const std::string &atom_id,
					  const mmdb::realtype &dist_esd) {

   if (false)
      std::cout << "adding plane " << comp_id <<  " " << plane_id
		<< " " << atom_id << std::endl;

   bool ifound = false;

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) { 
      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
	 if (dict_res_restraints[i].first == imol_enc) {
	    for (unsigned int ip=0; ip<dict_res_restraints[i].second.plane_restraint.size(); ip++) {
	       if (dict_res_restraints[i].second.plane_restraint[ip].plane_id == plane_id) { 
		  ifound = true;
		  dict_res_restraints[i].second.plane_restraint[ip].push_back_atom(atom_id, dist_esd);
		  break;
	       }
	    }
	    if (! ifound) {
	       // we have comp_id, but no planes of that id.
	       coot::dict_plane_restraint_t res(plane_id, atom_id, dist_esd);
	       dict_res_restraints[i].second.plane_restraint.push_back(res);
	       ifound = true;
	    }
	 }
      }
   }

   // It was not there.  This should only happen when plane restraints
   // are declared in the dictionary *before* the bonds (or angles).
   // Which is never, thanks to Alexei.
   //
   // So, there was no comp_id found 
   if (! ifound) { 
      // add the plane to the newly created dictionary_residue_restraints_t
      dictionary_residue_restraints_t rest(comp_id, read_number);
      std::pair<int, dictionary_residue_restraints_t> p(imol_enc, rest);
      dict_res_restraints.push_back(p);
      coot::dict_plane_restraint_t res(plane_id, atom_id, dist_esd);
      dict_res_restraints[dict_res_restraints.size()-1].second.plane_restraint.push_back(res);
   }
}
