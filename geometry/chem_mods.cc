/* geometry/protein-geometry.cc
 * 
 * Copyright 2011 The University of Oxford
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <string.h>
#include <iostream>
#include <map>
#include <algorithm>  // needed for sort? Yes.
#include <stdexcept>  // Thow execption.

#include "protein-geometry.hh"
#include "utils/coot-utils.hh"

#include <sys/types.h> // for stating
#include <sys/stat.h>

#include "clipper/core/clipper_util.h"

#include "lbg-graph.hh"

#include "compat/coot-sysdep.h"

#ifdef COOT_ENABLE_WINAPI_SUSPENSION
# undef AddAtom
#endif // COOT_ENABLE_WINAPI_SUSPENSION

std::ostream&
coot::operator<<(std::ostream &s, const coot::chem_mod_atom &a) {

   s << "[chem_mod_atom "
     << a.function << " atom_id :"
     << a.atom_id << ": new_atom_id :"
     << a.new_atom_id << ": new_type_symbol :"
     << a.new_type_symbol << ": new_type_energy :"
     << a.new_type_energy << ": new_partial_charge "
     << a.new_partial_charge << "]";
   return s;
} 

std::ostream&
coot::operator<<(std::ostream &s, const coot::chem_mod_tree &a) {

   s << "[chem_mod_tree "
     << a.function << " :"
     << a.atom_id << ": :"
     << a.atom_back << ": :"
     << a.back_type << ": :"
     << a.atom_forward << ": "
     << a.connect_type << ":]";
   return s;
} 

std::ostream&
coot::operator<<(std::ostream &s, const coot::chem_mod_bond &a) {

   s << "[chem_mod_bond "
     << a.function << " :"
     << a.atom_id_1 << ": :"
     << a.atom_id_2 << ": :"
     << a.new_type << ": "
     << a.new_value_dist << " "
     << a.new_value_dist_esd << "]";
   return s;
} 

std::ostream&
coot::operator<<(std::ostream &s, const coot::chem_mod_angle &a) {

   s << "[chem_mod_angle "
     << a.function << " "
     << a.atom_id_1 << " "
     << a.atom_id_2 << " "
     << a.atom_id_3 << " "
     << a.new_value_angle << " "
     << a.new_value_angle_esd << "]";
   return s;
} 

std::ostream&
coot::operator<<(std::ostream &s, const coot::chem_mod_tor &a) {

   s << "[chem_mod_tor "
     << a.function << " "
     << a.atom_id_1 << " "
     << a.atom_id_2 << " "
     << a.atom_id_3 << " "
     << a.atom_id_4 << " "
     << a.new_value_angle << " "
     << a.new_value_angle_esd << " "
     << a.new_period << "]";
   return s;
}

std::ostream&
coot::operator<<(std::ostream &s, const coot::chem_mod_plane &a) {

   s << "[chem_mod_plane function="
     << a.function << " "
     << a.plane_id << " ";
   s << " n_atoms=" << a.atom_id_esd.size();
   for (unsigned int i=0; i<a.atom_id_esd.size(); i++)
      s << "  " << a.atom_id_esd[i].first << " "
	<< a.atom_id_esd[i].second;
   s << "]";
   return s;
}

std::ostream&
coot::operator<<(std::ostream &s, const coot::chem_mod_chir &a) {

   s << "[chem_mod_chir "
     << a.function << " "
     << a.atom_id_centre << " "
     << a.atom_id_1 << " "
     << a.atom_id_2 << " "
     << a.atom_id_3 << " "
     << a.new_volume_sign << "]";
   return s;
}




int
coot::protein_geometry::add_mod(mmdb::mmcif::PData data) {

   int status = 0;

   int r = 0; 
   for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 
      
      mmdb::mmcif::PCategory cat = data->GetCategory(icat);
      std::string cat_name(cat->GetCategoryName());
      mmdb::mmcif::PLoop mmCIFLoop = data->GetLoop(cat_name.c_str()) ;
            
      if (mmCIFLoop == NULL) { 
	 std::cout << "null loop" << std::endl; 
      } else {
	 int n_chiral = 0;
	 if (cat_name == "_chem_mod_atom")
	    add_chem_mod_atom(mmCIFLoop);
	 if (cat_name == "_chem_mod_bond")
	    add_chem_mod_bond(mmCIFLoop);
	 if (cat_name == "_chem_mod_tree")
	    add_chem_mod_tree(mmCIFLoop);
	 if (cat_name == "_chem_mod_angle")
	    add_chem_mod_angle(mmCIFLoop);
	 if (cat_name == "_chem_mod_tor")
	    add_chem_mod_tor(mmCIFLoop);
	 if (cat_name == "_chem_mod_chir")
	    add_chem_mod_chir(mmCIFLoop);
	 if (cat_name == "_chem_mod_plane_atom")
	    add_chem_mod_plane(mmCIFLoop);
      }
   }
   return status;
}

void
coot::protein_geometry::add_chem_mod_atom( mmdb::mmcif::PLoop mmCIFLoop) {

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

      int ierr_tot = 0;
      int ierr;
      
      std::string mod_id;
      std::string function;
      std::string atom_id;
      std::string new_atom_id;
      std::string new_type_symbol;
      std::string new_type_energy;
      mmdb::realtype new_partial_charge;
      
      char *s;

      s = mmCIFLoop->GetString("mod_id", j, ierr);
      ierr_tot += ierr;
      if (ierr) std::cout << "   oops getting mod_id " << std::endl;
      if (s) mod_id = s;

      s = mmCIFLoop->GetString("function", j, ierr);
      if (ierr) std::cout << "   oops getting function " << std::endl;
      ierr_tot += ierr;
      if (s) function = s;

      s = mmCIFLoop->GetString("atom_id", j, ierr);
      if (ierr) std::cout << "   oops getting atom_id " << std::endl;
      ierr_tot += ierr;
      if (s) atom_id = s;

      s = mmCIFLoop->GetString("new_atom_id", j, ierr);
      if (ierr) std::cout << "   oops getting new_atom_id " << std::endl;
      ierr_tot += ierr;
      if (s) new_atom_id = s;

      s = mmCIFLoop->GetString("new_type_symbol", j, ierr);
      if (ierr) std::cout << "   oops getting new_type_symbol " << std::endl;
      ierr_tot += ierr;
      if (s) new_type_symbol = s;

      s = mmCIFLoop->GetString("new_type_energy", j, ierr);
      if (ierr) std::cout << "   oops getting new_type_energy " << std::endl;
      ierr_tot += ierr;
      if (s) new_type_energy = s;

      ierr = mmCIFLoop->GetReal(new_partial_charge, "new_partial_charge", j);
      // OK, with the old dictionary can fail here (the format of the
      // dictionary is wrong - you can't have default value "." for
      // numbers in cif.
      if (ierr) {
	 new_partial_charge = 0.0; // dummy value, for atoms that
				   // don't matter - e.g. function is
				   // "delete"
      } 
      
      if (ierr_tot == 0) {
	 coot::chem_mod_atom cma(function,
				 atom_id, new_atom_id,
				 new_type_symbol, new_type_energy,
				 new_partial_charge);
	 mods[mod_id].add_mod_atom(cma);
      } else {
	 std::cout << "oops in add_chem_mod_atom ierr_tot = "
		   << ierr_tot << std::endl;
	 std::cout << "   mod_id: \"" << mod_id
		   << "\"    function: \"" << function << "\" atom_id: \""
		   << atom_id_mmdb_expand(atom_id) << "\" new_atom_id: \"" 
		   << atom_id_mmdb_expand(new_atom_id) << "\" new_type_symbol: \"" 
		   << new_type_symbol << "\" new_type_energy: \"" 
		   << new_type_energy << "\" new_partial_charge: \"" 
		   << new_partial_charge << "\"" << std::endl;
      } 
   } 
      
} 

void
coot::protein_geometry::add_chem_mod_bond( mmdb::mmcif::PLoop mmCIFLoop) {

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

      int ierr_tot = 0;
      int ierr;
      
      std::string mod_id;
      std::string function;
      std::string atom_id_1;
      std::string atom_id_2;
      std::string new_type;
      mmdb::realtype new_value_dist;
      mmdb::realtype new_value_dist_esd;
      
      char *s;

      s = mmCIFLoop->GetString("mod_id", j, ierr);
      ierr_tot += ierr;
      if (s) mod_id = s;

      s = mmCIFLoop->GetString("function", j, ierr);
      ierr_tot += ierr;
      if (s) function = s;

      s = mmCIFLoop->GetString("atom_id_1", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id_1 = s;

      s = mmCIFLoop->GetString("atom_id_2", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id_2 = s;

      s = mmCIFLoop->GetString("new_type", j, ierr);
      ierr_tot += ierr;
      if (s) new_type = s;

      ierr = mmCIFLoop->GetReal(new_value_dist, "new_value_dist", j);
      ierr_tot += ierr; 
      
      ierr = mmCIFLoop->GetReal(new_value_dist_esd, "new_value_dist_esd", j);
      ierr_tot += ierr;

      if (ierr_tot == 0) {
	 coot::chem_mod_bond cmb(function,
				 atom_id_mmdb_expand(atom_id_1),
				 atom_id_mmdb_expand(atom_id_2), new_type,
				 new_value_dist, new_value_dist_esd);
	 mods[mod_id].add_mod_bond(cmb);
      } 
   } 
}

void
coot::protein_geometry::add_chem_mod_tree( mmdb::mmcif::PLoop mmCIFLoop) {

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

      int ierr_tot = 0;
      int ierr;
      
      std::string mod_id;
      std::string function;
      std::string atom_id;
      std::string atom_back;
      std::string back_type;
      std::string atom_forward;
      std::string connect_type;
      
      char *s;

      s = mmCIFLoop->GetString("mod_id", j, ierr);
      ierr_tot += ierr;
      if (s) mod_id = s;

      s = mmCIFLoop->GetString("function", j, ierr);
      ierr_tot += ierr;
      if (s) function = s;

      s = mmCIFLoop->GetString("atom_id", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id = s;

      s = mmCIFLoop->GetString("atom_back", j, ierr);
      ierr_tot += ierr;
      if (s) atom_back = s;

      s = mmCIFLoop->GetString("back_type", j, ierr);
      ierr_tot += ierr;
      if (s) back_type = s;

      s = mmCIFLoop->GetString("atom_forward", j, ierr);
      ierr_tot += ierr;
      if (s) atom_forward = s;

      s = mmCIFLoop->GetString("connect_type", j, ierr);
      ierr_tot += ierr;
      if (s) connect_type = s;

      if (ierr_tot == 0) {
	 coot::chem_mod_tree cmt(function, atom_id_mmdb_expand(atom_id),
				 atom_id_mmdb_expand(atom_back),
				 back_type, atom_id_mmdb_expand(atom_forward),
				 connect_type);
	 mods[mod_id].add_mod_tree(cmt);
      } 
   } 
}

void
coot::protein_geometry::add_chem_mod_angle(mmdb::mmcif::PLoop mmCIFLoop) {

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

      int ierr_tot = 0;
      int ierr;
      
      std::string mod_id;
      std::string function;
      std::string atom_id_1;
      std::string atom_id_2;
      std::string atom_id_3;
      mmdb::realtype new_value_angle;
      mmdb::realtype new_value_angle_esd;
      
      char *s;

      s = mmCIFLoop->GetString("mod_id", j, ierr);
      ierr_tot += ierr;
      if (s) mod_id = s;

      s = mmCIFLoop->GetString("function", j, ierr);
      ierr_tot += ierr;
      if (s) function = s;

      s = mmCIFLoop->GetString("atom_id_1", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id_1 = s;

      s = mmCIFLoop->GetString("atom_id_2", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id_2 = s;

      s = mmCIFLoop->GetString("atom_id_3", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id_3 = s;

      ierr = mmCIFLoop->GetReal(new_value_angle, "new_value_angle", j);
      ierr_tot += ierr; 
      
      ierr = mmCIFLoop->GetReal(new_value_angle_esd, "new_value_angle_esd", j);
      ierr_tot += ierr;

      if (ierr_tot == 0) {
	 coot::chem_mod_angle cma(function,
				  atom_id_mmdb_expand(atom_id_1),
				  atom_id_mmdb_expand(atom_id_2),
				  atom_id_mmdb_expand(atom_id_3),
				  new_value_angle,
				  new_value_angle_esd);
	 mods[mod_id].add_mod_angle(cma);
      } 
   } 
}

void
coot::protein_geometry::add_chem_mod_tor(mmdb::mmcif::PLoop mmCIFLoop) {

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

      int ierr_tot = 0;
      int ierr;
      
      std::string mod_id;
      std::string function;
      std::string atom_id_1;
      std::string atom_id_2;
      std::string atom_id_3;
      std::string atom_id_4;
      mmdb::realtype new_value_angle;
      mmdb::realtype new_value_angle_esd;
      int new_period;
      
      char *s;

      s = mmCIFLoop->GetString("mod_id", j, ierr);
      ierr_tot += ierr;
      if (s) mod_id = s;

      s = mmCIFLoop->GetString("function", j, ierr);
      ierr_tot += ierr;
      if (s) function = s;

      s = mmCIFLoop->GetString("atom_id_1", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id_1 = s;

      s = mmCIFLoop->GetString("atom_id_2", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id_2 = s;

      s = mmCIFLoop->GetString("atom_id_3", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id_3 = s;

      s = mmCIFLoop->GetString("atom_id_4", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id_4 = s;

      ierr = mmCIFLoop->GetReal(new_value_angle, "new_value_angle", j);
      ierr_tot += ierr; 
      
      ierr = mmCIFLoop->GetReal(new_value_angle_esd, "new_value_angle_esd", j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetInteger(new_period, "new_period", j);
      ierr_tot += ierr;

      if (ierr_tot == 0) {
	 coot::chem_mod_tor cmt(function,
				atom_id_mmdb_expand(atom_id_1),
				atom_id_mmdb_expand(atom_id_2), 
				atom_id_mmdb_expand(atom_id_3),  
				atom_id_mmdb_expand(atom_id_4), 
				new_value_angle,
				new_value_angle_esd,
				new_period);
	 mods[mod_id].add_mod_tor(cmt);
      } 
   } 
}

void
coot::protein_geometry::add_chem_mod_chir(mmdb::mmcif::PLoop mmCIFLoop) {

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

      int ierr_tot = 0;
      int ierr;
      
      std::string mod_id;
      std::string function;
      std::string atom_id_centre;
      std::string atom_id_1;
      std::string atom_id_2;
      std::string atom_id_3;
      std::string new_volume_sign;
      
      char *s;

      s = mmCIFLoop->GetString("mod_id", j, ierr);
      ierr_tot += ierr;
      if (s) mod_id = s;

      s = mmCIFLoop->GetString("function", j, ierr);
      ierr_tot += ierr;
      if (s) function = s;

      s = mmCIFLoop->GetString("atom_id_centre", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id_centre = s;

      s = mmCIFLoop->GetString("atom_id_1", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id_1 = s;

      s = mmCIFLoop->GetString("atom_id_2", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id_2 = s;

      s = mmCIFLoop->GetString("atom_id_3", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id_3 = s;

      s = mmCIFLoop->GetString("new_volume_sign", j, ierr);
      ierr_tot += ierr;
      if (s) new_volume_sign = s;

      int volume_sign_int =
	 coot::protein_geometry::chiral_volume_string_to_chiral_sign(new_volume_sign);
      

      if (ierr_tot == 0) {
	 coot::chem_mod_chir cmc(function,
				 atom_id_mmdb_expand(atom_id_centre),
				 atom_id_mmdb_expand(atom_id_1),
				 atom_id_mmdb_expand(atom_id_2),
				 atom_id_mmdb_expand(atom_id_3), 
				 volume_sign_int);
	 mods[mod_id].add_mod_chir(cmc);
      } else {
	 std::cout << "oops in add_chem_mod_chir ierr_tot is "
		   << ierr_tot << std::endl;
      } 
   }
}


void
coot::protein_geometry::add_chem_mod_plane(mmdb::mmcif::PLoop mmCIFLoop) {

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

      int ierr_tot = 0;
      int ierr;
      
      std::string mod_id;
      std::string plane_id;
      std::string function;
      std::string atom_id;
      mmdb::realtype new_dist_esd;
      
      char *s;

      s = mmCIFLoop->GetString("mod_id", j, ierr);
      ierr_tot += ierr;
      if (s) mod_id = s;

      s = mmCIFLoop->GetString("plane_id", j, ierr);
      ierr_tot += ierr;
      if (s) plane_id = s;

      s = mmCIFLoop->GetString("function", j, ierr);
      ierr_tot += ierr;
      if (s) function = s;

      s = mmCIFLoop->GetString("atom_id", j, ierr);
      ierr_tot += ierr;
      if (s) atom_id = s;

      ierr = mmCIFLoop->GetReal(new_dist_esd, "new_dist_esd", j);
      if (ierr)
	 new_dist_esd = 0.0; // some ignored dummy value
      ierr_tot += ierr;

      if (ierr_tot == 0 || function == "delete") {
	 // coot::chem_mod_plane plane(plane_id, function);
	 // mods[mod_id][plane].add_atom(atom_id_mmdb_expand(atom_id), new_dist_esd);
	 std::string atom_name = atom_id_mmdb_expand(atom_id);
	 mods[mod_id].add_plane_atom(plane_id, function, atom_name, new_dist_esd);
      } else {
	 std::cout << "oops in add_chem_mod_plane ierr_tot is "
		   << ierr_tot << std::endl;
      }
   }
}

// can throw a std::runtime_error
// 
std::pair<coot::chem_mod, coot::chem_mod>
coot::protein_geometry::get_chem_mods_for_link(const std::string &link_id) const {

   bool found = false; 

//    for (unsigned int ilink=0; ilink<chem_link_vec.size(); ilink++) {
//       const chem_link &chem_link = chem_link_vec[ilink];
//       if (chem_link.Id() == link_id) {
// 	 std::pair<std::string, std::string> mod_names = chem_link.chem_mod_names();
// 	 std::map<std::string, chem_mod>::const_iterator it_1 = mods.find(mod_names.first);
// 	 std::map<std::string, chem_mod>::const_iterator it_2 = mods.find(mod_names.second);
// 	 if (it_1 != mods.end()) {
// 	    if (it_2 != mods.end()) {
// 	       // we found both mods.
// 	       // std::cout << "we found both mods! " << std::endl;
// 	       return std::pair<chem_mod, chem_mod>(it_1->second, it_2->second);
// 	    } else {
// 	       std::cout << "DEBUG:: oops no " << mod_names.second
// 			 << " in mods" << std::endl;
// 	    }
// 	 } else {
// 	    std::cout << "DEBUG:: oops no " << mod_names.first << " in mods" << std::endl;
// 	 }
//    }

   std::map<unsigned int, std::vector<chem_link> >::const_iterator it;
   for (it=chem_link_map.begin(); it!=chem_link_map.end(); it++) {
      const std::vector<chem_link> &v = it->second;
      std::vector<chem_link>::const_iterator itv;
      for (itv=v.begin(); itv!=v.end(); itv++) {
	 const chem_link &cl = *itv;
	 if (cl.Id() == link_id) {
	    std::pair<std::string, std::string> mod_names = cl.chem_mod_names();
	    std::map<std::string, chem_mod>::const_iterator it_1 = mods.find(mod_names.first);
	    std::map<std::string, chem_mod>::const_iterator it_2 = mods.find(mod_names.second);
	    if (it_1 != mods.end()) {
	       if (it_2 != mods.end()) {
		  // we found both mods.
		  // std::cout << "we found both mods! " << std::endl;
		  return std::pair<chem_mod, chem_mod>(it_1->second, it_2->second);
	       } else {
		  std::cout << "DEBUG:: oops no " << mod_names.second
			    << " in mods" << std::endl;
	       }
	    } else {
	       std::cout << "DEBUG:: oops no " << mod_names.first << " in mods" << std::endl;
	    }
	 }
      }
   }

   throw std::runtime_error("No link found");
   
//    if (link.link_id == "") {
//       throw std::runtime_error("No link found");
//    } else {
      
//       if (it != mods.end())
// 	 return it->second;
//       else
// 	 throw std::runtime_error("No chem_mod found");
//    }

}



void
coot::protein_geometry::debug_mods() const {

   std::map<std::string, chem_mod>::const_iterator it;

   for (it=mods.begin(); it!=mods.end(); it++) {

      std::cout << "----- mod: " << it->first
		<< " --------------" << std::endl;
      std::cout << "::: " << it->second.atom_mods.size()
		<< " atom mods" << std::endl;
      for (unsigned int i=0; i<it->second.atom_mods.size(); i++) {
	 std::cout << "   " << it->second.atom_mods[i] << std::endl;
      } 
      std::cout << "::: " << it->second.tree_mods.size()
		<< " tree mods" << std::endl;
      for (unsigned int i=0; i<it->second.tree_mods.size(); i++) {
	 std::cout << "   " << it->second.tree_mods[i] << std::endl;
      } 
      std::cout << "::: " << it->second.bond_mods.size()
		<< " bond mods" << std::endl;
      for (unsigned int i=0; i<it->second.bond_mods.size(); i++) {
	 std::cout << "   " << it->second.bond_mods[i] << std::endl;
      } 
      std::cout << "::: " << it->second.bond_mods.size()
		<< " angle mods" << std::endl;
      for (unsigned int i=0; i<it->second.angle_mods.size(); i++) {
	 std::cout << "   " << it->second.angle_mods[i] << std::endl;
      } 
      std::cout << "::: " << it->second.tor_mods.size()
		<< " tor mods" << std::endl;
      for (unsigned int i=0; i<it->second.tor_mods.size(); i++) {
	 std::cout << "   " << it->second.tor_mods[i] << std::endl;
      } 
      std::cout << "::: " << it->second.plane_mods.size()
		<< " plane mods" << std::endl;
      for (unsigned int i=0; i<it->second.plane_mods.size(); i++) {
	 std::cout << "   " << it->second.plane_mods[i] << std::endl;
      } 
      std::cout << "::: " << it->second.chir_mods.size()
		<< " chir mods" << std::endl;
      for (unsigned int i=0; i<it->second.chir_mods.size(); i++) {
	 std::cout << "   " << it->second.chir_mods[i] << std::endl;
      } 
   }
}
