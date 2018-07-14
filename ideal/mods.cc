/* ideal/simple-restraint.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2008, 2009, 2010  by The University of Oxford
 * Copyright 2016 by Medical Research Council
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

#include <string.h> // for strcmp


// we don't want to compile anything if we don't have gsl
#ifdef HAVE_GSL


#include <fstream>
#include <algorithm> // for sort
#include <stdexcept>

#include "simple-restraint.hh"

coot::restraints_container_t::restraint_counts_t 
coot::restraints_container_t::apply_mods(int idr, mmdb::PPAtom res_selection,
					  int i_no_res_atoms,
					  mmdb::PResidue residue_p,
					 const coot::protein_geometry &geom) {

   coot::restraints_container_t::restraint_counts_t mod_counts;

   // does this residue have an OXT? (pre-cached).  If yes, add a mod_COO
   //
   if (residues_with_OXTs.size()) {
      if (std::find(residues_with_OXTs.begin(),
		    residues_with_OXTs.end(),
		    residue_p) != residues_with_OXTs.end()) {
	 apply_mod("COO", geom, idr, residue_p);
      }
   }

   return mod_counts;
}

void
coot::restraints_container_t::apply_mod(const std::string &mod_name,
					const coot::protein_geometry &geom,
					int idr,
					mmdb::PResidue residue_p) {

   // We crash here when linked with CCP4srs, geom.mods has been corrupted.
   //
   if (false) {
      std::map<std::string, coot::chem_mod>::const_iterator iit;
      for (iit=geom.mods.begin(); iit!=geom.mods.end(); iit++) 
	 std::cout << "  " << iit->first << std::endl;
   }

   std::map<std::string, coot::chem_mod>::const_iterator it = geom.mods.find(mod_name);

   if (it != geom.mods.end()) {

      for (unsigned int i=0; i<it->second.bond_mods.size(); i++) {
	 apply_mod_bond(it->second.bond_mods[i], residue_p);
      }
      for (unsigned int i=0; i<it->second.angle_mods.size(); i++) {
	 apply_mod_angle(it->second.angle_mods[i], residue_p);
      }
      for (unsigned int i=0; i<it->second.plane_mods.size(); i++) {
	 apply_mod_plane(it->second.plane_mods[i], residue_p);
      }
   } else {
      std::cout << "WARNING:: mod name \"" << mod_name << "\" not found in dictionary "
		<< std::endl;
   } 
}

void
coot::restraints_container_t::apply_mod_bond(const coot::chem_mod_bond &mod_bond,
					     mmdb::PResidue residue_p) {

   if (mod_bond.function == coot::CHEM_MOD_FUNCTION_ADD) {
      mod_bond_add(mod_bond, residue_p);
   }
   if (mod_bond.function == coot::CHEM_MOD_FUNCTION_CHANGE) {
      mod_bond_change(mod_bond, residue_p);
   }
   if (mod_bond.function == coot::CHEM_MOD_FUNCTION_DELETE) {
      mod_bond_delete(mod_bond, residue_p);
   }
}

void
coot::restraints_container_t::mod_bond_add(const coot::chem_mod_bond &mod_bond,
					   mmdb::PResidue residue_p) {

   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   
   int index_1 = -1, index_2 = -1;
   for (int iat_1=0; iat_1<n_residue_atoms; iat_1++) {
      std::string pdb_atom_name_1(residue_atoms[iat_1]->name);
      // std::cout << "comparing :" << pdb_atom_name_1 << ": with :" << mod_bond.atom_id_1
      // << ":" << std::endl;
      if (pdb_atom_name_1 == mod_bond.atom_id_1) {
	 for (int iat_2=0; iat_2<n_residue_atoms; iat_2++) {
	    std::string pdb_atom_name_2(residue_atoms[iat_2]->name);
	    if (pdb_atom_name_2 == mod_bond.atom_id_2) {
	       // check that they have the same alt conf
	       std::string alt_1(residue_atoms[iat_1]->altLoc);
	       std::string alt_2(residue_atoms[iat_2]->altLoc);
	       if (alt_1 == "" || alt_2 == "" || alt_1 == alt_2) {
		  residue_atoms[iat_1]->GetUDData(udd_atom_index_handle, index_1);
		  residue_atoms[iat_2]->GetUDData(udd_atom_index_handle, index_2);
		  bonded_atom_indices[index_1].push_back(index_2);
		  bonded_atom_indices[index_2].push_back(index_1);
		  std::vector<bool> fixed_flags = make_fixed_flags(index_1, index_2);

		  add(BOND_RESTRAINT, index_1, index_2,
		      fixed_flags,
		      mod_bond.new_value_dist,
		      mod_bond.new_value_dist_esd,
		      1.2);  // junk value
	       }
	    }
	 }
      }
   }
}

void
coot::restraints_container_t::mod_bond_change(const coot::chem_mod_bond &mod_bond,
					      mmdb::PResidue residue_p) {

   for (unsigned int i=0; i<restraints_vec.size(); i++) {
      if (restraints_vec[i].restraint_type == coot::BOND_RESTRAINT) {
	 const coot::simple_restraint &rest = restraints_vec[i];
	 if (atom[restraints_vec[i].atom_index_1]->residue == residue_p) {
	    if (atom[restraints_vec[i].atom_index_2]->residue == residue_p) {
	       std::string name_1 = atom[rest.atom_index_1]->name;
	       std::string name_2 = atom[rest.atom_index_2]->name;
	       if (name_1 == mod_bond.atom_id_1) {
		  if (name_2 == mod_bond.atom_id_2) {
		     restraints_vec[i].target_value = mod_bond.new_value_dist;
		     restraints_vec[i].sigma = mod_bond.new_value_dist_esd;

		     if (0) 
			std::cout << "DEBUG:: mod_bond_change() changed bond "
				  << coot::atom_spec_t(atom[restraints_vec[i].atom_index_1])
				  << " to " 
				  << coot::atom_spec_t(atom[restraints_vec[i].atom_index_2])
				  << " dist " <<  mod_bond.new_value_dist
				  << " esd " <<  mod_bond.new_value_dist_esd
				  << std::endl;
		  }
	       }
	    }
	 }
      }
   }
}

void
coot::restraints_container_t::mod_bond_delete(const coot::chem_mod_bond &mod_bond,
					      mmdb::PResidue residue_p) {


   std::vector<coot::simple_restraint>::iterator it;
   
   for (it=restraints_vec.begin(); it!=restraints_vec.end(); it++) { 
      if (it->restraint_type == coot::BOND_RESTRAINT) {
	 if (atom[it->atom_index_1]->residue == residue_p) {
	    if (atom[it->atom_index_2]->residue == residue_p) {
	       std::string name_1 = atom[it->atom_index_1]->name;
	       std::string name_2 = atom[it->atom_index_2]->name;
	       if (name_1 == mod_bond.atom_id_1) {
		  if (name_2 == mod_bond.atom_id_2) {
		     if (0) 
			std::cout << "DEBUG:: mod_bond_delete() delete bond "
				  << coot::atom_spec_t(atom[it->atom_index_1])
				  << " to " 
				  << coot::atom_spec_t(atom[it->atom_index_2])
				  << std::endl;
		     restraints_vec.erase(it);
		  }
	       }
	    }
	 }
      }
   }

}

void
coot::restraints_container_t::apply_mod_angle(const coot::chem_mod_angle &mod_angle,
					     mmdb::PResidue residue_p) {

   if (mod_angle.function == coot::CHEM_MOD_FUNCTION_ADD) {
      mod_angle_add(mod_angle, residue_p);
   }
   if (mod_angle.function == coot::CHEM_MOD_FUNCTION_CHANGE) {
      mod_angle_change(mod_angle, residue_p);
   }
   if (mod_angle.function == coot::CHEM_MOD_FUNCTION_DELETE) {
      mod_angle_delete(mod_angle, residue_p);
   }
}

void
coot::restraints_container_t::mod_angle_add(const coot::chem_mod_angle &mod_angle,
					    mmdb::PResidue residue_p) {

   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   
   int index_1 = -1, index_2 = -1, index_3 = -1;
   for (int iat_1=0; iat_1<n_residue_atoms; iat_1++) {
      std::string pdb_atom_name_1(residue_atoms[iat_1]->name);
      if (pdb_atom_name_1 == mod_angle.atom_id_1) {
	 for (int iat_2=0; iat_2<n_residue_atoms; iat_2++) {
	    std::string pdb_atom_name_2(residue_atoms[iat_2]->name);
	    if (pdb_atom_name_2 == mod_angle.atom_id_2) {
	       for (int iat_3=0; iat_3<n_residue_atoms; iat_3++) {
		  std::string pdb_atom_name_3(residue_atoms[iat_3]->name);
		  if (pdb_atom_name_3 == mod_angle.atom_id_3) {
		     
		     // check that they have the same alt conf
		     std::string alt_1(residue_atoms[iat_1]->altLoc);
		     std::string alt_2(residue_atoms[iat_2]->altLoc);
		     std::string alt_3(residue_atoms[iat_3]->altLoc);
		     if (((alt_1 == alt_2) && (alt_1 == alt_3)) ||
			 ((alt_1 == ""   ) && (alt_2 == alt_3)) ||
			 ((alt_2 == ""   ) && (alt_1 == alt_3)) ||
			 ((alt_3 == ""   ) && (alt_1 == alt_2)))
			{
			   
			   residue_atoms[iat_1]->GetUDData(udd_atom_index_handle, index_1);
			   residue_atoms[iat_2]->GetUDData(udd_atom_index_handle, index_2);
			   residue_atoms[iat_3]->GetUDData(udd_atom_index_handle, index_3);
			   std::vector<bool> fixed_flags =
			      make_fixed_flags(index_1, index_2, index_3);
			   bool is_single_Hydrogen_atom_angle_restraint = false;
			   unsigned int nH = 0;
			   if (is_hydrogen(residue_atoms[iat_1])) nH++;
			   if (is_hydrogen(residue_atoms[iat_3])) nH++;
			   if (nH == 1) is_single_Hydrogen_atom_angle_restraint = true;

			   add(ANGLE_RESTRAINT, index_1, index_2, index_3,
			       fixed_flags,
			       mod_angle.new_value_angle,
			       mod_angle.new_value_angle_esd,
			       is_single_Hydrogen_atom_angle_restraint);  // junk value
			}
		  }
	       }
	    }
	 }
      }
   }
}


void
coot::restraints_container_t::mod_angle_change(const coot::chem_mod_angle &mod_angle,
					       mmdb::PResidue residue_p) {

   for (unsigned int i=0; i<restraints_vec.size(); i++) {
      if (restraints_vec[i].restraint_type == coot::ANGLE_RESTRAINT) {
	 const coot::simple_restraint &rest = restraints_vec[i];
	 if (atom[restraints_vec[i].atom_index_1]->residue == residue_p) {
	    if (atom[restraints_vec[i].atom_index_2]->residue == residue_p) {
	       std::string name_1 = atom[rest.atom_index_1]->name;
	       std::string name_2 = atom[rest.atom_index_2]->name;
	       std::string name_3 = atom[rest.atom_index_3]->name;
	       if (name_1 == mod_angle.atom_id_1) {
		  if (name_2 == mod_angle.atom_id_2) {
		     if (name_3 == mod_angle.atom_id_3) {
			restraints_vec[i].target_value = mod_angle.new_value_angle;
			restraints_vec[i].sigma = mod_angle.new_value_angle_esd;
			if (0) 
			   std::cout << "DEBUG:: mod_angle_change() changed angle "
				     << coot::atom_spec_t(atom[restraints_vec[i].atom_index_1])
				     << " to " 
				     << coot::atom_spec_t(atom[restraints_vec[i].atom_index_2])
				     << " to " 
				     << coot::atom_spec_t(atom[restraints_vec[i].atom_index_3])
				     << " angle " <<  mod_angle.new_value_angle
				     << " esd " <<  mod_angle.new_value_angle_esd
				     << std::endl;
		     }
		  }
	       }
	    }
	 }
      }
   }
}



void
coot::restraints_container_t::mod_angle_delete(const coot::chem_mod_angle &mod_angle,
					       mmdb::PResidue residue_p) {


   std::vector<coot::simple_restraint>::iterator it;
   
   for (it=restraints_vec.begin(); it!=restraints_vec.end(); it++) { 
      if (it->restraint_type == coot::ANGLE_RESTRAINT) {
	 if (atom[it->atom_index_1]->residue == residue_p) {
	    if (atom[it->atom_index_2]->residue == residue_p) {
	       std::string name_1 = atom[it->atom_index_1]->name;
	       std::string name_2 = atom[it->atom_index_2]->name;
	       std::string name_3 = atom[it->atom_index_3]->name;
	       if (name_1 == mod_angle.atom_id_1) {
		  if (name_2 == mod_angle.atom_id_2) {
		     if (name_2 == mod_angle.atom_id_3) {
			if (0) 
			   std::cout << "DEBUG:: mod_angle_delete() delete angle "
				     << coot::atom_spec_t(atom[it->atom_index_1])
				     << " to " 
				  << coot::atom_spec_t(atom[it->atom_index_2])
				     << " to " 
				     << coot::atom_spec_t(atom[it->atom_index_3])
				     << std::endl;
			restraints_vec.erase(it);
		     }
		  }
	       }
	    }
	 }
      }
   }
}

void
coot::restraints_container_t::apply_mod_plane(const coot::chem_mod_plane &mod_plane,
					      mmdb::PResidue residue_p) {

   if (mod_plane.function == coot::CHEM_MOD_FUNCTION_ADD) {
      mod_plane_add(mod_plane, residue_p);
   }
   if (mod_plane.function == coot::CHEM_MOD_FUNCTION_DELETE) {
      mod_plane_delete(mod_plane, residue_p);
   }
}


void
coot::restraints_container_t::mod_plane_add(const coot::chem_mod_plane &mod_plane,
					    mmdb::PResidue residue_p) {
   
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   
   std::map<std::string, std::vector <int> > pos; // we worry about alt confs.

   for (unsigned int i=0; i<mod_plane.atom_id_esd.size(); i++) {
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 std::string atom_name(residue_atoms[iat]->name);
	 if (atom_name == mod_plane.atom_id_esd[i].first) {
	    int atom_index;
	    residue_atoms[iat]->GetUDData(udd_atom_index_handle, atom_index);
	    std::string altconf = residue_atoms[iat]->altLoc;
	    pos[altconf].push_back(atom_index);
	 }
      }
   }

   // iterate through all the alt confs (almost certainly only one)
   std::map<std::string, std::vector <int> >::const_iterator it;
   for (it=pos.begin(); it!=pos.end(); it++) {
      const std::vector<int> &position_indices = it->second;

      if (position_indices.size() > 3) {
	 double esd = 0.02;

	 std::vector<std::pair<int, double> > position_sigma_indices;
	 for (unsigned int ii=0; ii<position_indices.size(); ii++)
	    position_sigma_indices.push_back(std::pair<int, double> (position_indices[ii], 0.02));
	 
	 std::vector<bool> fixed_flags = make_fixed_flags(position_indices);
	 add_plane(position_sigma_indices, fixed_flags);
	 if (false) {
	    std::cout << "DEBUG:: mod_plane_add() adding plane\n";
	    for (unsigned int i=0; i<position_indices.size(); i++)
	       std::cout << "   " << coot::atom_spec_t(atom[position_indices[i]]) << "\n";
	 }
      }
   }
}

void
coot::restraints_container_t::mod_plane_delete(const coot::chem_mod_plane &mod_plane,
					       mmdb::PResidue residue_p) {

   std::vector<coot::simple_restraint>::iterator it;
   
   for (it=restraints_vec.begin(); it!=restraints_vec.end(); it++) { 
      if (it->restraint_type == coot::PLANE_RESTRAINT) {
	 bool in_same_residue = 1;
	 unsigned int n_found = 0;
	 // do the atoms of the mod_plane match the atoms of the restraint?
	 for (unsigned int iat=0; iat<it->plane_atom_index.size(); iat++) { 
	    for (unsigned int iat_mod=0; iat_mod<mod_plane.atom_id_esd.size(); iat_mod++) {
	       std::string atom_name = atom[it->plane_atom_index[iat].first]->name;
	       if (atom_name == mod_plane.atom_id_esd[iat_mod].first) {
		  if (atom[it->plane_atom_index[iat].first]->GetResidue() == residue_p) {
		     n_found++;
		     break;
		  }
	       }
	    }
	 }
	 if (n_found == it->plane_atom_index.size()) {

	    if (0) { 
	       std::cout << "DEBUG:: mod_plane_delete() delete plane ";
	       for (unsigned int iat=0; iat<it->plane_atom_index.size(); iat++)
		  std::cout << "   " << coot::atom_spec_t(atom[it->plane_atom_index[iat].first])
			    << "\n";
	    }
	    restraints_vec.erase(it);
	 } 
      }
   }
}

#endif // HAVE_GSL
