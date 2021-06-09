/* geometry/protein-geometry.cc
 * 
 * Copyright 2003, 2004, 2005, 2006 The University of York
 * Author: Paul Emsley
 * Copyright 2007, 2008, 2009, 2010, 2011, 2012 The University of Oxford
 * Copyright 2014, 2015, 2016 by Medical Research Council
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

#include <string.h>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>  // needed for sort? Yes.
#include <stdexcept>  // Thow execption.

#include "mini-mol/atom-quads.hh"
#include "geometry/protein-geometry.hh"
#include "utils/coot-utils.hh"

// #include <sys/types.h> // for stating
// #include <sys/stat.h>

#if !defined _MSC_VER
#include <unistd.h>
#else
#define DATADIR "C:/coot/share"
#define PKGDATADIR DATADIR
// stop using these, use win-compat functions
// #define S_ISDIR(m)  (((m) & S_IFMT) == S_IFDIR)
// #define S_ISREG(m)  (((m) & S_IFMT) == S_IFREG)
#endif

#include "clipper/core/clipper_util.h"

#include "compat/coot-sysdep.h"
#include "utils/win-compat.hh"

#include "lbg-graph.hh"


// Normally, we would use a const pointer (or a reference). But this
// is mmdb.
//
// return then number of links read (to pass on to
// init_refmac_mon_lib) so it doesn't return 0 (ie. fail) when reading
// links (no new atoms in a link, you see).
// 
int
coot::protein_geometry::init_links(mmdb::mmcif::PData data) {

   int r = 0; 
   for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 
      
      mmdb::mmcif::PCategory cat = data->GetCategory(icat);
      std::string cat_name(cat->GetCategoryName());

      // std::cout << "DEBUG:: init_link is handling " << cat_name << std::endl;

      mmdb::mmcif::PLoop mmCIFLoop = data->GetLoop(cat_name.c_str() );
            
      if (mmCIFLoop == NULL) { 
	 std::cout << "null loop" << std::endl; 
      } else {
	 int n_chiral = 0;
	 if (cat_name == "_chem_link")
	    add_chem_links(mmCIFLoop);
	 if (cat_name == "_chem_link_bond")
	    r += link_bond(mmCIFLoop);
	 if (cat_name == "_chem_link_angle")
	    link_angle(mmCIFLoop);
	 if (cat_name == "_chem_link_tor")
	    link_torsion(mmCIFLoop);
	 if (cat_name == "_chem_link_plane")
	    link_plane(mmCIFLoop);
	 if (cat_name == "_chem_link_chiral")
	    n_chiral = link_chiral(mmCIFLoop);
	 if (n_chiral) {
	    assign_link_chiral_volume_targets();
	 }
      }
   }
   return r;
} 


// References to the modifications
// to the link groups (the modifications
// themselves are in data_mod_list).
void
coot::protein_geometry::add_chem_links(mmdb::mmcif::PLoop mmCIFLoop) {


   char *s;
   int ierr;
   int ierr_tot = 0;

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {
      std::string chem_link_id;
      std::string chem_link_comp_id_1;
      std::string chem_link_mod_id_1;
      std::string chem_link_group_comp_1;
      std::string chem_link_comp_id_2;
      std::string chem_link_mod_id_2;
      std::string chem_link_group_comp_2;
      std::string chem_link_name;
      
      s = mmCIFLoop->GetString("id", j, ierr);
      ierr_tot += ierr;
      if (ierr)
	 std::cout << "WARNING add_chem_links error getting field \"id\"" << std::endl;
      if (s) chem_link_id = s;

      s = mmCIFLoop->GetString("comp_id_1", j, ierr);
      if (ierr)
	 std::cout << "WARNING add_chem_links error getting field \"comp_id_1\"" << std::endl;
      ierr_tot += ierr;
      if (s) chem_link_comp_id_1 = s;

      s = mmCIFLoop->GetString("mod_id_1", j, ierr);
      if (ierr)
	 std::cout << "WARNING add_chem_links error getting field \"mod_id_1\"" << std::endl;
      ierr_tot += ierr;
      if (s) chem_link_mod_id_1 = s;

      s = mmCIFLoop->GetString("group_comp_1", j, ierr);
      if (ierr)
	 std::cout << "WARNING add_chem_links error getting field \"group_comp_1\"" << std::endl;
      ierr_tot += ierr;
      if (s) chem_link_group_comp_1 = s;

      s = mmCIFLoop->GetString("comp_id_2", j, ierr);
      if (ierr)
	 std::cout << "WARNING add_chem_links error getting field \"comp_id_2\"" << std::endl;
      ierr_tot += ierr;
      if (s) chem_link_comp_id_2 = s;

      s = mmCIFLoop->GetString("mod_id_2", j, ierr);
      if (ierr)
	 std::cout << "WARNING add_chem_links error getting field \"mod_id_2\"" << std::endl;
      ierr_tot += ierr;
      if (s) chem_link_mod_id_2 = s;

      s = mmCIFLoop->GetString("group_comp_2", j, ierr);
      if (ierr)
	 std::cout << "WARNING add_chem_links error getting field \"group_comp_2\"" << std::endl;
      ierr_tot += ierr;
      if (s) chem_link_group_comp_2 = s;

      s = mmCIFLoop->GetString("name", j, ierr);
      if (ierr)
	 std::cout << "WARNING add_chem_links error getting field \"name\"" << std::endl;
      ierr_tot += ierr;
      if (s) chem_link_name = s;

      if (ierr_tot == 0) {
	 coot::chem_link clink(chem_link_id,
			       chem_link_comp_id_1, chem_link_mod_id_1, chem_link_group_comp_1,
			       chem_link_comp_id_2, chem_link_mod_id_2, chem_link_group_comp_2,
			       chem_link_name);
	 // std::cout << "Adding to chem_link_map: " << clink << std::endl;
	 // chem_link_vec.push_back(clink);
	 chem_link_map[clink.get_hash_code()].push_back(clink);
      } else {
	 std::cout << "WARNING:: an error occurred when trying to add link: "
		   << "\"" << chem_link_id << "\" "
		   << "\"" << chem_link_comp_id_1 << "\" "
		   << "\"" << chem_link_mod_id_1 << "\" "
		   << "\"" << chem_link_group_comp_1 << "\" "
		   << "\"" << chem_link_comp_id_2 << "\" "
		   << "\"" << chem_link_mod_id_2 << "\" "
		   << "\"" << chem_link_group_comp_2 << "\" "
		   << "\"" << chem_link_name << "\" "
		   << std::endl;
      }
   }
}

int
coot::protein_geometry::link_bond(mmdb::mmcif::PLoop mmCIFLoop) {
   std::string link_id;
   std::string atom_id_1, atom_id_2;
   std::string type;
   mmdb::realtype value_dist, value_dist_esd;
   int atom_1_comp_id, atom_2_comp_id;

   char *s;
   int ierr;
   int ierr_tot = 0;
   int n_link_bonds = 0;

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) { 
      s = mmCIFLoop->GetString("link_id",j,ierr);
      ierr_tot = ierr;
      if (s) link_id = s;

      s = mmCIFLoop->GetString("atom_id_1",j,ierr);
      ierr_tot += ierr;
      if (s) atom_id_1 = s;

      s = mmCIFLoop->GetString("atom_id_2",j,ierr);
      ierr_tot += ierr;
      if (s) atom_id_2 = s;

      ierr = mmCIFLoop->GetInteger(atom_1_comp_id, "atom_1_comp_id",j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetInteger(atom_2_comp_id, "atom_2_comp_id",j);
      ierr_tot += ierr;
   
      ierr = mmCIFLoop->GetReal(value_dist, "value_dist",j);
      ierr_tot += ierr; 
      
      ierr = mmCIFLoop->GetReal(value_dist_esd, "value_dist_esd",j);
      ierr_tot += ierr;

      if (ierr_tot == 0) {
	 link_add_bond(link_id,
		       atom_1_comp_id, atom_2_comp_id,
		       atom_id_1, atom_id_2,
		       value_dist, value_dist_esd);
	 n_link_bonds++;
      } else {
	 std::cout << "problem reading bond mmCIFLoop" << std::endl;
      } 
   }
   return n_link_bonds;
}

void
coot::protein_geometry::link_angle(mmdb::mmcif::PLoop mmCIFLoop) {

   std::string link_id;
   std::string atom_id_1, atom_id_2, atom_id_3;
   std::string type;
   mmdb::realtype value_angle, value_angle_esd;
   int atom_1_comp_id, atom_2_comp_id, atom_3_comp_id;

   char *s;
   int ierr;
   int ierr_tot = 0;

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) { 
      s = mmCIFLoop->GetString("link_id",j,ierr);
      ierr_tot = ierr;
      if (s) link_id = s;

      s = mmCIFLoop->GetString("atom_id_1",j,ierr);
      ierr_tot += ierr;
      if (s) atom_id_1 = s;

      s = mmCIFLoop->GetString("atom_id_2",j,ierr);
      ierr_tot += ierr;
      if (s) atom_id_2 = s;

      s = mmCIFLoop->GetString("atom_id_3",j,ierr);
      ierr_tot += ierr;
      if (s) atom_id_3 = s;

      ierr = mmCIFLoop->GetInteger(atom_1_comp_id, "atom_1_comp_id",j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetInteger(atom_2_comp_id, "atom_2_comp_id",j);
      ierr_tot += ierr;
   
      ierr = mmCIFLoop->GetInteger(atom_3_comp_id, "atom_3_comp_id",j);
      ierr_tot += ierr;
   
      ierr = mmCIFLoop->GetReal(value_angle, "value_angle",j);
      ierr_tot += ierr;
      
      ierr = mmCIFLoop->GetReal(value_angle_esd, "value_angle_esd",j);
      ierr_tot += ierr;

      if (ierr_tot == 0) {
	 link_add_angle(link_id,
		       atom_1_comp_id, atom_2_comp_id, atom_3_comp_id,
		       atom_id_1, atom_id_2, atom_id_3,
		       value_angle, value_angle_esd);
      } else {
	 std::cout << "problem reading link angle mmCIFLoop" << std::endl;
      } 
   }
}


void
coot::protein_geometry::link_torsion(mmdb::mmcif::PLoop mmCIFLoop) {

   std::string link_id;
   std::string  atom_id_1, atom_id_2, atom_id_3, atom_id_4;
   mmdb::realtype value_angle, value_angle_esd;
   int atom_1_comp_id, atom_2_comp_id, atom_3_comp_id, atom_4_comp_id;
   char *s;
   int ierr;
   int ierr_tot = 0;
   int period;
   std::string id("unknown"); // gets get to phi psi, etc

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

      s = mmCIFLoop->GetString("link_id",j,ierr);
      ierr_tot = ierr;
      if (s) link_id = s;
      
      id = "unknown"; // no need to error if "id" was not given
      s = mmCIFLoop->GetString("id",j,ierr);
      if (s) id = s; 
      
      s = mmCIFLoop->GetString("atom_id_1",j,ierr);
      ierr_tot += ierr;
      if (s) atom_id_1 = s;

      s = mmCIFLoop->GetString("atom_id_2",j,ierr);
      ierr_tot += ierr;
      if (s) atom_id_2 = s;

      s = mmCIFLoop->GetString("atom_id_3",j,ierr);
      ierr_tot += ierr;
      if (s) atom_id_3 = s;

      s = mmCIFLoop->GetString("atom_id_4",j,ierr);
      ierr_tot += ierr;
      if (s) atom_id_4 = s;

      ierr = mmCIFLoop->GetInteger(atom_1_comp_id, "atom_1_comp_id",j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetInteger(atom_2_comp_id, "atom_2_comp_id",j);
      ierr_tot += ierr;
   
      ierr = mmCIFLoop->GetInteger(atom_3_comp_id, "atom_3_comp_id",j);
      ierr_tot += ierr;
   
      ierr = mmCIFLoop->GetInteger(atom_4_comp_id, "atom_4_comp_id",j);
      ierr_tot += ierr;
   
      ierr = mmCIFLoop->GetReal(value_angle, "value_angle",j);
      ierr_tot += ierr;
      
      ierr = mmCIFLoop->GetReal(value_angle_esd, "value_angle_esd",j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetInteger(period, "period",j);
      ierr_tot += ierr;

      if (ierr_tot == 0) {
	 link_add_torsion(link_id,
			  atom_1_comp_id, atom_2_comp_id, atom_3_comp_id, atom_4_comp_id,
			  atom_id_1, atom_id_2, atom_id_3, atom_id_4,
			  value_angle, value_angle_esd, period, id);
      } else {
	 std::cout << "problem reading link torsion mmCIFLoop" << std::endl;
      } 
   }
}

int 
coot::protein_geometry::link_chiral(mmdb::mmcif::PLoop mmCIFLoop) {

   int n_chiral = 0;
   std::string chiral_id;
   std::string atom_id_c, atom_id_1, atom_id_2, atom_id_3;
   int atom_c_comp_id, atom_1_comp_id, atom_2_comp_id, atom_3_comp_id;
   int volume_sign;

   char *s; 
   int ierr;
   int ierr_tot = 0;
   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) { 

      // modify a reference (ierr)
      //
      s = mmCIFLoop->GetString("chiral_id",j,ierr);
      ierr_tot += ierr;
      if (s)
	 chiral_id = s; 

      ierr = mmCIFLoop->GetInteger(volume_sign, "volume_sign", j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetInteger(atom_c_comp_id, "atom_centre_comp_id",j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetInteger(atom_1_comp_id, "atom_1_comp_id",j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetInteger(atom_2_comp_id, "atom_2_comp_id",j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetInteger(atom_3_comp_id, "atom_3_comp_id",j);
      ierr_tot += ierr;

      s = mmCIFLoop->GetString("atom_id_centre",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_id_c = s;
      
      s = mmCIFLoop->GetString("atom_id_1",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_id_1 = s;

      if (ierr_tot == 0) {
	 link_add_chiral(chiral_id,
			 atom_c_comp_id,
			 atom_1_comp_id, atom_2_comp_id, atom_3_comp_id,
			 atom_id_c, atom_id_1, atom_id_2, atom_id_3,
			 volume_sign);
	 n_chiral++;
      } else {
	 std::cout << "problem reading link torsion mmCIFLoop" << std::endl;
      } 
   }
   return n_chiral;
}

void
coot::protein_geometry::link_plane(mmdb::mmcif::PLoop mmCIFLoop) {

   std::string link_id;
   std::string atom_id, plane_id;
   mmdb::realtype dist_esd;
   int atom_comp_id;

   char *s; 
   int ierr;
   int ierr_tot = 0;
   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) { 

      // modify a reference (ierr)
      //
      s = mmCIFLoop->GetString("link_id",j,ierr);
      ierr_tot += ierr;
      if (s)
	 link_id = s; 

      s = mmCIFLoop->GetString("atom_id",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_id = s; 

      ierr = mmCIFLoop->GetInteger(atom_comp_id, "atom_comp_id",j);
      ierr_tot += ierr;

      s = mmCIFLoop->GetString("plane_id",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 plane_id = s;

      ierr = mmCIFLoop->GetReal(dist_esd, "dist_esd",j);
      ierr_tot += ierr;

      if (ierr_tot == 0) {
	 link_add_plane(link_id, atom_id, plane_id, atom_comp_id, dist_esd); 
      } else {
	 std::cout << "problem reading link plane mmCIFLoop" << std::endl;
      }
   }
}

void
coot::protein_geometry::link_add_bond(const std::string &link_id,
				      int atom_1_comp_id,
				      int atom_2_comp_id,
				      const std::string &atom_id_1,
				      const std::string &atom_id_2,
				      mmdb::realtype value_dist,
				      mmdb::realtype value_dist_esd) {

   dict_link_bond_restraint_t lbr(atom_1_comp_id,
				  atom_2_comp_id,
				  atom_id_1,
				  atom_id_2,
				  value_dist,
				  value_dist_esd);

//    std::cout << "attempting to add link bond "
// 	     << link_id << " "
// 	     << atom_1_comp_id << " "
// 	     << atom_1_comp_id << " "
// 	     << value_dist << " "
// 	     << value_dist_esd << " dict_link_res_restraints size = "
// 	     << dict_link_res_restraints.size() << std::endl;
   
   short int ifound = 0;
   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {
      if (dict_link_res_restraints[i].link_id == link_id) {
	 ifound = 1;
	 dict_link_res_restraints[i].link_bond_restraint.push_back(lbr);
      }
   }

   if (! ifound) {
      dict_link_res_restraints.push_back(dictionary_residue_link_restraints_t(link_id));
      dict_link_res_restraints[dict_link_res_restraints.size()-1].link_bond_restraint.push_back(lbr);
   }
}

void
coot::protein_geometry::link_add_angle(const std::string &link_id,
				       int atom_1_comp_id,
				       int atom_2_comp_id,
				       int atom_3_comp_id,
				       const std::string &atom_id_1,
				       const std::string &atom_id_2,
				       const std::string &atom_id_3,
				       mmdb::realtype value_dist,
				       mmdb::realtype value_dist_esd) {
   
   dict_link_angle_restraint_t lar(atom_1_comp_id,
				   atom_2_comp_id,
				   atom_3_comp_id,
				   atom_id_1,
				   atom_id_2,
				   atom_id_3,
				   value_dist,
				   value_dist_esd);

   short int ifound = 0;
   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {
      if (dict_link_res_restraints[i].link_id == link_id) {
	 ifound = 1;
	 dict_link_res_restraints[i].link_angle_restraint.push_back(lar);
      }
   }

   if (! ifound) {
      dict_link_res_restraints.push_back(dictionary_residue_link_restraints_t(link_id));
      dict_link_res_restraints[dict_link_res_restraints.size()-1].link_angle_restraint.push_back(lar);
   }
}

void
coot::protein_geometry::link_add_torsion(const std::string &link_id,
					 int atom_1_comp_id,
					 int atom_2_comp_id,
					 int atom_3_comp_id,
					 int atom_4_comp_id,
					 const std::string &atom_id_1,
					 const std::string &atom_id_2,
					 const std::string &atom_id_3,
					 const std::string &atom_id_4,
					 mmdb::realtype value_dist,
					 mmdb::realtype value_dist_esd,
					 int period,
					 const std::string &id) {  // e.g. phi, psi, omega

   dict_link_torsion_restraint_t ltr(atom_1_comp_id,
				     atom_2_comp_id,
				     atom_3_comp_id,
				     atom_4_comp_id,
				     atom_id_1,
				     atom_id_2,
				     atom_id_3,
				     atom_id_4,
				     value_dist,
				     value_dist_esd,
				     period,
				     id);

   short int ifound = 0;
   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {
      if (dict_link_res_restraints[i].link_id == link_id) {
	 ifound = 1;
	 dict_link_res_restraints[i].link_torsion_restraint.push_back(ltr);
      }
   }

   if (! ifound) {
      dict_link_res_restraints.push_back(dictionary_residue_link_restraints_t(link_id));
      dict_link_res_restraints[dict_link_res_restraints.size()-1].link_torsion_restraint.push_back(ltr);
   }
}

void
coot::protein_geometry::link_add_chiral(const std::string &link_id,
					int atom_c_comp_id,
					int atom_1_comp_id,
					int atom_2_comp_id,
					int atom_3_comp_id,
					const std::string &atom_id_c,
					const std::string &atom_id_1,
					const std::string &atom_id_2,
					const std::string &atom_id_3,
					int volume_sign) {

}

void
coot::protein_geometry::link_add_plane(const std::string &link_id,
				       const std::string &atom_id,
				       const std::string &plane_id,
				       int atom_comp_id,
				       double dist_esd) {


   dict_link_plane_restraint_t lpr(atom_id, plane_id, atom_comp_id, dist_esd);
   bool ifound = 0;
   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {

      // std::cout << "comparing  link_ids: " << i << " " << dict_link_res_restraints[i].link_id
      // << " vs " << link_id << std::endl;
      
      if (dict_link_res_restraints[i].link_id == link_id) { // e.g "TRANS"
	 for (unsigned int ip=0; ip<dict_link_res_restraints[i].link_plane_restraint.size(); ip++) {
	    if (dict_link_res_restraints[i].link_plane_restraint[ip].plane_id == plane_id) {
	       ifound = true;

	       // now check the atom_id and comp_id of that atom are not already in this restraint.
	       // Only add this atom to the plane restraint if there is not already an atom there.
	       // If there is, then replace the atom in the restraint.
	       //
	       
	       int found_atom_index = coot::protein_geometry::UNSET_NUMBER;

	       for (unsigned int i_rest_at=0; i_rest_at<dict_link_res_restraints[i].link_plane_restraint[ip].n_atoms(); i_rest_at++) {
		  if (dict_link_res_restraints[i].link_plane_restraint[ip].atom_ids[i_rest_at] == atom_id) { 
		     if (dict_link_res_restraints[i].link_plane_restraint[ip].atom_comp_ids[i_rest_at] == atom_comp_id) {
			found_atom_index = i_rest_at;
			break;
		     }
		  }
	       }

	       if (found_atom_index != coot::protein_geometry::UNSET_NUMBER) {
		  // replace it then
		  dict_link_res_restraints[i].link_plane_restraint[ip].atom_ids[found_atom_index] = atom_id;
		  dict_link_res_restraints[i].link_plane_restraint[ip].atom_comp_ids[found_atom_index] = atom_comp_id;
	       } else {

		  dict_link_res_restraints[i].link_plane_restraint[ip].atom_ids.push_back(atom_id);
		  dict_link_res_restraints[i].link_plane_restraint[ip].atom_comp_ids.push_back(atom_comp_id);
	       }
	       
	       break;
	    }
	 } 
	 if (! ifound) {
	    // we have link_id, but no planes of that id
	    coot::dict_link_plane_restraint_t res(atom_id, plane_id, atom_comp_id, dist_esd);
	    dict_link_res_restraints[i].link_plane_restraint.push_back(res);
	    ifound = true;
	 }
      }
   }
   
   // It was not there.  This should only happen when plane restraints
   // are declared in the dictionary *before* the bonds (or angles).
   // Which is never, thanks to Alexei.
   //
   // So, there was no comp_id found 
   if (! ifound) {
      // add the plae to the newly created dictionary_residue_link_restraints_t
      dict_link_res_restraints.push_back(dictionary_residue_link_restraints_t(link_id));
      coot::dict_link_plane_restraint_t res(atom_id, plane_id, atom_comp_id, dist_esd);
      dict_link_res_restraints[dict_link_res_restraints.size()-1].link_plane_restraint.push_back(res);
   }
}

// throw an error on no such chem_link
// 
std::vector<std::pair<coot::chem_link, bool> >
coot::protein_geometry::matching_chem_link(const std::string &comp_id_1,
					   const std::string &group_1,
					   const std::string &comp_id_2,
					   const std::string &group_2) const {
   bool allow_peptide_link_flag = true;
   return matching_chem_link(comp_id_1, group_1, comp_id_2, group_2, allow_peptide_link_flag);
}

// throw an error on no chem links at all. (20100420, not sure why an
// exception is thrown, why not just return an empty vector?)
//
// Perhaps the throwing of the exception is a hang-over from when this
// function returned a pair, not a vector of chem_link (and
// assocciated order_switch_flags).
// 
std::vector<std::pair<coot::chem_link, bool> > 
coot::protein_geometry::matching_chem_link(const std::string &comp_id_1,
					   const std::string &group_1,
					   const std::string &comp_id_2,
					   const std::string &group_2,
					   bool allow_peptide_link_flag) const {

   bool switch_order_flag = false;
   bool found = false;
   bool debug = false;

   if (debug) {
      std::cout << "---------------------- Here are the chem_links: -----------------"
		<< std::endl;
      print_chem_links(); // prints the chem_link_map
   }

   // This needs to iterate to make the count now that we use a map
   // std::cout << "Testing vs " << chem_link_vec.size() << " chem links\n";

   unsigned int search_hash_code_f = chem_link::make_hash_code(comp_id_1, comp_id_2, group_1, group_2);
   unsigned int search_hash_code_b = chem_link::make_hash_code(comp_id_2, comp_id_1, group_2, group_1);

   if (debug)
      std::cout << "DEBUG:: matching_chem_link() " << search_hash_code_f << " " << search_hash_code_b << " "
		<< comp_id_1 << " " << comp_id_2 << " groups: " << group_1 << " " << group_2 << std::endl;

   // Is this link a TRANS peptide or a CIS?  Both have same group and
   // comp_ids.  Similarly, is is BETA-1-2 or BETA1-4 (etc).  We need
   // to decide later, don't just pick the first one that matches
   // (keep the order switch flag too).
   //

   // "gap" and "symmetry" have hash code 0 (blank strings)

   std::vector<std::pair<coot::chem_link, bool> > matching_chem_links;
   std::map<unsigned int, std::vector<chem_link> >::const_iterator it = chem_link_map.find(search_hash_code_f);
   if (it == chem_link_map.end()) {
      it = chem_link_map.find(search_hash_code_b);
      if (it != chem_link_map.end())
	 switch_order_flag = true; // used?
   }

   if (debug) {
      if (it != chem_link_map.end()) {
	 std::cout << "DEBUG:: matching_chem_link() found the hash at least! " << std::endl;
         std::cout << "DEBUG:: matching_chem_link() Here is the vector of chem links in the map:" << std::endl;
         const std::vector<chem_link> &v = it->second;
         std::vector<chem_link>::const_iterator itv;
         for (itv=v.begin(); itv!=v.end(); itv++) {
            const chem_link &cl = *itv;
            std::cout << "                 " << cl << std::endl;
         }
      } else {
	 std::cout << "DEBUG:: matching_chem_link() failed to find hash " << search_hash_code_f << " "
		   << search_hash_code_b << std::endl;
      }
   }

   if (it == chem_link_map.end()) {

      // NAG-ASN for example

      if (debug)
	 std::cout << "DEBUG:: matching_chem_link() iterator hit the chem_link_map end" << std::endl;

      unsigned int search_bl_1_f  = chem_link::make_hash_code(comp_id_1, comp_id_2, "", group_2);
      unsigned int search_bl_1_b  = chem_link::make_hash_code(comp_id_2, comp_id_1, group_2, "");
      unsigned int search_bl_2_f  = chem_link::make_hash_code(comp_id_1, comp_id_2, group_1, "");
      unsigned int search_bl_2_b  = chem_link::make_hash_code(comp_id_2, comp_id_1, "", group_1);
      unsigned int search_bl_bl_f = chem_link::make_hash_code(comp_id_1, comp_id_2, "", "");
      unsigned int search_bl_bl_b = chem_link::make_hash_code(comp_id_2, comp_id_1, "", "");

      std::set<chem_link> candidate_chem_links;
      std::vector<chem_link>::const_iterator itv;

      it = chem_link_map.find(search_bl_1_f);
      if (it != chem_link_map.end()) {
	 const std::vector<chem_link> &v = it->second;
	 for (itv=v.begin(); itv!=v.end(); itv++)
	    candidate_chem_links.insert(*itv);
      }

      it = chem_link_map.find(search_bl_1_b);
      if (it != chem_link_map.end()) {
	 const std::vector<chem_link> &v = it->second;
	 for (itv=v.begin(); itv!=v.end(); itv++)
	    candidate_chem_links.insert(*itv);
      }
      it = chem_link_map.find(search_bl_2_f);
      if (it != chem_link_map.end()) {
	 const std::vector<chem_link> &v = it->second;
	 for (itv=v.begin(); itv!=v.end(); itv++)
	    candidate_chem_links.insert(*itv);
      }
      it = chem_link_map.find(search_bl_2_b);
      if (it != chem_link_map.end()) {
	 const std::vector<chem_link> &v = it->second;
	 for (itv=v.begin(); itv!=v.end(); itv++)
	    candidate_chem_links.insert(*itv);
      }
      it = chem_link_map.find(search_bl_bl_f);
      if (it != chem_link_map.end()) {
	 const std::vector<chem_link> &v = it->second;
	 for (itv=v.begin(); itv!=v.end(); itv++)
	    candidate_chem_links.insert(*itv);
      }
      it = chem_link_map.find(search_bl_bl_b);
      if (it != chem_link_map.end()) {
	 const std::vector<chem_link> &v = it->second;
	 for (itv=v.begin(); itv!=v.end(); itv++)
	    candidate_chem_links.insert(*itv);
      }

      if (debug) {
	 std::cout << "DEBUG:: matching_chem_link() -------- here with candidate_chem_links size ------- "
		   << candidate_chem_links.size() << std::endl;
	 std::set<chem_link>::const_iterator it;
	 for(it=candidate_chem_links.begin(); it!=candidate_chem_links.end(); it++)
	    std::cout << "   " << *it << std::endl;
      }

      if (candidate_chem_links.size() > 0) {
	 const std::set<chem_link> &v = candidate_chem_links;

	 std::set<chem_link>::const_iterator itv;
	 for (itv=v.begin(); itv!=v.end(); itv++) {
	    const chem_link &cl = *itv;

	    std::pair<bool, bool> match_res =
	       cl.matches_comp_ids_and_groups(comp_id_1, group_1, comp_id_2, group_2);

	    if (match_res.first) {
	       if (cl.Id() != "gap" && cl.Id() != "symmetry") {
		  switch_order_flag = match_res.second;
		  found = true;
		  std::pair<coot::chem_link, bool> p(cl, switch_order_flag);

                  if (debug)
                     std::cout << "DEBUG:: matching_chem_link() ::::::::: pushing back matching chem link " << cl << std::endl;
		  matching_chem_links.push_back(p);
	       }
	    } else {
               if (debug)
                  std::cout << "DEBUG:: matching_chem_links() test for matches_comp_ids_and_groups failed " << std::endl;
            }
	 }
      }

   } else {

      // normal (say, peptide link) hit

      // v: the set of chem links that have this matching hash code
      //
      const std::vector<chem_link> &v = it->second;

      // on a match, set found and add the std::pair<coot::chem_link, bool> to matching_chem_links

      std::vector<chem_link>::const_iterator itv;
      for (itv=v.begin(); itv!=v.end(); itv++) {
	 const chem_link &cl = *itv;

         // std::cout << ":::::::::::::::::::::::::: testing comp_id and group match for " << cl << std::endl;

	 std::pair<bool, bool> match_res =
	    cl.matches_comp_ids_and_groups(comp_id_1, group_1, comp_id_2, group_2);

	 if (debug)
	    std::cout << "DEBUG:: matching_chem_links() ... matching_chem_link: found matching link "
		      << comp_id_1 << " " << comp_id_2 << " " 
		      << cl << std::endl;

	 if (false) // was debug but TMI ATM - we are looking for a bug in above code
	    std::cout << "    checking chem link:                             " << cl << " -> matched: "
		      << match_res.first << " need order-switch: " << match_res.second << std::endl;

	 if (match_res.first) {
	    if (cl.Id() != "gap" && cl.Id() != "symmetry") {
	       if (!cl.is_peptide_link_p() || allow_peptide_link_flag) {
		  switch_order_flag = match_res.second;
		  found = true;
		  std::pair<coot::chem_link, bool> p(cl, switch_order_flag);
                  if (debug)
                     std::cout << "DEBUG:: matching_chem_link(): pushing back chem link " << cl << " " << switch_order_flag
                               << std::endl;
		  matching_chem_links.push_back(p);

	       } else {
		  if (debug)
		     std::cout << "    reject link " << cl.Id() << " on peptide/allow-peptide test " << std::endl;
	       }
	    } else {
	       if (debug) {
		  std::cout << "    reject link \"" << cl.Id() << "\"" << std::endl;
	       }
	    }
	 } else {
            if (debug)
               std::cout << "    reject link " << cl.Id() << " matches_comp_ids_and_groups() found no match" << std::endl;
         }
      }
   }

   // When allow_peptide_link_flag is FALSE, we don't want to hear
   // about not making a link between ASP and VAL etc (otherwise we
   // do).
   // 
   if ( (!found) && (allow_peptide_link_flag)) {
      std::string rte = "INFO:: matching_chem_links() No chem link for groups \"";
      rte += group_1;
      rte += "\" \"";
      rte += group_2;
      rte += "\" and comp_ids \"";
      rte += comp_id_1;
      rte += "\" \"";
      rte += comp_id_2;
      rte += "\"";
      throw std::runtime_error(rte);
   }
   if (debug)
      std::cout << "DEBUG:: matching_chem_link() returns " << matching_chem_links.size()
		<< " matching chem links" << std::endl;
   return matching_chem_links;
}

// throw an error on no such chem_link
// 
std::vector<std::pair<coot::chem_link, bool> >
coot::protein_geometry::matching_chem_link_non_peptide(const std::string &comp_id_1,
						       const std::string &group_1,
						       const std::string &comp_id_2,
						       const std::string &group_2) const {
   return matching_chem_link(comp_id_1, group_1, comp_id_2, group_2, 0);
}

// throw an error on no such chem_link
// 
std::vector<std::pair<coot::chem_link, bool> >
coot::protein_geometry::matching_chem_link_non_peptide(const std::string &comp_id_1,
						       const std::string &group_1,
						       const std::string &comp_id_2,
						       const std::string &group_2,
						       mmdb::Manager *mol) const {
   // HACK FIXME 20150714
   return matching_chem_link(comp_id_1, group_1, comp_id_2, group_2, 0);
}

