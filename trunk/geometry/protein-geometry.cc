/* geometry/protein-geometry.cc
 * 
 * Copyright 2003, 2004, 2005, 2006 The University of York
 * Author: Paul Emsley
 * Copyright 2007, 2008, 2009, 2010 The University of Oxford
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

#include "atom-quads.hh"
#include "protein-geometry.hh"
#include "coot-utils.hh"

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if !defined _MSC_VER
#include <unistd.h>
#else
#define DATADIR "C:/coot/share"
#define PKGDATADIR DATADIR
#define S_ISDIR(m)  (((m) & S_IFMT) == S_IFDIR)
#define S_ISREG(m)  (((m) & S_IFMT) == S_IFREG)
#endif

#include "clipper/core/clipper_util.h"

#include "coot-sysdep.h"
#include "Cartesian.h"

#include "lbg-graph.hh"

// std::string 
// coot::basic_dict_restraint_t::atom_id_1_4c() const {
//    return atom_id_mmdb_expand(atom_id_1_); 
// }

// std::string 
// coot::basic_dict_restraint_t::atom_id_2_4c() const {
//    return atom_id_mmdb_expand(atom_id_2_);
// }

// std::string 
// coot::basic_dict_restraint_t::atom_id_mmdb_expand(const std::string &atomname) const {

//    std::string r;
//    int ilen = atomname.length();

//    if (ilen == 4) return atomname;
   
//    if (ilen == 1) {
//       r = " ";
//       r += atomname;
//       r += "  ";
//    } else {
//       if (ilen == 2) { 
// 	 r = " ";
// 	 r += atomname;
// 	 r += " ";
//       } else {
// 	 if (ilen == 3) {
// 	    r = " ";
// 	    r += atomname;
// 	 } else {
// 	    r = atomname;
// 	 }
//       }
//    }
//    return r;
// }

std::string
coot::atom_id_mmdb_expand(const std::string &atomname) { 
   std::string r;
   int ilen = atomname.length();
      
   if (ilen == 4) return atomname;
      
   if (ilen == 1) {
      r = " ";
      r += atomname;
      r += "  ";
   } else {
      if (ilen == 2) { 
	 r = " ";
	 r += atomname;
	 r += " ";
      } else {
	 if (ilen == 3) {
	    r = " ";
	    r += atomname;
	 } else {
	    r = atomname;
	 }
      }
   }
   return r;
}

std::string
coot::atom_id_mmdb_expand(const std::string &atomname, const std::string &element) {

   std::string r = coot::atom_id_mmdb_expand(atomname);

   if (element.length() == 2 && element[0] != ' ') {
      if (atomname.length() == 1) { // unlikely
	 r = " ";
	 r += atomname;
	 r += "  ";
      } else {
	 if (atomname.length() == 2) {
	    r = atomname;
	    r += "  ";
	 } else {
	    if (atomname.length() == 3) {
	       r = atomname;
	       r += " ";
	    } else {
	       r = atomname;
	    }
	 }
      }
   }
   if (0)  // debug
      std::cout << "Given :" << atomname << ": and element :" <<
	 element << ": returning :" << r << ":" << std::endl;
   return r;
}



coot::basic_dict_restraint_t::basic_dict_restraint_t(const std::string &at1,
						     const std::string &at2) {

   atom_id_1_ = at1;
   atom_id_2_ = at2;
}



// return the number of atoms read (not the number of bonds (because
// that is not a good measure of having read the file properly for
// (for example) CL)).
// 
int
coot::protein_geometry::init_refmac_mon_lib(std::string ciffilename, int read_number_in) {

   int ret_val = 0; 
   CMMCIFFile ciffile;

   // Here we would want to check through dict_res_restraints for the
   // existance of this restraint.  If it does exist, kill it by
   // changing its name to blank.  However, we don't know the name of
   // the restraint yet!  We only know that at the add_atom(),
   // add_bond() [etc] stage.
   

   struct stat buf;
   int istat = stat(ciffilename.c_str(), &buf);
   // Thanks Ezra Peisach for this this bug report

   if (istat != 0) {
      std::cout << "WARNING: in init_refmac_mon_lib file \"" << ciffilename
		<< "\" not found."
		<< "\n";
      return 0; // failure
   }

   if (! S_ISREG(buf.st_mode)) {
      std::cout << "WARNING: in init_refmac_mon_lib \"" << ciffilename
		<< "\" not read.  It is not a regular file.\n";
   } else { 

      int ierr = ciffile.ReadMMCIFFile(ciffilename.c_str());
   
      if (ierr!=CIFRC_Ok) {
	 std::cout << "dirty mmCIF file? " << ciffilename.c_str() << std::endl;
	 std::cout << "    Bad CIFRC_Ok on ReadMMCIFFile" << std::endl;
	 std::cout << "    " << GetErrorDescription(ierr) << std::endl;
	 char        err_buff[1000];
	 std::cout <<  "CIF error rc=" << ierr << " reason:" << 
	    GetCIFMessage (err_buff,ierr) << std::endl;


      } else {
	 if (verbose_mode)
	    std::cout << "There are " << ciffile.GetNofData() << " data in "
		      << ciffilename << std::endl; 
      
	 for(int idata=0; idata<ciffile.GetNofData(); idata++) { 
         
	    PCMMCIFData data = ciffile.GetCIFData(idata);
	    
	    // note that chem_link goes here to:
	    // 
	    if (std::string(data->GetDataName()).substr(0,5) == "link_") {
	       ret_val += init_links(data);
	    }

	    if (std::string(data->GetDataName()).length() > 7) {
	       if (std::string(data->GetDataName()).substr(0,8) == "mod_list") {
		  // this handles all "list_chem_mod"s in the file (just one call)
		  ret_val += add_chem_mods(data);

		  if (0) // debug
		     for (unsigned int i=0; i<chem_mod_vec.size(); i++)
			std::cout << "     " << chem_mod_vec[i] << std::endl;
	       }
	    }
	    
	    if (std::string(data->GetDataName()).length() > 4) {
	       // e.g. mod_NH3 ? 
	       if (std::string(data->GetDataName()).substr(0,4) == "mod_") {
		  ret_val += add_mod(data);
	       }
	    }
	       

	    if (std::string(data->GetDataName()).length() > 16) { 
	       if (std::string(data->GetDataName()).substr(0,17) == "comp_synonym_list") {
		  add_synonyms(data);
	       }
	    }
	    
	    int n_chiral = 0;
	    std::vector<std::string> comp_ids_for_chirals;
	    for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 

	       PCMMCIFCategory cat = data->GetCategory(icat);
	       std::string cat_name(cat->GetCategoryName());
	       
	       // All catagories have loops (AFAICS). 
	       // std::cout << "debug got catagory: " << cat_name << std::endl; 

	       PCMMCIFLoop mmCIFLoop = data->GetLoop(cat_name.c_str() );

	       int n_loop_time = 0;
	       if (mmCIFLoop == NULL) {

		  bool handled = false;
		  if (cat_name == "_chem_comp") {
		     // read the chemical component library which does
		     // not have a loop (the refmac files do) for the
		     // chem_comp info.
		     handled = 1;
		     PCMMCIFStruct structure = data->GetStructure(cat_name.c_str());
		     if (structure) {
			chem_comp_component(structure);
		     }
		  }

		  if (cat_name == "_chem_comp_chir") {
		     handled = 1;
		     PCMMCIFStruct structure = data->GetStructure(cat_name.c_str());
		     if (structure) {
			chem_comp_chir_structure(structure);
		     }
		  } 
		     
		  if (cat_name == "_chem_comp_tor") {
		     handled = 1;
		     PCMMCIFStruct structure = data->GetStructure(cat_name.c_str());
		     if (structure) {
			chem_comp_tor_structure(structure);
		     }
		  }
		  
		  if (! handled) 
		     std::cout << "in init_refmac_mon_lib() null loop for catagory "
			       << cat_name << std::endl; 
		  
	       } else {
               
		  n_loop_time++;

		  // We currently want to stop adding chem comp info
		  // if the chem_comp info comes from mon_lib_list.cif:
		  if (cat_name == "_chem_comp") { 
		     if (read_number_in != coot::protein_geometry::MON_LIB_LIST_CIF)
			chem_comp(mmCIFLoop);
		     else
			simple_mon_lib_chem_comp(mmCIFLoop);
		  }
		  
		  // monomer info, name, number of atoms etc.
		  if (cat_name == "_chem_comp_atom")
		     ret_val += comp_atom(mmCIFLoop); // and at the end pad up the atom names

		  // tree
		  if (cat_name == "_chem_comp_tree")
		     comp_tree(mmCIFLoop);

		  // bond
		  if (cat_name == "_chem_comp_bond")
		     comp_bond(mmCIFLoop);

		  // angle
		  if (cat_name == "_chem_comp_angle")
		     comp_angle(mmCIFLoop);

		  // tor
		  if (cat_name == "_chem_comp_tor")
		     comp_torsion(mmCIFLoop);

		  // chiral
		  if (cat_name == "_chem_comp_chir") {
		     std::pair<int, std::vector<std::string> > chirals = 
			comp_chiral(mmCIFLoop);
		     n_chiral += chirals.first;
		     for (unsigned int ichir=0; ichir<chirals.second.size(); ichir++)
			comp_ids_for_chirals.push_back(chirals.second[ichir]);
		  }
               
		  // plane
		  if (cat_name == "_chem_comp_plane_atom")
		     comp_plane(mmCIFLoop);

	       }
	    }
	    if (n_chiral) {
	       assign_chiral_volume_targets();
	       filter_chiral_centres(comp_ids_for_chirals);
	    }
	 } // for idata
      } // cif file is OK test
   } // is regular file test

   // debug_mods();
   
   return ret_val; // the number of atoms read.
}

void
coot::protein_geometry::chem_comp_component(PCMMCIFStruct structure) {

   int n_tags = structure->GetNofTags();
   std::string cat_name = structure->GetCategoryName();

//     std::cout << "DEBUG: ================= by structure: in category " << cat_name << " there are "
//  	     << n_tags << " tags" << std::endl;

   std::pair<bool, std::string> comp_id(0, "");
   std::pair<bool, std::string> three_letter_code(0, "");
   std::pair<bool, std::string> name(0, "");
   std::pair<bool, std::string> type(0, ""); // aka group?
   int number_of_atoms_all = coot::protein_geometry::UNSET_NUMBER;
   int number_of_atoms_nh  = coot::protein_geometry::UNSET_NUMBER;
   std::pair<bool, std::string> description_level(0, "");

   for (int itag=0; itag<n_tags; itag++) {
      std::string tag = structure->GetTag(itag);
      std::string field = structure->GetField(itag);
      // std::cout << " by structure got tag " << itag << " "
      // << tag << " field: " << f << std::endl;
      if (tag == "id")
	 comp_id = std::pair<bool, std::string> (1,field);
      if (tag == "three_letter_code")
	 three_letter_code = std::pair<bool, std::string> (1,field);
      if (tag == "name")
	 name = std::pair<bool, std::string> (1,field);
      if (tag == "type")
	 type = std::pair<bool, std::string> (1,field);
      if (tag == "descr_level")
	 description_level = std::pair<bool, std::string> (1,field);
      if (tag == "description_level")
	 description_level = std::pair<bool, std::string> (1,field);
      // number of atoms here too.

      if (tag == "number_atoms_all") { 
	 try {
	    number_of_atoms_all = coot::util::string_to_int(field);
	 }
	 catch (std::runtime_error rte) {
	    std::cout << rte.what() << std::endl;
	 }
      }
      if (tag == "number_atoms_nh") { 
	 try {
	    number_of_atoms_nh = coot::util::string_to_int(field);
	 }
	 catch (std::runtime_error rte) {
	    std::cout << rte.what() << std::endl;
	 }
      }
   }

   if (0) 
      std::cout
	 << "comp_id :" << comp_id.first << " :" << comp_id.second << ": "
	 << "three_letter_code :" << three_letter_code.first << " :" << three_letter_code.second << ": "
	 << "name :" << name.first << " :" << name.second << ": "
	 << "type :" << type.first << " :" << type.second << ": "
	 << "description_level :" << description_level.first << " :" << description_level.second << ": "
	 << std::endl;

   if (comp_id.first && three_letter_code.first && name.first) {
      mon_lib_add_chem_comp(comp_id.second, three_letter_code.second,
			    name.second, type.second,
			    number_of_atoms_all, number_of_atoms_nh,
			    description_level.second);
   } else { 
      // std::cout << "oooppps - something missing, not adding that" << std::endl;
   } 

}


// non-looping (single) tor
void
coot::protein_geometry::chem_comp_tor_structure(PCMMCIFStruct structure) {
   
   int n_tags = structure->GetNofTags();
   std::string cat_name = structure->GetCategoryName();

   if (0)
      std::cout << "DEBUG: ================= by chem_comp_tor by structure: in category "
		<< cat_name << " there are "
		<< n_tags << " tags" << std::endl;

   std::pair<bool, std::string> comp_id(0, "");
   std::pair<bool, std::string> torsion_id(0, "");
   std::pair<bool, std::string> atom_id_1(0, "");
   std::pair<bool, std::string> atom_id_2(0, "");
   std::pair<bool, std::string> atom_id_3(0, "");
   std::pair<bool, std::string> atom_id_4(0, "");
   std::pair<bool, int> period(0, 0);
   std::pair<bool, realtype> value_angle(0, 0);
   std::pair<bool, realtype> value_angle_esd(0, 0);
   
   for (int itag=0; itag<n_tags; itag++) {
      std::string tag = structure->GetTag(itag);
      std::string field = structure->GetField(itag);
      // std::cout << " by structure got tag " << itag << " \""
      // << tag << "\" field: \"" << field << "\"" << std::endl;
      if (tag == "comp_id")
	 comp_id = std::pair<bool, std::string> (1,field);
      if (tag == "torsion_id")
	 torsion_id = std::pair<bool, std::string> (1,field);
      if (tag == "atom_id_1")
	 atom_id_1 = std::pair<bool, std::string> (1,field);
      if (tag == "atom_id_2")
	 atom_id_2 = std::pair<bool, std::string> (1,field);
      if (tag == "atom_id_3")
	 atom_id_3 = std::pair<bool, std::string> (1,field);
      if (tag == "atom_id_4")
	 atom_id_4 = std::pair<bool, std::string> (1,field);
      if (tag == "period") { 
	 try { 
	    period = std::pair<bool, int> (1,coot::util::string_to_int(field));
	 }
	 catch (std::runtime_error rte) {
	    std::cout << "WARNING:: not an integer: " << field << std::endl;
	 }
      }
      if (tag == "value_angle") { 
	 try { 
	    value_angle = std::pair<bool, float> (1,coot::util::string_to_float(field));
	 }
	 catch (std::runtime_error rte) {
	    std::cout << "WARNING:: value_angle not an number: " << field << std::endl;
	 }
      }
      if (tag == "value_angle_esd") { 
	 try { 
	    value_angle_esd = std::pair<bool, float> (1,coot::util::string_to_float(field));
	 }
	 catch (std::runtime_error rte) {
	    std::cout << "WARNING:: value_angle_esd not an number: " << field << std::endl;
	 }
      }
   }

   if (comp_id.first && 
       atom_id_1.first && atom_id_2.first && atom_id_3.first && atom_id_4.first &&
       value_angle.first && value_angle_esd.first && 
       period.first) {
      mon_lib_add_torsion(comp_id.second,
			  torsion_id.second,
			  atom_id_1.second,
			  atom_id_2.second,
			  atom_id_3.second,
			  atom_id_4.second,
			  value_angle.second, value_angle_esd.second,
			  period.second);
   } else {
      std::cout << "WARNING:: chem_comp_tor_structure() something bad" << std::endl;
   } 
}

// non-looping (single) chir
void
coot::protein_geometry::chem_comp_chir_structure(PCMMCIFStruct structure) {

   int n_tags = structure->GetNofTags();
   std::string cat_name = structure->GetCategoryName();

   if (0)
      std::cout << "DEBUG: ================= by chem_comp_dot by structure: in category "
		<< cat_name << " there are "
		<< n_tags << " tags" << std::endl;

   std::pair<bool, std::string> comp_id(0, "");
   std::pair<bool, std::string>      id(0, "");
   std::pair<bool, std::string> atom_id_centre(0, "");
   std::pair<bool, std::string> atom_id_1(0, "");
   std::pair<bool, std::string> atom_id_2(0, "");
   std::pair<bool, std::string> atom_id_3(0, "");
   std::pair<bool, std::string> volume_sign(0, "");
   
   for (int itag=0; itag<n_tags; itag++) {
      std::string tag = structure->GetTag(itag);
      std::string field = structure->GetField(itag);
      // std::cout << " by structure got tag " << itag << " \""
      // << tag << "\" field: \"" << field << "\"" << std::endl;
      if (tag == "comp_id")
	 comp_id = std::pair<bool, std::string> (1,field);
      if (tag == "id")
	 id = std::pair<bool, std::string> (1,field);
      if (tag == "atom_id_centre")
	 atom_id_centre = std::pair<bool, std::string> (1,field);
      if (tag == "atom_id_1")
	 atom_id_1 = std::pair<bool, std::string> (1,field);
      if (tag == "atom_id_2")
	 atom_id_2 = std::pair<bool, std::string> (1,field);
      if (tag == "atom_id_3")
	 atom_id_3 = std::pair<bool, std::string> (1,field);
      if (tag == "volume_sign")
	 volume_sign = std::pair<bool, std::string> (1,field);
   }

   if (comp_id.first && atom_id_centre.first &&
       atom_id_1.first && atom_id_2.first && atom_id_3.first &&
       volume_sign.first) {
      mon_lib_add_chiral(comp_id.second,
			 id.second,
			 atom_id_centre.second,
			 atom_id_1.second,
			 atom_id_2.second,
			 atom_id_3.second,
			 volume_sign.second);
   } else {
      std::cout << "WARNING:: chem_comp_chir_structure() something bad" << std::endl;
   } 
}



void
coot::protein_geometry::assign_chiral_volume_targets() {

   for (unsigned int idict=0; idict<dict_res_restraints.size(); idict++) {
      if (dict_res_restraints[idict].has_unassigned_chiral_volumes()) {
	 if (0) 
	    std::cout << "DEBUG:: assign_chiral_volume_targets for dict_res_restraints entry: "
		      << idict << " " << dict_res_restraints[idict].residue_info.comp_id
		      << " has unassigned chiral volumes" << std::endl;
	 dict_res_restraints[idict].assign_chiral_volume_targets();
      }
   }

   assign_link_chiral_volume_targets();
}

// the chiral restraint for this comp_id(s) may need filtering
// (i.e. removing some of them if they are not real chiral centres
// (e.g. from prodrg restraints)).
void
coot::protein_geometry::filter_chiral_centres(const std::vector<std::string> &comp_ids_for_filtering) {

   for (unsigned int ichir=0; ichir<comp_ids_for_filtering.size(); ichir++) {
      int idx = get_monomer_restraints_index(comp_ids_for_filtering[ichir], 0);
      if (idx != -1) { 
	 const coot::dictionary_residue_restraints_t &restraints =
	    dict_res_restraints[idx];
	 std::vector<coot::dict_chiral_restraint_t> new_chirals =
	    filter_chiral_centres(restraints);
	 dict_res_restraints[idx].chiral_restraint = new_chirals;
      }
   }
} 

// Return a filtered list, that is don't include chiral centers that
// are connected to more than one hydrogen.
// 
std::vector<coot::dict_chiral_restraint_t>
coot::protein_geometry::filter_chiral_centres(const dictionary_residue_restraints_t &restraints) {

   std::vector<coot::dict_chiral_restraint_t> v;
   for (unsigned int ichir=0; ichir<restraints.chiral_restraint.size(); ichir++) { 
      int n_H=0;
      for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) {
	 if (restraints.bond_restraint[ib].atom_id_1_4c() ==
	     restraints.chiral_restraint[ichir].atom_id_c_4c()) {
	    if (restraints.element(restraints.bond_restraint[ib].atom_id_2_4c()) == " H")
	       n_H++;
	 }
	 if (restraints.bond_restraint[ib].atom_id_2_4c() ==
	     restraints.chiral_restraint[ichir].atom_id_c_4c()) {
	    if (restraints.element(restraints.bond_restraint[ib].atom_id_1_4c()) == " H")
	       n_H++;
	 }
      }
      if (n_H <= 1)
	 v.push_back(restraints.chiral_restraint[ichir]);
   }
   return v;
} 

void
coot::protein_geometry::assign_link_chiral_volume_targets() {

   for (unsigned int idict=0; idict<dict_link_res_restraints.size(); idict++) {
      if (dict_link_res_restraints[idict].has_unassigned_chiral_volumes()) {
	 dict_link_res_restraints[idict].assign_link_chiral_volume_targets();
      }
   }
}

void
coot::dictionary_residue_restraints_t::clear_dictionary_residue() {

   residue_info = coot::dict_chem_comp_t("", "", "", "", 0, 0, "");
   has_partial_charges_flag = 0;

   // need different constructors.
//    atom_info.resize(0);
//    bond_restraint.resize(0);
//    angle_restraint.resize(0);
//    torsion_restraint.resize(0);
//    chiral_restraint.resize(0);
//    plane_restraint.resize(0);
}

void
coot::protein_geometry::mon_lib_add_chem_comp(const std::string &comp_id,
					      const std::string &three_letter_code,
					      const std::string &name,
					      const std::string &group,
					      int number_atoms_all, int number_atoms_nh,
					      const std::string &description_level) {


//    std::cout << "DEBUG:: in mon_lib_add_chem_comp :"
// 	     << comp_id << ": :" << three_letter_code << ": :"
// 	     << name << ": :" << group << ": :"
// 	     << description_level << ": :" << number_atoms_all << ": :"
// 	     << number_atoms_nh << std::endl;

   coot::dict_chem_comp_t ri(comp_id, three_letter_code, name, group, number_atoms_all, number_atoms_nh,
			     description_level);
   short int ifound = 0;

       for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
       if (dict_res_restraints[i].residue_info.comp_id == comp_id) {
	  if (dict_res_restraints[i].read_number == read_number) { 
	     ifound = 1;
	     dict_res_restraints[i].residue_info = ri;
	     break;
	  } else {
	     // trash the old one then
	     dict_res_restraints[i].clear_dictionary_residue();
	  }
       }
    }

    if (! ifound) {
       // std::cout << "residue not found in mon_lib_add_chem_comp" << std::endl;
       dict_res_restraints.push_back(dictionary_residue_restraints_t(comp_id, read_number));
       dict_res_restraints[dict_res_restraints.size()-1].residue_info = ri;
    }
}


// add to simple_monomer_descriptions not dict_res_restraints.
void
coot::protein_geometry::simple_mon_lib_add_chem_comp(const std::string &comp_id,
					      const std::string &three_letter_code,
					      const std::string &name,
					      const std::string &group,
					      int number_atoms_all, int number_atoms_nh,
					      const std::string &description_level) {


   if (0)
      std::cout << "------- DEBUG:: in simple_mon_lib_add_chem_comp comp_id :"
		<< comp_id << ": three-letter-code :" << three_letter_code << ": name :"
		<< name << ": :" << group << ": descr-lev :"
		<< description_level << ": :" << number_atoms_all << ": :"
		<< number_atoms_nh << std::endl;

   // notice that we also pass the comp_id here (a different constructor needed);
   coot::dict_chem_comp_t ri(comp_id, three_letter_code, name, group, number_atoms_all,
			     number_atoms_nh, description_level);


   std::map<std::string,coot::dictionary_residue_restraints_t>::const_iterator it =
      simple_monomer_descriptions.find(comp_id);

   coot::dictionary_residue_restraints_t blank_res_rest;
   blank_res_rest.residue_info = ri; // now mostly blank
   simple_monomer_descriptions[comp_id] = blank_res_rest;
   
//    std::cout << "Added [residue info :"
// 	     << ri.comp_id << ": :"
// 	     << ri.three_letter_code << ": " 
// 	     << ri.group << "] " 
// 	     << " for key :"
// 	     << comp_id << ":" << " simple_monomer_descriptions size() "
// 	     << simple_monomer_descriptions.size() 
// 	     << std::endl;

}


void
coot::protein_geometry::mon_lib_add_atom(const std::string &comp_id,
					 const std::string &atom_id,
					 const std::string &atom_id_4c,
					 const std::string &type_symbol,
					 const std::string &type_energy,
					 const std::pair<bool, realtype> &partial_charge,
					 const std::pair<bool, clipper::Coord_orth> &model_pos,
					 const std::pair<bool, clipper::Coord_orth> &model_pos_ideal) { 

   // debugging
   bool debug = 0;
   if (debug) { 
      std::cout << "   mon_lib_add_atom  " << comp_id << " :" << atom_id << ": :"
		<< atom_id_4c << ": " << type_symbol << " " << type_energy << " ("
		<< partial_charge.first << "," << partial_charge.second << ")" << std::endl;
   } 

   coot::dict_atom at_info(atom_id, atom_id_4c, type_symbol, type_energy, partial_charge);
//    std::cout << "mon_lib_add_atom " << model_pos.first << " "
// 	     << model_pos.second.format() << std::endl;
//    std::cout << "mon_lib_add_atom " << model_pos_ideal.first << " "
// 	     << model_pos_ideal.second.format() << std::endl;
   
   if (model_pos.first)
      at_info.add_pos(coot::dict_atom::REAL_MODEL_POS, model_pos);
   if (model_pos_ideal.first)
      at_info.add_pos(coot::dict_atom::IDEAL_MODEL_POS, model_pos_ideal);

   bool ifound = 0;
   int this_index = -1; // unset

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
//       std::cout << "comparing comp_ids: :" << dict_res_restraints[i].comp_id
// 		<< ":  :" << comp_id << ":" << std::endl;
      
      if (dict_res_restraints[i].residue_info.comp_id == comp_id) {

// 	 std::cout << "comparing read numbers: "
// 		   << dict_res_restraints[i].read_number << " and "
// 		   << read_number << std::endl;
	 if (dict_res_restraints[i].read_number == read_number) { 
	    ifound = 1;
	    this_index = i;
	    dict_res_restraints[i].atom_info.push_back(at_info);
	    break;
	 } else {
	    // trash the old one then
	    dict_res_restraints[i].clear_dictionary_residue();
	 }
      }
   }

   if (! ifound) {
      // std::cout << "residue not found in mon_lib_add_atom" << std::endl;
      dict_res_restraints.push_back(dictionary_residue_restraints_t(comp_id, read_number));
      this_index = dict_res_restraints.size()-1;
      dict_res_restraints[this_index].atom_info.push_back(at_info);
   }

   if (debug) {
      std::cout << "   dictionary for " << dict_res_restraints[this_index].residue_info.comp_id
		<< " now contains " << dict_res_restraints[this_index].atom_info.size()
		<< " atoms" << std::endl;
      for (unsigned int i=0; i<dict_res_restraints[this_index].atom_info.size(); i++) { 
	 // 	  std::cout << "  " << i << "  " << dict_res_restraints[this_index].atom_info[i]
	 // 		    << std::endl;
      }
   }
}

// Because they were all at origin, for example.
// 
void
coot::protein_geometry::delete_atom_positions(const std::string &comp_id, int pos_type) {
   
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].residue_info.comp_id == comp_id) {
	 for (unsigned int iat=0; iat<dict_res_restraints[i].atom_info.size(); iat++) { 
	    if (pos_type == coot::dict_atom::IDEAL_MODEL_POS)
	       dict_res_restraints[i].atom_info[iat].pdbx_model_Cartn_ideal.first = false;
	    if (pos_type == coot::dict_atom::REAL_MODEL_POS)
	       dict_res_restraints[i].atom_info[iat].model_Cartn.first = false;
	 }
      }
   }
} 

void
coot::dict_atom::add_pos(int pos_type,
			 const std::pair<bool, clipper::Coord_orth> &model_pos) {

   if (pos_type == coot::dict_atom::IDEAL_MODEL_POS)
      pdbx_model_Cartn_ideal = model_pos;
   if (pos_type == coot::dict_atom::REAL_MODEL_POS)
      model_Cartn = model_pos;

}

void
coot::protein_geometry::mon_lib_add_tree(std::string comp_id,
					 std::string atom_id,
					 std::string atom_back,
					 std::string atom_forward,
					 std::string connect_type) {
   bool ifound = 0;
   coot::dict_chem_comp_tree_t ac(atom_id, atom_back, atom_forward, connect_type);
   // std::cout << "mon_lib_add_tree atom :" << atom_id << ":" << std::endl;
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].residue_info.comp_id == comp_id) {
	 ifound = 1;
	 dict_res_restraints[i].tree.push_back(ac);
	 break;
      }
   }
}

void
coot::protein_geometry::mon_lib_add_bond(std::string comp_id,
				   std::string atom_id_1,
				   std::string atom_id_2,
				   std::string type,
				   realtype value_dist, realtype value_dist_esd) {

//     std::cout << "adding bond for " << comp_id << " " << atom_id_1
// 	      << " " << atom_id_2 << " " << type << " " << value_dist
// 	      << " " << value_dist_esd << std::endl;

   // add a bond restraint to the list for comp_id.
   // The list container for comp_id is a dictionary_residue_restraints_t

   add_restraint(comp_id, dict_bond_restraint_t(atom_id_1,
						atom_id_2,
						type,
						value_dist,
						value_dist_esd));
}

void
coot::protein_geometry::mon_lib_add_bond_no_target_geom(std::string comp_id,
							std::string atom_id_1,
							std::string atom_id_2,
							std::string type) { 

   add_restraint(comp_id, dict_bond_restraint_t(atom_id_1, atom_id_2, type));
}


void 
coot::protein_geometry::add_restraint(std::string comp_id, const dict_bond_restraint_t &restr) { 

   // if comp is in the list, simply push back restr, 
   // 
   // if not, then push back a dict_bond_restraint_t for it, passing
   // the comp_id. 

   bool ifound = 0;

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) { 
      if (dict_res_restraints[i].residue_info.comp_id == comp_id) { 
	 ifound = 1;
	 dict_res_restraints[i].bond_restraint.push_back(restr); 
	 break;
      }
   } 

   // it was not there
   if (! ifound) { 
      dict_res_restraints.push_back(dictionary_residue_restraints_t(comp_id, read_number));
      // add the bond to the newly created dictionary_residue_restraints_t
      dict_res_restraints[dict_res_restraints.size()-1].bond_restraint.push_back(restr);
   }
}

void 
coot::protein_geometry::add_restraint(std::string comp_id, const dict_angle_restraint_t &restr) { 

   // if comp is in the list, simply push back restr, 
   // 
   // if not, then push back a dict_bond_restraint_t for it, passing
   // the comp_id. 

   short int ifound = 0;

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) { 
      if (dict_res_restraints[i].residue_info.comp_id == comp_id) { 
	 ifound = 1;
	 dict_res_restraints[i].angle_restraint.push_back(restr); 
	 break;
      }
   } 

   // it was not there
   if (! ifound) { 
      dict_res_restraints.push_back(dictionary_residue_restraints_t(comp_id, read_number));
      // add the angle to the newly created dictionary_residue_restraints_t
      dict_res_restraints[dict_res_restraints.size()-1].angle_restraint.push_back(restr);
   }
}


std::vector<coot::dict_torsion_restraint_t>
coot::dictionary_residue_restraints_t::get_non_const_torsions(bool include_hydrogen_torsions_flag) const {

   std::vector<coot::dict_torsion_restraint_t> v;
   for (unsigned int i=0; i<torsion_restraint.size(); i++) {
      if (! torsion_restraint[i].is_const()) {
	 if (include_hydrogen_torsions_flag) { 
	    v.push_back(torsion_restraint[i]);
	 } else {
	    // only add this torsion if neither of end atoms of the torsion are hydrogen.
	    if (!is_hydrogen(torsion_restraint[i].atom_id_1()))
	       if (!is_hydrogen(torsion_restraint[i].atom_id_4()))
		  v.push_back(torsion_restraint[i]);
	 }
      }
   }
   return v;
}

bool
coot::dict_torsion_restraint_t::is_pyranose_ring_torsion() const {

   // Needs fixup for PDBv3
   bool status = false;
   std::string ring_atoms[6] = { " C1 ", " C2 ", " C3 ", " C4 ", " C5 ", " O5 " };

   int n_matches = 0;
   for (unsigned int i=0; i<6; i++) { 
      if (atom_id_2_4c() == ring_atoms[i])
	 n_matches++;
      if (atom_id_3_4c() == ring_atoms[i])
	 n_matches++;
   }
   if (n_matches == 2)
      status = true;
   return status;
}

bool
coot::dict_link_torsion_restraint_t::is_pyranose_ring_torsion() const {

   // Needs fixup for PDBv3
   bool status = false;
   std::string ring_atoms[6] = { " C1 ", " C2 ", " C3 ", " C4 ", " C5 ", " O5 " };

   int n_matches = 0;
   for (unsigned int i=0; i<6; i++) { 
      if (atom_id_2_4c() == ring_atoms[i])
	 n_matches++;
      if (atom_id_3_4c() == ring_atoms[i])
	 n_matches++;
   }
   if (n_matches == 2)
      status = true;
   return status;
} 




bool
coot::dict_torsion_restraint_t::is_const() const {

   bool const_flag = 0;
   if (id_.length() > 5) {
      std::string bit = id_.substr(0,5);
      if (bit == "CONST")
	 const_flag = 1;
      if (bit == "const")
	 const_flag = 1;
   }
   return const_flag;
}


std::string
coot::dictionary_residue_restraints_t::atom_name_for_tree_4c(const std::string &atom_id) const {

   std::string r = atom_id;
   for (unsigned int iat=0; iat<atom_info.size(); iat++) {
      if (atom_info[iat].atom_id == atom_id) {
	 r = atom_info[iat].atom_id_4c;
      }
   }
   return r;
}


// look up the atom id in the atom_info (dict_atom vector).
// Return "" on no atom found with name atom_name;
// 
std::string
coot::dictionary_residue_restraints_t::element(const std::string &atom_name) const {

   std::string r = "";
   for (unsigned int iat=0; iat<atom_info.size(); iat++) {
      if (atom_info[iat].atom_id_4c == atom_name) {
	 r = atom_info[iat].type_symbol;
	 break;
      }
   }
   if (r.length() == 1)
      r = " " + r;
   
   // std::cout << " dictionary_residue_restraints_t::element()"
   // 	     << " on atom name :" << atom_name << ": returns :" << r << ":" << std::endl;
   return r;
}

// likewise look up the energy type.  Return "" on no atom fould
// with that atom_name.
// 
std::string
coot::dictionary_residue_restraints_t::type_energy(const std::string &atom_name) const {

   std::string r = "";
   for (unsigned int iat=0; iat<atom_info.size(); iat++) {
      if (atom_info[iat].atom_id_4c == atom_name) {
	 r = atom_info[iat].type_energy;
	 break;
      }
   }
   return r;
}

std::vector<std::string>
coot::dictionary_residue_restraints_t::neighbours(const std::string &atom_name, bool allow_hydrogen_neighbours_flag) const {

   std::vector<std::string> n;
   for (unsigned int i=0; i<bond_restraint.size(); i++) { 
      if (bond_restraint[i].atom_id_1() == atom_name) {
	 if (allow_hydrogen_neighbours_flag || ! is_hydrogen(bond_restraint[i].atom_id_2())) {
	    n.push_back(bond_restraint[i].atom_id_2());
	 }
      }
      if (bond_restraint[i].atom_id_2() == atom_name) {
	 if (allow_hydrogen_neighbours_flag || ! is_hydrogen(bond_restraint[i].atom_id_1())) {
	    n.push_back(bond_restraint[i].atom_id_1());
	 }
      }
   }
   return n;
} 





std::vector<std::string>
coot::dictionary_residue_restraints_t::get_attached_H_names(const std::string &atom_name) const {

   std::vector<std::string> v;
   for (unsigned int i=0; i<bond_restraint.size(); i++) { 
      if (bond_restraint[i].atom_id_1() == atom_name) {
	 // if (element(bond_restraint[i].atom_id_2()) == " H")
	 if (is_hydrogen(bond_restraint[i].atom_id_2()))
	    v.push_back(bond_restraint[i].atom_id_2());
      }
      if (bond_restraint[i].atom_id_2() == atom_name) {
	 // if (element(bond_restraint[i].atom_id_1()) == " H")
	 if (is_hydrogen(bond_restraint[i].atom_id_1()))
	    v.push_back(bond_restraint[i].atom_id_1());
      }
   }
   return v;
}



void 
coot::protein_geometry::add_restraint(std::string comp_id,
				      const dict_torsion_restraint_t &restr) { 

   // if comp is in the list, simply push back restr, 
   // 
   // if not, then push back a dict_bond_restraint_t for it, passing
   // the comp_id. 

   short int ifound = 0;

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) { 
      if (dict_res_restraints[i].residue_info.comp_id == comp_id) { 
	 ifound = 1;
	 dict_res_restraints[i].torsion_restraint.push_back(restr); 
	 break;
      }
   } 

   // it was not there
   if (! ifound) { 
      dict_res_restraints.push_back(dictionary_residue_restraints_t(comp_id, read_number));
      // add the angle to the newly created dictionary_residue_restraints_t
      dict_res_restraints[dict_res_restraints.size()-1].torsion_restraint.push_back(restr);
   }
}

void 
coot::protein_geometry::add_restraint(std::string comp_id,
				      const dict_chiral_restraint_t &restr) { 

   // if comp is in the list, simply push back restr, 
   // 
   // if not, then push back a dict_bond_restraint_t for it, passing
   // the comp_id. 

   short int ifound = 0;

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) { 
      if (dict_res_restraints[i].residue_info.comp_id == comp_id) { 
	 ifound = 1;
	 dict_res_restraints[i].chiral_restraint.push_back(restr); 
	 break;
      }
   } 

   // it was not there
   if (! ifound) { 
      dict_res_restraints.push_back(dictionary_residue_restraints_t(comp_id, read_number));
      // add the angle to the newly created dictionary_residue_restraints_t
      dict_res_restraints[dict_res_restraints.size()-1].chiral_restraint.push_back(restr);
   }
}

void
coot::protein_geometry::mon_lib_add_angle(std::string comp_id,
				    std::string atom_id_1,
				    std::string atom_id_2,
				    std::string atom_id_3,
				    realtype value_angle,
				    realtype value_angle_esd) {

//    std::cout << "adding angle " << comp_id <<  " " << atom_id_1
// 	     << " " << atom_id_2 << " " << atom_id_3 << " "
// 	     << value_angle << std::endl;

   add_restraint(comp_id, dict_angle_restraint_t(atom_id_1,
						atom_id_2,
						atom_id_3,
						value_angle,
						value_angle_esd));
}

void
coot::protein_geometry::mon_lib_add_torsion(std::string comp_id,
					    std::string torsion_id,
					    std::string atom_id_1,
					    std::string atom_id_2,
					    std::string atom_id_3,
					    std::string atom_id_4,
					    realtype value_angle,
					    realtype value_angle_esd,
					    int period) {

   if (0)
      std::cout << "adding torsion " << comp_id <<  " " << atom_id_1
		<< " " << atom_id_2 << " " << atom_id_3 << " "
		<< atom_id_4 << " "
		<< "value_angle: " << value_angle
		<< ", value_angle_esd: " << value_angle_esd
		<< ", period: " << period << std::endl;
   
   add_restraint(comp_id, dict_torsion_restraint_t(torsion_id,
						   atom_id_1,
						   atom_id_2,
						   atom_id_3,
						   atom_id_4,
						   value_angle,
						   value_angle_esd,
						   period));

}

// static int
int
coot::protein_geometry::chiral_volume_string_to_chiral_sign(const std::string &volume_sign) {

   int volume_sign_int = coot::dict_chiral_restraint_t::CHIRAL_VOLUME_RESTRAINT_VOLUME_SIGN_UNASSIGNED;
   if (volume_sign.length() > 3) { 

       if (volume_sign.substr(0,3) == "pos") { 
	  volume_sign_int = 1;
       }
       if (volume_sign.substr(0,3) == "neg") { 
	  volume_sign_int = -1;
       }
       if (volume_sign.substr(0,3) == "POS") { 
	  volume_sign_int = 1;
       }
       if (volume_sign.substr(0,3) == "NEG") { 
	  volume_sign_int = -1;
       }
       if (volume_sign == "both" || volume_sign == "BOTH") { 
	  volume_sign_int = coot::dict_chiral_restraint_t::CHIRAL_RESTRAINT_BOTH;
       }
    }
    return volume_sign_int;
}


// and the reverse 
//
// static 
std::string 
coot::protein_geometry::make_chiral_volume_string(int chiral_sign) { 

  std::string s;

  if (chiral_sign == -1) 
    s = "negative";
  if (chiral_sign == 1) 
    s = "positive";
  if (chiral_sign == coot::dict_chiral_restraint_t::CHIRAL_RESTRAINT_BOTH)
    s = "both";
  return s;
}



void
coot::protein_geometry::mon_lib_add_chiral(std::string comp_id,
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
	  add_restraint(comp_id, dict_chiral_restraint_t(id, atom_id_centre,
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
					  const std::string &plane_id,
					  const std::string &atom_id,
					  const realtype &dist_esd) {
   
//    std::cout << "adding plane " << comp_id <<  " " << plane_id
// 	     << " " << atom_id << endl;

   short int ifound = 0;

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) { 
      if (dict_res_restraints[i].residue_info.comp_id == comp_id) {
	 for (unsigned int ip=0; ip<dict_res_restraints[i].plane_restraint.size(); ip++) {
	    if (dict_res_restraints[i].plane_restraint[ip].plane_id == plane_id) { 
	       ifound = 1;
	       dict_res_restraints[i].plane_restraint[ip].push_back_atom(atom_id);
	       break;
	    }
	 }
	 if (! ifound ) {
	    // we have comp_id, but no planes of that id.
	    coot::dict_plane_restraint_t res(plane_id, atom_id, dist_esd);
	    dict_res_restraints[i].plane_restraint.push_back(res);
	    ifound = 1;
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
      dict_res_restraints.push_back(dictionary_residue_restraints_t(comp_id, read_number));
      coot::dict_plane_restraint_t res(plane_id, atom_id, dist_esd);
      dict_res_restraints[dict_res_restraints.size()-1].plane_restraint.push_back(res);
   }
}

std::ostream& coot::operator<<(std::ostream&s, coot::dict_plane_restraint_t rest) {

   s << "[plane-restraint: " << rest.plane_id << " " << rest.dist_esd() << " {"
     << rest.n_atoms() << " atoms} ";
   for (unsigned int iatom=0; iatom<rest.n_atoms(); iatom++) {
      s << ":" << rest[iatom] << ": ";
   }
   s << "]";
   return s;
}

std::ostream&
coot::operator<<(std::ostream &s, const coot::dict_torsion_restraint_t &rest) {
   s << "[torsion-restraint: " << rest.id() << " "
     << rest.atom_id_1_4c() << " "
     << rest.atom_id_2_4c() <<  " "
     << rest.atom_id_3_4c() <<  " "
     << rest.atom_id_4_4c() <<  " "
     << rest.angle() << " " 
     << rest.esd() << " " 
     << rest.periodicity();
   if (rest.is_const())
      s << " CONST ";
   s << "]";
   return s;
}

// hack for mac, ostream problems
std::string
coot::dict_torsion_restraint_t::format() const {

   std::string s = "[torsion-restraint: ";
   s +=  id();
   s += " ";
   s +=  atom_id_1_4c();
   s +=  " ";
   s +=  atom_id_2_4c();
   s +=   " ";
   s +=  atom_id_3_4c();
   s +=   " ";
   s +=  atom_id_4_4c();
   s +=   " ";
   s +=  coot::util::float_to_string(angle());
   s +=  " " ;
   s +=  coot::util::float_to_string(esd());
   s +=  " ";
   s +=  coot::util::int_to_string(periodicity());
   if (is_const())
      s +=  " CONST ";
   s +=  "]";
   return s;
}




// We currently want to stop adding chem comp info
// if the chem_comp info comes from mon_lib_list.cif:
//
// This is the function that we use read things other than
// MON_LIB_LIST_CIF.  i.e. bespoke ligand dictionaries.
// 
void
coot::protein_geometry::chem_comp(PCMMCIFLoop mmCIFLoop) {

   int ierr = 0;

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

      // modify a reference (ierr)
      // 
      std::string three_letter_code;
      std::string name;
      std::string group; // e.g. "L-peptide"
      int number_atoms_all;
      int number_atoms_nh;
      std::string description_level = "None";
      
      std::string id; // gets stored as comp_id, but is labelled "id"
		      // in the cif file chem_comp block.

      int ierr_tot = 0;
      char *s = mmCIFLoop->GetString("id", j, ierr);
      ierr_tot += ierr;
      if (s) 
	 id = s;
      
      s = mmCIFLoop->GetString("three_letter_code", j, ierr);
      ierr_tot += ierr;
      if (s)
	 three_letter_code = s;

      s = mmCIFLoop->GetString("name", j, ierr);
      ierr_tot += ierr;
      if (s)
	 name = s;

      s = mmCIFLoop->GetString("group", j, ierr);
      ierr_tot += ierr;
      if (s)
	 group = s; // e.g. "L-peptide"

      ierr = mmCIFLoop->GetInteger(number_atoms_all, "number_atoms_all", j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetInteger(number_atoms_nh, "number_atoms_nh", j);
      ierr_tot += ierr;

      // If desc_level is in the file, extract it, otherwise set it to "None"
      // 
      s = mmCIFLoop->GetString("desc_level", j, ierr);
      if (! ierr) {
	 if (s) { 
	    description_level = s;  // e.g. "." for full, I think
	 } else {
	    // if the description_level is "." in the cif file, then
	    // GetString() does not fail, but s is set to NULL
	    // (slightly surprising).
	    description_level = ".";
	 } 
      } else {
	 std::cout << "WARNING:: desc_level was not set " << j << std::endl;
      } 


      if (ierr_tot != 0) {
	 std::cout << "oops:: ierr_tot was " << ierr_tot << std::endl;
      } else { 

	 if (0) {
	    std::cout << "in chem_comp, adding chem_comp, description_level is :"
		      << description_level << ":" << std::endl;

	    std::cout << "Adding :" << id << ": :" << three_letter_code << ": :" << name << ": :"
		      << group << ": " << number_atoms_all << " "
		      << number_atoms_nh << " :" << description_level << ":" << std::endl;
	 }

	 delete_mon_lib(id); // delete it if it exists already.

	 mon_lib_add_chem_comp(id, three_letter_code, name,
			       group, number_atoms_all, number_atoms_nh,
			       description_level);
      }
   }
}

// We currently want to stop adding chem comp info
// if the chem_comp info comes from mon_lib_list.cif:
//
// This is the function that we use read MON_LIB_LIST_CIF
// 
void
coot::protein_geometry::simple_mon_lib_chem_comp(PCMMCIFLoop mmCIFLoop) {

   int ierr = 0;
   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) { 
      // modify a reference (ierr)
      // 
      char *s = mmCIFLoop->GetString("id", j, ierr);
      std::string three_letter_code;
      std::string name;
      std::string group; // e.g. "L-peptide"
      int number_atoms_all;
      int number_atoms_nh;
      std::string description_level = "None";

      if (ierr == 0) {
	 int ierr_tot = 0;
	 std::string comp_id(s);
	 s = mmCIFLoop->GetString("three_letter_code", j, ierr);
	 ierr_tot += ierr;
	 if (s)
	    three_letter_code = s;
	 else {
	    three_letter_code = "";
// 	    std::cout << "WARNING:: failed to get 3-letter code for comp_id: "
// 		      << comp_id << " error: " << ierr << std::endl;
	 }

	 s = mmCIFLoop->GetString("name", j, ierr);
	 ierr_tot += ierr;
	 if (s)
	    name = s;

	 s = mmCIFLoop->GetString("group", j, ierr);
	 ierr_tot += ierr;
	 if (s)
	    group = s; // e.g. "L-peptide"

	 ierr = mmCIFLoop->GetInteger(number_atoms_all, "number_atoms_all", j);
	 ierr_tot += ierr;

	 ierr = mmCIFLoop->GetInteger(number_atoms_nh, "number_atoms_nh", j);
	 ierr_tot += ierr;

	 s = mmCIFLoop->GetString("desc_level", j, ierr);
	 if (! ierr)
	    if (s)
	       description_level = s;  // e.g. "." for full, I think

	 if (ierr_tot == 0) {
	    simple_mon_lib_add_chem_comp(comp_id, three_letter_code, name,
					 group, number_atoms_all, number_atoms_nh,
					 description_level);

	 }
      }
   }
}

// return the number of atoms.
int 
coot::protein_geometry::comp_atom(PCMMCIFLoop mmCIFLoop) {

   // If the number of atoms with partial charge matches the number of
   // atoms, then set a flag in the residue that this monomer has
   // partial charges.

   int ierr = 0;
   int n_atoms = 0;
   int n_atoms_with_partial_charge = 0;
   // count the following to see if we need to delete the model/ideal
   // atom coordinates because there were all at the origin.
   int n_origin_ideal_atoms = 0; 
   int n_origin_model_atoms = 0;
   std::string comp_id; // used to delete atoms (if needed).
   
   std::string comp_id_for_partial_charges = "unset"; // unassigned.

   

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) { 

      // modify a reference (ierr)
      // 
      char *s = mmCIFLoop->GetString("comp_id",j,ierr);
      std::string atom_id;
      std::string type_symbol; 
      std::string type_energy = "unset";
      std::pair<bool, realtype> partial_charge(0,0);

      std::pair<bool, int> pdbx_align(0, 0);
      int xalign;
      int ierr_optional;
      
      std::pair<bool, std::string> pdbx_aromatic_flag(0,"");
      std::pair<bool, std::string> pdbx_leaving_atom_flag(0,"");
      std::pair<bool, std::string> pdbx_stereo_config_flag(0,"");
      std::pair<bool, clipper::Coord_orth> pdbx_model_Cartn_ideal;
      std::pair<bool, clipper::Coord_orth> model_Cartn;

      if (ierr == 0) {
	 int ierr_tot = 0; 
	 comp_id = std::string(s); // e.g. "ALA"

	 s = mmCIFLoop->GetString("atom_id",j,ierr);
	 ierr_tot += ierr;
	 if (s) {
	    atom_id = s; 
	 }
			
	 s = mmCIFLoop->GetString("type_symbol",j,ierr);
	 ierr_tot += ierr;
	 if (s) {
	    type_symbol = s; 
	 }

	 int ierr_optional = 0;
	 s = mmCIFLoop->GetString("type_energy",j,ierr_optional);
	 if (s) {
	    type_energy = s; 
	 }

	 ierr_optional = mmCIFLoop->GetInteger(xalign, "pdbx_align", j);
	 if (! ierr_optional)
	    pdbx_align = std::pair<bool, int> (1, xalign);

	 s = mmCIFLoop->GetString("pdbx_aromatic_flag", j, ierr_optional);
	 if (s) {
	    if (! ierr_optional) 
	       pdbx_aromatic_flag = std::pair<bool, std::string> (1, s);
	 } 

	 s = mmCIFLoop->GetString("pdbx_leaving_atom_flag", j, ierr_optional);
	 if (s) {
	    if (! ierr_optional) 
	       pdbx_leaving_atom_flag = std::pair<bool, std::string> (1, s);
	 } 

	 s = mmCIFLoop->GetString("pdbx_stereo_config_flag", j, ierr_optional);
	 if (s) {
	    if (! ierr_optional) 
	       pdbx_stereo_config_flag = std::pair<bool, std::string> (1, s);
	 }

	 realtype x,y,z;
	 pdbx_model_Cartn_ideal.first = 0;
	 int ierr_optional_x = mmCIFLoop->GetReal(x, "pdbx_model_Cartn_x_ideal", j);
	 int ierr_optional_y = mmCIFLoop->GetReal(y, "pdbx_model_Cartn_y_ideal", j);
	 int ierr_optional_z = mmCIFLoop->GetReal(z, "pdbx_model_Cartn_z_ideal", j);
	 if (ierr_optional_x == 0)
	    if (ierr_optional_y == 0)
	       if (ierr_optional_z == 0) { 
		  if (close_float_p(x, 0.0))
		     if (close_float_p(z, 0.0))
			if (close_float_p(z, 0.0))
			   n_origin_ideal_atoms++;
		  pdbx_model_Cartn_ideal = std::pair<bool, clipper::Coord_orth>(1, clipper::Coord_orth(x,y,z));
	       }

	 model_Cartn.first = 0;
	 ierr_optional_x = mmCIFLoop->GetReal(x, "model_Cartn_x", j);
	 ierr_optional_y = mmCIFLoop->GetReal(y, "model_Cartn_y", j);
	 ierr_optional_z = mmCIFLoop->GetReal(z, "model_Cartn_z", j);
	 if (ierr_optional_x == 0)
	    if (ierr_optional_y == 0)
	       if (ierr_optional_z == 0) { 
		  if (close_float_p(x, 0.0))
		     if (close_float_p(z, 0.0))
			if (close_float_p(z, 0.0))
			   n_origin_model_atoms++;
		  model_Cartn = std::pair<bool, clipper::Coord_orth>(1, clipper::Coord_orth(x,y,z));
	       }

	 // Try simple x, y, z (like the refmac dictionary that Garib sent has)
	 // 
	 if (model_Cartn.first == 0) {
	    ierr_optional_x = mmCIFLoop->GetReal(x, "x", j);
	    ierr_optional_y = mmCIFLoop->GetReal(y, "y", j);
	    ierr_optional_z = mmCIFLoop->GetReal(z, "z", j);
	    if (ierr_optional_x == 0)
	       if (ierr_optional_y == 0)
		  if (ierr_optional_z == 0) {
		     if (close_float_p(x, 0.0))
			if (close_float_p(z, 0.0))
			   if (close_float_p(z, 0.0))
			      n_origin_model_atoms++;
		     model_Cartn = std::pair<bool, clipper::Coord_orth>(1, clipper::Coord_orth(x,y,z));
		  }
	 }

	 // It's possible that this data type is not in the cif file,
	 // so don't fail if we can't read it.

	 realtype tmp_var;
	 ierr = mmCIFLoop->GetReal(tmp_var, "partial_charge", j);
	 if (ierr == 0) {
	    partial_charge = std::pair<bool, float>(1, tmp_var);
	    n_atoms_with_partial_charge++;
	 }

	 if (ierr_tot == 0) {

	    std::string padded_name = comp_atom_pad_atom_name(atom_id, type_symbol);
// 	    std::cout << "comp_atom_pad_atom_name: in :" << atom_id << ": out :"
// 		      << padded_name << ":" << std::endl;
	    n_atoms++;
	    if (comp_id_for_partial_charges != "bad match") { 
	       if (comp_id_for_partial_charges == "unset") {
		  comp_id_for_partial_charges = comp_id;
	       } else {
		  if (comp_id != comp_id_for_partial_charges) {
		     comp_id_for_partial_charges == "bad match";
		  }
	       }
	    }

	    if (0) 
	       std::cout << "debug:: calling mon_lib_add_atom: "
			 << ":" << comp_id << ":  "
			 << ":" << atom_id << ":  "
			 << ":" << padded_name << ":  "
			 << ":" << type_symbol << ":  "
			 << std::endl;
	    mon_lib_add_atom(comp_id, atom_id, padded_name, type_symbol, type_energy,
			     partial_charge, model_Cartn, pdbx_model_Cartn_ideal);

	 } else {
	    std::cout << " error on read " << ierr_tot << std::endl;
	 } 
      }
   }

   if (n_atoms_with_partial_charge == n_atoms) {
      if (comp_id_for_partial_charges != "unset") {
	 if (comp_id_for_partial_charges != "bad match") {
	    for (unsigned int id=0; id<dict_res_restraints.size(); id++) {
	       if (dict_res_restraints[id].residue_info.comp_id == comp_id_for_partial_charges) {
		  dict_res_restraints[id].set_has_partial_charges(1);
	       } 
	    }
	 }
      }
   }

   // Now we need to check that the atom ideal or model coordinates were not at the origin.
   // 
   if (n_origin_ideal_atoms > 2) // trash all ideal atoms
      delete_atom_positions(comp_id, coot::dict_atom::IDEAL_MODEL_POS);
   if (n_origin_model_atoms > 2) // trash all model/real atoms
      delete_atom_positions(comp_id, coot::dict_atom::REAL_MODEL_POS);
   
   return n_atoms;
}


void
coot::protein_geometry::comp_tree(PCMMCIFLoop mmCIFLoop) {

   std::string comp_id;
   std::string atom_id;
   std::string atom_back;
   std::string atom_forward; 
   std::string connect_type; 
   char *s; 
   
   int ierr;
   int ierr_tot = 0; 

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

      // reset them for next round, we don't want the to keep the old values if
      // they were not set.
      
      comp_id = "";
      atom_id = "";
      atom_back = "";
      atom_forward = "";
      connect_type = "";
      
      // modify a reference (ierr)
      // 
      s = mmCIFLoop->GetString("comp_id",j,ierr);
      ierr_tot += ierr;
      if (s)
	 comp_id = s; 

      s = mmCIFLoop->GetString("atom_id",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_id = s; 

      s = mmCIFLoop->GetString("atom_back",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_back = s; 

      s = mmCIFLoop->GetString("atom_forward",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_forward = s; 

      s = mmCIFLoop->GetString("connect_type",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 connect_type = s;

      if (ierr == 0) {
	 std::string padded_name_atom_id = atom_name_for_tree_4c(comp_id, atom_id);
	 std::string padded_name_atom_back = atom_name_for_tree_4c(comp_id, atom_back);
	 std::string padded_name_atom_forward = atom_name_for_tree_4c(comp_id, atom_forward);
	 mon_lib_add_tree(comp_id, padded_name_atom_id, padded_name_atom_back,
			  padded_name_atom_forward, connect_type); 
      } 
   }
}

// look up atom_id in the atom atom_info (dict_atom vector) of the comp_id restraint
// 
std::string
coot::protein_geometry::atom_name_for_tree_4c(const std::string &comp_id, const std::string &atom_id) const {

   std::string r = atom_id;
   for (int id=(dict_res_restraints.size()-1); id >=0; id--) {
      if (dict_res_restraints[id].residue_info.comp_id == comp_id) {
	 r = dict_res_restraints[id].atom_name_for_tree_4c(atom_id);
	 break;
      }
   }
   return r;
}


int
coot::protein_geometry::comp_bond(PCMMCIFLoop mmCIFLoop) {

   bool verbose_output = 0; // can be passed, perhaps.
   std::string comp_id;
   std::string atom_id_1, atom_id_2;
   std::string type;
   realtype value_dist = -1.0, value_dist_esd = -1.0;

   char *s; 
   int nbond = 0;
   int comp_id_index = -1; // not found initially

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) { 

      int ierr;
      int ierr_tot = 0; 
      int ierr_tot_for_ccd = 0;  // CCD comp_bond entries don't have
			         // target geometry (dist and esd) but we
			         // want to be able to read them in
			         // anyway (to get the bond orders for
			         // drawing).

   
      // modify a reference (ierr)
      //

      s = mmCIFLoop->GetString("comp_id",j,ierr);
      ierr_tot += ierr;
      ierr_tot_for_ccd =+ ierr;
      if (s) { 
	 comp_id = s;
	 for (int id=(dict_res_restraints.size()-1); id >=0; id--) {
	    if (dict_res_restraints[id].residue_info.comp_id == comp_id) {
	       comp_id_index = id;
	       break;
	    }
	 }
      }

      s = mmCIFLoop->GetString("atom_id_1", j, ierr);
      ierr_tot += ierr;
      ierr_tot_for_ccd =+ ierr;
      if (s) 
	 atom_id_1 = get_padded_name(s, comp_id_index);

      s = mmCIFLoop->GetString("atom_id_2", j, ierr);
      ierr_tot += ierr;
      ierr_tot_for_ccd =+ ierr;
      if (s) 
	 atom_id_2 = get_padded_name(s, comp_id_index);

      s = mmCIFLoop->GetString("type", j, ierr);
      if (s) 
	 type = s;
      // perhaps it was in the dictionary as "value_order"?
      if (ierr) {
	 s = mmCIFLoop->GetString("value_order", j, ierr);
	 if (! ierr) {
	    if (s) { // just in case (should not be needed).
	       std::string ss(s);
	       // convert from Chemical Component Dictionary value_order to
	       // refmac monomer library chem_comp_bond types. 
	       if (ss == "SING")
		  type = "single";
	       if (ss == "DOUB")
		  type = "double";
	       if (ss == "TRIP")
		  type = "triple";
	       
	       // Chemical Chemical Dictionary also has an aromatic flag, so
	       // we can have bonds that are (for example) "double"
	       // "aromatic" Hmm!  Food for thought.
	       
	       // Metal bonds are "SING" (i.e. CCD doesn't have metal
	       // bonds).
	    }
	 } 
      } 
      ierr_tot += ierr;
      ierr_tot_for_ccd =+ ierr;


      ierr = mmCIFLoop->GetReal(value_dist, "value_dist", j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetReal(value_dist_esd, "value_dist_esd", j);
      ierr_tot += ierr;

      if (ierr_tot == 0) {
	 mon_lib_add_bond(comp_id, atom_id_1, atom_id_2,
			  type, value_dist, value_dist_esd); 
	 nbond++;
      } else {

	 if (! ierr_tot_for_ccd) {

	    mon_lib_add_bond_no_target_geom(comp_id, atom_id_1, atom_id_2, type);
	    
	 } else {
	    // Hopeless - nothing worked...
	    
	    // std::cout << "DEBUG::  ierr_tot " << ierr_tot << std::endl;
	    if (verbose_output) { 
	       std::cout << "Fail on read " << atom_id_1 << ": :" << atom_id_2 << ": :"
			 << type << ": :" << value_dist << ": :" << value_dist_esd
			 << ":" << std::endl;
	    }
	 }
      } 
   }
   return nbond;
}

std::string
coot::protein_geometry::get_padded_name(const std::string &atom_id,
					const int &comp_id_index) const {
   std::string s;
   if (comp_id_index < 0) {
      std::cout << "ERROR:: disaster in get_padded_name for comp_id_index "
		<< comp_id_index << " and atom name " << atom_id << std::endl;
      return s;
   } else {
      for (unsigned int iat=0; iat<dict_res_restraints[comp_id_index].atom_info.size(); iat++) {
	 if (dict_res_restraints[comp_id_index].atom_info[iat].atom_id == atom_id) {
	    s = dict_res_restraints[comp_id_index].atom_info[iat].atom_id_4c;
	    break;
	 }
      }
   }
   return s;
}


void
coot::protein_geometry::comp_angle(PCMMCIFLoop mmCIFLoop) {

   std::string comp_id;
   std::string atom_id_1, atom_id_2, atom_id_3;
   realtype value_angle, value_angle_esd;
   int comp_id_index = -1;

   char *s; 
   int ierr;
   int ierr_tot = 0;

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) { 

      // modify a reference (ierr)
      //

      s = mmCIFLoop->GetString("comp_id",j,ierr);
      ierr_tot += ierr;
      if (s) { 
	 comp_id = s;
	 for (int id=(dict_res_restraints.size()-1); id >=0; id--) {
	    if (dict_res_restraints[id].residue_info.comp_id == comp_id) {
	       comp_id_index = id;
	       break;
	    }
	 }
      }

      s = mmCIFLoop->GetString("atom_id_1",j,ierr);
      ierr_tot += ierr;
      if (s) {
	 atom_id_1 = get_padded_name(std::string(s), comp_id_index);
      }

      s = mmCIFLoop->GetString("atom_id_2",j,ierr);
      ierr_tot += ierr;
      if (s) { 
	 atom_id_2 = get_padded_name(std::string(s), comp_id_index);
      }

      s = mmCIFLoop->GetString("atom_id_3",j,ierr);
      ierr_tot += ierr;
      if (s) {
	 atom_id_3 = get_padded_name(std::string(s), comp_id_index);
      }

      ierr = mmCIFLoop->GetReal(value_angle, "value_angle",j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetReal(value_angle_esd, "value_angle_esd",j);
      ierr_tot += ierr;

      if (ierr_tot == 0) {
	 mon_lib_add_angle(comp_id, atom_id_1, atom_id_2, atom_id_3,
			   value_angle, value_angle_esd);
// 	 std::cout << "added angle " << comp_id << " "
// 		   << atom_id_1 << " " 
// 		   << atom_id_2 << " " 
// 		   << atom_id_3 << " " 
// 		   << value_angle << " " 
// 		   << value_angle_esd << std::endl;
      }
   }
} 

void
coot::protein_geometry::comp_torsion(PCMMCIFLoop mmCIFLoop) {

   std::string comp_id, id;
   std::string atom_id_1, atom_id_2, atom_id_3, atom_id_4;
   realtype value_angle, value_angle_esd;
   int period;

   char *s; 
   int ierr;
   int ierr_tot = 0;
   int comp_id_index = -1;

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) { 

      // modify a reference (ierr)
      //

      s = mmCIFLoop->GetString("comp_id",j,ierr);
      ierr_tot += ierr;
      if (s) { 
	 comp_id = s;
	 for (int id=(dict_res_restraints.size()-1); id >=0; id--) {
	    if (dict_res_restraints[id].residue_info.comp_id == comp_id) {
	       comp_id_index = id;
	       break;
	    }
	 }
      }

      s = mmCIFLoop->GetString("id",j,ierr);
      ierr_tot += ierr;
      if (s)
	 id = s; 

      s = mmCIFLoop->GetString("atom_id_1",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_id_1 = get_padded_name(s, comp_id_index);

      s = mmCIFLoop->GetString("atom_id_2",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_id_2 = get_padded_name(s, comp_id_index);

      s = mmCIFLoop->GetString("atom_id_3",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_id_3 = get_padded_name(s, comp_id_index);

      s = mmCIFLoop->GetString("atom_id_4",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_id_4 = get_padded_name(s, comp_id_index);

      ierr = mmCIFLoop->GetReal(value_angle, "value_angle",j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetReal(value_angle_esd, "value_angle_esd",j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetInteger(period, "period",j);
      ierr_tot += ierr;

      if (ierr_tot == 0) { 
	 // we don't want to add restraints for non-variable (CONST) torsion
	 // angles (e.g. in Roberto's DAC).
	 // So, reject if comp_id starts with "CONST" or "const".

	 
// 	 short int add_it = 0;
// 	 if (id.length() > 5) {
// 	    std::string bit = id.substr(0, 5);
// 	    if ((bit != "CONST" && bit != "const")) {
// 	       add_it = 1;
// 	    }
// 	 } else {
// 	    add_it = 1;
// 	 }
	 

// 	 if (add_it == 0)
// 	    std::cout << "DEBUG:: rejecting torsion " << comp_id << " " << atom_id_1 << " " << atom_id_2
// 		      << " " << atom_id_3 << " " << atom_id_4 << " " << value_angle << std::endl;
// 	 else 
// 	    std::cout << "DEBUG:: adding torsion " << id << " " << comp_id << " " << atom_id_1
// 		      << " " << atom_id_2
// 		      << " " << atom_id_3 << " " << atom_id_4 << " " << value_angle << std::endl;
   
	 mon_lib_add_torsion(comp_id, id,
			     atom_id_1, atom_id_2,
			     atom_id_3, atom_id_4,
			     value_angle, value_angle_esd, period); 
      }
   }
} 


std::pair<int, std::vector<std::string> >
coot::protein_geometry::comp_chiral(PCMMCIFLoop mmCIFLoop) {

   std::string comp_id;
   std::string id, atom_id_centre;
   std::string atom_id_1, atom_id_2, atom_id_3;
   std::string volume_sign;
   int n_chiral = 0;
   std::vector<std::string> comp_id_vector;

   char *s; 
   int ierr;
   int ierr_tot = 0;
   int comp_id_index = -1;

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) { 

      // modify a reference (ierr)
      //

      s = mmCIFLoop->GetString("comp_id",j,ierr);
      ierr_tot += ierr;
      if (s) { 
	 comp_id = s;
	 for (int id=(dict_res_restraints.size()-1); id >=0; id--) {
	    if (dict_res_restraints[id].residue_info.comp_id == comp_id) {
	       comp_id_index = id;
	       break;
	    }
	 }
      }

      s = mmCIFLoop->GetString("id",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 id = s; 

      s = mmCIFLoop->GetString("atom_id_centre",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_id_centre = get_padded_name(s, comp_id_index);

      s = mmCIFLoop->GetString("atom_id_1",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_id_1 = get_padded_name(s, comp_id_index);

      s = mmCIFLoop->GetString("atom_id_2",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_id_2 = get_padded_name(s, comp_id_index);

      s = mmCIFLoop->GetString("atom_id_3",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_id_3 = get_padded_name(s, comp_id_index);

      s = mmCIFLoop->GetString("volume_sign",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 volume_sign = s;

      if (ierr_tot == 0) { 
	 mon_lib_add_chiral(comp_id, id, atom_id_centre,
			    atom_id_1,
			    atom_id_2,
			    atom_id_3,
			    volume_sign);
	 // add comp_id to comp_id_vector if it is not there already.
	 if (std::find(comp_id_vector.begin(), comp_id_vector.end(), comp_id) ==
	     comp_id_vector.end())
	    comp_id_vector.push_back(comp_id);
	 n_chiral++;
      }
   }

   return std::pair<int, std::vector<std::string> > (n_chiral, comp_id_vector);
}

void
coot::protein_geometry::comp_plane(PCMMCIFLoop mmCIFLoop) {

   std::string comp_id;
   std::string atom_id, plane_id;
   realtype dist_esd; 

   char *s; 
   int ierr;
   int ierr_tot = 0;
   int comp_id_index = -1;

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) { 

      // modify a reference (ierr)
      //
      s = mmCIFLoop->GetString("comp_id",j,ierr);
      ierr_tot += ierr;
      if (s) { 
	 comp_id = s;
	 for (int id=(dict_res_restraints.size()-1); id >=0; id--) {
	    if (dict_res_restraints[id].residue_info.comp_id == comp_id) {
	       comp_id_index = id;
	       break;
	    }
	 }
      }

      s = mmCIFLoop->GetString("atom_id",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 atom_id = get_padded_name(s, comp_id_index);

      s = mmCIFLoop->GetString("plane_id",j,ierr);
      ierr_tot += ierr;
      if (s) 
	 plane_id = s;

      ierr = mmCIFLoop->GetReal(dist_esd, "dist_esd",j);
      ierr_tot += ierr;

      if (ierr_tot == 0) {
	 mon_lib_add_plane(comp_id, plane_id, atom_id, dist_esd);
      } else {
	 std::cout << "problem reading comp plane" << std::endl;
      } 
   }
} 


void
coot::protein_geometry::info() const {


   std::cout << "::::: MONOMER GEOMETRY:" << std::endl;
   for (int idr=0; idr<size(); idr++) {
      // ejd says that "restraints" has an "n" in it.  Fixed.
      std::cout << dict_res_restraints[idr].residue_info.comp_id << std::endl;
      std::cout << "   " << dict_res_restraints[idr].bond_restraint.size()
		<< " bond restraints " << std::endl;
      std::cout << "   " << dict_res_restraints[idr].angle_restraint.size()
		<< " angle restraints " << std::endl;
      std::cout << "   " << dict_res_restraints[idr].torsion_restraint.size()
		<< " torsion restraints " << std::endl;
      std::cout << "   " << dict_res_restraints[idr].plane_restraint.size()
		<< " plane restraints " << std::endl;
//       for (int i=0; i<dict_res_restraints[idr].plane_restraint.size(); i++) {
// 	 for (int j=0; j<dict_res_restraints[idr].plane_restraint[i].atom_ids.size(); j++) {
// 	    std::cout << "      " << dict_res_restraints[idr].plane_restraint[i].plane_id
// 		      << " " << dict_res_restraints[idr].plane_restraint[i].atom_ids[j]
// 		      << std::endl;
// 	 }
//       }
   }

   std::cout << "::::: LINK GEOMETRY:" << std::endl;
   for (unsigned int idr=0; idr<dict_link_res_restraints.size(); idr++) {
      std::cout << dict_link_res_restraints[idr].link_id << std::endl;
      std::cout << "   " << dict_link_res_restraints[idr].link_bond_restraint.size()
		<< " link bond restraits " << std::endl;
      std::cout << "   " << dict_link_res_restraints[idr].link_angle_restraint.size()
		<< " link angle restraits " << std::endl;
      std::cout << "   " << dict_link_res_restraints[idr].link_torsion_restraint.size()
		<< " link torsion restraits " << std::endl;
      std::cout << "   " << dict_link_res_restraints[idr].link_plane_restraint.size()
		<< " link plane restraits " << std::endl;
   }
   
} 

// 
int
coot::protein_geometry::add_chem_mods(PCMMCIFData data) {

   int n_mods = 0;
   for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 
      
      PCMMCIFCategory cat = data->GetCategory(icat);
      std::string cat_name(cat->GetCategoryName());
      PCMMCIFLoop mmCIFLoop = data->GetLoop(cat_name.c_str() );
            
      if (mmCIFLoop == NULL) { 
	 std::cout << "null loop" << std::endl; 
      } else {
	 int n_chiral = 0;
	 if (cat_name == "_chem_mod")
	    n_mods += add_chem_mod(mmCIFLoop);
      }
   }
   return n_mods;
}


// 
void
coot::protein_geometry::add_synonyms(PCMMCIFData data) {

   for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 
      
      PCMMCIFCategory cat = data->GetCategory(icat);
      std::string cat_name(cat->GetCategoryName());
      PCMMCIFLoop mmCIFLoop = data->GetLoop(cat_name.c_str() );
            
      if (mmCIFLoop == NULL) { 
	 std::cout << "null loop" << std::endl; 
      } else {
	 int n_chiral = 0;
	 if (cat_name == "_chem_comp_synonym") {
	    add_chem_comp_synonym(mmCIFLoop);
	 }
      }
   }
}

void 
coot::protein_geometry::add_chem_comp_synonym(PCMMCIFLoop mmCIFLoop) {

   int ierr = 0;
   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

      int ierr_tot = 0;
      std::string comp_id;
      std::string comp_alternative_id;
      std::string mod_id;

      char *s = mmCIFLoop->GetString("comp_id", j, ierr);
      ierr_tot += ierr;
      if (s) comp_id = s;

      s = mmCIFLoop->GetString("comp_alternative_id", j, ierr);
      ierr_tot += ierr;
      if (s) comp_alternative_id = s;

      s = mmCIFLoop->GetString("mod_id", j, ierr);
      ierr_tot += ierr;
      if (s) mod_id = s;

      if (ierr_tot == 0) {
	 coot::protein_geometry::residue_name_synonym rns(comp_id,
							  comp_alternative_id,
							  mod_id);
	 residue_name_synonyms.push_back(rns);
      } 
   }
} 




int 
coot::protein_geometry::add_chem_mod(PCMMCIFLoop mmCIFLoop) {

   int n_chem_mods = 0;

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

      int ierr_tot = 0;
      int ierr;
      
      std::string id;
      std::string name;
      std::string comp_id; // often "." (default)
      std::string group_id;
      char *s;

      s = mmCIFLoop->GetString("id", j, ierr);
      ierr_tot += ierr;
      if (s) id = s;
      
      s = mmCIFLoop->GetString("name", j, ierr);
      ierr_tot += ierr;
      if (s) name = s;
      
      s = mmCIFLoop->GetString("comp_id", j, ierr);
      ierr_tot += ierr;
      if (s) comp_id = s;
      
      s = mmCIFLoop->GetString("group_id", j, ierr);
      ierr_tot += ierr;
      if (s) group_id = s;

      if (ierr_tot == 0) {
	 coot::list_chem_mod mod(id, name, comp_id, group_id);
	 chem_mod_vec.push_back(mod);
	 n_chem_mods++;
      } 
   }
   return n_chem_mods;
}


// Normally, we would use a const pointer (or a reference). But this
// is mmdb.
//
// return then number of links read (to pass on to
// init_refmac_mon_lib) so it doesn't return 0 (ie. fail) when reading
// links (no new atoms in a link, you see).
// 
int
coot::protein_geometry::init_links(PCMMCIFData data) {

   int r = 0; 
   for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 
      
      PCMMCIFCategory cat = data->GetCategory(icat);
      std::string cat_name(cat->GetCategoryName());

      // std::cout << "DEBUG:: init_link is handling " << cat_name << std::endl;

      PCMMCIFLoop mmCIFLoop = data->GetLoop(cat_name.c_str() );
            
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
coot::protein_geometry::add_chem_links(PCMMCIFLoop mmCIFLoop) {


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
	 // std::cout << "Adding to chem_link_vec: " << clink << std::endl;
	 chem_link_vec.push_back(clink);
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

std::pair<bool, bool>
coot::chem_link::matches_comp_ids_and_groups(const std::string &comp_id_1,
					     const std::string &group_1,
					     const std::string &comp_id_2,
					     const std::string &group_2) const {

   bool debug = 0;
   if (debug) { 
      std::cout << "   ------ DEBUG:: in matches_comp_ids_and_groups "
		<< id << " " << chem_link_name << ": input comp_ids "
		<< comp_id_1 << " and " << comp_id_2 << " vs :"
		<< chem_link_comp_id_1 << ": :"
		<< chem_link_comp_id_2 << ":" << std::endl; 
      std::cout << "      " << chem_link_name << ": input groups "
		<< group_1 << " and " << group_2 << " vs :"
		<< chem_link_group_comp_1 << ": :"
		<< chem_link_group_comp_2 << ":" << std::endl;
   }

   bool match = 0;
   bool order_switch = 0;

   std::string local_group_1 = group_1; 
   std::string local_group_2 = group_2;

   // chem_links specify "peptide" or "pyranose", but comp_groups are "L-peptide"/"D-pyranose".
   // So allow them to match.
   // 
   if (local_group_1 == "L-peptide")
      local_group_1 = "peptide";
   if (local_group_2 == "L-peptide")
      local_group_2 = "peptide";
   if (local_group_1 == "D-pyranose")
      local_group_1 = "pyranose";
   if (local_group_2 == "D-pyranose")
      local_group_2 = "pyranose";
   
   if (((chem_link_group_comp_1 == "") || (chem_link_group_comp_1 == local_group_1)) &&
       ((chem_link_group_comp_2 == "") || (chem_link_group_comp_2 == local_group_2)))
      if (((chem_link_comp_id_1 == "") || (chem_link_comp_id_1 == comp_id_1)) &&
	  ((chem_link_comp_id_2 == "") || (chem_link_comp_id_2 == comp_id_2)))
	 match = 1;

   if (((chem_link_group_comp_1 == "DNA/RNA") && (local_group_1 == "RNA") && 
	(chem_link_group_comp_1 == "DNA/RNA") && (local_group_2 == "RNA")))
      match = 1;

   if (((chem_link_group_comp_1 == "DNA/RNA") && (local_group_1 == "DNA") && 
	(chem_link_group_comp_1 == "DNA/RNA") && (local_group_2 == "DNA")))
      match = 1;

   if (debug) 
      std::cout << "    matches_comp_ids_and_groups given comp_id_1: \""
		<< comp_id_1 << "\" and group_1: \"" << group_1 << "\" and comp_id_2: \"" 
		<< comp_id_2 << "\" and group_2: \"" << group_2 
		<< "\" returns " << match << std::endl;
	
   if (match == 1)
      return std::pair<bool, bool>(match, order_switch);
      
   // And what about if the residues come here backward? We should
   // report a match and that they should be reversed to the calling
   // function?  

   // reverse index 
   if (((chem_link_group_comp_1 == "") || (chem_link_group_comp_1 == local_group_2)) &&
       ((chem_link_group_comp_2 == "") || (chem_link_group_comp_2 == local_group_1)))
      if (((chem_link_comp_id_1 == "") || (chem_link_comp_id_1 == comp_id_2)) &&
	  ((chem_link_comp_id_2 == "") || (chem_link_comp_id_2 == comp_id_1))) { 
	 match = 1;
	 order_switch = 1;
      }
   
   return std::pair<bool, bool>(match, order_switch);
}

std::ostream& coot::operator<<(std::ostream &s, coot::chem_link lnk) {

   s << "[chem_link: id: " << lnk.id
     << " [comp: " << lnk.chem_link_comp_id_1 << " group: " << lnk.chem_link_group_comp_1
     << " mod: " << lnk.chem_link_mod_id_1 << "] to "
     << " [comp: " << lnk.chem_link_comp_id_2 << " group: " << lnk.chem_link_group_comp_2
     << " mod: " << lnk.chem_link_mod_id_2 << "] " << lnk.chem_link_name << "]";
   return s; 
}

std::ostream& coot::operator<<(std::ostream &s, coot::list_chem_mod mod) {

   s << "[list_chem_mod: id: " << mod.id << " " 
     << "name: " << mod.name 
     << " comp_id :" << mod.comp_id 
     << ": group_id: " << mod.group_id
     << "]";
   return s;
} 


int
coot::protein_geometry::link_bond(PCMMCIFLoop mmCIFLoop) {
   std::string link_id;
   std::string atom_id_1, atom_id_2;
   std::string type;
   realtype value_dist, value_dist_esd;
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
coot::protein_geometry::link_angle(PCMMCIFLoop mmCIFLoop) {

   std::string link_id;
   std::string atom_id_1, atom_id_2, atom_id_3;
   std::string type;
   realtype value_angle, value_angle_esd;
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
coot::protein_geometry::link_torsion(PCMMCIFLoop mmCIFLoop) {

   std::string link_id;
   std::string  atom_id_1, atom_id_2, atom_id_3, atom_id_4;
   realtype value_angle, value_angle_esd;
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
coot::protein_geometry::link_chiral(PCMMCIFLoop mmCIFLoop) {

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
coot::protein_geometry::link_plane(PCMMCIFLoop mmCIFLoop) {

   std::string link_id;
   std::string atom_id, plane_id;
   realtype dist_esd;
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
				      realtype value_dist,
				      realtype value_dist_esd) {

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
				       realtype value_dist,
				       realtype value_dist_esd) {
   
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
					 realtype value_dist,
					 realtype value_dist_esd,
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


   dict_link_plane_restraint_t lpr(atom_id,plane_id, atom_comp_id, dist_esd);
   short int ifound = 0;
   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {
      if (dict_link_res_restraints[i].link_id == link_id) { // e.g "TRANS"
	 for (unsigned int ip=0; ip<dict_link_res_restraints[i].link_plane_restraint.size(); ip++) {
	    if (dict_link_res_restraints[i].link_plane_restraint[ip].plane_id == plane_id) {
	       ifound = 1;

	       // now check the atom_id and comp_id of that atom are not already in this restraint.
	       // Only add this atom to the plane restraint if there is not already an atom there.
	       // If there is, then replace the atom in the restraint.
	       //
	       
	       int found_atom_index = coot::protein_geometry::UNSET_NUMBER;

	       for (int i_rest_at=0; i_rest_at<dict_link_res_restraints[i].link_plane_restraint[ip].n_atoms(); i_rest_at++) {
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
	    ifound = 1;
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
   return matching_chem_link(comp_id_1, group_1, comp_id_2, group_2, 1);
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

   bool switch_order_flag = 0;
   bool found = 0;
   bool debug = 0;

//    if (debug) { 
//       std::cout << "---------------------- Here are the chem_links: -----------------"
// 		<< std::endl;
//       print_chem_links();
//    }
   
   // Is this link a TRANS peptide or a CIS?  Both have same group and
   // comp_ids.  Similarly, is is BETA-1-2 or BETA1-4 (etc).  We need
   // to decide later, don't just pick the first one that matches
   // (keep the order switch flag too).
   // 
   std::vector<std::pair<coot::chem_link, bool> > matching_chem_links;
   for (unsigned int i_chem_link=0; i_chem_link<chem_link_vec.size(); i_chem_link++) {
      std::pair<bool, bool> match_res =
	 chem_link_vec[i_chem_link].matches_comp_ids_and_groups(comp_id_1, group_1,
								comp_id_2, group_2);
      if (match_res.first) {

	 if (debug)
	    std::cout << "... matching link " << comp_id_1 << " " << comp_id_2 << " " 
		      << chem_link_vec[i_chem_link] << std::endl;
	 
	 // make sure that this link id is not a (currently) useless one.
	 if (chem_link_vec[i_chem_link].Id() != "gap" &&
	     chem_link_vec[i_chem_link].Id() != "symmetry") { 
	    coot::chem_link clt = chem_link_vec[i_chem_link];
	    if (!clt.is_peptide_link_p() || allow_peptide_link_flag) {
	       switch_order_flag = match_res.second;
	       found = 1;
	       std::pair<coot::chem_link, bool> p(clt, switch_order_flag);
	       matching_chem_links.push_back(p);

	       // no! We want all of them - not just the first glycosidic bond that matches
	       // i.e. don't return just BETA1-2 when we have a BETA1-4.
	       // break; // we only want to find one chem link for this comp_id pair.
	       
	    } else {
	       if (debug)
		  std::cout << "reject link on peptide/allow-peptide test " << std::endl;
	    } 
	 } else {
	    if (debug) { 
	       std::cout << "reject link \"" << chem_link_vec[i_chem_link].Id() << "\""
			 << std::endl;
	    } 
	 } 
      }
   }

   // When allow_peptide_link_flag is FALSE, we don't want to hear
   // about not making a link between ASP and VAL etc (otherwise we
   // do).
   // 
   if ( (!found) && (allow_peptide_link_flag)) {
      std::string rte = "INFO:: No chem link for groups \"";
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

// return "" on failure.
// no order switch is considered.
// 
std::string
coot::protein_geometry::find_glycosidic_linkage_type(CResidue *first, CResidue *second) const {

   // Fixup needed for PDBv3

   double critical_dist = 3.0; // A, less than that and Coot should
			       // try to make the bond.
   PPCAtom res_selection_1 = NULL;
   PPCAtom res_selection_2 = NULL;
   int i_no_res_atoms_1;
   int i_no_res_atoms_2;
   double d;
   std::vector<coot::glycosidic_distance> close;
   
   first->GetAtomTable( res_selection_1, i_no_res_atoms_1);
   second->GetAtomTable(res_selection_2, i_no_res_atoms_2);

   for (int i1=0; i1<i_no_res_atoms_1; i1++) { 
      clipper::Coord_orth a1(res_selection_1[i1]->x,
			     res_selection_1[i1]->y,
			     res_selection_1[i1]->z);
      for (int i2=0; i2<i_no_res_atoms_2; i2++) {
	 clipper::Coord_orth a2(res_selection_2[i2]->x,
				res_selection_2[i2]->y,
				res_selection_2[i2]->z);
	 d = (a1-a2).lengthsq();
	 if (d < critical_dist*critical_dist) {
	    close.push_back(coot::glycosidic_distance(res_selection_1[i1],
						      res_selection_2[i2],
						      sqrt(d)));
	 }
      }
   }

   std::sort(close.begin(), close.end());

   // if you consider to uncomment this to debug a repulsion instead
   // the forming of a glycosidic bond, consider the residue numbering
   // of the residues involved: the "residue 1" should have the O4 and
   // the "residue 2" (+1 residue number) should have the C1.
   // 
   if (0) { 
      std::cout << "DEBUG:: number of sorted distances in glycosidic_linkage: "
		<< close.size() << std::endl;
      for (unsigned int i=0; i<close.size(); i++) {
	 std::cout << "#### glyco close: " << close[i].distance << "  "
		   << close[i].at1->GetChainID() << " " 
		   << close[i].at1->GetSeqNum() << " " 
		   << close[i].at1->GetAtomName() << " " 
		   << " to "
		   << close[i].at2->GetChainID() << " " 
		   << close[i].at2->GetSeqNum() << " " 
		   << close[i].at2->GetAtomName() << " " 
		   << std::endl;
      }
   }

   std::string link_type("");

   // glyco_chiral constructor can throw an exception
   try { 
   
      float smallest_link_dist = 99999.9;
      for (unsigned int i=0; i<close.size(); i++) {
	 std::string name_1(close[i].at1->name);
	 std::string name_2(close[i].at2->name);


	 // First test the NAG-ASN link (that order - as per dictionary)
	 //
	 if (name_1 == " C1 ")
	    if (name_2 == " ND2")
	       if (close[i].distance < smallest_link_dist) {
		  smallest_link_dist = close[i].distance;
		  link_type = "NAG-ASN";
	       }
      
	 if (name_1 == " O4 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "BETA1-4");
		  if (glyco_chiral_quad.chiral_volume() > 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "BETA1-4";
		  }
	       }
      
	 if (name_1 == " O2 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "BETA1-2");
		  if (glyco_chiral_quad.chiral_volume() > 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "BETA1-2";
		  }
	       }
      
	 if (name_1 == " O3 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "BETA1-3");
		  if (glyco_chiral_quad.chiral_volume() > 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "BETA1-3";
		  }
	       }
      
	       
	 // This should never happen :-)
	 // There are no biosynthetic pathways to make an BETA2-3 link for a SIA.
	 // (SIA BETA2-3 would be axial if it existed)
	 // 
	 if (name_1 == " C2 " )
	    if (name_2 == " O3 ")
	       if (std::string(close[i].at1->GetResName()) == "SIA") { 
		  if (close[i].distance < smallest_link_dist) {
		     coot::atom_quad glyco_chiral_quad(first, second, "BETA2-3");
		     std::cout << "   glyco_chiral BETA2-3 "
			       << close[i].at1->GetResName() << " "
			       << close[i].at2->GetResName() << " "
			       << glyco_chiral_quad.chiral_volume() << std::endl;
		     if (glyco_chiral_quad.chiral_volume() > 0.0) { 
			smallest_link_dist = close[i].distance;
			link_type = "BETA2-3";
		     }
		  }
	       }
	       
	 if (name_1 == " O6 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "BETA1-6");
		  if (glyco_chiral_quad.chiral_volume() > 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "BETA1-6";
		  }
	       }
	       
	 if (name_1 == " O2 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "ALPHA1-2");
		  if (glyco_chiral_quad.chiral_volume() < 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "ALPHA1-2";
		  }
	       }
	       
	 if (name_1 == " O3 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "ALPHA1-3");
		  if (glyco_chiral_quad.chiral_volume() < 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "ALPHA1-3";
		  }
	       }

	 if (name_1 == " C2 " )
	    if (name_2 == " O3 ")
	       if (std::string(close[i].at1->GetResName()) == "SIA") { 
		  if (close[i].distance < smallest_link_dist) {
		     coot::atom_quad glyco_chiral_quad(first, second, "ALPHA2-3");
		     std::cout << "   glyco_chiral ALPHA2-3 "
			       << close[i].at1->GetResName() << " "
			       << close[i].at2->GetResName() << " "
			       << glyco_chiral_quad.chiral_volume() << std::endl;
		     if (glyco_chiral_quad.chiral_volume() < 0.0) { 
			smallest_link_dist = close[i].distance;
			link_type = "ALPHA2-3";
		     }
		  }
	       }
      
	 if (name_1 == " O4 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "ALPHA1-4");
		  if (glyco_chiral_quad.chiral_volume() < 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "ALPHA1-4";
		  }
	       }
      
	 if (name_1 == " O6 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "ALPHA1-6");
		  if (glyco_chiral_quad.chiral_volume() < 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "ALPHA1-6";
		  }
	       }
      }
   }
   catch (std::runtime_error rte) {
      std::cout << "WARNING::" << rte.what() << std::endl;
   }

   if (0) 
      std::cout << "   debug:: find_glycosidic_linkage_type() for "
		<< first->GetChainID() << " " << first->GetSeqNum() << " " << first->GetInsCode()
		<< first->GetResName() << ","
		<< second->GetChainID() << " " << second->GetSeqNum() << " " << second->GetInsCode()
		<< second->GetResName() 
		<< " returns \"" << link_type << "\""
		<< std::endl;
   
   return link_type;
}

std::pair<std::string, bool>
coot::protein_geometry::find_glycosidic_linkage_type_with_order_switch(CResidue *first, CResidue *second) const {

   std::pair<std::string, bool> r("", false);

   std::string l = find_glycosidic_linkage_type(first, second);

   if (l == "") { 
      l = find_glycosidic_linkage_type(second,first);
      if (l != "") {
	 r.first = l;
	 r.second = true;
      } 
   } else {
      r.first = l;
      r.second = false;
   } 
   return r;
} 



      
void
coot::protein_geometry::set_verbose(bool verbose_flag) {
   verbose_mode = verbose_flag;
}

int 
coot::protein_geometry::init_standard() {

   // So, first we check if COOT_REFMAC_LIB_DIR has been set.  If it
   // try to use it.  If the directory fails to exist, try next...
   //
   // Next, try CLIBD_MON (which is e.g." ~/ccp4/ccp4-6.0/lib/data/monomers")
   //
   // Next, try to use the CCP4 lib/data dir.
   //
   // Then we check if the files are where should be given a proper
   // install (PKGDATADIR).  If so, use that...  [This is the prefered
   // option].
   //
   
   std::string dir = DATADIR; 
   std::string hardwired_default_place;
   hardwired_default_place = util::append_dir_dir(dir, "coot");
   hardwired_default_place = util::append_dir_dir(hardwired_default_place, "lib");
   bool using_clibd_mon = 0; 

   std::string mon_lib_dir; 
   short int env_dir_fails = 0;
   int istat;
   struct stat buf;
   char *cmld = NULL;
   
   char *s = getenv("COOT_REFMAC_LIB_DIR");
   if (s) {
      istat = stat(s, &buf);
      if (istat) { 
	 env_dir_fails = 1;
	 std::cout << "WARNING:: Coot REFMAC dictionary override COOT_REFMAC_LIB_DIR"
		   << "failed to find a dictionary " << s << std::endl;
      } else {
	 mon_lib_dir = s;
      }
   }

   if (!s || env_dir_fails) {
      cmld = getenv("COOT_MONOMER_LIB_DIR"); // for phenix.
      // we find $COOT_MONOMER_LIB_DIR/a/ALA.cif
      if (cmld) {
	 mon_lib_dir = cmld;
      }
   } 
      
   if (!s || env_dir_fails) {

      // OK, so COOT_REFMAC_LIB_DIR didn't provide a library.

      // Next try CLIBD_MON:
      s = getenv("CLIBD_MON");
      if (s) {
	 istat = stat(s, &buf);
	 if (istat) { 
	    env_dir_fails = 1;
	 } else {
	    env_dir_fails = 0;
	    std::cout << "INFO:: Using Standard CCP4 Refmac dictionary from"
		      << " CLIBD_MON: " << s << std::endl;
	    mon_lib_dir = s;
	    using_clibd_mon = 1;
	 }
      }

      
      if (!s || env_dir_fails) {
	 // Next, try CCP4_LIB

	 s = getenv("CCP4_LIB");
	 if (s) {
	    std::cout << "INFO:: Using Standard CCP4 Refmac dictionary: "
		      << s << std::endl;
	    mon_lib_dir = s;

	 } else {

	    // OK, CCP4 failed to give us a dictionary, now try the
	    // version that comes with Coot:

	    int istat = stat(hardwired_default_place.c_str(), &buf);
	    if (istat == 0) {
	       mon_lib_dir = hardwired_default_place;
	    } else {

	       // OK, let's look for $COOT_PREFIX/share/coot/lib (as you
	       // would with the binary distros)

	       s = getenv("COOT_PREFIX");
	       if (s) {
		  std::string lib_dir = util::append_dir_dir(s, "share");
		  lib_dir = util::append_dir_dir(lib_dir, "coot");
		  lib_dir = util::append_dir_dir(lib_dir, "lib");
		  istat = stat(lib_dir.c_str(), &buf);
		  if (istat == 0) {
		     mon_lib_dir = lib_dir;
		  } else {
		     std::cout << "WARNING:: COOT_PREFIX set, but no dictionary lib found\n";
		  }
	       } else {
		  std::cout << "WARNING:: COOT_PREFIX not set, all attemps to "
			    << "find dictionary lib failed\n";
	       }
	    }
	 }
      }
   }
   
   if (mon_lib_dir.length() > 0) {
      std::string filename = mon_lib_dir;
      // contains the linkages:
      filename += "/data/monomers/list/mon_lib_list.cif";
      if (using_clibd_mon) {
	 filename = mon_lib_dir;
	 filename += "list/mon_lib_list.cif";
      }
      // now check that that file is there:
      struct stat buf;
      int istat = stat(filename.c_str(), &buf);
      if ((istat == 0) && (! S_ISREG(buf.st_mode))) {
	 std::cout << "ERROR: dictionary " << filename
		   << " is not a regular file" << std::endl;
      } else {
	 // OK 
	 
      }
   } else { 
      std::cout << "WARNING: Failed to read restraints dictionary. "
		<< std::endl; 
   }

   // setting up CCP4 sets mon_lib_cif to
   // $CCP4_MASTER/lib/data/monomers (by using $CLIBD_MON).
   // 
   std::string mon_lib_cif = mon_lib_dir + "/data/monomers/list/mon_lib_list.cif";
   std::string energy_cif_file_name = mon_lib_dir + "/data/monomers/ener_lib.cif";
   if (using_clibd_mon) { 
      mon_lib_cif = mon_lib_dir + "/list/mon_lib_list.cif";
      energy_cif_file_name = mon_lib_dir + "/ener_lib.cif";
   }
   if (cmld) { 
      mon_lib_cif = cmld;
      mon_lib_cif += "/list/mon_lib_list.cif";
      energy_cif_file_name = std::string(cmld) + "/ener_lib.cif";
   }
   
   init_refmac_mon_lib(mon_lib_cif, coot::protein_geometry::MON_LIB_LIST_CIF);
   // now the protein monomers:
   read_number = 1;
   std::vector <std::string> protein_mono = standard_protein_monomer_files();
   for (unsigned int i=0; i<protein_mono.size(); i++) { 
      std::string monomer_cif_file = protein_mono[i];
      if (!cmld && !using_clibd_mon) {
	 monomer_cif_file = "data/monomers/" + monomer_cif_file;
      }
      refmac_monomer(mon_lib_dir, monomer_cif_file); // update read_number too :)
   }


   read_energy_lib(energy_cif_file_name);

   return read_number;
}


int 
coot::protein_geometry::refmac_monomer(const std::string &s, // dir
				       const std::string &protein_mono) { // extra path to file
   
   std::string filename = s; 
   filename = coot::util::append_dir_file(s, protein_mono);
   struct stat buf;
   int istat = stat(filename.c_str(), &buf);
   if (istat == 0) { 
      if (S_ISREG(buf.st_mode)) {
	 init_refmac_mon_lib(filename, read_number);
	 read_number++;
      } else {

#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
	 
	 // indenting goes strange if you factor out the else code here.
	 
	 if (! S_ISDIR(buf.st_mode) ) {
	    std::cout << "WARNING: " << filename 
		      << ": no such file (or directory)\n";
	 } else { 
	    std::cout << "ERROR: file " << filename
		      << " is not a regular file" << std::endl;
	 }
#else
	 if (! S_ISDIR(buf.st_mode) && 
	     ! S_ISLNK(buf.st_mode) ) {
	    std::cout << "WARNING: " << filename 
		      << ": no such file (or directory)\n";
	 } else { 
	    std::cout << "ERROR: file " << filename
		      << " is not a regular file" << std::endl;
	 }
#endif
      }
   }
   return read_number;
}

// Return 0 on failure to do a dynamic add (actually, the number of
// atoms read).
// 
int
coot::protein_geometry::try_dynamic_add(const std::string &resname, int read_number) {

   int success = 0;  // fail initially and is set to the number of
		     // atoms read from the mmcif dictionary file in
		     // init_refmac_mon_lib().

   // So what is happening here?
   //
   // It is a heirachy of setting
   //
   // The highest priority is COOT_REFMAC_LIB_DIR, if that is set we use it.
   // If COOT_REFMAC_LIB_DIR is then try CLIB (a CCP4 setting).
   // If that is not set, then we fall back to the default directory:
   // $prefix/share/coot onto which we tag a "lib" dir.
   // 

   char *s  = getenv("COOT_REFMAC_LIB_DIR");
   char *cmld = getenv("COOT_MONOMER_LIB_DIR");
   if (! s) {
      s  = getenv("CLIB");
      if (! s) {
	 std::string tmp_string(PKGDATADIR);
	 tmp_string = util::append_dir_dir(tmp_string, "lib");
	 s = new char[tmp_string.length() + 1];
	 strcpy(s, tmp_string.c_str());
      } else {
	 std::cout << "INFO:: using standard CCP4 Refmac dictionary"
		   << " to search for " << resname << std::endl; 
      } 
   }

   if (s) {
      std::string filename(s);
      std::string beta_anomer_name; 
      std::string alpha_anomer_name; 
      std::string alt_beta_anomer_name; 
      std::string alt_alpha_anomer_name; 
      if (cmld) { 
	 filename = cmld;
	 filename += "/";
      } else {
	 filename += "/data/monomers/";
      }
      if (resname.length() > 0) {
	 const char rs = resname[0];
	 const char v = tolower(rs); // get the sub directory name
         char v1[2];
	 v1[0] = v;
	 v1[1] = '\0';
	 std::string letter(v1);
	 filename += letter;
	 filename += "/";
	 std::string upcased_resname_filename = filename;
	 if (resname.length() > 2 ) { 
	    if (resname[2] != ' ') { 
	       filename += resname;
	       upcased_resname_filename += coot::util::upcase(resname);
	    } else {
	       filename += resname.substr(0,2);
	       upcased_resname_filename += coot::util::upcase(resname.substr(0,2));
	    }
	 } else {
	    filename += resname;
	    upcased_resname_filename += coot::util::upcase(resname);	    
	 }
	 beta_anomer_name = filename;
	 beta_anomer_name += "-b-D.cif";
	 alt_beta_anomer_name = filename;
	 alt_beta_anomer_name += "-B-D.cif";
	 alpha_anomer_name  = filename;
	 alpha_anomer_name += "-a-L.cif";
	 alt_alpha_anomer_name  = filename;
	 alt_alpha_anomer_name += "-A-L.cif";
	 filename += ".cif";
	 upcased_resname_filename += ".cif";
	 
	 struct stat buf;
	 int istat = stat(filename.c_str(), &buf);
	 if (istat == 0) { 
	    if (S_ISREG(buf.st_mode)) {
	       success = init_refmac_mon_lib(filename, read_number);
	    } else {
	       
	       // continue with regular file code
#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
	       if (! S_ISDIR(buf.st_mode) ) {
		  std::cout << "WARNING: " << filename 
			    << ": no such file (or directory)\n";
	       } else { 
		  std::cout << "ERROR: dictionary " << filename
			    << " is not a regular file" << std::endl;
	       }
#else
	       if (! S_ISDIR(buf.st_mode) && 
		   ! S_ISLNK(buf.st_mode) ) {
		  std::cout << "WARNING: " << filename 
			    << ": no such file (or directory)\n";
	       } else { 
		  std::cout << "ERROR: dictionary " << filename
			    << " is not a regular file" << std::endl;
	       }
#endif
	    }
	 } else { 
	    
	    // Regular file doesn't exist,
	    
	    // Try the upcased filename
	    // 
	    istat = stat(upcased_resname_filename.c_str(), &buf);
	    if (istat == 0) {
		  success = init_refmac_mon_lib(upcased_resname_filename, read_number);
	    } else { 

	       // try the beta anomer version
	       istat = stat(beta_anomer_name.c_str(), &buf);
	       if (istat == 0) {
		  success = init_refmac_mon_lib(beta_anomer_name, read_number);
	       } else {
		  // try the upcased file name e.g. xxx/NAG-B-D.cif
		  istat = stat(alt_beta_anomer_name.c_str(), &buf);
		  if (istat == 0) {
		     success = init_refmac_mon_lib(alt_beta_anomer_name, read_number);
		  } else {
		     // alpha?
		     istat = stat(alpha_anomer_name.c_str(), &buf);
		     if (istat == 0) {
			success = init_refmac_mon_lib(alpha_anomer_name, read_number);
		     } else {
			istat = stat(alt_alpha_anomer_name.c_str(), &buf);
			if (istat == 0) {
			   success = init_refmac_mon_lib(alt_alpha_anomer_name, read_number);
			}
		     }
		  }
	       }
	    }
	 }
      } 
   }

   // now, did we get something with minimal description? If so,
   // delete it, it was a fail.
   // 
   std::pair<bool, dictionary_residue_restraints_t> p = get_monomer_restraints(resname);
   if (! p.first) {
      success = 0;
   } else {
      // elide minimal description restraints.
      if (p.second.residue_info.description_level == "M") { 
	 success = 0;
	 delete_mon_lib(resname);
      }
   }
   return success;
} 

std::vector <coot::dict_torsion_restraint_t>
coot::protein_geometry::get_monomer_torsions_from_geometry(const std::string &monomer_type) {

   short int ifound = 0;
   std::vector <coot::dict_torsion_restraint_t> rv;
   
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].residue_info.comp_id == monomer_type) {
	 ifound = 1;
// 	 std::cout << "DEBUG:: found " << monomer_type << " in restraints dictionary" 
// 		   << std::endl;
	 return dict_res_restraints[i].torsion_restraint;
      }
   }

   // check synonyms before 3-letter codes.
   // 
   if (ifound == 0) { 
      // OK, that failed to, perhaps there is a synonym?
      for (unsigned int i=0; i<residue_name_synonyms.size(); i++) { 
	 if (residue_name_synonyms[i].comp_alternative_id == monomer_type) {
	    int ndict = dict_res_restraints.size();
	    for (int i=0; i<ndict; i++) {
	       if (dict_res_restraints[i].residue_info.comp_id == residue_name_synonyms[i].comp_id) {
		  ifound = 1;
		  rv = dict_res_restraints[i].torsion_restraint;
		  break;
	       }
	    }
	 }
	 if (ifound)
	    break;
      }
   }

   if (ifound == 0) {
      int read_number = 40; // just a filler, FIXME.
      int nbonds = try_dynamic_add(monomer_type, read_number);
      if (nbonds > 0) {
	 // OK, we got it, what are the torsions?
	 for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
	    if (dict_res_restraints[i].residue_info.comp_id == monomer_type) {
	       ifound = 1;
	       rv = dict_res_restraints[i].torsion_restraint;
	       break;
	    }
	 }
      }
   }

   

   if (ifound == 0) { // which it should be if we get here with that
		      // return in the loop
      std::cout << "WARNING: residue type " << monomer_type << " not found "
		<< "in restraints dictionary (torsion)" << std::endl;
   }
   rv = filter_torsion_restraints(rv);
   return rv;
}  

std::vector <coot::dict_torsion_restraint_t>
coot::protein_geometry::get_monomer_torsions_from_geometry(const std::string &monomer_type,
							   bool find_hydrogen_torsions_flag) const {

   bool ifound = 0;
   int ii = -1; // unset, set when torsion found, ii is the index in
		// dict_res_restraints for the torsion restraints,
		// needed so that we can ask about hydrogens.
   std::vector <coot::dict_torsion_restraint_t> unfiltered_torsion_restraints;
   std::vector <coot::dict_torsion_restraint_t> filtered_torsion_restraints;
   
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].residue_info.comp_id == monomer_type) {
	 ifound = 1;
	 unfiltered_torsion_restraints = dict_res_restraints[i].torsion_restraint;
	 ii = i; // used for filtering out hydrogens 
      }
   }

   if (ifound == 0) { // which it should be if we get here with that
		      // return in the loop
      // try dynamic add?
      std::cout << "WARNING: residue type " << monomer_type << " not found "
		<< "in restraints dictionary (in get_monomer_torsions_from_geometry(mon, hy)" << std::endl;
   } else {
      if (find_hydrogen_torsions_flag) {
	 filtered_torsion_restraints = unfiltered_torsion_restraints;
      } else {
	 // we don't want torsions that move Hydrogens
	 int nt = dict_res_restraints[ii].torsion_restraint.size();
	 for (int it=0; it<nt; it++) { 
	    if (!dict_res_restraints[ii].is_hydrogen(dict_res_restraints[ii].torsion_restraint[it].atom_id_1())) {
	       if (!dict_res_restraints[ii].is_hydrogen(dict_res_restraints[ii].torsion_restraint[it].atom_id_4())) {
		  filtered_torsion_restraints.push_back(dict_res_restraints[ii].torsion_restraint[it]);
	       }
	    }
	 }
      }
   }

   // more filtering (only one version of a torsion_restraint that
   // have the same middle atoms).
   // 
   filtered_torsion_restraints = filter_torsion_restraints(filtered_torsion_restraints);
   return filtered_torsion_restraints;
}

// Return only one version of a torsions restraint that have the same
// middle atoms.
// 
std::vector <coot::dict_torsion_restraint_t>
coot::protein_geometry::filter_torsion_restraints(const std::vector <coot::dict_torsion_restraint_t> &restraints_in) const {

   std::vector <coot::dict_torsion_restraint_t> r;

   for (unsigned int i=0; i<restraints_in.size(); i++) {
      std::string a2 = restraints_in[i].atom_id_2_4c();
      std::string a3 = restraints_in[i].atom_id_3_4c();
      bool match = 0;
      for (unsigned int j=0; j<r.size(); j++) {
	 if (r[j].atom_id_2_4c() == a2)
	    if (r[j].atom_id_3_4c() == a3)
	       match = 1;
      }
      if (match == 0)
	 r.push_back(restraints_in[i]);
   }

   std::sort(r.begin(), r.end(), torsion_restraints_comparer);
   return r;
}

// static
bool
coot::protein_geometry::torsion_restraints_comparer(const coot::dict_torsion_restraint_t &a, const coot::dict_torsion_restraint_t &b) {
   
      std::string a2 = a.atom_id_2_4c();
      std::string a3 = a.atom_id_3_4c();
      std::string b2 = b.atom_id_2_4c();
      std::string b3 = b.atom_id_3_4c();

      if (a2 < b2)
	 return 0;
      else
	 if (a2 > b2)
	    return 1;
	 else
	    if (a3 < b3)
	       return 0;
      
      return 1;
}


std::vector <coot::dict_chiral_restraint_t>
coot::protein_geometry::get_monomer_chiral_volumes(const std::string monomer_type) const { 

   bool ifound = 0;
   std::vector <coot::dict_chiral_restraint_t> rv;
   
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].residue_info.comp_id == monomer_type) {
	 ifound = 1;
	 rv = dict_res_restraints[i].chiral_restraint;
	 break;
      }
   }
   // OK so the monomer_type did not match the comp_id.  Perhaps the
   // comp_id was not the same as the three letter code, so let's
   // check the monomer_type against the three_letter_codes.
   // 
   if (ifound == 0) {
      for (int i=0; i<dict_res_restraints.size(); i++) {
	 if (dict_res_restraints[i].residue_info.three_letter_code == monomer_type) {
	    ifound = 1;
	    rv = dict_res_restraints[i].chiral_restraint;
	    break;
	 }
      }
   } 
   if (ifound == 0) {
      // try dynamic add?
      std::cout << "WARNING: residue type " << monomer_type << " not found "
		<< "in restraints dictionary (chiral)" << std::endl;
   }
   return rv;
}



std::vector <std::string>
coot::protein_geometry::standard_protein_monomer_files() const {

   std::vector <std::string> s;

   s.push_back("a/ALA.cif");
   s.push_back("a/ASP.cif");
   s.push_back("a/ASN.cif");
   s.push_back("c/CYS.cif");
   s.push_back("g/GLN.cif");
   s.push_back("g/GLY.cif");
   s.push_back("g/GLU.cif");
   s.push_back("p/PHE.cif");
   s.push_back("h/HIS.cif");
   s.push_back("i/ILE.cif");
   s.push_back("l/LYS.cif");
   s.push_back("l/LEU.cif");
   s.push_back("m/MET.cif");
   s.push_back("m/MSE.cif");
   s.push_back("p/PRO.cif");
   s.push_back("a/ARG.cif");
   s.push_back("s/SER.cif");
   s.push_back("t/THR.cif");
   s.push_back("v/VAL.cif");
   s.push_back("t/TRP.cif");
   s.push_back("t/TYR.cif");

   s.push_back("p/PO4.cif");
   s.push_back("s/SO4.cif");
   s.push_back("g/GOL.cif");
   s.push_back("e/ETH.cif");
   s.push_back("c/CIT.cif");

   s.push_back("a/AR.cif");
   s.push_back("a/AD.cif");
   s.push_back("c/CR.cif");
   s.push_back("c/CD.cif");
   s.push_back("g/GR.cif");
   s.push_back("g/GD.cif");
   s.push_back("t/TD.cif");
   s.push_back("u/UR.cif");

   // new-style (CCP4 6.2) RNA names
   s.push_back("a/A.cif");
   s.push_back("c/C.cif");
   s.push_back("g/G.cif");
   s.push_back("u/U.cif");

   // new-style (CCP4 6.2) DNA names
   s.push_back("d/DA.cif");
   s.push_back("d/DC.cif");
   s.push_back("d/DG.cif");
   s.push_back("d/DT.cif");
   
   s.push_back("h/HOH.cif");

   return s;
}


bool
coot::protein_geometry::have_dictionary_for_residue_type(const std::string &monomer_type,
							 int read_number_in) { 

   bool ifound = 0;
   int ndict = dict_res_restraints.size();
   read_number = read_number_in;
   for (int i=0; i<ndict; i++) {
      if (dict_res_restraints[i].residue_info.comp_id == monomer_type) {
	 if (! dict_res_restraints[i].is_from_sbase_data()) {
	    ifound = 1;
	    break;
	 }
      }
   }

   // check synonyms before checking three-letter-codes
   
   if (ifound == 0) { 
      // OK, that failed to, perhaps there is a synonym?
      for (unsigned int i=0; i<residue_name_synonyms.size(); i++) { 
	 if (residue_name_synonyms[i].comp_alternative_id == monomer_type) {
	    for (int j=0; j<ndict; j++) {
	       if (dict_res_restraints[j].residue_info.comp_id == residue_name_synonyms[i].comp_id) {
		  ifound = 1;
		  break;
	       }
	    }
	 }
	 if (ifound)
	    break;
      }
   }
	 
   // OK so the monomer_type did not match the comp_id.  Perhaps the
   // comp_id was not the same as the three letter code, so let's
   // check the monomer_type against the three_letter_codes.
   // 
   if (ifound == 0) {
      for (int i=0; i<ndict; i++) {
	 if (dict_res_restraints[i].residue_info.three_letter_code == monomer_type) {
	    if (! dict_res_restraints[i].is_from_sbase_data()) {
	       ifound = 1;
	       break;
	    }
	 }
      }
   }

   if (ifound == 0) {
      ifound = try_dynamic_add(monomer_type, read_number);
   }
   return ifound;
}

// this is const because there is no dynamic add
bool
coot::protein_geometry::have_dictionary_for_residue_type_no_dynamic_add(const std::string &monomer_type) const {

   bool ifound = 0;
   int ndict = dict_res_restraints.size();
   for (int i=0; i<ndict; i++) {
      if (dict_res_restraints[i].residue_info.comp_id == monomer_type) {
	 if (! dict_res_restraints[i].is_from_sbase_data()) {
	    ifound = 1;
	    break;
	 }
      }
   }
   return ifound;
} 


bool
coot::protein_geometry::have_dictionary_for_residue_types(const std::vector<std::string> &residue_types) {

   bool have_all = 1;
   int read_number = 30; // hack dummy thing.
   for (unsigned int i=0; i<residue_types.size(); i++) {
      int ifound = have_dictionary_for_residue_type(residue_types[i], read_number);
      if (ifound == 0) {
	 have_all = 0;
      } 
      read_number++;
   }
   return have_all;
}

// this is const because there is no dynamic add
bool
coot::protein_geometry::have_at_least_minimal_dictionary_for_residue_type(const std::string &monomer_type) const {

   bool ifound = 0;
   int ndict = dict_res_restraints.size();
   for (int i=0; i<ndict; i++) {
      if (dict_res_restraints[i].residue_info.comp_id == monomer_type) {
	 ifound = 1;
	 break;
      }
   }
   return ifound;
} 



// Check that the atom names in the residue match the atom names in
// the dictionary.  Atom " OXT" is treated as a special case (it does
// not cause a failure when " OXT" is not in the dictionary for the
// residue).
//
// Similarly, we can not return problematic status (0) if the
// hydrogens do not match, if caller wishes.
//
// This does only monomer by monomer testing.
//
// There is no test of DEL-O1 (for example) atoms when making links
// between monomers.
// 
// Return in pair.first the state of the match and in .second, the
// list of atoms that do not match the dictionary.
// 
std::pair<bool, std::vector<std::string> >
coot::protein_geometry::atoms_match_dictionary(CResidue *residue_p,
					       bool check_hydrogens_too_flag,
					       const coot::dictionary_residue_restraints_t &restraints) {

   std::vector<std::string> atom_name_vec;
   bool status = 1;

   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

   bool debug = 0;
   if (debug) {
      std::cout << "=== atoms_match_dictionary() with these residue atom names ======= " << std::endl;
      for (int i=0; i<n_residue_atoms; i++) {
	 std::cout << i << "  :" << residue_atoms[i]->name << ":" << std::endl;
      } 
      std::cout << "=== atoms_match_dictionary() with these residue atom names ======= " << std::endl;
      for (unsigned int irat=0; irat<restraints.atom_info.size(); irat++) {
	 std::cout << irat << "  :" << restraints.atom_info[irat].atom_id_4c
		   << ":" << std::endl;
      } 
      

   } 

   
   for (unsigned int i=0; i<n_residue_atoms; i++) {

      if (! residue_atoms[i]->isTer()) { 
	 std::string residue_atom_name(residue_atoms[i]->name);
	 std::string ele(residue_atoms[i]->element);

	 bool found = 0;
	 if (ele == " H")
	    if (check_hydrogens_too_flag == 0)
	       found = 1;

	 if (! found) { 
	    for (unsigned int irestraint_atom_name=0; irestraint_atom_name<restraints.atom_info.size(); irestraint_atom_name++) {
	       if (restraints.atom_info[irestraint_atom_name].atom_id_4c == residue_atom_name) {
		  found = 1;
		  break;
	       }
	    }
	 }
	 if (! found) {
	    if (residue_atom_name != " OXT") { 
	       atom_name_vec.push_back(residue_atom_name);
	       status = 0;
	    }
	 }
      }
   }
   
   return std::pair<bool, std::vector<std::string> > (status, atom_name_vec);
}



// return a pair, overall status, and pair of residue names and
// atom names that dont't match.
//
std::pair<bool, std::vector<std::pair<std::string, std::vector<std::string> > > >
coot::protein_geometry::atoms_match_dictionary(const std::vector<CResidue *> residues,
		       bool check_hydrogens_too_flag) {

   bool status = 1;
   std::vector<std::pair<std::string, std::vector<std::string> > > p;
   
   
   for (unsigned int ires=0; ires<residues.size(); ires++) { 
      std::string res_name(residues[ires]->GetResName());
      std::pair<bool, coot::dictionary_residue_restraints_t> restraints =
	 get_monomer_restraints(res_name);
      if (restraints.first) { 
	 std::pair<bool, std::vector<std::string> > r =
	    atoms_match_dictionary(residues[ires], check_hydrogens_too_flag,
				   restraints.second);
	 if (r.first == 0) {
	    std::pair<std::string, std::vector<std::string> > p_bad(res_name, r.second);
	    p.push_back(p_bad);
	    status = 0;
	 }
      }
   }

   return std::pair<bool, std::vector<std::pair<std::string, std::vector<std::string> > > > (status, p);
}


// return a pair, the first is status (1 if the name was found, 0 if not)
// 
std::pair<bool, std::string>
coot::protein_geometry::get_monomer_name(const std::string &comp_id) const {

   std::pair<bool, std::string> r(0,"");

   std::map<std::string,coot::dictionary_residue_restraints_t>::const_iterator it = 
      simple_monomer_descriptions.find(comp_id);

   if (it != simple_monomer_descriptions.end()) {
      r.first = 1;
      std::string s = it->second.residue_info.name;
      r.second = coot::util::remove_trailing_whitespace(s);
   } 

   return r;
} 




// Try comparing vs the comp_id first, if that fails compare the
// three_letter_code to the monomer_type.
//
// In future, try to come here only with the monomer_type adjusted to
// the comp_id, for example, monomer_type should be "NAG-b-D", not
// "NAG".
// 
std::pair<bool, coot::dictionary_residue_restraints_t>
coot::protein_geometry::get_monomer_restraints(const std::string &monomer_type) const {

   return get_monomer_restraints_internal(monomer_type, 0);
   
}

std::pair<bool, coot::dictionary_residue_restraints_t>
coot::protein_geometry::get_monomer_restraints_at_least_minimal(const std::string &monomer_type) const {

   return get_monomer_restraints_internal(monomer_type, 1);
}

std::pair<bool, coot::dictionary_residue_restraints_t>
coot::protein_geometry::get_monomer_restraints_internal(const std::string &monomer_type, bool allow_minimal_flag) const {

   coot::dictionary_residue_restraints_t t(std::string("(null)"), 0);
   std::pair<bool, coot::dictionary_residue_restraints_t> r(0,t);

   unsigned int nrest = dict_res_restraints.size();
   for (unsigned int i=0; i<nrest; i++) {
      if (dict_res_restraints[i].residue_info.comp_id  == monomer_type) {
	 if ((allow_minimal_flag == 1) || (! dict_res_restraints[i].is_from_sbase_data())) { 
	    r.second = dict_res_restraints[i];
	    r.first = 1;
	    break;
	 }
      }
   }

   if (!r.first) {
      // OK, that failed to, perhaps there is a synonym?
      for (unsigned int i=0; i<residue_name_synonyms.size(); i++) { 
	 if (residue_name_synonyms[i].comp_alternative_id == monomer_type) {
	    int ndict = dict_res_restraints.size();
	    for (int j=0; j<ndict; j++) {
	       if (dict_res_restraints[j].residue_info.comp_id == residue_name_synonyms[i].comp_id) {
		  r.first = 1;
		  r.second = dict_res_restraints[j];
		  break;
	       }
	    }
	 }
	 if (r.first)
	    break;
      }
   }
   
   if (!r.first) {
      for (unsigned int i=0; i<nrest; i++) {
	 if (dict_res_restraints[i].residue_info.three_letter_code  == monomer_type) {
	    if ((allow_minimal_flag == 1) || (! dict_res_restraints[i].is_from_sbase_data())) { 
	       r.second = dict_res_restraints[i];
	       r.first = 1;
	       break;
	    }
	 }
      }
   }
   return r;
}

// return -1 on monomer not found.
int
coot::protein_geometry::get_monomer_restraints_index(const std::string &monomer_type, bool allow_minimal_flag) const {

   int r = -1;

   unsigned int nrest = dict_res_restraints.size();
   for (unsigned int i=0; i<nrest; i++) {
      if (dict_res_restraints[i].residue_info.comp_id  == monomer_type) {
	 if ((allow_minimal_flag == 1) || (! dict_res_restraints[i].is_from_sbase_data())) { 
	    r = i;
	    break;
	 }
      }
   }

   if (r == -1) {
      // OK, that failed to, perhaps there is a synonym?
      for (unsigned int i=0; i<residue_name_synonyms.size(); i++) { 
	 if (residue_name_synonyms[i].comp_alternative_id == monomer_type) {
	    int ndict = dict_res_restraints.size();
	    for (int j=0; j<ndict; j++) {
	       if (dict_res_restraints[j].residue_info.comp_id == residue_name_synonyms[i].comp_id) {
		  r = 1;
		  break;
	       }
	    }
	 }
	 if (r != -1)
	    break;
      }
   }
   
   if (r == -1) {
      for (unsigned int i=0; i<nrest; i++) {
	 if (dict_res_restraints[i].residue_info.three_letter_code  == monomer_type) {
	    if ((allow_minimal_flag == 1) || (! dict_res_restraints[i].is_from_sbase_data())) { 
	       r = i;
	       break;
	    }
	 }
      }
   }
   
   return r;
}


std::string
coot::protein_geometry::get_type_energy(const std::string &atom_name,
					const std::string &residue_name) const { 
   // return "" if not found, else return the energy type found in ener_lib.cif
   //
   std::string r;
   int indx = get_monomer_restraints_index(residue_name, 1);
   if (indx != -1) {
      const coot::dictionary_residue_restraints_t &restraints = dict_res_restraints[indx];
      r = restraints.type_energy(atom_name);
   } 
   return r;
   
}

// return -1.1 on failure to look up.
// 
double
coot::protein_geometry::get_vdw_radius(const std::string &atom_name,
				       const std::string &residue_name,
				       bool use_vdwH_flag) const {
   double r = -1.1;
   int indx = get_monomer_restraints_index(residue_name, 1);
   if (indx != -1) {
      const coot::dictionary_residue_restraints_t &restraints = dict_res_restraints[indx];
      std::string et = restraints.type_energy(atom_name);
      if (et != "") {
	 std::map<std::string, energy_lib_atom>::const_iterator it =
	    energy_lib.atom_map.find(et);
	 if (it != energy_lib.atom_map.end()) {
	    if (use_vdwH_flag)
	       r = it->second.vdwh_radius;
	    else
	       r = it->second.vdw_radius;
	 }
      }
   } else {
      std::cout << "  no restraints for type " << residue_name << std::endl;
   }
   return r;
}

// expand to 4c, the atom_id, give that it should match an atom of the type comp_id.
// Used in chem mods, when we don't know the comp_id until residue modification time.
// 
std::string
coot::protein_geometry::atom_id_expand(const std::string &atom_id,
				       const std::string &comp_id) const {

   std::string s = atom_id;
   int idx = get_monomer_restraints_index(comp_id, 1);
   if (idx != -1) {
      const coot::dictionary_residue_restraints_t &restraints =
	 dict_res_restraints[idx];
      const std::vector<dict_atom> &atoms = restraints.atom_info;
      for (unsigned int iat=0; iat<atoms.size(); iat++) { 
	 if (atoms[iat].atom_id == atom_id) { 
	    s = atoms[iat].atom_id_4c;
	    break;
	 }
      }
   }
//    std::cout << "atom_id_expand() \"" << atom_id << "\" expanded to \""
// 	     << s << "\"" << std::endl;
   return s; 
} 






// Hmmm... empty function, needs examining.
// 
// Use dynamic add if necessary.
// 
// Return -1 if the thing was not found or added.
int
coot::protein_geometry::get_monomer_type_index(const std::string &monomer_type) { 

   int i = -1;

   return i;
}

// Return 1 for hydrogen or deuterium, 0 for not found or not a hydrogen.
//
bool
coot::dictionary_residue_restraints_t::is_hydrogen(const std::string &atom_name) const {

   bool r = 0;
   for (unsigned int i=0; i<atom_info.size(); i++) {
      if (atom_info[i].atom_id_4c == atom_name) {
	 if (atom_info[i].type_symbol == "H" ||
	     atom_info[i].type_symbol == " H" ||
	     atom_info[i].type_symbol == "D") {
	    r = 1;
	    break;
	 }
      }
   }
   return r;
}

bool
coot::dictionary_residue_restraints_t::has_unassigned_chiral_volumes() const {
   bool r = 0;
   for (unsigned int ic=0; ic<chiral_restraint.size(); ic++) {
      if (chiral_restraint[ic].has_unassigned_chiral_volume()) {
	 r = 1;
	 break;
      }
   }
   return r;
}

bool
coot::dictionary_residue_link_restraints_t::has_unassigned_chiral_volumes() const {
   bool r = 0;
   for (unsigned int ic=0; ic<link_chiral_restraint.size(); ic++) {
      if (link_chiral_restraint[ic].has_unassigned_chiral_volume()) {
	 r = 1;
	 break;
      }
   }
   return r;
}


int
coot::dictionary_residue_restraints_t::assign_chiral_volume_targets() {

   int ich = 0;

   if (0) 
      std::cout << "DEBUG:: in dictionary_residue_restraints_t::assign_chiral_volume_targets "
		<< "there are " << chiral_restraint.size() << " chiral restraints for "
		<< residue_info.comp_id << " \n";
    
   for (unsigned int i=0; i<chiral_restraint.size(); i++) {
      chiral_restraint[i].assign_chiral_volume_target(bond_restraint, angle_restraint);
      ich++;
   }
   return ich;
}

int
coot::dictionary_residue_link_restraints_t::assign_link_chiral_volume_targets() {

   int ic = 0;
   for (unsigned int i=0; i<link_chiral_restraint.size(); i++) {
      std::vector <coot::dict_bond_restraint_t> bond_restraints_1;
      std::vector <coot::dict_bond_restraint_t> bond_restraints_2;
      std::vector <coot::dict_angle_restraint_t> angle_restraints_1;
      std::vector <coot::dict_angle_restraint_t> angle_restraints_2;
      std::vector <coot::dict_link_bond_restraint_t> link_bonds;
      std::vector <coot::dict_link_angle_restraint_t> link_angles;
      
      link_chiral_restraint[i].assign_chiral_volume_target(bond_restraints_1,
							   angle_restraints_1,
							   bond_restraints_2,
							   angle_restraints_2,
							   link_bonds,
							   link_angles);
      ic++;
   }
   return ic;
}

double
coot::dict_chiral_restraint_t::assign_chiral_volume_target(const std::vector <dict_bond_restraint_t> &bonds,
							   const std::vector <dict_angle_restraint_t> &angles) {

   double vol = -1;
   double a = -1, b = -1, c = -1;
   double alpha = -1, beta = -1, gamma = -1;
   std::string mmdb_centre_atom =  atom_id_mmdb_expand(local_atom_id_centre);
   std::string mmdb_local_atom_id_1 = atom_id_mmdb_expand(local_atom_id_1);
   std::string mmdb_local_atom_id_2 = atom_id_mmdb_expand(local_atom_id_2);
   std::string mmdb_local_atom_id_3 = atom_id_mmdb_expand(local_atom_id_3);
   
   // local_atom_id_centre to local_atom_id_1 bond length
   for (unsigned int i=0; i<bonds.size(); i++) {
      try { 
	 if (bonds[i].atom_id_1_4c() == mmdb_centre_atom) {
	    if (bonds[i].atom_id_2_4c() == atom_id_mmdb_expand(local_atom_id_1)) { 
	       a = bonds[i].value_dist();
	    }
	    if (bonds[i].atom_id_2_4c() == atom_id_mmdb_expand(local_atom_id_2)) { 
	       b = bonds[i].value_dist();
	    }
	    if (bonds[i].atom_id_2_4c() == atom_id_mmdb_expand(local_atom_id_3)) { 
	       c = bonds[i].value_dist();
	    }
	 }
	 if (bonds[i].atom_id_2_4c() == atom_id_mmdb_expand(local_atom_id_centre)) {
	    if (bonds[i].atom_id_1_4c() == atom_id_mmdb_expand(local_atom_id_1)) { 
	       a = bonds[i].value_dist();
	    }
	    if (bonds[i].atom_id_1_4c() == atom_id_mmdb_expand(local_atom_id_2)) { 
	       b = bonds[i].value_dist();
	    }
	    if (bonds[i].atom_id_1_4c() == atom_id_mmdb_expand(local_atom_id_3)) { 
	       c = bonds[i].value_dist();
	    }
	 }
      }
      catch (std::runtime_error rte) {
	 // do nothing, it's not really an error if the dictionary
	 // doesn't have target geometry (the bonding description came
	 // from a Chemical Component Dictionary entry for example).
      } 
   }

   for (unsigned int i=0; i<angles.size(); i++) {
      if (angles[i].atom_id_2_4c() == mmdb_centre_atom) {
	 if ((angles[i].atom_id_1_4c() == mmdb_local_atom_id_2 &&
	      angles[i].atom_id_3_4c() == mmdb_local_atom_id_3) ||
	     (angles[i].atom_id_3_4c() == mmdb_local_atom_id_2 &&
	      angles[i].atom_id_1_4c() == mmdb_local_atom_id_3))  {
	    alpha = clipper::Util::d2rad(angles[i].angle());
	 }
	 if ((angles[i].atom_id_1_4c() == mmdb_local_atom_id_1 &&
	      angles[i].atom_id_3_4c() == mmdb_local_atom_id_3) ||
	     (angles[i].atom_id_3_4c() == mmdb_local_atom_id_1 &&
	      angles[i].atom_id_1_4c() == mmdb_local_atom_id_3))  {
	    beta = clipper::Util::d2rad(angles[i].angle());
	 }
	 if ((angles[i].atom_id_1_4c() == mmdb_local_atom_id_1 &&
	      angles[i].atom_id_3_4c() == mmdb_local_atom_id_2) ||
	     (angles[i].atom_id_3_4c() == mmdb_local_atom_id_1 &&
	      angles[i].atom_id_1_4c() == mmdb_local_atom_id_2))  {
	    gamma = clipper::Util::d2rad(angles[i].angle());
	 }
      }
   }

   
   if (a > 0 && b > 0 && c > 0) {
//       std::cout << "DEBUG:: found all distances in chiral restraint" << std::endl;
      if (alpha > 0 && beta > 0 && gamma > 0) {
// 	 std::cout << "DEBUG:: found all angles in chiral restraint" << std::endl;
	 vol = assign_chiral_volume_target_internal(a, b, c, alpha, beta, gamma);
      } else {
// 	 std::cout << "DEBUG:: failed to find all angles in chiral restraint"
// 		   << alpha << " " << beta << " " << gamma << std::endl;
      }
   } else {
//       std::cout << "DEBUG:: failed to find all distances in chiral restraint"
// 		<< a << " " << b << " " << c << std::endl;
   }
   return vol;
}

double
coot::dict_link_chiral_restraint_t::assign_chiral_volume_target(const std::vector <dict_bond_restraint_t> &bonds_1,
								const std::vector <dict_angle_restraint_t> &angles_1,
								const std::vector <dict_bond_restraint_t> &bonds_2,
								const std::vector <dict_angle_restraint_t> &angles_2,
								const std::vector <dict_link_bond_restraint_t> &link_bonds,
								const std::vector <dict_link_angle_restraint_t> &link_angles) {

   double d = 0;

   return d;
} 


// angles in radians.
double
coot::dict_chiral_restraint_t::assign_chiral_volume_target_internal(double a, double b, double c,
								    double alpha, double beta, double gamma) {

   // set target_volume
   // from: abc ( 1 - cos^2(alpha) - cos^2(beta) - cos^2(gamma) + 2(cos(alpha) + cos(beta) + cos(gamma)))^0.5

   double cos_alpha = cos(alpha);
   double cos_beta  = cos(beta);
   double cos_gamma = cos(gamma);
   
   double cos_2_alpha = cos_alpha * cos_alpha;
   double cos_2_beta  = cos_beta  * cos_beta;
   double cos_2_gamma = cos_gamma * cos_gamma;

//    std::cout << "input a, b, c, alpha, beta, gamma " << a << " "
// 	     << b << " " << c << " "
// 	     << clipper::Util::rad2d(alpha) << " "
// 	     << clipper::Util::rad2d(beta) << " "
// 	     << clipper::Util::rad2d(gamma) << " " << std::endl;

//    double a_bit = 1-cos_2_alpha-cos_2_beta-cos_2_gamma;
//    double b_bit = 2 * cos_alpha * cos_beta * cos_gamma;
//    double c_bit = a_bit + b_bit;

//    std::cout << "bits: " << a_bit << " " << b_bit << " " << c_bit << std::endl;

   target_volume_ = volume_sign * a*b*c*sqrt(1-cos_2_alpha-cos_2_beta-cos_2_gamma + 2*cos_alpha*cos_beta*cos_gamma);

   volume_sigma_ = 0.2;  // test value

   if (0)
      std::cout << "DEBUG:: assign_chiral_volume_target_internal() target_volume chiral: "
		<< target_volume_ << std::endl;
   
   return target_volume_;
}


std::string
coot::protein_geometry::three_letter_code(const unsigned int &i) const {

   std::string r = dict_res_restraints[i].residue_info.three_letter_code;
   if (r == "")
      r = dict_res_restraints[i].residue_info.comp_id;
   return r;
}

// add "synthetic" 5 atom planar peptide restraint
void
coot::protein_geometry::add_planar_peptide_restraint() {

   std::string link_id = "TRANS";
   std::string plane_id = "plane3";
   realtype dist_esd = 0.05;

   std::string atom_id; 
   std::vector<std::pair<int, std::string> > v;
   v.push_back(std::pair<int, std::string> (1, "CA"));
   v.push_back(std::pair<int, std::string> (1, "C"));
   v.push_back(std::pair<int, std::string> (1, "O"));
   v.push_back(std::pair<int, std::string> (2, "N"));
   v.push_back(std::pair<int, std::string> (2, "CA"));

   for (unsigned int i=0; i<v.size(); i++) 
      // link_add_plane(link_id, atom_id,     plane_id, atom_comp_id, dist_esd); 
         link_add_plane(link_id, v[i].second, plane_id, v[i].first,   dist_esd); 
}


void
coot::protein_geometry::remove_planar_peptide_restraint() {

   std::string link_id = "TRANS";
   std::string plane_id = "plane3";
   bool ifound = 0;

   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {
      if (dict_link_res_restraints[i].link_id == link_id) { // e.g "TRANS"

	 std::vector<coot::dict_link_plane_restraint_t>::iterator it;
	 for (it = dict_link_res_restraints[i].link_plane_restraint.begin();
	      it != dict_link_res_restraints[i].link_plane_restraint.end(); it++) {
	    if (it->plane_id == plane_id) {
	       ifound = 1;
	       // let's remove it
	       std::cout << "INFO:: before removal of plane3 TRANS has " 
			 << dict_link_res_restraints[i].link_plane_restraint.size()
			 << " plane restraints\n";
 	       dict_link_res_restraints[i].link_plane_restraint.erase(it);
	       std::cout << "INFO::  after removal of plane3 TRANS has " 
			 << dict_link_res_restraints[i].link_plane_restraint.size()
			 << " plane restraints\n";
	       break;
	    }
	 }
      }
      if (ifound)
	 break;
   }
}

// Do the link restraints contain a planar peptide restraint?
bool
coot::protein_geometry::planar_peptide_restraint_state() const {

   std::string link_id = "TRANS";
   std::string plane_id = "plane3";
   bool ifound = 0;

   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {
      if (dict_link_res_restraints[i].link_id == link_id) { // e.g "TRANS"
	 
	 std::vector<coot::dict_link_plane_restraint_t>::const_iterator it;
	 for (it = dict_link_res_restraints[i].link_plane_restraint.begin();
	      it != dict_link_res_restraints[i].link_plane_restraint.end(); it++) {
	    if (it->plane_id == plane_id) {
	       ifound = 1;
	       break;
	    }
	 }
      }
   }
   return ifound;
} 


// restraints for omega for both CIS and TRANS links (and
// PTRANS)
void
coot::protein_geometry::add_omega_peptide_restraints() {

   double esd = 5.0; // perhaps this should be passed?
   std::vector<std::pair<std::string, double> > v;
   v.push_back(std::pair<std::string, double> ("TRANS",  180.0));
   v.push_back(std::pair<std::string, double> ("PTRANS", 180.0));
   v.push_back(std::pair<std::string, double> ("CIS",    0.0));
   v.push_back(std::pair<std::string, double> ("PCIS",   0.0));

   for (unsigned int iv=0; iv<v.size(); iv++) {
      std::string link_id = v[iv].first;
      // period is 0 (like the dictionary).  A good thing?
      link_add_torsion(link_id, 1, 1, 2, 2, "CA", "C", "N", "CA", v[iv].second, esd, 0, "omega");
   }

}



void
coot::protein_geometry::remove_omega_peptide_restraints() {

   std::vector<std::string> v;
   v.push_back("TRANS");
   v.push_back("PTRANS");
   v.push_back("CIS");
   v.push_back("PCIS");

   bool ifound = 0;
   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {
      if (dict_link_res_restraints[i].link_id == v[i]) { // is TRANS, say

	 std::vector<coot::dict_link_torsion_restraint_t>::iterator it;
	 for (it = dict_link_res_restraints[i].link_torsion_restraint.begin();
	      it != dict_link_res_restraints[i].link_torsion_restraint.end(); it++) {
	    if (it->id() == "omega") {
	       ifound = 1;
 	       dict_link_res_restraints[i].link_torsion_restraint.erase(it);
	       break;
	    }
	 }
      }
   }
}



// a list of three-letter-codes (should that be comp_ids?) that match
// the string in the chem_comp name using the simple_monomer_descriptions
// 
std::vector<std::pair<std::string, std::string> >
coot::protein_geometry::matching_names(const std::string &test_string,
				       short int allow_minimal_descriptions) const {

   std::vector<std::pair<std::string, std::string> > v;
   std::string test_string_dc = coot::util::downcase(test_string);
   std::map<std::string,coot::dictionary_residue_restraints_t>::const_iterator it;


   
   for (it=simple_monomer_descriptions.begin();
	it!=simple_monomer_descriptions.end();
	it++) {
      std::string name_dc = coot::util::downcase(it->second.residue_info.name);
      std::string::size_type ifound = name_dc.find(test_string_dc);
      if (ifound != std::string::npos) {
// 	 std::cout << "test_string :" << test_string << ": matched :"
// 		   << it->second.residue_info.comp_id << ": :"
// 		   << it->second.residue_info.name
// 		   << ":" << std::endl;
 	 std::pair<std::string, std::string> p(it->second.residue_info.comp_id,
					       it->second.residue_info.name);
	 v.push_back(p);
      }
   }

   return v;
}


void
coot::dictionary_residue_restraints_t::write_cif(const std::string &filename) const {

   PCMMCIFFile mmCIFFile = new CMMCIFFile();
      
   PCMMCIFData   mmCIFData = NULL;
   PCMMCIFStruct mmCIFStruct;
   char S[2000];
   
   //  2.1  Example 1: add a structure into mmCIF object

   int rc;

   rc = mmCIFFile->AddMMCIFData("comp_list");
   mmCIFData = mmCIFFile->GetCIFData("comp_list");
   rc = mmCIFData->AddStructure ("_chem_comp", mmCIFStruct);
   // std::cout << "rc on AddStructure returned " << rc << std::endl;
   if (rc!=CIFRC_Ok && rc!=CIFRC_Created)  {
      // badness!
      std::cout << "rc not CIFRC_Ok " << rc << std::endl;
      printf ( " **** error: attempt to retrieve Loop as a Structure.\n" );
      if (!mmCIFStruct)  {
	 printf ( " **** error: mmCIFStruct is NULL - report as a bug\n" );
      }
   } else {
      if (rc==CIFRC_Created) 
	 printf ( " -- new structure created\n" );
      else 
	 printf(" -- structure was already in mmCIF, it will be extended\n");
      std::cout << "SUMMARY:: rc CIFRC_Ok or newly created. " << std::endl;

      PCMMCIFLoop mmCIFLoop = new CMMCIFLoop; // 20100212
      // data_comp_list, id, three_letter_code, name group etc:

      rc = mmCIFData->AddLoop("_chem_comp", mmCIFLoop);
      int i=0;
      const char *s = residue_info.comp_id.c_str();
      mmCIFLoop->PutString(s, "id", i);
      s = residue_info.three_letter_code.c_str();
      mmCIFLoop->PutString(s, "three_letter_code", i);
      s = residue_info.name.c_str();
      mmCIFLoop->PutString(s, "name", i);
      s =  residue_info.group.c_str();
      mmCIFLoop->PutString(s, "group", i);
      int nat = residue_info.number_atoms_all;
      mmCIFLoop->PutInteger(nat, "number_atoms_all", i);
      nat = residue_info.number_atoms_nh;
      mmCIFLoop->PutInteger(nat, "number_atoms_nh", i);
      s = residue_info.description_level.c_str();
      mmCIFLoop->PutString(s, "desc_level", i);
      
      std::string comp_record = "comp_list";
      mmCIFData->PutDataName(comp_record.c_str()); // 'data_' record

      // atom loop

      std::string comp_monomer_name = "comp_";
      comp_monomer_name += residue_info.comp_id.c_str(); 
      rc = mmCIFFile->AddMMCIFData(comp_monomer_name.c_str());
      mmCIFData = mmCIFFile->GetCIFData(comp_monomer_name.c_str());
      rc = mmCIFData->AddLoop("_chem_comp_atom", mmCIFLoop);
      
      if (rc == CIFRC_Ok || rc == CIFRC_Created) {
	 for (int i=0; i<atom_info.size(); i++) {
	    const char *ss =  residue_info.comp_id.c_str();
	    mmCIFLoop->PutString(ss, "comp_id", i);
	    ss = atom_info[i].atom_id.c_str();
	    mmCIFLoop->PutString(ss, "atom_id", i);
	    ss = atom_info[i].type_symbol.c_str();
	    mmCIFLoop->PutString(ss, "type_symbol", i);
	    ss = atom_info[i].type_energy.c_str();
	    mmCIFLoop->PutString(ss, "type_energy", i);
	    if (atom_info[i].partial_charge.first) {
	       float v = atom_info[i].partial_charge.second;
	       mmCIFLoop->PutReal(v, "partial_charge", i);
	    }
	 }
      }

      // bond loop

      rc = mmCIFData->AddLoop("_chem_comp_bond", mmCIFLoop);
      if (rc == CIFRC_Ok || rc == CIFRC_Created) {
	 // std::cout << " number of bonds: " << bond_restraint.size() << std::endl;
	 for (int i=0; i<bond_restraint.size(); i++) {
	    // std::cout << "ading bond number " << i << std::endl;
	    const char *ss = residue_info.comp_id.c_str();
	    mmCIFLoop->PutString(ss, "comp_id", i);
	    ss = bond_restraint[i].atom_id_1_4c().c_str();
	    mmCIFLoop->PutString(ss, "atom_id_1", i);
	    ss = bond_restraint[i].atom_id_2_4c().c_str();
	    mmCIFLoop->PutString(ss, "atom_id_2", i);
	    ss = bond_restraint[i].type().c_str();
	    mmCIFLoop->PutString(ss, "type", i);
	    float v = bond_restraint[i].value_dist();
	    try { 
	       mmCIFLoop->PutReal(v, "value_dist", i);
	       v = bond_restraint[i].value_esd();
	       mmCIFLoop->PutReal(v, "value_dist_esd", i);
	    }
	    catch (std::runtime_error rte) {
	       // do nothing, it's not really an error if the dictionary
	       // doesn't have target geometry (the bonding description came
	       // from a Chemical Component Dictionary entry for example).
	    } 
	 }
      }

      // angle loop

      rc = mmCIFData->AddLoop("_chem_comp_angle", mmCIFLoop);
      if (rc == CIFRC_Ok || rc == CIFRC_Created) {
	 // std::cout << " number of angles: " << angle_restraint.size() << std::endl;
	 for (int i=0; i<angle_restraint.size(); i++) {
	    // std::cout << "ading angle number " << i << std::endl;
	    char *ss = (char *) residue_info.comp_id.c_str();
	    mmCIFLoop->PutString(ss, "comp_id", i);
	    ss = (char *) angle_restraint[i].atom_id_1_4c().c_str();
	    mmCIFLoop->PutString(ss, "atom_id_1", i);
	    ss = (char *) angle_restraint[i].atom_id_2_4c().c_str();
	    mmCIFLoop->PutString(ss, "atom_id_2", i);
	    ss = (char *) angle_restraint[i].atom_id_3_4c().c_str();
	    mmCIFLoop->PutString(ss, "atom_id_3", i);
	    float v = angle_restraint[i].angle();
	    mmCIFLoop->PutReal(v, "angle", i);
	    v = angle_restraint[i].esd();
	    mmCIFLoop->PutReal(v, "esd", i);
	 }
      }

      // torsion loop

      if (torsion_restraint.size() > 0) { 
	 rc = mmCIFData->AddLoop("_chem_comp_tor", mmCIFLoop);
	 if (rc == CIFRC_Ok || rc == CIFRC_Created) {
	    // std::cout << " number of torsions: " << torsion_restraint.size() << std::endl;
	    for (int i=0; i<torsion_restraint.size(); i++) {
	       // std::cout << "ading torsion number " << i << std::endl;
	       char *ss = (char *) residue_info.comp_id.c_str();
	       mmCIFLoop->PutString(ss, "comp_id", i);
	       ss = (char *) torsion_restraint[i].id().c_str();
	       mmCIFLoop->PutString(ss, "id", i);
	       ss = (char *) torsion_restraint[i].atom_id_1_4c().c_str();
	       mmCIFLoop->PutString(ss, "atom_id_1", i);
	       ss = (char *) torsion_restraint[i].atom_id_2_4c().c_str();
	       mmCIFLoop->PutString(ss, "atom_id_2", i);
	       ss = (char *) torsion_restraint[i].atom_id_3_4c().c_str();
	       mmCIFLoop->PutString(ss, "atom_id_3", i);
	       ss = (char *) torsion_restraint[i].atom_id_4_4c().c_str();
	       mmCIFLoop->PutString(ss, "atom_id_4", i);
	       float v = torsion_restraint[i].angle();
	       mmCIFLoop->PutReal(v, "angle", i);
	       v = torsion_restraint[i].esd();
	       mmCIFLoop->PutReal(v, "angle_esd", i);
	       int p = torsion_restraint[i].periodicity();
	       mmCIFLoop->PutInteger(p, "period", i);
	    }
	 }
      }

      // chiral loop
      // 
      if (chiral_restraint.size() > 0) { 
	 rc = mmCIFData->AddLoop("_chem_comp_chir", mmCIFLoop);
	 if (rc == CIFRC_Ok || rc == CIFRC_Created) {
	    // std::cout << " number of chirals: " << chiral_restraint.size() << std::endl;
	    for (int i=0; i<chiral_restraint.size(); i++) {
	       // std::cout << "ading chiral number " << i << std::endl;
	       const char *ss = residue_info.comp_id.c_str();
	       mmCIFLoop->PutString(ss, "comp_id", i);
	       ss = chiral_restraint[i].Chiral_Id().c_str();
	       mmCIFLoop->PutString(ss, "id", i);
	       ss = chiral_restraint[i].atom_id_c_4c().c_str();
	       mmCIFLoop->PutString(ss, "atom_id_centre", i);
	       ss = chiral_restraint[i].atom_id_1_4c().c_str();
	       mmCIFLoop->PutString(ss, "atom_id_1", i);
	       ss = chiral_restraint[i].atom_id_2_4c().c_str();
	       mmCIFLoop->PutString(ss, "atom_id_2", i);
	       ss = chiral_restraint[i].atom_id_3_4c().c_str();
	       mmCIFLoop->PutString(ss, "atom_id_3", i);
	       int sign = chiral_restraint[i].volume_sign;
	       ss = "both";
	       if (sign == 1)
		  ss = "positiv";
	       if (sign == -1)
		  ss = "negativ";
	       mmCIFLoop->PutString(ss, "volume_sign", i);
	    }
	 }
      }

      // plane loop
      if (plane_restraint.size() > 0) { 
	 rc = mmCIFData->AddLoop("_chem_comp_plane_atom", mmCIFLoop);
	 if (rc == CIFRC_Ok || rc == CIFRC_Created) {
	    // std::cout << " number of planes: " << plane_restraint.size() << std::endl;
	    int icount = 0;
	    for (int i=0; i<plane_restraint.size(); i++) {
	       // std::cout << "DEBUG:: adding plane number " << i << std::endl;
	       for (int iat=0; iat<plane_restraint[i].n_atoms(); iat++) {
		  char *ss = (char *) residue_info.comp_id.c_str();
		  mmCIFLoop->PutString(ss, "comp_id", icount);
		  ss = (char *) plane_restraint[i].plane_id.c_str();
		  mmCIFLoop->PutString(ss, "plane_id", icount);
		  ss = (char *) plane_restraint[i].atom_id(iat).c_str();
		  mmCIFLoop->PutString(ss, "atom_id", icount);
		  float v = plane_restraint[i].dist_esd();
		  mmCIFLoop->PutReal(v, "dist_esd", icount);
		  icount++;
	       }
	    }
	 }
      }
      
      mmCIFFile->WriteMMCIFFile(filename.c_str());

   }
   delete mmCIFFile; // deletes all its attributes too.
}



// make a connect file specifying the bonds to Hydrogens
bool
coot::protein_geometry::hydrogens_connect_file(const std::string &resname,
					       const std::string &filename) const {

   bool r = 0;
   std::pair<short int, coot::dictionary_residue_restraints_t> p =
      get_monomer_restraints(resname);

   if (p.first) {
      std::vector<coot::dict_bond_restraint_t> bv = p.second.bond_restraint;
      if (bv.size() > 0) {
	 // try to open the file then:
	 std::ofstream connect_stream(filename.c_str());
	 if (connect_stream) {
	    int n_atoms = p.second.atom_info.size();
	    connect_stream << "# Generated by Coot" << std::endl;
	    connect_stream << "RESIDUE   " << resname << "   " << n_atoms << std::endl;
	    std::vector<std::pair<std::string, std::vector<std::string> > > assoc; 
	    for (unsigned int i=0; i<bv.size(); i++) {
	       std::string atom1 = bv[i].atom_id_1();
	       std::string atom2 = bv[i].atom_id_2();
	       // find atom1
	       bool found = 0;
	       int index_1 = -1;
	       int index_2 = -1;
	       for (unsigned int j=0; j<assoc.size(); j++) {
		  if (atom1 == assoc[j].first) {
		     found = 1;
		     index_1 = j;
		     break;
		  } 
	       }
	       if (found == 1) {
		  assoc[index_1].second.push_back(atom2);
	       } else {
		  // we need to add a new atom:
		  std::vector<std::string> vt;
		  vt.push_back(atom2);
		  std::pair<std::string, std::vector<std::string> > p(atom1, vt);
		  assoc.push_back(p);
	       }
	       // find atom2
	       found = 0;
	       for (unsigned int j=0; j<assoc.size(); j++) {
		  if (atom2 == assoc[j].first) {
		     found = 1;
		     index_2 = j;
		     break;
		  }
	       }
	       if (found == 1) {
		  assoc[index_2].second.push_back(atom1);
	       } else {
		  // we need to add a new atom:
		  std::vector<std::string> vt;
		  vt.push_back(atom1);
		  std::pair<std::string, std::vector<std::string> > p(atom2, vt);
		  assoc.push_back(p);
	       } 
	    }

	    r = 1;
	    // for each atom in assoc
	    for (unsigned int i=0; i<assoc.size(); i++) {
	       
	       connect_stream << "CONECT     " << assoc[i].first << "    "
			      << assoc[i].second.size();
	       for (unsigned int ii=0; ii<assoc[i].second.size(); ii++) {
		  connect_stream << assoc[i].second[ii] << " ";
	       }
	       connect_stream << std::endl;
	    }
	 }
      }
   }
   return r;
} 
      
// constructor
coot::simple_cif_reader::simple_cif_reader(const std::string &cif_dictionary_file_name) {
   
   CMMCIFFile ciffile;
   struct stat buf;
   int istat = stat(cif_dictionary_file_name.c_str(), &buf);
   if (istat != 0) {
      std::cout << "WARNIG:: cif dictionary " << cif_dictionary_file_name
		<< " not found" << std::endl;
   } else {
      int ierr = ciffile.ReadMMCIFFile((char *)cif_dictionary_file_name.c_str());
      if (ierr != CIFRC_Ok) {
	 std::cout << "Dirty mmCIF file? " << cif_dictionary_file_name
		   << std::endl;
      } else {
	 for(int idata=0; idata<ciffile.GetNofData(); idata++) { 
         
	    PCMMCIFData data = ciffile.GetCIFData(idata);
	    for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 
	       PCMMCIFCategory cat = data->GetCategory(icat);
	       std::string cat_name(cat->GetCategoryName());
	       PCMMCIFLoop mmCIFLoop =
		  data->GetLoop(cat_name.c_str() );
	       if (mmCIFLoop == NULL) { 
		  std::cout << "null loop" << std::endl; 
	       } else {
		  if (cat_name == "_chem_comp") {
		     int ierr = 0;
		     for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {
			char *n = mmCIFLoop->GetString("name", j, ierr);
			char *t = mmCIFLoop->GetString("three_letter_code", j, ierr);
			if (n && t) {
			   names.push_back(n);
			   three_letter_codes.push_back(t);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
}

bool
coot::simple_cif_reader::has_restraints_for(const std::string &res_type) {

   bool r = 0;
   for (unsigned int i=0; i<three_letter_codes.size(); i++) {
      if (three_letter_codes[i] == res_type) {
	 r = 1;
	 break;
      }
   }
   return r;
}

// replace (return 1)  or add (if not replacable) (return 0).
// 
bool
coot::protein_geometry::replace_monomer_restraints(std::string monomer_type,
						   const coot::dictionary_residue_restraints_t &mon_res_in) {
   bool s = 0;

   coot::dictionary_residue_restraints_t mon_res = mon_res_in;
   mon_res.assign_chiral_volume_targets();
   
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].residue_info.comp_id == monomer_type) {
	 dict_res_restraints[i] = mon_res;
	 s = 1;
	 break;
      }
   }

   if (s == 0) {
      dict_res_restraints.push_back(mon_res);
   } 
   return s;
}

std::vector<std::string>
coot::protein_geometry::monomer_types() const {
   std::vector<std::string> v;
   for (int i=0; i<dict_res_restraints.size(); i++) {
      v.push_back(dict_res_restraints[i].residue_info.comp_id);
   }
   return v;
}

// Thow an exception if we can't get the group of r.
//
// 20100420 First compare against the three_letter_code and if not
// found try testing against the comp_id (added 20100420).  This is
// needed so that this function returns the group from a Ur/U (that's
// the comp_id/tlc) [and similarly for other bases].  This is needed
// so that find_link_type_rigourous() works for link type "p" (RNA/DNA
// stuff).  (Problem found when trying to sphere refine a on RNA -
// 3l0u).
// 
std::string
coot::protein_geometry::get_group(CResidue *r) const {

   std::string res_name = r->GetResName();
   return get_group(res_name);
}


std::string
coot::protein_geometry::get_group(const std::string &res_name_in) const {
   
   bool found = 0;
   std::string group;
   std::string res_name = res_name_in;
   if (res_name.length() > 3)
      res_name = res_name.substr(0,2);
   for (unsigned int i=0; i<size(); i++) {
      if (three_letter_code(i) == res_name) {
	 found = 1;
	 group = (*this)[i].residue_info.group;
	 break;
      }
   }

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) { 
      if (dict_res_restraints[i].residue_info.comp_id == res_name) {
	 found = 1;
	 group = dict_res_restraints[i].residue_info.group;
	 break;
      }
   }

   if (! found) {
      std::string s = "No dictionary group found for residue type :";
      s += res_name;
      s += ":";
      throw std::runtime_error(s);
   }
   return group;
}

CResidue *
coot::protein_geometry::get_residue(const std::string &comp_id, bool idealised_flag) {

   CResidue *residue_p = NULL;

   // force use of try_dynamic_add (if needed).
   bool r = have_dictionary_for_residue_type(comp_id, 42);

   bool make_hetatoms = ! coot::util::is_standard_residue_name(comp_id);
   
   std::vector<CAtom *> atoms;
   for (int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].residue_info.comp_id == comp_id) {

	 std::vector<coot::dict_atom> atom_info = dict_res_restraints[i].atom_info;
	 int atom_index = 0;
	 for (unsigned int iat=0; iat<atom_info.size(); iat++) {
	    

	    clipper::Coord_orth p(0,0,0);
	    bool flag_and_have_coords = 0;

	    if (idealised_flag && atom_info[iat].pdbx_model_Cartn_ideal.first) {
	       p = atom_info[iat].pdbx_model_Cartn_ideal.second;
	       flag_and_have_coords = 1;
	    }

	    if (! flag_and_have_coords) { 
	       // OK, try model_Cartn (and that is idealised if the dictionary was refmac)
	       // (better than nothing).
	       // 
	       if (atom_info[iat].model_Cartn.first) {
		  p = atom_info[iat].model_Cartn.second;
		  flag_and_have_coords = 1;
	       }
	    }
	    
	    if (flag_and_have_coords) { 
	       CAtom *atom = new CAtom;
	       realtype occ = 1.0;
	       realtype b = 20.0;
	       std::string ele = atom_info[iat].type_symbol; // element
	       atom->SetCoordinates(p.x(), p.y(), p.z(), occ, b);
	       atom->SetAtomName(atom_info[iat].atom_id_4c.c_str());

//	       Strange things happen when I use this...
// 	       atom->SetAtomName(atom_index, -1,
// 				 atom_info[iat].atom_id_4c.c_str(),
// 				 "", "", ele.c_str());
	       
	       atom->SetElementName(ele.c_str());
	       if (make_hetatoms)
		  atom->Het = 1;
	       atoms.push_back(atom);
	       atom_index++; // for next round
	    }
	 }
      }
   }
   if (atoms.size() > 0) {
      residue_p = new CResidue;
      residue_p->SetResID(comp_id.c_str(), 1, "");
      for (unsigned int iat=0; iat<atoms.size(); iat++) 
	 residue_p->AddAtom(atoms[iat]);
   }
   return residue_p;
} 



CMMDBManager *
coot::protein_geometry::mol_from_dictionary(const std::string &three_letter_code,
					    bool idealised_flag) {

   CMMDBManager *mol = NULL;
   CResidue *residue_p = get_residue(three_letter_code, idealised_flag);
   if (residue_p) { 
      CChain *chain_p = new CChain;
      chain_p->SetChainID("A");
      chain_p->AddResidue(residue_p);
      CModel *model_p = new CModel;
      model_p->AddChain(chain_p);
      mol = new CMMDBManager;
      mol->AddModel(model_p);
   }
   return mol;
}

void
coot::protein_geometry::print_chem_links() const {

   for (unsigned int i_chem_link=0; i_chem_link<chem_link_vec.size(); i_chem_link++) {
      std::cout<< i_chem_link << " " << chem_link_vec[i_chem_link] << "\n";
   } 

} 

// delete comp_id from dict_res_restraints (if it exists there).
bool
coot::protein_geometry::delete_mon_lib(std::string comp_id) {

   bool deleted = false; 
   std::vector<coot::dictionary_residue_restraints_t>::iterator it;
   for (it=dict_res_restraints.begin(); it!=dict_res_restraints.end(); it++) {
      if (it->residue_info.comp_id == comp_id) { 
	 dict_res_restraints.erase(it);
	 deleted = 1;
	 break;
      }
   }
   return deleted;
} 

bool
coot::protein_geometry::linkable_residue_types_p(const std::string &this_res_type,
						 const std::string &env_res_type) {

   std::pair<short int, coot::dictionary_residue_restraints_t> r1 = get_monomer_restraints(this_res_type);
   std::pair<short int, coot::dictionary_residue_restraints_t> r2 = get_monomer_restraints(env_res_type);

   bool r = 0;
   if (r1.first) {
      if (r1.second.residue_info.group != "non-polymer")
	 r = 1;
   }
   if (r2.first) {
      if (r2.second.residue_info.group != "non-polymer")
	 r = 1;
   }
   return r;
} 

bool
coot::protein_geometry::OXT_in_residue_restraints_p(const std::string &residue_type) const {

   bool r = 0;
   std::pair<bool, coot::dictionary_residue_restraints_t> p = get_monomer_restraints(residue_type);
   if (p.first) {
      for (unsigned int i=0; i<p.second.atom_info.size(); i++) {
	 if (p.second.atom_info[i].atom_id_4c == " OXT") {
	    r = 1;
	    break;
	 }
      }
   } else {
      if (0) 
	 std::cout << "INFO:: residue type :" << residue_type << ": not found in dictionary"
		   << std::endl;
   } 
   return r;
}


std::vector<std::vector<std::string> >
coot::dictionary_residue_restraints_t::get_ligand_ring_list() const {

   // get a list of aromatic bonds, so that they can be used to find
   // aromatic rings.
   // 
   std::vector<std::pair<std::string, std::string> > bonds;
   for (unsigned int irest=0; irest<bond_restraint.size(); irest++) {
      std::pair<std::string, std::string> p(bond_restraint[irest].atom_id_1_4c(),
					    bond_restraint[irest].atom_id_2_4c());
      bonds.push_back(p);
   }
   
   coot::aromatic_graph_t bond_list(bonds);
   std::vector<std::vector<std::string> > ring_list = bond_list.ring_list();

   if (0) {
      std::cout << "----------- " << ring_list.size() << " rings ---------- " << std::endl;
      for (unsigned int i=0; i<ring_list.size(); i++) {
	 std::cout << "ring " << i << "\n   ";
	 for (unsigned int j=0; j<ring_list[i].size(); j++) { 
	    std::cout << ring_list[i][j] << "  ";
	 }
	 std::cout << std::endl;
      }
   }
   return ring_list;
}

std::vector<std::vector<std::string> >
coot::dictionary_residue_restraints_t::get_ligand_aromatic_ring_list() const {

   // get a list of aromatic bonds, so that they can be used to find
   // aromatic rings.
   // 
   std::vector<std::pair<std::string, std::string> > bonds;
   for (unsigned int irest=0; irest<bond_restraint.size(); irest++) {
      if (bond_restraint[irest].type() == "aromatic") {
	 std::pair<std::string, std::string> p(bond_restraint[irest].atom_id_1_4c(),
					       bond_restraint[irest].atom_id_2_4c());
	 bonds.push_back(p);
      }
   }
   
   coot::aromatic_graph_t arom(bonds);
   std::vector<std::vector<std::string> > ring_list = arom.ring_list();

   if (0) {
      std::cout << "----------- " << ring_list.size() << " rings ---------- " << std::endl;
      for (unsigned int i=0; i<ring_list.size(); i++) {
	 std::cout << "ring " << i << "\n   ";
	 for (unsigned int j=0; j<ring_list[i].size(); j++) { 
	    std::cout << ring_list[i][j] << "  ";
	 }
	 std::cout << std::endl;
      }
   }
   return ring_list;
}



// Find the bonded neighbours of the given atoms - throw an
// exception if residue name is not in dictionary.
// 
std::vector<std::string>
coot::protein_geometry::get_bonded_neighbours(const std::string &residue_name,
					      const std::string &atom_name_1,
					      const std::string &atom_name_2, 
					      bool also_2nd_order_neighbs_flag) const {

   std::vector<std::string> v;

   std::vector<std::string> v_2nd_order; // only filled for triple bonds

   std::pair<bool, coot::dictionary_residue_restraints_t> restraints =
      get_monomer_restraints_at_least_minimal(residue_name);

   if (restraints.first) { 
      for (unsigned int i=0; i<restraints.second.bond_restraint.size(); i++) {
// 	 std::cout << "-- comparing " << atom_name_1 << " " << atom_name_2 << " -- to -- "
// 		   << restraints.second.bond_restraint[i].atom_id_1_4c() << " "
// 		   << restraints.second.bond_restraint[i].atom_id_2_4c() << std::endl;
	 if (restraints.second.bond_restraint[i].atom_id_1_4c() == atom_name_1)
	    if (restraints.second.bond_restraint[i].atom_id_2_4c() != atom_name_2) {
	       // std::cout << " adding a " << restraints.second.bond_restraint[i].atom_id_2_4c() << std::endl;
	       std::string at_name = restraints.second.bond_restraint[i].atom_id_2_4c();
	       v.push_back(at_name);
	       if (also_2nd_order_neighbs_flag) { 
		  std::vector<std::string> nv = 
		     coot::protein_geometry::get_bonded_neighbours(residue_name, 
								   atom_name_1, at_name);
		  for (unsigned int in=0; in<nv.size(); in++)
		     v_2nd_order.push_back(nv[in]);
	       }
	    }
	 if (restraints.second.bond_restraint[i].atom_id_1_4c() == atom_name_2)
	    if (restraints.second.bond_restraint[i].atom_id_2_4c() != atom_name_1) { 
	       // std::cout << " adding b " << restraints.second.bond_restraint[i].atom_id_2_4c() << std::endl;
	       std::string at_name = restraints.second.bond_restraint[i].atom_id_2_4c();
	       v.push_back(at_name);
	       if (also_2nd_order_neighbs_flag) { 
		  std::vector<std::string> nv = 
		     coot::protein_geometry::get_bonded_neighbours(residue_name, 
								   atom_name_1, at_name);
		  for (unsigned int in=0; in<nv.size(); in++)
		     v_2nd_order.push_back(nv[in]);
	       }
	    }
	 if (restraints.second.bond_restraint[i].atom_id_2_4c() == atom_name_1)
	    if (restraints.second.bond_restraint[i].atom_id_1_4c() != atom_name_2) {
	       // std::cout << " adding c " << restraints.second.bond_restraint[i].atom_id_1_4c() << std::endl;
	       std::string at_name = restraints.second.bond_restraint[i].atom_id_1_4c();
	       v.push_back(at_name);
	       if (also_2nd_order_neighbs_flag) { 
		  std::vector<std::string> nv = 
		     coot::protein_geometry::get_bonded_neighbours(residue_name, 
								   atom_name_1, at_name);
		  for (unsigned int in=0; in<nv.size(); in++)
		     v_2nd_order.push_back(nv[in]);
	       }
	    }
	 if (restraints.second.bond_restraint[i].atom_id_2_4c() == atom_name_2)
	    if (restraints.second.bond_restraint[i].atom_id_1_4c() != atom_name_1) {
	       // std::cout << " adding d " << restraints.second.bond_restraint[i].atom_id_1_4c() << std::endl;
	       std::string at_name = restraints.second.bond_restraint[i].atom_id_1_4c();
	       v.push_back(at_name);
	       if (also_2nd_order_neighbs_flag) { 
		  std::vector<std::string> nv = 
		     coot::protein_geometry::get_bonded_neighbours(residue_name, 
								   atom_name_1, at_name);
		  for (unsigned int in=0; in<nv.size(); in++)
		     v_2nd_order.push_back(nv[in]);
	       }
	    }
      }

      // add the neighbour neighbours to v (if needed):
      if (also_2nd_order_neighbs_flag) { 
	 for(unsigned int in=0; in<v_2nd_order.size(); in++)
	    if (std::find(v.begin(), v.end(), v_2nd_order[in]) == v.end())
	       v.push_back(v_2nd_order[in]);
      } 

      if (v.size()) {
	 // add the initial atom names if they are not already there
	 if (std::find(v.begin(), v.end(), atom_name_1) == v.end())
	    v.push_back(atom_name_1);
	 if (std::find(v.begin(), v.end(), atom_name_2) == v.end())
	    v.push_back(atom_name_2);
      } 
   } else {
      std::string m = "No dictionary for ";
      m += residue_name;
      throw std::runtime_error(m);
   } 
   return v;
} 
