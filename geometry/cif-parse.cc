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

// return the number of atoms read (not the number of bonds (because
// that is not a good measure of having read the file properly for
// (for example) CL)).
// 
coot::read_refmac_mon_lib_info_t
coot::protein_geometry::init_refmac_mon_lib(std::string ciffilename, int read_number_in,
					    int imol_enc_in) {

   if (false)
      std::cout << "DEBUG:: init_refmac_mon_lib() " << ciffilename << " "
		<< read_number_in << " imol_enc_in: " << imol_enc_in << std::endl;

   int imol_enc = imol_enc_in;
   if (imol_enc_in == IMOL_ENC_UNSET) {
      imol_enc = IMOL_ENC_ANY;

   } else {
      
      // it's too late (too deep) by the time we get here, to use IMOL_ENC_AUTO - that should
      // have been done by the calling function (where, if auto, then set imol_enc_in to imol_no)
      // for a specific model
      // or
      // IMOL_ENC_ANY for any model.
      //
      // so we need to provide outside access to something like this, by which the outside
      // function can decide whether to use imol_no or IMOL_ENC_ANY
      // if (std::find(non_auto_load_residue_names, comp_id) == non_auto_load_residue_names.end())..
      // we have it:  is_non_auto_load_ligand()

   }
      

   read_refmac_mon_lib_info_t rmit;
   mmdb::mmcif::File ciffile;
   std::vector<std::string> comp_ids; // found in the ciffile

   // Here we would want to check through dict_res_restraints for the
   // existance of this restraint.  If it does exist, kill it by
   // changing its name to blank.  However, we don't know the name of
   // the restraint yet!  We only know that at the add_atom(),
   // add_bond() [etc] stage.

   // added 20120121, to fix Jack Sack and Kevin Madauss bug (so that
   // in mon_lib_add_chem_comp(), the read numbers between new file
   // and previous are different, hence trash the old restraints).
   // 
   read_number = read_number_in;

   //    struct stat buf;
   //    int istat = stat(ciffilename.c_str(), &buf);
   // Thanks Ezra Peisach for this this bug report

   if (! is_regular_file(ciffilename)) {
      std::string s = "WARNING: in init_refmac_mon_lib, file \"";
      s += ciffilename;
      s += "\" not found.";
      std::cout <<  s << "\n";
      rmit.error_messages.push_back(s);
      rmit.success = false;
      return rmit;

   } else {

      int ierr = ciffile.ReadMMCIFFile(ciffilename.c_str());
      std::string comp_id_1; 
      std::string comp_id_2;  // initially unset
   
      if (ierr!=mmdb::mmcif::CIFRC_Ok) {
	 std::cout << "dirty mmCIF file? " << ciffilename << std::endl;
	 std::cout << "    Bad mmdb::mmcif::CIFRC_Ok on ReadMMCIFFile" << std::endl;

	 std::cout << "    " << mmdb::GetErrorDescription(mmdb::ERROR_CODE(ierr))
		   << std::endl;
	 
 	 char        err_buff[1000];
 	 std::cout <<  "CIF error rc=" << ierr << " reason:" << 
 	    mmdb::mmcif::GetCIFMessage (err_buff, ierr) << std::endl;

	 rmit.success = false;
	 std::string s = "Dirty mmCIF file? ";
	 s += ciffilename;
	 rmit.error_messages.push_back(s);
	 s = "Bad mmdb::mmcif::CIFRC_Ok on ReadMMCIFFile";
	 rmit.error_messages.push_back(s);
	 s = mmdb::GetErrorDescription(mmdb::ERROR_CODE(ierr));
	 rmit.error_messages.push_back(s);
	 clipper::String cs = "CIF error rc=";
	 cs += ierr;
	 cs += " reason:";
	 cs += mmdb::mmcif::GetCIFMessage (err_buff, ierr);
	 rmit.error_messages.push_back(cs);

      } else {
	 if (verbose_mode)
	    std::cout << "There are " << ciffile.GetNofData() << " data in "
		      << ciffilename << std::endl; 
      
	 for(int idata=0; idata<ciffile.GetNofData(); idata++) {
	    
	    mmdb::mmcif::PData data = ciffile.GetCIFData(idata);

	    // note that chem_link goes here to:
	    // 
	    if (std::string(data->GetDataName()).substr(0,5) == "link_") {
	       rmit.n_links += init_links(data);
	    }

	    if (std::string(data->GetDataName()).length() > 7) {
	       if (std::string(data->GetDataName()).substr(0,8) == "mod_list") {
		  // this handles all "list_chem_mod"s in the file (just one call)
		  rmit.n_links += add_chem_mods(data); // check this

		  if (0) // debug
		     for (unsigned int i=0; i<chem_mod_vec.size(); i++)
			std::cout << "     " << chem_mod_vec[i] << std::endl;
	       }
	    }

	    if (std::string(data->GetDataName()).length() > 4) {
	       // e.g. mod_NH3 ? 
	       if (std::string(data->GetDataName()).substr(0,4) == "mod_") {
		  rmit.n_links += add_mod(data);
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

	       mmdb::mmcif::PCategory cat = data->GetCategory(icat);
	       std::string cat_name(cat->GetCategoryName());
	       
	       // All catagories have loops (AFAICS). 
	       // std::cout << "DEBUG:: got catagory: " << cat_name << std::endl; 

	       mmdb::mmcif::PLoop mmCIFLoop = data->GetLoop(cat_name.c_str() );

	       int n_loop_time = 0;
	       if (mmCIFLoop == NULL) {

		  bool handled = false;
		  if (cat_name == "_chem_comp") {
		     // read the chemical component library which does
		     // not have a loop (the refmac files do) for the
		     // chem_comp info.
		     handled = 1;
		     mmdb::mmcif::PStruct structure = data->GetStructure(cat_name.c_str());
		     if (structure) {
			comp_id_1 = chem_comp_component(structure, imol_enc);
		     }
		  }

		  if (cat_name == "_pdbx_chem_comp_model") {
		     handled = 1;
		     mmdb::mmcif::PStruct structure = data->GetStructure(cat_name.c_str());
		     if (structure) {
			comp_id_1 = pdbx_chem_comp_model(structure, imol_enc);
		     }
		  }

		  if (cat_name == "_chem_comp_chir") {
		     handled = 1;
		     mmdb::mmcif::PStruct structure = data->GetStructure(cat_name.c_str());
		     if (structure) {
			chem_comp_chir_structure(structure, imol_enc);
		     }
		  }

		  if (cat_name == "_chem_comp_tor") {
		     handled = 1;
		     mmdb::mmcif::PStruct structure = data->GetStructure(cat_name.c_str());
		     if (structure) {
			chem_comp_tor_structure(structure, imol_enc);
		     }
		  }
		  
		  if (! handled)   // this can happen if there is not an atom loop, e.g. dictionary
		                   // with one atom e.g. AM.cif (Americium ion)
		     std::cout << "WARNING:: in init_refmac_mon_lib() unhandled catagory \""
			       << cat_name << "\" file: " << ciffilename << std::endl; 
		  
	       } else {
               
		  n_loop_time++;

		  // We currently want to stop adding chem comp info
		  // if the chem_comp info comes from mon_lib_list.cif:
		  if (cat_name == "_chem_comp") {
		     if (read_number_in != coot::protein_geometry::MON_LIB_LIST_CIF)
			comp_id_2 = chem_comp(mmCIFLoop, imol_enc);
		     else
			comp_id_2 = simple_mon_lib_chem_comp(mmCIFLoop, imol_enc);
		  }

		  // monomer info, name, number of atoms etc.
		  if (cat_name == "_chem_comp_atom")
		     rmit.n_atoms += comp_atom(mmCIFLoop, imol_enc); // and at the end pad up the atom names

		  // tree
		  if (cat_name == "_chem_comp_tree")
		     comp_tree(mmCIFLoop, imol_enc);

		  // bond
		  if (cat_name == "_chem_comp_bond")
		     rmit.n_bonds += comp_bond(mmCIFLoop, imol_enc);

		  // angle
		  if (cat_name == "_chem_comp_angle")
		     comp_angle(mmCIFLoop, imol_enc);

		  // tor
		  if (cat_name == "_chem_comp_tor")
		     comp_torsion(mmCIFLoop, imol_enc);

		  // chiral
		  if (cat_name == "_chem_comp_chir") {
		     std::pair<int, std::vector<std::string> > chirals = 
			comp_chiral(mmCIFLoop, imol_enc);
		     n_chiral += chirals.first;
		     for (unsigned int ichir=0; ichir<chirals.second.size(); ichir++)
			comp_ids_for_chirals.push_back(chirals.second[ichir]);
		  }

		  // plane
		  if (cat_name == "_chem_comp_plane_atom")
		     comp_plane(mmCIFLoop, imol_enc);

		  // PDBx stuff
		  if (cat_name == "_pdbx_chem_comp_descriptor")
		     pdbx_chem_comp_descriptor(mmCIFLoop, imol_enc);

		  // PDBx model molecule
		  if (cat_name == "_pdbx_chem_comp_model_atom")
		     rmit.n_atoms += comp_atom(mmCIFLoop, imol_enc, true);
		  if (cat_name == "_pdbx_chem_comp_model_bond")
		     rmit.n_atoms += comp_bond(mmCIFLoop, imol_enc, true);
	       }
	    }
	    if (n_chiral) {
	       assign_chiral_volume_targets();
	       filter_chiral_centres(imol_enc, comp_ids_for_chirals);
	    }
	 } // for idata
	 add_cif_file_name(ciffilename, comp_id_1, comp_id_2, imol_enc);
      } // cif file is OK test

      if (false)
	 std::cout << "returning with comp_id_1 " << comp_id_1
		   << " comp_id_2 " << comp_id_2 << std::endl;

      bool allow_minimal_flag = true;
      if (! comp_id_1.empty()) {
	 comp_ids.push_back(comp_id_1);
	 rmit.monomer_idx = get_monomer_restraints_index(comp_id_1, imol_enc, allow_minimal_flag);
      }
      if (! comp_id_2.empty()) {
	 comp_ids.push_back(comp_id_2);
	 rmit.monomer_idx = get_monomer_restraints_index(comp_id_2, imol_enc, allow_minimal_flag);
      }
   } // is regular file test

   // debug_mods();

   if (comp_ids.size() > 0)
      rmit.comp_id = comp_ids[0];

   return rmit;
}


// 
void
coot::protein_geometry::add_cif_file_name(const std::string &cif_file_name,
					  const std::string &comp_id1,
					  const std::string &comp_id2,
					  int imol_enc) {

   std::string comp_id = comp_id1;
   if (comp_id == "")
      comp_id = comp_id2;
   if (comp_id != "") {
      int idx = get_monomer_restraints_index(comp_id2, imol_enc, true);
      if (idx != -1) {
	 dict_res_restraints[idx].second.cif_file_name = cif_file_name;
      }
   }
}


std::string
coot::protein_geometry::get_cif_file_name(const std::string &comp_id,
					  int imol_enc) const {

   std::string r; 
   int idx = get_monomer_restraints_index(comp_id, imol_enc, true);
   if (idx != -1)
      r = dict_res_restraints[idx].second.cif_file_name;
   return r;
}
      


// return the comp id (so that later we can associate the file name with the comp_id).
// 
std::string 
coot::protein_geometry::chem_comp_component(mmdb::mmcif::PStruct structure, int imol_enc) {

   int n_tags = structure->GetNofTags();
   std::string cat_name = structure->GetCategoryName();

   if (false)
      std::cout << "DEBUG: ================= chem_comp_component() by structure: in category "
		<< cat_name << " there are "
		<< n_tags << " tags" << std::endl;
    
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
      if (tag == "desc_level") // not descr_level
	 description_level = std::pair<bool, std::string> (1,field);
      if (tag == "description_level")
	 description_level = std::pair<bool, std::string> (1,field);
      // number of atoms here too.

      if (tag == "number_atoms_all") { 
	 try {
	    number_of_atoms_all = coot::util::string_to_int(field);
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << rte.what() << std::endl;
	 }
      }
      if (tag == "number_atoms_nh") { 
	 try {
	    number_of_atoms_nh = coot::util::string_to_int(field);
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << rte.what() << std::endl;
	 }
      }
   }

   if (false) 
      std::cout
	 << "chem_comp_component() comp_id :" << comp_id.first << " :" << comp_id.second << ": "
	 << "three_letter_code :" << three_letter_code.first << " :" << three_letter_code.second
	 << ": "
	 << "name :" << name.first << " :" << name.second << ": "
	 << "type :" << type.first << " :" << type.second << ": "
	 << "description_level :" << description_level.first << " :" << description_level.second
	 << ": "
	 << std::endl;

   if (comp_id.first && three_letter_code.first && name.first) {

      mon_lib_add_chem_comp(comp_id.second, imol_enc,
			    three_letter_code.second,
			    name.second, type.second,
			    number_of_atoms_all, number_of_atoms_nh,
			    description_level.second);
   } else { 
      // std::cout << "oooppps - something missing, not adding that" << std::endl;
   }

   if (comp_id.first)
      return comp_id.second;
   else
      return "";
}

std::string
coot::protein_geometry::pdbx_chem_comp_model(mmdb::mmcif::PStruct structure, int imol_enc) {

   std::string id;
   int n_tags = structure->GetNofTags();
   for (int itag=0; itag<n_tags; itag++) {
      std::string tag = structure->GetTag(itag);
      std::string field = structure->GetField(itag);
      if (tag == "id")
	 id = field;
   }
   return id;
}




// non-looping (single) tor
void
coot::protein_geometry::chem_comp_tor_structure(mmdb::mmcif::PStruct structure, int imol_enc) {
   
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
   std::pair<bool, mmdb::realtype> value_angle(0, 0);
   std::pair<bool, mmdb::realtype> value_angle_esd(0, 0);
   
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
	 catch (const std::runtime_error &rte) {
	    std::cout << "WARNING:: not an integer: " << field << std::endl;
	 }
      }
      if (tag == "value_angle") { 
	 try { 
	    value_angle = std::pair<bool, float> (1,coot::util::string_to_float(field));
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "WARNING:: value_angle not an number: " << field << std::endl;
	 }
      }
      if (tag == "value_angle_esd") { 
	 try { 
	    value_angle_esd = std::pair<bool, float> (1,coot::util::string_to_float(field));
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "WARNING:: value_angle_esd not an number: " << field << std::endl;
	 }
      }
   }

   if (comp_id.first && 
       atom_id_1.first && atom_id_2.first && atom_id_3.first && atom_id_4.first &&
       value_angle.first && value_angle_esd.first && 
       period.first) {
      mon_lib_add_torsion(comp_id.second,
			  imol_enc,
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
coot::protein_geometry::chem_comp_chir_structure(mmdb::mmcif::PStruct structure, int imol_enc) {

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
			 imol_enc,
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


// add to simple_monomer_descriptions not dict_res_restraints.
void
coot::protein_geometry::simple_mon_lib_add_chem_comp(const std::string &comp_id,
						     int imol_enc,
						     const std::string &three_letter_code,
						     const std::string &name,
						     const std::string &group,
						     int number_atoms_all, int number_atoms_nh,
						     const std::string &description_level) {


   if (false)
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
coot::dict_atom::add_pos(int pos_type,
			 const std::pair<bool, clipper::Coord_orth> &model_pos) {

   if (pos_type == coot::dict_atom::IDEAL_MODEL_POS)
      pdbx_model_Cartn_ideal = model_pos;
   if (pos_type == coot::dict_atom::REAL_MODEL_POS) {
      model_Cartn = model_pos;
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
					 dict_bond_restraint_t::aromaticity_t arom_in) {

   if (false)
      std::cout << "adding bond for " << comp_id << " " << atom_id_1
		<< " " << atom_id_2 << " " << type << " " << value_dist
		<< " " << value_dist_esd << std::endl;

   // add a bond restraint to the list for comp_id.
   // The list container for comp_id is a dictionary_residue_restraints_t

   add_restraint(comp_id, imol_enc, dict_bond_restraint_t(atom_id_1,
							  atom_id_2,
							  type,
							  value_dist,
							  value_dist_esd,
							  arom_in));
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
coot::protein_geometry::add_restraint(std::string comp_id,
				      int imol_enc,
				      const dict_bond_restraint_t &restr) { 

   // if comp is in the list, simply push back restr, 
   // 
   // if not, then push back a dict_bond_restraint_t for it, passing
   // the comp_id. 

   bool ifound = 0;

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) { 
      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
	 if (dict_res_restraints[i].first == imol_enc) {
	    ifound = 1;
	    dict_res_restraints[i].second.bond_restraint.push_back(restr); 
	    break;
	 }
      }
   } 

   // it was not there
   if (! ifound) {
      dictionary_residue_restraints_t rest(comp_id, read_number);
      std::pair<int, dictionary_residue_restraints_t> p(imol_enc, rest);
      dict_res_restraints.push_back(p);
      // add the bond to the newly created dictionary_residue_restraints_t
      dict_res_restraints[dict_res_restraints.size()-1].second.bond_restraint.push_back(restr);
   }
}

void
coot::protein_geometry::add_restraint(std::string comp_id,
				      int imol_enc,
				      const dict_angle_restraint_t &restr) {

   // if comp is in the list, simply push back restr,
   // 
   // if not, then push back a dict_bond_restraint_t for it, passing
   // the comp_id.

   short int ifound = 0;

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
	 if (dict_res_restraints[i].first == imol_enc) {
	    ifound = 1;
	    dict_res_restraints[i].second.angle_restraint.push_back(restr);
	    break;
	 }
      }
   }

   // it was not there
   if (! ifound) {
      dictionary_residue_restraints_t rest(comp_id, read_number);
      std::pair<int, dictionary_residue_restraints_t> p(imol_enc, rest);
      dict_res_restraints.push_back(p);
      // add the angle to the newly created dictionary_residue_restraints_t
      dict_res_restraints[dict_res_restraints.size()-1].second.angle_restraint.push_back(restr);
   }
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

// We currently want to stop adding chem comp info
// if the chem_comp info comes from mon_lib_list.cif:
//
// This is the function that we use read things other than
// MON_LIB_LIST_CIF.  i.e. bespoke ligand dictionaries.
//
// return the chem_comp
// 
std::string
coot::protein_geometry::chem_comp(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc) {

   int ierr = 0;
   std::string returned_chem_comp; 

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

      std::string id; // gets stored as comp_id, but is labelled "id"
		      // in the cif file chem_comp block.

      // modify a reference (ierr)
      // 
      std::string three_letter_code;
      std::string name;
      std::string group; // e.g. "L-peptide"
      int number_atoms_all;
      int number_atoms_nh;
      std::string description_level = "None";
      
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
      if (s) {
	 group = s; // e.g. "L-peptide"
	 if (group == "L-PEPTIDE") // fix acedrg output
	    group = "L-peptide";
      }

      ierr = mmCIFLoop->GetInteger(number_atoms_all, "number_atoms_all", j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetInteger(number_atoms_nh, "number_atoms_nh", j);
      ierr_tot += ierr;

      char *release_status_cs = mmCIFLoop->GetString("release_status", j, ierr);
      std::string release_status;
      if (release_status_cs)
	 release_status = release_status_cs; // can be "OBS" or "REL"

      // If desc_level is in the file, extract it, otherwise set it to "None"
      //
      int ierr_description = 0;
      s = mmCIFLoop->GetString("desc_level", j, ierr_description);
      if (! ierr_description) {
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
	 description_level = "."; // full
      }

      if (ierr_tot != 0) {
	 std::cout << "oops:: ierr_tot was " << ierr_tot << std::endl;
      } else {
	 // std::cout << "--------- chem_comp() calls delete_mon_lib() " << id << " " << imol_enc
	 // << std::endl;
	 delete_mon_lib(id, imol_enc); // delete it if it exists already.
	 mon_lib_add_chem_comp(id, imol_enc,
			       three_letter_code, name,
			       group, number_atoms_all, number_atoms_nh,
			       description_level);
	 returned_chem_comp = id;
      }
   }
   return returned_chem_comp;
}

// We currently want to stop adding chem comp info
// if the chem_comp info comes from mon_lib_list.cif:
//
// This is the function that we use read MON_LIB_LIST_CIF
// 
// return the comp_id
std::string
coot::protein_geometry::simple_mon_lib_chem_comp(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc) {

   int ierr = 0;
   std::string comp_id;
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
	 comp_id = s;
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
	    simple_mon_lib_add_chem_comp(comp_id, imol_enc,
					 three_letter_code, name,
					 group, number_atoms_all, number_atoms_nh,
					 description_level);

	 }
      }
   }
   return comp_id;
}

// is_from_pdbx_model_atom is a optional argument bool false default
//
// return the number of atoms.
int 
coot::protein_geometry::comp_atom(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc,
				  bool is_from_pdbx_model_atom) {

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
   //
   std::string model_id; // pdbx_model_atom
   int ordinal_id; // pdbx_model_atom
   
   std::string comp_id_for_partial_charges = "unset"; // unassigned.

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

      // modify a reference (ierr)
      // 
      char *s = mmCIFLoop->GetString("comp_id",j,ierr);
      std::string atom_id;
      std::string type_symbol; 
      std::string type_energy = "unset";
      std::pair<bool, mmdb::realtype> partial_charge(0,0);

      std::string model_id; // pdbx_model_atom
      int ierr_pdbx = 0;
      char *pdbx_s = mmCIFLoop->GetString("model_id",j,ierr_pdbx);
      if (pdbx_s) {
	 model_id = pdbx_s;
      }
      ordinal_id = -1; // unset
      int ierr_pdbx_2 = mmCIFLoop->GetInteger(ordinal_id, "ordinal_id", j);

      std::pair<bool, int> pdbx_align(0, 0);
      int xalign;
      int pdbx_charge;
      int ierr_optional;
      int ierr_stereo_config;
      
      std::pair<bool, std::string> pdbx_leaving_atom_flag(false, "");
      std::pair<bool, std::string> pdbx_stereo_config_flag(false, "");
      std::pair<bool, int> formal_charge(false, 0); // read from PDB cif _chem_comp_atom.charge
      std::pair<bool, clipper::Coord_orth> pdbx_model_Cartn_ideal;
      std::pair<bool, clipper::Coord_orth> model_Cartn;
      // for cleanliness in debugging output
      pdbx_model_Cartn_ideal.second = clipper::Coord_orth(-1, -1, -1);
      model_Cartn.second            = clipper::Coord_orth(-1, -1, -1);
      dict_atom::aromaticity_t aromaticity = dict_atom::UNASSIGNED;

      if (ierr == 0 || (is_from_pdbx_model_atom && (ierr_pdbx_2 == 0))) {
	 int ierr_tot = 0;
	 if (s)
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

	 ierr_optional = mmCIFLoop->GetInteger(pdbx_charge, "charge", j);
	 if (! ierr_optional)
	    formal_charge = std::pair<bool, int> (true, pdbx_charge);

	 ierr_optional = mmCIFLoop->GetInteger(xalign, "pdbx_align", j);
	 if (! ierr_optional)
	    pdbx_align = std::pair<bool, int> (1, xalign);

	 s = mmCIFLoop->GetString("pdbx_aromatic_flag", j, ierr_optional);
	 if (s) {
	    if (! ierr_optional) {
	       std::string ss(s);
	       if (ss == "Y" || ss == "y")
		  aromaticity = dict_atom::AROMATIC;
	       if (ss == "N" || ss == "n")
		  aromaticity = dict_atom::NON_AROMATIC;
	    }
	 }

	 // just in case someone marks their aromatic atoms in this way.
	 s = mmCIFLoop->GetString("aromaticity", j, ierr_optional);
	 if (s) {
	    if (! ierr_optional) {
	       std::string ss(s);
	       if (ss == "Y" || ss == "y")
		  aromaticity = dict_atom::AROMATIC;
	    }
	 }
	 
	 s = mmCIFLoop->GetString("pdbx_leaving_atom_flag", j, ierr_optional);
	 if (s) {
	    if (! ierr_optional) 
	       pdbx_leaving_atom_flag = std::pair<bool, std::string> (true, s);
	 }


	 s = mmCIFLoop->GetString("pdbx_stereo_config", j, ierr_stereo_config);
	 if (s) {
	    if (! ierr_stereo_config)
	       pdbx_stereo_config_flag = std::pair<bool, std::string> (true, s);
	 }

	 mmdb::realtype x,y,z;
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
		  model_Cartn = std::pair<bool, clipper::Coord_orth>(1, clipper::Coord_orth(x,y,z));
		  if (close_float_p(x, 0.0))
		     if (close_float_p(z, 0.0))
			if (close_float_p(z, 0.0))
			   n_origin_model_atoms++;
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
		     // model_Cartn = std::pair<bool, clipper::Coord_orth>(true, clipper::Coord_orth(x,y,z));
		     model_Cartn.first = true;
		     model_Cartn.second = clipper::Coord_orth(x,y,z);
		     if (close_float_p(x, 0.0))
			if (close_float_p(z, 0.0))
			   if (close_float_p(z, 0.0))
			      n_origin_model_atoms++;
		  }
	 }

	 // It's possible that this data type is not in the cif file,
	 // so don't fail if we can't read it.

	 mmdb::realtype tmp_var;
	 ierr = mmCIFLoop->GetReal(tmp_var, "partial_charge", j);
	 if (ierr == 0) {
	    partial_charge = std::pair<bool, float>(1, tmp_var);
	    n_atoms_with_partial_charge++;
	 }

	 // ierr_tot will not be 0 for pdbx model atoms
	 // ...
	 if (ierr_tot == 0 || is_from_pdbx_model_atom) {

	    std::string padded_name = comp_atom_pad_atom_name(atom_id, type_symbol);
//  	    std::cout << "comp_atom_pad_atom_name: in :" << atom_id << ": out :"
//  		      << padded_name << ":" << std::endl;
	    n_atoms++;
	    if (comp_id_for_partial_charges != "bad match") { 
	       if (comp_id_for_partial_charges == "unset") {
		  comp_id_for_partial_charges = comp_id;
	       } else {
		  if (comp_id != comp_id_for_partial_charges) {
		     comp_id_for_partial_charges = "bad match";
		  }
	       }
	    }

	    if (is_from_pdbx_model_atom)
	       if (! model_id.empty())
		  comp_id = model_id;  // e.g. M_010_00001

	    if (false)
	       std::cout << "debug:: calling mon_lib_add_atom: "
			 << ":" << comp_id << ":  "
			 << ":" << atom_id << ":  "
			 << ":" << padded_name << ":  "
			 << ":" << type_symbol << ":  "
			 << "stereo-config: " << pdbx_stereo_config_flag.first
			 << " " << pdbx_stereo_config_flag.second << " "
			 << "model-pos " << model_Cartn.first << " " << model_Cartn.second.format() << " "
			 << "ideal-pos " << pdbx_model_Cartn_ideal.first << " "
			 << pdbx_model_Cartn_ideal.second.format()
			 << std::endl;

	    dict_atom atom(atom_id, padded_name, type_symbol, type_energy, partial_charge);
	    atom.aromaticity = aromaticity;
	    atom.formal_charge = formal_charge;
	    atom.pdbx_stereo_config = pdbx_stereo_config_flag;
	    if (model_Cartn.first)
	       atom.add_pos(dict_atom::REAL_MODEL_POS, model_Cartn);

	    if (pdbx_model_Cartn_ideal.first)
	       atom.add_pos(dict_atom::IDEAL_MODEL_POS, pdbx_model_Cartn_ideal);

	    atom.formal_charge      = formal_charge;
	    atom.aromaticity        = aromaticity;
	    atom.pdbx_stereo_config = pdbx_stereo_config_flag;

	    if (is_from_pdbx_model_atom)
	       if (ierr_pdbx == 0)
		  if (ierr_pdbx_2 == 0)
		     atom.add_ordinal_id(ordinal_id);

	    mon_lib_add_atom(comp_id, imol_enc, atom);

	 } else {
	    std::cout << " error on read " << ierr_tot << std::endl;
	 } 
      }
   }

   if (n_atoms_with_partial_charge == n_atoms) {
      if (comp_id_for_partial_charges != "unset") {
	 if (comp_id_for_partial_charges != "bad match") {
	    for (unsigned int id=0; id<dict_res_restraints.size(); id++) {
	       if (dict_res_restraints[id].second.residue_info.comp_id == comp_id_for_partial_charges) {
		  if (dict_res_restraints[id].first == imol_enc) {
		     dict_res_restraints[id].second.set_has_partial_charges(1);
		  }
	       } 
	    }
	 }
      }
   }

   // Now we need to check that the atom ideal or model coordinates were not at the origin.
   // 
   if (n_origin_ideal_atoms > 2) // trash all ideal atoms
      delete_atom_positions(comp_id, imol_enc, coot::dict_atom::IDEAL_MODEL_POS);
   if (n_origin_model_atoms > 2) // trash all model/real atoms
      delete_atom_positions(comp_id, imol_enc, coot::dict_atom::REAL_MODEL_POS);
   
   return n_atoms;
}


void
coot::protein_geometry::comp_tree(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc) {

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
	 mon_lib_add_tree(comp_id, imol_enc, padded_name_atom_id, padded_name_atom_back,
			  padded_name_atom_forward, connect_type); 
      }
   }
}

// look up atom_id in the atom atom_info (dict_atom vector) of the comp_id restraint
// 
std::string
coot::protein_geometry::atom_name_for_tree_4c(const std::string &comp_id, const std::string &atom_id) const {

   std::string r = atom_id;

   if (dict_res_restraints.size() > 0) { 
      for (int id=(dict_res_restraints.size()-1); id >=0; id--) {
	 if (dict_res_restraints[id].second.residue_info.comp_id == comp_id) {
	    r = dict_res_restraints[id].second.atom_name_for_tree_4c(atom_id);
	    break;
	 }
      }
   }
   return r;
}


int
coot::protein_geometry::comp_bond(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc, bool is_from_pdbx_model_bond) {

   bool verbose_output = 0; // can be passed, perhaps.
   std::string comp_id;
   std::string atom_id_1, atom_id_2;
   std::string type;
   mmdb::realtype value_dist = -1.0, value_dist_esd = -1.0;
   std::string model_id; // pdbx_chem_comp_model_bond

   char *s; 
   int nbond = 0;
   int comp_id_index = -1; // not found initially

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) { 

      int ierr;
      int ierr_tot = 0;
      int ierr_optional = 0;
      int ierr_tot_for_ccd = 0;  // CCD comp_bond entries don't have
			         // target geometry (dist and esd) but we
			         // want to be able to read them in
			         // anyway (to get the bond orders for
			         // drawing).
      dict_bond_restraint_t::aromaticity_t aromaticity(dict_bond_restraint_t::UNASSIGNED);

   
      // modify a reference (ierr)
      //

      s = mmCIFLoop->GetString("comp_id",j,ierr);
      ierr_tot += ierr;
      ierr_tot_for_ccd += ierr;
      if (s) { 
	 comp_id = s;
	 for (int id=(dict_res_restraints.size()-1); id >=0; id--) {
	    if (dict_res_restraints[id].second.residue_info.comp_id == comp_id) {
	       comp_id_index = id;
	       break;
	    }
	 }
      }

      if (is_from_pdbx_model_bond) {
	 int ierr_pdbx = 0;
	 s = mmCIFLoop->GetString("model_id",j,ierr_pdbx);
	 if (! ierr_pdbx) {
	    if (s) {
	       model_id = s;
	       for (int id=(dict_res_restraints.size()-1); id >=0; id--) {
		  if (dict_res_restraints[id].second.residue_info.comp_id == model_id) {
		     comp_id_index = id;
		     break;
		  }
	       }
	    } else {
	       std::cout << "oops! null s" << std::endl;
	    }
	 } else {
	    std::cout << "oops ierr_pdbx is not 0" << std::endl;
	 }
      }

      if (comp_id_index == -1) {

	 std::cout << "WARNING:: failed to find dictionary entry index for "
		   << comp_id << std::endl;

      } else {

	 s = mmCIFLoop->GetString("atom_id_1", j, ierr);
	 ierr_tot += ierr;
	 ierr_tot_for_ccd += ierr;
	 if (s) 
	    atom_id_1 = get_padded_name(s, comp_id_index);

	 s = mmCIFLoop->GetString("atom_id_2", j, ierr);
	 ierr_tot += ierr;
	 ierr_tot_for_ccd += ierr;
	 if (s) 
	    atom_id_2 = get_padded_name(s, comp_id_index);

	 s = mmCIFLoop->GetString("type", j, ierr);
	 if (s) 
	    type = s;
	 // perhaps it was in the dictionary as "value_order"?
	 if (ierr) {
	    s = mmCIFLoop->GetString("value_order", j, ierr);
	 }

	 if (! ierr) {
	    if (s) { // just in case (should not be needed).
	       std::string ss(s);
	       // convert from Chemical Component Dictionary value_order to
	       // refmac monomer library chem_comp_bond types - 
	       // or FeiDrg output.
	       if (ss == "SING")
		  type = "single";
	       if (ss == "DOUB")
		  type = "double";
	       if (ss == "TRIP")
		  type = "triple";
	       if (ss == "trip") // acedrg
		  type = "triple";
	       if (ss == "TRIPLE")
		  type = "triple";

	       if (ss == "SINGLE")
		  type = "single";
	       if (ss == "DOUBLE")
		  type = "double";
	       if (ss == "DELOC")
		  type = "deloc";
	       
	       // Chemical Chemical Dictionary also has an aromatic flag, so
	       // we can have bonds that are (for example) "double"
	       // "aromatic" Hmm!  Food for thought.
	       
	       // Metal bonds are "SING" (i.e. CCD doesn't have metal
	       // bonds).
	    }
	 }

	 // Acedrg marks aromatic bonds in this way.
	 // 
	 s = mmCIFLoop->GetString("aromaticity", j, ierr_optional);
	 if (s) {
	    if (! ierr_optional) {
	       std::string ss(s);
	       if (ss == "Y" || ss == "y")
		  aromaticity = dict_bond_restraint_t::AROMATIC;
	       if (ss == "N" || ss == "n")
		  aromaticity = dict_bond_restraint_t::NON_AROMATIC;
	    }
	 }

	 // PDB marks aromaticity this way:
	 // 
	 s = mmCIFLoop->GetString("pdbx_aromatic_flag", j, ierr_optional);
	 if (s) {
	    if (! ierr_optional) {
	       std::string ss(s);
	       if (ss == "Y" || ss == "y")
		  aromaticity = dict_bond_restraint_t::AROMATIC;
	       if (ss == "N" || ss == "n")
		  aromaticity = dict_bond_restraint_t::NON_AROMATIC;
	    }
	 }

	 
	 ierr_tot += ierr;
	 ierr_tot_for_ccd += ierr;

	 ierr = mmCIFLoop->GetReal(value_dist, "value_dist", j);
	 ierr_tot += ierr;

	 ierr = mmCIFLoop->GetReal(value_dist_esd, "value_dist_esd", j);
	 ierr_tot += ierr;

	 if (ierr_tot == 0) {

	    mon_lib_add_bond(comp_id, imol_enc, atom_id_1, atom_id_2,
			     type, value_dist, value_dist_esd, aromaticity); 
	    nbond++;
	 } else {

	    if (! ierr_tot_for_ccd || is_from_pdbx_model_bond) {

	       if (is_from_pdbx_model_bond)
		  comp_id = model_id;

	       mon_lib_add_bond_no_target_geom(comp_id, imol_enc, atom_id_1, atom_id_2, type, aromaticity);
	       nbond++;
	    
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
   }

   return nbond;
}
void
coot::protein_geometry::comp_angle(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc) {

   std::string comp_id;
   std::string atom_id_1, atom_id_2, atom_id_3;
   mmdb::realtype value_angle, value_angle_esd;
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
	    if (dict_res_restraints[id].second.residue_info.comp_id == comp_id) {
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
	 mon_lib_add_angle(comp_id, imol_enc,
			   atom_id_1, atom_id_2, atom_id_3,
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
coot::protein_geometry::comp_torsion(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc) {

   std::string comp_id, id;
   std::string atom_id_1, atom_id_2, atom_id_3, atom_id_4;
   mmdb::realtype value_angle, value_angle_esd;
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
	    if (dict_res_restraints[id].second.residue_info.comp_id == comp_id) {
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
   
	 mon_lib_add_torsion(comp_id, imol_enc,
			     id,
			     atom_id_1, atom_id_2,
			     atom_id_3, atom_id_4,
			     value_angle, value_angle_esd, period); 
      }
   }
} 


std::pair<int, std::vector<std::string> >
coot::protein_geometry::comp_chiral(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc) {

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
	    if (dict_res_restraints[id].second.residue_info.comp_id == comp_id) {
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
	 mon_lib_add_chiral(comp_id, imol_enc,
			    id, atom_id_centre,
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
coot::protein_geometry::comp_plane(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc) {

   std::string comp_id;
   std::string atom_id, plane_id;
   mmdb::realtype dist_esd; 

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
	    if (dict_res_restraints[id].second.residue_info.comp_id == comp_id) {
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
	 mon_lib_add_plane(comp_id, imol_enc, plane_id, atom_id, dist_esd);
      } else {
	 std::cout << "problem reading comp plane" << std::endl;
      } 
   }
} 
// 
int
coot::protein_geometry::add_chem_mods(mmdb::mmcif::PData data) {

   int n_mods = 0;
   for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 
      
      mmdb::mmcif::PCategory cat = data->GetCategory(icat);
      std::string cat_name(cat->GetCategoryName());
      mmdb::mmcif::PLoop mmCIFLoop = data->GetLoop(cat_name.c_str() );
            
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
coot::protein_geometry::add_synonyms(mmdb::mmcif::PData data) {

   for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 
      
      mmdb::mmcif::PCategory cat = data->GetCategory(icat);
      std::string cat_name(cat->GetCategoryName());
      mmdb::mmcif::PLoop mmCIFLoop = data->GetLoop(cat_name.c_str() );
            
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
coot::protein_geometry::add_chem_comp_synonym(mmdb::mmcif::PLoop mmCIFLoop) {

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
coot::protein_geometry::add_chem_mod(mmdb::mmcif::PLoop mmCIFLoop) {

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
   bool found = false;
   bool debug = false;
   
   if (debug) {
      std::cout << "---------------------- Here are the chem_links: -----------------"
		<< std::endl;
      print_chem_links();
   }

   // This needs to iterate to make the count now that we use a map
   // std::cout << "Testing vs " << chem_link_vec.size() << " chem links\n";

   unsigned int search_hash_code_f = chem_link::make_hash_code(comp_id_1, comp_id_2, group_1, group_2);
   unsigned int search_hash_code_b = chem_link::make_hash_code(comp_id_2, comp_id_1, group_2, group_1);

   if (debug)
      std::cout << "here in matching_chem_link() " << search_hash_code_f << " " << search_hash_code_b << " "
		<< comp_id_1 << " " << comp_id_2 << " " << group_1 << " " << group_2 << std::endl;

   // Is this link a TRANS peptide or a CIS?  Both have same group and
   // comp_ids.  Similarly, is is BETA-1-2 or BETA1-4 (etc).  We need
   // to decide later, don't just pick the first one that matches
   // (keep the order switch flag too).
   //

   // "gap" and "symmetry" have hash code 0 (blank strings)

   std::vector<std::pair<coot::chem_link, bool> > matching_chem_links;
   std::map<unsigned int, std::vector<chem_link> >::const_iterator it =
      chem_link_map.find(search_hash_code_f);
   if (it == chem_link_map.end()) {
      it = chem_link_map.find(search_hash_code_b);
      if (it != chem_link_map.end())
	 switch_order_flag = true; // used?
   }

   if (debug) {
      if (it != chem_link_map.end())
	 std::cout << "matching_chem_link() found the hash at least! " << std::endl;
      else
	 std::cout << "matching_chem_link() failed to find hash " << search_hash_code_f << " "
		   << search_hash_code_b << std::endl;
   }

   if (it == chem_link_map.end()) {

      // NAG-ASN for example

      unsigned int search_bl_1_f = chem_link::make_hash_code(comp_id_1, comp_id_2, "", group_2);
      unsigned int search_bl_1_b = chem_link::make_hash_code(comp_id_2, comp_id_1, group_2, "");
      unsigned int search_bl_2_f = chem_link::make_hash_code(comp_id_1, comp_id_2, group_1, "");
      unsigned int search_bl_2_b = chem_link::make_hash_code(comp_id_2, comp_id_1, "", group_1);

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

      // std::cout << "-------- here with candidate_chem_links size ------- "
      // << candidate_chem_links.size() << std::endl;

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

		  // std::cout << "::::::::: adding matching chem link " << cl << std::endl;
		  matching_chem_links.push_back(p);
	       }
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

	 std::pair<bool, bool> match_res =
	    cl.matches_comp_ids_and_groups(comp_id_1, group_1, comp_id_2, group_2);

	 if (debug)
	    std::cout << "... matching_chem_link: found matching link "
		      << comp_id_1 << " " << comp_id_2 << " " 
		      << cl << std::endl;

	 if (debug)
	    std::cout << "   checking chem link: " << cl << " -> "
		      << match_res.first << " " << match_res.second << std::endl;

	 if (match_res.first) {
	    if (cl.Id() != "gap" && cl.Id() != "symmetry") {
	       if (!cl.is_peptide_link_p() || allow_peptide_link_flag) {
		  switch_order_flag = match_res.second;
		  found = true;
		  std::pair<coot::chem_link, bool> p(cl, switch_order_flag);
		  matching_chem_links.push_back(p);

	       } else {
		  if (debug)
		     std::cout << "reject link on peptide/allow-peptide test " << std::endl;
	       }
	    } else {
	       if (debug) {
		  std::cout << "reject link \"" << cl.Id() << "\"" << std::endl;
	       }
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
   if (debug)
      std::cout << "matching_chem_link() returns " << matching_chem_links.size()
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

   // struct stat buf;
   char *cmld = NULL;
   
   char *s = getenv("COOT_REFMAC_LIB_DIR");
   if (s) {
      if (! is_dir_or_link(s)) { 
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

	 std::string ss(s); // might have trailing "/"
	 ss = coot::util::remove_trailing_slash(ss);
	 if (! is_dir_or_link(ss)) { 
	    env_dir_fails = 1;
	 } else {
	    env_dir_fails = 0;
	    std::cout << "INFO:: Using Standard CCP4 Refmac dictionary from"
		      << " CLIBD_MON: " << s << std::endl;
	    mon_lib_dir = s;
	    using_clibd_mon = 1;
	    // strip any trailing / from mon_lib_dir
	    if (mon_lib_dir.length() > 0) {
	       if (mon_lib_dir.at(mon_lib_dir.length()-1) == '/')
		  mon_lib_dir = mon_lib_dir.substr(0,mon_lib_dir.length()-1);
	    }
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

	    if (is_dir_or_link(hardwired_default_place)) {
	       mon_lib_dir = hardwired_default_place;
	    } else {

	       // OK, let's look for $COOT_PREFIX/share/coot/lib (as you
	       // would with the binary distros)

	       s = getenv("COOT_PREFIX");
	       if (s) {
		  std::string lib_dir = util::append_dir_dir(s, "share");
		  lib_dir = util::append_dir_dir(lib_dir, "coot");
		  lib_dir = util::append_dir_dir(lib_dir, "lib");
		  if (is_dir_or_link(lib_dir)) {
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
      mon_lib_dir =  coot::util::intelligent_debackslash(mon_lib_dir);
      std::string filename = mon_lib_dir;
      // contains the linkages:
      filename += "/data/monomers/list/mon_lib_list.cif";
      if (using_clibd_mon) {
	 filename = util::remove_trailing_slash(mon_lib_dir);
	 filename += "/list/mon_lib_list.cif";
      }
      // now check that that file is there:
      if (! is_regular_file(filename)) {
	 std::cout << "ERROR: dictionary " << filename << " is not a regular file"
		   << std::endl;
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
   
   init_refmac_mon_lib(mon_lib_cif, protein_geometry::MON_LIB_LIST_CIF);
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
   
   int imol_enc = IMOL_ENC_ANY; // maybe pass this?
   
   std::string filename = util::append_dir_file(s, protein_mono);
   if (is_regular_file(filename)) {
      init_refmac_mon_lib(filename, read_number, imol_enc);
      read_number++;
   } else {
      std::cout << "WARNING: file " << filename << " is not a regular file"
		<< std::endl;
   }
   return read_number;
}

// quote atom name as needed - i.e. CA -> CA, CA' -> "CA'"
std::string 
coot::dictionary_residue_restraints_t::quoted_atom_name(const std::string &an) const {

   std::string n = an;
   bool has_quotes = false;

   for (unsigned int i=0; i<an.size(); i++) {
      if (an[i] == '\'') {
	 has_quotes = true;
	 break;
      } 
   }
   if (has_quotes)
      n = "\"" + an + "\"";

   return n;
}

void
coot::dictionary_residue_restraints_t::write_cif(const std::string &filename) const {

   mmdb::mmcif::File *mmCIFFile = new mmdb::mmcif::File(); // d
      
   mmdb::mmcif::Data   *mmCIFData = NULL;
   mmdb::mmcif::Struct *mmCIFStruct;
   char S[2000];
   
   //  2.1  Example 1: add a structure into mmCIF object

   int rc;

   rc = mmCIFFile->AddCIFData("comp_list");
   mmCIFData = mmCIFFile->GetCIFData("comp_list");
   rc = mmCIFData->AddStructure ("_chem_comp", mmCIFStruct);

   if (rc!=mmdb::mmcif::CIFRC_Ok && rc!=mmdb::mmcif::CIFRC_Created)  {
      // badness!
      std::cout << "rc not mmdb::mmcif::CIFRC_Ok " << rc << std::endl;
      printf ( " **** error: attempt to retrieve Loop as a Structure.\n" );
      if (!mmCIFStruct)  {
	 printf ( " **** error: mmCIFStruct is NULL - report as a bug\n" );
      }
   } else {
      if (rc == mmdb::mmcif::CIFRC_Created) {
	 // printf ( " -- new structure created\n" );
      } else { 
	 printf(" -- structure was already in mmCIF, it will be extended\n");
      }
      // std::cout << "SUMMARY:: rc mmdb::mmcif::CIFRC_Ok or newly created. " << std::endl;

      mmdb::mmcif::Loop *mmCIFLoop = new mmdb::mmcif::Loop; // 20100212
      // data_comp_list, id, three_letter_code, name group etc:

      rc = mmCIFData->AddLoop("_chem_comp", mmCIFLoop);
      int i=0;
      const char *s = residue_info.comp_id.c_str();
      mmCIFLoop->PutString(s, "id", i);
      s = residue_info.three_letter_code.c_str();
      mmCIFLoop->PutString(s, "three_letter_code", i);
      std::string raw_name = residue_info.name.c_str();
      std::string quoted_name = util::single_quote(raw_name, "\"");
      mmCIFLoop->PutString(raw_name.c_str(), "name", i);
      s =  residue_info.group.c_str();
      mmCIFLoop->PutString(s, "group", i);
      int nat = residue_info.number_atoms_all;
      mmCIFLoop->PutInteger(nat, "number_atoms_all", i);
      nat = number_of_non_hydrogen_atoms();
      mmCIFLoop->PutInteger(nat, "number_atoms_nh", i);
      s = residue_info.description_level.c_str();
      mmCIFLoop->PutString(s, "desc_level", i);
      
      std::string comp_record = "comp_list";
      mmCIFData->PutDataName(comp_record.c_str()); // 'data_' record

      // atom loop

      std::string comp_monomer_name = "comp_";
      comp_monomer_name += residue_info.comp_id; 
      rc = mmCIFFile->AddCIFData(comp_monomer_name.c_str());
      mmCIFData = mmCIFFile->GetCIFData(comp_monomer_name.c_str());

      // shall we add coordinates too?
      //
      bool add_coordinates = true; // if no atoms have coords, this gets set to false.
      int n_atoms_with_coords = 0;
      for (unsigned int i=0; i<atom_info.size(); i++) {
	 if (atom_info[i].model_Cartn.first)
	    n_atoms_with_coords++;
      }
      
      if (n_atoms_with_coords == 0)
	 add_coordinates = false;
      
      if (atom_info.size()) { 
	 rc = mmCIFData->AddLoop("_chem_comp_atom", mmCIFLoop);
	 if (rc == mmdb::mmcif::CIFRC_Ok || rc == mmdb::mmcif::CIFRC_Created) {
	    for (unsigned int i=0; i<atom_info.size(); i++) {
	       const dict_atom &ai = atom_info[i];
	       const char *ss =  residue_info.comp_id.c_str();
	       mmCIFLoop->PutString(ss, "comp_id", i);
	       std::string annw = util::remove_whitespace(ai.atom_id).c_str();
	       std::string qan = quoted_atom_name(annw);
	       mmCIFLoop->PutString(annw.c_str(), "atom_id", i);
	       std::string up_type_symbol = util::upcase(atom_info[i].type_symbol);
	       ss = up_type_symbol.c_str();
	       mmCIFLoop->PutString(ss, "type_symbol", i);
	       // std::cout << "up_type_symbol: " << up_type_symbol << std::endl;
	       std::string up_type_energy = util::upcase(atom_info[i].type_energy);
	       // std::cout << "up_type_energy: " << up_type_energy << std::endl;
	       ss = up_type_energy.c_str();
	       mmCIFLoop->PutString(ss, "type_energy", i);
	       if (atom_info[i].partial_charge.first) {
		  float v = atom_info[i].partial_charge.second;
		  mmCIFLoop->PutReal(v, "partial_charge", i, 4);
	       }
	       if (add_coordinates) {
		  if (atom_info[i].model_Cartn.first) { 
		     float x = atom_info[i].model_Cartn.second.x();
		     float y = atom_info[i].model_Cartn.second.y();
		     float z = atom_info[i].model_Cartn.second.z();
		     mmCIFLoop->PutReal(x, "x", i, 6);
		     mmCIFLoop->PutReal(y, "y", i, 6);
		     mmCIFLoop->PutReal(z, "z", i, 6);
		  }
	       } 
	    }
	 }
      }

      // bond loop

      if (bond_restraint.size()) { 
	 rc = mmCIFData->AddLoop("_chem_comp_bond", mmCIFLoop);
	 if (rc == mmdb::mmcif::CIFRC_Ok || rc == mmdb::mmcif::CIFRC_Created) {
	    // std::cout << " number of bonds: " << bond_restraint.size() << std::endl;
	    for (unsigned int i=0; i<bond_restraint.size(); i++) {
	       // std::cout << "ading bond number " << i << std::endl;
	       const char *ss = residue_info.comp_id.c_str();
	       mmCIFLoop->PutString(ss, "comp_id", i);
	       std::string id_1 = util::remove_whitespace(bond_restraint[i].atom_id_1_4c());
	       std::string id_2 = util::remove_whitespace(bond_restraint[i].atom_id_2_4c());
	       std::string qan_1 = quoted_atom_name(id_1);
	       std::string qan_2 = quoted_atom_name(id_2);
	       ss = id_1.c_str();
	       mmCIFLoop->PutString(id_1.c_str(), "atom_id_1", i);
	       ss = id_2.c_str();
	       mmCIFLoop->PutString(id_2.c_str(), "atom_id_2", i);
	       ss = bond_restraint[i].type().c_str();
	       mmCIFLoop->PutString(ss, "type", i);
	       try {
	          float v = bond_restraint[i].value_dist();
		  mmCIFLoop->PutReal(v, "value_dist", i, 5);
		  v = bond_restraint[i].value_esd();
		  mmCIFLoop->PutReal(v, "value_dist_esd", i, 3);
	       }
	       catch (const std::runtime_error &rte) {
		  // do nothing, it's not really an error if the dictionary
		  // doesn't have target geometry (the bonding description came
		  // from a Chemical Component Dictionary entry for example).
	       }
	    }
	 }
      }

      // angle loop

      if (angle_restraint.size()) { 
	 rc = mmCIFData->AddLoop("_chem_comp_angle", mmCIFLoop);
	 if (rc == mmdb::mmcif::CIFRC_Ok || rc == mmdb::mmcif::CIFRC_Created) {
	    // std::cout << " number of angles: " << angle_restraint.size() << std::endl;
	    for (unsigned int i=0; i<angle_restraint.size(); i++) {
	       // std::cout << "ading angle number " << i << std::endl;
	       std::string id_1 = util::remove_whitespace(angle_restraint[i].atom_id_1_4c());
	       std::string id_2 = util::remove_whitespace(angle_restraint[i].atom_id_2_4c());
	       std::string qan_1 = quoted_atom_name(id_1);
	       std::string qan_2 = quoted_atom_name(id_2);
	       const char *ss = residue_info.comp_id.c_str();
	       mmCIFLoop->PutString(ss, "comp_id", i);
	       ss = id_1.c_str();
	       mmCIFLoop->PutString(id_1.c_str(), "atom_id_1", i);
	       ss = id_2.c_str();
	       mmCIFLoop->PutString(id_2.c_str(), "atom_id_2", i);

	       // bug fix(!) intermediate value my_ss clears up casting
	       // problem (inheritance-related?) on writing.
	       std::string id_3 = util::remove_whitespace(angle_restraint[i].atom_id_3_4c());
	       // 20170305 std::string qan_3 = quoted_atom_name(my_ss);
	       // ss = qan_3.c_str();
	       mmCIFLoop->PutString(id_3.c_str(), "atom_id_3", i);

	       float v = angle_restraint[i].angle();
	       mmCIFLoop->PutReal(v, "value_angle", i, 5);
	       v = angle_restraint[i].esd();
	       mmCIFLoop->PutReal(v, "value_angle_esd", i, 3);
	    }
	 }
      }

      // torsion loop

      if (torsion_restraint.size() > 0) { 
	 rc = mmCIFData->AddLoop("_chem_comp_tor", mmCIFLoop);
	 if (rc == mmdb::mmcif::CIFRC_Ok || rc == mmdb::mmcif::CIFRC_Created) {
	    // std::cout << " number of torsions: " << torsion_restraint.size() << std::endl;
	    for (unsigned int i=0; i<torsion_restraint.size(); i++) {
	       // std::cout << "ading torsion number " << i << std::endl;
	       std::string id_1 = util::remove_whitespace(torsion_restraint[i].atom_id_1_4c());
	       std::string id_2 = util::remove_whitespace(torsion_restraint[i].atom_id_2_4c());
	       std::string id_3 = util::remove_whitespace(torsion_restraint[i].atom_id_3_4c());
	       std::string id_4 = util::remove_whitespace(torsion_restraint[i].atom_id_4_4c());
	       std::string qan_1 = quoted_atom_name(id_1);
	       std::string qan_2 = quoted_atom_name(id_2);
	       std::string qan_3 = quoted_atom_name(id_3);
	       std::string qan_4 = quoted_atom_name(id_4);
	       const char *ss = residue_info.comp_id.c_str();
	       mmCIFLoop->PutString(ss, "comp_id", i);
	       ss = torsion_restraint[i].id().c_str();
	       mmCIFLoop->PutString(ss, "id", i);
	       ss = id_1.c_str();
	       mmCIFLoop->PutString(id_1.c_str(), "atom_id_1", i);
	       ss = id_2.c_str();
	       mmCIFLoop->PutString(id_2.c_str(), "atom_id_2", i);
	       ss = id_3.c_str();
	       mmCIFLoop->PutString(id_3.c_str(), "atom_id_3", i);
	       ss = id_4.c_str();
	       mmCIFLoop->PutString(id_4.c_str(), "atom_id_4", i);
	       float v = torsion_restraint[i].angle();
	       mmCIFLoop->PutReal(v, "value_angle", i, 5);
	       v = torsion_restraint[i].esd();
	       mmCIFLoop->PutReal(v, "value_angle_esd", i, 3);
	       int p = torsion_restraint[i].periodicity();
	       mmCIFLoop->PutInteger(p, "period", i);
	    }
	 }
      }

      // chiral loop
      // 
      if (chiral_restraint.size() > 0) { 
	 rc = mmCIFData->AddLoop("_chem_comp_chir", mmCIFLoop);
	 if (rc == mmdb::mmcif::CIFRC_Ok || rc == mmdb::mmcif::CIFRC_Created) {
	    // std::cout << " number of chirals: " << chiral_restraint.size() << std::endl;
	    for (unsigned int i=0; i<chiral_restraint.size(); i++) {
	       // std::cout << "ading chiral number " << i << std::endl;
	       const char *ss = residue_info.comp_id.c_str();
	       std::string id_c = util::remove_whitespace(chiral_restraint[i].atom_id_c_4c());
	       std::string id_1 = util::remove_whitespace(chiral_restraint[i].atom_id_1_4c());
	       std::string id_2 = util::remove_whitespace(chiral_restraint[i].atom_id_2_4c());
	       std::string id_3 = util::remove_whitespace(chiral_restraint[i].atom_id_3_4c());
	       std::string qan_c = quoted_atom_name(id_c);
	       std::string qan_1 = quoted_atom_name(id_1);
	       std::string qan_2 = quoted_atom_name(id_2);
	       std::string qan_3 = quoted_atom_name(id_3);
	       mmCIFLoop->PutString(ss, "comp_id", i);
	       ss = chiral_restraint[i].Chiral_Id().c_str();
	       mmCIFLoop->PutString(ss, "id", i);
	       ss = id_c.c_str();
	       mmCIFLoop->PutString(id_c.c_str(), "atom_id_centre", i);
	       ss = id_1.c_str();
	       mmCIFLoop->PutString(id_1.c_str(), "atom_id_1", i);
	       ss = id_2.c_str();
	       mmCIFLoop->PutString(id_2.c_str(), "atom_id_2", i);
	       ss = id_3.c_str();
	       mmCIFLoop->PutString(id_3.c_str(), "atom_id_3", i);
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
	 if (rc == mmdb::mmcif::CIFRC_Ok || rc == mmdb::mmcif::CIFRC_Created) {
	    // std::cout << " number of planes: " << plane_restraint.size() << std::endl;
	    int icount = 0;
	    for (unsigned int i=0; i<plane_restraint.size(); i++) {
	       // std::cout << "DEBUG:: adding plane number " << i << std::endl;
	       for (int iat=0; iat<plane_restraint[i].n_atoms(); iat++) {
		  const char *ss = residue_info.comp_id.c_str();
		  mmCIFLoop->PutString(ss, "comp_id", icount);
		  ss = plane_restraint[i].plane_id.c_str();
		  mmCIFLoop->PutString(ss, "plane_id", icount);
		  
		  std::string id = util::remove_whitespace(plane_restraint[i].atom_id(iat));
		  std::string qan = quoted_atom_name(id);
		  ss = id.c_str();
		  mmCIFLoop->PutString(id.c_str(), "atom_id", icount);

		  float v = plane_restraint[i].dist_esd(iat);
		  mmCIFLoop->PutReal(v, "dist_esd", icount, 4);
		  icount++;
	       }
	    }
	 }
      }


      write_cif_pdbx_chem_comp_descriptor(mmCIFData);

      // delete mmCIFLoop; // crashed when enabled?
      
      int status = mmCIFFile->WriteMMCIFFile(filename.c_str());
      if (status == 0) 
	 std::cout << "INFO:: wrote mmCIF \"" << filename << "\"" << std::endl;
      else 
	 std::cout << "INFO:: on write mmCIF \"" << filename << "\" status: "
		   << status << std::endl;
   }
   delete mmCIFFile; // deletes all its attributes too.
}

void
coot::dictionary_residue_restraints_t::write_cif_pdbx_chem_comp_descriptor(mmdb::mmcif::Data *mmCIFData) const {

   

} 


// make a connect file specifying the bonds to Hydrogens
bool
coot::protein_geometry::hydrogens_connect_file(const std::string &resname,
					       const std::string &filename) const {

   bool r = 0;
   std::pair<short int, dictionary_residue_restraints_t> p =
      get_monomer_restraints(resname, IMOL_ENC_ANY);

   if (p.first) {
      std::vector<dict_bond_restraint_t> bv = p.second.bond_restraint;
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
   
   mmdb::mmcif::File ciffile;
   if (! is_regular_file(cif_dictionary_file_name)) {
      std::cout << "WARNIG:: cif dictionary " << cif_dictionary_file_name
		<< " not found" << std::endl;
   } else {
      int ierr = ciffile.ReadMMCIFFile(cif_dictionary_file_name.c_str());
      if (ierr != mmdb::mmcif::CIFRC_Ok) {
	 std::cout << "Dirty mmCIF file? " << cif_dictionary_file_name
		   << std::endl;
      } else {
	 for(int idata=0; idata<ciffile.GetNofData(); idata++) { 
         
	    mmdb::mmcif::PData data = ciffile.GetCIFData(idata);
	    for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 
	       mmdb::mmcif::PCategory cat = data->GetCategory(icat);
	       std::string cat_name(cat->GetCategoryName());
	       mmdb::mmcif::PLoop mmCIFLoop =
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

// maybe these function need their own file.  For now they can go here.
// 
void
coot::protein_geometry::pdbx_chem_comp_descriptor(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc) {

   std::string comp_id;
   std::string type;
   std::string program;
   std::string program_version;
   std::string descriptor;
   
   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {
      int ierr;
      int ierr_tot = 0;
      char *s;
      s = mmCIFLoop->GetString("comp_id",j,ierr);
      ierr_tot += ierr;
      if (s) comp_id = s;
      s = mmCIFLoop->GetString("program",j,ierr);
      if (s) program = s;
      s = mmCIFLoop->GetString("program_version",j,ierr);
      if (s) program_version = s;
      s = mmCIFLoop->GetString("descriptor",j,ierr);
      ierr_tot += ierr;
      if (s) descriptor = s;
      s = mmCIFLoop->GetString("type",j,ierr);
      ierr_tot += ierr;
      if (s) type = s;
      
      if (ierr_tot == 0) {
	 pdbx_chem_comp_descriptor_item descr(type, program, program_version, descriptor);
	 add_pdbx_descriptor(comp_id, imol_enc, descr);
      }
   }
}



// imol_enc will always be IMOL_ENC_ANY, surely
void
coot::protein_geometry::add_pdbx_descriptor(const std::string &comp_id,
					    int imol_enc,
					    pdbx_chem_comp_descriptor_item &descr) {

   // like the others, not using iterators because we don't have a
   // comparitor using a string (the comp_id).
   // 
   bool found = false;
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
	 found = true;
	 dict_res_restraints[i].second.descriptors.descriptors.push_back(descr);
	 break;
      }
   }
   if (! found) {
      dictionary_residue_restraints_t rest(comp_id, read_number);
      rest.descriptors.descriptors.push_back(descr);
      std::pair<int, dictionary_residue_restraints_t> p(imol_enc, rest);
      dict_res_restraints.push_back(p);
   } 
}
