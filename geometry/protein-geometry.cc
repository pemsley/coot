/* geometry/protein-geometry.cc
 * 
 * Copyright 2003, 2004, 2005, 2006 The University of York
 * Author: Paul Emsley
 * Copyright 2007, 2008, 2009, 2010, 2011, 2012 The University of Oxford
 * Copyright 2014, 2015, 2016 by Medical Research Council
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <algorithm>  // needed for sort? Yes.
#include <stdexcept>  // Thow execption.

#include "utils/win-compat.hh"
#include "mini-mol/atom-quads.hh"
#include "geometry/protein-geometry.hh"
#include "utils/coot-utils.hh"

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if !defined _MSC_VER
#include <unistd.h>
#else
#define XDATADIR "C:/coot/share"
#define PKGDATADIR XDATADIR
#define S_ISDIR(m)  (((m) & S_IFMT) == S_IFDIR)
#define S_ISREG(m)  (((m) & S_IFMT) == S_IFREG)
#endif

#include "clipper/core/clipper_util.h"

#include "compat/coot-sysdep.h"

#include "lbg-graph.hh"
#include "utils/xdg-base.hh"

#include "utils/logging.hh"
extern logging logger;

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

bool
coot::protein_geometry::matches_imol(int imol_dict, int imol_enc) const {

   if (imol_dict == IMOL_ENC_ANY)
      return true;
   else
      return (imol_dict == imol_enc);
}


void
coot::protein_geometry::assign_chiral_volume_targets() {

   for (unsigned int idict=0; idict<dict_res_restraints.size(); idict++) {
      if (dict_res_restraints[idict].second.has_unassigned_chiral_volumes()) {
	 if (false)
	    std::cout << "DEBUG:: assign_chiral_volume_targets for dict_res_restraints entry: "
		      << idict << " "
		      << dict_res_restraints[idict].second.residue_info.comp_id
		      << " has unassigned chiral volumes" << std::endl;
	 dict_res_restraints[idict].second.assign_chiral_volume_targets();
      }
   }

   assign_link_chiral_volume_targets();
}

// the chiral restraint for this comp_id(s) may need filtering
// (i.e. removing some of them if they are not real chiral centres
// (e.g. from prodrg restraints)).
void
coot::protein_geometry::filter_chiral_centres(int imol, const std::vector<std::string> &comp_ids_for_filtering) {

   for (unsigned int ichir=0; ichir<comp_ids_for_filtering.size(); ichir++) {
      bool minimal = false;
      int idx = get_monomer_restraints_index(comp_ids_for_filtering[ichir], imol, minimal);
      if (idx != -1) {
	 const coot::dictionary_residue_restraints_t &restraints =
	    dict_res_restraints[idx].second;
	 std::vector<coot::dict_chiral_restraint_t> new_chirals =
	    filter_chiral_centres(restraints);
	 dict_res_restraints[idx].second.chiral_restraint = new_chirals;
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

// Because they were all at origin, for example.
// 
void
coot::protein_geometry::delete_atom_positions(const std::string &comp_id,
					      int imol_enc,
					      int pos_type) {
   
   // for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
   // if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {

   int idr = get_monomer_restraints_index(comp_id, imol_enc, false);
   if (idr >= 0) {
      for (unsigned int iat=0; iat<dict_res_restraints[idr].second.atom_info.size(); iat++) { 
	 if (pos_type == dict_atom::IDEAL_MODEL_POS)
	    dict_res_restraints[idr].second.atom_info[iat].pdbx_model_Cartn_ideal.first = false;
	 if (pos_type == dict_atom::REAL_MODEL_POS)
	    dict_res_restraints[idr].second.atom_info[iat].model_Cartn.first = false;
      }
   }
}



void 
coot::protein_geometry::add_restraint(std::string comp_id,
				      int imol_enc,
				      const dict_torsion_restraint_t &restr) { 

   // if comp is in the list, simply push back restr, 
   // 
   // if not, then push back a dict_bond_restraint_t for it, passing
   // the comp_id. 

   bool ifound = false;

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) { 
      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
	 if (dict_res_restraints[i].first == imol_enc) {
	    ifound = true;
	    dict_res_restraints[i].second.torsion_restraint.push_back(restr); 
	    break;
	 }
      }
   }

   // it was not there
   if (! ifound) {
      dictionary_residue_restraints_t drr(comp_id, read_number);
      drr.torsion_restraint.push_back(restr);
      std::pair<int, dictionary_residue_restraints_t> p(imol_enc, drr);
      dict_res_restraints.push_back(p);
      // add the angle to the newly created dictionary_residue_restraints_t

      // I don't think that we need to do this - they can be pre-added.
      // use back? or reverse iterator? FIXME
      // dict_res_restraints[dict_res_restraints.size()-1].second.torsion_restraint.push_back(restr);
   }
}

void 
coot::protein_geometry::add_restraint(std::string comp_id,
				      int imol_enc,
				      const dict_chiral_restraint_t &restr) { 

   // if comp is in the list, simply push back restr, 
   // 
   // if not, then push back a dict_bond_restraint_t for it, passing
   // the comp_id. 

   bool ifound = false;

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
	 if (dict_res_restraints[i].first == imol_enc) {
	    ifound = true;
	    dict_res_restraints[i].second.chiral_restraint.push_back(restr); 
	    break;
	 }
      }
   }


   // it was not there
   if (! ifound) {
      std::cout << "---------------------------- oops missing in add_restraint() chiral " << std::endl;
      dictionary_residue_restraints_t drr(comp_id, read_number);
      drr.chiral_restraint.push_back(restr);
      std::pair<int, dictionary_residue_restraints_t> p(imol_enc, drr);
      dict_res_restraints.push_back(p);
      // add the angle to the newly created dictionary_residue_restraints_t
      // dict_res_restraints[dict_res_restraints.size()-1].
   }
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

std::ostream&
coot::operator<<(std::ostream &s, const dict_atom &at) {

   s << "dict_atom: "
     << "atom_id :" << at.atom_id << ":  "
     << "atom-id-4c :" << at.atom_id_4c << ":  "
     << "type-symbol :" << at.type_symbol << ":  "
     << "pdbx_stereo_config: " << at.pdbx_stereo_config.first
     << " \"" << at.pdbx_stereo_config.second << "\" ";
   if (at.formal_charge.first)
      s << "formal-charge " << at.formal_charge.second << " ";
   else
      s << "no-formal-charge ";
   if (at.partial_charge.first)
      s << "partial-charge " << at.partial_charge.second << " ";
   else
      s << "no-partial-charge ";
   s << "\n      model-pos-flag: " << at.model_Cartn.first << " ";
   if (at.model_Cartn.first)
      s << at.model_Cartn.second.format() << " ";
   else
      s << "no model_Cartn ";
   s << "ideal-pos-flag: " << at.pdbx_model_Cartn_ideal.first << " ";
   if (at.pdbx_model_Cartn_ideal.first)
      s << at.pdbx_model_Cartn_ideal.second.format();
   else
      s << "no model_Cartn_ideal";
   return s;
}


std::ostream&
coot::operator<<(std::ostream &s, const dict_bond_restraint_t &rest) {

   s << "[bond-restraint: " 
     << rest.atom_id_1_4c() << " "
     << rest.atom_id_2_4c() << " "
     << rest.type() << " " << std::setw(7) << rest.value_dist() << " " << rest.value_esd()
     <<  "]";
   return s;
}

std::ostream& coot::operator<<(std::ostream &s, const coot::dict_chem_comp_t &rest) {

   s << "[dict_chem_comp comp_id: \"" << rest.comp_id << "\" 3-letter-code: \""
     << rest.three_letter_code << "\" name: \"" << rest.name << "\" group: \"" << rest.group
     << "\" descr-level: \"" << rest.description_level << "\" "  << rest.number_atoms_all
     << " " << rest.number_atoms_nh << "]";
   return s;
}

std::ostream&
coot::operator<<(std::ostream &s, const dict_chiral_restraint_t &rest) {

   s << "[chiral: " << rest.Chiral_Id() << " "
     << rest.atom_id_c_4c() << " "
     << rest.atom_id_1_4c() << " "
     << rest.atom_id_2_4c() << " "
     << rest.atom_id_3_4c() << " "
     << rest.volume_sign << "]";
   return s;
}


std::ostream& coot::operator<<(std::ostream&s, coot::dict_plane_restraint_t rest) {

   s << "[plane-restraint: " << rest.plane_id << " " << " {"
     << rest.n_atoms() << " atoms} ";
   for (int iatom=0; iatom<rest.n_atoms(); iatom++) {
      s << ":" << rest[iatom].first << " " << rest[iatom].second << ": ";
   }
   s << "]";
   return s;
}

std::ostream&
coot::operator<<(std::ostream &s, const dict_angle_restraint_t &rest) {

   s << "[angle-restraint: " 
     << rest.atom_id_1_4c() << " "
     << rest.atom_id_2_4c() << " "
     << rest.atom_id_3_4c() << " "
     << rest.angle_ << " " << rest.angle_esd_
     <<  "]";
   return s;
} 

std::ostream&
coot::operator<<(std::ostream &s, const coot::dict_torsion_restraint_t &rest) {
   s << "[torsion-restraint: " << rest.id() << " "
     << "\"" << rest.atom_id_1_4c() << "\"" << " "
     << "\"" << rest.atom_id_2_4c() << "\"" << " "
     << "\"" << rest.atom_id_3_4c() << "\"" << " "
     << "\"" << rest.atom_id_4_4c() << "\"" << " "
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




std::string
coot::protein_geometry::get_padded_name(const std::string &atom_id,
					const int &comp_id_index) const {
   std::string s;
   if (comp_id_index < 0) {
      std::cout << "ERROR:: disaster in get_padded_name for comp_id_index "
		<< comp_id_index << " and atom name \"" << atom_id << "\"" << std::endl;
      return s;
   } else {
      for (unsigned int iat=0; iat<dict_res_restraints[comp_id_index].second.atom_info.size(); iat++) {
	 if (dict_res_restraints[comp_id_index].second.atom_info[iat].atom_id == atom_id) {
	    s = dict_res_restraints[comp_id_index].second.atom_info[iat].atom_id_4c;
	    break;
	 }
      }
   }
   return s;
}




void
coot::protein_geometry::info() const {

   std::cout << "::::: MONOMER GEOMETRY:" << std::endl;
   for (unsigned int idr=0; idr<size(); idr++) {
      // ejd says that "restraints" has an "n" in it.  Fixed.
      std::cout << dict_res_restraints[idr].second.residue_info.comp_id << std::endl;
      std::cout << "   " << dict_res_restraints[idr].second.bond_restraint.size()
                << " bond restraints " << std::endl;
      std::cout << "   " << dict_res_restraints[idr].second.angle_restraint.size()
                << " angle restraints " << std::endl;
      std::cout << "   " << dict_res_restraints[idr].second.torsion_restraint.size()
                << " torsion restraints " << std::endl;
      std::cout << "   " << dict_res_restraints[idr].second.plane_restraint.size()
                << " plane restraints " << std::endl;
//       for (int i=0; i<dict_res_restraints[idr].plane_restraint.size(); i++) {
//          for (int j=0; j<dict_res_restraints[idr].plane_restraint[i].atom_ids.size(); j++) {
//             std::cout << "      " << dict_res_restraints[idr].plane_restraint[i].plane_id
//                       << " " << dict_res_restraints[idr].plane_restraint[i].atom_ids[j]
//                       << std::endl;
//          }
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

void
coot::protein_geometry::set_verbose(bool verbose_flag) {
   verbose_mode = verbose_flag;
}


// return empty file name on failure.
std::string
coot::protein_geometry::comp_id_to_file_name(const std::string &comp_id) const {

   std::string file_name;

   if (comp_id.length() > 0) {
      char *cmld = getenv("COOT_MONOMER_LIB_DIR");
      std::string d;
      if (! cmld) {
         d = PKGDATADIR;
         d = util::append_dir_dir(d, "lib");
         d = util::append_dir_dir(d, "data");
         d = util::append_dir_dir(d, "monomers");
      } else {
         d = cmld;
      }

      if (! d.empty()) {
         std::string c0(1,comp_id[0]);
         std::string lc0 = util::downcase(c0);
         d = util::append_dir_dir(d, lc0);
         std::string lfn = comp_id + ".cif";
         file_name = util::append_dir_file(d, lfn);
      }
   }
   return file_name;
}

int
coot::protein_geometry::check_and_try_dynamic_add(const std::string &resname, int imol_enc, int read_number) {

   bool try_autoload = true;
   return have_dictionary_for_residue_type(resname, imol_enc, read_number, try_autoload);
}

// Return 0 on failure to do a dynamic add (actually, the number of
// atoms read).
//
// try_dynamic_add() will add with an imol of IMOL_ENC_ANY
//
int
coot::protein_geometry::try_dynamic_add(const std::string &resname, int read_number) {

   bool debug = false;
   // if (verbose_mode) debug = true;

   int success = 0;  // fail initially and is set to the number of
                     // atoms read from the mmcif dictionary file in
                     // init_refmac_mon_lib().

   // If this is INH, DRG etc, don't try to auto-add
   //
   if (is_non_auto_load_ligand(resname)) {
      std::cout << "INFO:: comp-id: " << resname
                << " is marked for non-autoloading - stopping dynamic_add() now "
                << std::endl;
      return success;
   }

   // So what is happening here?
   //
   // It is a heirachy of setting
   //
   // The highest priority is COOT_REFMAC_LIB_DIR, if that is set we use it.
   // If COOT_REFMAC_LIB_DIR is then try CLIB (a CCP4 setting).
   // If that is not set, then we fall back to the default directory:
   // $prefix/share/coot onto which we tag a "lib" dir.
   //

   char *s         = getenv("COOT_REFMAC_LIB_DIR");
   char *cmld      = getenv("COOT_MONOMER_LIB_DIR");
   char *clibd     = getenv("CLIBD");     // a CCP4 setting
   char *clibd_mon = getenv("CLIBD_MON"); // a CCP4 setting

   if (s) {

      if (debug)
         std::cout << "PATH p-A" << std::endl;

   } else {

      if (debug)
         std::cout << "PATH p-B" << std::endl;

      s = clibd;

      if (! s) {
         if (debug)
            std::cout << "PATH p-C" << std::endl;
         if (debug)
            std::cout << "DEBUG:: try_dynamic_add() using package_data_dir(): " << package_data_dir()
                      << std::endl;
         std::string dir_1 = package_data_dir();
         std::string dir_2 = util::append_dir_dir(dir_1, "lib");
         std::string dir_3 = util::append_dir_dir(dir_2, "data");
         s = new char[dir_3.length() + 1];
         strcpy(s, dir_3.c_str());
      } else {
         if (verbose_mode)
            std::cout << "INFO:: using standard CCP4 Refmac dictionary"
                      << " to search for \"" << resname << "\"" << std::endl;
      }
   }

   if (!s) {

      std::string pref_dir = coot::prefix_dir(); // checks COOT_PREFIX env var
      std::string ss = pref_dir + "/share/coot/lib/data";
      if (ss.length() < 2048)
         strcpy(cmld, ss.c_str());
   }

   if (debug) {
      std::cout << ":::::::::::::::: debug in try_dynamic_add() here A with s " << s << std::endl;
      if (cmld)
         std::cout << ":::::::::::::::: debug in try_dynamic_add() here A with cmld " << cmld << std::endl;
      else
         std::cout << ":::::::::::::::: debug in try_dynamic_add() here A with null cmld " << std::endl;
   }

   // 20241001-PE hostage to fortune...
   // when will I be back here?

   {
      std::string filename(s);
      std::string beta_anomer_name;
      std::string alpha_anomer_name;
      std::string alt_beta_anomer_name;
      std::string alt_alpha_anomer_name;

      if (cmld) {
         filename = cmld;
         filename += "/";
      } else {
         filename += "/monomers/";
      }

      filename = coot::util::intelligent_debackslash(filename);

      if (debug)
         std::cout << "debug:: try_dynamic_add(): here C filename " << filename << std::endl;

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
            if (coot::is_regular_file(filename)) {

               coot::read_refmac_mon_lib_info_t rmit = init_refmac_mon_lib(filename, read_number);
               success = rmit.success;
            } else {

               // error/warning
               if (! coot::is_dir_or_link(filename)) {
                  std::cout << "WARNING: " << filename << ": no such file (or directory)\n";
               } else {
                  std::cout << "ERROR: dictionary " << filename << " is not a regular file" << std::endl;
               }
            }
         } else {

            // Regular file doesn't exist,

            // Try the upcased filename
            //
            istat = stat(upcased_resname_filename.c_str(), &buf);
            if (istat == 0) {
               read_refmac_mon_lib_info_t rmit = init_refmac_mon_lib(upcased_resname_filename, read_number);
               success = rmit.success;
            } else {

               // try the beta anomer version
               istat = stat(beta_anomer_name.c_str(), &buf);
               if (istat == 0) {
                  read_refmac_mon_lib_info_t rmit = init_refmac_mon_lib(beta_anomer_name, read_number);
                  success = rmit.success;
               } else {
                  // try the upcased file name e.g. xxx/NAG-B-D.cif
                  istat = stat(alt_beta_anomer_name.c_str(), &buf);
                  if (istat == 0) {
                     read_refmac_mon_lib_info_t rmit = init_refmac_mon_lib(alt_beta_anomer_name, read_number);
                     success = rmit.success;
                  } else {
                     // alpha?
                     istat = stat(alpha_anomer_name.c_str(), &buf);
                     if (istat == 0) {
                        read_refmac_mon_lib_info_t rmit = init_refmac_mon_lib(alpha_anomer_name, read_number);
                        success = rmit.success;
                     } else {
                        istat = stat(alt_alpha_anomer_name.c_str(), &buf);
                        if (istat == 0) {
                           read_refmac_mon_lib_info_t rmit = init_refmac_mon_lib(alt_alpha_anomer_name, read_number);
                           success = rmit.success;
                        }
                     }
                  }
               }
            }
         }
      }
   }

   if (!success) {
      // try the XDG Base Directory Protocol cache
      xdg_t xdg;
      bool debug = false;
      std::filesystem::path ch = xdg.get_cache_home();
      if (std::filesystem::exists(ch)) {
         std::filesystem::path monomers_path = ch / "monomers";
         if (std::filesystem::exists(monomers_path)) {
            const char rs = resname[0];
            const char v = tolower(rs); // get the sub directory name
            std::string letter(1, v);
            std::filesystem::path sub_dir = monomers_path / letter;
            if (std::filesystem::exists(sub_dir)) {
               std::string cif_file_name = resname + ".cif";
               std::filesystem::path cif_file_path = sub_dir / cif_file_name;
               if (std::filesystem::exists(cif_file_path)) {
                  // read it
                  read_refmac_mon_lib_info_t rmit = init_refmac_mon_lib(cif_file_path.string(), read_number);
                  success = rmit.success;
               } else {
                  // we will need to download it then
                  // and put it in the above directory)
                  if (debug)
                     std::cout << "DEBUG:: try_dynamic_add(): " << cif_file_path.string() << " does not exist" << std::endl;
               }
            } else {
               if (debug)
                  std::cout << "DEBUG:: try_dynamic_add(): " << sub_dir.string() << " does not exist" << std::endl;
            }
         } else {
            if (debug)
               std::cout << "DEBUG:: try_dynamic_add(): " << monomers_path.string() << " does not exist" << std::endl;
         }
      } else {
         if (debug)
            std::cout << "DEBUG:: try_dynamic_add(): " << ch.string() << " does not exist" << std::endl;
      }
   }

   // now, did we get something with minimal description? If so,
   // delete it, it was a fail.
   //
   std::pair<bool, dictionary_residue_restraints_t> p = get_monomer_restraints(resname, IMOL_ENC_ANY);

   if (p.first) {
      if (resname == "3GP") {
         for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
            if (dict_res_restraints[i].second.residue_info.comp_id == resname) {
               auto &dict = dict_res_restraints[i].second;
               dict.move_3GP_atoms();
               break;
            }
         }
      }
   }

   if (! p.first) {
      success = 0;
   } else {
      // elide minimal description restraints.
      if (p.second.residue_info.description_level == "M") {
         success = 0;
         delete_mon_lib(resname, IMOL_ENC_ANY);
      }
   }
   return success;
}


bool
coot::protein_geometry::is_non_auto_load_ligand(const std::string &resname) const {

   bool r = false;
   std::vector<std::string>::const_iterator it;
   for (it=non_auto_load_residue_names.begin(); it!=non_auto_load_residue_names.end(); ++it) {
      if (*it == resname) {
         r = true;
         break;
      }
   }
   return r;
}

void
coot::protein_geometry::add_non_auto_load_residue_name(const std::string &res_name) {

   bool found = false;
   std::vector<std::string>::const_iterator it;
   for (it=non_auto_load_residue_names.begin(); it!=non_auto_load_residue_names.end(); it++) {
      if (*it == res_name) {
         found = true;
         break;
      }
      if (found)
         break;
   }
   if (! found)
      non_auto_load_residue_names.push_back(res_name);
}

void
coot::protein_geometry::remove_non_auto_load_residue_name(const std::string &res_name) {

   std::vector<std::string>::iterator it;
   for (it=non_auto_load_residue_names.begin(); it!=non_auto_load_residue_names.end(); it++) {
      if (*it == res_name) {
         non_auto_load_residue_names.erase(it);
         break;
      }
   }
}



void
coot::protein_geometry::fill_default_non_auto_load_residue_names() { // build-it default
   non_auto_load_residue_names.push_back("LIG");
   non_auto_load_residue_names.push_back("DRG");
   non_auto_load_residue_names.push_back("INH");
   non_auto_load_residue_names.push_back("01");
   non_auto_load_residue_names.push_back("02");
   non_auto_load_residue_names.push_back("03");
   non_auto_load_residue_names.push_back("04");
   non_auto_load_residue_names.push_back("05");
   non_auto_load_residue_names.push_back("06");
   non_auto_load_residue_names.push_back("07");
   non_auto_load_residue_names.push_back("08");
   non_auto_load_residue_names.push_back("09");
   for (unsigned int i=1; i<10; i++) {
      for (unsigned int j=0; j<10; j++) {
         std::string ss = std::to_string(i) + std::to_string(j);
         non_auto_load_residue_names.push_back(ss);
      }
   }
}


std::vector <coot::dict_torsion_restraint_t>
coot::protein_geometry::get_monomer_torsions_from_geometry(const std::string &monomer_type,
                                                           int imol_enc) {

   bool ifound = false;
   std::vector <coot::dict_torsion_restraint_t> rv;

   int idr = get_monomer_restraints_index(monomer_type, imol_enc, false);

   if (idr >= 0 ) {
      ifound = true;
      return dict_res_restraints[idr].second.torsion_restraint;
   }

   // check synonyms before 3-letter codes.  Maybe needs matches_imol() test
   //
   if (! ifound) {
      // OK, that failed to, perhaps there is a synonym?
      for (unsigned int i=0; i<residue_name_synonyms.size(); i++) {
	 if (residue_name_synonyms[i].comp_alternative_id == monomer_type) {
	    if (matches_imol(dict_res_restraints[i].first, imol_enc)) { // is this the right test?
	       int ndict = dict_res_restraints.size();
	       for (int i=0; i<ndict; i++) {
		  if (dict_res_restraints[i].second.residue_info.comp_id == residue_name_synonyms[i].comp_id) {
		     ifound = 1;
		     rv = dict_res_restraints[i].second.torsion_restraint;
		     break;
		  }
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
	    if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	       ifound = 1;
	       rv = dict_res_restraints[i].second.torsion_restraint;
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
							   int imol_enc,
							   bool find_hydrogen_torsions_flag) const {

   bool ifound = 0;
   int ii = -1; // unset, set when torsion found, ii is the index in
		// dict_res_restraints for the torsion restraints,
		// needed so that we can ask about hydrogens.
   std::vector <coot::dict_torsion_restraint_t> unfiltered_torsion_restraints;
   std::vector <coot::dict_torsion_restraint_t> filtered_torsion_restraints;

   int idr = get_monomer_restraints_index(monomer_type, imol_enc, false);

   if (idr >= 0) {
      ifound = 1;
      unfiltered_torsion_restraints = dict_res_restraints[idr].second.torsion_restraint;
      ii = idr; // used for filtering out hydrogens
   }

   if (! ifound) { // which it should be if we get here with that
		      // return in the loop
      // try dynamic add?
      std::cout << "WARNING: residue type " << monomer_type << " not found "
		<< "in restraints dictionary (in get_monomer_torsions_from_geometry(mon, hy)" << std::endl;
   } else {
      if (find_hydrogen_torsions_flag) {
	 filtered_torsion_restraints = unfiltered_torsion_restraints;
      } else {
	 // we don't want torsions that move Hydrogens
	 int nt = dict_res_restraints[ii].second.torsion_restraint.size();
	 for (int it=0; it<nt; it++) {
            // std::cout << "testing for hydrogen: this \"" << dict_res_restraints[ii].second << "\" and this \""
            //           << dict_res_restraints[ii].second << std::endl;
	    if (!dict_res_restraints[ii].second.is_hydrogen(dict_res_restraints[ii].second.torsion_restraint[it].atom_id_1())) {
	       if (!dict_res_restraints[ii].second.is_hydrogen(dict_res_restraints[ii].second.torsion_restraint[it].atom_id_4())) {
		  filtered_torsion_restraints.push_back(dict_res_restraints[ii].second.torsion_restraint[it]);
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
      bool match = false;
      for (unsigned int j=0; j<r.size(); j++) {
	 if (r[j].atom_id_2_4c() == a2)
	    if (r[j].atom_id_3_4c() == a3)
	       match = true;
	 if (r[j].atom_id_2_4c() == a3)
	    if (r[j].atom_id_3_4c() == a2)
	       match = true;
      }
      if (! match) {
	 const dict_torsion_restraint_t &tr = restraints_in[i];
	 r.push_back(tr);
      }
   }


   // why do we crash (on the mac) in the following sort?
   //
   // for (unsigned int i=0; i<restraints_in.size(); i++)
   //    std::cout << " filter_torsion_restraints(): " << restraints_in[i] << std::endl;
   // 20180227 still crashing.
   // Because it is running the sort even if r is of length 0? No, not that.
   // std::cout << " r.size() " << r.size() << std::endl;
   try {
      std::sort(r.begin(), r.end(), torsion_restraints_comparer);
   }
   catch (const std::bad_alloc &e) {
      std::cout << "ERROR caught bad alloc when sorting in filter_torsion_restraints()" << std::endl;
   }
   return r;
}

// static
bool
coot::protein_geometry::torsion_restraints_comparer(const coot::dict_torsion_restraint_t &a, const coot::dict_torsion_restraint_t &b) {

//    const std::string &a2 = a.atom_id_2_4c();
//    const std::string &a3 = a.atom_id_3_4c();
//    const std::string &b2 = b.atom_id_2_4c();
//    const std::string &b3 = b.atom_id_3_4c();

//    const std::string &a2 = a.atom_id_2();
//    const std::string &a3 = a.atom_id_3();
//    const std::string &b2 = b.atom_id_2();
//    const std::string &b3 = b.atom_id_3();

   // we leave these (debugging/test) copies in for now (was testing if reference was the problem
   // (it wasn't).
   //
   std::string a2 = a.atom_id_2();
   std::string a3 = a.atom_id_3();
   std::string b2 = b.atom_id_2();
   std::string b3 = b.atom_id_3();

   if (a2 < b2) {
      return false;
   } else {
      if (a2 > b2) {
 	 return true;
      } else {
	 // a2 and b2 are equal

	 if (a3 <= b3) { // 20180227-PE changed to <=, stops crash

	    // std::cout << "a3 " << a3 << std::endl;
	    // std::cout << "b3 " << b3 << std::endl;

	    return false;

	 }
      }
   }

   //       else
   // 	 if (a3 < b3)
   // 	    return false;

   return true;
}


std::vector <coot::dict_chiral_restraint_t>
coot::protein_geometry::get_monomer_chiral_volumes(const std::string &monomer_type,
						   int imol_enc) const { 

   bool ifound = 0;
   std::vector <coot::dict_chiral_restraint_t> rv;
   
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].first == imol_enc) {
	 if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	    ifound = true;
	    rv = dict_res_restraints[i].second.chiral_restraint;
	    break;
	 }
      }
   }
   if (! ifound) {
      for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
	 if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
	    if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	       ifound = true;
	       rv = dict_res_restraints[i].second.chiral_restraint;
	       break;
	    }
	 }
      }
   }
   // OK so the monomer_type did not match the comp_id.  Perhaps the
   // comp_id was not the same as the three letter code, so let's
   // check the monomer_type against the three_letter_codes.
   // 
   if (ifound == 0) {
      for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
	 if (dict_res_restraints[i].second.residue_info.three_letter_code == monomer_type) {
	    ifound = 1;
	    rv = dict_res_restraints[i].second.chiral_restraint;
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

std::pair<bool, coot::dict_atom>
coot::protein_geometry::get_monomer_atom_info(const std::string &monomer_name,
					      const std::string &atom_name,
					      int imol_enc) const {

   bool status = false;
   dict_atom da;

   std::pair<bool, dictionary_residue_restraints_t> dr = get_monomer_restraints(monomer_name, imol_enc);
   if (dr.first) {
      const std::vector<dict_atom> &atom_info = dr.second.atom_info;
      for (std::size_t i=0; i<atom_info.size(); i++) {
	 dict_atom a = atom_info[i];
	 // std::cout << "comparing atom names \"" << atom_name << "\" and \""
	 // << a.atom_id_4c << "\"" << std::endl;
	 if (atom_name == a.atom_id_4c) {
	    da = a;
	    status = true;
	 }
      }
   }

   return std::pair<bool, dict_atom> (status, da);
}

std::vector<std::pair<int, std::string> >
coot::protein_geometry::get_monomer_names() const {

   std::vector<std::pair<int, std::string> > v;
   unsigned int nrest = dict_res_restraints.size();
   for (unsigned int i=0; i<nrest; i++) {
      const auto &r = dict_res_restraints[i];
      std::pair<int, std::string> p(r.first, r.second.residue_info.comp_id);
      v.push_back(p);
   }
   return v;
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
   s.push_back("c/CIT.cif");
   s.push_back("e/EDO.cif");
   // s.push_back("e/ETH.cif");

   // s.push_back("a/AR.cif"); // argon
   // s.push_back("c/CR.cif"); // Chromium
   // s.push_back("g/GR.cif"); old

   // s.push_back("a/AD.cif");
   // s.push_back("c/CD.cif");
   // s.push_back("g/GD.cif");
   // s.push_back("t/TD.cif");
   // s.push_back("u/UR.cif");

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

   // Metals
   s.push_back("n/NA.cif");
   s.push_back("k/K.cif");
   s.push_back("c/CA.cif");
   s.push_back("m/MG.cif");
   s.push_back("n/NI.cif");
   s.push_back("z/ZN.cif");

   return s;
}


// optional arg: bool try_autoload_if_needed=true
bool
coot::protein_geometry::have_dictionary_for_residue_type(const std::string &monomer_type,
							 int imol_enc,
							 int read_number_in,
							 bool try_autoload_if_needed) {

   bool debug = false;

   bool ifound = false;
   int ndict = dict_res_restraints.size();
   std::string path = "--start--";
   read_number = read_number_in;

   // ---------------- FIXME ----------------------------------------------------
   // ---------------- FIXME ----------------------------------------------------
   // ---------------- FIXME ----------------------------------------------------
   int idr = get_monomer_restraints_index(monomer_type, imol_enc, true);
   if (idr >= 0) {
      ifound = true;
   }

   if (debug)
      std::cout << "INFO:: pg::have_dictionary_for_residue_type() idr here is " << idr << std::endl;

   // check synonyms before checking three-letter-codes

   if (! ifound) {
      // OK, that failed to, perhaps there is a synonym?
      for (unsigned int i=0; i<residue_name_synonyms.size(); i++) {
	 if (residue_name_synonyms[i].comp_alternative_id == monomer_type) {
	    for (int j=0; j<ndict; j++) {
	       if (dict_res_restraints[j].second.residue_info.comp_id == residue_name_synonyms[i].comp_id) {
		  ifound = true;
		  break;
	       }
	    }
	 }
	 if (ifound) {
            path = "path-1";
	    break;
         }
      }
   }

   // OK so the monomer_type did not match the comp_id.  Perhaps the
   // comp_id was not the same as the three letter code, so let's
   // check the monomer_type against the three_letter_codes.
   //
   if (ifound == 0) {
      for (int i=0; i<ndict; i++) {
	 if (dict_res_restraints[i].second.residue_info.three_letter_code == monomer_type) {
	    if (! dict_res_restraints[i].second.is_bond_order_data_only()) {
	       ifound = 1;
               path = "path-2";
	       break;
	    }
	 }
      }
   }

   if (! ifound) {
      if (try_autoload_if_needed) {
	 ifound = try_dynamic_add(monomer_type, read_number);
         path = "path-3";
	 // std::cout << "DEBUG:: here in have_dictionary_for_residue_type() try_dynamic_add returned "
         // << ifound << std::endl;
      }
   }

   if (debug)
      std::cout << "INFO:: pg::have_dictionary_for_residue_type() " << monomer_type
                << " " << imol_enc << " path " << path << " returns " << ifound << std::endl;

   return ifound;
}

// this is const because there is no dynamic add
// This is a test for "shall we add a new dictionary?"
bool
coot::protein_geometry::have_dictionary_for_residue_type_no_dynamic_add(const std::string &monomer_type, int imol_no) const {

   bool ifound = false;

   int ndict = dict_res_restraints.size();
   for (int i=0; i<ndict; i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	 if (matches_imol(dict_res_restraints[i].first, imol_no)) {
	    if (! dict_res_restraints[i].second.is_bond_order_data_only()) {
	       ifound = 1;
	       break;
	    }
	 }
      }
   }
   return ifound;
}


bool
coot::protein_geometry::have_dictionary_for_residue_types(const std::vector<std::string> &residue_types,
							  int imol_enc,
							  int read_number) {

   bool have_all = 1;
   for (unsigned int i=0; i<residue_types.size(); i++) {
      int ifound = have_dictionary_for_residue_type(residue_types[i], imol_enc, read_number);
      if (ifound == 0) {
	 have_all = 0;
      }
      read_number++;
   }
   return have_all;
}

// Return false if there are no bond restraints
bool
coot::protein_geometry::have_restraints_dictionary_for_residue_types(const std::vector<std::string> &residue_types,
                                                                     int imol_enc,
                                                                     int read_number) {

   if (false) {
      std::cout << "debug:: in have_restraints_dictionary_for_residue_types() --- start --- " << std::endl;
      for (const auto &r : residue_types)
	 std::cout << "debug:: in have_restraints_dictionary_for_residue_types() type :" << r << ":" << std::endl;
   }

   bool have_all = true;
   for (unsigned int i=0; i<residue_types.size(); i++) {
      if (! have_all) continue;
      const std::string &rt = residue_types[i];
      int idx = get_monomer_restraints_index(rt, imol_enc, false);  // const function
      if (idx != -1) {
         const coot::dictionary_residue_restraints_t &restraints = dict_res_restraints[idx].second;
	 // this test does not make sense for MG, or other single atoms
         // if (restraints.bond_restraint.empty()) {
	 bool has_bonds = false;
	 if (restraints.atom_info.size() > 1) has_bonds = true;
         if (has_bonds && restraints.bond_restraint.empty()) {
            have_all = false;
            break;
         } else {
            for (auto it=restraints.bond_restraint.begin(); it!=restraints.bond_restraint.end(); ++it) {
               try {
                  // this will throw an exception if there are no bond *restraints*
                  float v = it->value_dist();
               }
               catch (const std::runtime_error &e) {
                  have_all = false;
               }
            }
         }
      } else {
         have_all = false;
         break;
      }
      read_number++;
   }
   return have_all;
}


// this is const because there is no dynamic add
//
// maybe imol should be passed.
bool
coot::protein_geometry::have_at_least_minimal_dictionary_for_residue_type(const std::string &monomer_type,
									  int imol) const {

   //std::cout << "debug:: in have_at_least_minimal_dictionary_for_residue_type() " << monomer_type
   // << " " << imol << std::endl;

   bool ifound = false;
   std::size_t ndict = dict_res_restraints.size();

   for (std::size_t i=0; i<ndict; i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	 if (matches_imol(dict_res_restraints[i].first, IMOL_ENC_ANY)) {
	    ifound = true;
	    break;
	 }
	 if (matches_imol(dict_res_restraints[i].first, imol)) {
	    ifound = true;
	    break;
	 }
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
//  function name: do_the_atom_names_match_the_dictionary?() - it is a question
//  not an imperative.
// 
std::pair<bool, std::vector<std::string> >
coot::protein_geometry::atoms_match_dictionary(mmdb::Residue *residue_p,
					       bool check_hydrogens_too_flag,
					       bool apply_bond_distance_check,
					       const coot::dictionary_residue_restraints_t &restraints) const {

   std::vector<std::string> atom_name_vec;
   bool status = 1; // nothing fails to match (so far).

   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

   bool debug = false;
   if (debug) {
      std::cout << "=== atoms_match_dictionary() with these residue atom names ======= " << std::endl;
      for (int i=0; i<n_residue_atoms; i++) {
	 std::cout << i << "  :" << residue_atoms[i]->name << ": and ele "
		   << residue_atoms[i]->element << std::endl;
      } 
      std::cout << "=== atoms_match_dictionary() with these residue atom names ======= " << std::endl;
      for (unsigned int irat=0; irat<restraints.atom_info.size(); irat++) {
	 std::cout << irat << "  :" << restraints.atom_info[irat].atom_id_4c
		   << ":" << std::endl;
      } 
   }

   for (int i=0; i<n_residue_atoms; i++) {

      if (! residue_atoms[i]->isTer()) { 
	 std::string residue_atom_name(residue_atoms[i]->name);
	 std::string ele(residue_atoms[i]->element);

	 bool found = 0;
	 // PDBv3 FIXME
	 if (ele == " H" || ele == " D")
	    if (check_hydrogens_too_flag == false)
	       found = true;

	 if (! found) {
	    for (unsigned int irestraint_atom_name=0; irestraint_atom_name<restraints.atom_info.size(); irestraint_atom_name++) {
	       if (restraints.atom_info[irestraint_atom_name].atom_id_4c == residue_atom_name) {
		  found = 1;
        atom_name_vec.push_back(residue_atom_name);

		  break;
	       }
	    }
	 }
	 if (! found) {
	    if (residue_atom_name != " OXT") {
         if (debug) std::cout << "here d" << std::endl;
               if (std::find(atom_name_vec.begin(), atom_name_vec.end(), residue_atom_name) == atom_name_vec.end()) {
                  if (debug) std::cout << "here e" << std::endl;
                  atom_name_vec.push_back(residue_atom_name);
               }
               status = false;
	    }
	 }
      }
   }

   // We can finally fail to match because we have some very long
   // bonds (but no atom name mismatches, of course)
   // 
   if (status && apply_bond_distance_check)
      status = atoms_match_dictionary_bond_distance_check(residue_p, check_hydrogens_too_flag, restraints);

   return std::pair<bool, std::vector<std::string> > (status, atom_name_vec);
}


std::pair<bool, std::vector<std::string> >
coot::protein_geometry::atoms_match_dictionary(int imol,
					       mmdb::Residue *residue_p,
					       bool check_hydrogens_too_flag,
					       bool apply_bond_distance_check) const {

   std::string res_name(residue_p->GetResName());
   std::pair<bool, dictionary_residue_restraints_t> restraints = get_monomer_restraints(res_name, imol);

   if (restraints.first) {
      return atoms_match_dictionary(residue_p, check_hydrogens_too_flag,
				    apply_bond_distance_check, restraints.second);
   } else { 
      std::vector<std::string> atom_name_vec;
      return std::pair<bool, std::vector<std::string> > (false, atom_name_vec);
   }
}


// return a pair, overall status, and vector of pairs of residue names and
// atom names that dont't match.
//
std::pair<bool, std::vector<std::pair<mmdb::Residue *, std::vector<std::string> > > >
coot::protein_geometry::atoms_match_dictionary(int imol,
					       const std::vector<mmdb::Residue *> &residues,
					       bool check_hydrogens_too_flag,
					       bool apply_bond_distance_check) const {

   bool status = true;
   std::vector<std::pair<mmdb::Residue *, std::vector<std::string> > > p;
   
   for (unsigned int ires=0; ires<residues.size(); ires++) { 
      std::string res_name(residues[ires]->GetResName());
      std::pair<bool, coot::dictionary_residue_restraints_t> restraints =
	 get_monomer_restraints(res_name, imol);
      if (restraints.first) {
	 std::pair<bool, std::vector<std::string> > r =
	    atoms_match_dictionary(residues[ires],
				   check_hydrogens_too_flag,
				   apply_bond_distance_check,
				   restraints.second);
	 if (r.first == false) {
	    std::pair<mmdb::Residue *, std::vector<std::string> > p_bad(residues[ires], r.second);
	    p.push_back(p_bad);
	    status = 0;
	 }
      } else {
	 std::cout << "ERROR:: atoms_match_dictionary() --- no restraints" << std::endl;
      }
   }
   return std::pair<bool, std::vector<std::pair<mmdb::Residue *, std::vector<std::string> > > > (status, p);
}

bool
coot::protein_geometry::atoms_match_dictionary_bond_distance_check(mmdb::Residue *residue_p,
								   bool check_hydrogens_too_flag,
								   const coot::dictionary_residue_restraints_t &restraints) const { 

   bool status = true; // good
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   if (n_residue_atoms > 2) { 
      for (unsigned int ibond=0; ibond<restraints.bond_restraint.size(); ibond++) {
	 for (int iat=0; iat<(n_residue_atoms-1); iat++) {
	    const mmdb::Atom *at_1 = residue_atoms[iat];
	    std::string atom_name_1(at_1->name);
	    if (restraints.bond_restraint[ibond].atom_id_1_4c() == atom_name_1) { 
	       for (int jat=iat+1; jat<n_residue_atoms; jat++) {
		  const mmdb::Atom *at_2 = residue_atoms[jat];
		  std::string atom_name_2(at_2->name);
		  if (restraints.bond_restraint[ibond].atom_id_2_4c() == atom_name_2) {
		     std::string alt_conf_1(at_1->altLoc);
		     std::string alt_conf_2(at_2->altLoc);
		     if (alt_conf_1 == alt_conf_2) {
			clipper::Coord_orth pt1(at_1->x, at_1->y, at_1->z);
			clipper::Coord_orth pt2(at_2->x, at_2->y, at_2->z);
			double d = (pt1-pt2).lengthsq();
			if (d > 10) {
			   status = false;
			   break;
			}
		     }
		  }
		  if (status == false)
		     break;
	       }
	    }
	    if (status == false)
	       break;
	 }
	 if (status == false)
	    break;
      }
   }
   return status;
}

// return a pair, the first is status (1 if the name was found, 0 if not)
// 
std::pair<bool, std::string>
coot::protein_geometry::get_monomer_name(const std::string &comp_id, int imol_enc) const {

   std::pair<bool, std::string> r(false,"");

   bool allow_minimal_flag = true;
   std::pair<bool, dictionary_residue_restraints_t> rest =
      get_monomer_restraints_internal(comp_id, imol_enc, allow_minimal_flag);

   if (rest.first) { 
      r.first = true;
      std::string s = rest.second.residue_info.name;
      r.second = coot::util::remove_trailing_whitespace(s);
   }
   return r;
} 



bool
coot::protein_geometry::copy_monomer_restraints(const std::string &monomer_type, int imol_enc_current, int imol_enc_new) {

   bool status = false;
   std::pair<bool, coot::dictionary_residue_restraints_t> r =
      get_monomer_restraints_internal(monomer_type, imol_enc_current, false);

   if (r.first) {
      status = true;
      dict_res_restraints.push_back(std::make_pair(imol_enc_new, r.second));
   }

   return status;

}


// Try comparing vs the comp_id first, if that fails compare the
// three_letter_code to the monomer_type.
//
// In future, try to come here only with the monomer_type adjusted to
// the comp_id, for example, monomer_type should be "NAG-b-D", not
// "NAG".
//
std::pair<bool, coot::dictionary_residue_restraints_t>
coot::protein_geometry::get_monomer_restraints(const std::string &monomer_type,
					       int imol_enc) const {

   return get_monomer_restraints_internal(monomer_type, imol_enc, 0);

}

std::pair<bool, coot::dictionary_residue_restraints_t>
coot::protein_geometry::get_monomer_restraints_at_least_minimal(const std::string &monomer_type,
								int imol_enc) const {

   return get_monomer_restraints_internal(monomer_type, imol_enc, 1);
}

std::pair<bool, coot::dictionary_residue_restraints_t>
coot::protein_geometry::get_monomer_restraints_internal(const std::string &monomer_type,
							int imol_enc,
							bool allow_minimal_flag) const {

   if (false)
      std::cout << "debug:: get_monomer_restraints_internal() called with type: "
                << monomer_type << " imol_enc: " << imol_enc
                << " allow-minimal: " << allow_minimal_flag << std::endl;

   // 20161028
   // Compiling with SRS causes a crash when we access dict_res_restraints.
   // Needs more testing.

   coot::dictionary_residue_restraints_t t(std::string("(null)"), 0);
   std::pair<bool, dictionary_residue_restraints_t> r(0,t);
   unsigned int nrest = dict_res_restraints.size();

   // This is how it used to be - starting from the beginning.
   // 20161004
   // Now we want to start from the end.
   //
//    for (unsigned int i=0; i<nrest; i++) {
//       if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
// 	 if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
// 	    r.second = dict_res_restraints[i].second;
// 	    r.first = true;
// 	    break;
// 	 }
//       }
//    }

   // Compiler error.  Not sure why
   // std::vector<std::pair<int, dictionary_residue_restraints_t> >::reverse_iterator rit;
   // for (rit=dict_res_restraints.rbegin(); rit!=dict_res_restraints.rend(); rit--) {
   // }

   // Try to match most recent exact match by model number
   //
   for (int i=(nrest-1); i>=0; i--) {
      if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
 	 if (dict_res_restraints[i].first == imol_enc) {
 	    r.second = dict_res_restraints[i].second;
 	    r.first = true;
 	    break;
 	 }
      }
   }

   // Try to match any dictionary (this is what will match in most cases)
   //
   if (!r.first) {
      for (int i=(nrest-1); i>=0; i--) {
	 // std::cout << "Matching :" << dict_res_restraints[i].second.residue_info.comp_id
	 // << ": :" << monomer_type << ":" << std::endl;
	 if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	    // std::cout << "   comparing " << dict_res_restraints[i].first << " " << imol_enc << std::endl;
	    if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
	       r.second = dict_res_restraints[i].second;
	       r.first = true;
	       break;
	    }
            // try "any" then...
	    if (matches_imol(dict_res_restraints[i].first, IMOL_ENC_ANY)) {
	       r.second = dict_res_restraints[i].second;
	       r.first = true;
	       break;
	    }
	 }
      }
   }

   if (!r.first) {
      // OK, that failed, perhaps there is a synonym?
      for (unsigned int i=0; i<residue_name_synonyms.size(); i++) {
	 if (residue_name_synonyms[i].comp_alternative_id == monomer_type) {
	    int ndict = dict_res_restraints.size();
	    for (int j=0; j<ndict; j++) {
	       if (dict_res_restraints[j].second.residue_info.comp_id == residue_name_synonyms[i].comp_id) {
		  r.first = true;
		  r.second = dict_res_restraints[j].second;
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
	 if (dict_res_restraints[i].second.residue_info.three_letter_code == monomer_type) {
	    if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
	       if (allow_minimal_flag || (! dict_res_restraints[i].second.is_bond_order_data_only())) {
		  r.second = dict_res_restraints[i].second;
		  r.first = true;
		  break;
	       }
	    }
	 }
      }
   }
   return r;
}

// return -1 on monomer not found.
int
coot::protein_geometry::get_monomer_restraints_index(const std::string &monomer_type,
						     int imol_enc,
						     bool allow_minimal_flag) const {

   int r = -1;
   bool debug = false; // hello again

   unsigned int nrest = dict_res_restraints.size();
   for (unsigned int i=0; i<nrest; i++) {

      if (debug)
	 std::cout << "in get_monomer_restraints_index() comparing dict: \""
		   << dict_res_restraints[i].second.residue_info.comp_id << "\" vs mine: \"" << monomer_type
		   << "\" and dict: " << dict_res_restraints[i].first << " vs mine: " <<  imol_enc
		   << "     with allow_minimal_flag " << allow_minimal_flag << std::endl;

      if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	 if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
	    // if (dict_res_restraints[i].first == imol_enc) {
	    if ((allow_minimal_flag == 1) || (! dict_res_restraints[i].second.is_bond_order_data_only())) {
	       r = i;
	       break;
	    }
	 }
      }
   }

   if (r == -1) {
      for (unsigned int i=0; i<nrest; i++) {
	 if (dict_res_restraints[i].second.residue_info.comp_id  == monomer_type) {
	    if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
	       if ((allow_minimal_flag == 1) || (! dict_res_restraints[i].second.is_bond_order_data_only())) {
		  r = i;
		  break;
	       }
	    }
	 }
      }
   }

   if (r == -1) {
      // OK, that failed to, perhaps there is a synonym?
      for (unsigned int i=0; i<residue_name_synonyms.size(); i++) {
	 if (residue_name_synonyms[i].comp_alternative_id == monomer_type) {
	    if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
	       int ndict = dict_res_restraints.size();
	       for (int j=0; j<ndict; j++) {
		  if (dict_res_restraints[j].second.residue_info.comp_id == residue_name_synonyms[i].comp_id) {
		     r = j;
		     break;
		  }
	       }
	    }
	    if (r != -1)
	       break;
	 }
      }
   }

   if (r == -1) {
      for (unsigned int i=0; i<nrest; i++) {
	 if (dict_res_restraints[i].second.residue_info.three_letter_code == monomer_type) {
	    if (matches_imol(dict_res_restraints[i].first, imol_enc)) {
	       if ((allow_minimal_flag == 1) || (! dict_res_restraints[i].second.is_bond_order_data_only())) {
		  r = i;
		  break;
	       }
	    }
	 }
      }
   }

   // std::cout << "get_monomer_restraints_index() for " << monomer_type << " returns r " << r << std::endl;

   return r;
}

// 20250124-PE and the reverse
//
// return the second blank on lookup failure
std::pair<int, std::string>
coot::protein_geometry::get_monomer_name(int monomer_index) const {

   int imol = -1;
   std::string comp_id;
   int n_rest = dict_res_restraints.size();
   if (monomer_index >= 0) {
      if (monomer_index < n_rest) {
         imol = dict_res_restraints[monomer_index].first;
         comp_id = dict_res_restraints[monomer_index].second.residue_info.comp_id;
      }
   }
   return std::make_pair(imol, comp_id);
}



std::string
coot::protein_geometry::get_type_energy(const std::string &atom_name,
					const std::string &residue_name,
					int imol) const {
   // return "" if not found, else return the energy type found in ener_lib.cif
   //
   std::string r;
   int indx = get_monomer_restraints_index(residue_name, imol, 1); // in get_type_energy()
   if (indx != -1) {
      const coot::dictionary_residue_restraints_t &restraints = dict_res_restraints[indx].second;
      r = restraints.type_energy(atom_name);
   }
   return r;
}

// return -1.1 on failure to look up.
//
double
coot::protein_geometry::get_vdw_radius(const std::string &atom_name,
				       const std::string &residue_name,
				       int imol,
				       bool use_vdwH_flag) const {
   double r = -1.1;
   int indx = get_monomer_restraints_index(residue_name, imol, 1); // in get_vdw_radius
   if (indx != -1) {
      const coot::dictionary_residue_restraints_t &restraints = dict_res_restraints[indx].second;
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

// calculated once and then stored
bool
coot::protein_geometry::atom_is_metal(mmdb::Atom *atom) const {

   // PDBv3 FIXME

   bool status = false;
   std::string atom_name(atom->GetAtomName());
   if (atom_name == "NA" || atom_name == "CA" || atom_name == "LI") {
      return true;
   } else {
      if (atom_name == "BE" || atom_name == "K" || atom_name == "RB") {
	 return true;
      } else {
	 if (atom_name == "SR" || atom_name == "CS" || atom_name == "BA") {
	    return true;
	 } else {
	    if (atom_name == "SC" || atom_name == "TI" || atom_name == "V" || atom_name == "CR") {
	       return true;
	    } else {
	       if (atom_name == "MN" || atom_name == "FE" || atom_name == "CO" || atom_name == "NI") {
		  return true;
	       } else {
		  if (atom_name == "CU" || atom_name == "ZN" || atom_name == "ZR" || atom_name == "MO") {
		     return true;
		  } else {
		     if (atom_name == "AG" || atom_name == "AU" || atom_name == "PT" || atom_name == "HG") {
			return true;
		     } else {
			if (atom_name == "OS" || atom_name == "PB" || atom_name == " K" || atom_name == " W") {
			   return true;
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   return status;
}


// expand to 4c, the atom_id, give that it should match an atom of the type comp_id.
// Used in chem mods, when we don't know the comp_id until residue modification time.
// 
std::string
coot::protein_geometry::atom_id_expand(const std::string &atom_id,
				       const std::string &comp_id,
				       int imol_enc) const {

   std::string s = atom_id;
   bool allow_minimal = true;
   int idx = get_monomer_restraints_index(comp_id, imol_enc, allow_minimal);
   if (idx != -1) {
      const coot::dictionary_residue_restraints_t &restraints =
	 dict_res_restraints[idx].second;
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

std::string
coot::protein_geometry::three_letter_code(const unsigned int &i) const {

   std::string r;
   if (i < dict_res_restraints.size()) {
      r = dict_res_restraints[i].second.residue_info.three_letter_code;
      if (r == "")
         r = dict_res_restraints[i].second.residue_info.comp_id;
   }
   return r;
}

// add "synthetic" 5 atom planar peptide restraint
void
coot::protein_geometry::add_planar_peptide_restraint() {

   std::string plane_id = "plane-5-atoms";
   mmdb::realtype dist_esd = 0.08; // was 0.11 (why?)

   std::string atom_id; 
   std::vector<std::pair<int, std::string> > v;
   v.push_back(std::pair<int, std::string> (1, "CA"));
   v.push_back(std::pair<int, std::string> (1, "C"));
   v.push_back(std::pair<int, std::string> (1, "O"));
   v.push_back(std::pair<int, std::string> (2, "N"));
   v.push_back(std::pair<int, std::string> (2, "CA"));

   for (unsigned int i=0; i<v.size(); i++) {
      link_add_plane("TRANS",  v[i].second, plane_id, v[i].first, dist_esd);
      link_add_plane("PTRANS", v[i].second, plane_id, v[i].first, dist_esd);
   }
}

bool
coot::protein_geometry::make_tight_planar_peptide_restraint() {

   std::string link_id("TRANS");
   std::string plane_id("plane-5-atoms");
   bool ifound = false;

   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {
      if (dict_link_res_restraints[i].link_id == link_id) { // e.g "TRANS"
	 std::vector<dict_link_plane_restraint_t>::iterator it;
	 for (it = dict_link_res_restraints[i].link_plane_restraint.begin();
	      it != dict_link_res_restraints[i].link_plane_restraint.end(); it++) {
	    if (it->plane_id == plane_id) {
	       it->set_dist_esd(0.03); // guess value
	       ifound = true;
	       break;
	    }
	 }
      }
   }
   return ifound;
}


void
coot::protein_geometry::remove_planar_peptide_restraint() {

   std::string link_id = "TRANS";
   std::string plane_id = "plane-5-atoms";
   unsigned int ifound = 0;

   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {
      if ((dict_link_res_restraints[i].link_id ==  "TRANS")  ||
          (dict_link_res_restraints[i].link_id == "PTRANS")) {

	 std::vector<coot::dict_link_plane_restraint_t>::iterator it;
	 for (it = dict_link_res_restraints[i].link_plane_restraint.begin();
	      it != dict_link_res_restraints[i].link_plane_restraint.end(); it++) {
	    if (it->plane_id == plane_id) {
	       ifound++;
	       if (0)
		  std::cout << "INFO:: before removal of plane3 TRANS has " 
			    << dict_link_res_restraints[i].link_plane_restraint.size()
			    << " plane restraints\n";
	       
	       // let's remove it
 	       dict_link_res_restraints[i].link_plane_restraint.erase(it);
	       
	       if (0)
		  std::cout << "INFO::  after removal of plane3 TRANS has " 
			    << dict_link_res_restraints[i].link_plane_restraint.size()
			    << " plane restraints\n";
	       break;
	    }
	 }
      }
      if (ifound == 2)
	 break;
   }
}

// Do the link restraints contain a planar peptide restraint?
bool
coot::protein_geometry::planar_peptide_restraint_state() const {

   std::string link_id = "TRANS";
   std::string plane_id = "plane-5-atoms";
   bool ifound = false;

   for (unsigned int i=0; i<dict_link_res_restraints.size(); i++) {
      if (dict_link_res_restraints[i].link_id == link_id) { // e.g "TRANS"
	 
	 std::vector<coot::dict_link_plane_restraint_t>::const_iterator it;
	 for (it = dict_link_res_restraints[i].link_plane_restraint.begin();
	      it != dict_link_res_restraints[i].link_plane_restraint.end(); it++) {
	    if (it->plane_id == plane_id) {
	       ifound = true;
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
   std::map<std::string,coot::dictionary_residue_restraints_t>::const_iterator it;

   std::vector<std::string> test_name_fragments =
      util::split_string(test_string, " ");

   for (it=simple_monomer_descriptions.begin();
	it!=simple_monomer_descriptions.end();
	it++) {
      std::string name_dc = util::downcase(it->second.residue_info.name);

      std::string::size_type ifound = std::string::npos;      

      for (unsigned int i=0; i<test_name_fragments.size(); i++) {

	 const std::string &test_string = test_name_fragments[i];
	 std::string test_string_dc = util::downcase(test_string);

	 ifound = name_dc.find(test_string_dc);
	 if (ifound == std::string::npos)
	    break;
      }
      if (ifound != std::string::npos) {

	 // 	 std::cout << "test_string :" << test_string << ": matched :"
	 // 		   << it->second.residue_info.comp_id << ": :"
	 // 		   << it->second.residue_info.name
	 // 		   << ":" << std::endl;

	 std::pair<std::string, std::string> p(it->second.residue_info.comp_id,
					       it->second.residue_info.name);
	 v.push_back(p);
      } else {
         // std::cout << "No match for " << it->first << " " << it->second.residue_info.name << std::endl;
      }
   }
   return v;
}


// replace (return 1)  or add (if not replacable) (return 0).
// 
bool
coot::protein_geometry::replace_monomer_restraints(std::string monomer_type,
						   int imol_enc,
						   const coot::dictionary_residue_restraints_t &mon_res_in) {
   bool s = false;

   dictionary_residue_restraints_t mon_res = mon_res_in;
   mon_res.assign_chiral_volume_targets();
   
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	 if (dict_res_restraints[i].first == imol_enc) {
	    dict_res_restraints[i].second = mon_res;
	    s = true;
	    break;
	 }
      }
   }

   if (!s) {
      std::pair<int, dictionary_residue_restraints_t> p(imol_enc, mon_res_in);
      dict_res_restraints.push_back(p);
   }
   return s;
}


// Keep everything that we have already, replace only those
// parts that are in mon_res_in.
// Used to update bond and angle restraints from Mogul.
// 
bool
coot::protein_geometry::replace_monomer_restraints_conservatively(std::string monomer_type, 
								  const dictionary_residue_restraints_t &mon_res) {

   bool s = false;
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == monomer_type) {
	 replace_monomer_restraints_conservatively_bonds( i, mon_res);
	 replace_monomer_restraints_conservatively_angles(i, mon_res);
	 s = true;
	 break;
      }
   }
   return s;
}

void
coot::protein_geometry::replace_monomer_restraints_conservatively_bonds(int irest,
									const dictionary_residue_restraints_t &mon_res) {

   // replace bonds
   // 
   for (unsigned int ibond=0; ibond<dict_res_restraints[irest].second.bond_restraint.size(); ibond++) { 
      for (unsigned int jbond=0; jbond<mon_res.bond_restraint.size(); jbond++) {

	 if (dict_res_restraints[irest].second.bond_restraint[ibond].atom_id_1_4c() ==
	     mon_res.bond_restraint[jbond].atom_id_1_4c()) {
	    if (dict_res_restraints[irest].second.bond_restraint[ibond].atom_id_2_4c() ==
		mon_res.bond_restraint[jbond].atom_id_2_4c()) {
	       dict_res_restraints[irest].second.bond_restraint[ibond] =
		  mon_res.bond_restraint[jbond];
	       break;
	    }
	 }

	 if (dict_res_restraints[irest].second.bond_restraint[ibond].atom_id_1_4c() ==
	     mon_res.bond_restraint[jbond].atom_id_2_4c()) {
	    if (dict_res_restraints[irest].second.bond_restraint[ibond].atom_id_2_4c() ==
		mon_res.bond_restraint[jbond].atom_id_1_4c()) {
	       dict_res_restraints[irest].second.bond_restraint[ibond] =
		  mon_res.bond_restraint[jbond];
	       break;
	    }
	 }
      }
   }
}




void
coot::protein_geometry::replace_monomer_restraints_conservatively_angles(int irest,
									 const dictionary_residue_restraints_t &mon_res) {

   for (unsigned int iangle=0; iangle<dict_res_restraints[irest].second.angle_restraint.size(); iangle++) { 
      for (unsigned int jangle=0; jangle<mon_res.angle_restraint.size(); jangle++) {

	 // middle atom the same
	 // 
	 if (dict_res_restraints[irest].second.angle_restraint[iangle].atom_id_2_4c() ==
	     mon_res.angle_restraint[jangle].atom_id_2_4c()) {

	    // check for either way round of 1 and 3:

	    if (dict_res_restraints[irest].second.angle_restraint[iangle].atom_id_1_4c() ==
		mon_res.angle_restraint[jangle].atom_id_1_4c()) {
	       if (dict_res_restraints[irest].second.angle_restraint[iangle].atom_id_3_4c() ==
		   mon_res.angle_restraint[jangle].atom_id_3_4c()) {
		  dict_res_restraints[irest].second.angle_restraint[iangle] =
		     mon_res.angle_restraint[jangle];
	       }
	    }

	    if (dict_res_restraints[irest].second.angle_restraint[iangle].atom_id_1_4c() ==
		mon_res.angle_restraint[jangle].atom_id_3_4c()) {
	       if (dict_res_restraints[irest].second.angle_restraint[iangle].atom_id_3_4c() ==
		   mon_res.angle_restraint[jangle].atom_id_1_4c()) {
		  dict_res_restraints[irest].second.angle_restraint[iangle] =
		     mon_res.angle_restraint[jangle];
	       }
	    }
	 }
      }
   }
}


std::vector<std::string>
coot::protein_geometry::monomer_types() const {
   std::vector<std::string> v;
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      v.push_back(dict_res_restraints[i].second.residue_info.comp_id);
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
coot::protein_geometry::get_group(mmdb::Residue *r) const {

   std::string res_name = r->GetResName();
   return get_group(res_name);
}


std::string
coot::protein_geometry::get_group(const std::string &res_name_in) const {

   bool found = false;
   std::string group;
   std::string res_name = res_name_in;

   // 20250331-PE Why would I do this?
   //             Comment it out.
   // if (res_name.length() > 3)
   //    res_name = res_name.substr(0,2);

   unsigned int s = size(); // fails if the protein_geometry pointer was not valid
   for (unsigned int i=0; i<s; i++) {
      if (three_letter_code(i) == res_name) {
	 found = true;
	 group = (*this)[i].second.residue_info.group;
	 break;
      }
   }

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      if (dict_res_restraints[i].second.residue_info.comp_id == res_name) {
	 found = true;
	 group = dict_res_restraints[i].second.residue_info.group;
	 break;
      }
   }

   if (! found) {
      std::string s = "WARNING:: get_group(): No dictionary group found for residue type :";
      s += res_name;
      s += ":";
      throw std::runtime_error(s);
   }
   return group;
}


std::vector<std::string>
coot::protein_geometry::residue_names_with_no_dictionary(mmdb::Manager *mol, int imol_no) const {

   std::vector<std::string> v;

   if (mol) {
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         std::set<std::string> already_tested_names;
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) {
	       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	       if (residue_p) {
		  std::string residue_name = residue_p->GetResName();
                  if (already_tested_names.find(residue_name) == already_tested_names.end()) {
                     if (! have_dictionary_for_residue_type_no_dynamic_add(residue_name, imol_no))
                        if (std::find(v.begin(), v.end(), residue_name) == v.end())
                           v.push_back(residue_name);
                     already_tested_names.insert(residue_name);
                  }
	       }
	    }
	 }
      }
   }
   return v;
}

bool
coot::protein_geometry::read_extra_dictionaries_for_molecule(mmdb::Manager *mol, int imol, int *read_number_p) {

   if (! mol) return false;

   std::vector<std::string> v = residue_names_with_no_dictionary(mol, imol);
   if (false) {
      std::cout << "-------------- debug residue names with no dictionary " << std::endl;
      for (unsigned int i=0; i<v.size(); i++) {
         std::cout << "           " << v[i] << std::endl;
      }
   }
   bool success = true;
   for (std::size_t i=0; i<v.size(); i++) {
      const std::string &rn = v[i];
      int success_for_residue = try_dynamic_add(rn, *read_number_p);
      if (success_for_residue == 0)
	 success = false;
      *read_number_p += 1;
   }

   return success;
}



// optional arg: bool try_autoload_if_needed=true.
// optional arg: float b_factor.
mmdb::Residue *
coot::protein_geometry::get_residue(const std::string &comp_id, int imol_enc,
				    bool idealised_flag,
				    bool try_autoload_if_needed, // default true
                                    float b_factor // default 20
                                    ) {

   auto debug_residue = [] (mmdb::Residue *residue_p) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::cout << "debug:: in get_residue(): atom " << iat << " " << at-> name
                      << " at " << at->x << " " << at->y << " " << at->z << std::endl;
         }
      }
   };

   // If the coordinates for the model are (0,0,0) then this function
   // returns a null.

   bool debug = false;

   mmdb::Residue *residue_p = NULL;

   // might use try_dynamic_add (if needed).
   bool r = have_dictionary_for_residue_type(comp_id, imol_enc, try_autoload_if_needed);

   if (debug)
      std::cout << "------------------ in pr::get_residue() have_dictionary_for_residue_type() returns  "
                << r << std::endl;
   if (r) {
      for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
         const dictionary_residue_restraints_t &rest = dict_res_restraints[i].second;
         if (debug)
            std::cout << "   testing comp_id " << rest.residue_info.comp_id << std::endl;
         if (rest.residue_info.comp_id == comp_id) {
            int imol_for_dict = dict_res_restraints[i].first;
            if (imol_for_dict == imol_enc) {
               residue_p = rest.GetResidue(idealised_flag, b_factor);
               // debug_residue(residue_p);
               break;
            }
         }
      }
   }
   if (debug)
      std::cout << "------------------ in pr::get_residue() returns " << residue_p << std::endl;

   return residue_p;
}



mmdb::Manager *
coot::protein_geometry::mol_from_dictionary(const std::string &three_letter_code,
					    int imol_enc,
					    bool idealised_flag) {

   auto debug_mol = [] (mmdb::Manager *mol) {
      for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     for (int iat=0; iat<n_atoms; iat++) {
                        mmdb::Atom *at = residue_p->GetAtom(iat);
                        if (! at->isTer()) {
                           std::cout << "pg::mol_from_dictionary(): atom " << iat << " " << at->name
                                     << " at " << at->x << " "  << at->y << " " << at->z
                                     << std::endl;
                        }
                     }
                  }
               }
            }
         }
      }
   };

   mmdb::Manager *mol = NULL;
   mmdb::Residue *residue_p = get_residue(three_letter_code, imol_enc, idealised_flag);

   if (residue_p) {
      mmdb::Chain *chain_p = new mmdb::Chain;
      chain_p->SetChainID("A");
      chain_p->AddResidue(residue_p);
      mmdb::Model *model_p = new mmdb::Model;
      model_p->AddChain(chain_p);
      mol = new mmdb::Manager;
      mol->AddModel(model_p);
   } else {
      std::cout << "WARNING:: protein_geometry::mol_from_dictionary(): "
                << "Null residue in mol_from_dictionary() for "
                << three_letter_code << std::endl;
   }
   // debug_mol(mol);
   return mol;
}

mmdb::Manager *
coot::protein_geometry::mol_from_dictionary(int monomer_index,
					    int imol_enc,
					    bool idealised_flag) {

   mmdb::Manager *mol = NULL;
   mmdb::Residue *residue_p = 0;
   float b_factor = 30.0;
   int r_size = dict_res_restraints.size();

   if (monomer_index >= 0)
      if (monomer_index < r_size)
	 residue_p = dict_res_restraints[monomer_index].second.GetResidue(idealised_flag, b_factor);

   if (residue_p) {
      mmdb::Chain *chain_p = new mmdb::Chain;
      chain_p->SetChainID("A");
      chain_p->AddResidue(residue_p);
      mmdb::Model *model_p = new mmdb::Model;
      model_p->AddChain(chain_p);
      mol = new mmdb::Manager;
      mol->AddModel(model_p);
   } else {
      std::cout << "WARNING:: Null residue in mol_from_dictionary() for idx "
		<< monomer_index << std::endl;
   }

   std::cout << "DEBUG:: mol_from_dictionary() returns " << mol << std::endl;
   return mol;
}


// delete comp_id from dict_res_restraints (if it exists there).
bool
coot::protein_geometry::delete_mon_lib(const std::string &comp_id, int imol_enc) {

   bool deleted = false;
   std::vector<std::pair<int, coot::dictionary_residue_restraints_t> >::iterator it;
   for (it=dict_res_restraints.begin(); it!=dict_res_restraints.end(); ++it) {
      if (it->second.residue_info.comp_id == comp_id) {
	 if (it->first == imol_enc) {
	    dict_res_restraints.erase(it);
	    deleted = 1;
	    break;
	 }
      }
   }
   
   return deleted;
} 

bool
coot::protein_geometry::OXT_in_residue_restraints_p(const std::string &residue_type) const {

   bool r = 0;
   std::pair<bool, coot::dictionary_residue_restraints_t> p =
      get_monomer_restraints(residue_type, IMOL_ENC_ANY);
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



// Find the bonded neighbours of the given atoms - throw an
// exception if residue name is not in dictionary.
// 
std::vector<std::string>
coot::protein_geometry::get_bonded_neighbours(const std::string &residue_name,
					      int imol_enc,
					      const std::string &atom_name_1,
					      const std::string &atom_name_2, 
					      bool also_2nd_order_neighbs_flag) const {

   std::vector<std::string> v;

   std::vector<std::string> v_2nd_order; // only filled for triple bonds

   std::pair<bool, coot::dictionary_residue_restraints_t> restraints =
      get_monomer_restraints_at_least_minimal(residue_name, imol_enc);

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
		     coot::protein_geometry::get_bonded_neighbours(residue_name,  imol_enc,
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
		     coot::protein_geometry::get_bonded_neighbours(residue_name, imol_enc,
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
		     coot::protein_geometry::get_bonded_neighbours(residue_name, imol_enc,
								   atom_name_1, at_name);
		  for (unsigned int in=0; in<nv.size(); in++)
		     v_2nd_order.push_back(nv[in]);
	       }
	    }
	 if (restraints.second.bond_restraint[i].atom_id_2_4c() == atom_name_2) {
	    if (restraints.second.bond_restraint[i].atom_id_1_4c() != atom_name_1) {
	       // std::cout << " adding d " << restraints.second.bond_restraint[i].atom_id_1_4c() << std::endl;
	       std::string at_name = restraints.second.bond_restraint[i].atom_id_1_4c();
	       v.push_back(at_name);
	       if (also_2nd_order_neighbs_flag) { 
		  std::vector<std::string> nv = 
		     coot::protein_geometry::get_bonded_neighbours(residue_name, imol_enc,
								   atom_name_1, at_name);
		  for (unsigned int in=0; in<nv.size(); in++)
		     v_2nd_order.push_back(nv[in]);
	       }
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

std::vector<std::pair<std::string, std::string> >
coot::protein_geometry::get_bonded_and_1_3_angles(const std::string &comp_id, int imol) const {

   std::vector<std::pair<std::string, std::string> > v;
   int idx = get_monomer_restraints_index(comp_id, imol, true);
   if (idx != -1) {
      const dictionary_residue_restraints_t &rest = dict_res_restraints[idx].second;
      for (unsigned int i=0; i<rest.bond_restraint.size(); i++) {
	 std::pair<std::string, std::string> p(rest.bond_restraint[i].atom_id_1_4c(),
					       rest.bond_restraint[i].atom_id_2_4c());
	 v.push_back(p);
      }
      for (unsigned int i=0; i<rest.angle_restraint.size(); i++) {
	 std::pair<std::string, std::string> p(rest.angle_restraint[i].atom_id_1_4c(),
					       rest.angle_restraint[i].atom_id_3_4c());
	 v.push_back(p);
      }
   }
   return v;
} 



// should this return a set - or filter the results for duplicates?
//
std::vector<std::string>
coot::protein_geometry::monomer_restraints_comp_ids() const {

   std::vector<std::string> v;
   for (unsigned int i=0; i<dict_res_restraints.size(); i++)
      v.push_back(dict_res_restraints[i].second.residue_info.comp_id);
   return v;
}


// can throw a std::runtime_error
std::string
coot::protein_geometry::Get_SMILES_for_comp_id(const std::string &comp_id, int imol_enc) const {

   bool found = false;
   std::string s;
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {

      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {
         if (dict_res_restraints[i].first == imol_enc) {

            unsigned int nd = dict_res_restraints[i].second.descriptors.descriptors.size();
            for (unsigned int idesc=0; idesc<nd; idesc++) {
               if (dict_res_restraints[i].second.descriptors.descriptors[idesc].type == "SMILES_CANONICAL") {
                  s = dict_res_restraints[i].second.descriptors.descriptors[idesc].descriptor;
                  found = true;
                  break;
               }
            }
         }
         if (found)
            break;
      }
   }

   if (! found){
      // check non-canonical
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {

      if (dict_res_restraints[i].second.residue_info.comp_id == comp_id) {

	 unsigned int nd = dict_res_restraints[i].second.descriptors.descriptors.size();
	 for (unsigned int idesc=0; idesc<nd; idesc++) {
	    if (dict_res_restraints[i].second.descriptors.descriptors[idesc].type == "SMILES") {
	       s = dict_res_restraints[i].second.descriptors.descriptors[idesc].descriptor;
	       found = true;
	       break;
	    }
	 }
      }
      if (found)
	 break;
   }
}

   if (! found){
      std::string mess = "No SMILES in dictionary for ";
      mess += comp_id;
      throw (std::runtime_error(mess));
   }
   return s;
} 

      



// This uses have_dictionary_for_residue_type() (and thus
// try_dynamic_add() if needed).
// 
// Return -1 if residue type not found.
// 
int
coot::protein_geometry::n_hydrogens(const std::string &residue_type) {

   int n_hydrogens = -1;
   
   std::pair<bool, dictionary_residue_restraints_t> r =
      get_monomer_restraints(residue_type, IMOL_ENC_ANY);

   if (r.first) {
      n_hydrogens = 0; // not not-found
      for (unsigned int iat=0; iat<r.second.atom_info.size(); iat++) { 
	 if (r.second.atom_info[iat].type_symbol == " H")
	    n_hydrogens++;
	 else 
	 if (r.second.atom_info[iat].type_symbol == "H")
	    n_hydrogens++;
      }
   }
   return n_hydrogens;
} 

// This uses have_dictionary_for_residue_type() (and thus
// try_dynamic_add() if needed).
// 
// Return -1 if residue type not found.
// 
int
coot::protein_geometry::n_non_hydrogen_atoms(const std::string &residue_type) {

   int n_atoms = -1;
   
   std::pair<bool, dictionary_residue_restraints_t> r =
      get_monomer_restraints(residue_type, IMOL_ENC_ANY);

   if (r.first) {
      n_atoms = 0; // not not-found
      for (unsigned int iat=0; iat<r.second.atom_info.size(); iat++) { 
	 if (r.second.atom_info[iat].type_symbol != " H")
	    if (r.second.atom_info[iat].type_symbol != "H")
	       n_atoms++;
      }
   }
   return n_atoms;
}




// debug
void
coot::protein_geometry::debug() const {

   std::cout << "### debug(): " << dict_res_restraints.size() << " entries " << std::endl;
   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      int imol = dict_res_restraints[i].first;
      std::string imol_str = std::string("          ") + util::int_to_string(imol);
      if (imol == IMOL_ENC_ANY)
	 imol_str = "IMOL_ENC_ANY";
      if (imol == IMOL_ENC_AUTO)
	 imol_str = "IMOL_ENC_AUTO";
      if (imol == IMOL_ENC_UNSET)
	 imol_str = "IMOL_ENC_UNSET";
      std::cout << "     " << i << " imol: " << imol_str << " \""
		<< dict_res_restraints[i].second.residue_info << "\"" << std::endl;
   }
}

void coot::protein_geometry::all_plane_restraints_to_improper_dihedrals() {

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      dict_res_restraints[i].second.improper_dihedral_restraint.clear();
      for (unsigned int j=0; j<dict_res_restraints[i].second.plane_restraint.size(); j++) {
         std::vector<atom_name_quad> qs = dict_res_restraints[i].second.plane_restraint_to_improper_dihedrals(j);
         for(unsigned int k=0; k<qs.size(); k++) {
            dict_improper_dihedral_restraint_t r(qs[k].atom_name(0), qs[k].atom_name(1), qs[k].atom_name(2), qs[k].atom_name(3));
            dict_res_restraints[i].second.improper_dihedral_restraint.push_back(r);
         }
      }
   }

}

void coot::protein_geometry::delete_plane_restraints() {

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      dict_res_restraints[i].second.plane_restraint.clear();
   }

}



void
coot::protein_geometry::print_dictionary_store() const {

   for (unsigned int i=0; i<dict_res_restraints.size(); i++) {
      int imol_enc = dict_res_restraints[i].first;
      const auto &rest = dict_res_restraints[i].second;
      std::cout << i << " " << rest.residue_info << " for imol-enc: " << imol_enc << std::endl;
   }

}


std::vector<std::pair<std::string, std::string> >
coot::protein_geometry::get_acedrg_atom_types(const std::string &comp_id,
                                              int imol_enc) const {

   std::vector<std::pair<std::string, std::string> > v;

   std::pair<bool, dictionary_residue_restraints_t> r =
      get_monomer_restraints_internal(comp_id, imol_enc, false);
   if (r.first) {
      const auto &restraints = r.second;
      std::vector<dict_atom>::const_iterator it;
      for (it=restraints.atom_info.begin(); it!=restraints.atom_info.end(); ++it) {
         if (! it->acedrg_atom_type.empty()) {
            auto pair = std::make_pair(it->atom_id, it->acedrg_atom_type);
            v.push_back(pair);
         }
      }
   }
   return v;

}
