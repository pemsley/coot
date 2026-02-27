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
#define XDATADIR "C:/coot/share"
#define PKGDATADIR XDATADIR
// stop using these, use win-compat functions
// #define S_ISDIR(m)  (((m) & S_IFMT) == S_IFDIR)
// #define S_ISREG(m)  (((m) & S_IFMT) == S_IFREG)
#endif

#include "clipper/core/clipper_util.h"

#include "compat/coot-sysdep.h"
#include "utils/win-compat.hh"

#include "lbg-graph.hh"

#include "utils/logging.hh"
extern logging logger;

void
coot::protein_geometry::set_only_bonds(int dict_idx) {

   if (dict_idx != -1) {
      int drr_size = dict_res_restraints.size(); // type change
      if (dict_idx < drr_size) {
         dictionary_residue_restraints_t &dict = dict_res_restraints[dict_idx].second;
         std::map<std::string, std::vector<std::pair<unsigned int, std::string> > > atom_name_map;
         for (unsigned int i=0; i<dict.bond_restraint.size(); i++) {
            const auto &bond = dict.bond_restraint[i];
            if (bond.type() == "single") {
               const std::string &atom_name_1 = bond.atom_id_1();
               const std::string &atom_name_2 = bond.atom_id_2();
               bool is_H_1 = dict.is_hydrogen(atom_name_1);
               bool is_H_2 = dict.is_hydrogen(atom_name_2);
               if (! is_H_1 && !is_H_2) {
                  std::pair<unsigned int, std::string> p1(i, "first");
                  std::pair<unsigned int, std::string> p2(i, "second");
                  atom_name_map[bond.atom_id_1()].push_back(p1);
                  atom_name_map[bond.atom_id_2()].push_back(p2);
               }
            }
         }
         std::map<std::string, std::vector<std::pair<unsigned int, std::string> > >::const_iterator it;
         for (it=atom_name_map.begin(); it!=atom_name_map.end(); ++it) {
            const std::vector<std::pair<unsigned int, std::string> > &v = it->second;
            if (v.size() == 1) {
               const auto &atom_name = it->first;
               unsigned int bond_index = v[0].first;
               const std::string &pos  = v[0].second;
               dict.bond_restraint[bond_index].set_only_bond(pos, true);
               if (false)
                  std::cout << "debug:: set_only_bond() in " << dict.residue_info.comp_id
                            << " atom " << atom_name << " " << pos
                            << " " << dict.bond_restraint[bond_index].type()
                            << " has only one non-Hydrogen bond" << std::endl;
            }
         }
      }
   }
}

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
            // std::cout << "There are " << ciffile.GetNofData() << " data in "
            //                       << ciffilename << std::endl;
            logger.log(log_t::INFO, "There are ", std::to_string(ciffile.GetNofData()),
                       " data in ", ciffilename);

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

               // All categories have loops (AFAICS).
               // std::cout << "DEBUG:: got category: " << cat_name << std::endl;

               mmdb::mmcif::PLoop mmCIFLoop = data->GetLoop(cat_name.c_str() );

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

                  if (cat_name == "_lib") {
                     mmdb::mmcif::Struct *structure = data->GetStructure(cat_name.c_str());
                     if (structure) {
                        parse_lib_info(structure);
                        handled = true; // hack for now - so that we don't get the warning message
                     }
                  }

                  if (cat_name == "_gphl_chem_comp_info") {
                     mmdb::mmcif::Struct *structure = data->GetStructure(cat_name.c_str());
                     if (structure) {
                        gphl_chem_comp_info(structure, imol_enc);
                        handled = true;
                     }
                  }

                  if (! handled)   // this can happen if there is not an atom loop, e.g. dictionary
                                   // with one atom e.g. AM.cif (Americium ion)
                     std::cout << "WARNING:: in init_refmac_mon_lib() unhandled category \""
                               << cat_name << "\" file: " << ciffilename << std::endl;

               } else {


                  // We currently want to stop adding chem comp info
                  // if the chem_comp info comes from mon_lib_list.cif:
                  if (cat_name == "_chem_comp") {
                     if (read_number_in != coot::protein_geometry::MON_LIB_LIST_CIF) {
                        comp_id_2 = chem_comp(mmCIFLoop, imol_enc);
                     } else {
                        comp_id_2 = simple_mon_lib_chem_comp(mmCIFLoop, imol_enc);
                     }
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

                  // PDBe depiction -  needs _pdbe_chem_comp_bond_depiction parser also
                  if (cat_name == "_pdbe_chem_comp_atom_depiction")
                     pdbe_chem_comp_atom_depiction(mmCIFLoop, imol_enc);

                  if (cat_name == "_pdbx_chem_comp_description_generator")
                     pdbx_chem_comp_description_generator(mmCIFLoop, imol_enc);

                  if (cat_name == "_chem_comp_acedrg")
                     chem_comp_acedrg(mmCIFLoop, imol_enc);
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

   set_only_bonds(rmit.monomer_idx);

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
coot::protein_geometry::parse_lib_info(mmdb::mmcif::PStruct structure) {

   // 20220223-PE fill this at some stage

   // example:
   // data_lib
   // _lib.name mon_lib
   // _lib.version 5.60
   //  _lib.update 16/02/22

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
//              << ri.comp_id << ": :"
//              << ri.three_letter_code << ": " 
//              << ri.group << "] " 
//              << " for key :"
//              << comp_id << ":" << " simple_monomer_descriptions size() "
//              << simple_monomer_descriptions.size() 
//              << std::endl;

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
      char *s = mmCIFLoop->GetString("id", j, ierr); // 20220223-PE this might return null
      std::string three_letter_code;
      std::string name;
      std::string group; // e.g. "L-peptide"
      int number_atoms_all;
      int number_atoms_nh;
      std::string description_level = "None";

      if (ierr == 0) {
         int ierr_tot = 0;
         if (s) { // 20220223-PE add protection for null id extraction.
            comp_id = s;
            s = mmCIFLoop->GetString("three_letter_code", j, ierr);
            ierr_tot += ierr;
            if (s)
               three_letter_code = s;
            else {
               three_letter_code = "";
               //             std::cout << "WARNING:: failed to get 3-letter code for comp_id: "
               //                       << comp_id << " error: " << ierr << std::endl;
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
   }
   return comp_id;
}

void
coot::protein_geometry::pdbe_chem_comp_atom_depiction(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc) {

   int ierr = 0;
   std::vector<depiction_atom_t> dav;
   std::set<std::string> comp_id_set;
   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {
      std::string atom_id, element;
      double model_Cartn_x, model_Cartn_y;
      int pdbx_ordinal;
      char *s = mmCIFLoop->GetString("comp_id", j, ierr);
      if (s) {
         std::string comp_id = std::string(s);
         comp_id_set.insert(comp_id);
      }
      s = mmCIFLoop->GetString("atom_id", j, ierr);
      if (s) {
         atom_id = std::string(s);
      }
      s = mmCIFLoop->GetString("element", j, ierr);
      if (s) {
         element = std::string(s);
      }
      int ierr_x   = mmCIFLoop->GetReal(model_Cartn_x, "model_Cartn_x", j);
      int ierr_y   = mmCIFLoop->GetReal(model_Cartn_y, "model_Cartn_y", j);
      int ierr_ord = mmCIFLoop->GetInteger(pdbx_ordinal, "pdbx_ordinal", j);

      if (ierr == 0 && ierr_x == 0 && ierr_y == 0 && ierr_ord == 0) {
         depiction_atom_t da(atom_id, element, model_Cartn_x, model_Cartn_y, pdbx_ordinal);
         dav.push_back(da);
      }
   }
   if (! dav.empty()) {
      if (comp_id_set.size() == 1) {
         std::string comp_id = *comp_id_set.begin();
         chem_comp_atom_depiction_t d(comp_id, dav);
         int idx = get_monomer_restraints_index(comp_id, imol_enc, true);
         if (idx >= 0) {
            dict_res_restraints[idx].second.depiction = d;
            std::cout << "debug:: pdbe_chem_comp_atom_depiction() added depiction of "
                      << dav.size() << " atoms " << std::endl;
         }
      }
   }
}

void
coot::protein_geometry::pdbx_chem_comp_description_generator(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc) {

   int ierr = 0;
   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {
      std::string comp_id;
      std::string program_name;
      std::string program_version;
      std::string descriptor;
      char *s = mmCIFLoop->GetString("comp_id", j, ierr);
      if (s) {
         comp_id = std::string(s);
      }
      s = mmCIFLoop->GetString("program_name", j, ierr);
      if (s) {
         program_name = std::string(s);
      }
      s = mmCIFLoop->GetString("program_version", j, ierr);
      if (s) {
         program_version = std::string(s);
      }
      s = mmCIFLoop->GetString("descriptor", j, ierr);
      if (s) {
         descriptor = std::string(s);
      }
      if (ierr == 0) {
         pdbx_chem_comp_description_generator_t dg(program_name, program_version, descriptor);
         int idx = get_monomer_restraints_index(comp_id, imol_enc, true);
         if (idx >= 0) {
            dict_res_restraints[idx].second.description_generation = dg;
         }
      }
   }
}

void
coot::protein_geometry::gphl_chem_comp_info(mmdb::mmcif::PStruct structure, int imol_enc) {

   gphl_chem_comp_info_t gphl_chem_comp_info;
   std::vector<std::string> keys = {
      "comp_id",
      "arguments",
      "run_date",
      "grade2_version",
      "grade2_date",
      "rdkit_version",
      "input_from",
      "input_data",
      "input_inchi",
      "input_inchikey",
      "input_inchi_match",
      "rdkit_sanitization_problem",
      "partial_charges_source",
      "force_field",
      "initial_energy",
      "final_energy",
      "mogul_version",
      "mogul_data_libraries",
      "csd_version",
      "csd_python_api",
      "coordinates_source",
      "geometry_optimize_program",
      "geometry_optimize_steps",
      "geometry_optimize_initial_function",
      "geometry_optimize_final_function",
      "geometry_optimize_initial_rms_gradient",
      "geometry_optimize_final_rms_gradient",
      "geometry_optimize_initial_rms_bond_deviation",
      "geometry_optimize_final_rms_bond_deviation",
      "elapsed_seconds"};

   int n_tags = structure->GetNofTags();
   for (int itag=0; itag<n_tags; itag++) {
      std::string tag   = structure->GetTag(itag);
      std::string field = structure->GetField(itag);
      gphl_chem_comp_info.add(tag, field);
   }

   int idx = gphl_chem_comp_info.get_index("comp_id");
   if (idx >= 0) {
      const std::string &comp_id = gphl_chem_comp_info[idx].second;
      int idx_rest = get_monomer_restraints_index(comp_id, imol_enc, true);
      if (idx_rest >= 0) {
         dict_res_restraints[idx_rest].second.gphl_chem_comp_info = gphl_chem_comp_info;
         std::cout << "debug:: adding a gphl info for " << comp_id << " of size " << gphl_chem_comp_info.info.size() << std::endl;
      }
   }
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
//              std::cout << "comp_atom_pad_atom_name: in :" << atom_id << ": out :"
//                        << padded_name << ":" << std::endl;
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
      dict_bond_restraint_t::bond_length_type_t blt(dict_bond_restraint_t::UNKNOWN); // nuclear or electron

   
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
	 s = mmCIFLoop->GetString("aromatic", j, ierr_optional);
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

         mmdb::realtype value_dist_nuclear     = -1.0;
         mmdb::realtype value_dist_nuclear_esd = -1.0;
         int ierr_nuc     = mmCIFLoop->GetReal(value_dist_nuclear,     "value_dist_nucleus",     j);
         int ierr_nuc_esd = mmCIFLoop->GetReal(value_dist_nuclear_esd, "value_dist_nucleus_esd", j);

	 if (ierr_tot == 0) {

            // if value_dist_nucleus_esd is negative, it's treated as unset
            //
	    mon_lib_add_bond(comp_id, imol_enc, atom_id_1, atom_id_2,
			     type,
                             value_dist, value_dist_esd,
                             value_dist_nuclear, value_dist_nuclear_esd,
                             aromaticity, blt);
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

void
coot::protein_geometry::chem_comp_acedrg(mmdb::mmcif::PLoop mmCIFLoop, int imol_enc) {

   std::string comp_id;
   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {
      int ierr = 0;
      int ierr_tot = 0;
      std::string atom_id;
      std::string atom_type;
      char *s = mmCIFLoop->GetString("comp_id", j, ierr);
      if (! ierr)
         if (s)
            comp_id = s;
      ierr_tot += ierr;
      s = mmCIFLoop->GetString("atom_id", j, ierr);
      if (! ierr) {
         atom_id = s;
      }
      ierr_tot += ierr;
      s = mmCIFLoop->GetString("atom_type", j, ierr);
      if (! ierr) {
         atom_type = s;
      }
      ierr_tot += ierr;

      if (ierr_tot == 0) {
         mon_lib_add_acedrg_atom_type(comp_id, imol_enc, atom_id, atom_type);
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
              list_chem_mod mod(id, name, comp_id, group_id);
              chem_mod_vec.push_back(mod);
              n_chem_mods++;
      }
   }
   return n_chem_mods;
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

   bool debug = false;

   // std::string dir = DATADIR;
   std::string pkg_data_dir = package_data_dir(); // xxx/share/coot
   std::string hardwired_default_place = util::append_dir_dir(pkg_data_dir, "lib");
   bool using_clibd_mon = false;

   if (debug) {
      std::cout << "DEBUG:: init_standard(): pkg_data_dir: " << pkg_data_dir << std::endl;
      std::cout << "DEBUG:: init_standard(): hardwired_default_place: " << hardwired_default_place << std::endl;
   }

   std::string mon_lib_dir;
   short int env_dir_fails = 0;
   int istat;

   // struct stat buf;
   char *cmld = NULL;

   char *s = getenv("COOT_REFMAC_LIB_DIR");
   // it is not refmac, it's CCP4
   if (! s)
      s = getenv("COOT_MONOMER_LIB_DIR");
   if (! s)
      s = getenv("COOT_CCP4_LIB_DIR");
   if (s) {
      if (is_dir_or_link(s)) {
          mon_lib_dir = s;
      } else {
          env_dir_fails = 1;
          std::cout << "WARNING:: Coot REFMAC dictionary override COOT_REFMAC_LIB_DIR "
                    << s << " " << "failed to find the monomer library " << std::endl;
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
            if (verbose_mode)
               // std::cout << "INFO:: Using Standard CCP4 Refmac dictionary from"
               //           << " CLIBD_MON: " << s << std::endl;
               logger.log(log_t::INFO, "Using Standard CCP4 Refmac dictionary from CLIBD_MON: " + std::string(s));
            mon_lib_dir = s;
            using_clibd_mon = true;
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
            if (verbose_mode)
               // std::cout << "INFO:: Using Standard CCP4 Refmac dictionary: "
               //           << s << std::endl;
               logger.log(log_t::INFO, "Using Standard CCP4 Refmac dictionary: " + std::string(s));
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
                     std::cout << "WARNING:: init_standard(): env var COOT_PREFIX set, but no dictionary lib found\n";
                     std::cout << "WARNING:: init_standard(): env var COOT_PREFIX was set to \"" << s << "\"\n";
                  }
               } else {
                  std::cout << "WARNING:: env var COOT_PREFIX not set" << std::endl;
                  mon_lib_dir.clear();
               }
            }
         }
      }
   }

   if (debug)
      std::cout << "DEBUG:: Here with mon_lib_dir set to " << mon_lib_dir << std::endl;

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

      if (false)
         std::cout << "calling init_refmac_mon_lib() on" << mon_lib_cif << std::endl;
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
   } else {
      std::cout << "WARNING: Failed to read restraints dictionary. "
                << std::endl;
   }


   return read_number;
}


int
coot::protein_geometry::refmac_monomer(const std::string &dir, // dir
                                       const std::string &protein_mono) { // extra path to file

   int imol_enc = IMOL_ENC_ANY; // maybe pass this?

   std::string filename = util::append_dir_file(dir, protein_mono);
   if (is_regular_file(filename)) {
      init_refmac_mon_lib(filename, read_number, imol_enc);
      read_number++;
   } else {
      if (coot::file_exists(filename)) {
         std::cout << "WARNING:: file " << filename << " is not a regular file"
                   << std::endl;
      } else {
         std::cout << "WARNING:: file " << filename << " does not exist"
                   << std::endl;
      }
   }
   return read_number;
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
               std::string up_type_symbol = util::upcase(util::remove_whitespace(atom_info[i].type_symbol));
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

      if (bond_restraint.size() > 0) {
         // bool nuclear_distances_flag = false;
         rc = mmCIFData->AddLoop("_chem_comp_bond", mmCIFLoop);
         if (rc == mmdb::mmcif::CIFRC_Ok || rc == mmdb::mmcif::CIFRC_Created) {
            // std::cout << " number of bonds: " << bond_restraint.size() << std::endl;

            // nuclear_distances_flag means that we only have one distance - and it's the
            // nuclear distance. So we need to "invent" non-nuclear distance for bonds
            // to hydrogen atoms
            for (unsigned int i=0; i<bond_restraint.size(); i++) {

               const dict_bond_restraint_t &br = bond_restraint[i];
               std::string value_dist("value_dist");
               std::string value_dist_esd("value_dist_esd");
               // std::cout << "ading bond number " << i << std::endl;
               const char *ss = residue_info.comp_id.c_str();
               mmCIFLoop->PutString(ss, "comp_id", i);
               std::string id_1 = util::remove_whitespace(br.atom_id_1_4c());
               std::string id_2 = util::remove_whitespace(br.atom_id_2_4c());
               std::string qan_1 = quoted_atom_name(id_1);
               std::string qan_2 = quoted_atom_name(id_2);
               ss = id_1.c_str();
               mmCIFLoop->PutString(id_1.c_str(), "atom_id_1", i);
               ss = id_2.c_str();
               mmCIFLoop->PutString(id_2.c_str(), "atom_id_2", i);
               std::string bond_type = bond_restraint[i].type();
               mmCIFLoop->PutString(bond_type.c_str(), "type", i);
               try {

                  if (nuclear_distances_flag) {

                     // Non-nuclear

                     float v = bond_restraint[i].value_dist();
                     if (is_bond_to_hydrogen_atom(br))
                        v /= 1.08;
                     mmCIFLoop->PutReal(v, value_dist.c_str(), i, 5);
                     v = bond_restraint[i].value_esd();
                     mmCIFLoop->PutReal(v, value_dist_esd.c_str(), i, 3);

                     // Nuclear

                     value_dist     = "value_dist_nucleus";
                     value_dist_esd = "value_dist_nucleus_esd";

                     v = bond_restraint[i].value_dist();
                     mmCIFLoop->PutReal(v, value_dist.c_str(), i, 5);
                     v = bond_restraint[i].value_esd();
                     mmCIFLoop->PutReal(v, value_dist_esd.c_str(), i, 3);

                  } else {

                     float v = bond_restraint[i].value_dist();
                     mmCIFLoop->PutReal(v, value_dist.c_str(), i, 5);
                     v = bond_restraint[i].value_esd();
                     mmCIFLoop->PutReal(v, value_dist_esd.c_str(), i, 3);

                  }

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

      if (angle_restraint.size() > 0) {
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

               // note strange string usage is due to static analysis use-after-free error report
               // don't reuse ss
               std::string id_1 = util::remove_whitespace(torsion_restraint[i].atom_id_1_4c());
               std::string id_2 = util::remove_whitespace(torsion_restraint[i].atom_id_2_4c());
               std::string id_3 = util::remove_whitespace(torsion_restraint[i].atom_id_3_4c());
               std::string id_4 = util::remove_whitespace(torsion_restraint[i].atom_id_4_4c());
               std::string qan_1 = quoted_atom_name(id_1);
               std::string qan_2 = quoted_atom_name(id_2);
               std::string qan_3 = quoted_atom_name(id_3);
               std::string qan_4 = quoted_atom_name(id_4);
               std::string residue_id = residue_info.comp_id;
               mmCIFLoop->PutString(residue_id.c_str(), "comp_id", i);
               // ss = torsion_restraint[i].id().c_str();
               std::string torsion_id = torsion_restraint[i].id();
               mmCIFLoop->PutString(torsion_id.c_str(), "id", i);
               // ss = id_1.c_str();
               std::string torsion_atom_id_1 = id_1;
               mmCIFLoop->PutString(torsion_atom_id_1.c_str(), "atom_id_1", i);
               // ss = id_2.c_str();
               std::string torsion_atom_id_2 = id_2;
               mmCIFLoop->PutString(torsion_atom_id_2.c_str(), "atom_id_2", i);
               // ss = id_3.c_str();
               std::string torsion_atom_id_3 = id_3;
               mmCIFLoop->PutString(torsion_atom_id_3.c_str(), "atom_id_3", i);
               // ss = id_4.c_str();
               std::string torsion_atom_id_4 = id_4;
               mmCIFLoop->PutString(torsion_atom_id_4.c_str(), "atom_id_4", i);
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
               // const char *ss = residue_info.comp_id.c_str();
               std::string id_c = util::remove_whitespace(chiral_restraint[i].atom_id_c_4c());
               std::string id_1 = util::remove_whitespace(chiral_restraint[i].atom_id_1_4c());
               std::string id_2 = util::remove_whitespace(chiral_restraint[i].atom_id_2_4c());
               std::string id_3 = util::remove_whitespace(chiral_restraint[i].atom_id_3_4c());
               std::string qan_c = quoted_atom_name(id_c);
               std::string qan_1 = quoted_atom_name(id_1);
               std::string qan_2 = quoted_atom_name(id_2);
               std::string qan_3 = quoted_atom_name(id_3);
               std::string residue_id = residue_info.comp_id;
               mmCIFLoop->PutString(residue_id.c_str(), "comp_id", i);
               // ss = chiral_restraint[i].Chiral_Id().c_str();
               std::string chiral_id = chiral_restraint[i].Chiral_Id();
               mmCIFLoop->PutString(chiral_id.c_str(), "id", i);
               // ss = id_c.c_str();
               mmCIFLoop->PutString(id_c.c_str(), "atom_id_centre", i);
               mmCIFLoop->PutString(id_1.c_str(), "atom_id_1", i);
               mmCIFLoop->PutString(id_2.c_str(), "atom_id_2", i);
               mmCIFLoop->PutString(id_3.c_str(), "atom_id_3", i);
               int sign = chiral_restraint[i].volume_sign;
               std::string sign_string = "both";
               if (sign == 1)
                  sign_string = "positiv";
               if (sign == -1)
                  sign_string = "negativ";
               mmCIFLoop->PutString(sign_string.c_str(), "volume_sign", i);
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
      if (status == 0) {
              // std::cout << "INFO:: wrote mmCIF \"" << filename << "\"" << std::endl;
              logger.log(log_t::INFO, "wrote mmCIF", filename);
      } else {
              // std::cout << "INFO:: on write mmCIF \"" << filename << "\" status: "
              //                << status << std::endl;
              logger.log(log_t::INFO, "on write mmCIF", filename, "status:", status);
      }
   }
   delete mmCIFFile; // deletes all its attributes too.
}

void
coot::dictionary_residue_restraints_t::write_cif_pdbx_chem_comp_descriptor(mmdb::mmcif::Data *mmCIFData) const {
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
