/*
 * api/coot-molecule.cc
 * 
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */


#include <iostream>
#include <sstream>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-coord-extras.hh"
#include "coords/mmdb.hh"
#include "coot-molecule.hh"
#include "ideal/pepflip.hh"
#include "rama-plot-phi-psi.hh"

#include "pli/sdf-interface-for-export.hh"
#include "ligand/side-chain-densities.hh"

#include "add-terminal-residue.hh"
#include "molecules-container.hh"

bool
coot::molecule_t::is_valid_model_molecule() const {

   bool status = false;
   if (atom_sel.mol)
      status = true;
   return status;

}

int
coot::molecule_t::close_yourself() {

   int status = 0;
   if (is_closed_flag)
      return status;

   if (is_valid_model_molecule()) {
      atom_sel.clear_up();
      status = 1;
   }
   if (is_valid_map_molecule()) {
      clipper::Xmap<float> xmap_empty;
      std::swap(xmap, xmap_empty);
      status = 1;
   }
   is_closed_flag = true;
   return status;
}

bool
coot::molecule_t::is_valid_map_molecule() const {

   bool status = false;
   if (! xmap.is_null()) {
      status = true;
   }
   return status;
}


mmdb::Atom *
coot::molecule_t::cid_to_atom(const std::string &cid) const {

   mmdb::Atom *atom_p = 0;
   if (atom_sel.mol) {
      int selHnd = atom_sel.mol->NewSelection(); // d
      mmdb::Atom **SelAtoms = nullptr;
      int nSelAtoms = 0;
      atom_sel.mol->Select(selHnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selHnd, SelAtoms, nSelAtoms);
      if (nSelAtoms > 0) {
         atom_p = SelAtoms[0];
      }
      atom_sel.mol->DeleteSelection(selHnd);
   }
   return atom_p;
}

mmdb::Residue *
coot::molecule_t::cid_to_residue(const std::string &cid) const {

   mmdb::Residue *residue_p = 0;
   if (atom_sel.mol) {
      int selHnd = atom_sel.mol->NewSelection(); // d
      mmdb::Residue **SelResidues;
      int nSelResidues = 0;
      atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, cid.c_str(), mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
      if (nSelResidues > 0) {
         residue_p = SelResidues[0];
      }
      atom_sel.mol->DeleteSelection(selHnd);
   }
   return residue_p;
}

std::pair<bool, coot::residue_spec_t>
coot::molecule_t::cid_to_residue_spec(const std::string &cid) const {

   bool status = false;
   coot::residue_spec_t rs;
   if (atom_sel.mol) {
      int selHnd = atom_sel.mol->NewSelection(); // d
      mmdb::Residue **SelResidues;
      int nSelResidues = 0;
      atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, cid.c_str(), mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
      if (nSelResidues > 0) {
         mmdb::Residue *residue_p = SelResidues[0];
         coot::residue_spec_t rs_inner(residue_p);
         rs = rs_inner;
         status = true;
      }
      atom_sel.mol->DeleteSelection(selHnd);
   }
   return std::make_pair(status, rs);
}

std::pair<bool, coot::atom_spec_t>
coot::molecule_t::cid_to_atom_spec(const std::string &cid) const {

   bool status = false;
   coot::atom_spec_t atom_spec;
   if (atom_sel.mol) {
      int selHnd = atom_sel.mol->NewSelection(); // d
      mmdb::Atom **SelAtoms;
      int nSelAtoms = 0;
      atom_sel.mol->Select(selHnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selHnd, SelAtoms, nSelAtoms);
      if (nSelAtoms > 0) {
         mmdb::Atom *atom_p = SelAtoms[0];
         coot::atom_spec_t atom_spec_inner(atom_p);
         atom_spec = atom_spec_inner;
         status = true;
      }
      atom_sel.mol->DeleteSelection(selHnd);
   }
   return std::make_pair(status, atom_spec);
}

std::vector<mmdb::Residue *>
coot::molecule_t::cid_to_residues(const std::string &atom_selection_cids) const {

   std::vector<mmdb::Residue *> v;
   if (! atom_sel.mol) return v;

   std::set<mmdb::Residue *> residue_set;
   std::vector<std::string> cid_v = coot::util::split_string(atom_selection_cids, "||");
   int selHnd = atom_sel.mol->NewSelection(); // d
   if (! cid_v.empty()) {
      int nSelResidues = 0;
      mmdb::Residue **SelResidues = 0;
      atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
      for (const auto &cid : cid_v) {
         atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, cid.c_str(), mmdb::SKEY_OR);
         atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
         for (int i=0; i<nSelResidues; i++) {
            residue_set.insert(SelResidues[i]);
         }
      }
   }

   std::set<mmdb::Residue *>::const_iterator it;
   v.reserve(residue_set.size()); // micro-optimization
   for (it=residue_set.begin(); it!=residue_set.end(); ++it)
      v.push_back(*it);

   return v;
}

// can return null
mmdb::Residue *
coot::molecule_t::get_residue(const std::string &residue_cid) const {

   mmdb::Residue *residue_p = cid_to_residue(residue_cid);
   return residue_p;
}




// restore from (previous) backup
void
coot::molecule_t::restore_from_backup(int mod_index, const std::string &cwd) {

#if 0

   auto print_history = [] (const std::vector<std::string> &history_filename_vec) {
      for (unsigned int i=0; i<history_filename_vec.size(); i++) {
         std::cout << "  " << i << " " << history_filename_vec[i] << std::endl;
      }
   };

   if (true)
      std::cout << "debug:: restore_from_backup() requested mod_index: " << mod_index
                << " history size: " << history_filename_vec.size() << std::endl;

   if (history_filename_vec.empty()) {
      std::cout << "ERROR:: in restore_from_backup(): empty history_filename_vec " << std::endl;
      return;
   }

   if (mod_index >= int(history_filename_vec.size())) {
      std::cout << "ERROR:: in restore_from_backup(): bad mod_index " << mod_index << std::endl;
      print_history(history_filename_vec);
      return;
   }
   if (mod_index < 0) {
      std::cout << "ERROR:: in restore_from_backup(): bad mod_index " << mod_index << std::endl;
      print_history(history_filename_vec);
      return;
   }

   std::string file_name = history_filename_vec[mod_index];
   // hostage to forture here?
   // bool use_gemmi = false; // 20240112-PE now use the class data item (which is usually true)
   atom_selection_container_t asc = get_atom_selection(file_name, use_gemmi);
   if (asc.read_success) {
      save_info.set_modification_index(mod_index);
      atom_sel.clear_up();
      atom_sel = asc;
      // 20221018-PE no bond regeneration, maybe there should be?
   }

#endif

}


void
coot::molecule_t::replace_molecule_by_model_from_file(const std::string &pdb_file_name) {

   // bool use_gemmi = false; // 20240112-PE now use the class data item (which is usually true)
   atom_selection_container_t asc = get_atom_selection(pdb_file_name, use_gemmi, true, false);
   if (asc.read_success) {
      atom_sel.clear_up();
      atom_sel = asc;
   }
}


int
coot::molecule_t::write_coordinates(const std::string &file_name) const {

   int err = 1;
   if (atom_sel.n_selected_atoms > 0) {
      std::string ext = coot::util::file_name_extension(file_name);
      if (coot::util::extension_is_for_shelx_coords(ext)) {
         write_shelx_ins_file(file_name);
      } else {
         if (ext == ".cif") {
            mmdb::byte bz = mmdb::io::GZM_NONE; // 20221018-PE  this should be used
            err = coot::write_coords_cif(atom_sel.mol, file_name);
         } else {
            mmdb::byte bz = mmdb::io::GZM_NONE; // 20221018-PE  this should be used too
            err = coot::write_coords_pdb(atom_sel.mol, file_name);
         }
      }
   }
   return err; // the return value of WritePDBASCII() or WriteCIFASCII(). mmdb return type
}

unsigned int
coot::molecule_t::get_number_of_atoms() const {

   unsigned int n = 0;
   // for(int imod = 1; imod<=asc.mol->GetNumberOfModels(); imod++) {
   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
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
                     n++;
                  }
               }
            }
         }
      }
   }
   return n;
}

int
coot::molecule_t::get_number_of_hydrogen_atoms() const {

   int n = 0;
   // for(int imod = 1; imod<=asc.mol->GetNumberOfModels(); imod++) {
   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
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
		  std::string ele(at->element);
		  if (ele == " H") {
		     if (! at->isTer()) {
			n++;
		     }
                  }
               }
            }
         }
      }
   }
   return n;
}

#include "coot-utils/atom-selection-container.hh"

float
coot::molecule_t::get_molecule_diameter() const {

   // sample atom pairs

   float f = -1;
   if (atom_sel.mol) {
      f = coot::get_molecule_diameter(atom_sel);
   }
   return f;
}

//! Get Radius of Gyration
//!
//! @param imol is the model molecule index
//!
//! @return the molecule centre. If the number is less than zero, there
//! was a problem finding the molecule or atoms.
double
coot::molecule_t::get_radius_of_gyration() const {

   double d = -1.0;
   if (is_valid_model_molecule()) {
      std::pair<bool, double> rgp = coot::radius_of_gyration(atom_sel.mol);
      if (rgp.first) {
	 d = rgp.second;
      }
   }
   return d;
}



std::string
coot::molecule_t::name_for_display_manager() const {

   std::string s("");

   bool show_paths_in_display_manager_flag = false;
   if (show_paths_in_display_manager_flag) {
      s = name;
   } else {
      if (is_valid_model_molecule()) {
         std::string::size_type islash = name.find_last_of("/");
         if (islash == std::string::npos) {
            s = name;
         } else {
            s = name.substr(islash+1, name.length());
         }
      } else {
         // This is a map, so we want to strip of xxx/ from each of
         // the (space separated) strings.
         // e.g.:
         // thing/other.mtz -> other.mtz
         // but
         // Averaged -> Averaged

         std::vector<std::string> v = coot::util::split_string(name, " ");
         for (unsigned int i=0; i<v.size(); i++) {
            if (i > 0)
               s += " ";
            std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(v[i]);
            if (p.second == "")
               s += v[i];
            else
               s += p.second;
         }
      }
   }
   return s;
}

std::string
coot::molecule_t::dotted_chopped_name() const {

   int go_to_atom_menu_label_n_chars_max = 80;
   std::string ss = coot::util::int_to_string(imol_no);
   ss += " " ;
   int ilen = name.length();
   int left_size = ilen-go_to_atom_menu_label_n_chars_max;
   if (left_size <= 0) {
      // no chop
      left_size = 0;
   } else {
      // chop
      ss += "...";
   }
   ss += name.substr(left_size, ilen);
   return ss;
}


#if 0 // old
// backups:

// Backup filename: return a stub.
//
std::string
coot::molecule_t::get_save_molecule_filename(const std::string &dir) {

   auto replace_char = [] (const std::string &s, char a) {
                          std::string r = s;
                          int slen = s.length();
                          for (int i=0; i<slen; i++) {
                             if (r[i] == a)
                                r[i] = '_';
                          }
                          return r;
                       };


   bool decolonify = true;
   bool backup_compress_files_flag = false;
   bool unpathed_backup_file_names_flag = false;
   std::string t_name_1 = name;

   if (unpathed_backup_file_names_flag)
      t_name_1 = name_for_display_manager();
   std::string t_name_2 = replace_char(t_name_1, '/');
   std::string t_name_3 = replace_char(t_name_2, ' ');

   if (save_time_string.empty()) {
      time_t t;
      time(&t);
      char *chars_time = ctime(&t);
      int l = strlen(chars_time);
      save_time_string = chars_time;
      if (! save_time_string.empty()) {
         std::string::size_type l = save_time_string.length();
         save_time_string = save_time_string.substr(0, l-1);
      }
      save_time_string = replace_char(save_time_string, ' ');
      save_time_string = replace_char(save_time_string, '/');
      if (decolonify)
         save_time_string = replace_char(save_time_string, ':');
   }
   std::string time_string = save_time_string;
   std::string t_name_4 = t_name_3 + "_" + time_string;

   //    std::string index_string = coot::util::int_to_string(history_index);
   std::string index_string = modification_info.index_string();
   std::string t_name_5 = t_name_4 + "_modification_" + index_string;

   std::string extension = ".pdb";
   if (coot::is_mmcif_filename(name))
      extension = ".cif";
   if (is_from_shelx_ins_flag)
      extension = ".res";
   if (backup_compress_files_flag)
      extension += ".gz";

   std::string t_name_6 = t_name_5 + extension;

   std::string save_file_name = coot::util::append_dir_file(dir, t_name_6);
   return save_file_name;

}

#endif

#if 0 // 20240221-PE this had a space reminder for me to look at it before commit
void
coot::molecule_t::save_history_file_name(const std::string &file) {

   // 20221016-PE fold this function into save_info

   // First, history_index is zero and the vec is zero,
   // normal service, then another backup: history_index is 1 and vec is 1.
   //
   int history_filename_vec_size = history_filename_vec.size();
   if (modification_info.modification_index == history_filename_vec_size) {
      history_filename_vec.push_back(file);
      std::cout << "debug:: in save_history_file_name() added to history vec" << file << " " << history_filename_vec.size() << std::endl;
   } else {
      // we have gone back in history.
      //
      if (modification_info.modification_index < history_filename_vec_size) {
         history_filename_vec[modification_info.modification_index] = file;
      }
   }
}

#endif

void
coot::molecule_t::transform_by(mmdb::mat44 mat) {

   bool verbose = false;

   if (is_valid_model_molecule()) {
      clipper::Coord_orth co;
      clipper::Coord_orth trans_pos;
      make_backup("transform_by");
      clipper::Mat33<double> clipper_mat(mat[0][0], mat[0][1], mat[0][2],
                                         mat[1][0], mat[1][1], mat[1][2],
                                         mat[2][0], mat[2][1], mat[2][2]);
      clipper::Coord_orth cco(mat[0][3], mat[1][3], mat[2][3]);
      clipper::RTop_orth rtop(clipper_mat, cco);
      if (verbose)
         std::cout << "INFO:: coordinates transformed by orthogonal matrix: \n"
                   << rtop.format() << std::endl;
      clipper::Rotation rtn( clipper_mat );
      clipper::Polar_ccp4 polar = rtn.polar_ccp4();
      clipper::Euler_ccp4 euler = rtn.euler_ccp4();
      if (verbose) {
         std::cout << "  Rotation - polar (omega,phi,kappa)  " << clipper::Util::rad2d(polar.omega()) << " " << clipper::Util::rad2d(polar.phi())  << " " << clipper::Util::rad2d(polar.kappa()) << std::endl;
         std::cout << "  Rotation - euler (alpha,beta,gamma) " << clipper::Util::rad2d(euler.alpha()) << " " << clipper::Util::rad2d(euler.beta()) << " " << clipper::Util::rad2d(euler.gamma()) << std::endl;
         std::cout << "  Translation - Angstroms             " << cco.x() << " " << cco.y() << " " << cco.z() << " " << std::endl;
      }
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
         mmdb::Atom *at = atom_sel.atom_selection[i];
         co = clipper::Coord_orth(at->x, at->y, at->z);
         trans_pos = co.transform(rtop);
         at->x = trans_pos.x();
         at->y = trans_pos.y();
         at->z = trans_pos.z();
         if (false) // debugging
            if (co.x() < 0.0 && co.x() > -10.0)
               if (co.y() < 20.0 && co.y() > 10.0)
                  if (co.z() < 30.0 && co.z() > 20.0)
                     std::cout << i << " from " << co.format() << " to " << at->x << " " << at->y << " " << at->z << std::endl;
      }
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
   }
}



// shelx stuff
//
std::pair<int, std::string>
coot::molecule_t::write_shelx_ins_file(const std::string &filename) const {

   // std::cout << "DEBUG:: starting write_shelx_ins_file in molecule "<< std::endl;
   // shelxins.debug();

   std::pair<int, std::string> p(1, "");

   if (atom_sel.n_selected_atoms > 0) {
      // 20221018-PE  restore this when write_ins_file() is const.
      //              This function needs to be const because write_coordinates() is const.
      // p = shelxins.write_ins_file(atom_sel.mol, filename, is_from_shelx_ins_flag);
   } else {
      p.second = "WARNING:: No atoms to write!";
   }
   return p;
}


std::vector<mmdb::Residue *>
coot::molecule_t::select_residues(const residue_spec_t &residue_spec, const std::string &mode_in) const {

   // why is this not in utils? Make it so.

   auto all_residues = [] (mmdb::Manager *mol) {

      std::vector<mmdb::Residue *> rv;
      if (!mol) return rv;
      int n_models = mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            mmdb::Chain *chain_p;
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               chain_p = model_p->GetChain(ichain);
               int nres = chain_p->GetNumberOfResidues();
               mmdb::Residue *residue_p;
               for (int ires=0; ires<nres; ires++) {
                  residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     rv.push_back(residue_p);
                  }
               }
            }
         }
      }
      return rv;
   };

   std::string mode = mode_in;
   if (mode == "LITERAL") mode = "SINGLE";
   std::vector<mmdb::Residue *> rv;
   mmdb::Manager *mol = atom_sel.mol;
   mmdb::Residue *residue_p = coot::util::get_residue(residue_spec, mol);
   if (residue_p) {
      if (mode == "SINGLE") {
         rv.push_back(residue_p);
      }
      if (mode == "TRIPLE") {
         mmdb::Residue *r_p_1 = coot::util::get_following_residue(residue_spec, mol);
         mmdb::Residue *r_m_1 = coot::util::get_previous_residue(residue_spec, mol);
         if (r_m_1) rv.push_back(r_m_1);
         if (true ) rv.push_back(residue_p);
         if (r_p_1) rv.push_back(r_p_1);
      }
      if (mode == "QUINTUPLE") {
         mmdb::Residue *r_p_1 = coot::util::get_following_residue(residue_spec, mol);
         mmdb::Residue *r_p_2 = coot::util::get_following_residue(coot::residue_spec_t(r_p_1), mol);
         mmdb::Residue *r_m_1 = coot::util::get_previous_residue(residue_spec, mol);
         mmdb::Residue *r_m_2 = coot::util::get_previous_residue(coot::residue_spec_t(r_m_1), mol);
         if (r_m_2) rv.push_back(r_m_2);
         if (r_m_1) rv.push_back(r_m_1);
         if (true ) rv.push_back(residue_p);
         if (r_p_1) rv.push_back(r_p_1);
         if (r_p_2) rv.push_back(r_p_2);
      }
      if (mode == "HEPTUPLE") {
         mmdb::Residue *r_p_1 = coot::util::get_following_residue(residue_spec, mol);
         mmdb::Residue *r_p_2 = coot::util::get_following_residue(coot::residue_spec_t(r_p_1), mol);
         mmdb::Residue *r_p_3 = coot::util::get_following_residue(coot::residue_spec_t(r_p_2), mol);
         mmdb::Residue *r_m_1 = coot::util::get_previous_residue(residue_spec, mol);
         mmdb::Residue *r_m_2 = coot::util::get_previous_residue(coot::residue_spec_t(r_m_1), mol);
         mmdb::Residue *r_m_3 = coot::util::get_previous_residue(coot::residue_spec_t(r_m_2), mol);
         if (r_m_3) rv.push_back(r_m_3);
         if (r_m_2) rv.push_back(r_m_2);
         if (r_m_1) rv.push_back(r_m_1);
         if (true ) rv.push_back(residue_p);
         if (r_p_1) rv.push_back(r_p_1);
         if (r_p_2) rv.push_back(r_p_2);
         if (r_p_3) rv.push_back(r_p_3);
      }
      if (mode == "CHAIN") {
         mmdb::Chain *chain_p = residue_p->GetChain();
         rv = coot::util::residues_in_chain(chain_p);
      }
      if (mode == "ALL") {
         rv = all_residues(mol);
      }
      if (mode == "SPHERE") {
         float radius = 4.2;
         // do these need to be sorted here?
         auto v = coot::residues_near_residue(residue_p, mol, radius);
         rv.push_back(residue_p);
         std::move(v.begin(), v.end(), std::back_inserter(rv));
      }
      if (mode == "BIG_SPHERE") {
         float radius = 8.0;
         // do these need to be sorted here?
         auto v = coot::residues_near_residue(residue_p, mol, radius);
         rv.push_back(residue_p);
         // std::move(v.begin(), v.end(), std::back_inserter(rv));
         for (const auto &r : v)
            rv.push_back(r);
      }
   }

   return rv;
}

//!
std::vector<mmdb::Residue *>
coot::molecule_t::select_residues(const std::string &multi_cids, const std::string &mode) const {

   auto set_to_vec = [] (std::set<mmdb::Residue *> rs) {
      std::vector<mmdb::Residue *> rv;
      std::set<mmdb::Residue *>::const_iterator it;
      for (it=rs.begin(); it!=rs.end(); ++it)
         rv.push_back(*it);
      return rv;
   };

   std::vector<mmdb::Residue *> rv;
   std::set<mmdb::Residue *> rs;

   std::vector<std::string> v = coot::util::split_string(multi_cids, "||");
   if (! v.empty()) {
      for (const auto &cid : v) {
         int selHnd = atom_sel.mol->NewSelection(); // d
         mmdb::Residue **SelResidues;
         int nSelResidues = 0;
         atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, cid.c_str(), mmdb::SKEY_OR);
         atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
         for (int i=0; i<nSelResidues; i++) {
            mmdb::Residue *residue_p = SelResidues[i];
            coot::residue_spec_t residue_spec(residue_p);
            std::vector<mmdb::Residue *> neighbs = select_residues(residue_spec, mode);
            for (unsigned int j=0; j<neighbs.size(); j++) {
               rs.insert(neighbs[j]);
            }
         }
         atom_sel.mol->DeleteSelection(selHnd);
      }
   }

   rv = set_to_vec(rs);
   return rv;
}


//! resno_start and resno_end are inclusive
std::vector<mmdb::Residue *>
coot::molecule_t::select_residues(const std::string &chain_id, int resno_start, int resno_end) const {

   std::vector<mmdb::Residue *> rv;

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         std::string chain_id_this(chain_p->GetChainID());
         if (chain_id_this == chain_id) {
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int res_no_this = residue_p->GetSeqNum();
                  if (res_no_this >= resno_start) {
                     if (res_no_this <= resno_end) {
                        rv.push_back(residue_p);
                     }
                  }
               }
            }
         }
      }
   }
   return rv;
}




bool
coot::molecule_t::moving_atom_matches(mmdb::Atom *at, int this_mol_index_maybe) const {

   bool matches = false;
   if (atom_sel.n_selected_atoms > 0) {
      if (this_mol_index_maybe >= atom_sel.n_selected_atoms) {
         return false;
      } else {
         std::string atom_name_mov = at->name;
         std::string ins_code_mov  = at->GetInsCode();
         std::string alt_conf_mov  = at->altLoc;
         std::string chain_id_mov  = at->GetChainID();
         int resno_mov = at->GetSeqNum();

         std::string atom_name_ref = atom_sel.atom_selection[this_mol_index_maybe]->name;
         std::string ins_code_ref  = atom_sel.atom_selection[this_mol_index_maybe]->GetInsCode();
         std::string alt_conf_ref  = atom_sel.atom_selection[this_mol_index_maybe]->altLoc;
         std::string chain_id_ref  = atom_sel.atom_selection[this_mol_index_maybe]->GetChainID();
         int resno_ref = atom_sel.atom_selection[this_mol_index_maybe]->GetSeqNum();

         if (atom_name_ref == atom_name_mov) {
            if (ins_code_ref == ins_code_mov) {
               if (resno_ref == resno_mov) {
                  if (alt_conf_ref == alt_conf_mov) {
                     if (chain_id_mov == chain_id_ref) { // 20170612 extra condition added, Oliver Clarke bug
                        matches = true;
                     }
                  }
               }
            }
         }
      }
   }
   return matches;
}


int
coot::molecule_t::full_atom_spec_to_atom_index(const coot::atom_spec_t &atom_spec) const {

   return full_atom_spec_to_atom_index(atom_spec.chain_id,
                                       atom_spec.res_no,
                                       atom_spec.ins_code,
                                       atom_spec.atom_name,
                                       atom_spec.alt_conf);

}

// return -1 on no atom found.
int
coot::molecule_t::full_atom_spec_to_atom_index(const std::string &chain,
                                               int resno,
                                               const std::string &insertion_code,
                                               const std::string &atom_name,
                                               const std::string &alt_conf) const {

   int iatom_index = -1;

   // some protection for null molecule.
   if (! atom_sel.mol) {
      std::cout << "ERROR:: null molecule for molecule number "
                << imol_no << " pointer: " << atom_sel.mol
                << " (in full_atom_spec_to_atom_index)" << std::endl;
      return -1;
   }

   int selHnd = atom_sel.mol->NewSelection(); // d

   atom_sel.mol->SelectAtoms(selHnd, 0, chain.c_str(),
                            resno, insertion_code.c_str(), // start, insertion code
                            resno, insertion_code.c_str(), // end, insertion code
                            "*", // residue name
                            atom_name.c_str(),
                            "*", // elements
                            alt_conf.c_str()); // alt locs

   int nSelAtoms;
   mmdb::PPAtom local_SelAtom;
   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

   if (true)
      std::cout << "DEBUG:: full_atom_spec_to_atom_index() for :" << chain << ": "
                << resno << " :" << insertion_code << ": :"
                << atom_name << ": :" << alt_conf << ": finds " << nSelAtoms <<  " atoms\n";

   if (nSelAtoms == 0) {

      std::cout << "WARNING:: full_atom_spec_to_atom_index() Could not find "
                << "\"" << atom_name << "\"," << "\"" << alt_conf  << "\"" << "/"
                << resno << insertion_code << "/" << chain << " in this molecule: ("
                <<  imol_no << ") " << name << std::endl;

      int selHnd2 = atom_sel.mol->NewSelection(); // d

      atom_sel.mol->SelectAtoms(selHnd2, 0,
                                chain.c_str(),
                                resno, "*", // start, insertion code
                                resno, "*", // end, insertion code
                                "*", // residue name
                                "*", // atom name
                                "*", // elements
                                "*"); // alt locs

      atom_sel.mol->GetSelIndex(selHnd2, local_SelAtom, nSelAtoms);

      if (false) { // debugging.
         std::cout << "There were " << nSelAtoms << " atoms in that residue:\n";
         std::cout << "debug:: full_atom_spec_to_atom_index() resno " << resno
                   << " (cf MinInt4) " << mmdb::MinInt4 << "\n";
         if (resno == mmdb::MinInt4) {
            std::cout << "      residue with resno MinInt4\n";
         } else {
            for (int i=0; i<nSelAtoms; i++) {
               std::cout << "      " << local_SelAtom[i] << "\n";
            }
         }
      }

      atom_sel.mol->DeleteSelection(selHnd2);

   } else {

      int idx = 0;
      if (nSelAtoms != 1) {
         // the wildcard atom selection case "*HO2"
         for (int i=0; i<nSelAtoms; i++) {
            if (std::string(local_SelAtom[i]->GetChainID()) == chain) {
               if (local_SelAtom[i]->residue->seqNum == resno) {
                  if (std::string(local_SelAtom[i]->GetInsCode()) == insertion_code) {
                     if (std::string(local_SelAtom[i]->name) == atom_name) {
                        if (std::string(local_SelAtom[i]->altLoc) == alt_conf) {
                           idx = i;
                           break;
                        }
                     }
                  }
               }
            }
         }
      }

      int iatom_index_udd = -1;
      int ic;
      if (local_SelAtom[idx]->GetUDData(atom_sel.UDDAtomIndexHandle, ic) == mmdb::UDDATA_Ok) {
         iatom_index_udd = ic;
      }
      iatom_index = iatom_index_udd;
   }
   atom_sel.mol->DeleteSelection(selHnd);

   return iatom_index;
}




// helper function for above function
bool
coot::molecule_t::movable_atom(mmdb::Atom *mol_atom, bool replace_coords_with_zero_occ_flag) const {

   // std::cout << "debug:: movable_atom() called with atom " << mol_atom << std::endl;

   if (! mol_atom) {
      std::cout << "ERROR:: null mol_atom in movable_atom()" << std::endl;
      return false;
   }

   bool m = true;

   if ((mol_atom->occupancy < 0.0001) &&
       (mol_atom->occupancy > -0.0001))
      if (replace_coords_with_zero_occ_flag == 0)
         m = 0; // zero occupancy and "dont move zero occ atoms is set"
   return m;
}


// We just added a new atom to a residue, now we need to adjust the
// occupancy of the other atoms (so that we don't get residues with
// atoms whose occupancy is greater than 1.0 (Care for SHELX molecule?)).
// at doesn't have to be in residue.
//
// Perhaps this can be a utility function?
//
void
coot::molecule_t::adjust_occupancy_other_residue_atoms(mmdb::Atom *at,
                                                       mmdb::Residue *residue,
                                                       short int force_sum_1_flag) {

   if (!is_from_shelx_ins_flag) {
      int nResidueAtoms;
      mmdb::PPAtom ResidueAtoms = 0;
      residue->GetAtomTable(ResidueAtoms, nResidueAtoms);
      float new_atom_occ = at->occupancy;
      std::string new_atom_name(at->name);
      std::string new_atom_altconf(at->altLoc);
      std::vector<mmdb::Atom *> same_name_atoms;
      float sum_occ = 0;
      for (int i=0; i<nResidueAtoms; i++) {
         std::string this_atom_name(ResidueAtoms[i]->name);
         std::string this_atom_altloc(ResidueAtoms[i]->altLoc);
         if (this_atom_name == new_atom_name) {
            if (this_atom_altloc != new_atom_altconf) {
               same_name_atoms.push_back(ResidueAtoms[i]);
               sum_occ += ResidueAtoms[i]->occupancy;
            }
         }
      }

      //
      if (sum_occ > 0.01) {
         if (same_name_atoms.size() > 0) {
            float other_atom_occ_sum = 0.0;
            for (unsigned int i=0; i<same_name_atoms.size(); i++)
               other_atom_occ_sum += same_name_atoms[i]->occupancy;

            float remainder = 1.0 - new_atom_occ;
            float f = remainder/other_atom_occ_sum;
            for (unsigned int i=0; i<same_name_atoms.size(); i++) {
               if (0) // debug
                  std::cout << "debug " << same_name_atoms[i]
                            << " mulitplying occ " << same_name_atoms[i]->occupancy
                            << " by " << remainder << "/" << other_atom_occ_sum << "\n";
               same_name_atoms[i]->occupancy *= f;
            }
         }
      }
   }
}


// Put the regularization results back into the molecule:
//
//// Recall that regularized_asc contains an atom_selection_container_t
// with the new coordinates in.  the mol contains all the molecule and
// the atom_selection contains just the moving parts.
//
// This function no longer does a backup or updates the save_info!
// The calling function should do that.
//
void
coot::molecule_t::replace_coords(const atom_selection_container_t &asc,
                                 bool change_altconf_occs_flag,
                                 bool replace_coords_with_zero_occ_flag) {

   // 20221122-PE convert this to fast indexing one day!

   int n_atom = 0;
   int tmp_index;
   bool debug = false;
   float add_alt_conf_new_atoms_occupancy = 0.5; // was a static in graphics_info_t

   // make_backup(); // Is replace_coords() the right place for make_backup()?
                     // Perhaps the calling function should make the backup?
                     // Let's presume so.

   if (true) {
      std::cout << "DEBUG:: --------------- replace_coords replacing "
                << asc.n_selected_atoms << " atoms " << std::endl;
      for (int i=0; i<asc.n_selected_atoms; i++) {
         mmdb::Atom *atom = asc.atom_selection[i];
         bool is_ter_state = atom->isTer();
         std::cout << "DEBUG:: in replace_coords, intermediate atom: " << i << " " << atom << " "
                   << "chain-id: "
                   << atom->residue->GetChainID() <<  ": "
                   << atom->residue->seqNum << " inscode \""
                   << atom->GetInsCode() << "\" name \""
                   << atom->name << "\" altloc \""
                   << atom->altLoc << "\" occupancy: "
                   << atom->occupancy << " :"
                   << " TER state: " << is_ter_state << std::endl;
      }
   }

   // For each atom in the new set of atoms:
   //
   for (int i=0; i<asc.n_selected_atoms; i++) {
      int idx = -1;
      mmdb::Atom *atom = asc.atom_selection[i];

      // std::cout << "replace_coords(): atom at " << i << " " << atom << " " << coot::atom_spec_t(atom) << std::endl;
      if (! atom->isTer()) {

         if (debug) { // debug
            std::cout << "considering replacement for selected atom " << coot::atom_spec_t(atom) << std::endl;

            //
            // idx = atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
            // atom->residue->seqNum, std::string(atom->name));

         }

         // std::cout << "------------------ replace_coords() with UDDOldAtomIndexHandle() " << asc.UDDOldAtomIndexHandle << std::endl;

         if (asc.UDDOldAtomIndexHandle >= 0) { // OK for fast atom indexing

            std::cout << "------------ replace_coords() path A" << std::endl;

            if (debug)
               std::cout << "... OK for fast atom indexing, asc.UDDOldAtomIndexHandle: "
                         << asc.UDDOldAtomIndexHandle
                         << " for atom " << coot::atom_spec_t(atom) << std::endl;

            if (atom->GetUDData(asc.UDDOldAtomIndexHandle, tmp_index) == mmdb::UDDATA_Ok) {

               if (debug)
                  std::cout << "OK, good GetUDData() for atom " << coot::atom_spec_t(atom) << std::endl;
               if (tmp_index >= 0) {
                  if (moving_atom_matches(atom, tmp_index)) {
                     // std::cout << "      DEBUG:: successfully found old atom index" << std::endl;
                     idx = tmp_index;
                  } else {
                     // std::cout << "DEBUG:: atom index mismatch" << std::endl;
                     idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
                                                        atom->residue->seqNum,
                                                        std::string(atom->GetInsCode()),
                                                        std::string(atom->name),
                                                        std::string(atom->altLoc));
                     // std::cout << "DEBUG:: full_atom_spec_to_atom_index gives index: " << idx << std::endl;
                  }
               } else {
                  // This shouldn't happen.
                  std::cout << "Good Handle, bad index found for old atom: specing" << std::endl;
                  idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
                                                     atom->residue->seqNum,
                                                     std::string(atom->GetInsCode()),
                                                     std::string(atom->name),
                                                     std::string(atom->altLoc));
               }
            } else {

               std::cout << "ERROR:: non-bad handle (" << asc.UDDOldAtomIndexHandle
                         <<  "), but bad GetUDData() for atom " << coot::atom_spec_t(atom) << std::endl;

            }
         } else {

            if (true)
               std::cout << "DEBUG:: asc.UDDOldAtomIndexHandle is "
                         << asc.UDDOldAtomIndexHandle << " using full atom spec to atom index..."
                         << std::endl;

            idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
                                               atom->residue->seqNum,
                                               std::string(atom->GetInsCode()),
                                               std::string(atom->name),
                                               std::string(atom->altLoc));

            std::cout << "full_atom_spec_to_atom_index() returned " << idx << " for " << coot::atom_spec_t(atom) << std::endl;
            if (idx != -1) {
               mmdb::Atom *mol_atom = atom_sel.atom_selection[idx];
               std::cout << "mol_atom " << coot::atom_spec_t(mol_atom) << std::endl;
            }

            if (idx == -1) {
               std::cout << "DEBUG:: idx: " << idx << "\n";
               std::cout << "ERROR:: failed to find atom in molecule: chain-id :"
                         << std::string(atom->residue->GetChainID()) <<  ": res_no "
                         << atom->residue->seqNum << " inscode :"
                         << std::string(atom->GetInsCode()) << ": name :"
                         << std::string(atom->name) << ": altloc :"
                         << std::string(atom->altLoc) << ":" << std::endl;
            }
         }

         // std::cout << "----- replace_coords() with change_altconf_occs_flag " << change_altconf_occs_flag << std::endl;

         if (change_altconf_occs_flag) {
            if (idx >= 0) {
               n_atom++;
               mmdb::Atom *mol_atom = atom_sel.atom_selection[idx];
               float atom_occ = atom->occupancy;
               // if this is a shelx molecule, then we don't change
               // occupancies this way.  We do it by changing the FVAR
               if (is_from_shelx_ins_flag) {
                  atom_occ = mol_atom->occupancy;

                  // OK, one more go.  We have an occupancy of 31 or -31
                  // say.  Now, the alt conf atoms has been immediately
                  // added with the old occupancy for the actual FVAR number
                  // - this happens before we get to twiddle the occupancy
                  // slider.  So here we have to find out the index of the
                  // replaced atom and set it's fvar to whatever the slider
                  // value had been set to.

                  int fvar_number = coot::ShelxIns::shelx_occ_to_fvar(atom_occ);
                  if (fvar_number > 1) {
                     //                std::cout << "DEBUG:: replace_coords: setting fvar number "
                     //                          <<  fvar_number << " (generated from occ " << atom_occ << ") to "
                     //                          << graphics_info_t::add_alt_conf_new_atoms_occupancy << std::endl;
                     shelxins.set_fvar(fvar_number, add_alt_conf_new_atoms_occupancy);
                  }

                  if (true) {
                     coot::Cartesian old_pos(mol_atom->x, mol_atom->y, mol_atom->z);
                     coot::Cartesian new_pos(atom->x, atom->y, atom->z);
                     double d = (new_pos - old_pos).amplitude();
                     if (false) {
                        std::cout << "    changing coords for atom with idx " << idx << " "
                                  << coot::atom_spec_t(mol_atom) << std::endl;
                        std::cout << "   " << old_pos << " " << new_pos << " moved-by " << d << std::endl;
                     }
                  }

                  if (movable_atom(mol_atom, replace_coords_with_zero_occ_flag))
                     mol_atom->SetCoordinates(atom->x,
                                              atom->y,
                                              atom->z,
                                              atom_occ,
                                              mol_atom->tempFactor);
               } else {
                  if (movable_atom(mol_atom, replace_coords_with_zero_occ_flag))
                     mol_atom->SetCoordinates(atom->x,
                                              atom->y,
                                              atom->z,
                                              atom_occ,
                                              mol_atom->tempFactor);
               }

               // similarly we adjust occupancy if this is not a shelx molecule
               if (! is_from_shelx_ins_flag) {
                  adjust_occupancy_other_residue_atoms(mol_atom, mol_atom->residue, 0);
               }
               // std::cout << atom << " coords replace " << idx << " " << mol_atom << std::endl;
            } else {
               std::cout << "ERROR:: bad atom index in replace_coords replacing atom: "
                         << atom << std::endl;
            }
         } else {

            // "don't change alt confs" mode

            if (idx != -1 ) {
               if (idx < atom_sel.n_selected_atoms) { // 20240724-PE was <= !
                  mmdb::Atom *mol_atom = atom_sel.atom_selection[idx];
                  if (! mol_atom) {
                     std::cout << "ooops:: mol_atom is null in replace_coords()" << std::endl;
                  }
                  bool is_movable_atom = movable_atom(mol_atom, replace_coords_with_zero_occ_flag);
                  if (is_movable_atom) {
                     if (debug) { // debug
                        coot::Cartesian old_pos(mol_atom->x, mol_atom->y, mol_atom->z);
                        coot::Cartesian new_pos(atom->x, atom->y, atom->z);
                        double d = (new_pos - old_pos).amplitude();
                        if (false) { // debug
                           std::cout << "    changing coords for atom with idx " << idx << " " << coot::atom_spec_t(mol_atom)
                                     << std::endl;
                           std::cout << "   " << old_pos << " " << new_pos << " moved-by " << d << std::endl;
                        }
                     }
                     mol_atom->SetCoordinates(atom->x,
                                              atom->y,
                                              atom->z,
                                              mol_atom->occupancy,
                                              mol_atom->tempFactor);
                     n_atom++;
                  }
               } else {
                  std::cout << "ERROR:: Trapped error! in replace_coords() late block: idx "
                            << idx << " but atom_sel.n_selected_atoms " << atom_sel.n_selected_atoms
                            << std::endl;
               }
            } else {
               std::cout << "WARNING:: bad atom idx -1" << std::endl;
            }
         }
      }
   }
   // std::cout << "INFO:: replace_coords: " << n_atom << " atoms updated." << std::endl;

   // have_unsaved_changes_flag = 1;
   // save_info.new_modification("replace_coords()");

   if (show_symmetry) {  // internal
      update_symmetry();
   }

   // make_bonds_type_checked(__FUNCTION__);

}


int coot::molecule_t::flip_peptide(const coot::atom_spec_t &as_in, const std::string &alt_conf) {

   make_backup("flip_peptide");
   coot::atom_spec_t as = as_in;
   if (as.atom_name == " N  ")
      as.res_no--;
   int result = coot::pepflip(atom_sel.mol, as.chain_id, as.res_no, as.ins_code, alt_conf);
   // save_info.new_modification("flip_peptide");
   return result;

}


std::vector<coot::phi_psi_prob_t>
coot::molecule_t::ramachandran_validation(const ramachandrans_container_t &rc) const {

   auto have_close_peptide_bond = [] (mmdb::Residue *residue_1, mmdb::Residue *residue_2) {

      bool status = false;
      mmdb::Atom **residue_atoms_1 = 0;
      int n_residue_atoms_1 = 0;
      // I should iterate over all alt confs
      residue_1->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
      for (int iat=0; iat<n_residue_atoms_1; iat++) {
         mmdb::Atom *at_1 = residue_atoms_1[iat];
         if (! at_1->isTer()) {
            std::string atom_name_1(at_1->GetAtomName());
            if (atom_name_1 == " C  ") {
               coot::Cartesian pt_c(at_1->x, at_1->y, at_1->z);
               mmdb::Atom **residue_atoms_2 = 0;
               int n_residue_atoms_2 = 0;
               // I should iterate over all alt confs
               residue_2->GetAtomTable(residue_atoms_2, n_residue_atoms_2);
               for (int iat=0; iat<n_residue_atoms_2; iat++) {
                  mmdb::Atom *at_2 = residue_atoms_2[iat];
                  if (! at_2->isTer()) {
                     std::string atom_name_2(at_2->GetAtomName());
                     if (atom_name_2 == " N  ") {
                        coot::Cartesian pt_n(at_1->x, at_1->y, at_1->z);
                        double dd = coot::Cartesian::lengthsq(pt_c, pt_n);
                        double d = std::sqrt(dd);
                        if (d < 3.0) {
                           status = true;
                           break;
                        }
                     }
                  }
               }
            }
         }
      }
      return status;
   };

   auto get_HA_unit_vector = [] (mmdb::Residue *r) {
      bool status = false;
      coot::Cartesian dir;
      mmdb::Atom *CA = r->GetAtom(" CA ");
      mmdb::Atom *C  = r->GetAtom(" C  ");
      mmdb::Atom *N  = r->GetAtom(" N  ");
      mmdb::Atom *CB = r->GetAtom(" CB ");

      if (CA && C && N && CB) {
         coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
         coot::Cartesian  c_pos( C->x,  C->y,  C->z);
         coot::Cartesian  n_pos( N->x,  N->y,  N->z);
         coot::Cartesian cb_pos(CB->x, CB->y, CB->z);
         coot::Cartesian dir_1 = ca_pos - c_pos;
         coot::Cartesian dir_2 = ca_pos - n_pos;
         coot::Cartesian dir_3 = ca_pos - cb_pos;
         coot::Cartesian r = dir_1 + dir_2 + dir_3;
         dir = r.unit();
         status = true;
      } else {
         if (CA && C && N) {
            coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
            coot::Cartesian  c_pos( C->x,  C->y,  C->z);
            coot::Cartesian  n_pos( N->x,  N->y,  N->z);
            coot::Cartesian dir_1 = ca_pos - c_pos;
            coot::Cartesian dir_2 = ca_pos - n_pos;
            coot::Cartesian r = dir_1 + dir_2;
            dir = r.unit();
            status = true;
         }
      }
      return std::make_pair(status, dir);
   };

   std::vector<coot::phi_psi_prob_t> v;

   float rama_ball_pos_offset_scale = 0.6;

   rama_plot::phi_psis_for_model_t ppm(atom_sel.mol);
   // This: std::map<coot::residue_spec_t, phi_psi_t> ppm.phi_psi  is now filled

   std::map<residue_spec_t, rama_plot::phi_psi_t>::const_iterator it;
   for (it=ppm.phi_psi.begin(); it!=ppm.phi_psi.end(); ++it) {
      const auto &phi_psi(it->second);
      mmdb::Residue *rp = phi_psi.residue_prev;
      mmdb::Residue *rt = phi_psi.residue_this;
      mmdb::Residue *rn = phi_psi.residue_next;
      if (rp && rt && rn) {
         if (have_close_peptide_bond(rp, rt)) {
            if (have_close_peptide_bond(rt, rn)) {
               mmdb::Atom *at = rt->GetAtom(" CA "); // 20221006-PE alt-confs another day
               if (at) {
                  coot::Cartesian pos(at->x, at->y, at->z);
                  std::pair<bool, coot::Cartesian> hav = get_HA_unit_vector(rt);
                  coot::Cartesian offset(0,0,rama_ball_pos_offset_scale);
                  if (hav.first) offset = hav.second * rama_ball_pos_offset_scale;
                  coot::util::phi_psi_t cupp(rp, rt, rn);
                  coot::phi_psi_prob_t ppp(cupp, pos + offset, rc);
                  v.push_back(ppp);
               }
            }
         }
      }
   }
   return v;
}


// returns either the specified atom or null if not found
mmdb::Atom *
coot::molecule_t::get_atom(const coot::atom_spec_t &atom_spec) const {

   mmdb::Atom *at = coot::util::get_atom(atom_spec, atom_sel.mol);
   return at;
}

glm::vec4
coot::molecule_t::colour_holder_to_glm(const coot::colour_holder &ch) const {
   return glm::vec4(ch.red, ch.green, ch.blue, ch.alpha);
}

#include "utils/dodec.hh"

// returns either the specified residue or null if not found
mmdb::Residue *
coot::molecule_t::get_residue(const coot::residue_spec_t &residue_spec) const {

   mmdb::Residue *r = coot::util::get_residue(residue_spec, atom_sel.mol);
   return r;

}

std::string
coot::molecule_t::get_residue_name(const residue_spec_t &residue_spec) const {

   std::string n;
   mmdb::Residue *r = coot::util::get_residue(residue_spec, atom_sel.mol);
   if (r) n = r->GetResName();
   return n;
}



std::pair<bool, coot::Cartesian>
coot::molecule_t::get_HA_unit_vector(mmdb::Residue *r) const {

   bool status = false;
   coot::Cartesian dir;
   mmdb::Atom *CA = r->GetAtom(" CA ");
   mmdb::Atom *C  = r->GetAtom(" C  ");
   mmdb::Atom *N  = r->GetAtom(" N  ");
   mmdb::Atom *CB = r->GetAtom(" CB ");

   if (CA && C && N && CB) {
      coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
      coot::Cartesian  c_pos( C->x,  C->y,  C->z);
      coot::Cartesian  n_pos( N->x,  N->y,  N->z);
      coot::Cartesian cb_pos(CB->x, CB->y, CB->z);
      coot::Cartesian dir_1 = ca_pos - c_pos;
      coot::Cartesian dir_2 = ca_pos - n_pos;
      coot::Cartesian dir_3 = ca_pos - cb_pos;
      coot::Cartesian r = dir_1 + dir_2 + dir_3;
      dir = r.unit();
      status = true;
   } else {
      if (CA && C && N) {
         coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
         coot::Cartesian  c_pos( C->x,  C->y,  C->z);
         coot::Cartesian  n_pos( N->x,  N->y,  N->z);
         coot::Cartesian dir_1 = ca_pos - c_pos;
         coot::Cartesian dir_2 = ca_pos - n_pos;
         coot::Cartesian r = dir_1 + dir_2;
         dir = r.unit();
         status = true;
      }
   }
   return std::make_pair(status, dir);
};

coot::simple_mesh_t
coot::molecule_t::get_rotamer_dodecs(coot::protein_geometry *geom_p,
                                     coot::rotamer_probability_tables *rpt) {

   // THis function is an API version of:
   //
   // void
   // Mesh::make_graphical_bonds_rotamer_dodecs(const graphical_bonds_container &gbc,
   // const glm::vec3 &screen_up_dir)

   auto cartesian_to_glm = [] (const coot::Cartesian &c) {
      return glm::vec3(c.x(), c.y(), c.z()); };

   auto clipper_to_glm = [] (const clipper::Coord_orth &c) {
      return glm::vec3(c.x(), c.y(), c.z()); };

   auto clipper_to_cartesian = [] (const clipper::Coord_orth &c) {
      return Cartesian(c.x(), c.y(), c.z()); };

   simple_mesh_t m;

   // use bonds_box

   std::set<int> dummy;
   bool do_rota_markup = true;
   bool change_c_only_flag = true;
   bool goodsell_mode = false;
   bool draw_hydrogen_atoms_flag = true;
   bool draw_missing_loops_flag = true;
   bool force_rebonding = true;
   make_colour_by_chain_bonds(geom_p, dummy, change_c_only_flag, goodsell_mode, draw_hydrogen_atoms_flag, draw_missing_loops_flag, do_rota_markup, rpt, force_rebonding);

   if (false)
      std::cout << "DEBUG:: in get_rotamer_dodecs() bonds_box.n_rotamer_markups " << bonds_box.n_rotamer_markups
                << std::endl;

   if (bonds_box.n_rotamer_markups > 0) {

      auto &vertices = m.vertices;
      auto &triangles = m.triangles;

      glm::vec4 col(0.6, 0.2, 0.8, 1.0); // starting colour
      dodec d;
      std::vector<clipper::Coord_orth> coords = d.coords();
      std::vector<glm::vec3> dodec_postions(coords.size());
      for (unsigned int i=0; i<coords.size(); i++)
         dodec_postions[i] = clipper_to_glm(coords[i]);

      std::vector<coot::api::vnc_vertex> dodec_vertices;
      std::vector<g_triangle> dodec_triangles;
      dodec_triangles.reserve(36);

      for (unsigned int iface=0; iface<12; iface++) {

         std::vector<coot::api::vnc_vertex> face_verts;
         std::vector<g_triangle> face_triangles;
         face_triangles.reserve(3);

         std::vector<unsigned int> indices_for_face = d.face(iface);
         glm::vec3 ns(0,0,0);
         for (unsigned int j=0; j<5; j++)
            ns += dodec_postions[indices_for_face[j]];
         glm::vec3 normal = glm::normalize(ns);

         for (unsigned int j=0; j<5; j++) {
            glm::vec3 &pos = dodec_postions[indices_for_face[j]];
            coot::api::vnc_vertex v(0.5f * pos, normal, col);
            face_verts.push_back(v);
         }

         face_triangles.push_back(g_triangle(0,1,2));
         face_triangles.push_back(g_triangle(0,2,3));
         face_triangles.push_back(g_triangle(0,3,4));

         unsigned int idx_base = dodec_vertices.size();
         unsigned int idx_tri_base = dodec_triangles.size();
         dodec_vertices.insert(dodec_vertices.end(), face_verts.begin(), face_verts.end());
         dodec_triangles.insert(dodec_triangles.end(), face_triangles.begin(), face_triangles.end());
         for (unsigned int jj=idx_tri_base; jj<dodec_triangles.size(); jj++)
            dodec_triangles[jj].rebase(idx_base);
      }

      // now there is a dodec at the origin, dodec_vertices and dodec_triangle

      // let's make copies of that and move them around to the residues

      double rama_ball_pos_offset_scale = 1.5; // may need tweaking, (was 1.2)

      if (false)
         std::cout << "DEBUG:: in get_rotamer_dodecs() there were " << bonds_box.n_rotamer_markups
                   << " rotamer markups " << std::endl;

      for (int i=0; i<bonds_box.n_rotamer_markups; i++) {
         const rotamer_markup_container_t &rm = bonds_box.rotamer_markups[i];
         const residue_spec_t &residue_spec = rm.spec;
         mmdb::Residue *residue_p = get_residue(residue_spec);
         Cartesian offset(0,0,rama_ball_pos_offset_scale);
         if (residue_p) {
            std::pair<bool, coot::Cartesian> hav = get_HA_unit_vector(residue_p);
            if (hav.first) offset = hav.second * 1.6;
         }

         glm::vec3 atom_pos = clipper_to_glm(rm.pos) + cartesian_to_glm(offset);
         // 20221126-PE Calm down the ultra-bright rota dodec:
         auto rm_col = rm.col;
         rm_col.scale_intensity(0.75); // was 0.6 in Mesh-from-graphical-bonds.cc
         auto this_dodec_colour = colour_holder_to_glm(rm_col);

         std::vector<coot::api::vnc_vertex> this_dodec_vertices = dodec_vertices; // at the origin to start

         // now move it.
         for (unsigned int j=0; j<this_dodec_vertices.size(); j++) {
            auto &vertex = this_dodec_vertices[j];
            vertex.pos  += atom_pos;
            vertex.normal = -vertex.normal;  // 20221018-PE reverse the normal - hmm.
            vertex.color = this_dodec_colour;
            if (false)
               std::cout << "DEBUG:: in get_rotamer_dodecs() atom_pos " << glm::to_string(vertex.pos)
                         << " rama_markup_col " << rm.col
                         << " color " << glm::to_string(vertex.color) << std::endl;
         }

         // fill the colour map with the colour for this dodec
         m.colour_index_to_colour_map[i] = this_dodec_colour;

         // 20230110-PE  Hmmm colour_index is no longer a member of g_triangle (there is another class:
         // g_triangle_with_colour_index).
         //
         // change the colour index of dodec_triangles
         // for (auto &tri : dodec_triangles)
         // tri.colour_index = i;

         unsigned int idx_base = vertices.size();
         unsigned int idx_tri_base = triangles.size();
         vertices.insert(vertices.end(), this_dodec_vertices.begin(), this_dodec_vertices.end());
         triangles.insert(triangles.end(), dodec_triangles.begin(), dodec_triangles.end());
         for (unsigned int jj=idx_tri_base; jj<triangles.size(); jj++)
            triangles[jj].rebase(idx_base);
      }
   }
   return m;
}

coot::instanced_mesh_t
coot::molecule_t::get_rotamer_dodecs_instanced(protein_geometry *geom_p, rotamer_probability_tables *rpt) {

   auto cartesian_to_glm = [] (const coot::Cartesian &c) {
      return glm::vec3(c.x(), c.y(), c.z()); };

   auto clipper_to_glm = [] (const clipper::Coord_orth &c) {
      return glm::vec3(c.x(), c.y(), c.z()); };

   instanced_mesh_t m;

   std::set<int> dummy;
   bool do_rota_markup = true;
   bool change_c_only_flag = true;
   bool goodsell_mode = false;
   bool draw_hydrogen_atoms_flag = true;
   bool draw_missing_loops_flag = true;
   bool force_rebonding = true;
   make_colour_by_chain_bonds(geom_p, dummy, change_c_only_flag, goodsell_mode, draw_hydrogen_atoms_flag, draw_missing_loops_flag, do_rota_markup, rpt, force_rebonding);

   if (false)
      std::cout << "DEBUG:: in get_rotamer_dodecs_instanced() bonds_box.n_rotamer_markups "
                << bonds_box.n_rotamer_markups << std::endl;

   dodec d;
   std::vector<clipper::Coord_orth> coords = d.coords();
   std::vector<glm::vec3> dodec_postions(coords.size());
   for (unsigned int i=0; i<coords.size(); i++)
      dodec_postions[i] = clipper_to_glm(coords[i]);

   std::vector<coot::api::vn_vertex> dodec_vertices;
   std::vector<g_triangle> dodec_triangles;
   dodec_triangles.reserve(36);

   for (unsigned int iface=0; iface<12; iface++) {

      std::vector<coot::api::vn_vertex> face_verts;
      std::vector<g_triangle> face_triangles;
      face_triangles.reserve(3);

      std::vector<unsigned int> indices_for_face = d.face(iface);
      glm::vec3 ns(0,0,0);
      for (unsigned int j=0; j<5; j++)
         ns += dodec_postions[indices_for_face[j]];
      glm::vec3 normal = glm::normalize(ns);

      for (unsigned int j=0; j<5; j++) {
         glm::vec3 &pos = dodec_postions[indices_for_face[j]];
         coot::api::vn_vertex v(0.5f * pos, normal);
         face_verts.push_back(v);
      }

      face_triangles.push_back(g_triangle(0,1,2));
      face_triangles.push_back(g_triangle(0,2,3));
      face_triangles.push_back(g_triangle(0,3,4));

      unsigned int idx_base = dodec_vertices.size();
      unsigned int idx_tri_base = dodec_triangles.size();
      dodec_vertices.insert(dodec_vertices.end(), face_verts.begin(), face_verts.end());
      dodec_triangles.insert(dodec_triangles.end(), face_triangles.begin(), face_triangles.end());
      for (unsigned int jj=idx_tri_base; jj<dodec_triangles.size(); jj++)
         dodec_triangles[jj].rebase(idx_base);
   }

   instanced_geometry_t ig(dodec_vertices, dodec_triangles);
   double rama_ball_pos_offset_scale = 1.5;
   float size = 1.0;
   glm::vec3 size_3(size, size, size);

   // std::cout << "DEBUG:: in coot::molecule_t::get_rotamer_dodecs_instanced(): n_rotamer_markups: "
   //           << bonds_box.n_rotamer_markups << std::endl;

   for (int i=0; i<bonds_box.n_rotamer_markups; i++) {
      // because of the way they are sized, the rotamer markups can contain non-valid rotamers.
      const rotamer_markup_container_t &rm = bonds_box.rotamer_markups[i];
      if (rm.rpi.state == coot::rotamer_probability_info_t::OK) { // should be coot::rotamer_probability_info_t::OK
         const residue_spec_t &residue_spec = rm.spec;
         mmdb::Residue *residue_p = get_residue(residue_spec);
         Cartesian offset(0,0,rama_ball_pos_offset_scale);
         if (residue_p) {
            std::pair<bool, coot::Cartesian> hav = get_HA_unit_vector(residue_p);
            if (hav.first) offset = hav.second * 1.6;
         }
         glm::vec3 atom_pos = clipper_to_glm(rm.pos) + cartesian_to_glm(offset);
         auto rm_col = rm.col;
         rm_col.scale_intensity(0.75); // was 0.6 in Mesh-from-graphical-bonds.cc
         auto this_dodec_colour = colour_holder_to_glm(rm_col);
         instancing_data_type_A_t id(atom_pos, this_dodec_colour, size_3);
         ig.instancing_data_A.push_back(id);
      }
   }

   m.add(ig);
   // std::cout << "in cm::get_rotamer_dodecs_instanced() geom size is " << m.geom.size()
   // << " sending back " << m.geom[0].instancing_data_A.size() << " in 0th" << std::endl;
   return m;
}


#include "ligand/backrub-rotamer.hh"

std::pair<bool,float>
coot::molecule_t::backrub_rotamer(const std::string &chain_id, int res_no,
                                  const std::string &ins_code, const std::string &alt_conf,
                                  const clipper::Xmap<float> &xmap_in,
                                  const coot::protein_geometry &pg) {

   // this doesn't check alt conf - perhaps it should.
   auto move_atoms = [] (const minimol::residue &mres, mmdb::Residue *r) {
      for(const auto &atom : mres.atoms) {
         int n_atoms = 0;
         mmdb::Atom **residue_atoms = 0;
         r->GetAtomTable(residue_atoms, n_atoms);
         for (int i=0; i<n_atoms; i++) {
            mmdb::Atom *at = residue_atoms[i];
            std::string atom_name(at->GetAtomName());
            if (atom_name == atom.name) {
               at->x = atom.pos.x();
               at->y = atom.pos.y();
               at->z = atom.pos.z();
            }
         }
      }
   };

   auto move_backrub_atoms = [move_atoms] (mmdb::Residue *prev_res,
                                 mmdb::Residue *this_res,
                                 mmdb::Residue *next_res,
                                 const minimol::fragment &frag) {
      for (mmdb::Residue *r : {prev_res, this_res, next_res}) {
         if (r) {
            std::string chain_id;
            int res_no = r->GetSeqNum();
            const minimol::residue &mres = frag[res_no];
            int n_atoms = mres.atoms.size();
            // std::cout << "found " << n_atoms << " atoms in minimol residue for " << res_no << std::endl;
            if (n_atoms > 0) {
               move_atoms(mres, r);
            }
         }
      }
   };

   bool status = false;
   float score = -1;
   bool refinement_move_atoms_with_zero_occupancy_flag = true; // pass this?

   // std::cout << "debug:: molecule_t::backrub_rotamer() at " << chain_id << " " << res_no << std::endl;

   residue_spec_t res_spec(chain_id, res_no, ins_code);
   mmdb::Residue *res = get_residue(res_spec);
   if (! res) {
      std::cout << "   WARNING:: residue in molecule :" << chain_id << ": "
                << res_no << " inscode :" << ins_code << ": altconf :"
                << alt_conf << ":" << std::endl;
   } else {
      std::string monomer_type = res->GetResName();
      std::pair<short int, coot::dictionary_residue_restraints_t> p =
         pg.get_monomer_restraints(monomer_type, imol_no);
      coot::dictionary_residue_restraints_t restraints = p.second;

      if (p.first) {

         make_backup("backrub_rotamer");

         try {

            mmdb::Residue *prev_res = coot::util::previous_residue(res);
            mmdb::Residue *next_res = coot::util::next_residue(res);
            mmdb::Manager *mol = atom_sel.mol;
            coot::backrub br(chain_id, res, prev_res, next_res, alt_conf, mol,
                             &xmap_in); // use a const pointer for the map
            // std::cout << "------------ done making a backrub" << std::endl;
            // std::cout << "------------ calling br.search()" << std::endl;
            std::pair<coot::minimol::molecule,float> m = br.search(restraints);
            std::vector<coot::atom_spec_t> baddie_waters = br.waters_for_deletion();
            int ich = 0;

            move_backrub_atoms(prev_res, res, next_res, m.first[ich]);
            status = true; // 20240727-PE oops I deleted this when I wrote move_backrub_atoms()
                           // and forgot to restore it

            if (baddie_waters.size())
               delete_atoms(baddie_waters);

            atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
            atom_sel.mol->FinishStructEdit();
         }
         catch (const std::runtime_error &rte) {
            std::cout << "WARNING:: in backrub_rotamer(): thrown " << rte.what()
                      << " with status " << status << std::endl;
         }
         // if we make a backup, then we also make a new modification
         // save_info.new_modification("backrub_rotamer()");
      } else {
         std::cout << " No restraints found for " << monomer_type << std::endl;
      }
   }
   return std::pair<bool,float> (status, score);
}

std::pair<bool,float>
coot::molecule_t::backrub_rotamer(mmdb::Residue *residue_p,
                                  const clipper::Xmap<float> &xmap,
                                  const coot::protein_geometry &pg) {

   std::string alt_conf;
   return backrub_rotamer(residue_p->GetChainID(),
                          residue_p->GetSeqNum(),
                          residue_p->GetInsCode(),
                          alt_conf,
                          xmap, pg);

}




int
coot::molecule_t::auto_fit_rotamer(const std::string &chain_id, int res_no, const std::string &ins_code,
                                   const std::string &alt_conf,
                                   const clipper::Xmap<float> &xmap, const coot::protein_geometry &geom) {

   int status = 0;

   // let's just use backrub rotamer for now.
   std::pair<bool, float> r = backrub_rotamer(chain_id, res_no, ins_code, alt_conf, xmap, geom);
   status = r.first;

   return status;
}


int
coot::molecule_t::delete_side_chain(const residue_spec_t &residue_spec) {

   int status = 0;
   mmdb::Residue *residue_p = get_residue(residue_spec);
   if (residue_p) {

      make_backup("delete_side_chain");
      bool was_deleted = false;
      // do we include CB? I forget.
      std::vector<std::string> main_chain_atoms_list = { " C  ", " N  ", " H  ", " O  ", " CA ",  " HA ", " CB " };
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::vector<mmdb::Atom *> atoms_to_be_deleted;
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         std::string atom_name(at->GetAtomName());
         if (std::find(main_chain_atoms_list.begin(), main_chain_atoms_list.end(), atom_name) == main_chain_atoms_list.end()) {
            atoms_to_be_deleted.push_back(at);
         }
      }

      for (auto &at : atoms_to_be_deleted) {
         delete at;
         was_deleted = true;
      }

      if (was_deleted) {
         status = true;
         atom_sel.mol->FinishStructEdit();
         atom_sel = make_asc(atom_sel.mol);
         // make_bonds_type_checked(__FUNCTION__);
         // have_unsaved_changes_flag = 1;
         // save_info.new_modification("delete_side_chain()");
         // unlikely to be necessary:
         trim_atom_label_table();
      }
   }
   return status;
}


int
coot::molecule_t::delete_atoms(const std::vector<coot::atom_spec_t> &atom_specs) {

   short int was_deleted = 0;
   int n_deleleted_atoms = 0;

   if (atom_sel.n_selected_atoms > 0) {
      if (atom_specs.size() > 0) {
         make_backup("delete_atoms");
         for (unsigned int i=0; i<atom_specs.size(); i++) {
            int SelHnd = atom_sel.mol->NewSelection();
            // how about a function that calls this:
            // select_atomspec_atoms(atom_sel.mol, atom_specs[i])
            // or  a member function of an atom_spec_t:
            //    atom_specs[i].select_atoms(mol)
            //
            mmdb::PAtom *atoms = NULL;
            int n_atoms;
            atom_sel.mol->SelectAtoms(SelHnd, 0, atom_specs[i].chain_id.c_str(),
                                      atom_specs[i].res_no, atom_specs[i].ins_code.c_str(),
                                      atom_specs[i].res_no, atom_specs[i].ins_code.c_str(),
                                      "*",
                                      atom_specs[i].atom_name.c_str(),
                                      "*",
                                      atom_specs[i].alt_conf.c_str()
                                      );
            atom_sel.mol->GetSelIndex(SelHnd, atoms, n_atoms);
            if (n_atoms) {
               delete atoms[0];
               atoms[0] = NULL;
               n_deleleted_atoms++;
               was_deleted = 1;
            }
            atom_sel.mol->DeleteSelection(SelHnd);
         }
      }

      // potentially
      if (was_deleted) {
         atom_sel.mol->FinishStructEdit();
         atom_sel = make_asc(atom_sel.mol);
         // make_bonds_type_checked(__FUNCTION__);
         // have_unsaved_changes_flag = 1;
         // save_info.new_modification("delete_atoms");
         // unlikely to be necessary:
         trim_atom_label_table();

         // update_symmetry();
      }
   }
   return n_deleleted_atoms;
}

void
coot::molecule_t::update_symmetry() {

   // 20221013-PE no symmetry yet
}

void
coot::molecule_t::trim_atom_label_table() {

   // 20221013-PE no atom labels yet
}

void
coot::molecule_t::delete_ghost_selections() {

   // 20221013-PE no ghosts at the moment
}

int
coot::molecule_t::delete_atom(coot::atom_spec_t &atom_spec) {

   int was_deleted = 0;
   mmdb::Chain *chain;
   mmdb::Residue *residue_of_deleted_atom = NULL;

   std::string chain_id = atom_spec.chain_id;
   int resno            = atom_spec.res_no;
   std::string ins_code = atom_spec.ins_code;
   std::string atname   = atom_spec.atom_name;
   std::string altconf  = atom_spec.alt_conf;

   // run over chains of the existing mol
   int nchains = atom_sel.mol->GetNumberOfChains(1);
   for (int ichain=0; ichain<nchains; ichain++) {

      chain = atom_sel.mol->GetChain(1,ichain);
      std::string mol_chain_id(chain->GetChainID());

      // Note, if in the PDB file, the chain id is not set to
      // something, A, B, C etc, then the chain id from mmdb is ""
      // (not " "!)

//       std::cout << "debug:: delete_atom comparing chain_ids :"
//                 << chain_id << ": vs :" << mol_chain_id << ":"
//                 << std::endl;

      if (chain_id == mol_chain_id) {

         int nres = chain->GetNumberOfResidues();
         for (int ires=0; ires<nres; ires++) {
            mmdb::PResidue res = chain->GetResidue(ires);
            std::string ins_code_local = res->GetInsCode();

            if (res) {
//                std::cout << "debug:: delete_atom: comparing residue "
//                          << res->GetSeqNum() << " :" << ins_code_local
//                          << ":" << " to " << resno << " :" << ins_code << ":"
//                          << std::endl;
               if (res->GetSeqNum() == resno) {
                  if (ins_code_local == ins_code) {

                     // so we have a matching residue:
                     // std::cout << "debug:: delete_atom: we have a matching residue "
                     // << resno << " :" << ins_code << ":" << std::endl;

                     mmdb::PPAtom residue_atoms;
                     int nResidueAtoms;
                     std::string mol_atom_name;
                     res->GetAtomTable(residue_atoms, nResidueAtoms);
                     for (int iat=0; iat<nResidueAtoms; iat++) {

                        mol_atom_name = residue_atoms[iat]->name;
                        if (atname == mol_atom_name) {

                           if (std::string(residue_atoms[iat]->altLoc) == altconf) {

                              make_backup("delete_atom");
                              atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
                              delete_ghost_selections();
                              res->DeleteAtom(iat);
                              was_deleted = 1;
                              residue_of_deleted_atom = res;
                              break;
                           }
                        }
                     }
                  }
               }
            }
            if (was_deleted)
               break;
         }
      }
      if (was_deleted)
         break;
   }

   // potentially
   if (was_deleted) {
      atom_sel.mol->FinishStructEdit();


      // Now reset/recalculate the occupancy and the altLoc of the
      // remaining atoms in the residue with the same atom name.
      //
      mmdb::PPAtom atoms = NULL;
      int n_atoms;
      mmdb::Atom *at = 0;
      int n_matching_name = 0;
      residue_of_deleted_atom->GetAtomTable(atoms, n_atoms);
      for (int iat=0; iat<n_atoms; iat++) {
         std::string res_atom_name = atoms[iat]->name;
         if (res_atom_name == atname) {
            at = atoms[iat];
            n_matching_name++;
         }
      }
      if (n_matching_name == 1) { // one atom of this name left in the residue, so
                                    // remove its altconf string
         if (at) {
            strncpy(at->altLoc, "", 2);
            // set the occupancy to 1.0 of the remaining atom if it was not zero.
            if (at->occupancy > 0.009)
               at->occupancy = 1.0;
         }
      }

      atom_sel = make_asc(atom_sel.mol);
      // make_bonds_type_checked(__FUNCTION__);
      // have_unsaved_changes_flag = 1;
      // save_info.new_modification("delete_atom");
      // unlikely to be necessary:
      trim_atom_label_table();
      update_symmetry();
   }
   return was_deleted;
}


void
coot::molecule_t::delete_any_link_containing_residue(const coot::residue_spec_t &res_spec) {

   // this doesn't do a backup - the calling function is in charge of that

   if (atom_sel.mol) {
      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {
         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
         if ((res_spec.model_number == imod) || (res_spec.model_number == mmdb::MinInt4)) {
            mmdb::LinkContainer *links = model_p->GetLinks();
            int n_links = model_p->GetNumberOfLinks();
            for (int ilink=1; ilink<=n_links; ilink++) {
               mmdb::Link *link_p = model_p->GetLink(ilink);
               // mmdb::Link *link = static_cast<mmdb::Link *>(links->GetContainerClass(ilink));

               if (link_p) {

                  // must pass a valid link_p
                  std::pair<coot::atom_spec_t, coot::atom_spec_t> link_atoms = coot::link_atoms(link_p, model_p);
                  coot::residue_spec_t res_1(link_atoms.first);
                  coot::residue_spec_t res_2(link_atoms.second);
                  // std::cout << "found link " << res_1 << " to "  << res_2 << std::endl;
                  if (res_spec == res_1) {
                     delete_link(link_p, model_p);
                  }
                  if (res_spec == res_2) {
                     delete_link(link_p, model_p);
                  }
               } else {
                  std::cout << "ERROR:: Null link_p for link " << ilink << " of " << n_links << std::endl;
               }
            }
         }
      }
   }
}

void
coot::molecule_t::delete_link(mmdb::Link *link, mmdb::Model *model_p) {

   // Copy out the links, delete all links and add the saved links back

   std::vector<mmdb::Link *> saved_links;
   int n_links = model_p->GetNumberOfLinks();
   for (int ilink=1; ilink<=n_links; ilink++) {
      mmdb::Link *model_link = model_p->GetLink(ilink);
      if (model_link != link) {
         mmdb::Link *copy_link = new mmdb::Link(*model_link);
         saved_links.push_back(copy_link);
      }
   }

   model_p->RemoveLinks();
   for (unsigned int i=0; i<saved_links.size(); i++) {
      model_p->AddLink(saved_links[i]);
   }
}

int
coot::molecule_t::delete_residue(coot::residue_spec_t &residue_spec) {

   // this body copied from
   // short int
   // molecule_class_info_t::delete_residue(int model_number,
   //              const std::string &chain_id, int resno,
   //              const std::string &ins_code)

   int was_deleted = 0;
   if (atom_sel.mol) {

      mmdb::Chain *chain;

      // run over chains of the existing mol
      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {

         if (0)
            std::cout << "debug:: delete_residue() comparing imod: "
                      << imod << " and model_number "
                      << residue_spec.model_number << std::endl;

         if ((imod == residue_spec.model_number) || (residue_spec.model_number == mmdb::MinInt4)) {

            int nchains = atom_sel.mol->GetNumberOfChains(imod);
            for (int ichain=0; ichain<nchains; ichain++) {

               chain = atom_sel.mol->GetChain(imod,ichain);
               std::string mol_chain_id(chain->GetChainID());

               if (residue_spec.chain_id == mol_chain_id) {

                  // std::cout << "debug:: matching chain_ids on  " << chain_id << std::endl;

                  int nres = chain->GetNumberOfResidues();
                  for (int ires=0; ires<nres; ires++) {
                     mmdb::Residue *res = chain->GetResidue(ires);
                     if (res) {
                        if (res->GetSeqNum() == residue_spec.res_no) {

                           // so we have a matching residue:
                           int iseqno = res->GetSeqNum();
                           mmdb::pstr inscode = res->GetInsCode();
                           std::string inscodestr(inscode);
                           if (residue_spec.ins_code == inscodestr) {
                              make_backup("delete_residue");
                              atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
                              delete_ghost_selections();
                              chain->DeleteResidue(iseqno, inscode);
                              was_deleted = 1;
                              res = NULL;
                              break;
                           }
                        }
                     }
                  }
               }
               if (was_deleted)
                  break;
            }
         }
      }
   }

   if (was_deleted) {

      // we can't do this after the modification: it has to be done before
      // atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);

      atom_sel.atom_selection = NULL;
      coot::residue_spec_t spec(residue_spec.model_number, residue_spec.chain_id, residue_spec.res_no, residue_spec.ins_code);
      delete_any_link_containing_residue(spec);
      atom_sel.mol->FinishStructEdit();
      coot::util::pdbcleanup_serial_residue_numbers(atom_sel.mol);
      atom_sel = make_asc(atom_sel.mol);

      // save_info.new_modification(__FUNCTION__);
      trim_atom_label_table();
      update_symmetry();
   }
   return was_deleted;
}

int
coot::molecule_t::delete_residue_atoms_with_alt_conf(coot::residue_spec_t &residue_spec, const std::string &alt_conf) {

   int status = 0;
   std::vector<mmdb::Atom *> atoms_to_be_deleted;

   mmdb::Residue *residue_p = get_residue(residue_spec);
   if (residue_p) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         std::string al(at->altLoc);
         if (al == alt_conf)
            atoms_to_be_deleted.push_back(at);
      }
   }
   bool was_deleted = false;
   for (auto &at : atoms_to_be_deleted) {
      delete at;
      was_deleted = true;
   }
   if (was_deleted) {
      status = true;
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
   }
   return status;
}


int
coot::molecule_t::delete_chain_using_atom_cid(const std::string &cid) {

   int done = 0;

   std::pair<bool, atom_spec_t> spec_pair = cid_to_atom_spec(cid);
   if (spec_pair.first) {
      const auto &spec = spec_pair.second;
      const std::string &chain_id =  spec.chain_id;
      for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               if (chain_p) {
                  std::string this_chain_id(chain_p->GetChainID());
                  if (this_chain_id == chain_id) {
                     make_backup("delete_chain_using_atom_cid");
                     model_p->DeleteChain(ichain);
                     done = true;
                  }
               }
            }
         }
      }
   }

   if (done) {
      atom_sel.mol->FinishStructEdit();
      // update_molecule_after_additions();
   }

   return done;
}

int
coot::molecule_t::delete_literal_using_cid(const std::string &atom_selection_cids) {

   // cid is an atom selection, e.g. containing a residue range

   int status = 0;
   std::vector<mmdb::Atom *> atoms_to_be_deleted;

   mmdb::Atom **selection_atoms = 0;
   int n_selection_atoms = 0;
   int selHnd = atom_sel.mol->NewSelection(); // d
   std::vector<std::string> v = coot::util::split_string(atom_selection_cids, "||");
   if (! v.empty())
      for (const auto &cid : v)
         atom_sel.mol->Select(selHnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_OR);
   atom_sel.mol->GetSelIndex(selHnd, selection_atoms, n_selection_atoms);

   if (selection_atoms) {
      for (int iat=0; iat<n_selection_atoms; iat++) {
         mmdb::Atom *at = selection_atoms[iat];
         if (at)
            atoms_to_be_deleted.push_back(at);
      }
   }

   if (! atoms_to_be_deleted.empty()) {
      std::string s = std::string("delete-literal-using-cid ") + atom_selection_cids;
      make_backup(s);

      for (unsigned int iat=0; iat<atoms_to_be_deleted.size(); iat++) {
         delete atoms_to_be_deleted[iat];
         atoms_to_be_deleted[iat] = NULL;
      }
      status = 1;
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      coot::util::pdbcleanup_serial_residue_numbers(atom_sel.mol);

            // save_info.new_modification(s);
   }
   return status;
}

#include "geometry/main-chain.hh"

int
coot::molecule_t::change_alt_locs(const std::string &cid, const std::string &change_mode) {

   auto change_residue_alt_locs = [] (mmdb::Residue *residue_p) {
      int status = 0;
      std::vector<mmdb::Atom *> atoms_with_alt_loc_A;
      std::vector<mmdb::Atom *> atoms_with_alt_loc_B;
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         std::string alt_loc = at->altLoc;
         if (alt_loc == "A") atoms_with_alt_loc_A.push_back(at);
         if (alt_loc == "B") atoms_with_alt_loc_B.push_back(at);
      }
      for (auto atom : atoms_with_alt_loc_A) strncpy(atom->altLoc, "B", 2);
      for (auto atom : atoms_with_alt_loc_B) strncpy(atom->altLoc, "A", 2);
      if (! atoms_with_alt_loc_A.empty()) status = 1;
      if (! atoms_with_alt_loc_B.empty()) status = 1;
      return status;
   };

   auto change_side_chain_alt_locs = [] (mmdb::Residue *residue_p) {
      int status = 0;
      std::vector<mmdb::Atom *> atoms_with_alt_loc_A;
      std::vector<mmdb::Atom *> atoms_with_alt_loc_B;
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (is_main_chain_p(at)) {
         } else {
            std::string alt_loc = at->altLoc;
            if (alt_loc == "A") atoms_with_alt_loc_A.push_back(at);
            if (alt_loc == "B") atoms_with_alt_loc_B.push_back(at);
         }
      }
      for (auto atom : atoms_with_alt_loc_A) strncpy(atom->altLoc, "B", 2);
      for (auto atom : atoms_with_alt_loc_B) strncpy(atom->altLoc, "A", 2);
      if (! atoms_with_alt_loc_A.empty()) status = 1;
      if (! atoms_with_alt_loc_B.empty()) status = 1;
      return status;
   };

   auto change_main_chain_alt_locs = [] (mmdb::Residue *residue_p) {
      int status = 0;
      std::vector<mmdb::Atom *> atoms_with_alt_loc_A;
      std::vector<mmdb::Atom *> atoms_with_alt_loc_B;
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (is_main_chain_p(at)) {
            std::string alt_loc = at->altLoc;
            if (alt_loc == "A") atoms_with_alt_loc_A.push_back(at);
            if (alt_loc == "B") atoms_with_alt_loc_B.push_back(at);
         }
      }
      for (auto atom : atoms_with_alt_loc_A) strncpy(atom->altLoc, "B", 2);
      for (auto atom : atoms_with_alt_loc_B) strncpy(atom->altLoc, "A", 2);
      if (! atoms_with_alt_loc_A.empty()) status = 1;
      if (! atoms_with_alt_loc_B.empty()) status = 1;
      return status;
   };

   auto change_atom_list_alt_locs = [] (mmdb::Residue *residue_p, const std::string &change_mode) {
      int status = 0;
      std::vector<std::string> names = util::split_string(change_mode, ",");
      if (! names.empty()) {
         for (const auto &n : names) {
            std::string name = n;
            // hacky
            if (n.length() == 3) name = std::string(" ") + n;
            if (n.length() == 2) name = std::string(" ") + n + std::string(" ");
            if (n.length() == 1) name = std::string(" ") + n + std::string("  ");
            mmdb::Atom **residue_atoms = 0;
            int n_residue_atoms = 0;
            residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
            for (int iat=0; iat<n_residue_atoms; iat++) {
               mmdb::Atom *at = residue_atoms[iat];
               std::string atom_name(at->name);
               if (atom_name == name) {
                  std::string alt_loc = at->altLoc;
                  if (alt_loc == "A") {
                     strncpy(at->altLoc, "B", 2);
                     status = 1;
                  } else {
                     if (alt_loc == "B") {
                        strncpy(at->altLoc, "A", 2);
                        status = 1;
                     }
                  }
               }
            }
         }
      }
      return status;
   };

   int status = 0;
   std::vector<mmdb::Residue *> residues = cid_to_residues(cid);
   if (! residues.empty()) {
      for(auto residue : residues) {
         if (change_mode == "residue") {
            status = change_residue_alt_locs(residue);
         } else {
            if (change_mode == "side-chain") {
               status = change_side_chain_alt_locs(residue);
            } else {
               if (change_mode == "main-chain") {
                  status = change_main_chain_alt_locs(residue);
               } else {
                  status = change_atom_list_alt_locs(residue, change_mode);
               }
            }
         }
      }
   }
   return status;
}

std::vector<std::string>
coot::molecule_t::non_standard_residue_types_in_model() const {

   std::vector<std::string> v;
   if (atom_sel.mol) {
      v = coot::util::non_standard_residue_types_in_molecule(atom_sel.mol);
   }
   return v;
}

int
coot::molecule_t::sfcalc_genmap(const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                                const clipper::HKL_data<clipper::data32::Flag> &free,
                                clipper::Xmap<float> *xmap_p) {

   bool sane = sanity_check_atoms(atom_sel.mol);

   if (sane) {
      coot::util::sfcalc_genmap(atom_sel.mol, fobs, free, xmap_p);
   } else {
      std::cout << "ERROR:: coordinates were not sane" << std::endl;
   }
   return 0;
}

coot::util::sfcalc_genmap_stats_t
coot::molecule_t::sfcalc_genmaps_using_bulk_solvent(const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                                                         const clipper::HKL_data<clipper::data32::Flag> &free,
                                                         clipper::Xmap<float> *xmap_2fofc_p,
                                                         clipper::Xmap<float> *xmap_fofc_p) {

   coot::util::sfcalc_genmap_stats_t stats;
   bool sane = sanity_check_atoms(atom_sel.mol);
   if (sane) {

      clipper::Cell cell = xmap_2fofc_p->cell();
      if (true) {
         // sanity check data
         const clipper::HKL_info &hkls_check = fobs.base_hkl_info();
         const clipper::Spacegroup &spgr_check = hkls_check.spacegroup();

         if (false) {
            std::cout << "DEBUG:: Sanity check A in molecule_t::sfcalc_genmaps_using_bulk_solvent(): HKL_info: "
                      << "cell: " << hkls_check.cell().format() << " "
                      << "spacegroup: " << spgr_check.symbol_xhm() << " "
                      << "resolution: " << hkls_check.resolution().limit() << " "
                      << "invsqreslim: " << hkls_check.resolution().invresolsq_limit() << " "
                      << std::endl;
            std::cout << "DEBUG:: Sanity check B in molecule_t::sfcalc_genmaps_using_bulk_solvent(): Cell fofc-map"
                      << xmap_fofc_p->cell().format() << std::endl;
         }
      }

      stats = coot::util::sfcalc_genmaps_using_bulk_solvent(atom_sel.mol, fobs, free, cell, xmap_2fofc_p, xmap_fofc_p);

      // maybe format() should be inside coot::util::sfcalc_genmap_stats_t
      if (false) {
         std::cout << "\n R-factor      : " << stats.r_factor << "\n Free R-factor : " << stats.free_r_factor << "\n";
         std::cout << "\n Bulk Correction Volume: " << stats.bulk_solvent_volume;
         std::cout << "\n Bulk Correction Factor: " << stats.bulk_correction << "\n";
         std::cout << "\nNumber of spline params: " << stats.n_splines << "\n";
      }

   } else {
      std::cout << "ERROR:: coordinates were not sane" << std::endl;
   }
   return stats;
}

bool
coot::molecule_t::sanity_check_atoms(mmdb::Manager *mol) const {

   bool sane = true;
   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (! model_p) {
         std::cout << "ERROR:: Bad model " << imod << std::endl;
         sane = false;
      } else {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            if (! chain_p) {
               std::cout << "ERROR:: Bad chain with index " << ichain << "  in model "
                         << imod << std::endl;
               sane = false;
            } else {
               int nres = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<nres; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (! residue_p) {
                     std::cout << "ERROR:: Bad residue with index " << ires << "  in chain "
                               << chain_p->GetChainID() << std::endl;
                     sane = false;
                  } else {
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     for (int iat=0; iat<n_atoms; iat++) {
                        mmdb::Atom *at = residue_p->GetAtom(iat);
                        if (! at) {
                           std::cout << "ERROR:: Bad atom with index " << iat << "  in residue "
                                     << coot::residue_spec_t(residue_p) << std::endl;
                           sane = false;
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return sane;
}

// refinement is done in the container
// int
// coot::molecule_t::refine_residues(const residue_spec_t &residue_spec, refine_residues_mode mode) {

//    int status = 0;

//    return status;
// }

std::vector<coot::atom_spec_t>
coot::molecule_t::get_fixed_atoms() const {

   return fixed_atom_specs;
}

#include "ideal/simple-restraint.hh"

int
coot::molecule_t::refine_direct(std::vector<mmdb::Residue *> rv, const std::string &alt_loc, const clipper::Xmap<float> &xmap,
                                unsigned int max_number_of_threads,
                                float map_weight, int n_cycles, const coot::protein_geometry &geom,
                                bool do_rama_plot_restraints, float rama_plot_weight,
                                bool do_torsion_restraints, float torsion_weight,
                                bool refinement_is_quiet) {

   bool make_trans_peptide_restraints = true;

   int status =  0;
   std::vector<coot::atom_spec_t> fixed_atom_specs;
   std::vector<std::pair<bool,mmdb::Residue *> > local_residues;
   for (const auto &r : rv)
      local_residues.push_back(std::make_pair(false, r));

   if (false) { // debugging
      std::cout << "---------- local_residues " << local_residues.size() << " --------" << std::endl;
      for (unsigned int i=0; i<local_residues.size(); i++) {
         std::cout << "                 " << i << " " << local_residues[i].second << std::endl;
      }
   }

   make_backup("refine_direct");
   mmdb::Manager *mol = atom_sel.mol;
   std::vector<mmdb::Link> links;
   coot::restraints_container_t restraints(local_residues,
                                           links,
                                           geom,
                                           mol,
                                           fixed_atom_specs, &xmap);

   if (refinement_is_quiet)
      restraints.set_quiet_reporting();

   if (do_rama_plot_restraints) {
      restraints.set_rama_type(RAMA_TYPE_ZO);
      restraints.set_rama_plot_weight(rama_plot_weight);
   }

   if (do_torsion_restraints) {
      restraints.set_torsion_restraints_weight(torsion_weight);
   }

   restraints.set_map_weight(map_weight);
   restraints.add_map(map_weight);

   restraint_usage_Flags flags = TYPICAL_RESTRAINTS;
   if (do_torsion_restraints) flags = TYPICAL_RESTRAINTS_WITH_TORSIONS;
   pseudo_restraint_bond_type pseudos = NO_PSEUDO_BONDS;

   int n_threads = max_number_of_threads;
   ctpl::thread_pool thread_pool(n_threads);
   restraints.thread_pool(&thread_pool, n_threads);

   int imol = imol_no;
   bool do_auto_helix_restraints = true;
   bool do_auto_strand_restraints = false;
   bool do_h_bond_restraints = false;
   bool do_residue_internal_torsions = do_torsion_restraints;
   restraints.make_restraints(imol, geom, flags,
                              do_residue_internal_torsions,
                              make_trans_peptide_restraints,
                              rama_plot_weight, do_rama_plot_restraints,
                              do_auto_helix_restraints,
                              do_auto_strand_restraints,
                              do_h_bond_restraints,
                              pseudos);
   restraints.add_extra_restraints(imol, "stored extra retraints called from refine_direct()",
                                   extra_restraints, geom);
   int nsteps_max = n_cycles;
   short int print_chi_sq_flag = 1;
   restraints.minimize(flags, nsteps_max, print_chi_sq_flag);
   geometry_distortion_info_container_t gd = restraints.geometric_distortions();
   if (! refinement_is_quiet)
      gd.print();
   restraints.unset_fixed_during_refinement_udd();
   status = 1;

   // save_info.new_modification("refine_direct");

   return status;
}

int
coot::molecule_t::refine_using_last_restraints(int n_steps) {

   if (! last_restraints) return 0; // "success" (hmm)

   restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
   flags = coot::TYPICAL_RESTRAINTS; // Add GM here?

   // check here that the refinement hasn't already finished.
   //
   coot::refinement_results_t rr = last_restraints->minimize(flags, n_steps, 0);
   return rr.progress;
}


std::pair<int, std::string>
coot::molecule_t::add_terminal_residue_directly(const residue_spec_t &spec, const std::string &new_res_type,
                                                const coot::protein_geometry &geom,
                                                const clipper::Xmap<float> &xmap,
                                                mmdb::Manager *standard_residues_asc_mol, // for RNA
                                                ctpl::thread_pool &static_thread_pool) {

   std::pair<int, std::string> r;
   mmdb::Residue *residue_p = util::get_residue(spec, atom_sel.mol);
   if (residue_p) {
      if (util::is_nucleotide_by_dict(residue_p, geom)) {
         execute_simple_nucleotide_addition(residue_p, standard_residues_asc_mol);
      } else {
         std::string terminus_type = get_term_type(residue_p, atom_sel.mol);
         float bf_new = default_temperature_factor_for_new_atoms;
         make_backup("add_terminal_residue_directly");
         r = add_terminal_residue(imol_no, terminus_type, residue_p,
                                  atom_sel.mol, atom_sel.UDDAtomIndexHandle,
                                  spec.chain_id, new_res_type,
                                  bf_new, xmap, geom, static_thread_pool);
         atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
         atom_sel.mol->FinishStructEdit();
         util::pdbcleanup_serial_residue_numbers(atom_sel.mol);
         atom_sel = make_asc(atom_sel.mol);
         // save_info.new_modification("add-terminal-residue");
      }
   } else {
      std::cout << "WARNING:: in add_terminal_residue_directly() null residue_p " << std::endl;
   }
   return r;
}


int
coot::molecule_t::mutate(const coot::residue_spec_t &spec, const std::string &new_res_type) {

   make_backup("mutate");
   atom_sel.delete_atom_selection();
   mmdb::Residue *residue_p = coot::util::get_residue(spec, atom_sel.mol);
   int status = coot::util::mutate(residue_p, new_res_type);
   // std::cout << "mutate status " << status << std::endl;

   // 20221121-PE should thise function calls be in coot::util::mutate()? I think so.
   atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol); // regen the atom indices
   // save_info.new_modification("mutate() " + new_res_type);
   return status;

}

#include "ligand/side-chain.hh"

int
coot::molecule_t::side_chain_180(const coot::residue_spec_t &residue_spec, const std::string &alt_conf,
                                 coot::protein_geometry *geom_p) {

   int status = 0;
   mmdb::Residue *residue_p = coot::util::get_residue(residue_spec, atom_sel.mol);
   if (residue_p) {
      // sub functions us coot::protein_geometry *geom_p. Maybe they change it?
      coot::do_180_degree_side_chain_flip(residue_spec, alt_conf, atom_sel.mol, geom_p); // void
      status = 1;
   }
   return status;
}



int
coot::molecule_t::move_molecule_to_new_centre(const coot::Cartesian &new_centre) {

   int status = 0;
   if (is_valid_model_molecule()) {
      std::pair<bool, clipper::Coord_orth> cm = coot::centre_of_molecule(atom_sel.mol);
      if (cm.first) {
         make_backup("move_molecule_to_new_centre");
         coot::Cartesian delta = new_centre - coot::Cartesian(cm.second.x(), cm.second.y(), cm.second.z());
         for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
            mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
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
                              at->x += delta.x();
                              at->y += delta.y();
                              at->z += delta.z();
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   // save_info.new_modification("move_molecule_to_new_centre");
   return status;
}


coot::Cartesian
coot::molecule_t::get_molecule_centre() const {

   auto mmdb_to_cartesian = [] (mmdb::Atom *at) {
      return Cartesian(at->x, at->y, at->z);
   };

   coot::Cartesian c(0,0,0);
   unsigned int n_atoms_total = 0;
   if (is_valid_model_molecule()) {
      std::pair<bool, clipper::Coord_orth> cm = coot::centre_of_molecule(atom_sel.mol);
      if (cm.first) {
         for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
            mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
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
                              coot::Cartesian pt = mmdb_to_cartesian(at);
                              c += pt;
                              n_atoms_total++;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   if (n_atoms_total > 0) {
      float sf = 1.0/static_cast<float>(n_atoms_total);
      c = c * sf;
   }
   return c;
}


// return a non-empty string on a problem
//
std::string
coot::molecule_t::jed_flip_internal(coot::atom_tree_t &tree,
                                    const std::vector<coot::dict_torsion_restraint_t> &interesting_torsions,
                                    const std::string &atom_name,
                                    bool invert_selection) {

   // This function was copied to coot-utils - don't edit this, edit the coot-utils
   // version and call it from here - possibly delete this.
   //
   // But this does have_unsaved_changes_flag and make_backup.

   std::string problem_string;
   unsigned int selected_idx = 0;

   if (interesting_torsions.size() > 0) {

      unsigned int best_fragment_size = 9999;
      if (interesting_torsions.size() > 1) {
         // select the best torsion based on fragment size.
         for (unsigned int i=0; i<interesting_torsions.size(); i++) {

            std::string atn_1 = interesting_torsions[i].atom_id_2_4c();
            std::string atn_2 = interesting_torsions[i].atom_id_3_4c();
            bool reverse = false; // dummy value

            std::pair<unsigned int, unsigned int> p = tree.fragment_sizes(atn_1, atn_2, reverse);

            if (p.first < best_fragment_size) {
               best_fragment_size = p.first;
               selected_idx = i;
            }
            if (p.second < best_fragment_size) {
               best_fragment_size = p.second;
               selected_idx = i;
            }
         }
      }

      const auto &int_tor = interesting_torsions[selected_idx];
      problem_string = jed_flip_internal(tree, int_tor, atom_name, invert_selection); // does a backup
   }
   return problem_string;
}

// return a non-null string on a problem
//
std::string
coot::molecule_t::jed_flip_internal(coot::atom_tree_t &tree,
                                    const coot::dict_torsion_restraint_t &torsion,
                                    const std::string &atom_name,
                                    bool invert_selection) {

   // This function was copied to coot-utils - don't edit this, edit the coot-utils
   // version and call it from here - possibly delete this.
   //
   // But this does have_unsaved_changes_flag and make_backup.

   std::string problem_string;

   make_backup("jed_flip_internal");

   bool reverse = false; // reverse the moving dog<->tail fragment?

   if (invert_selection)
      reverse = true;

   std::string atn_1 = torsion.atom_id_2_4c();
   std::string atn_2 = torsion.atom_id_3_4c();

   if (torsion.atom_id_3_4c() == atom_name) {
      atn_1 = torsion.atom_id_3_4c();
      atn_2 = torsion.atom_id_2_4c();
   }

   int period = torsion.periodicity();

   if (period > 1) {

      double angle = 360/double(period);
      std::pair<unsigned int, unsigned int> p = tree.fragment_sizes(atn_1, atn_2, false);

      if (false) {  // debug
         std::cout << "flip this torsion: " << torsion << std::endl;
         std::cout << "DEBUG:: jed_flip_internal() fragment sizes: " << p.first << " " << p.second
                   << std::endl;
      }

      if (p.first > p.second)
         reverse = !reverse;

      tree.rotate_about(atn_1, atn_2, angle, reverse);
      // have_unsaved_changes_flag = 1;
      // save_info.new_modification(__FUNCTION__);
   } else {
      problem_string = "Selected torsion had a periodicity of ";
      problem_string += clipper::String(period);
   }
   return problem_string;
}




std::string
coot::molecule_t::jed_flip(coot::residue_spec_t &spec,
                           const std::string &atom_name,
                           const std::string &alt_conf,
                           bool invert_selection,
                           coot::protein_geometry *geom) {

   // std::cout << "########## jed_flip() " << spec << " " << atom_name << " " << invert_selection << std::endl;

   // This function was copied to coot-utils - don't edit this, edit the coot-utils
   // version and call it from here - possibly delete this.
   //
   // But this does have_unsaved_changes_flag and make_backup.

   std::string problem_string;

   mmdb::Residue *residue = get_residue(spec);
   if (! residue) {
      std::cout << "WARNING:: no residue " << spec << " found in molecule" << std::endl;
   } else {

      // Does atom atom_name with given alt_conf exist in this residue?
      mmdb::Atom *clicked_atom = 0;
      int clicked_atom_idx = -1;
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      residue->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         std::string an(residue_atoms[iat]->name);
         if (an == atom_name) {
            std::string ac(residue_atoms[iat]->altLoc);
            if (ac == alt_conf) {
               clicked_atom = residue_atoms[iat];
               clicked_atom_idx = iat;
               break;
            }
         }
      }

      if (! clicked_atom) {
         std::cout << "WARNING:: atom \"" << atom_name << "\" not found in residue " << std::endl;
      } else {

         // std::cout << "OK, found clicked atom " << atom_spec_t(clicked_atom) << std::endl;

         // make_backup(); // backup is done by jed_flip_internal

         std::string monomer_type = residue->GetResName();

         std::pair<bool, coot::dictionary_residue_restraints_t> p =
            geom->get_monomer_restraints(monomer_type, imol_no);

         if (! p.first) {
            std::cout << "WARNING residue type " << monomer_type << " not found in dictionary" << std::endl;
         } else {
            bool iht = false;  // include_hydrogen_torsions_flag
            std::vector<coot::dict_torsion_restraint_t> all_torsions = p.second.get_non_const_torsions(iht);

            if (all_torsions.size() == 0) {
               problem_string = "There are no non-CONST torsions for this residue type";
            } else {
               std::vector<std::vector<std::string> > ring_atoms_sets = p.second.get_ligand_ring_list();
               std::vector<coot::dict_torsion_restraint_t> interesting_torsions;
               for (unsigned int it=0; it<all_torsions.size(); it++) {

                  bool is_ring_torsion_flag = all_torsions[it].is_ring_torsion(ring_atoms_sets);
                  if (! all_torsions[it].is_ring_torsion(ring_atoms_sets)) {
                     if (all_torsions[it].atom_id_2_4c() == atom_name)
                        interesting_torsions.push_back(all_torsions[it]);
                     if (all_torsions[it].atom_id_3_4c() == atom_name)
                        interesting_torsions.push_back(all_torsions[it]);
                  }
               }

               if (interesting_torsions.size() == 0) {
                  problem_string = "There are no non-CONST non-ring torsions for this atom";
               } else {

                  // make a constructor?
                  atom_selection_container_t residue_asc;
                  residue_asc.n_selected_atoms = n_residue_atoms;
                  residue_asc.atom_selection = residue_atoms;
                  residue_asc.mol = 0;

                  coot::contact_info contact = coot::getcontacts(residue_asc, alt_conf, monomer_type, imol_no, geom);
                  std::vector<std::vector<int> > contact_indices = contact.get_contact_indices_with_reverse_contacts();

                    if (false) {
                       std::cout << " debug:: =========== in jed_flip() contact indices ======= " << std::endl;
                       for (unsigned int ic1 = 0; ic1 < contact_indices.size(); ic1++) {
                          std::cout << "in jed_flip(): index " << ic1 << " : ";
                          for (unsigned int ic2 = 0; ic2 < contact_indices[ic1].size(); ic2++)
                             std::cout << contact_indices[ic1][ic2] << " ";
                          std::cout << std::endl;
                       }
                    }

                  try {
                     coot::atom_tree_t tree(contact_indices, clicked_atom_idx, residue, alt_conf);
                     problem_string = jed_flip_internal(tree, interesting_torsions, atom_name, invert_selection);
                     atom_sel.mol->FinishStructEdit();

                     // save_info.new_modification("jed_flip");

                  }
                  catch (const std::runtime_error &rte) {
                     std::cout << "ERROR:: run-time-error " << rte.what() << " - giving up" << std::endl;
                  }
                  // make_bonds_type_checked(__FUNCTION__);
               }
            }
         }
      }
   }
   return problem_string;
}


#include "ligand/ligand.hh"

coot::minimol::molecule
coot::molecule_t::eigen_flip_residue(const coot::residue_spec_t &residue_spec) {

   coot::minimol::molecule m;

   mmdb::Residue *res = get_residue(residue_spec);
   if (!res) {
      std::cout << "DEBUG:: residue not found " << residue_spec
                << " in molecule number " << imol_no << std::endl;
   } else {

      make_backup("eigen_flip_residue");
      coot::ligand lig;
      coot::minimol::residue r(res);
      coot::minimol::fragment f(res->GetChainID());
      f.residues.push_back(coot::minimol::residue(res));
      coot::minimol::molecule ligand;
      ligand.fragments.push_back(f);

      ligand_flip_number++;
      if (ligand_flip_number == 4)
         ligand_flip_number = 0;

      lig.install_ligand(ligand);
      m = lig.flip_ligand(ligand_flip_number);

      // have_unsaved_changes_flag = 1;
      // save_info.new_modification(__FUNCTION__);

      replace_coords(make_asc(m.pcmmdbmanager()), 0, 1);
   }
   return m;
}



std::vector<std::string>
coot::molecule_t::chains_in_model() const {

   std::vector<std::string> v;
   if (is_valid_model_molecule())
      v = util::chains_in_molecule(atom_sel.mol);

   return v;
}


std::vector<std::pair<coot::residue_spec_t, std::string> >
coot::molecule_t::get_single_letter_codes_for_chain(const std::string &chain_id) const {

   std::vector<std::pair<coot::residue_spec_t, std::string> > v;
   if (is_valid_model_molecule()) {

      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            std::string chain_id_this(chain_p->GetChainID());
            if (chain_id_this == chain_id) {
               int nres;
               mmdb::PResidue *residue_table = 0;
               chain_p->GetResidueTable(residue_table, nres);
               std::vector<std::pair<mmdb::Residue *, int> > sorted = util::sort_residues_by_seqno(residue_table, nres);

               for (unsigned int ires=0; ires<sorted.size(); ires++) {
                  const auto &r(sorted[ires]);
                  residue_spec_t res_spec(r.first);
                  std::string res_name(r.first->GetResName());
                  std::string s = util::three_letter_to_one_letter_with_specials(res_name);
                  v.push_back(std::make_pair(res_spec, s));
               }
            }
         }
      }
   }
   return v;
}

std::vector<std::string>
coot::molecule_t::get_residue_names_with_no_dictionary(const coot::protein_geometry &geom) const {

   std::vector<std::string> v = geom.residue_names_with_no_dictionary(atom_sel.mol, imol_no);
   return v;
}


int
coot::molecule_t::apply_transformation_to_atom_selection(const std::string &atom_selection_cid,
                                                         int n_atoms,
                                                         clipper::Coord_orth &rotation_centre,
                                                         clipper::RTop_orth &rtop) {

   auto mmdb_to_clipper = [] (mmdb::Atom *at) {
      return clipper::Coord_orth(at->x, at->y, at->z);
   };

   int n_atoms_moved = 0;

   if (is_valid_model_molecule()) {

      mmdb::Atom **selection_atoms = 0;
      int n_selection_atoms = 0;
      int selHnd = atom_sel.mol->NewSelection(); // d
      atom_sel.mol->Select(selHnd, mmdb::STYPE_ATOM, atom_selection_cid.c_str(), mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selHnd, selection_atoms, n_selection_atoms);
      // does the number of atoms that we selected actually match the number of atoms that the caller thinks that
      // we should have selected?

      if (selection_atoms) {
         if (n_selection_atoms == n_atoms) {
            make_backup("apply_transformation_to_atom_selection");
            for (int iat=0; iat<n_selection_atoms; iat++) {
               mmdb::Atom *at = selection_atoms[iat];
               if (! at->isTer()) {
                  clipper::Coord_orth pt = mmdb_to_clipper(at);
                  clipper::Coord_orth p1 = pt - rotation_centre;
                  clipper::Coord_orth p2 = rtop * p1;
                  clipper::Coord_orth p3 = p2 - rotation_centre;
                  at->x = p3.x();
                  at->y = p3.y();
                  at->z = p3.z();
                  n_atoms_moved++;
               }
            }
         } else {
            std::cout << "ERROR in apply_transformation_to_atom_selection() mismatch atom in selection "
                      << n_atoms << " " << n_selection_atoms << std::endl;
         }
      }
      atom_sel.mol->DeleteSelection(selHnd);
      // save_info.new_modification(__FUNCTION__);
   }
   return n_atoms_moved;
}


int
coot::molecule_t::new_positions_for_residue_atoms(const std::string &residue_cid, const std::vector<coot::api::moved_atom_t> &moved_atoms) {

   mmdb::Residue *residue_p = cid_to_residue(residue_cid);
   if (residue_p) {
      return new_positions_for_residue_atoms(residue_p, moved_atoms, true);
   } else {
      std::cout << "ERROR:: in new_positions_for_residue_atoms() failed to find residue " << residue_cid << std::endl;
      return -1;
   }
}

int
coot::molecule_t::new_positions_for_residue_atoms(mmdb::Residue *residue_p, const std::vector<coot::api::moved_atom_t> &moved_atoms,
                                                  bool do_backup) {

   // calling function does the make_backup() - otherwise this function could be called 100s of
   // times for one moved chain.

   int n_atoms_moved = 0;
   if (residue_p) {
      if (do_backup)
         make_backup("new_positions_for_residue_atoms");
      for (unsigned int i=0; i<moved_atoms.size(); i++) {
         const api::moved_atom_t &mva = moved_atoms[i];

         mmdb::Atom **residue_atoms = 0;
         int n_residue_atoms = 0;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               std::string atom_name(at->GetAtomName());
               if (atom_name == mva.atom_name) {
                  std::string alt_conf(at->altLoc);
                  if (alt_conf == mva.alt_conf) {
                     at->x = mva.x;
                     at->y = mva.y;
                     at->z = mva.z;
                     n_atoms_moved++;
                  }
               }
            }
         }
      }
   } else {
      std::cout << "ERROR:: in new_positions_for_residue_atoms() failed to find residue " << std::endl;
   }
   return n_atoms_moved;
}

int
coot::molecule_t::new_positions_for_atoms_in_residues(const std::vector<coot::api::moved_residue_t> &moved_residues) {

   int status = 0;
   if (!moved_residues.empty()) {
      make_backup("new_positions_for_atoms_in_residues");
      for (unsigned int i=0; i<moved_residues.size(); i++) {
         const api::moved_residue_t &mvr = moved_residues[i];
         coot::residue_spec_t res_spec(mvr.chain_id, mvr.res_no, mvr.ins_code);
         mmdb::Residue *residue_p = get_residue(res_spec);
         new_positions_for_residue_atoms(residue_p, mvr.moved_atoms, false); // we make_backup() above.
      }
      // save_info.new_modification(__FUNCTION__);
   }
   return status;
}

// we might want to be looking at another molecule's model (e.g. in difference maps)
coot::residue_spec_t
coot::molecule_t::get_residue_closest_to(mmdb::Manager *mol, const clipper::Coord_orth &co) const {

   residue_spec_t spec;

   double d_best = 999999999999.0;
   int imod = 1;
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
                     double dx = at->x - co.x();
                     double dy = at->y - co.y();
                     double dz = at->z - co.z();
                     double dd = dx * dx + dy * dy + dz * dz;
                     if (dd < d_best) {
                        d_best = dd;
                        spec = residue_spec_t(residue_p);
                     }
                  }
               }
            }
         }
      }
   }
   return spec;
}


#include "coot-utils/merge-molecules.hh"

//! merge molecules - copy the atom of mols into this molecule
int
coot::molecule_t::merge_molecules(const std::vector<mmdb::Manager *> &mols) {

  make_backup("merge_molecules");
  int n_new_atoms = 0;
  int n_pre = atom_sel.n_selected_atoms;
  mmdb::Manager *mol = atom_sel.mol;
  atom_sel.delete_atom_selection();
  coot::merge_molecules(mol, mols);

  if (false) { // delete this when fixed
     int imod = 1;
     mmdb::Model *model_p = mol->GetModel(imod);
     if (model_p) {
        int n_chains = model_p->GetNumberOfChains();
        for (int ichain=0; ichain<n_chains; ichain++) {
           mmdb::Chain *chain_p = model_p->GetChain(ichain);
           std::string chain_id(chain_p->GetChainID());
           std::cout << "found chain with chain-id " << chain_id << std::endl;
        }
     }
  }

  atom_sel = make_asc(mol);
  int n_post = atom_sel.n_selected_atoms;
  n_new_atoms = n_post - n_pre;
  // save_info.new_modification(__FUNCTION__);
  return n_new_atoms;

}



// return status [1 means "usable"] and a chain id [status = 0 when
// there are 2*26 chains...]
//
std::pair<bool, std::string>
coot::molecule_t::unused_chain_id() const {

   std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
   std::pair<bool, std::string> s(false,"");
   mmdb::Chain *chain_p;
   if (atom_sel.n_selected_atoms > 0) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(1);
      int nchains = model_p->GetNumberOfChains();

      for (int ich=0; ich<nchains; ich++) {
         chain_p = model_p->GetChain(ich);
         std::string this_chain_id = chain_p->GetChainID();
         std::string::size_type idx = r.find(this_chain_id);
         if (idx != std::string::npos) {
            r.erase(idx,1);
         }
      }
      if (r.length()) {
         s.first = true;
         s.second = r.substr(0,1);
      }
   } else {
      s.first = true;
      s.second = "A";
   }
   return s;
}


void
coot::molecule_t::remove_TER_on_last_residue(mmdb::Chain *chain_p) {

   int n_residues = chain_p->GetNumberOfResidues();
   if (n_residues > 0) {
      mmdb::Residue *r = chain_p->GetResidue(n_residues-1); // last residue
      if (r)
         remove_TER_internal(r);
   }
}

// remove TER record from residue
//
void
coot::molecule_t::remove_TER_internal(mmdb::Residue *res_p) {

   int n_residue_atoms;
   mmdb::PPAtom residue_atoms;
   bool deleted = 0;
   if (res_p) {
      res_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
         if (residue_atoms[i]->isTer()) {
            res_p->DeleteAtom(i);
            deleted = 1;
         }
      }
   }
   if (deleted) {
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
   }
}


int
coot::molecule_t::insert_waters_into_molecule(const coot::minimol::molecule &water_mol, const std::string &res_name) {

   int istat = 0;  // set to failure initially

   // So run over the the chains of the existing molecule looking for
   // a solvent chain.  If there isn't one we simply use
   // append_to_molecule()
   //
   int nchains = atom_sel.mol->GetNumberOfChains(1);
   mmdb::Chain *chain_p = NULL;
   mmdb::Chain *solvent_chain_p = NULL;
   short int i_have_solvent_chain_flag = 0;
   for (int ichain=0; ichain<nchains; ichain++) {

      chain_p = atom_sel.mol->GetChain(1,ichain);
      if (chain_p->isSolventChain()) {
         solvent_chain_p = chain_p;
         std::string mol_chain_id(chain_p->GetChainID());
         i_have_solvent_chain_flag = 1;
      }
   }


   // For every atom in water_mol, create a new atom and a new residue
   // for it. Add the residue to our model's solvent chain and the
   // atom the the residue (of course).
   //
   if (i_have_solvent_chain_flag == 0) {

      // We didn't manage to find a solvent chain.
      // We need to create a new chain.
      chain_p = new mmdb::Chain;
      atom_sel.mol->GetModel(1)->AddChain(chain_p);
      std::pair<bool, std::string> u = unused_chain_id();
      if (u.first)
         chain_p->SetChainID(u.second.c_str());
      else
         chain_p->SetChainID("Z");
   } else {
      chain_p = solvent_chain_p; // put it back, (kludgey, should use
                                 // solvent_chain_p from here, not chain_p).
      // OK, we also need to remove any TER cards that are in that chain_p
      remove_TER_on_last_residue(solvent_chain_p);
   }

//    std::cout << "Debug:: choose chain " << chain_p->GetChainID()
//              << " with have_solvent flag: " << i_have_solvent_chain_flag
//              << std::endl;
//    std::cout << "Debug:: isSolvent for each residue of chain: " << std::endl;
//    for (int tmp_r=0; tmp_r<chain_p->GetNumberOfResidues(); tmp_r++) {
//       mmdb::Residue *rtmp = chain_p->GetResidue(tmp_r);
//       short int flag = isSolvent(rtmp->name);
//          std::cout << rtmp->name << " is solvent? " << flag << std::endl;
//    }

   std::pair<short int, int> p = coot::util::max_resno_in_chain(chain_p);
   // float bf = graphics_info_t::default_new_atoms_b_factor; // 20.0 by default
   float bf = 20.0;
   int max_resno;
   if (p.first) {
      max_resno = p.second;
   } else {
      max_resno = 0;
   }
   if (p.first || (i_have_solvent_chain_flag == 0)) {
      make_backup("insert_waters_into_molecule");
      std::cout << "INFO:: Adding to solvent chain: " << chain_p->GetChainID() << std::endl;
      atom_sel.delete_atom_selection();
      int prev_max_resno = max_resno;
      mmdb::Residue *new_residue_p = NULL;
      mmdb::Atom    *new_atom_p = NULL;
      int water_count = 0;
      float occ = 1.0;
      if (is_from_shelx_ins_flag)
         occ = 11.0;
      for (unsigned int ifrag=0; ifrag<water_mol.fragments.size(); ifrag++) {
         for (int ires=water_mol[ifrag].min_res_no();
              ires<=water_mol[ifrag].max_residue_number();
              ires++) {
            for (unsigned int iatom=0; iatom<water_mol[ifrag][ires].atoms.size(); iatom++) {
               new_residue_p = new mmdb::Residue;
               new_residue_p->SetResName(res_name.c_str());
               new_residue_p->seqNum = prev_max_resno + 1 + water_count;
               water_count++;
               bf = water_mol[ifrag][ires][iatom].temperature_factor;
               new_atom_p = new mmdb::Atom;
               new_atom_p->SetCoordinates(water_mol[ifrag][ires][iatom].pos.x(),
                                          water_mol[ifrag][ires][iatom].pos.y(),
                                          water_mol[ifrag][ires][iatom].pos.z(), occ, bf);
               if (false)
                  std::cout << "debug:: add water " << ires << " " << water_mol[ifrag][ires][iatom].pos.format()
                            << " with b " << bf << std::endl;
               new_atom_p->SetAtomName(water_mol[ifrag][ires][iatom].name.c_str());
               new_atom_p->Het = 1; // waters are now HETATMs
               strncpy(new_atom_p->element, water_mol[ifrag][ires][iatom].element.c_str(), 3);
               strncpy(new_atom_p->altLoc, water_mol[ifrag][ires][iatom].altLoc.c_str(), 2);

               // residue number, atom name, occ, coords, b factor

               // add the atom to the residue and the residue to the chain
               new_residue_p->AddAtom(new_atom_p);
               chain_p->AddResidue(new_residue_p);
            }
         }
      }
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      coot::util::pdbcleanup_serial_residue_numbers(atom_sel.mol);
      atom_sel = make_asc(atom_sel.mol);
      // update_molecule_after_additions(); // sets unsaved changes flag
      update_symmetry();
      // save_info.new_modification("insert_waters_into_molecule");
   }
   return istat;
}


int
coot::molecule_t::append_to_molecule(const coot::minimol::molecule &water_mol) {

   int istat = 0; // fail status initially.
   int n_atom = 0;  // 0 new atoms added initially.
   float default_new_atoms_b_factor = 20.0;

   if (atom_sel.n_selected_atoms > 0) {

      make_backup("append_to_molecule(water_mol)");

      // run over the chains in water_mol (there is only one for waters)
      //
      for (unsigned int ifrag=0; ifrag<water_mol.fragments.size(); ifrag++) {

//          std::cout << "DEBUG:: append_to_molecule: fragment id id for frag " << ifrag
//                    << " is " << water_mol[ifrag].fragment_id << std::endl;

         short int imatch = 0;

         // Run over chains of the existing mol, to see if there
         // already exists a chain with the same chain id as the
         // waters we want to add.  Only if imatch is 0 does this
         // function do anything.
         //
         int nchains = atom_sel.mol->GetNumberOfChains(1);
         mmdb::Chain *chain;
         for (int ichain=0; ichain<nchains; ichain++) {

            chain = atom_sel.mol->GetChain(1,ichain);
            std::string mol_chain_id(chain->GetChainID());

            if (water_mol.fragments[ifrag].fragment_id == mol_chain_id) {
               //
               imatch = 1;
               istat = 1;
               std::cout << "INFO:: Can't add waters from additional molecule "
                         << "chain id = " << mol_chain_id << std::endl
                         << "INFO:: That chain id already exists in this molecule"
                         << std::endl;
               break;
            }
         }

         mmdb::Model *model_p = atom_sel.mol->GetModel(1);
         if (imatch == 0) {
            // There was not already a chain in this molecule of that name.

            mmdb::Chain *new_chain_p;
            mmdb::Atom *new_atom_p;
            mmdb::Residue *new_residue_p;

            new_chain_p = new mmdb::Chain;
            std::cout << "DEBUG INFO:: chain id of new chain :"
                      << water_mol[ifrag].fragment_id << ":" << std::endl;
            new_chain_p->SetChainID(water_mol[ifrag].fragment_id.c_str());
            model_p->AddChain(new_chain_p);

            for (int ires=water_mol[ifrag].min_res_no();
                 ires<=water_mol[ifrag].max_residue_number();
                 ires++) {

               if (water_mol[ifrag][ires].atoms.size() > 0) {
                  new_residue_p = new mmdb::Residue;
                  new_residue_p->seqNum = ires;
                  strcpy(new_residue_p->name, water_mol[ifrag][ires].name.c_str());
                  new_chain_p->AddResidue(new_residue_p);
                  for (unsigned int iatom=0; iatom<water_mol[ifrag][ires].atoms.size(); iatom++) {

                     new_atom_p = new mmdb::Atom;
                     new_atom_p->SetAtomName(water_mol[ifrag][ires][iatom].name.c_str());
                     new_atom_p->SetElementName(water_mol[ifrag][ires][iatom].element.c_str());
                     new_atom_p->SetCoordinates(water_mol[ifrag][ires][iatom].pos.x(),
                                                water_mol[ifrag][ires][iatom].pos.y(),
                                                water_mol[ifrag][ires][iatom].pos.z(),
                                                1.0, default_new_atoms_b_factor);
                     new_residue_p->AddAtom(new_atom_p);
                     n_atom++;
                  }
               }
            }
         }
      }

      std::cout << "INFO:: " << n_atom << " atoms added to molecule." << std::endl;
      if (n_atom > 0) {
         atom_sel.mol->FinishStructEdit();
         // update_molecule_after_additions(); // sets unsaved changes flag
         // save_info.new_modification(__FUNCTION__);
      }
   }

   return istat;
}


coot::omega_distortion_info_container_t
coot::molecule_t::peptide_omega_analysis(const protein_geometry &geom, const std::string &chain_id,
                                         bool mark_cis_peptides_as_bad_flag) const {


   restraints_container_t rc(atom_sel, chain_id, nullptr);
   omega_distortion_info_container_t odi = rc.omega_trans_distortions(geom, mark_cis_peptides_as_bad_flag);
   return odi;
}


int
coot::molecule_t::cis_trans_conversion(const std::string &atom_cid, mmdb::Manager *standard_residues_mol) {

   if (! is_valid_model_molecule()) return 0;

   int status = 0;
   bool is_N_flag = false;
   mmdb::Atom *at = cid_to_atom(atom_cid);
   std::string atom_name(at->GetAtomName());
   if (atom_name == " N  ") is_N_flag = true;
   if (at) {
      status = coot::util::cis_trans_conversion(at, is_N_flag, atom_sel.mol, standard_residues_mol);
   }
   return status;
}


int
coot::molecule_t::add_compound(const coot::dictionary_residue_restraints_t &restraints,
                               const coot::Cartesian &position, const clipper::Xmap<float> &xmap,
                               float map_rmsd) {

#if 0 // this is in molecules_container_modelling because of merge_molecules().

   int status = 0;

   auto move_residue = [] (mmdb::Residue *residue_p, const coot::Cartesian &position) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      coot::Cartesian position_sum;
      unsigned int n_atoms = 0;
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            n_atoms++;
            position_sum + coot::Cartesian(at->x, at->y, at->z);
         }
      }
      if (n_atoms > 0) {
         float inv = 1.0 / static_cast<float>(n_atoms);
         coot::Cartesian current_position(position_sum * inv);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               at->x += position.x() - current_position.x();
               at->y += position.y() - current_position.y();
               at->z += position.z() - current_position.z();
            }
         }
      }
   };

   bool idealised_flag = true;
   float b_factor = coot::util::median_temperature_factor(atom_sel.atom_selection,
                                                          atom_sel.n_selected_atoms,
                                                          0, 100, false, false);

   mmdb::Residue *residue_p = restraints.GetResidue(idealised_flag, b_factor);
   if (residue_p) {
      move_residue(residue_p, position);
      mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(residue_p);
      if (mol) {
         atom_selection_container_t asc = make_asc(mol);
         std::vector<atom_selection_container_t> ascs = { asc };
         auto merge_results = merge_molecules(ascs);
         if (merge_results.first == 1) {
            if (merge_results.second.size() == 1) {
               coot::residue_spec_t res_spec = merge_results.second[0].spec;
               int n_trials = 100;
               float translation_scale_factor = 1.0;
               float d = fit_to_map_by_random_jiggle(res_spec, xmap, map_rmsd, n_trials, translation_scale_factor);
               std::cout << "d: " << d << std::endl;
            }
         }
      }
   }
   return status;
#endif

   return 0;
}


#include "coot-utils/reduce.hh"

int
coot::molecule_t::add_hydrogen_atoms(protein_geometry *geom_p) {

   make_backup("add_hydrogen_atoms");
   atom_sel.delete_atom_selection();
   coot::reduce r(atom_sel.mol, imol_no);
   r.add_geometry(geom_p);
   bool go_nuclear_flag = false;
   r.add_hydrogen_atoms(go_nuclear_flag);
   coot::util::pdbcleanup_serial_residue_numbers(atom_sel.mol);
   atom_sel = make_asc(atom_sel.mol); // it would be better if there was a member function to do this.

   // save_info.new_modification("add-hydrogen-atoms");
   return 1;

}

int
coot::molecule_t::delete_hydrogen_atoms() {

   make_backup("delete_hydrogen_atoms");
   atom_sel.delete_atom_selection();
   coot::reduce r(atom_sel.mol, imol_no);
   r.delete_hydrogen_atoms();
   atom_sel = make_asc(atom_sel.mol); // is there an asc function for this?
   // save_info.new_modification("add-hydrogen-atoms");
   return 1;

}


std::vector<coot::residue_spec_t>
coot::molecule_t::get_non_standard_residues_in_molecule() const {

   std::vector<coot::residue_spec_t> rv;
   std::vector<std::string> v = util::non_standard_residue_types_in_molecule(atom_sel.mol);

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<n_res; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               std::string rn(residue_p->GetResName());
               const std::vector<std::string>::const_iterator it = std::find(v.begin(), v.end(), rn);
               if (it != v.end()) {
                  if (rn == "HOH") {
                     // bypass
                  } else {
                     coot::residue_spec_t spec(residue_p);
                     spec.string_user_data = rn;
                     rv.push_back(spec);
                  }
               }
            }
         }
      }
   }
   return rv;
}

std::vector<std::string>
coot::molecule_t::get_residue_types_without_dictionaries(const protein_geometry &geom) const {

   std::vector<std::string> v;
   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<n_res; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               std::string rn(residue_p->GetResName());
               if (geom.have_dictionary_for_residue_type_no_dynamic_add(rn, imol_no)) {
                  // don't add
               } else {
                  v.push_back(rn);
               }
            }
         }
      }
   }
   return v;
}



std::vector<std::string>
coot::molecule_t::get_chain_ids() const {

   std::vector<std::string> chain_ids;
   mmdb::Manager *mol = atom_sel.mol;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         if (n_res > 0) {
            chain_ids.push_back(std::string(chain_p->GetChainID()));
         }
      }
   }
   return chain_ids;
}

//! Get the chains that are related by NCS:
std::vector<std::vector<std::string> >
coot::molecule_t::get_ncs_related_chains() const {

   std::vector<std::vector<std::string> > v;
   int model_number = 1;
   std::vector<std::vector<mmdb::Chain *> > ncs_related_chains = coot::ncs_related_chains(atom_sel.mol, model_number);
   std::cout << "found ncs_related_chains size " << ncs_related_chains.size() << std::endl;
   for (const auto &vv : ncs_related_chains) {
      std::cout << "vv size " << vv.size() << std::endl;
      std::vector<std::string> vc;
      for (const auto &c : vv) {
         std::string chid = c->GetChainID();
         std::cout << " " << chid;
         vc.push_back(chid);
      }
      std::cout << std::endl;
      v.push_back(vc);
   }
   return v;
}



#include "density-contour/gaussian-surface.hh"

// 20230206-PE maybe pass the colour map later - std::vector<std::string, std::string>
// which would override the built-in colour rules
//
coot::simple_mesh_t
coot::molecule_t::get_gaussian_surface(float sigma, float contour_level,
                                       float box_radius, float grid_scale, float b_factor) const {

   auto colour_holder_to_glm = [] (const coot::colour_holder &ch) {
      return glm::vec4(ch.red, ch.green, ch.blue, ch.alpha);
   };

   coot::simple_mesh_t mesh;

   if (is_valid_model_molecule()) {

      std::vector<std::string> chain_ids = get_chain_ids();
      mmdb::Manager *mol = atom_sel.mol;

      for (unsigned int i_ch=0; i_ch<chain_ids.size(); i_ch++) {
         const auto &chain_id = chain_ids[i_ch];
         coot::gaussian_surface_t gauss_surf(mol, chain_id, sigma, contour_level, box_radius, grid_scale, b_factor);
         coot::simple_mesh_t gs_mesh = gauss_surf.get_surface();

         // do we have a colour rule to change the colour of that surface?
         //
         // This is hacky because in general colour rules are for selection (even down to the atom)
         // but here we are interested only in colour rules that apply to just a chain
         // So... I am looking only for colour rules that are for this chain
         //
         for (unsigned int icr=0; icr<colour_rules.size(); icr++) {
            const std::string &colour_rule_cid = colour_rules[icr].first;
            const std::string &colour          = colour_rules[icr].second;
            if (std::string("//" + chain_id) == colour_rule_cid ||
                colour_rule_cid == "/"    ||
                colour_rule_cid == "//"   ||
                colour_rule_cid == "/*/*/*/*") {
               coot::colour_holder ch(colour);
               glm::vec4 col = colour_holder_to_glm(ch);
               gs_mesh.change_colour(col);
            }
         }
         mesh.add_submesh(gs_mesh);
      }
   }

   return mesh;
}


coot::simple_mesh_t
coot::molecule_t::get_chemical_features_mesh(const std::string &cid,
                                             const coot::protein_geometry &geom) const {

   coot::simple_mesh_t mesh;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   mmdb::Residue *residue_p = cid_to_residue(cid);
   if (residue_p) {
      std::vector<simple_mesh_t> meshes = chemical_features::generate_meshes(imol_no, residue_p, geom);
      for (const auto &cf_mesh : meshes)
         mesh.add_submesh(cf_mesh);
   }
#endif

   return mesh;
}


//! add an alternative conformation for the specified residue
int
coot::molecule_t::add_alternative_conformation(const std::string &cid) {

   int status = 0;

   auto move = [] (mmdb::Atom *at, const clipper::Coord_orth &o) {
      at->x += o.x();
      at->y += o.y();
      at->z += o.z();
   };

   auto set_offset = [] (clipper::Coord_orth &offset, mmdb::Residue *residue_p) {

      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::vector<mmdb::Atom *> new_atoms;
      mmdb::Atom *c_at = 0;
      mmdb::Atom *n_at = 0;
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         std::string atom_name(at->GetAtomName());
         if (atom_name == " N  ") n_at = at;
         if (atom_name == " C  ") c_at = at;
      }
      if (c_at && n_at) {
         clipper::Coord_orth c_pos(c_at->x, c_at->y, c_at->z);
         clipper::Coord_orth n_pos(n_at->x, n_at->y, n_at->z);
         clipper::Coord_orth cn_unit((c_pos-n_pos).unit());
         clipper::Coord_orth arb(1,2,3);
         clipper::Coord_orth arb_uv(arb.unit());
         clipper::Coord_orth cp = clipper::Coord_orth(clipper::Coord_orth::cross(cn_unit, arb_uv));
         offset = 0.2 * cp;
      }
   };

   mmdb::Residue *residue_p = cid_to_residue(cid);
   if (residue_p) {
      atom_sel.delete_atom_selection();
      clipper::Coord_orth offset(0.0, 0.0, 0.2);
      set_offset(offset, residue_p); // modify reference
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::vector<mmdb::Atom *> new_atoms;
      make_backup("add_alternative_conformation");
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {

            // 20230116-PE this is the trivial/basic/common alt conf split method.
            //             If the residue already has an alt-conf then things will
            //             be more complex. Let's just use this for now - it will
            //             handle the vast majority of cases.
            //
            std::string current_alt_conf(at->altLoc);
            if (current_alt_conf.empty()) {
               mmdb::Atom *at_new = new mmdb::Atom;
               at_new->Copy(at);
               move(at_new, -offset);
               strcpy(at_new->altLoc, "B");
               at_new->occupancy = 0.5;
               new_atoms.push_back(at_new);

               move(at, offset);
               at->occupancy = 0.5;
               strcpy(at->altLoc, "A");
            }
         }
      }

      for (unsigned int j=0; j<new_atoms.size(); j++)
         residue_p->AddAtom(new_atoms[j]);

      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      // save_info.new_modification("add-alt-conf");

   } else {
      std::cout << "WARNING:: add_alternative_conformation() Residue " << cid << " not found in molecule"
                << std::endl;
   }

   return status;

}


//! add atoms to a partially-filled side chain
int
coot::molecule_t::fill_partial_residue(const residue_spec_t &res_spec, const std::string &alt_conf,
                                       const clipper::Xmap<float> &xmap, const coot::protein_geometry &geom) {

   int status = 0;

   // these functions do their own make_backup().

   mmdb::Residue *residue_p = get_residue(res_spec);
   if (residue_p) {
      std::string residue_type = residue_p->GetResName();
      int success_a = mutate(res_spec, residue_type); // fill missing atoms
      if (success_a) {
         int success_b = auto_fit_rotamer(res_spec.chain_id, res_spec.res_no, res_spec.ins_code, alt_conf, xmap, geom);
         if (success_b)
            status = 1;
      }
   }
   return status;
}


#include "coot-utils/coot-coord-extras.hh"

//! add atoms to a partially-filled side chain
int
coot::molecule_t::fill_partial_residues(const clipper::Xmap<float> &xmap, protein_geometry *geom_p) {

   int status = 0;
   bool do_missing_hydrogen_atoms_flag = false;
   coot::util::missing_atom_info mai = coot::util::missing_atoms(atom_sel.mol, do_missing_hydrogen_atoms_flag, geom_p);

   if (! mai.residues_with_missing_atoms.empty()) {
      for (unsigned int i=0; i<mai.residues_with_missing_atoms.size(); i++) {
         mmdb::Residue *residue_p = mai.residues_with_missing_atoms[i];
         int res_no =  mai.residues_with_missing_atoms[i]->GetSeqNum();
         std::string chain_id = mai.residues_with_missing_atoms[i]->GetChainID();
         std::string residue_type = mai.residues_with_missing_atoms[i]->GetResName();
         std::string inscode = mai.residues_with_missing_atoms[i]->GetInsCode();
         std::string altloc("");
         coot::residue_spec_t res_spec(residue_p);
         mutate(res_spec, residue_type); // fill missing atoms, does make_backup()
         int success_b = auto_fit_rotamer(chain_id, res_no, inscode, altloc, xmap, *geom_p);
         if (success_b)
            status = 1;
      }
   }
   return status;
}


#include "ideal/add-linked-cho.hh"
void
coot::molecule_t::add_named_glyco_tree(const std::string &glycosylation_name, const std::string &chain_id, int res_no,
                                       const clipper::Xmap<float> &xmap,
                                       coot::protein_geometry *geom) {

   // the atom selection gets updated.
   float mt = get_median_temperature_factor();
   float new_atoms_b_factor = 1.55 * mt;
   coot::cho::add_named_glyco_tree(glycosylation_name, &atom_sel, imol_no,
                                   new_atoms_b_factor,
                                   xmap, geom, chain_id, res_no);

}



// --------------- rigid body fit
#include "rigid-body-fit.hh"
int
coot::molecule_t::rigid_body_fit(const std::string &multi_cids, const clipper::Xmap<float> &xmap) {

   int status = 0;

   std::vector<std::string> v = coot::util::split_string(multi_cids, "||");

   if (! v.empty()) {

      // udd_atom_selection is (just) the selection for the moving atoms.
      // we make it and delete it. the atom_sel selection is not touched.
      int udd_atom_selection = atom_sel.mol->NewSelection(); // d

      for (const auto &cid : v) {
         atom_sel.mol->Select(udd_atom_selection, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_OR);
         if (true) { // debugging the selection
            mmdb::PAtom *atoms = NULL;
            int n_atoms;
            atom_sel.mol->GetSelIndex(udd_atom_selection, atoms, n_atoms);
            std::cout << "----------- debug:: in rigid_body_fit() we selected " << n_atoms << " atoms " << std::endl;
            std::cout << "----------- after selection " << cid << " n_atoms " << n_atoms << std::endl;
         }
      }

      make_backup("rigid_body_fit " + multi_cids);
      // update the atoms of atom-sel.mol
      coot::api::rigid_body_fit(atom_sel.mol, udd_atom_selection, xmap);
      status = 1;
      atom_sel.mol->DeleteSelection(udd_atom_selection);
   }
   // save_info.new_modification("rigid-body-fit " + multi_cids);
   return status;
}


// ----------------------------------- symmetry ----------------------------------

coot::symmetry_info_t
coot::molecule_t::get_symmetry(float symmetry_search_radius, const coot::Cartesian &rotation_centre) const {

   int symmetry_shift_search_size = 2; // is 1 in graphics-info-statics
   molecule_extents_t extents(atom_sel, symmetry_search_radius);
   if (true) // turn this off later.
      std::cout << "extents: " << extents << std::endl;
   std::vector<std::pair<symm_trans_t, Cell_Translation> > symm_trans_boxes =
      extents.which_boxes(rotation_centre, atom_sel, symmetry_shift_search_size);
   mmdb::realtype a[6];
   mmdb::realtype vol;
   int orthcode;
   atom_sel.mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
   Cell cell(a[0], a[1], a[2], clipper::Util::d2rad(a[3]), clipper::Util::d2rad(a[4]), clipper::Util::d2rad(a[5]));
   return symmetry_info_t(symm_trans_boxes, cell);
}


// should the be in coot_molecule_rotamers?

#include "ligand/richardson-rotamer.hh"

//
//! change rotamers
coot::molecule_t::rotamer_change_info_t
coot::molecule_t::change_to_next_rotamer(const coot::residue_spec_t &res_spec, const std::string &alt_conf,
                                         const coot::protein_geometry &pg) {

   auto crni = change_rotamer_number(res_spec, alt_conf, 1, pg);
   return crni;
}

//
//! change rotamers
coot::molecule_t::rotamer_change_info_t
coot::molecule_t::change_to_previous_rotamer(const coot::residue_spec_t &res_spec, const std::string &alt_conf,
                                         const coot::protein_geometry &pg) {
   return change_rotamer_number(res_spec, alt_conf, -1, pg);
}
//
//! change rotamers
coot::molecule_t::rotamer_change_info_t
coot::molecule_t::change_to_first_rotamer(const coot::residue_spec_t &res_spec, const std::string &alt_conf,
                                         const coot::protein_geometry &pg) {
   return change_rotamer_number(res_spec, alt_conf, 0, pg);
}

//
//! change rotamers
coot::molecule_t::rotamer_change_info_t
coot::molecule_t::change_rotamer_number(const coot::residue_spec_t &res_spec, const std::string &alt_conf,
                                           int rotamer_change_direction,
                                           const coot::protein_geometry &pg) {

   coot::molecule_t::rotamer_change_info_t rci; // initially "fail" values
   int i_done = 0;
   mmdb::Residue *residue_p = util::get_residue(res_spec, atom_sel.mol);
   if (residue_p) {
      int current_rotamer = -1; // not found
      std::map<residue_spec_t, int>::const_iterator it = current_rotamer_map.find(res_spec);
      if (it != current_rotamer_map.end()) {
         current_rotamer = it->second;
      }

      coot::richardson_rotamer d(residue_p, alt_conf, atom_sel.mol, 0.01, 0);
      std::string res_type = residue_p->GetResName();
      if (res_type == "GLY") return rci;
      if (res_type == "ALA") return rci;
      std::pair<short int, coot::dictionary_residue_restraints_t> p =
         pg.get_monomer_restraints(res_type, imol_no);

      if (p.first) {
         float prob_cut = 0.01;
         std::vector<simple_rotamer> rotamers = d.get_rotamers(res_type, prob_cut);
         int rotamer_number = 0;
         if (rotamer_change_direction == 1)
            rotamer_number = current_rotamer + 1;
         if (rotamer_change_direction == -1)
            rotamer_number = current_rotamer - 1;
         int n_rotamers = rotamers.size();
         if (rotamer_number >= n_rotamers)
            rotamer_number = 0;
         if (rotamer_number < 0)
            rotamer_number = n_rotamers - 1;
         if (rotamer_number < 0) // protection
            return rci;

         coot::dictionary_residue_restraints_t rest = p.second;
         mmdb::Residue *moving_res = d.GetResidue(rest, rotamer_number);
         if (moving_res) {
            make_backup("change_rotamer_number()");
            i_done = set_residue_to_rotamer_move_atoms(residue_p, moving_res);
            delete moving_res; // or moving_res->chain?
            current_rotamer_map[res_spec] = rotamer_number;
            rci.status = i_done;
            rci.rank = rotamer_number;
            rci.name = rotamers[rotamer_number].rotamer_name();
            rci.richardson_probability = rotamers[rotamer_number].Probability_rich();

         }
      } else {
         std::cout << "WARNING:: change_rotamer_number() Failed to get monomer restraints for " << res_type << std::endl;
      }
   } else {
      std::cout << "WARNING:: change_rotamer_number no residue found" << res_spec << std::endl;
   }
   if (! i_done)
      std::cout << "WARNING:: change_rotamer_number(): set rotamer number failed" << std::endl;

   // save_info.new_modification("change_rotamer_number()");
   return rci;

}


int
coot::molecule_t::set_residue_to_rotamer_move_atoms(mmdb::Residue *res, mmdb::Residue *moving_res) {

   // note: calling function makes the backup

   int i_done = 0;
   int n_ref_atoms;
   mmdb::PPAtom ref_residue_atoms = 0;
   int n_mov_atoms;
   mmdb::PPAtom mov_residue_atoms= 0;

   res->GetAtomTable(ref_residue_atoms, n_ref_atoms);
   moving_res->GetAtomTable(mov_residue_atoms, n_mov_atoms);

   int n_atoms = 0;
   for (int imov=0; imov<n_mov_atoms; imov++) {
      std::string atom_name_mov(mov_residue_atoms[imov]->name);
      std::string alt_loc_mov(mov_residue_atoms[imov]->altLoc);
      for (int iref=0; iref<n_ref_atoms; iref++) {
         std::string atom_name_ref(ref_residue_atoms[iref]->name);
         std::string alt_loc_ref(ref_residue_atoms[iref]->altLoc);
         if (atom_name_mov == atom_name_ref) {
            if (alt_loc_mov == alt_loc_ref) {
               ref_residue_atoms[iref]->x = mov_residue_atoms[imov]->x;
               ref_residue_atoms[iref]->y = mov_residue_atoms[imov]->y;
               ref_residue_atoms[iref]->z = mov_residue_atoms[imov]->z;
               n_atoms++;
               i_done = 1;
            }
         }
      }
   }

   // 20240516-PE I aom only moving atoms - so I don't need any of this:
   // (it destroys atom_sel and atom_sel.UDDOldAtomIndexHandle)
   //
   // if (i_done) {
      // atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      // atom_sel.mol->FinishStructEdit();
      // atom_sel = make_asc(atom_sel.mol);
   // }
   return i_done;
}


// add or update.
void
coot::molecule_t::add_target_position_restraint(const std::string &atom_cid, float pos_x, float pos_y, float pos_z) {

   // make this a class member - and clear it when refinement starts.


   mmdb::Atom *at = cid_to_atom(atom_cid);
   if (at) {
      // try to find it...
      bool done = false;
      for (unsigned int i=0; i<atoms_with_position_restraints.size(); i++) {
         if (atoms_with_position_restraints[i].first == at) {
            clipper::Coord_orth p(pos_x, pos_y, pos_z);
            atoms_with_position_restraints[i].second = p;
            done = true;
         }
      }
      if (!done) {
         clipper::Coord_orth p(pos_x, pos_y, pos_z);
         auto pp = std::make_pair(at, p);
         atoms_with_position_restraints.push_back(pp);
      }
   }
}


// add or update.
coot::instanced_mesh_t
coot::molecule_t::add_target_position_restraint_and_refine(const std::string &atom_cid, float pos_x, float pos_y, float pos_z,
                                                           int n_cycles, coot::protein_geometry *geom_p) {

   coot::instanced_mesh_t m;
   add_target_position_restraint(atom_cid, pos_x, pos_y, pos_z);

   for (unsigned int i=0; i<atoms_with_position_restraints.size(); i++) {
      const auto &pp = atoms_with_position_restraints[i];
      clipper::Coord_orth p = pp.second;
      mmdb::Atom *at = pp.first;
      at->x = p.x(); at->y = p.y(); at->z = p.z();
   }

   if (n_cycles < 0) {
      // simple communication test.
   } else {
      // actually do some refinement then
      if (last_restraints) {
         clipper::Coord_orth pos(pos_x, pos_y, pos_z);
         mmdb::Atom *at = cid_to_atom(atom_cid);
         if (at) {
            coot::atom_spec_t spec(at);
            last_restraints->add_atom_pull_restraint(spec, pos);
            std::cout << "debug:: in wrapped_add_target_position_restraint() calling refine_using_last_restraints() "
                      << n_cycles << " cycles " << std::endl;
            refine_using_last_restraints(n_cycles);
         } else {
            std::cout << "wrapped_add_target_position_restraint() failed to find atom given " << atom_cid << std::endl;
         }
      } else {
         std::cout << "DEBUG:: in wrapped_add_target_position_restraint() last_restraints was empty! " << std::endl;
      }
   }

   std::string mode = "COLOUR-BY-CHAIN-AND-DICTIONARY";
   unsigned int smoothness_factor = 1;
   bool show_atoms_as_aniso_flag = false;
   bool show_aniso_atoms_as_ortep_flag = false;
   float aniso_probability = 0.5f;
   m = get_bonds_mesh_instanced(mode, geom_p, true, 0.1, 1.4,
                                show_atoms_as_aniso_flag,
                                aniso_probability,
                                show_aniso_atoms_as_ortep_flag,
                                smoothness_factor, true, true);
   return m;

}

//! clear any and all drag-atom target position restraints
void
coot::molecule_t::clear_target_position_restraints() {
   atoms_with_position_restraints.clear();
   // now actually remove them from the refinement...
   if (last_restraints)
      last_restraints->clear_all_atom_pull_restraints();
}

//! clear
void
coot::molecule_t::clear_target_position_restraint(const std::string &atom_cid) {

   mmdb:: Atom *at = cid_to_atom(atom_cid);
   if (at) {
      coot::atom_spec_t spec(at);
      if (last_restraints) {
         last_restraints->clear_atom_pull_restraint(spec);
      }
   }
}

//
void
coot::molecule_t::turn_off_when_close_target_position_restraint() {

   if (last_restraints) {
      last_restraints->turn_off_when_close_target_position_restraint();
   }
}



//! call this after molecule refinement has finished (say when the molecule molecule is accepted into the
//! original molecule
void
coot::molecule_t::clear_refinement() {

   if (last_restraints) {
      std::cout << "debug:: ---------- clear_refinement() ---------- " << std::endl;
      delete last_restraints;
      last_restraints = nullptr; // clear_refinement();
   }
}



void
coot::molecule_t::fix_atom_selection_during_refinement(const std::string &atom_selection_cid) {

   // get atoms in atom selection.
   // match those with thatom in the refinement and mark them as fixed. c.f. fixed (anchored) atom
   int selHnd = atom_sel.mol->NewSelection();
   atom_sel.mol->Select(selHnd, mmdb::STYPE_ATOM, atom_selection_cid.c_str(), mmdb::SKEY_NEW);
   int nSelAtoms = 0;
   mmdb::Atom **SelAtom = nullptr;
   atom_sel.mol->GetSelIndex(selHnd, SelAtom, nSelAtoms);
   for (int i=0; i<nSelAtoms; i++) {
      mmdb::Atom *at = SelAtom[i];
   }
   atom_sel.mol->DeleteSelection(selHnd);
}

// refine all of this molecule - the links and non-bonded contacts will be determined from mol_ref;
void
coot::molecule_t::init_all_molecule_refinement(int imol_ref_mol, coot::protein_geometry &geom,
                                               const clipper::Xmap<float> &xmap_in, float map_weight,
                                               ctpl::thread_pool *thread_pool) {

   // maybe we can use mol_ref as the molecule to call in add_neighbor_residues_for_refinement_help(mmdb::Manager *mol)
   // then we don't need a special version of copy_fragment(). Or so it currently seems to me.

   bool make_trans_peptide_restraints = true;
   bool do_rama_plot_restraints = false;
   bool refinement_is_quiet = true;

   auto get_all_residues_in_molecule = [] (mmdb::Manager *mol) {
      std::vector<mmdb::Residue *> rv;
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  rv.push_back(residue_p);
               }
            }
         }
      }
      return rv;
   };

   std::vector<mmdb::Residue *> rv = get_all_residues_in_molecule(atom_sel.mol);
   std::vector<coot::atom_spec_t> fixed_atom_specs;
   std::vector<std::pair<bool,mmdb::Residue *> > local_residues;
   for (const auto &r : rv)
      local_residues.push_back(std::make_pair(false, r));

   // neighb_residues are from a different mmdb::Manager - will that work in refinement?
   if (! neighbouring_residues.empty())
      local_residues.insert(local_residues.end(), neighbouring_residues.begin(), neighbouring_residues.end());

   make_backup("init_all_molecule_refinement");
   mmdb::Manager *mol = atom_sel.mol;
   std::vector<mmdb::Link> links;
   last_restraints = new coot::restraints_container_t(local_residues, links, geom, mol, fixed_atom_specs, &xmap_in);

   if (refinement_is_quiet)
      last_restraints->set_quiet_reporting();

   // std::cout << "DEBUG:: using restraints with map_weight " << map_weight << std::endl;
   last_restraints->add_map(map_weight);
   coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
   flags = coot::TYPICAL_RESTRAINTS;
   coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;

   unsigned int n_threads = 8;
   last_restraints->thread_pool(thread_pool, n_threads);

   // user-defined LIG diction ahve been assigned to imol_ref_mol, not this one (this one
   // is a temporary molecule used only for refinement).
   last_restraints->make_restraints(imol_ref_mol, geom, flags, 1, make_trans_peptide_restraints,
                                    1.0, do_rama_plot_restraints, true, true, false, pseudos);

   if (last_restraints->size() == 0) {
      // Failure.
      // clear up
      delete last_restraints;
      last_restraints = nullptr; // failure to setup
   }
}

// ------------------------------- rsr utils - add in the environment of this fragment molecule
// from the reidue from which this fragment was copied
void
coot::molecule_t::add_neighbor_residues_for_refinement_help(mmdb::Manager *mol) {

   // cid may not be needed because we want the neighbour of all the residues of this new (fragment) molecule

   auto get_residues_in_selection = [] (mmdb::Manager *mol, int selHnd) {
      std::vector<std::pair<bool,mmdb::Residue *> > residues_vec;
      int nSelResidues = 0;
      mmdb::Residue **SelResidues = 0;
      mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
      if (nSelResidues > 0) {
         residues_vec.resize(nSelResidues);
         for (int i=0; i<nSelResidues; i++) {
            residues_vec.push_back(std::make_pair(true, SelResidues[i]));
         }
      }
      return residues_vec;
   };

   auto map_of_sets_to_residue_vec = [] (const std::map<mmdb::Residue *, std::set<mmdb::Residue *> > &rnr) {

      std::vector<std::pair<bool, mmdb::Residue *> > neighb_residues;
      std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it;

      std::set<mmdb::Residue *> keys; // the residues in the molten zone
      for (it=rnr.begin(); it!=rnr.end(); ++it) {
         const auto &key = it->first;
         keys.insert(key);
      }

      for (it=rnr.begin(); it!=rnr.end(); ++it) {
         std::set<mmdb::Residue *>::const_iterator it_s;
         const auto &s = it->second;
         for (it_s=s.begin(); it_s!=s.end(); ++it_s) {
            mmdb::Residue *r(*it_s);
            // because this is a pair, I can't use find()
            bool found = false;
            for (unsigned int i=0; i<neighb_residues.size(); i++) {
               if (neighb_residues[i].second == r) {
                  found = true;
                  break;
               }
            }
            if (! found) {
               // we don't want environment residue to be any of the residues that are in the molten zone
               if (keys.find(r) == keys.end()) {
                  auto rp = std::pair<bool, mmdb::Residue *> (false, *it_s);
                  neighb_residues.push_back(rp);
               }
            }
         }
      }
      return neighb_residues;
   };

   // now code to save the environment of the residues in the new fragment
   int selHnd_residues = mol->NewSelection(); // d
   mol->Select(selHnd_residues, mmdb::STYPE_RESIDUE, "//", mmdb::SKEY_NEW); // all residues. Maybe "/1/"?
   float dist_crit = 5.0;
   std::vector<std::pair<bool,mmdb::Residue *> > residues_vec = get_residues_in_selection(mol, selHnd_residues);
   std::map<mmdb::Residue *, std::set<mmdb::Residue *> > rnr =
      coot::residues_near_residues(residues_vec, mol, dist_crit);

   neighbouring_residues = map_of_sets_to_residue_vec(rnr);
}

// static
std::string
coot::molecule_t::file_to_string(const std::string &file_name) {

   std::string s;
   std::string line;
   std::ifstream f(file_name.c_str());
   if (!f) {
      std::cout << "Failed to open " << file_name << std::endl;
   } else {
      while (std::getline(f, line)) {
         s += line;
         s += "\n";
      }
   }
   return s;
}


//! @return a model molecule imol as a string. Return emtpy string on error
std::string
coot::molecule_t::molecule_to_PDB_string() const {


   std::string s;

   if (is_valid_model_molecule()) {
      atom_sel.mol->WritePDBASCII("tmp.pdb");
      s = file_to_string("tmp.pdb");
   }

   return s;

}

//! @return a model molecule imol as a string. Return emtpy string on error
std::string
coot::molecule_t::molecule_to_mmCIF_string() const {

   std::string s;
   if (is_valid_model_molecule()) {

      mmdb::Manager *mol_copy = new mmdb::Manager;
      mol_copy->Copy(atom_sel.mol, mmdb::MMDBFCM_All);
      mol_copy->WriteCIFASCII("tmp.cif");
      s = file_to_string("tmp.cif");
      delete mol_copy;
   }
   return s;
}

// ------------------------------ put these functions in coot_molecule_refine.cc --------------


//! @return a list of residues specs that have atoms within dist of the atoms of the specified residue
std::vector<coot::residue_spec_t>
coot::molecule_t::residues_near_residue(const std::string &residue_cid, float dist) const {

   std::vector<coot::residue_spec_t> v;
   mmdb::Residue *residue_p = cid_to_residue(residue_cid);
   if (residue_p) {
      coot::residue_spec_t res_spec_in(residue_p);
      v = coot::residues_near_residue(res_spec_in, atom_sel.mol, dist);
   }
   return v;
}




//! export map molecule as glTF
void
coot::molecule_t::export_map_molecule_as_gltf(clipper::Coord_orth &p, float radius, float contour_level,
                                              const std::string &file_name) {

   coot::simple_mesh_t map_mesh = get_map_contours_mesh(p, radius, contour_level, false, nullptr);
   bool as_binary = true; // test the extension of file_name
   map_mesh.export_to_gltf(file_name, gltf_pbr_roughness, gltf_pbr_metalicity, as_binary);

}

//! export model molecule as glTF - This API will change - we want to specify surfaces and ribbons too.
void
coot::molecule_t::export_model_molecule_as_gltf(const std::string &mode,
                                                const std::string &selection_cid,
                                                coot::protein_geometry *geom,
                                                bool against_a_dark_background,
                                                float bonds_width, float atom_radius_to_bond_width_ratio, int smoothness_factor,
                                                bool draw_hydrogen_atoms_flag, bool draw_missing_residue_loops,
                                                const std::string &file_name) {

   bool show_atoms_as_aniso_flag = true;
   bool show_aniso_atoms_as_ortep_flag = false; // pass these

   instanced_mesh_t im = get_bonds_mesh_for_selection_instanced(mode, selection_cid, geom, against_a_dark_background,
                                                                bonds_width, atom_radius_to_bond_width_ratio,
                                                                show_atoms_as_aniso_flag,
                                                                show_aniso_atoms_as_ortep_flag,
                                                                smoothness_factor,
                                                                draw_hydrogen_atoms_flag, draw_missing_residue_loops);

   coot::simple_mesh_t sm = coot::instanced_mesh_to_simple_mesh(im);
   bool as_binary = true; // test the extension of file_name
   sm.export_to_gltf(file_name, gltf_pbr_roughness, gltf_pbr_metalicity, as_binary);

}

void
coot::molecule_t::export_molecular_representation_as_gltf(const std::string &atom_selection_cid,
                                                         const std::string &colour_scheme, const std::string &style,
                                                         int secondary_structure_usage_flag,
                                                         const std::string &file_name) {

   coot::simple_mesh_t sm = get_molecular_representation_mesh(atom_selection_cid, colour_scheme, style, secondary_structure_usage_flag);
   bool as_binary = true; // test the extension of file_name
   sm.export_to_gltf(file_name, gltf_pbr_roughness, gltf_pbr_metalicity, as_binary);
}

void
coot::molecule_t::export_chemical_features_as_gltf(const std::string &cid,
                                                   const coot::protein_geometry &geom,
                                                   const std::string &file_name) const {

   coot::simple_mesh_t sm = get_chemical_features_mesh(cid, geom);
   bool as_binary = true; // test the extension of file_name
   sm.export_to_gltf(file_name, gltf_pbr_roughness, gltf_pbr_metalicity, as_binary);
}




//! Interactive B-factor refinement (fun).
//! "factor" might typically be say 0.9 or 1.1
void
coot::molecule_t::multiply_residue_temperature_factors(const std::string &cid, float factor) {

   if (atom_sel.mol) {
      int selHnd = atom_sel.mol->NewSelection(); // d
      mmdb::Atom **SelAtoms = nullptr;
      int nSelAtoms = 0;
      atom_sel.mol->Select(selHnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selHnd, SelAtoms, nSelAtoms);
      if (nSelAtoms > 0) {
         for (int i=0; i<nSelAtoms; i++) {
            mmdb:: Atom *at = SelAtoms[i];
            if (! at->isTer()) {
               float new_B = at->tempFactor * factor;
               at->tempFactor = new_B;
            }
         }
      }
      atom_sel.mol->DeleteSelection(selHnd);
   }
}


#include "coot-utils/coot-coord-extras.hh"

// match those of the passed (reference residue (from a different
// molecule, typically).
//
int
coot::molecule_t::match_torsions(mmdb::Residue *res_reference,
                                 const std::vector <coot::dict_torsion_restraint_t> &tr_ref_res,
                                 const coot::protein_geometry &geom) {

   int n_torsions_moved = 0;
   make_backup("match_torsions");

   mmdb::Residue *res_ligand = coot::util::get_first_residue(atom_sel.mol); // this could/should be replaced
                                                                       // by something that allows
                                                                       // any residue in the molecule
                                                                       // to move to match the
                                                                       // reference residue.

   if (res_ligand) { // the local (moving) residue is xxx_ligand
      std::string res_name_ligand(res_ligand->GetResName());
      std::pair<bool, coot::dictionary_residue_restraints_t> ligand_restraints_info =
         geom.get_monomer_restraints(res_name_ligand, imol_no);
      if (ligand_restraints_info.first) {
         std::vector <coot::dict_torsion_restraint_t> tr_ligand =
            geom.get_monomer_torsions_from_geometry(res_name_ligand, imol_no, 0);
         if (tr_ligand.size()) {

            // find the matching torsion between res_ligand and res_reference and then
            // set the torsions of res_ligand to match those of res_reference.
            //
            // moving the res_ligand
            coot::match_torsions mt(res_ligand, res_reference, ligand_restraints_info.second);
            n_torsions_moved = mt.match(tr_ligand, tr_ref_res);
            atom_sel.mol->FinishStructEdit();
         } else {
            std::cout << "WARNING torsion restraints of ligand: size 0" << std::endl;
         }
      } else {
         std::cout << "WARNING ligand_restraints_info.first failed " << std::endl;
      }
   } else {
      std::cout << "WARNING:: null ligand residue (trying to get first) " << std::endl;
   }
   return n_torsions_moved;
}



void
coot::molecule_t::transform_by(const clipper::RTop_orth &rtop, mmdb::Residue *residue_moving) {

   mmdb::Atom **residue_atoms = nullptr;
   int n_residue_atoms = 0;
   residue_moving->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iatom=0; iatom<n_residue_atoms; iatom++) {
      clipper::Coord_orth p(residue_atoms[iatom]->x,
                            residue_atoms[iatom]->y,
                            residue_atoms[iatom]->z);
      clipper::Coord_orth p2 = p.transform(rtop);
      residue_atoms[iatom]->x = p2.x();
      residue_atoms[iatom]->y = p2.y();
      residue_atoms[iatom]->z = p2.z();
   }

   atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();

}


void
coot::molecule_t::transform_by(const clipper::RTop_orth &rtop) {

   for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
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
                        clipper::Coord_orth pos = coot::co(at);
                        clipper::Coord_orth p2 = pos.transform(rtop);
                        at->x = p2.x();
                        at->y = p2.y();
                        at->z = p2.z();
                     }
                  }
               }
            }
         }
      }
   }
}

void
coot::molecule_t::print_secondary_structure_info() const {

   for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
         coot::util::print_secondary_structure_info(model_p);
      }
   }

}

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
//! if the ligand cid specifies more than one residue, only the first is returned.
//! @return nullptr on error or failure to specify a ligand.
RDKit::ROMol *
coot::molecule_t::rdkit_mol(const std::string &ligand_cid) {

   RDKit::ROMol *mol = nullptr;
   mmdb::Residue *residue_p = cid_to_residue(ligand_cid);
   return mol;
}
#endif


//! get the median temperature factor for the model
//! @return a negative number on failure.
float
coot::molecule_t::get_median_temperature_factor() const {

   float b = coot::util::median_temperature_factor(atom_sel.atom_selection, atom_sel.n_selected_atoms, 2.0, 2222.2, false, false);
   return b;

}

float
coot::molecule_t::get_temperature_factor_of_atom(const std::string &atom_cid) const {

   float b = -1.1f;
   mmdb:: Atom *at = cid_to_atom(atom_cid);
   if (at) {
      b = at->tempFactor;
   }
   return b;

}


#include "utils/coot-fasta.hh"

void
coot::molecule_t::associate_sequence_with_molecule(const std::string &chain_id, const std::string &sequence) {

   // input_sequences.push_back(std::make_pair(chain_id, sequence));
   fasta f(chain_id, sequence, fasta::SIMPLE_STRING);
   multi_fasta_seq.add(f);

}

//! try to fit all of the sequences to all of the chains
void
coot::molecule_t::assign_sequence(const clipper::Xmap<float> &xmap, const coot::protein_geometry &geom) {

   auto apply_sequence = [] (const std::vector<mmdb::Residue *> &residues,
                             const std::string &sequence) {
      // caller checks that the lengths match
      for (unsigned int ires=0; ires<residues.size(); ires++) {
         mmdb::Residue *residue_p = residues[ires];
         char letter = sequence[ires];
         std::string new_residue_type = coot::util::single_letter_to_3_letter_code(letter);
         coot::util::mutate(residue_p, new_residue_type);
      }
   };

   auto apply_fasta_multi_to_fragment = [apply_sequence, geom, this] (mmdb::Manager *mol,
                                            const std::string &chain_id,
                                             int resno_start,
                                             int resno_end,
                                             const clipper::Xmap<float> &xmap,
                                             const coot::fasta_multi &fam) {

      std::vector<mmdb::Residue *> residues;
      side_chain_densities scd;
      unsigned int n_sequences = fam.size();
      for (unsigned int idx=0; idx<n_sequences; idx++) {
         const std::string &sequence = fam[idx].sequence;
         const std::string &name = fam[idx].name;
         std::pair<std::string, std::vector<mmdb::Residue *> > a_run_of_residues =
         scd.setup_test_sequence(mol, chain_id, resno_start, resno_end, xmap);
         if (a_run_of_residues.first.empty()) {
            bool print_slider_results = true;
            scd.test_sequence(a_run_of_residues.second, xmap, name, sequence, print_slider_results);
            bool only_probable = false;
            bool print_sequencing_solutions = true;
            coot::side_chain_densities::results_t new_sequence_result = scd.get_result(only_probable, print_sequencing_solutions);
            std::string new_sequence = new_sequence_result.sequence;
            std::cout << "new sequence  " << new_sequence << std::endl;
            int offset = new_sequence_result.offset;
            if (! new_sequence.empty()) {
               int sl = new_sequence.length();
               int resno_count = resno_end - resno_start + 1;
               std::cout << "compare sl " << sl << " resno_count " << resno_count << std::endl;
               if (sl == resno_count) {
                  std::cout << "..... now apply the sequence" << std::endl;
                  residues = a_run_of_residues.second;
                  apply_sequence(residues, new_sequence);
               }
            }
         } else {
            std::cout << "Error when generating a run-of-residues" << std::endl;
            std::cout << " " << a_run_of_residues.first << std::endl;
         }
      }
      return residues;
   };

   coot::side_chain_densities scd;
   mmdb::Manager *mol = atom_sel.mol;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         std::string chain_id(chain_p->GetChainID());
         int nres = 0;
         mmdb::PResidue *residue_table = 0;
         chain_p->GetResidueTable(residue_table, nres);
         if (nres > 10) {
             int idx_end = nres - 1;
             int resno_start = residue_table[0]->GetSeqNum();
             int resno_end   = residue_table[idx_end]->GetSeqNum();

            {
               atom_sel.delete_atom_selection();
               std::vector<mmdb::Residue *> residues =
                  apply_fasta_multi_to_fragment(atom_sel.mol, chain_id, resno_start, resno_end, xmap, multi_fasta_seq);
               atom_sel.regen_atom_selection();
               // backrub rotamer (actully replace_coords()) uses UDDOldAtomIndexHandle. I don't understand why
               if (false)
                  std::cout << "#### after regen_atom_selection() n_selected_atoms " << atom_sel.n_selected_atoms
                            << " and UDDAtomIndexHandle is " << atom_sel.UDDAtomIndexHandle << std::endl;
               atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
               atom_sel.mol->FinishStructEdit();
               util::pdbcleanup_serial_residue_numbers(atom_sel.mol);

               bool debug_atom_indexing = false;
               if (debug_atom_indexing) {
                  for (int i = 0; i < atom_sel.n_selected_atoms; i++) {
                     int idx = -1;
                     mmdb::Atom *at = atom_sel.atom_selection[i];
                     if (at->GetUDData(atom_sel.UDDAtomIndexHandle, idx) == mmdb::UDDATA_Ok) {
                        std::cout << "OK " << i << " " << idx << std::endl;
                     } else {
                        std::cout << "udd lookup failure for i " << i << std::endl;
                     }
                  }
               }

               for (unsigned int ires=0; ires<residues.size(); ires++) {
                  if (false)
                     std::cout << "#### after regen_atom_selection()"
                              << " mol " << atom_sel.mol
                              << " n_selected_atoms " << atom_sel.n_selected_atoms
                              << " atom_selection " << atom_sel.atom_selection
                              << " and UDDOldAtomIndexHandle is " << atom_sel.UDDOldAtomIndexHandle << std::endl;
                  mmdb::Residue *residue_p = residues[ires];
                  this->backrub_rotamer(residue_p, xmap, geom);
               }
            }

         } else {
            std::cout << "Chain must have at least 10 residue" << std::endl;
         }
      }
   }
   write_coordinates("test-add-sc.pdb");
}




//! get the residue CA position
//!
//! @return a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values
std::vector<double>
coot::molecule_t::get_residue_CA_position(const std::string &cid) const {

   std::vector<double> v;
   mmdb::Residue *residue_p = cid_to_residue(cid);
   if (residue_p) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string name = at->GetAtomName();
            if (name == " CA ") {
               v.push_back(at->x);
               v.push_back(at->y);
               v.push_back(at->z);
               break;
            }
         }
      }
   }
   return v;
}

//! get the avarge residue position
//!
//! @return a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values
std::vector<double>
coot::molecule_t::get_residue_average_position(const std::string &cid) const {

   std::vector<double> v;
   mmdb::Residue *residue_p = cid_to_residue(cid);
   if (residue_p) {
      std::vector<clipper::Coord_orth> atom_positions;
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            clipper::Coord_orth p = co(at);
            atom_positions.push_back(p);
         }
      }
      if (! atom_positions.empty()) {
         clipper::Coord_orth sum(0,0,0);
         for (const auto &pos : atom_positions)
            sum += pos;
         double is = 1.0/static_cast<double>(atom_positions.size());
         v = {sum.x() * is, sum.y() * is, sum.z() * is};
      }
   }
   return v;
}

//! get the avarge residue side-chain position
//!
//! @return a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values
std::vector<double>
coot::molecule_t::get_residue_sidechain_average_position(const std::string &cid) const {

   std::vector<double> v;
   mmdb::Residue *residue_p = cid_to_residue(cid);
   if (residue_p) {
      std::vector<clipper::Coord_orth> side_chain_positions;
      std::set<std::string> main_chain_atoms;
      main_chain_atoms.insert(" CA ");
      main_chain_atoms.insert(" C  ");
      main_chain_atoms.insert(" N  ");
      main_chain_atoms.insert(" O  ");
      main_chain_atoms.insert(" H  ");
      main_chain_atoms.insert(" HA ");
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string atom_name(at->GetAtomName());
            if (main_chain_atoms.find(atom_name) == main_chain_atoms.end()) {
               clipper::Coord_orth p = co(at);
               side_chain_positions.push_back(p);
            }
         }
      }
      if (! side_chain_positions.empty()) {
         clipper::Coord_orth sum(0,0,0);
         for (const auto &pos : side_chain_positions)
            sum += pos;
         double is = 1.0/static_cast<double>(side_chain_positions.size());
         v = {sum.x() * is, sum.y() * is, sum.z() * is};
      }
   }
   return v;
}

std::pair<int, double>
coot::molecule_t::get_torsion(const std::string &cid, const std::vector<std::string> &atom_names) const {
   std::pair<int, double> p(0,0);
   mmdb::Residue * residue_p = cid_to_residue(cid);
   if (residue_p) {
      if (atom_names.size() == 4) {
         atom_name_quad quad(atom_names[0], atom_names[1], atom_names[2], atom_names[3]);
         double torsion = quad.torsion(residue_p);
         p.first =1;
         p.second = torsion;
      }

   }

   return p;
}

//! set occupancy
//!
//! set the occupancy for the given atom selection
//!
//! @param imol is the model molecule index
//! @param cod is the atom selection CID
void
coot::molecule_t::set_occupancy(const std::string &cid, float occ_new) {

   int selHnd = atom_sel.mol->NewSelection(); // d
   mmdb::Atom **SelAtoms = nullptr;
   int nSelAtoms = 0;
   atom_sel.mol->Select(selHnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_NEW);
   atom_sel.mol->GetSelIndex(selHnd, SelAtoms, nSelAtoms);

   for (int i=0; i<nSelAtoms; i++) {
      mmdb:: Atom *at = SelAtoms[i];
      if (! at->isTer()) {
         at->occupancy = occ_new;
      }
   }
   atom_sel.mol->DeleteSelection(selHnd);
}


std::vector<std::pair<std::string, std::string> >
coot::molecule_t::get_sequence_info() const {

   std::vector<std::pair<std::string, std::string> > v;
   return v;

}

#include <mmdb2/mmdb_math_align.h>
#include <mmdb2/mmdb_tables.h>
#include "coot-utils/coot-align.hh"

//! return the mismatches/mutations:
coot::chain_mutation_info_container_t
coot::molecule_t::get_mutation_info() const {

   auto get_chain_sequence = [] (mmdb::Chain *chain_p,
                                 mmdb::Residue **selResidues,
                                 int nSelResidues) {
      std::string s;
      for (int ires=0; ires<nSelResidues; ires++) {
         mmdb::Residue *residue_p = selResidues[ires];
         char r[10];
         if (residue_p) {
            auto rn = residue_p->GetResName();
            mmdb::Get1LetterCode(rn, r);
            s += r;
         }
      }
      return s;
   };

   auto align_on_chain = [get_chain_sequence] (mmdb::Chain *chain_p,
                                               mmdb::PResidue *SelResidues,
                                               int nSelResidues,
                                               const std::string &target) {

      chain_mutation_info_container_t ch_info(chain_p->GetChainID());
      std::string model = get_chain_sequence(chain_p, SelResidues, nSelResidues);
      mmdb::math::Alignment align;
      mmdb::realtype wgap = 0.0;
      mmdb::realtype wspace = -1.0;
      std::string stripped_target = util::remove_whitespace(target);
      align.Align(model.c_str(), stripped_target.c_str());
      std::string s = align.GetAlignedS();
      std::string t = align.GetAlignedT();
      std::cout << "model:  " << s << std::endl;
      std::cout << "target: " << t << std::endl;
      std::cout << "score:  " << align.GetScore() << std::endl;
      if (s.length() == t.length()) {
	 std::vector<int> resno_offsets(s.length(), 0);
	 for (unsigned int ires=0; ires<s.length(); ires++) {
	    if (s[ires] != t[ires]) {
	       // Case 1: simple mutation:
	       if (s[ires] != '-' && t[ires] != '-') {
		  std::string target_type = coot::util::single_letter_to_3_letter_code(t[ires]);
		  residue_spec_t res_spec(ires);
		  ch_info.add_mutation(res_spec, target_type);
	       }

	       // Case 2: model has an insertion
	       if (s[ires] != '-' && t[ires] == '-') {
		  for (unsigned int i=ires+1; i<s.length(); i++)
		     resno_offsets[i] -= 1;
		  residue_spec_t res_spec(ires);
		  ch_info.add_deletion(res_spec);
	       }

	       // Case 3: model has a deletion
	       if (s[ires] == '-' && t[ires] != '-') {
		  for (unsigned int i=ires+1; i<s.length(); i++)
		     resno_offsets[i] += 1;
		  residue_spec_t res_spec(ires);
		  std::string target_type = coot::util::single_letter_to_3_letter_code(t[ires]);
		  ch_info.add_insertion(res_spec, target_type);
	       }
	    }
         }
      }
      ch_info.rationalize_insertions();
      return ch_info;
   };

   chain_mutation_info_container_t ch_info;

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 mmdb::Chain *chain_p = model_p->GetChain(ichain);
	 std::string chain_id(chain_p->GetChainID());
	 std::pair<bool, std::string> target_pair = multi_fasta_seq.get_fasta_for_name(chain_id);
	 if (target_pair.first) {
	    std::string target = target_pair.second;
	    if (! target.empty()) {
	       int n_res = chain_p->GetNumberOfResidues();
	       mmdb::Residue **chain_residues = 0;
	       chain_p->GetResidueTable(chain_residues, n_res);
	       ch_info = align_on_chain(chain_p, chain_residues, n_res, target);
	    }
	 }
      }
   }

   return ch_info;
}

void
coot::molecule_t::set_temperature_factors_using_cid(const std::string &cid, float temp_fact) {

   if (atom_sel.mol) {
      int selHnd = atom_sel.mol->NewSelection(); // d
      mmdb::Atom **SelAtoms = nullptr;
      int nSelAtoms = 0;
      atom_sel.mol->Select(selHnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selHnd, SelAtoms, nSelAtoms);
      if (nSelAtoms > 0) {
         for (int i=0; i<nSelAtoms; i++) {
            mmdb::Atom *atom = SelAtoms[i];
            atom->tempFactor = temp_fact;
         }
      }
      atom_sel.mol->DeleteSelection(selHnd);
   }
}

#include "coot-utils/json.hpp"
using json = nlohmann::json;

//! get pucker info
//!
//! @return a json string or an empty string on failure
std::string
coot::molecule_t::get_pucker_analysis_info() const {

   std::string s;

   std::vector<std::pair<mmdb::Residue *, pucker_analysis_info_t> > puckers;
   std::string alt_conf = "";
   if (atom_sel.mol) {
      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            if (n_res > 1) {
               for (int ires=0; ires<(n_res-1); ires++) {
                  mmdb::Residue *residue_p      = chain_p->GetResidue(ires);
                  mmdb::Residue *residue_next_p = chain_p->GetResidue(ires+1);
                  if (residue_p) {
                     if (residue_p->GetNumberOfAtoms() > 14) {
                        try {
                           pucker_analysis_info_t pai(residue_p, alt_conf);
                           double d = pai.phosphate_distance_to_base_plane(residue_next_p);
                           // store d in the markup info
                           pai.markup_info.phosphate_distance_to_base_plane = d;
                           puckers.push_back(std::make_pair(residue_p, pai));
                        }
                        catch (const std::runtime_error &e) {
                           // it's OK.
                           // std::cout << "WARNING::" << e.what() << std::endl;
                        }
                     }
                  }
               }
            }
         }
      }
   }
   if (! puckers.empty()) {

      json j = json::array();
      for (unsigned int i=0; i<puckers.size(); i++) {
         const auto &pi = puckers[i].second;
         mmdb::Residue *residue_p  = puckers[i].first;
         json j_plane_distortion = pi.plane_distortion;
         json j_out_of_plane_distance = pi.out_of_plane_distance;
         json j_markup_info_phosphorus_distance_to_base_plane = pi.markup_info.phosphate_distance_to_base_plane;
         json j_puckered_atom = pi.puckered_atom();
         json j_res_name = residue_p->GetResName();
         json j_chain_id = residue_p->GetChainID();
         json j_res_no = residue_p->GetSeqNum();
         json j_markup_info_base_ring_centre;
         json j_markup_info_base_ring_normal;
         json j_markup_info_base_phosphorus_position;
         json j_markup_info_base_projected_point;
         j_markup_info_base_ring_centre["x"] = pi.markup_info.base_ring_centre.x();
         j_markup_info_base_ring_centre["y"] = pi.markup_info.base_ring_centre.y();
         j_markup_info_base_ring_centre["z"] = pi.markup_info.base_ring_centre.z();
         j_markup_info_base_ring_normal["x"] = pi.markup_info.base_ring_normal.x();
         j_markup_info_base_ring_normal["y"] = pi.markup_info.base_ring_normal.y();
         j_markup_info_base_ring_normal["z"] = pi.markup_info.base_ring_normal.z();
         j_markup_info_base_phosphorus_position["x"] = pi.markup_info.phosphorus_position.x();
         j_markup_info_base_phosphorus_position["y"] = pi.markup_info.phosphorus_position.y();
         j_markup_info_base_phosphorus_position["z"] = pi.markup_info.phosphorus_position.z();
         j_markup_info_base_projected_point["x"] = pi.markup_info.projected_point.x();
         j_markup_info_base_projected_point["y"] = pi.markup_info.projected_point.y();
         j_markup_info_base_projected_point["z"] = pi.markup_info.projected_point.z();
         json j_pucker;
         j_pucker["plane_distortion"]      = j_plane_distortion;
         j_pucker["out_of_plane_distance"] = j_out_of_plane_distance;
         j_pucker["puckered_atom"]         = j_puckered_atom;
         j_pucker["chain_id"]              = j_chain_id;
         j_pucker["res_no"]                = j_res_no;
         j_pucker["res_name"]              = j_res_name;
         j_pucker["base_ring_centre"]             = j_markup_info_base_ring_centre;
         j_pucker["base_ring_normal"]             = j_markup_info_base_ring_normal;
         j_pucker["phosphorus_position"]          = j_markup_info_base_phosphorus_position;
         j_pucker["projected_point"]             = j_markup_info_base_projected_point;
         j_pucker["phosphate_distance_to_base_plane"] = j_markup_info_phosphorus_distance_to_base_plane;
         j.push_back(j_pucker);
      }
      s = j.dump(4);
   }
   return s;
}
