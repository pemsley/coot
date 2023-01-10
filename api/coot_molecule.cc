
#include <iostream>
#include <sstream>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coords/mmdb.hh"
#include "coot_molecule.hh"
#include "ideal/pepflip.hh"
#include "rama-plot-phi-psi.hh"

bool
coot::molecule_t::is_valid_model_molecule() const {

   bool status = false;
   if (atom_sel.mol)
      status = true;
   return status;

}

int
coot::molecule_t::close_yourself() {

   int status = 1;
   if (is_valid_model_molecule()) {
      atom_sel.clear_up();
      status = 1;
   }
   if (is_valid_map_molecule()) {
      clipper::Xmap<float> xmap_empty;
      xmap = xmap_empty;
      status = 1;
   }
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
      mmdb::Atom **SelAtoms;
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


// restore from (previous) backup
void
coot::molecule_t::restore_from_backup(int mod_index, const std::string &cwd) {

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
   atom_selection_container_t asc = get_atom_selection(file_name);
   if (asc.read_success) {
      save_info.set_modification_index(mod_index);
      atom_sel.clear_up();
      atom_sel = asc;
      // 20221018-PE no bond regeneration, maybe there should be?
   }

}

int
coot::molecule_t::undo() {

   make_backup();
   int status = 0;
   std::string cwd = coot::util::current_working_dir();
   int prev_mod_index = save_info.get_previous_modification_index();
   std::cout << ":::::::::::::::: undo requests prev_mod_index " << prev_mod_index << std::endl;
   restore_from_backup(prev_mod_index, cwd);
   return status;
}

int
coot::molecule_t::redo() {

   int status = 0;
   std::string cwd = coot::util::current_working_dir();
   int mod_index = save_info.get_next_modification_index();
   std::cout << ":::::::::::::::: redo requests mod_index " << mod_index << std::endl;
   restore_from_backup(mod_index, cwd);
   return status;
}

int
coot::molecule_t::write_coordinates(const std::string &file_name) const {

   int err = 1;
   if (atom_sel.n_selected_atoms > 0) {
      std::string ext = coot::util::file_name_extension(file_name);
      if (coot::util::extension_is_for_shelx_coords(ext)) {
         write_shelx_ins_file(file_name);
      } else {
         mmdb::byte bz = mmdb::io::GZM_NONE; // 20221018-PE  this should be used too
         err = coot::write_coords_pdb(atom_sel.mol, file_name);
      }
   }
   return err;
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
   std::string index_string = save_info.index_string();
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

void
coot::molecule_t::save_history_file_name(const std::string &file) {

   // 20221016-PE fold this function into save_info

   // First, history_index is zero and the vec is zero,
   // normal service, then another backup: history_index is 1 and vec is 1.
   //
   int history_filename_vec_size = history_filename_vec.size();
   if (save_info.modification_index == history_filename_vec_size) {
      history_filename_vec.push_back(file);
   } else {
      // we have gone back in history.
      //
      if (save_info.modification_index < history_filename_vec_size) {
         history_filename_vec[save_info.modification_index] = file;
      }
   }
}


#include <sys/stat.h>

std::string
coot::molecule_t::make_backup() {

   if (false) {
      std::cout << "start make_backup() for molecule " << imol_no << std::endl;
      std::cout << "start make_backup() for molecule with " << atom_sel.n_selected_atoms << " atoms " << std::endl;
   }

   // nothing done yet.

   std::string info_message;

   auto make_maybe_backup_dir = [] (const std::string &backup_dir) {
      return util::create_directory(backup_dir);
   };

   bool backup_this_molecule = true; // if needed, give user control later
   bool backup_compress_files_flag = false;

   if (backup_this_molecule) {
      std::string backup_dir("coot-backup");

      //shall we use the environment variable instead?
      char *env_var = getenv("COOT_BACKUP_DIR");

#ifdef EMSCRIPTEN

#else

      if (env_var) {
         struct stat buf;

         // we better debackslash the directory (for windows)
         std::string tmp_dir = env_var;
         tmp_dir = coot::util::intelligent_debackslash(tmp_dir);
         int err = stat(tmp_dir.c_str(), &buf);

         if (!err) {
            if (! S_ISDIR(buf.st_mode)) {
               env_var = NULL;
            }
         } else {
            env_var = NULL;
         }
      }
#endif

      if (env_var)
         backup_dir = env_var;

      if (atom_sel.mol) {
         int dirstat = make_maybe_backup_dir(backup_dir);

         if (dirstat != 0) {
            // fallback to making a directory in $HOME
            std::string home_dir = coot::get_home_dir();
            if (! home_dir.empty()) {
               backup_dir = coot::util::append_dir_dir(home_dir, "coot-backup");
               dirstat = make_maybe_backup_dir(backup_dir);
               if (dirstat != 0) {
                  std::cout << "WARNING:: backup directory "<< backup_dir
                            << " failure to exist or create" << std::endl;
               } else {
                  std::cout << "INFO using backup directory " << backup_dir << std::endl;
               }
            } else {
               std::cout << "WARNING:: backup directory "<< backup_dir
                         << " failure to exist or create" << std::endl;
            }
         }

         if (dirstat == 0) {
            // all is hunkey-dorey.  Directory exists.

            std::string backup_file_name = get_save_molecule_filename(backup_dir);
            std::cout << "INFO:: make_backup() backup file name " << backup_file_name << std::endl;

            mmdb::byte gz;
            if (backup_compress_files_flag) {
               gz = mmdb::io::GZM_ENFORCE;
            } else {
               gz = mmdb::io::GZM_NONE;
            }

            // Writing out a modified binary mmdb like this results in the
            // file being unreadable (crash in mmdb read).
            //
            int istat;
            if (! is_from_shelx_ins_flag) {
               bool write_as_cif = false;
               if (coot::is_mmcif_filename(name))
                  write_as_cif = true;

               istat = write_atom_selection_file(atom_sel, backup_file_name, write_as_cif, gz);

               if (true) { // 20221021-PE this should not be needed
                  struct stat buf;
                  int err = stat(backup_file_name.c_str(), &buf);
                  if (err == 0) {
                     std::cout << "DEBUG:: in make_backup() " << backup_file_name << " confirmed as existing" << std::endl;
                  } else {
                     std::cout << "DEBUG:: in make_backup() " << backup_file_name << " does not exist!" << std::endl;
                  }
               }

               // WriteMMDBF returns 0 on success, else mmdb:Error_CantOpenFile (15)
               if (istat) {
                  std::string warn;
                  warn = "WARNING:: WritePDBASCII failed! Return status ";
                  warn += istat;
                  // g.info_dialog_and_text(warn);
                  info_message = warn;
               }
            } else {
               std::pair<int, std::string> p = write_shelx_ins_file(backup_file_name);
               istat = p.first;
            }

            save_history_file_name(backup_file_name);
            // 20221016-PE old history counting - now use save_info.
            // if (history_index == max_history_index)
            // max_history_index++;
            // history_index++;
         }
      } else {
         std::cout << "WARNING:: BACKUP:: Ooops - no atoms to backup for this empty molecule"
                   << std::endl;
      }
   } else {
      // Occasionally useful but mostly tedious...
      // std::cout << "INFO:: backups turned off on this molecule"
      // << std::endl;
   }
   return info_message;
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
coot::molecule_t::select_residues(const residue_spec_t &residue_spec, const std::string &mode) const {

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

   if (false)
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

   // make_backup();// 20221016-PE why was this here? Ah, for interactive use/intermdiate atoms
                    // let's remove it for now

   if (false) {
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

            if (false)
               std::cout << "DEBUG:: asc.UDDOldAtomIndexHandle is "
                         << asc.UDDOldAtomIndexHandle << " using full atom spec to atom index..."
                         << std::endl;

            idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
                                               atom->residue->seqNum,
                                               std::string(atom->GetInsCode()),
                                               std::string(atom->name),
                                               std::string(atom->altLoc));

            // std::cout << "full_atom_spec_to_atom_index() returned " << idx << " for " << coot::atom_spec_t(atom) << std::endl;

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
                  // say.  Now, the alt conf atoms has been immmediately
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
               if (idx <= atom_sel.n_selected_atoms) {
                  mmdb::Atom *mol_atom = atom_sel.atom_selection[idx];
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
                  std::cout << "Trapped error! idx " << idx << " but atom_sel.n_selected_atoms " << atom_sel.n_selected_atoms
                            << std::endl;
               }
            } else {
               std::cout << "WARNING:: bad atom idx -1" << std::endl;
            }
         }
      }
   }
   std::cout << "INFO:: replace_coords: " << n_atom << " atoms updated." << std::endl;

   // have_unsaved_changes_flag = 1;
   // save_info.new_modification("replace_coords()");

   if (show_symmetry) {  // internal
      update_symmetry();
   }

   // make_bonds_type_checked(__FUNCTION__);

}


int coot::molecule_t::flip_peptide(const coot::atom_spec_t &as_in, const std::string &alt_conf) {

   make_backup();
   coot::atom_spec_t as = as_in;
   if (as.atom_name == " N  ")
      as.res_no--;
   int result = coot::pepflip(atom_sel.mol, as.chain_id, as.res_no, as.ins_code, alt_conf);
   save_info.new_modification("flip_peptide");
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
   return glm::vec4(ch.red, ch.green, ch.blue, 1.0f);
}

#include "utils/dodec.hh"

// returns either the specified residue or null if not found
mmdb::Residue *
coot::molecule_t::get_residue(const coot::residue_spec_t &residue_spec) const {

   mmdb::Residue *r = coot::util::get_residue(residue_spec, atom_sel.mol);
   return r;

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

   if (true)
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

         glm::vec3 atom_pos = cartesian_to_glm(rm.pos) + cartesian_to_glm(offset);
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

         // change the colour index of dodec_triangles
         for (auto &tri : dodec_triangles)
            tri.colour_index = i;
         unsigned int idx_base = vertices.size();
         unsigned int idx_tri_base = triangles.size();
         vertices.insert(vertices.end(), this_dodec_vertices.begin(), this_dodec_vertices.end());
         triangles.insert(triangles.end(), dodec_triangles.begin(), dodec_triangles.end());
         for (unsigned int jj=idx_tri_base; jj<triangles.size(); jj++)
            triangles[jj].rebase(idx_base);
      }
   }
   std::cout << "DEBUG:: ending get_rotamer_dodecs() with mesh " << m.vertices.size() << " vertices and "
             << m.triangles.size() << " triangle." << std::endl;
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

   if (true)
      std::cout << "DEBUG:: in get_rotamer_dodecs_instanced() bonds_box.n_rotamer_markups " << bonds_box.n_rotamer_markups
                << std::endl;

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

   std::cout << "DEBUG:: in coot::molecule_t::get_rotamer_dodecs_instanced(): n_rotamer_markups: " << bonds_box.n_rotamer_markups << std::endl;

   for (int i=0; i<bonds_box.n_rotamer_markups; i++) {
      const rotamer_markup_container_t &rm = bonds_box.rotamer_markups[i];
      const residue_spec_t &residue_spec = rm.spec;
      mmdb::Residue *residue_p = get_residue(residue_spec);
      Cartesian offset(0,0,rama_ball_pos_offset_scale);
      if (residue_p) {
         std::pair<bool, coot::Cartesian> hav = get_HA_unit_vector(residue_p);
         if (hav.first) offset = hav.second * 1.6;
      }
      glm::vec3 atom_pos = cartesian_to_glm(rm.pos) + cartesian_to_glm(offset);
      auto rm_col = rm.col;
      rm_col.scale_intensity(0.75); // was 0.6 in Mesh-from-graphical-bonds.cc
      auto this_dodec_colour = colour_holder_to_glm(rm_col);
      instancing_data_type_A_t id(atom_pos, this_dodec_colour, size_3);
      // std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!! adding to data A!" << std::endl;
      ig.instancing_data_A.push_back(id);
   }

   m.add(ig);
   return m;
}


#include "ligand/backrub-rotamer.hh"

std::pair<bool,float>
coot::molecule_t::backrub_rotamer(const std::string &chain_id, int res_no,
                                  const std::string &ins_code, const std::string &alt_conf,
                                  const clipper::Xmap<float> &xmap_in,
                                  const coot::protein_geometry &pg) {

   bool status = false;
   float score = -1;
   bool refinement_move_atoms_with_zero_occupancy_flag = true; // pass this?

   std::cout << "debug:: molecule_t::backrub_rotamer() starts " << chain_id << " " << res_no << std::endl;

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
         try {

            make_backup();
            mmdb::Residue *prev_res = coot::util::previous_residue(res);
            mmdb::Residue *next_res = coot::util::next_residue(res);
            mmdb::Manager *mol = atom_sel.mol;
            coot::backrub br(chain_id, res, prev_res, next_res, alt_conf, mol,
                             &xmap_in); // use a pointer for the map
            std::pair<coot::minimol::molecule,float> m = br.search(restraints);
            std::vector<coot::atom_spec_t> baddie_waters = br.waters_for_deletion();
            score = m.second;
            status = true;
            atom_selection_container_t fragment_asc = make_asc(m.first.pcmmdbmanager());
            replace_coords(fragment_asc, 0, refinement_move_atoms_with_zero_occupancy_flag);
            if (baddie_waters.size())
               delete_atoms(baddie_waters);

            atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
            atom_sel.mol->FinishStructEdit();
         }
         catch (const std::runtime_error &rte) {
            std::cout << "WARNING:: thrown " << rte.what() << std::endl;
         }
      } else {
         std::cout << " No restraints found for " << monomer_type << std::endl;
      }
   }
   if (status) {
      write_coordinates("post_backrub_rotamer.pdb");
      save_info.new_modification("backrub_rotamer()");
   }
   return std::pair<bool,float> (status, score);
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

      make_backup();
      bool was_deleted = false;
      // do we include CB? I forget.
      std::vector<std::string> main_chain_atoms_list = { " C  ", " N  ", " H  ", " O  ", " CA ",  " HA ", " CB " };
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         std::string atom_name(at->GetAtomName());
         if (std::find(main_chain_atoms_list.begin(), main_chain_atoms_list.end(), atom_name) == main_chain_atoms_list.end()) {
            delete at;
            was_deleted = true;
         }
      }

      if (was_deleted) {
         status = true;
         atom_sel.mol->FinishStructEdit();
         atom_sel = make_asc(atom_sel.mol);
         // make_bonds_type_checked(__FUNCTION__);
         // have_unsaved_changes_flag = 1;
         save_info.new_modification("delete_side_chain()");
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
      if (atom_specs.size() > 0)
         make_backup();
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
      save_info.new_modification("delete_atoms");
      // unlikely to be necessary:
      trim_atom_label_table();

      // update_symmetry();
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

                              make_backup();
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
      save_info.new_modification("delete_atom");
      // unlikely to be necessary:
      trim_atom_label_table();
      update_symmetry();
   }
   return was_deleted;
}


void
coot::molecule_t::delete_any_link_containing_residue(const coot::residue_spec_t &res_spec) {

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
                              make_backup();
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

      save_info.new_modification(__FUNCTION__);
      trim_atom_label_table();
      update_symmetry();
   }
   return was_deleted;
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
                     make_backup();
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
coot::molecule_t::delete_literal_using_cid(const std::string &atom_selection_cid) {

   // cid is an atom selection, e.g. containing a residue range

   int status = 0;
   std::vector<mmdb::Atom *> atoms_to_be_deleted;

   mmdb::Atom **selection_atoms = 0;
   int n_selection_atoms = 0;
   int selHnd = atom_sel.mol->NewSelection(); // d
   atom_sel.mol->Select(selHnd, mmdb::STYPE_ATOM, atom_selection_cid.c_str(), mmdb::SKEY_NEW);
   atom_sel.mol->GetSelIndex(selHnd, selection_atoms, n_selection_atoms);

   if (selection_atoms) {
      for (int iat=0; iat<n_selection_atoms; iat++) {
         mmdb::Atom *at = selection_atoms[iat];
         if (at)
            atoms_to_be_deleted.push_back(at);
      }
   }

   if (! atoms_to_be_deleted.empty()) {
      make_backup();

      for (unsigned int iat=0; iat<atoms_to_be_deleted.size(); iat++) {
         delete atoms_to_be_deleted[iat];
         atoms_to_be_deleted[iat] = NULL;
      }
      status = 1;
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      coot::util::pdbcleanup_serial_residue_numbers(atom_sel.mol);

      std::string s = std::string("delete-literal-using-cid ") + atom_selection_cid;
      save_info.new_modification(s);
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

         std::cout << "DEBUG:: Sanity check A in molecule_t::sfcalc_genmaps_using_bulk_solvent(): HKL_info: "
                   << "cell: " << hkls_check.cell().format() << " "
                   << "spacegroup: " << spgr_check.symbol_xhm() << " "
                   << "resolution: " << hkls_check.resolution().limit() << " "
                   << "invsqreslim: " << hkls_check.resolution().invresolsq_limit() << " "
                   << std::endl;
         std::cout << "DEBUG:: Sanity check B in molecule_t::sfcalc_genmaps_using_bulk_solvent(): Cell fofc-map"
                   << xmap_fofc_p->cell().format() << std::endl;
      }

      stats = coot::util::sfcalc_genmaps_using_bulk_solvent(atom_sel.mol, fobs, free, cell, xmap_2fofc_p, xmap_fofc_p);

      // maybe format() should be inside coot::util::sfcalc_genmap_stats_t
      std::cout << "\n R-factor      : " << stats.r_factor << "\n Free R-factor : " << stats.free_r_factor << "\n";
      std::cout << "\n Bulk Correction Volume: " << stats.bulk_solvent_volume;
      std::cout << "\n Bulk Correction Factor: " << stats.bulk_correction << "\n";
      std::cout << "\nNumber of spline params: " << stats.n_splines << "\n";

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
                                float map_weight, const coot::protein_geometry &geom, bool refinement_is_quiet) {

   bool make_trans_peptide_restraints = true;
   bool do_rama_plot_restraints = false;

   int status =  0;
   std::vector<coot::atom_spec_t> fixed_atom_specs;
   std::vector<std::pair<bool,mmdb::Residue *> > local_residues;
   for (const auto &r : rv)
      local_residues.push_back(std::make_pair(false, r));

   if (true) {

      make_backup();
      mmdb::Manager *mol = atom_sel.mol;
      std::vector<mmdb::Link> links;
      coot::restraints_container_t restraints(local_residues,
                                              links,
                                              geom,
                                              mol,
                                              fixed_atom_specs, &xmap);

      if (refinement_is_quiet)
         restraints.set_quiet_reporting();

      std::cout << "DEBUG:: using restraints with map_weight " << map_weight << std::endl;
      restraints.add_map(map_weight);
      coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
      flags = coot::TYPICAL_RESTRAINTS;
      coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;

      int n_threads = 4; // coot::get_max_number_of_threads();
      ctpl::thread_pool thread_pool(n_threads);
      restraints.thread_pool(&thread_pool, n_threads);

      int imol = 0; // dummy
      restraints.make_restraints(imol, geom, flags, 1, make_trans_peptide_restraints,
                                 1.0, do_rama_plot_restraints, true, true, false, pseudos);
      int nsteps_max = 4000;
      short int print_chi_sq_flag = 1;
      restraints.minimize(flags, nsteps_max, print_chi_sq_flag);
      coot::geometry_distortion_info_container_t gd = restraints.geometric_distortions();
      if (! refinement_is_quiet)
         gd.print();

      save_info.new_modification("refine_direct");

   }
   return status;
}


#include "add-terminal-residue.hh"

std::pair<int, std::string>
coot::molecule_t::add_terminal_residue_directly(const residue_spec_t &spec, const std::string &new_res_type,
                                                const coot::protein_geometry &geom,
                                                const clipper::Xmap<float> &xmap) {

   std::pair<int, std::string> r;
   mmdb::Residue *residue_p = util::get_residue(spec, atom_sel.mol);
   if (residue_p) {
      std::string terminus_type = coot::get_term_type(residue_p, atom_sel.mol);
      float bf_new = default_temperature_factor_for_new_atoms;
      make_backup();
      r = add_terminal_residue(imol_no, terminus_type, residue_p,
                               atom_sel.mol, atom_sel.UDDAtomIndexHandle,
                               spec.chain_id, new_res_type,
                               bf_new, xmap, geom);
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      coot::util::pdbcleanup_serial_residue_numbers(atom_sel.mol);
      atom_sel = make_asc(atom_sel.mol);
      save_info.new_modification("add-terminal-residue");
   } else {
      std::cout << "WARNING:: in add_terminal_residue_directly() null residue_p " << std::endl;
   }
   return r;
}


int
coot::molecule_t::mutate(const coot::residue_spec_t &spec, const std::string &new_res_type) {

   atom_sel.delete_atom_selection();
   mmdb::Residue *residue_p = coot::util::get_residue(spec, atom_sel.mol);
   int status = coot::util::mutate(residue_p, new_res_type);
   // std::cout << "mutate status " << status << std::endl;

   // 20221121-PE should thise function calls be in coot::util::mutate()? I think so.
   atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol); // regen the atom indices
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

   make_backup();

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
      save_info.new_modification(__FUNCTION__);
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

         make_backup();

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

                  coot::contact_info contact = coot::getcontacts(residue_asc, monomer_type, imol_no, geom);
                  std::vector<std::vector<int> > contact_indices =
                     contact.get_contact_indices_with_reverse_contacts();

                  try {
                     coot::atom_tree_t tree(contact_indices, clicked_atom_idx, residue, alt_conf);
                     problem_string = jed_flip_internal(tree, interesting_torsions, atom_name, invert_selection);
                     atom_sel.mol->FinishStructEdit();
                  }
                  catch (const std::runtime_error &rte) {
                     std::cout << "RUNTIME ERROR:: " << rte.what() << " - giving up" << std::endl;
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

      make_backup();
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
      save_info.new_modification(__FUNCTION__);
      std::cout << "DEBUG:: eigen_flip_residue() now save_info modification_index " << save_info.modification_index
                << std::endl;

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
   }
   return n_atoms_moved;
}


int
coot::molecule_t::new_positions_for_residue_atoms(const std::string &residue_cid, const std::vector<moved_atom_t> &moved_atoms) {

   mmdb::Residue *residue_p = cid_to_residue(residue_cid);
   return new_positions_for_residue_atoms(residue_p, moved_atoms);
}

int
coot::molecule_t::new_positions_for_residue_atoms(mmdb::Residue *residue_p, const std::vector<moved_atom_t> &moved_atoms) {

   int n_atoms_moved = 0;
   if (residue_p) {
      for (unsigned int i=0; i<moved_atoms.size(); i++) {
         const moved_atom_t &mva = moved_atoms[i];

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
   }
   return n_atoms_moved;
}

int
coot::molecule_t::new_positions_for_atoms_in_residues(const std::vector<moved_residue_t> &moved_residues) {

   int status = 0;
   for (unsigned int i=0; i<moved_residues.size(); i++) {
      const moved_residue_t &mvr = moved_residues[i];
      coot::residue_spec_t res_spec(mvr.chain_id, mvr.res_no, mvr.ins_code);
      mmdb::Residue *residue_p = get_residue(res_spec);
      new_positions_for_residue_atoms(residue_p, mvr.moved_atoms);
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

  make_backup();
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
  save_info.new_modification(__FUNCTION__);
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
coot::molecule_t::insert_waters_into_molecule(const coot::minimol::molecule &water_mol) {

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
      make_backup();
      std::cout << "INFO:: Adding to solvent chain: " << chain_p->GetChainID()
                << std::endl;
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
               new_residue_p->SetResName("HOH");
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
      atom_sel.mol->FinishStructEdit();
      // update_molecule_after_additions(); // sets unsaved changes flag
      update_symmetry();
   }
   return istat;
}


int
coot::molecule_t::append_to_molecule(const coot::minimol::molecule &water_mol) {

   int istat = 0; // fail status initially.
   int n_atom = 0;  // 0 new atoms added initially.
   float default_new_atoms_b_factor = 20.0;

   if (atom_sel.n_selected_atoms > 0) {

      make_backup();

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
         save_info.new_modification(__FUNCTION__);
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
