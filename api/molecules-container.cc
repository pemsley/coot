/*
 * api/molecules-container.cc
 *
 * Copyright 2020, 2021, 2-22 by Medical Research Council
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

#include <iomanip>
#include <filesystem>

#include <sys/types.h> // for stating
#include <sys/stat.h>

#include "molecules-container.hh"
#include "geometry/protein-geometry.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "ideal/pepflip.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/secondary-structure-headers.hh"
#include "coot-utils/oct.hh"
#include "coot-utils/read-sm-cif.hh"

#include "coords/Bond_lines.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-extras.hh"

#include "mmdb2/mmdb_atom.h"
#include "utils/logging.hh"
extern logging logger;

// statics
std::atomic<bool> molecules_container_t::restraints_lock(false);
std::atomic<bool> molecules_container_t::on_going_updating_map_lock(false);
std::string molecules_container_t::restraints_locking_function_name; // I don't know why this needs to be static
std::vector<atom_pull_info_t> molecules_container_t::atom_pulls;
// 20221018-PE not sure that this needs to be static.
clipper::Xmap<float> *molecules_container_t::dummy_xmap = new clipper::Xmap<float>;
//bool molecules_container_t::make_backups_flag = true;

molecules_container_t::~molecules_container_t() {

   standard_residues_asc.clear_up();
}

//! the one and only constructor
// verbose is an optional arg, default true
molecules_container_t::molecules_container_t(bool verbose) :
   ramachandrans_container(ramachandrans_container_t()),
   thread_pool(8) {
   if (! verbose) geom.set_verbose(false);
   init();
}


//! init (private)
void
molecules_container_t::init() {

   use_gemmi = true;
   imol_refinement_map = -1;
   imol_difference_map = -1;
   setup_syminfo();
   mmdb::InitMatType();
   geometry_init_standard(); // do this by default now
   refinement_immediate_replacement_flag = true; // 20221018-PE for WebAssembly for the moment
   imol_moving_atoms = -1;
   refinement_is_quiet = true;
   show_timings = true;
   cif_dictionary_read_number = 40;
   // refinement
   continue_threaded_refinement_loop = false;
   particles_have_been_shown_already_for_this_round_flag = false;
   map_weight = 50.0;
   geman_mcclure_alpha = 0.01;

   map_sampling_rate = 1.8;
   draw_missing_residue_loops_flag = true;

   read_standard_residues();
   interrupt_long_term_job = false;
   contouring_time = 0;
   make_backups_flag = true;

   // thread_pool.resize(8); // now in the constructor

   use_rama_plot_restraints = false;
   rama_plot_restraints_weight = 1.0;

   use_torsion_restraints = false;
   torsion_restraints_weight = 1.0;

   map_is_contoured_using_thread_pool_flag = false;

   ligand_water_to_protein_distance_lim_max = 3.4;
   ligand_water_to_protein_distance_lim_min = 2.4;
   ligand_water_variance_limit = 0.1;
   ligand_water_sigma_cut_off = 1.75; // max moorhen points for tutorial 1.

   max_number_of_simple_mesh_vertices = 200000;

   // debug();
   //
   // size_t sss = sizeof(molecules_container_t);
   // std::cout << "::::::::::::::::: sizeof molecules_container_t " << sss << std::endl;
}

unsigned int
molecules_container_t::get_max_number_of_simple_mesh_vertices() const {
   return max_number_of_simple_mesh_vertices;
}

void
molecules_container_t::set_max_number_of_simple_mesh_vertices(unsigned int n) {
   max_number_of_simple_mesh_vertices = n;
}

//! don't use this in emscript
coot::molecule_t &
molecules_container_t::operator[] (unsigned int imol) {
   // maybe this should throw an exception on out-of-range?
   return molecules[imol];
}

//! don't use this in ecmascript
mmdb::Manager *
molecules_container_t::get_mol(unsigned int imol) const { // 20221018-PE function name change


   if (is_valid_model_molecule(imol)) {
      return molecules[imol].atom_sel.mol;
   } else {
      return nullptr;
   }
}

//! Fill the rotamer probability tables (currently not ARG and LYS)
void
molecules_container_t::fill_rotamer_probability_tables() {
   if (! rot_prob_tables.tried_and_failed()) {

      std::string tables_dir = coot::package_data_dir();
      char *data_dir = getenv("COOT_DATA_DIR");
      if (data_dir) {
         tables_dir = data_dir;
      }
      tables_dir += "/rama-data";
      rot_prob_tables.set_tables_dir(tables_dir);
      bool ignore_lys_and_arg_flag = true; // 20221018-PE remove this flag when rotamer probabiity
      // tables are read from a binary file (and is fast enough
      // to include lys and arg).
      rot_prob_tables.fill_tables(ignore_lys_and_arg_flag);
   }
}

bool
molecules_container_t::contains_unsaved_models() const {
   for (const auto &m : molecules) {
      if (m.have_unsaved_changes()) return true;
   }
   return false;
}

//! Save the unsaved model - this function has not yet been written!
void
molecules_container_t::save_unsaved_model_changes() {
   for (const auto &m : molecules) {
      if (m.have_unsaved_changes()) {
         // something fun here. - whatever it is though, don't put it in this header.
      }
   }
}



//! Get the package version
//!
//! @return the package version, e.g. "1.1.11" - if this is a not yet a release version
//! the version will end in a "+", such as "1.1.11+"
std::string
molecules_container_t::package_version() const {

   std::string s(PACKAGE_VERSION);
   return s;

}

//! get imol_enc_any
//!
//! @return the value of imol_enc_any (meaning "the molecule number for dictionary that
// can be used with any molecule")
int
molecules_container_t::get_imol_enc_any() const {

   return coot::protein_geometry::IMOL_ENC_ANY;
}


bool
molecules_container_t::is_valid_model_molecule(int imol) const {
   bool status = false;
   if (imol >= 0) {
      int ms = molecules.size();
      if (imol < ms) {
         status = molecules[imol].is_valid_model_molecule();
      }
   }
   return status;
}

bool
molecules_container_t::is_valid_map_molecule(int imol) const {

   bool status = false;
   if (imol >= 0) {
      int ms = molecules.size();
      if (imol < ms) {
         status = molecules[imol].is_valid_map_molecule();
      }
   }
   return status;
}

//! make the logging output go to a file
//!
//! @param file_name the looging file name
void
molecules_container_t::set_logging_file(const std::string &file_name) {

   logger.set_log_file(file_name);

}


//! Control the logging
//!
//! @param level is the logging level
void
molecules_container_t::set_logging_level(const std::string &level) {

   logging::output_t ot = logging::output_t::TERMINAL;
   if (level == "LOW")       ot = logging::output_t::INTERNAL;
   if (level == "HIGH")      ot = logging::output_t::TERMINAL;
   if (level == "DEBUGGING") ot = logging::output_t::TERMINAL_WITH_DEBUGGING;
   logger.set_output_type(ot);
}



//! @return is this a difference map?
bool
molecules_container_t::is_a_difference_map(int imol) const {

   bool status = 0;
   if (is_valid_map_molecule(imol)) {
      status = molecules[imol].is_difference_map_p();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}


//! create an empty molecule
//! @return the index of the new molecule
int
molecules_container_t::new_molecule(const std::string &name) {

   int n_mol = molecules.size();
   molecules.push_back(coot::molecule_t(name, n_mol));
   return n_mol;
}


int
molecules_container_t::close_molecule(int imol) {

   int status = 0;
   int ms = molecules.size(); // type change
   if (imol < ms) {
      if (imol >= 0) {
         molecules[imol].close_yourself();
         status = 1;
      }
   }
   return status;
}

//! The adds a number of empty molecules to the internal vector of molecules
//! Note that his is not like `reserve` as it will increase the molecule index
//! of the next added molecule by `n_empty`.
void
molecules_container_t::create_empty_molecules(unsigned int n_empty) {

   if (n_empty > 0) {
      unsigned int n_start = molecules.size();
      for (unsigned int i=0; i<n_empty; i++) {
         molecules.push_back(coot::molecule_t("--empty--", n_start+i));
      }
   }
}


void
molecules_container_t::clear() {

   molecules.clear();
   molecules.shrink_to_fit();
}

//! delete the most recent/last molecule in the molecule vector
void
molecules_container_t::pop_back() {

   if (! molecules.empty()) {
      molecules.pop_back();
   }
}

// delete the most recent/last closed molecule in the molecule vector, until the first
// non-closed molecule is found (working from the end)
void
molecules_container_t::end_delete_closed_molecules() {

   if (! molecules.empty()) {
      while (true) {
         if (molecules.back().is_closed()) {
            std::vector<coot::molecule_t>::iterator end_iter = molecules.end();
            std::vector<coot::molecule_t>::iterator prev_iter = end_iter - 1;
            molecules.erase(prev_iter);
         } else {
            break;
         }
         if (molecules.empty()) break;
      }
   }
}


void
molecules_container_t::debug() const {

   // debug:
   char *env_var = getenv("SYMINFO");
   if (! env_var) {
      std::cout << "ERROR:: SYMINFO was not set" << std::endl;
   } else {
      std::string s(env_var);
      std::cout << "DEBUG:: SYMINFO was set to " << s << std::endl;

      struct stat buf;
      int status = stat(s.c_str(), &buf);
      if (status != 0) { // standard-residues file was not found in
                         // default location either...
        std::cout << "ERROR:: syminfo file " << s << " was not found" << std::endl;
      } else {
        std::cout << "DEBUG:: syminfo file " << s << " was found" << std::endl;
      }
   }
}

void
molecules_container_t::set_map_is_contoured_with_thread_pool(bool state) {
   map_is_contoured_using_thread_pool_flag = state;
}


std::string
molecules_container_t::get_molecule_name(int imol) const {

   int ms = molecules.size();
   if (imol < ms)
      if (imol >= 0)
         return molecules[imol].get_name();
   return std::string("");
}

//! set the molecule name
void
molecules_container_t::set_molecule_name(int imol, const std::string &new_name) {

   if (is_valid_map_molecule(imol))
      molecules[imol].set_molecule_name(new_name);
   if (is_valid_model_molecule(imol))
      molecules[imol].set_molecule_name(new_name);
}

api::cell_t
molecules_container_t::get_cell(int imol) const {

   api::cell_t c;
   if (is_valid_map_molecule(imol)) {
      clipper::Cell cell = molecules[imol].xmap.cell();
      c = api::cell_t(cell.a(), cell.b(), cell.c(), cell.alpha(), cell.beta(), cell.gamma());
   }

   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = molecules[imol].atom_sel.mol;
      mmdb::realtype a[6];
      mmdb::realtype vol;
      int orthcode;
      mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
      c = api::cell_t(a[0], a[1], a[2],
                      clipper::Util::d2rad(a[3]),
                      clipper::Util::d2rad(a[4]),
                      clipper::Util::d2rad(a[5]));
   }
   return c;
}

//! Get the middle of the "molecule blob" in cryo-EM reconstruction maps
//! @return a `coot::util::map_molecule_centre_info_t`.
coot::util::map_molecule_centre_info_t
molecules_container_t::get_map_molecule_centre(int imol) const {

   coot::util::map_molecule_centre_info_t mc;
   if (is_valid_map_molecule(imol)) {
      mc = molecules[imol].get_map_molecule_centre();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid map molecule " << imol << std::endl;
   }
   return mc;
}


void
molecules_container_t::display_molecule_names_table() const {

   for (unsigned int imol=0; imol<molecules.size(); imol++) {
      if (molecules[imol].is_closed()) {
         std::cout << imol << " ---closed---" << std::endl;
      } else {
         std::cout << imol << " " << std::setw(40) << molecules[imol].get_name() << std::endl;
      }
   }
}

unsigned int
molecules_container_t::get_number_of_atoms(int imol) const {

   unsigned int n = 0;
   if (is_valid_model_molecule(imol)) {
      n = molecules[imol].get_number_of_atoms();
   }
   return n;
}

int
molecules_container_t::get_number_of_hydrogen_atoms(int imol) const {

   int n = -1;
   if (is_valid_model_molecule(imol)) {
      n = molecules[imol].get_number_of_hydrogen_atoms();
   }
   return n;
}

void
molecules_container_t::set_draw_missing_residue_loops(bool state) {
   // std::cout << "****** in set_draw_missing_residue_loops() with state " << state << std::endl;
   draw_missing_residue_loops_flag = state;
}


void
molecules_container_t::testing_start_long_term_job(unsigned int n_seconds) {

   if (interrupt_long_term_job) {
      interrupt_long_term_job = false;
      return;
   }
   unsigned int n_ms_count = 0;
   unsigned int n_ms_per_cycle = 300;
   while (true) {
      double d = long_term_job_stats.time_difference();
      long_term_job_stats.function_value = 0.01 * d * d;
      if (interrupt_long_term_job) {
         interrupt_long_term_job = false;
         break;
      }
      if (n_seconds > 0)
         if (n_ms_count > n_seconds * 1000)
            break;
      std::this_thread::sleep_for(std::chrono::milliseconds(n_ms_per_cycle));
      n_ms_count += n_ms_per_cycle;
   }

}

void
molecules_container_t::testing_stop_long_term_job() {

   interrupt_long_term_job = true;

}

void
molecules_container_t::read_standard_residues() {

   std::string standard_env_dir = "COOT_STANDARD_RESIDUES";

   const char *env_var_filename = getenv(standard_env_dir.c_str());
   if (! env_var_filename) {

      std::string dir = coot::package_data_dir();
      std::string standard_file_name = coot::util::append_dir_file(dir, "standard-residues.pdb");

      std::filesystem::path standard_residues_file_path(standard_file_name);

      // failure case first:
      if (! std::filesystem::is_regular_file(standard_residues_file_path)) { // standard-residues file was not found in
                                                                             // default location either...
         std::cout << "WARNING:: default location: " << standard_file_name << std::endl;
         std::cout << "WARNING:: Can't find standard residues file in the default location \n";
         std::cout << "         and environment variable for standard residues ";
         std::cout << standard_env_dir << "\n";
         std::cout << "         is not set.";
         std::cout << " Mutations will not be possible\n";
         // mark as not read then:
         standard_residues_asc.read_success = 0;
         standard_residues_asc.n_selected_atoms = 0;
         // std::cout << "DEBUG:: standard_residues_asc marked as empty" << std::endl;

      } else {

#if 0
         // unresolved (linking related?) startup bug here:
         std::cout << "------------------ read_standard_residues() C map_sampling_rate " << map_sampling_rate << std::endl;
         std::cout << "------------------ read_standard_residues() C " << std::endl;
         atom_selection_container_t t_asc = get_atom_selection(standard_file_name, true, true, false);
         std::cout << "------------------ read_standard_residues() D map_sampling_rate " << map_sampling_rate << std::endl;
         // std::cout << "------------------ read_standard_residues() D " << std::endl;
         // standard_residues_asc = t_asc; // Here's the problem
#else

         mmdb::ERROR_CODE err;
         mmdb::Manager *mol = new mmdb::Manager;
         err = mol->ReadCoorFile(standard_file_name.c_str());
         if (err) {
            std::cout << "There was an error reading " << standard_file_name << ". \n";
            std::cout << "ERROR " << err << " READ: " << mmdb::GetErrorDescription(err) << std::endl;
            delete mol;
         } else {

            mmdb::PPAtom atom_selection = 0;
            int n_selected_atoms = 0;
            int SelHnd = mol->NewSelection(); // d
            mol->SelectAtoms(SelHnd, 1, "*",
                             mmdb::ANY_RES, "*",
                             mmdb::ANY_RES, "*",
                             "*","*","!H","*", mmdb::SKEY_NEW);
            standard_residues_asc.mol              = mol;
            standard_residues_asc.n_selected_atoms = n_selected_atoms;
            standard_residues_asc.atom_selection   = atom_selection;
            standard_residues_asc.read_success     = 1;
            standard_residues_asc.SelectionHandle  = SelHnd;

         }
#endif

      }
   } else {
      bool use_gemmi = true;
      standard_residues_asc = get_atom_selection(env_var_filename, use_gemmi, true, false);
   }

}



//! get the active atom
std::pair<int, std::string>
molecules_container_t::get_active_atom(float x, float y, float z, const std::string &displayed_model_molecules_list) const {

   auto atom_to_cid = [] (mmdb::Atom *at) {
      if (at) {
         std::string s = "/";
         s += std::to_string(at->GetModelNum());
         s += "/";
         s += std::string(at->GetChainID());
         s += "/";
         s += std::to_string(at->GetSeqNum());
         s += std::string(at->GetInsCode());
         s += "/";
         s += std::string(at->GetAtomName());
         std::string a(at->altLoc);
         if (! a.empty()) {
            s += ":";
            s += std::string();
         }
         return s;
      } else {
         return std::string("");
      }
   };

   int imol_closest = -1;
   std::string cid;
   std::vector<std::string> number_strings = coot::util::split_string(displayed_model_molecules_list, ":");
   std::vector<int> mols;
   for (const auto &item : number_strings) {
      int idx = coot::util::string_to_int(item);
      if (is_valid_model_molecule(idx))
         mols.push_back(idx);
   }

   float best_distance_sqrd = 99999999999999999.9;
   int best_imol = -1;
   mmdb::Atom *best_atom = 0;
   coot::Cartesian screen_centre(x,y,z);
   for (unsigned int ii=0; ii<mols.size(); ii++) {
      int imol = mols[ii];
      mmdb::Manager *mol = molecules[imol].atom_sel.mol;
      if (mol) {
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
                              coot::Cartesian atom_pos(at->x, at->y, at->z);
                              float dd = coot::Cartesian::lengthsq(screen_centre, atom_pos);
                              if (dd < best_distance_sqrd) {
                                 best_distance_sqrd = dd;
                                 best_imol = imol;
                                 best_atom = at;
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   if (best_atom) {
      imol_closest = best_imol;
      cid = atom_to_cid(best_atom);
   }
   return std::make_pair(imol_closest, cid);
}


void
molecules_container_t::update_updating_maps(int imol) {

   // should I return the stats?

   if (updating_maps_info.imol_model == imol) {
      if (updating_maps_info.maps_need_an_update) {
         if (is_valid_model_molecule(imol)) {
            if (is_valid_map_molecule(updating_maps_info.imol_2fofc)) {
               if (is_valid_map_molecule(updating_maps_info.imol_fofc)) {
                  coot::util::sfcalc_genmap_stats_t stats =
                     sfcalc_genmaps_using_bulk_solvent(imol,
                                                       updating_maps_info.imol_2fofc,
                                                       updating_maps_info.imol_fofc,
                                                       updating_maps_info.imol_with_data_info_attached);
                  // sfcalc_genmaps_using_bulk_solvent() sets latest_sfcalc_stats
                  updating_maps_info.maps_need_an_update = false;
               }
            }
         }
      } else {
         // 20221106-PE add debugging for now
         std::cout << "in updating_maps_info() maps_need_an_update is false" << std::endl;
      }
   }
}

void
molecules_container_t::set_updating_maps_need_an_update(int imol) {

   if (updating_maps_info.imol_model == imol)
      updating_maps_info.maps_need_an_update = true;

}

molecules_container_t::r_factor_stats
molecules_container_t::get_r_factor_stats() {

   int rpn_8 = calculate_new_rail_points();
   int rpt_8 = rail_points_total();
   auto latest_r_factors = get_latest_sfcalc_stats();

   r_factor_stats stats;
   stats.r_factor = latest_r_factors.r_factor;
   stats.free_r_factor = latest_r_factors.free_r_factor;
   stats.rail_points_total = rpt_8;
   stats.rail_points_new   = rpn_8;

   // std::cout << ":::::: get_r_factor_stats() " << r_factor_stats_as_string(stats) << std::endl;
   return stats;

}

std::string
molecules_container_t::r_factor_stats_as_string(const molecules_container_t::r_factor_stats &rfs) const {

   std::string s;
   s += "R-factor " + std::to_string(rfs.r_factor);
   s += " Free-R-factor " + std::to_string(rfs.free_r_factor);
   s += " Moorhen-Points-Total  " + std::to_string(rfs.rail_points_total);
   s += " Moorhen-Points-New  " + std::to_string(rfs.rail_points_new);
   return s;
}


coot::atom_spec_t
molecules_container_t::atom_cid_to_atom_spec(int imol, const std::string &cid) const {

   coot::atom_spec_t spec;
   if (is_valid_model_molecule(imol)) {
      auto p = molecules[imol].cid_to_atom_spec(cid);
      if (p.first) {
         spec = p.second;
      } else {
         std::cout << "WARNING:: molecule_class_info_t::atom_cid_to_atom_spec() no matching atom " << cid << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return spec;
}


coot::residue_spec_t
molecules_container_t::residue_cid_to_residue_spec(int imol, const std::string &cid) const   {

   coot::residue_spec_t spec;

   if (is_valid_model_molecule(imol)) {
      auto p = molecules[imol].cid_to_residue_spec(cid);
      if (p.first) {
         spec = p.second;
      } else {
         std::cout << "WARNING:: molecule_class_info_t::residue_cid_to_residue_spec() no matching residue " << cid << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return spec;
}


void
molecules_container_t::geometry_init_standard() {
   geom.init_standard();
}


int
molecules_container_t::undo(int imol) {
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].undo();
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

int
molecules_container_t::redo(int imol) {
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].redo();
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}




int
molecules_container_t::flip_peptide(int imol, const coot::atom_spec_t &as, const std::string &alt_conf) {

   int result = 0;
   if (is_valid_model_molecule(imol)) {
      result = molecules[imol].flip_peptide(as, alt_conf);
      set_updating_maps_need_an_update(imol);
   }
   return result;
}

int
molecules_container_t::flip_peptide_using_cid(int imol, const std::string &atom_cid, const std::string &alt_conf) {

   int result = 0;
   if (is_valid_model_molecule(imol)) {
      auto &m = molecules[imol];
      std::pair<bool, coot::atom_spec_t> as = m.cid_to_atom_spec(atom_cid);
      if (as.first) {
         const auto &atom_spec = as.second;
         result = molecules[imol].flip_peptide(atom_spec, alt_conf); // N check in here
         set_updating_maps_need_an_update(imol);
      }
   }
   return result;
}

int
molecules_container_t::install_model(const coot::molecule_t &m) {

   int size = molecules.size();
   molecules.push_back(m);
   return size;
}


//! read a coordinates file (mmcif or PDB)
//! @return the new molecule index on success and -1 on failure
int
molecules_container_t::read_coordinates(const std::string &file_name) {

   auto print_the_SSE = [] (mmdb::Manager *mol) {

      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               mmdb::byte sse = residue_p->SSE;
               int sse_int = sse;
               std::string rn = residue_p->GetResName();
               std::cout << "   " << coot::residue_spec_t(residue_p) << " " << rn << " " << sse_int << std::endl;
            }
         }
      }
   };

   auto debug_model = [] (mmdb::Model *model_p, const std::string &tag) {

      int n_atoms = 0;
      if (model_p) {
         std::cout << tag << " model_p: " << model_p << std::endl;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            std::cout << tag << "   chain_p: " << chain_p << std::endl;
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               std::cout << tag << "      residue_p: " << residue_p << std::endl;
               if (residue_p) {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     if (at) {
                        std::cout << tag << "         at: " << at << " " << coot::atom_spec_t(at) << std::endl;
                        if (! at->isTer()) {
                           n_atoms++;
                        }
                     } else {
                        std::cout << tag << "         at " << " was null " <<  iat
                                  << " of " << n_atoms << std::endl;
                     }
                  }
               }
            }
         }
      }
      std::cout << "n_atoms in model_p " << model_p << " " << n_atoms << std::endl;
   };

   int status = -1;
   atom_selection_container_t asc = get_atom_selection(file_name, use_gemmi, true, false);
   if (asc.read_success) {

      if (false) { // debugging the model
         mmdb::Manager *mol = new mmdb::Manager;
         mol->ReadPDBASCII(file_name.c_str());
         for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
            mmdb::Model *model_p = mol->GetModel(imod);
            std::string tag = "read_coordinates_molecule_MODEL_" + std::to_string(imod);
            debug_model(model_p, tag);
         }
      }

      // print_the_SSE(asc.mol); make the function print_the_SSE() available in the API. It may be
      // useful there.

      // 20221011-PE this constructor doesn't call make_bonds().
      int imol = molecules.size();
      coot::molecule_t m = coot::molecule_t(asc, imol, file_name);
      // m.make_bonds(&geom, &rot_prob_tables); // where does this go? Here or as a container function?
      molecules.push_back(m);
      status = imol;
   } else {
      logger.log(log_t::DEBUG, logging::function_name_t(__FUNCTION__),
		 { "asc.read_success was", asc.read_success,
		      "for", file_name});
   }
   return status;
}


int
molecules_container_t::read_pdb(const std::string &file_name) {

   return read_coordinates(file_name);
}

//! Read a Small molecule CIF file
//!
//! @param file_name is the cif file-name
//!
//! @return the new molecule index on success and -1 on failure
int
molecules_container_t::read_small_molecule_cif(const std::string &file_name) {

   int imol = -1;

   coot::smcif cif;
   mmdb::Manager *mol = cif.read_sm_cif(file_name);
   if (mol) {
      imol = molecules.size();
      atom_selection_container_t asc = make_asc(mol);
      molecules.push_back(coot::molecule_t(asc, imol, file_name));
   }
   return imol;
}



//! read a PDB file (or mmcif coordinates file, despite the name) to
//! replace the current molecule. This will only work if the molecules
//! is already a model molecule
void
molecules_container_t::replace_molecule_by_model_from_file(int imol, const std::string &pdb_file_name) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].replace_molecule_by_model_from_file(pdb_file_name);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}


int
molecules_container_t::import_cif_dictionary(const std::string &cif_file_name, int imol_enc) {

   coot::read_refmac_mon_lib_info_t r = geom.init_refmac_mon_lib(cif_file_name,
                                                                 cif_dictionary_read_number, imol_enc);
   cif_dictionary_read_number++;

   if (false)
      std::cout << "debug:: import_cif_dictionary() cif_file_name: " << cif_file_name
                << " for imol_enc " << imol_enc << " success " << r.success << " with "
                << r.n_atoms << " atoms " << r.n_bonds << " bonds " << r.n_links << " links"
                << " and monomer_idx " << r.monomer_idx << std::endl;

   if (false)
      geom.print_dictionary_store();

   return r.success;

}


int
molecules_container_t::get_monomer_from_dictionary(const std::string &comp_id,
                                                   int imol_enc,
                                                   bool idealised_flag) {

   int istat = -1; // unfound molecule

   // int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
   mmdb::Manager *mol = geom.mol_from_dictionary(comp_id, imol_enc, idealised_flag);
   if (mol) {
      int imol = molecules.size();
      std::string name = comp_id;
      name += "_from_dict";
      // graphics_info_t::molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1);
      // move_molecule_to_screen_centre_internal(imol);
      atom_selection_container_t asc = make_asc(mol);
      coot::molecule_t m = coot::molecule_t(asc, imol, name);
      molecules.push_back(m);
      istat = imol;
   } else {
      std::cout << "WARNING:: Null mol from mol_from_dictionary() with comp_id " << comp_id << " "
                << idealised_flag << std::endl;
   }
   return istat;
}


int
molecules_container_t::get_monomer(const std::string &comp_id) {

   int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
   int imol = get_monomer_from_dictionary(comp_id, imol_enc, true); // idealized
   return imol;
}

int
molecules_container_t::get_monomer_and_position_at(const std::string &comp_id, int imol_enc, float x, float y, float z) {

   int imol = get_monomer_from_dictionary(comp_id, imol_enc, true);
   if (is_valid_model_molecule(imol)) {
      move_molecule_to_new_centre(imol, x, y, z);
   }
   return imol;
}


//! return the group for the give list of residue names
std::vector<std::string>
molecules_container_t::get_groups_for_monomers(const std::vector<std::string> &residue_names) const {

   std::vector<std::string> v;
   std::vector<std::string>::const_iterator it;
   for (it=residue_names.begin(); it!=residue_names.end(); ++it) {
      v.push_back(geom.get_group(*it));
   }
   return v;
}

//! return the group for the give residue name
std::string
molecules_container_t::get_group_for_monomer(const std::string &residue_name) const {

   std::string s = geom.get_group(residue_name);
   return s;
}


#if RDKIT_HAS_CAIRO_SUPPORT // 20231129-PE Cairo is not allowed in Moorhen.
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>
#include "lidia-core/rdkit-interface.hh"
#endif

//! write a PNG for the given compound_id
void
molecules_container_t::write_png(const std::string &compound_id, int imol_enc,
                                 const std::string &file_name) const {

#if RDKIT_HAS_CAIRO_SUPPORT // 20231129-PE Cairo is not allowed in Moorhen.
                            // 20231221-PE but is in Coot/libcootapi/chapi

   // For now, let's use RDKit PNG depiction, not lidia-core/pyrogen

   std::pair<bool, coot::dictionary_residue_restraints_t> r_p =
      geom.get_monomer_restraints(compound_id, imol_enc);

   if (r_p.first) {
      const auto &restraints = r_p.second;
      std::pair<int, RDKit::RWMol> mol_pair = coot::rdkit_mol_with_2d_depiction(restraints);
      int conf_id = mol_pair.first;
      if (conf_id >= 0) {
         const auto &rdkit_mol(mol_pair.second);
         RDKit::MolDraw2DCairo drawer(500, 500);
         drawer.drawMolecule(rdkit_mol, nullptr, nullptr, nullptr, conf_id);
         drawer.finishDrawing();
         std::string dt = drawer.getDrawingText();
         std::ofstream f(file_name.c_str());
         f << dt;
         f << "\n";
         f.close();
      }
   }
#endif
}



// 20221030-PE nice to have one day
// int
// molecules_container_t::get_monomer_molecule_by_network(const std::string &text) {

//    int imol = -1;
//    return imol;
// }


int
molecules_container_t::write_coordinates(int imol, const std::string &file_name) const {

   if (false) {
      mmdb::Manager *mol = get_mol(imol);
      mol->WriteCIFASCII("write_coords_molecules_container_fn_start.cif");
   }

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].write_coordinates(file_name);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}


int
molecules_container_t::read_mtz(const std::string &file_name,
                                const std::string &f, const std::string &phi, const std::string &weight,
                                bool use_weight, bool is_a_difference_map) {

   int imol = -1; // currently "failure"
   int imol_in_hope = molecules.size();

   std::string name_in = file_name + std::string(" ") + std::string(f) + std::string(" ") + std::string(phi);
   coot::molecule_t m(name_in, imol_in_hope);

   // 20230417-PE hack!
   // map_sampling_rate = 1.8; // because init problems. Compiler config related?
   // std::cout << "........... read_mtz() map_sampling_rate " << map_sampling_rate << std::endl;

   bool status = coot::util::map_fill_from_mtz(&m.xmap, file_name, f, phi, weight, use_weight, map_sampling_rate);
   if (is_a_difference_map)
      m.set_map_is_difference_map(true);
   if (status) {
      molecules.push_back(m);
      imol = imol_in_hope;
      if (false)
         std::cout << "DEBUG:: in read_mtz() " << file_name << " " << f << " " << phi << " imol map: " << imol
                   << " diff-map-status: " << is_a_difference_map << std::endl;
   }
   return imol;
}

int
molecules_container_t::replace_map_by_mtz_from_file(int imol,
                                                    const std::string &file_name, const std::string &f, const std::string &phi,
                                                    const std::string &weight, bool use_weight) {

   int status = 0;
   if (is_valid_map_molecule(imol)) {
      clipper::Xmap<float> &xmap = molecules[imol].xmap;
      status = coot::util::map_fill_from_mtz(&xmap, file_name, f, phi, weight, use_weight, map_sampling_rate);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;

}



#include "coot-utils/cmtz-interface.hh"
#include "coot-utils/mtz-column-auto-read.hh"

// utlity function
//
// We should allow labels that are simply "FWT" and "PHWT" without
// dataset and xtal info.
int
molecules_container_t::valid_labels(const std::string &mtz_file_name, const std::string &f_col, const std::string &phi_col,
                                    const std::string &weight_col, int use_weights) const {

   int valid = 0;

   short int have_f = 0;
   short int have_phi = 0;
   short int have_weight = 1; // later turn on test if we have weights.

   std::string f_col_str(f_col);
   std::string phi_col_str(phi_col);
   std::string weight_col_str("");

   if (use_weights)
      weight_col_str = weight_col;

   // These now return have 0 members on failure
   //
//    char **f_cols      = get_f_cols(mtz_file_name, &n_f);
//    char **phi_cols    = get_phi_cols(mtz_file_name, &n_phi);
//    char **weight_cols = get_weight_cols(mtz_file_name, &n_weight);
//    char **d_cols      = get_d_cols(mtz_file_name, &n_d); // anom

   coot::mtz_column_types_info_t r = coot::get_mtz_columns(mtz_file_name);

   // Check first the MTZ column labels that don't have a slash
   for (unsigned int i=0; i<r.f_cols.size(); i++) {
      std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(r.f_cols[i].column_label);
      if (p.second.length() > 0)
         if (p.second == f_col_str) {
            have_f = 1;
            break;
         }
   }
   for (unsigned int i=0; i<r.phi_cols.size(); i++) {
      std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(r.phi_cols[i].column_label);
      if (p.second.length() > 0)
	      if (p.second == phi_col_str) {
	         have_phi = 1;
	         break;
	      }
   }
   if (use_weights) {
      for (unsigned int i=0; i<r.weight_cols.size(); i++) {
	      std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(r.weight_cols[i].column_label);
	      if (p.second.length() > 0)
	         if (p.second == weight_col_str) {
	            have_weight = 1;
	            break;
	         }
      }
   }

   // Now check the MTZ column labels that *do* have a slash
   if (r.f_cols.size() > 0) {
      for (unsigned int i=0; i< r.f_cols.size(); i++) {
	      if (f_col_str == r.f_cols[i].column_label) {
	         have_f = 1;
	         break;
	      }
      }
   } else {
      std::cout << "ERROR: no f_cols! " << std::endl;
   }

   // We can be trying to make an anomalous fourier.
   if (! have_f) {
      if (r.d_cols.size() > 0) {
	 for (unsigned int i=0; i< r.d_cols.size(); i++) {
	    std::cout << "comparing " << f_col_str << " " << r.d_cols[i].column_label << std::endl;
	    if (f_col_str == r.d_cols[i].column_label) {
	       have_f = 1;
	       break;
	    }
	    std::pair<std::string, std::string> p =
	       coot::util::split_string_on_last_slash(r.d_cols[i].column_label);
	    if (p.second.length() > 0) {
	       if (f_col_str == p.second) {
		  have_f = 1;
		  break;
	       }
	    }
	 }
      }
   }

   if (r.phi_cols.size() > 0) {
      for (unsigned int i=0; i< r.phi_cols.size(); i++) {
	 if (phi_col_str == r.phi_cols[i].column_label) {
	    have_phi = 1;
	    break;
	 }
      }
   } else {
      std::cout << "ERROR: no phi_cols! " << std::endl;
   }

   if (use_weights) {
      have_weight = 0;
      weight_col_str = std::string(weight_col);
      if (r.weight_cols.size() > 0) {
	 for (unsigned int i=0; i< r.weight_cols.size(); i++) {
	    if (weight_col_str == r.weight_cols[i].column_label) {
	       have_weight = 1;
	       break;
	    }
	 }
      } else {
	 std::cout << "ERROR: bad (null) weight_cols! " << std::endl;
      }
   }

   if (have_f && have_phi && have_weight)
      valid = 1;

   if (false)  // debug
      std::cout << "DEBUG:: done checking for valid column labels... returning "
		<< valid << " have-f: " << have_f << " have_phi: " << have_phi << " "
		<< "have_weight: " << have_weight << std::endl;
   return valid;
}


//! Read the given mtz file.
//! @return a vector of the maps created from reading the file
std::vector<molecules_container_t::auto_read_mtz_info_t>
molecules_container_t::auto_read_mtz(const std::string &mtz_file_name) {

   std::vector<molecules_container_t::auto_read_mtz_info_t> mol_infos;

   std::vector<coot::mtz_column_trials_info_t> auto_mtz_pairs;

   // Built-ins
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("FWT",          "PHWT",      false));
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("2FOFCWT",      "PH2FOFCWT", false));
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("DELFWT",       "PHDELWT",   true ));
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("FOFCWT",       "PHFOFCWT",  true ));
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("FDM",          "PHIDM",     false));
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("FAN",          "PHAN",      true));
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("F_ano",        "PHI_ano",   true));
   auto_mtz_pairs.push_back(coot::mtz_column_trials_info_t("F_early-late", "PHI_early-late", true));

   auto add_r_free_column_label = [] (auto_read_mtz_info_t *a, const coot::mtz_column_types_info_t &r) {
      for (unsigned int i=0; i<r.r_free_cols.size(); i++) {
         const std::string &l = r.r_free_cols[i].column_label;
         if (! l.empty()) {
            a->Rfree = l;
            break;
         }
      }
   };

   auto add_Fobs = [] (auto_read_mtz_info_t *armi_p, const auto_read_mtz_info_t &aFobs) {
      if (!aFobs.F_obs.empty()) {
         armi_p->F_obs = aFobs.F_obs;
         armi_p->sigF_obs = aFobs.sigF_obs;
         armi_p->Rfree = aFobs.Rfree;
      }
   };

   // 20221001-PE if there is one F and one PHI col, read that also (and it is not a difference map)
   coot::mtz_column_types_info_t r = coot::get_mtz_columns(mtz_file_name);

   // and now the observed data, this relies on the column label being of the form
   // /crystal/dataset/label
   //
   auto_read_mtz_info_t armi_fobs;
   for (unsigned int i=0; i<r.f_cols.size(); i++) {
      const std::string &f = r.f_cols[i].column_label;
      // example f: "/2vtq/1/FP"
      std::string  nd_f = coot::util::file_name_non_directory(f);
      std::string dir_f = coot::util::file_name_directory(f);
      for (unsigned int j=0; j<r.sigf_cols.size(); j++) {
         const std::string &sf = r.sigf_cols[j].column_label;
         std::string test_string = std::string(dir_f + std::string("SIG") + nd_f);
         if (sf == test_string) {
            auto_read_mtz_info_t armi;
            armi.set_fobs_sigfobs(f, sf);
            add_r_free_column_label(&armi, r); // modify armi possibly
            mol_infos.push_back(armi);
            if (armi_fobs.F_obs.empty())
               armi_fobs = armi;
         }
      }
   }


   for (unsigned int i=0; i<auto_mtz_pairs.size(); i++) {
      const coot::mtz_column_trials_info_t &b = auto_mtz_pairs[i];
      if (valid_labels(mtz_file_name.c_str(), b.f_col.c_str(), b.phi_col.c_str(), "", 0)) {
         int imol = read_mtz(mtz_file_name, b.f_col, b.phi_col, "", 0, b.is_diff_map);
	      if (is_valid_map_molecule(imol)) {
            auto_read_mtz_info_t armi(imol, b.f_col, b.phi_col);
            add_Fobs(&armi, armi_fobs);
	         mol_infos.push_back(armi);
         }
      }
   }

   if (r.f_cols.size() == 1) {
      if (r.phi_cols.size() == 1) {
         int imol = read_mtz(mtz_file_name, r.f_cols[0].column_label, r.phi_cols[0].column_label, "", false, false);
         auto_read_mtz_info_t armi(imol, r.f_cols[0].column_label, r.phi_cols[0].column_label);
         add_Fobs(&armi, armi_fobs);
         mol_infos.push_back(auto_read_mtz_info_t(armi));
      }
   }

   for (unsigned int i=0; i<r.f_cols.size(); i++) {
      std::string s = r.f_cols[i].column_label;
      std::string::size_type idx = s.find(".F_phi.F");
      if (idx != std::string::npos) {
	      std::string prefix = s.substr(0, idx);
	      std::string trial_phi_col = prefix + ".F_phi.phi";
	      for (unsigned int j=0; j<r.phi_cols.size(); j++) {
	         if (r.phi_cols[j].column_label == trial_phi_col) {
	            std::string f_col   = r.f_cols[i].column_label;
	            std::string phi_col = r.phi_cols[j].column_label;
               int imol = read_mtz(mtz_file_name, f_col, phi_col, "", false, false);
               if (is_valid_map_molecule(imol)) {
                  auto_read_mtz_info_t armi(imol, f_col, phi_col);
                  add_Fobs(&armi, armi_fobs);
                  mol_infos.push_back(armi);
               }
	         }
	      }
      }
   }

   return mol_infos;
}

#include "clipper-ccp4-map-file-wrapper.hh"
#include "coot-utils/slurp-map.hh"



//! rotamer validation information
coot::validation_information_t
molecules_container_t::rotamer_analysis(int imol_model) const {

   coot::validation_information_t r;
   r.name = "Rotamer analysis";
#ifdef EMSCRIPTEN
   r.type = "PROBABILITY";
#else
   r.type = coot::PROBABILITY;
#endif

   if (is_valid_model_molecule(imol_model)) {

      mmdb::Manager *mol = molecules[imol_model].atom_sel.mol;

      // fill these
      mmdb::PResidue *SelResidues = 0;
      int nSelResidues = 0;

      int selHnd = mol->NewSelection(); // yes, it's deleted.
      int imod = 1; // multiple models don't work on validation graphs

      mol->Select(selHnd, mmdb::STYPE_RESIDUE, imod,
                           "*", // chain_id
                           mmdb::ANY_RES, "*",
                           mmdb::ANY_RES, "*",
                           "*",  // residue name
                           "*",  // Residue must contain this atom name?
                           "*",  // Residue must contain this Element?
                           "*",  // altLocs
                           mmdb::SKEY_NEW // selection key
                           );
      mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

      for (int ir=0; ir<nSelResidues; ir++) {
         mmdb::Residue *residue_p = SelResidues[ir];
         coot::residue_spec_t res_spec(residue_p);
         mmdb::PAtom *residue_atoms=0;
         int n_residue_atoms;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

         // double residue_density_score = coot::util::map_score(residue_atoms, n_residue_atoms, xmap, 1);

         if (n_residue_atoms > 5) {

            std::string res_name = residue_p->GetResName();
            if (true) {

               coot::rotamer rot(residue_p);
               coot::rotamer_probability_info_t rpi = rot.probability_of_this_rotamer();
               double prob = rpi.probability;

               std::string l = res_spec.label();
               std::string atom_name = coot::util::intelligent_this_residue_mmdb_atom(residue_p)->GetAtomName();
               const std::string &chain_id = res_spec.chain_id;
               int this_resno = res_spec.res_no;
               coot::atom_spec_t atom_spec(chain_id, this_resno, res_spec.ins_code, atom_name, "");
               coot::residue_validation_information_t rvi(res_spec, atom_spec, prob, l);
               r.add_residue_validation_information(rvi, chain_id);
            }
         }
      }
      mol->DeleteSelection(selHnd);
   }
   r.set_min_max();
   return r;
}

double
molecules_container_t::phi_psi_probability(const coot::util::phi_psi_t &phi_psi, const ramachandrans_container_t &rc) const {

      const clipper::Ramachandran *rama = &rc.rama;

      if (phi_psi.residue_name() == "PRO") rama = &rc.rama_pro;
      if (phi_psi.residue_name() == "GLY") rama = &rc.rama_gly;

      // if (phi_psi.residue_name() == "ILE" || phi_psi.residue_name() == "VAL" ) rama = &rc.rama_ileval;
      // if (phi_psi.is_pre_pro())
      // if (phi_psi.residue_name() != "GLY")
      // rama = &rc.rama_pre_pro;

      double rama_prob = rama->probability(clipper::Util::d2rad(phi_psi.phi()),
                                           clipper::Util::d2rad(phi_psi.psi()));
      return rama_prob;
}

//! ramachandran validation information (formatted for a graph, not 3d)
coot::validation_information_t
molecules_container_t::ramachandran_analysis(int imol_model) const {

   coot::validation_information_t vi;
   vi.name = "Ramachandran plot Probability";
#ifdef EMSCRIPTEN
   vi.type = "PROBABILITY";
#else
   vi.type = coot::PROBABILITY;
#endif
   std::vector<coot::phi_psi_prob_t> rv = ramachandran_validation(imol_model);

   for (unsigned int i=0; i<rv.size(); i++) {
      std::string chain_id = rv[i].phi_psi.chain_id;
      coot::residue_spec_t residue_spec(rv[i].phi_psi.chain_id, rv[i].phi_psi.residue_number, rv[i].phi_psi.ins_code);
      double pr = rv[i].probability;
      std::string label = rv[i].phi_psi.chain_id + std::string(" ") + std::to_string(rv[i].phi_psi.residue_number);
      if (! rv[i].phi_psi.ins_code.empty())
         label += std::string(" ") + rv[i].phi_psi.ins_code;
      coot::atom_spec_t atom_spec(residue_spec.chain_id, residue_spec.res_no, residue_spec.ins_code, " CA ", "");
      coot::residue_validation_information_t rvi(residue_spec, atom_spec, pr, label);
      if (false)
         std::cout << "         " << residue_spec << " " << rv[i].phi_psi.phi() << " " << rv[i].phi_psi.psi()
                   << " pr " << pr << " " << std::endl;
      vi.add_residue_validation_information(rvi, chain_id);
   }
   vi.set_min_max();
   return vi;
}

//! ramachandran validation information (formatted for a graph, not 3d) for a given chain in a given molecule
//! 20230127-PE This function does not exist yet.
//!
//! @returns a `coot::validation_information_t`
coot::validation_information_t
molecules_container_t::ramachandran_analysis_for_chain(int imol_model, const std::string &user_chain_id) const {

   coot::validation_information_t vi;
   vi.name = "Ramachandran plot Probability";
#ifdef EMSCRIPTEN
   vi.type = "PROBABILITY";
#else
   vi.type = coot::PROBABILITY;
#endif
   std::vector<coot::phi_psi_prob_t> rv = ramachandran_validation(imol_model);

   for (unsigned int i=0; i<rv.size(); i++) {
      std::string chain_id = rv[i].phi_psi.chain_id;
      if (chain_id != user_chain_id) continue;
      coot::residue_spec_t residue_spec(rv[i].phi_psi.chain_id, rv[i].phi_psi.residue_number, rv[i].phi_psi.ins_code);
      double pr = rv[i].probability;
      std::string label = rv[i].phi_psi.chain_id + std::string(" ") + std::to_string(rv[i].phi_psi.residue_number);
      if (! rv[i].phi_psi.ins_code.empty())
         label += std::string(" ") + rv[i].phi_psi.ins_code;
      coot::atom_spec_t atom_spec(residue_spec.chain_id, residue_spec.res_no, residue_spec.ins_code, " CA ", "");
      coot::residue_validation_information_t rvi(residue_spec, atom_spec, pr, label);
      if (false)
         std::cout << "         " << residue_spec << " " << rv[i].phi_psi.phi() << " " << rv[i].phi_psi.psi()
                   << " pr " << pr << " " << std::endl;
      vi.add_residue_validation_information(rvi, chain_id);
   }
   vi.set_min_max();
   return vi;
}


//! peptide omega validation information
//! @returns a `validation_information_t`
coot::validation_information_t
molecules_container_t::peptide_omega_analysis(int imol) const {

   coot::validation_information_t vi;
   vi.name = "Peptide Omega Deviation";
#ifdef EMSCRIPTEN
   vi.type = "TORSION_ANGLE";
#else
   vi.type = coot::TORSION_ANGLE;
#endif

   if (is_valid_model_molecule(imol)) {

      bool mark_cis_peptides_as_bad_flag = false;
      bool m = mark_cis_peptides_as_bad_flag;
      std::vector<std::string> chain_ids = molecules[imol].chains_in_model();
      for (const auto &chain_id : chain_ids) {
         coot::chain_validation_information_t cvi(chain_id);
         coot::omega_distortion_info_container_t odi = molecules.at(imol).peptide_omega_analysis(geom, chain_id, m);
         for (const auto &od : odi.omega_distortions) {
            // oops - we have forgotten about the insertion code.
            coot::residue_spec_t res_spec(chain_id, od.resno, "");
            coot::atom_spec_t atom_spec(chain_id, od.resno, "", " CA ", "");
            std::string label = od.info_string;
            coot::residue_validation_information_t rvi(res_spec, atom_spec, od.distortion, label);
            cvi.add_residue_validation_information(rvi);
         }
         vi.cviv.push_back(cvi);
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return vi;
}


// #include "vertex.hh" // neeeded?

coot::simple_mesh_t
molecules_container_t::test_origin_cube() const {

   coot::simple_mesh_t mesh;

   std::vector<coot::api::vnc_vertex> vertices;
   std::vector<g_triangle> triangles;

   glm::vec4 c(0.5, 0.2, 0.5, 1.0); // colour

   // bottom
   coot::api::vnc_vertex v0(glm::vec3(0, 0, 0), glm::vec3(0,0,-1), c); vertices.push_back(v0);
   coot::api::vnc_vertex v1(glm::vec3(1, 0, 0), glm::vec3(0,0,-1), c); vertices.push_back(v1);
   coot::api::vnc_vertex v2(glm::vec3(0, 1, 0), glm::vec3(0,0,-1), c); vertices.push_back(v2);
   coot::api::vnc_vertex v3(glm::vec3(1, 1, 0), glm::vec3(0,0,-1), c); vertices.push_back(v3);

   // top
   coot::api::vnc_vertex v4(glm::vec3(0, 0, 1), glm::vec3(0,0,1), c); vertices.push_back(v4);
   coot::api::vnc_vertex v5(glm::vec3(1, 0, 1), glm::vec3(0,0,1), c); vertices.push_back(v5);
   coot::api::vnc_vertex v6(glm::vec3(0, 1, 1), glm::vec3(0,0,1), c); vertices.push_back(v6);
   coot::api::vnc_vertex v7(glm::vec3(1, 1, 1), glm::vec3(0,0,1), c); vertices.push_back(v7);

   // left
   coot::api::vnc_vertex v8 (glm::vec3(0, 0, 0), glm::vec3(-1,0,0), c); vertices.push_back(v8);
   coot::api::vnc_vertex v9 (glm::vec3(0, 1, 0), glm::vec3(-1,0,0), c); vertices.push_back(v9);
   coot::api::vnc_vertex v10(glm::vec3(0, 0, 1), glm::vec3(-1,0,0), c); vertices.push_back(v10);
   coot::api::vnc_vertex v11(glm::vec3(0, 1, 1), glm::vec3(-1,0,0), c); vertices.push_back(v11);

   // right
   coot::api::vnc_vertex v12(glm::vec3(1, 0, 0), glm::vec3(1,0,0), c); vertices.push_back(v12);
   coot::api::vnc_vertex v13(glm::vec3(1, 1, 0), glm::vec3(1,0,0), c); vertices.push_back(v13);
   coot::api::vnc_vertex v14(glm::vec3(1, 0, 1), glm::vec3(1,0,0), c); vertices.push_back(v14);
   coot::api::vnc_vertex v15(glm::vec3(1, 1, 1), glm::vec3(1,0,0), c); vertices.push_back(v15);

   // front
   coot::api::vnc_vertex v16(glm::vec3(0, 0, 0), glm::vec3(0,-1,0), c); vertices.push_back(v16);
   coot::api::vnc_vertex v17(glm::vec3(1, 0, 0), glm::vec3(0,-1,0), c); vertices.push_back(v17);
   coot::api::vnc_vertex v18(glm::vec3(0, 0, 1), glm::vec3(0,-1,0), c); vertices.push_back(v18);
   coot::api::vnc_vertex v19(glm::vec3(1, 0, 1), glm::vec3(0,-1,0), c); vertices.push_back(v19);

   // back
   coot::api::vnc_vertex v20(glm::vec3(0, 1, 0), glm::vec3(0,1,0), c); vertices.push_back(v20);
   coot::api::vnc_vertex v21(glm::vec3(1, 1, 0), glm::vec3(0,1,0), c); vertices.push_back(v21);
   coot::api::vnc_vertex v22(glm::vec3(0, 1, 1), glm::vec3(0,1,0), c); vertices.push_back(v22);
   coot::api::vnc_vertex v23(glm::vec3(1, 1, 1), glm::vec3(0,1,0), c); vertices.push_back(v23);

   triangles.push_back(g_triangle( 0, 1, 2));
   triangles.push_back(g_triangle( 1, 3, 2));
   triangles.push_back(g_triangle( 4, 5, 6));
   triangles.push_back(g_triangle( 5, 7, 6));
   triangles.push_back(g_triangle( 8, 9,10));
   triangles.push_back(g_triangle( 9,11,10));
   triangles.push_back(g_triangle(12,13,14));
   triangles.push_back(g_triangle(13,15,14));
   triangles.push_back(g_triangle(16,17,18));
   triangles.push_back(g_triangle(17,19,18));
   triangles.push_back(g_triangle(20,21,22));
   triangles.push_back(g_triangle(21,23,22));

   for (auto &vertex : vertices) {
      vertex.pos *= 10.0;
      vertex.pos += glm::vec3(-5.0, -5.0, -5.0);
   }

   coot::simple_mesh_t m(vertices, triangles);
   // m.translate(glm::vec3(-0.5, -0.5, -0.5));
   return m;
}

std::vector<coot::phi_psi_prob_t>
molecules_container_t::ramachandran_validation(int imol) const {

   std::vector<coot::phi_psi_prob_t> v;
   if (is_valid_model_molecule(imol))
      v = molecules[imol].ramachandran_validation(ramachandrans_container);
   return v;
}

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

coot::simple_mesh_t
molecules_container_t::get_ramachandran_validation_markup_mesh(int imol) const {

   // this function should be pushed into the coot::molecule_t class
   // (which means that the mesh will be copied)

   unsigned int num_subdivisions = 2;  // pass this
   float rama_ball_radius = 0.5;

   auto prob_raw_to_colour_rotation = [] (float prob) {
                                         if (prob > 0.5) prob = 0.5; // 0.4 and 2.5 f(for q) might be better (not tested)
                                         // good probabilities have q = 0
                                         // bad probabilities have q 0.66
                                         double q = (1.0 - 2.0 * prob);
                                         q = pow(q, 20.0);
                                         return q;
   };

   auto cartesian_to_glm = [] (const coot::Cartesian &c) {
                              return glm::vec3(c.x(), c.y(), c.z());
   };

   auto phi_psi_probability = [] (const coot::util::phi_psi_t &phi_psi, const ramachandrans_container_t &rc) {

      const clipper::Ramachandran *rama = &rc.rama;

      if (phi_psi.residue_name() == "PRO") rama = &rc.rama_pro;
      if (phi_psi.residue_name() == "GLY") rama = &rc.rama_gly;

      // if (phi_psi.residue_name() == "ILE" || phi_psi.residue_name() == "VAL" ) rama = &rc.rama_ileval;
      // if (phi_psi.is_pre_pro())
      // if (phi_psi.residue_name() != "GLY")
      // rama = &rc.rama_pre_pro;

      double rama_prob = rama->probability(clipper::Util::d2rad(phi_psi.phi()),
                                           clipper::Util::d2rad(phi_psi.psi()));
      return rama_prob;
   };

   auto test_ramachandran_probabilities = [] (const ramachandrans_container_t &rc) {

      std::vector<const clipper::Ramachandran *> ramas = { &rc.rama, &rc.rama_gly, &rc.rama_pro, &rc.rama_non_gly_pro };

      for (unsigned int ir=0; ir<ramas.size(); ir++) {
         for (unsigned int i=0; i<10; i++) {
            for (unsigned int j=0; j<10; j++) {
               double phi = static_cast<double>(i * 36.0) - 180.0;
               double psi = static_cast<double>(j * 36.0) - 180.0;
               double p = rc.rama.probability(phi, psi);
               std::cout << ir << "   "
                         << std::setw(10) << phi << " " << std::setw(10) << psi << " "
                         << std::setw(10) << p << std::endl;
            }
         }
      }
   };

   // test_ramachandran_probabilities(ramachandrans_container); // don't use rama_pre_pro without CLIPPER_HAS_TOP8000

   coot::simple_mesh_t mesh;

   // 20221126-PE Calm down the ultra-bright rama balls:
   float sober_factor = 0.75;

   if (is_valid_model_molecule(imol)) {

      std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octaball = tessellate_octasphere(num_subdivisions);

      std::vector<coot::phi_psi_prob_t> ramachandran_goodness_spots = ramachandran_validation(imol);
      // now convert positions into meshes of balls
      int n_ramachandran_goodness_spots = ramachandran_goodness_spots.size();
      for (int i=0; i<n_ramachandran_goodness_spots; i++) {
         const coot::Cartesian &position = ramachandran_goodness_spots[i].position;
         // std::cout << "goodness spot " << i << " position " << position << std::endl;
         const coot::phi_psi_prob_t &phi_psi = ramachandran_goodness_spots[i];
         double prob_raw = phi_psi.probability;
         double q = prob_raw_to_colour_rotation(prob_raw);
         coot::colour_holder col = coot::colour_holder(q, 0.0, 1.0, false, std::string(""));
         glm::vec3 ball_position = cartesian_to_glm(position);
         unsigned int idx_base = mesh.vertices.size();
         unsigned int idx_tri_base = mesh.triangles.size();
         for (unsigned int ibv=0; ibv<octaball.first.size(); ibv++) {
            glm::vec4 col_v4(sober_factor * col.red, sober_factor * col.green, sober_factor * col.blue, 1.0f);
            const glm::vec3 &vertex_position = octaball.first[ibv];
            coot::api::vnc_vertex vertex(ball_position + rama_ball_radius * vertex_position, vertex_position, col_v4);
            mesh.vertices.push_back(vertex);
         }
         std::vector<g_triangle> octaball_triangles = octaball.second;
         mesh.triangles.insert(mesh.triangles.end(), octaball_triangles.begin(), octaball_triangles.end());

         for (unsigned int k=idx_tri_base; k<mesh.triangles.size(); k++)
            mesh.triangles[k].rebase(idx_base);
      }
   }
   return mesh;
}


mmdb::Atom *
molecules_container_t::get_atom(int imol, const coot::atom_spec_t &atom_spec) const {

   mmdb::Atom *r = nullptr;
   if (is_valid_model_molecule(imol)) {
      return molecules[imol].get_atom(atom_spec);
   }
   return r;
}

mmdb::Residue *
molecules_container_t::get_residue(int imol, const coot::residue_spec_t &residue_spec) const {

   mmdb::Residue *r = nullptr;
   if (is_valid_model_molecule(imol)) {
      return molecules[imol].get_residue(residue_spec);
   }
   return r;
}

// returns either the specified atom or null if not found
mmdb::Atom *
molecules_container_t::get_atom_using_cid(int imol, const std::string &cid) const {

   mmdb::Atom *at_p = nullptr;
   if (is_valid_model_molecule(imol)) {
      std::pair<bool, coot::atom_spec_t> p = molecules[imol].cid_to_atom_spec(cid);
      if (p.first) {
         mmdb::Atom *at_m = molecules[imol].get_atom(p.second);
         if (at_m) {
            at_p = new mmdb::Atom;
            at_p->Copy(at_m);
         }
      }
   }
   return at_p; // maybe a memory leak?
}

// returns either the specified residue or null if not found
mmdb::Residue *
molecules_container_t::get_residue_using_cid(int imol, const std::string &cid) const {

   auto deep_copy_residue_local = [] (mmdb::Residue *residue) {

      mmdb::Residue *rres = new mmdb::Residue;
      rres->SetResID(residue->GetResName(), residue->GetSeqNum(), residue->GetInsCode());

      mmdb::PPAtom residue_atoms;
      int nResidueAtoms;
      residue->GetAtomTable(residue_atoms, nResidueAtoms);
      for(int iat=0; iat<nResidueAtoms; iat++) {
         mmdb::Atom *atom_p = new mmdb::Atom;
         atom_p->Copy(residue_atoms[iat]);
         rres->AddAtom(atom_p);
      }
      return rres;
   };

   mmdb::Residue *residue_p = nullptr;

   if (is_valid_model_molecule(imol)) {
      std::pair<bool, coot::residue_spec_t> p = molecules[imol].cid_to_residue_spec(cid);
      if (p.first) {
         mmdb::Residue *rr = molecules[imol].get_residue(p.second);
         if (rr) {
            residue_p = deep_copy_residue_local(rr);
         }
      }
   }
   return residue_p;
}


//! get header info.
//! @return an object with header info. Sparce at the moment.
moorhen::header_info_t
molecules_container_t::get_header_info(int imol) const {

   auto get_author_info = [] (mmdb::Manager *mol) {

      std::vector<std::string> author_lines;
      access_mol *am = static_cast<access_mol *>(mol);
      const mmdb::Title *tt = am->GetTitle();
      mmdb::Title *ttmp = const_cast<mmdb::Title *>(tt);
      access_title *at = static_cast<access_title *> (ttmp);
      mmdb::TitleContainer *author_container = at->GetAuthor();
      unsigned int al = author_container->Length();
      for (unsigned int i=0; i<al; i++) {
         mmdb::Author *a_line = mmdb::PAuthor(author_container->GetContainerClass(i));
         if (a_line) {
            std::string line(a_line->Line);
            author_lines.push_back(line);
         }
      }
      return author_lines;
   };

   auto get_journal_info = [] (mmdb::Manager *mol) {

      std::vector<std::string> journal_lines;
      access_mol *am = static_cast<access_mol *>(mol);

      const mmdb::Title *tt = am->GetTitle();
      mmdb::Title *ttmp = const_cast<mmdb::Title *>(tt);
      access_title *at = static_cast<access_title *> (ttmp);
      mmdb::TitleContainer *journal_container = at->GetJournal();
      unsigned int al = journal_container->Length();
      for (unsigned int i=0; i<al; i++) {
         mmdb::Journal *j_line = mmdb::PJournal(journal_container->GetContainerClass(i));
         if (j_line) {
            std::string line(j_line->Line);
            journal_lines.push_back(line);
         }
      }
      return journal_lines;
   };

   bool screen_output = false;
   moorhen::header_info_t header;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = molecules[imol].atom_sel.mol;
      if (mol) {

         std::string title = coot::get_title(mol);
         std::vector<std::string> compound_lines = coot::get_compound_lines(mol);
         std::vector<std::string>   author_lines = get_author_info(mol);
         std::vector<std::string>  journal_lines = get_journal_info(mol);

         header.compound_lines = compound_lines;
         header.author_lines   = author_lines;
         header.journal_lines  = journal_lines;

         coot::secondary_structure_header_records sshr(mol, false);
         mmdb::Model *model_p = mol->GetModel(1);
         if (model_p) {
            coot::util::print_secondary_structure_info(model_p);
            int nhelix = model_p->GetNumberOfHelices();
            int nsheet = model_p->GetNumberOfSheets();
            std::cout << "INFO:: There are " << nhelix << " helices and " << nsheet << " sheets\n";
            for (int ih=1; ih<=nhelix; ih++) {
               mmdb:: Helix *helix_p = model_p->GetHelix(ih);
               if (helix_p) {
                  if (screen_output)
                     std::cout << helix_p->serNum      << " " << helix_p->helixID    << " "
                               << helix_p->initChainID << " " << helix_p->initSeqNum << " "
                               << helix_p->endChainID  << " " << helix_p->endSeqNum  << " "
                               << helix_p->length      << " " << helix_p->comment    << std::endl;
                  moorhen::helix_t helix(helix_p->serNum, helix_p->helixID,
                                         helix_p->initResName, helix_p->initChainID, helix_p->initSeqNum, helix_p->initICode,
                                         helix_p->endResName,  helix_p->endChainID,  helix_p->endSeqNum,  helix_p->endICode,
                                         helix_p->helixClass, helix_p->comment, helix_p->length);
                  header.helix_info.push_back(helix);
               } else {
                  std::cout << "ERROR: no helix!?" << std::endl;
               }
            }

            for (int isheet=0; isheet<nsheet; isheet++) {
               mmdb::Sheet *sheet_p = model_p->GetSheet(isheet);
                 if (sheet_p) {
                    int n_strand = sheet_p->nStrands;
                    for (int istrand=0; istrand<n_strand; istrand++) {
                       mmdb::Strand *strand_p = sheet_p->strand[istrand];
                       moorhen::strand_t strand(strand_p->strandNo,
                                                strand_p->initResName, strand_p->initChainID, strand_p->initSeqNum, strand_p->initICode,
                                                strand_p->endResName,  strand_p->endChainID,  strand_p->endSeqNum,  strand_p->endICode,
                                                strand_p->sense);
                       header.strand_info.push_back(strand);
                    }
                 }
            }
         }
      }
   }
   return header;
}



int
molecules_container_t::move_molecule_to_new_centre(int imol, float x, float y, float z) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::Cartesian new_centre(x,y,z);
      status = molecules[imol].move_molecule_to_new_centre(new_centre);
      set_updating_maps_need_an_update(imol); // weird thing to do usually
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

//! get the atom position - don't use this in emscript
std::pair<bool, coot::Cartesian>
molecules_container_t::get_atom_position(int imol, coot::atom_spec_t &atom_spec) {

   mmdb::Atom *at = get_atom(imol, atom_spec);
   if (at) {
      return std::pair<bool, coot::Cartesian> (true, coot::Cartesian(at->x, at->y, at->z));
   } else {
      return std::pair<bool, coot::Cartesian> (false, coot::Cartesian(0,0,0));
   }

}



coot::Cartesian
molecules_container_t::get_molecule_centre(int imol) const {

   coot::Cartesian c;
   if (is_valid_model_molecule(imol)) {
      c = molecules[imol].get_molecule_centre();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return c;
}


int
molecules_container_t::write_map(int imol, const std::string &file_name) const {

   int status= 0;
   if (is_valid_map_molecule(imol)) {
      status = molecules[imol].write_map(file_name);
   }
   return status;

}

// Mode is "COLOUR-BY-CHAIN-AND-DICTIONARY" or "CA+LIGANDS"
coot::simple_mesh_t
molecules_container_t::get_bonds_mesh(int imol, const std::string &mode,
                                      bool against_a_dark_background,
                                      float bonds_width, float atom_radius_to_bond_width_ratio,
                                      int smoothness_factor) {

   bool draw_hydrogen_atoms_flag = true; // pass this

   auto tp_0 = std::chrono::high_resolution_clock::now();

   coot::simple_mesh_t sm;
   if (is_valid_model_molecule(imol)) {
      sm = molecules[imol].get_bonds_mesh(mode, &geom, against_a_dark_background, bonds_width, atom_radius_to_bond_width_ratio,
                                          smoothness_factor, draw_hydrogen_atoms_flag, draw_missing_residue_loops_flag);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   auto tp_1 = std::chrono::high_resolution_clock::now();
   if (show_timings) {
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      std::cout << "---------- timings: for get_bonds_mesh(): : " << d10 << " milliseconds " << std::endl;
   }

   return sm;
}

void
molecules_container_t::add_to_non_drawn_bonds(int imol, const std::string &atom_selection_cid) {

   if (is_valid_model_molecule(imol)) {
       molecules[imol].add_to_non_drawn_bonds(atom_selection_cid);
   }
}

void
molecules_container_t::clear_non_drawn_bonds(int imol) {

   if (is_valid_model_molecule(imol)) {
       molecules[imol].clear_non_drawn_bonds();
   }
}


void
molecules_container_t::print_non_drawn_bonds(int imol) const {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].print_non_drawn_bonds();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}

//! @return an ``instanced_mesh_t``
coot::instanced_mesh_t
molecules_container_t::get_bonds_mesh_instanced(int imol, const std::string &mode,
                                                bool against_a_dark_background,
                                                float bond_width, float atom_radius_to_bond_width_ratio,
                                                bool show_atoms_as_aniso_flag,
                                                bool show_aniso_atoms_as_ortep_flag,
                                                bool show_aniso_atoms_as_empty_flag,
                                                bool draw_hydrogen_atoms_flag,
                                                int smoothness_factor) {

   float aniso_probability = 0.5f; // pass this in the api

   auto tp_0 = std::chrono::high_resolution_clock::now();

   coot::instanced_mesh_t im;
   if (is_valid_model_molecule(imol)) {

      // testing colours
      // set_use_bespoke_carbon_atom_colour(imol, true);
      // coot::colour_t col(0.0999, 0.0888, 0.0777);
      // set_bespoke_carbon_atom_colour(imol, col);
      im = molecules[imol].get_bonds_mesh_instanced(mode, &geom, against_a_dark_background, bond_width, atom_radius_to_bond_width_ratio,
                                                    show_atoms_as_aniso_flag,
                                                    aniso_probability,
                                                    show_aniso_atoms_as_ortep_flag,
                                                    show_aniso_atoms_as_empty_flag,
                                                    smoothness_factor, draw_hydrogen_atoms_flag,
                                                    draw_missing_residue_loops_flag);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   auto tp_1 = std::chrono::high_resolution_clock::now();
   if (show_timings) {
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      std::cout << "---------- timings: for get_bonds_mesh_instanced(): : " << d10 << " milliseconds " << std::endl;
   }

   return im;
}

//! As above, but only return the bonds for the atom selection.
//! Typically one would call this with a wider bond_with than one would use for standards atoms (all molecule)
//!
//! @return a ``coot::instanced_mesh_t``
coot::instanced_mesh_t
molecules_container_t::get_bonds_mesh_for_selection_instanced(int imol, const std::string &atom_selection_cid,
                                                              const std::string &mode,
                                                              bool against_a_dark_background,
                                                              float bond_width, float atom_radius_to_bond_width_ratio,
                                                              bool show_atoms_as_aniso_flag,
                                                              bool show_aniso_atoms_as_ortep_flag,
                                                              bool show_aniso_atoms_as_empty_flag,
                                                              bool draw_hydrogen_atoms_flag,
                                                              int smoothness_factor) {

   // auto tp_0 = std::chrono::high_resolution_clock::now();

   coot::instanced_mesh_t im;
   if (is_valid_model_molecule(imol)) {
      im = molecules[imol].get_bonds_mesh_for_selection_instanced(mode, atom_selection_cid,
                                                                  &geom, against_a_dark_background, bond_width, atom_radius_to_bond_width_ratio,
                                                                  show_atoms_as_aniso_flag,
                                                                  show_aniso_atoms_as_ortep_flag,
                                                                  show_aniso_atoms_as_empty_flag,
                                                                  smoothness_factor,
                                                                  draw_hydrogen_atoms_flag,
                                                                  draw_missing_residue_loops_flag);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return im;
}

coot::instanced_mesh_t
molecules_container_t::get_goodsell_style_mesh_instanced(int imol, float colour_wheel_rotation_step,
                                                         float saturation, float goodselliness) {

   coot::instanced_mesh_t im;
   if (is_valid_model_molecule(imol)) {
      im = molecules[imol].get_goodsell_style_mesh_instanced(&geom, colour_wheel_rotation_step, saturation, goodselliness);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return im;
}



//! return the colour table (for testing)
std::vector<glm::vec4>
molecules_container_t::get_colour_table(int imol, bool against_a_dark_background) const {

   std::vector<glm::vec4> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].make_colour_table(against_a_dark_background);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}



//! user-defined colour-index to colour
void
molecules_container_t::set_user_defined_bond_colours(int imol, const std::map<unsigned int, std::array<float, 4> > &colour_map) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].set_user_defined_bond_colours(colour_map);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}

//! user-defined atom selection to colour index
void
molecules_container_t::set_user_defined_atom_colour_by_selection(int imol, const std::vector<std::pair<std::string, unsigned int> > &indexed_residues_cids, bool colour_applies_to_non_carbon_atoms_also) {

   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = molecules[imol].atom_sel.mol; // mol in the following argument need not be this mol
      molecules[imol].set_user_defined_atom_colour_by_selections(indexed_residues_cids, colour_applies_to_non_carbon_atoms_also, mol);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}





coot::simple_mesh_t
molecules_container_t::get_map_contours_mesh(int imol, double position_x, double position_y, double position_z,
                                             float radius, float contour_level) {

   auto tp_0 = std::chrono::high_resolution_clock::now();
   coot::simple_mesh_t mesh;
   try {
      if (is_valid_map_molecule(imol)) {
         clipper::Coord_orth position(position_x, position_y, position_z);
         if (updating_maps_info.maps_need_an_update) {
            update_updating_maps(updating_maps_info.imol_model);
         }

         mesh = molecules[imol].get_map_contours_mesh(position, radius, contour_level, map_is_contoured_using_thread_pool_flag, &thread_pool);
      } else {
         std::cout << "WARNING:: get_map_contours_mesh() Not a valid map molecule " << imol << std::endl;
      }
   }
   catch (...) {
      std::cout << "An error occured in " << __FUNCTION__<< "() - this should not happen " << std::endl;
   }
   auto tp_1 = std::chrono::high_resolution_clock::now();
   contouring_time = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   if (mesh.vertices.size() > max_number_of_simple_mesh_vertices) {
      mesh.clear();
      mesh.status = 0;
      mesh.set_name("too-many-vertices-in-mesh");
   }
   return mesh;
}

//! get the mesh for the map contours using another map for colouring
//!
coot::simple_mesh_t
molecules_container_t::get_map_contours_mesh_using_other_map_for_colours(int imol_ref, int imol_map_for_colouring,
                                                                         double position_x, double position_y, double position_z,
                                                                         float radius, float contour_level,
                                                                         float other_map_for_colouring_min_value,
                                                                         float other_map_for_colouring_max_value,
                                                                         bool invert_colour_ramp) {
   coot::simple_mesh_t mesh;
   try {
      if (is_valid_map_molecule(imol_ref)) {
         if (is_valid_map_molecule(imol_map_for_colouring)) {
            clipper::Coord_orth position(position_x, position_y, position_z);
            molecules[imol_ref].set_other_map_for_colouring_min_max(other_map_for_colouring_min_value,
                                                                    other_map_for_colouring_max_value);
            molecules[imol_ref].set_other_map_for_colouring_invert_colour_ramp(invert_colour_ramp);
	    if (colour_map_by_other_map_user_defined_table.is_set()) {
	       mesh = molecules[imol_ref].get_map_contours_mesh_using_other_map_for_colours(position, radius, contour_level,
											    colour_map_by_other_map_user_defined_table,
											    molecules[imol_map_for_colouring].xmap);
	    } else {
	       mesh = molecules[imol_ref].get_map_contours_mesh_using_other_map_for_colours(position, radius, contour_level,
											    molecules[imol_map_for_colouring].xmap);
	    }
	 }
      }
   }
   catch (...) {
      std::cout << "An error occured in " << __FUNCTION__<< "() - this should not happen " << std::endl;
   }
   return mesh;
}


void
molecules_container_t::set_map_colour(int imol, float r, float g, float b) {

   if (is_valid_map_molecule(imol)) {
      coot::colour_holder ch(r,g,b);
      molecules[imol].set_map_colour(ch);
   }
}

void molecules_container_t::set_colour_map_for_map_coloured_by_other_map(std::vector<std::pair<double, std::vector<double> > > colour_table ) {

   bool debug = false;
   if (debug) {
      for (const auto &c : colour_table) {
	 double stop = c.first;
	 auto col = c.second;
	 std::cout << stop << " ";
	 for (auto p : col) {
	    std::cout << p << " ";
	 }
	 std::cout << std::endl;
      }
   }

   colour_map_by_other_map_user_defined_table.clear();

   for (const auto &c : colour_table) {
      float stop = c.first;
      auto col = c.second;
      if (col.size() > 2) {
	 glm::vec3 c(col[0], col[1], col[2]);
	 colour_map_by_other_map_user_defined_table.add_stop(stop, c);
      }
   }
}



// get the rotamer dodecs for the model
coot::simple_mesh_t
molecules_container_t::get_rotamer_dodecs(int imol) {
   coot::simple_mesh_t m;
   if (is_valid_model_molecule(imol)) {
      return molecules[imol].get_rotamer_dodecs(&geom, &rot_prob_tables);
   } else {
      std::cout << "WARNING:: in " << __FUNCTION__ << "() imol " << imol << " was not a valid model molecule " << std::endl;
   }
   return m;
}

//! get the rotamer dodecs for the model, not const because it regenerates the bonds.
//! @return an `instanced_mesh_t`
coot::instanced_mesh_t
molecules_container_t::get_rotamer_dodecs_instanced(int imol) {

   coot::instanced_mesh_t im;
   if (is_valid_model_molecule(imol)) {
      im = molecules[imol].get_rotamer_dodecs_instanced(&geom, &rot_prob_tables);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return im;
}


int
molecules_container_t::auto_fit_rotamer(int imol,
                                        const std::string &chain_id, int res_no, const std::string &ins_code,
                                        const std::string &alt_conf,
                                        int imol_map) {
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         const clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
         // std::cout << "debug:: mc::auto_fit_rotamer() calling the coot_molecule version with "
         //           << chain_id << " " << res_no << " " << alt_conf << std::endl;
         status = molecules[imol].auto_fit_rotamer(chain_id, res_no, ins_code, alt_conf, xmap, geom);
         set_updating_maps_need_an_update(imol);
      } else {
         std::cout << "ERROR:: mc::auto_fit_rotamer() not a valid map index " << imol_map << std::endl;
      }
   } else {
      std::cout << "ERROR:: mc::auto_fit_rotamer() not a valid model molecule " << imol << std::endl;
   }
   return status;
}

coot::molecule_t::rotamer_change_info_t
molecules_container_t::change_to_next_rotamer(int imol, const std::string &residue_cid, const std::string &alt_conf)  {

   coot::molecule_t::rotamer_change_info_t rci;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec = residue_cid_to_residue_spec(imol, residue_cid);
      rci = molecules[imol].change_to_next_rotamer(res_spec, alt_conf, geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return rci;
}

coot::molecule_t::rotamer_change_info_t
molecules_container_t::change_to_previous_rotamer(int imol, const std::string &residue_cid, const std::string &alt_conf)  {

   coot::molecule_t::rotamer_change_info_t rci;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec = residue_cid_to_residue_spec(imol, residue_cid);
      rci = molecules[imol].change_to_previous_rotamer(res_spec, alt_conf, geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return rci;
}


coot::molecule_t::rotamer_change_info_t
molecules_container_t::change_to_first_rotamer(int imol, const std::string &residue_cid, const std::string &alt_conf)  {

   coot::molecule_t::rotamer_change_info_t rci;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec = residue_cid_to_residue_spec(imol, residue_cid);
      rci = molecules[imol].change_to_first_rotamer(res_spec, alt_conf, geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return rci;
}


//! Change to the nth rotamer
//!
//! @param imol is the model molecule index
//! @param residue_cid is the atom selection CID e.g "//A/15" (all the atoms in residue 15 of chain A)
//! @param alt_conf is the alternate conformation, e.g. "A" or "B"
//!
//! @return the state of the change.
int molecules_container_t::set_residue_to_rotamer_number(int imol, const std::string &residue_cid,
                                                         const std::string &alt_conf, int rotamer_number) {

   int state = 0;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec = residue_cid_to_residue_spec(imol, residue_cid);
      state = molecules[imol].set_residue_to_rotamer_number(res_spec, alt_conf, rotamer_number, geom);
   }
   return state;
}


std::pair<int, unsigned int>
molecules_container_t::delete_atom(int imol,
                                   const std::string &chain_id, int res_no, const std::string &ins_code,
                                   const std::string &atom_name, const std::string &alt_conf) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec(chain_id, res_no, ins_code, atom_name, alt_conf);
      status = molecules[imol].delete_atom(atom_spec);
      set_updating_maps_need_an_update(imol);
   }
   unsigned int atom_count = get_number_of_atoms(imol);
   return std::make_pair(status, atom_count);
}

std::pair<int, unsigned int>
molecules_container_t::delete_atom_using_cid(int imol, const std::string &cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, cid);
      status = molecules[imol].delete_atom(atom_spec);
      set_updating_maps_need_an_update(imol);
   }
   unsigned int atom_count = get_number_of_atoms(imol);
   return std::make_pair(status, atom_count);
}



std::pair<int, unsigned int>
molecules_container_t::delete_residue(int imol,
                                      const std::string &chain_id, int res_no, const std::string &ins_code) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec(chain_id, res_no, ins_code);
      status = molecules[imol].delete_residue(residue_spec);
      set_updating_maps_need_an_update(imol);
   }
   unsigned int atom_count = get_number_of_atoms(imol);
   return std::make_pair(status, atom_count);
}


std::pair<int, unsigned int>
molecules_container_t::delete_residue_using_cid(int imol, const std::string &residue_cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec = residue_cid_to_residue_spec(imol, residue_cid);
      status = molecules[imol].delete_residue(residue_spec);
      set_updating_maps_need_an_update(imol);
   }
   unsigned int atom_count = get_number_of_atoms(imol);
   return std::make_pair(status, atom_count);
}

std::pair<int, unsigned int>
molecules_container_t::delete_residue_atoms_using_cid(int imol, const std::string &atom_cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, atom_cid);
      coot::residue_spec_t residue_spec(atom_spec);
      status = molecules[imol].delete_residue(residue_spec);
      set_updating_maps_need_an_update(imol);
   }
   unsigned int atom_count = get_number_of_atoms(imol);
   return std::make_pair(status, atom_count);
}

std::pair<int, unsigned int>
molecules_container_t::delete_residue_atoms_with_alt_conf(int imol, const std::string &chain_id,
                                                          int res_no, const std::string &ins_code,
                                                          const std::string &alt_conf) {
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      std::string atom_cid = std::string("//") + chain_id + std::string("/") + std::to_string(res_no) + ins_code;
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, atom_cid);
      coot::residue_spec_t residue_spec(atom_spec);
      status = molecules[imol].delete_residue_atoms_with_alt_conf(residue_spec, alt_conf);
      set_updating_maps_need_an_update(imol);
   }
   unsigned int atom_count = get_number_of_atoms(imol);
   return std::make_pair(status, atom_count);
}



std::pair<int, unsigned int>
molecules_container_t::delete_chain_using_cid(int imol, const std::string &cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].delete_chain_using_atom_cid(cid);
      set_updating_maps_need_an_update(imol);
   }
   unsigned int atom_count = get_number_of_atoms(imol);
   return std::make_pair(status, atom_count);
}


//! delete the atoms specified in the CID selection
//! @return 1 on successful deletion, return 0 on failure to delete.
std::pair<int, unsigned int>
molecules_container_t::delete_literal_using_cid(int imol, const std::string &cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].delete_literal_using_cid(cid);
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   unsigned int atom_count = get_number_of_atoms(imol);
   return std::make_pair(status, atom_count);
}



//where scope in ["ATOM","WATER", "RESIDUE","CHAIN","MOLECULE"]
std::pair<int, unsigned int>
molecules_container_t::delete_using_cid(int imol, const std::string &cid, const std::string &scope) {

   std::pair<int, unsigned int> r(0,0);
   if (scope == "ATOM") {
      r = delete_atom_using_cid(imol, cid);
      set_updating_maps_need_an_update(imol);
   }
   if (scope == "RESIDUE") {
      r = delete_residue_atoms_using_cid(imol, cid);
      set_updating_maps_need_an_update(imol);
   }
   if (scope == "CHAIN") {
      r = delete_chain_using_cid(imol, cid);
      set_updating_maps_need_an_update(imol);
   }
   if (scope == "LITERAL") {
      r = delete_literal_using_cid(imol, cid);
      set_updating_maps_need_an_update(imol);
   }
   if (scope == "MOLECULE") {
      int status = close_molecule(imol);
      if (status == 1) r.first = 1;
      set_updating_maps_need_an_update(imol);
   }
   return r;
}

// Old API
// int
// molecules_container_t::load_dictionary_file(const std::string &monomer_cif_file_name) {

//    int status = 0;

//    int read_number = 44;
//    geom.init_refmac_mon_lib(monomer_cif_file_name, read_number);
//    return status;
// }

std::vector<std::string>
molecules_container_t::non_standard_residue_types_in_model(int imol) const {
   std::vector<std::string> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].non_standard_residue_types_in_model();
   }
   return v;
}

//! @return the mean of the map or -1 is `imol_map` is not a map molecule index
float
molecules_container_t::get_map_mean(int imol) const {
   float m = -1.1;
   if (is_valid_map_molecule(imol)) {
      m = molecules[imol].get_map_mean();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid map molecule " << imol << std::endl;
   }
   return m;
}

float
molecules_container_t::get_map_rmsd_approx(int imol) const {
   float rmsd = -1.1;
   if (is_valid_map_molecule(imol)) {
      rmsd = molecules[imol].get_map_rmsd_approx();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid map molecule " << imol << std::endl;
   }
   return rmsd;
}


std::vector<coot::molecule_t::interesting_place_t>
molecules_container_t::difference_map_peaks(int imol_map, int imol_protein, float n_rmsd) const {

   std::vector<coot::molecule_t::interesting_place_t> v;
   if (is_valid_model_molecule(imol_protein)) {
      if (is_valid_map_molecule(imol_map)) {
         mmdb::Manager *m = get_mol(imol_protein);
         v = molecules[imol_map].difference_map_peaks(m, n_rmsd);
      } else {
         std::cout << "debug:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_map << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_protein << std::endl;
   }
   return v;
}



// return a useful message if the addition did not work
std::pair<int, std::string>
molecules_container_t::add_terminal_residue_directly(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) {

   std::string new_res_type = "ALA";
   int status = 0;
   std::string message;

   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t spec(chain_id, res_no, ins_code);
      mmdb::Residue *residue_p = molecules[imol].get_residue(spec);
      bool is_RNA = coot::util::is_nucleotide_by_dict(residue_p, geom);
      if (is_valid_map_molecule(imol_refinement_map) || is_RNA) {
         clipper::Xmap<float> &xmap = molecules[imol_refinement_map].xmap;
         coot::residue_spec_t residue_spec(chain_id, res_no, ins_code);
         std::pair<int, std::string> m = molecules[imol].add_terminal_residue_directly(residue_spec, new_res_type,
                                                                                       geom, xmap,
                                                                                       standard_residues_asc.mol,
                                                                                       thread_pool);
         status  = m.first;
         message = m.second;
         if (! message.empty())
            std::cout << "WARNING:: add_terminal_residue_directly(): " << message << std::endl;
         // write_coordinates(imol, "post-add-terminal-residue.pdb");
         set_updating_maps_need_an_update(imol);
      } else {
         std::cout << "debug:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_refinement_map << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return std::make_pair(status, message);
}

// 20221023-PE return an int for now so that I can write the binding
int
molecules_container_t::add_terminal_residue_directly_using_cid(int imol, const std::string &cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, cid);
      if (! atom_spec.empty()) {
         auto p = add_terminal_residue_directly(imol, atom_spec.chain_id, atom_spec.res_no, atom_spec.ins_code);
         status = p.first;
      }
   }
   return status;
}

//! buccaneer building, called by the above
int
molecules_container_t::add_terminal_residue_directly_using_bucca_ml_growing_using_cid(int imol, const std::string &cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, cid);
      coot::residue_spec_t res_spec(atom_spec);
      status = add_terminal_residue_directly_using_bucca_ml_growing(imol, res_spec);
   }
   return status;
}



// reset the rail_points (calls reset_the_rail_points()), updates the maps (using internal/clipper SFC)
// so, update your contour lines meshes after calling this function.
int
molecules_container_t::connect_updating_maps(int imol_model, int imol_with_data_info_attached, int imol_map_2fofc, int imol_map_fofc) {

   int status = 0;

   rail_point_history.clear();
   updating_maps_info.imol_model = imol_model;
   updating_maps_info.imol_2fofc = imol_map_2fofc;
   updating_maps_info.imol_fofc  = imol_map_fofc;
   updating_maps_info.imol_with_data_info_attached = imol_with_data_info_attached;
   imol_difference_map = imol_map_fofc;

   // Let's force a sfcalc_genmap here.
   updating_maps_info.maps_need_an_update = true;
   update_updating_maps(imol_model);

   return status;
}

void
molecules_container_t::associate_data_mtz_file_with_map(int imol, const std::string &data_mtz_file_name,
                                                        const std::string &f_col, const std::string &sigf_col,
                                                        const std::string &free_r_col) {

   if (is_valid_map_molecule(imol) || is_valid_model_molecule(imol)) {
      // 20221018-PE if free_r_col is not valid then Coot will (currently) crash on the structure factor calculation
      molecules[imol].associate_data_mtz_file_with_map(data_mtz_file_name, f_col, sigf_col, free_r_col);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid molecule " << imol << std::endl;
   }
}

/*! \brief Calculate structure factors from the model and update the given difference
           map accordingly */

// copied from:
// void
// graphics_info_t::sfcalc_genmap(int imol_model,
//                                int imol_map_with_data_attached,
//                                int imol_updating_difference_map) {
void
molecules_container_t::sfcalc_genmap(int imol_model,
                                     int imol_map_with_data_attached,
                                     int imol_updating_difference_map) {

   // I am keen for this function to be fast - so that it can be used with cryo-EM structures
   //
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map_with_data_attached)) {
         if (true) {
            if (is_valid_map_molecule(imol_updating_difference_map)) {
               if (molecules[imol_updating_difference_map].is_difference_map_p()) {
                  clipper::Xmap<float> *xmap_p = &molecules[imol_updating_difference_map].xmap;
                  try {
                     if (! on_going_updating_map_lock) {
                        on_going_updating_map_lock = true;
                        molecules[imol_map_with_data_attached].fill_fobs_sigfobs();
                        const clipper::HKL_data<clipper::data32::F_sigF> *fobs_data =
                           molecules[imol_map_with_data_attached].get_original_fobs_sigfobs();
                        const clipper::HKL_data<clipper::data32::Flag> *free_flag =
                           molecules[imol_map_with_data_attached].get_original_rfree_flags();
                        if (fobs_data && free_flag) {
                           molecules[imol_model].sfcalc_genmap(*fobs_data, *free_flag, xmap_p);
                        } else {
                           std::cout << "sfcalc_genmap() either fobs_data or free_flag were not set " << std::endl;
                        }
                        on_going_updating_map_lock = false;
                     } else {
                        std::cout << "DEBUG:: on_going_updating_map_lock was set! - aborting map update." << std::endl;
                     }
                  }
                  catch (const std::runtime_error &rte) {
                     std::cout << rte.what() << std::endl;
                  }
               } else {
                  std::cout << "sfcalc_genmap() not a valid difference map " << imol_updating_difference_map << std::endl;
               }
            } else {
               std::cout << "sfcalc_genmap() not a valid map (diff) " << imol_updating_difference_map << std::endl;
            }
         }
      } else {
         std::cout << "sfcalc_genmap() not a valid map " << imol_map_with_data_attached << std::endl;
      }
   } else {
      std::cout << "sfcalc_genmap() not a valid model " << imol_model << std::endl;
   }
}


#include "coot-utils/diff-diff-map-peaks.hh"

coot::util::sfcalc_genmap_stats_t
molecules_container_t::sfcalc_genmaps_using_bulk_solvent(int imol_model,
                                                         int imol_map_2fofc,  // this map should have the data attached.
                                                         int imol_map_fofc,
                                                         int imol_with_data_info_attached) {

   coot::util::sfcalc_genmap_stats_t stats;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map_2fofc)) {
         if (is_valid_map_molecule(imol_map_fofc)) {
            if (molecules[imol_map_fofc].is_difference_map_p()) {
               try {
                  if (! on_going_updating_map_lock) {
                     on_going_updating_map_lock = true;
                     molecules[imol_with_data_info_attached].fill_fobs_sigfobs();

                     // 20210815-PE used to be const reference (get_original_fobs_sigfobs() function changed too)
                     // const clipper::HKL_data<clipper::data32::F_sigF> &fobs_data = molecules[imol_map_with_data_attached].get_original_fobs_sigfobs();
                     // const clipper::HKL_data<clipper::data32::Flag> &free_flag = molecules[imol_map_with_data_attached].get_original_rfree_flags();
                     // now the full object (40us for RNAse test).
                     // 20210815-PE OK, the const reference was not the problem. But we will leave it as it is now, for now.
                     //
                     clipper::HKL_data<clipper::data32::F_sigF> *fobs_data_p = molecules[imol_with_data_info_attached].get_original_fobs_sigfobs();
                     clipper::HKL_data<clipper::data32::Flag>   *free_flag_p = molecules[imol_with_data_info_attached].get_original_rfree_flags();

                     if (fobs_data_p && free_flag_p) {

                        if (true) { // sanity check data

                           const clipper::HKL_info &hkls_check = fobs_data_p->base_hkl_info();
                           const clipper::Spacegroup &spgr_check = hkls_check.spacegroup();
                           const clipper::Cell &cell_check = fobs_data_p->base_cell();
                           const clipper::HKL_sampling &sampling_check = fobs_data_p->hkl_sampling();

                           if (false) {
                              std::cout << "DEBUG:: in sfcalc_genmaps_using_bulk_solvent() imol_map_with_data_attached "
                                        << imol_map_2fofc << std::endl;

                              std::cout << "DEBUG:: Sanity check in graphics_info_t:sfcalc_genmaps_using_bulk_solvent(): HKL_info: "
                                        << "base_cell: " << cell_check.format() << " "
                                        << "spacegroup: " << spgr_check.symbol_xhm() << " "
                                        << "sampling is null: " << sampling_check.is_null() << " "
                                        << "resolution: " << hkls_check.resolution().limit() << " "
                                        << "invsqreslim: " << hkls_check.resolution().invresolsq_limit() << " "
                                        << "num_reflections: " << hkls_check.num_reflections()
                                        << std::endl;
                           }
                        }

                        clipper::Xmap<float> &xmap_2fofc = molecules[imol_map_2fofc].xmap;
                        clipper::Xmap<float> &xmap_fofc  = molecules[imol_map_fofc].xmap;
                        molecules[imol_map_fofc].updating_maps_previous_difference_map = xmap_fofc;
                        stats = molecules[imol_model].sfcalc_genmaps_using_bulk_solvent(*fobs_data_p, *free_flag_p, &xmap_2fofc, &xmap_fofc);

                        { // diff difference map peaks
                           float rmsd = get_map_rmsd_approx(imol_map_fofc);
                           float base_level = 2.0 * rmsd;  // was 0.2 - this might need to be computed from the rmsd.
                           const clipper::Xmap<float> &m1 = molecules[imol_map_fofc].updating_maps_previous_difference_map;
                           const clipper::Xmap<float> &m2 = xmap_fofc;
                           std::vector<std::pair<clipper::Coord_orth, float> > v1 = coot::diff_diff_map_peaks(m1, m2, base_level);
                           // std::cout << "***************************** got " << v1.size() << " diff diff map peaks.... "
                           // << " using base level " << base_level << " with map rmsd " << rmsd << std::endl;
                           molecules[imol_map_fofc].set_updating_maps_diff_diff_map_peaks(v1);
                        }

                     } else {
                        std::cout << "ERROR:: null data pointer in graphics_info_t::sfcalc_genmaps_using_bulk_solvent() " << std::endl;
                     }
                     on_going_updating_map_lock = false;
                  }
               }
               catch (const std::runtime_error &rte) {
                  std::cout << rte.what() << std::endl;
               }
            }
         }
      }
   }
   latest_sfcalc_stats = stats;
   return stats;
}

//! shift_field B-factor refinement
//! @return success status
bool
molecules_container_t::shift_field_b_factor_refinement(int imol, int imol_with_data_attached) {

   bool status = false;
   int imol_map = imol_with_data_attached;
   try {
      if (is_valid_model_molecule(imol)) {
         if (is_valid_map_molecule(imol_map)) {
            molecules[imol_map].fill_fobs_sigfobs();
            const clipper::HKL_data<clipper::data32::F_sigF> *fobs_data = molecules[imol_map].get_original_fobs_sigfobs();
            const clipper::HKL_data<clipper::data32::Flag>  *rfree_flag = molecules[imol_map].get_original_rfree_flags();
            std::cout << "debug:: fobs_data" << fobs_data << " rfree " << rfree_flag << std::endl;
            if (fobs_data && rfree_flag) {
               status = molecules[imol].shiftfield_b_factor_refinement(*fobs_data, *rfree_flag);
               set_updating_maps_need_an_update(imol);
            }
         }
      }
   }
   catch(const std::runtime_error& rte) {
      std::cout << rte.what() << '\n';
   }
   return status;
}

//! @return a vector the position where the difference map has been flattened.
//! The associated float value is the ammount that the map has been flattened.
std::vector<std::pair<clipper::Coord_orth, float> >
molecules_container_t::get_diff_diff_map_peaks(int imol_map_fofc,
                                               float screen_centre_x, float screen_centre_y, float screen_centre_z) const {

   clipper::Coord_orth screen_centre(screen_centre_x, screen_centre_y, screen_centre_z); // also, is this used in this function?
   std::vector<std::pair<clipper::Coord_orth, float> > v;
   if (is_valid_map_molecule(imol_map_fofc)) {
      v = molecules[imol_map_fofc].get_updating_maps_diff_diff_map_peaks(screen_centre);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_map_fofc << std::endl;
   }
   return v;

}


int
molecules_container_t::rail_points_total() const { // the sum of all the rail ponts accumulated
   return rail_points_t::total(rail_point_history);
}

int
molecules_container_t::calculate_new_rail_points() {

   float rmsd = get_map_rmsd_approx(imol_difference_map);
   if (! rail_point_history.empty()) {
      const rail_points_t &prev = rail_point_history.back();
      rail_points_t new_points(rmsd, prev);
      rail_point_history.push_back(new_points);
      return new_points.map_rail_points_delta;
   } else {
      rail_points_t prev = rail_points_t(rmsd);
      rail_points_t new_points(rmsd, prev);
      rail_point_history.push_back(new_points);
      return new_points.map_rail_points_delta;
   }
}


// static
void
molecules_container_t::thread_for_refinement_loop_threaded() {

   // I think that there is a race condition here
   // check_and_warn_inverted_chirals_and_cis_peptides()
   // get called several times when the refine loop ends
   // (with success?).

   bool use_graphics_interface_flag = false;
   bool refinement_immediate_replacement_flag = true;

#if 0 // 20221018-PE this might not be the right thing

   if (restraints_lock) {
      if (false)
         std::cout << "debug:: thread_for_refinement_loop_threaded() restraints locked by "
                   << restraints_locking_function_name << std::endl;
      return;
   } else {

      if (use_graphics_interface_flag) {

         if (!refinement_immediate_replacement_flag) {

            // if there's not a refinement redraw function already running start up a new one.
            if (threaded_refinement_redraw_timeout_fn_id == -1) {
               GSourceFunc cb = GSourceFunc(regenerate_intermediate_atoms_bonds_timeout_function_and_draw);
               // int id = gtk_timeout_add(15, cb, NULL);

               int timeout_ms = 15;
               timeout_ms = 30; // 20220503-PE try this value
               int id = g_timeout_add(timeout_ms, cb, NULL);
               threaded_refinement_redraw_timeout_fn_id = id;
            }
         }
      }

      continue_threaded_refinement_loop = true;
      std::thread r(refinement_loop_threaded);
      r.detach();
   }
#endif

}



int
molecules_container_t::refine_direct(int imol, std::vector<mmdb::Residue *> rv, const std::string &alt_loc, int n_cycles) {

   // note to self: did you set imol_refinement_map?

   if (false)
      std::cout << "starting mc::refine_direct() with imol " << imol
                << " and imol_refinement_map " << imol_refinement_map
                << std::endl;

   // this is not stored in molecules_container!
   unsigned int max_number_of_threads = thread_pool.size();

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_refinement_map)) {
         const clipper::Xmap<float> &xmap = molecules[imol_refinement_map].xmap;
         status = molecules[imol].refine_direct(rv, alt_loc, xmap, max_number_of_threads,
                                                map_weight, n_cycles, geom,
                                       use_rama_plot_restraints, rama_plot_restraints_weight,
                                       use_torsion_restraints, torsion_restraints_weight,
                                       refinement_is_quiet);
         set_updating_maps_need_an_update(imol);
      } else {
	 logger.log(log_t::WARNING, logging::function_name_t(__FUNCTION__),
		    "not a valid map molecule, imol_refinement_map:", imol_refinement_map);
      }
   }
   return status;
}

int
molecules_container_t::refine_residues_using_atom_cid(int imol, const std::string &cid, const std::string &mode, int n_cycles) {

   auto debug_selected_residues = [cid] (const std::vector<mmdb::Residue *> &rv) {
      std::cout << "debug:: selection: refine_residues_using_atom_cid(): selected these " << rv.size() << " residues"
                << " from cid: \"" << cid << "\"" << std::endl;
      std::vector<mmdb::Residue *>::const_iterator it;
      for (it=rv.begin(); it!=rv.end(); ++it)
         std::cout << "   " << coot::residue_spec_t(*it) << std::endl;
   };

   if (false)
      std::cout << "starting refine_residues_using_atom_cid() with imol " << imol
                << " and imol_refinement_map " << imol_refinement_map
                << std::endl;

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_refinement_map)) {
         // coot::atom_spec_t spec = atom_cid_to_atom_spec(imol, cid);
         // status = refine_residues(imol, spec.chain_id, spec.res_no, spec.ins_code, spec.alt_conf, mode, n_cycles);
         std::vector<mmdb::Residue *> rv = molecules[imol].select_residues(cid, mode);

         // debug_selected_residues(rv);
         std::string alt_conf = "";
         status = refine_direct(imol, rv, alt_conf, n_cycles);
         set_updating_maps_need_an_update(imol);
      } else {
         std::cout << "WARNING:: " << __FUNCTION__ << " Not a valid map molecule " << imol_refinement_map << std::endl;
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << " Not a valid model molecule " << imol << std::endl;
   }
   return status;
}



int
molecules_container_t::refine_residues(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
                                       const std::string &alt_conf, const std::string &mode, int n_cycles) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec(chain_id, res_no, ins_code);
      std::vector<mmdb::Residue *> rv = molecules[imol].select_residues(residue_spec, mode);
      if (! rv.empty()) {
         status = refine_direct(imol, rv, alt_conf, n_cycles);
         set_updating_maps_need_an_update(imol);
      } else {
         std::cout << "WARNING:: in refine_residues() - empty residues." << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

int
molecules_container_t::refine_residue_range(int imol, const std::string &chain_id, int res_no_start, int res_no_end,
                                            int n_cycles) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      std::vector<mmdb::Residue *> rv = molecules[imol].select_residues(chain_id, res_no_start, res_no_end);
      if (! rv.empty()) {
         std::string alt_conf = "";
         status = refine_direct(imol, rv, alt_conf, n_cycles);
         set_updating_maps_need_an_update(imol);
      } else {
         std::cout << "WARNING:: in refine_residues() - empty residues." << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}



coot::refinement_results_t
molecules_container_t::refine_residues_vec(int imol,
                                           const std::vector<mmdb::Residue *> &residues,
                                           const std::string &alt_conf,
                                           mmdb::Manager *mol) {
   bool use_map_flag = true;
   if (false)
      std::cout << "INFO:: refine_residues_vec() with altconf \"" << alt_conf << "\"" << std::endl;

   coot::refinement_results_t rr = generate_molecule_and_refine(imol, residues, alt_conf, mol, use_map_flag);
   return rr;
}

// return -1 on failure to find a residue for insertion index
//
int
molecules_container_t::find_serial_number_for_insert(int seqnum_new,
                                                     const std::string &ins_code_for_new,
                                                     mmdb::Chain *chain_p) const {

   int iserial_no = -1;
   if (chain_p) {
      int current_diff = 999999;
      int nres = chain_p->GetNumberOfResidues();
      for (int ires=0; ires<nres; ires++) { // ires is a serial number
         mmdb::Residue *residue = chain_p->GetResidue(ires);

         // we are looking for the smallest negative diff:
         //
         int diff = residue->GetSeqNum() - seqnum_new;
         if ( (diff > 0) && (diff < current_diff) ) {
            iserial_no = ires;
            current_diff = diff;
         } else {
            if (diff == 0) {
               std::string ins_code_this = residue->GetInsCode();
               if (ins_code_this > ins_code_for_new) {
                  iserial_no = ires;
                  break;
               }
            }
         }
      }
   }
   return iserial_no;
}



std::pair<mmdb::Manager *, std::vector<mmdb::Residue *> >
molecules_container_t::create_mmdbmanager_from_res_vector(const std::vector<mmdb::Residue *> &residues,
                                                          int imol,
                                                          mmdb::Manager *mol_in,
                                                          std::string alt_conf) {

   // returned entities
   mmdb::Manager *new_mol = 0;
   std::vector<mmdb::Residue *> rv; // gets checked

   float dist_crit = 5.0;
   bool debug = false;

   if (debug) {
      std::cout << "############ starting create_mmdbmanager_from_res_vector() with these "
                << " residues " << std::endl;
      for (std::size_t ii=0; ii<residues.size(); ii++)
         std::cout << "   " << coot::residue_spec_t(residues[ii])  << std::endl;
      int udd_atom_index_handle = mol_in->GetUDDHandle(mmdb::UDR_ATOM, "atom index");
      std::cout << "############ udd for atom index from seeding molecule " << udd_atom_index_handle
                << std::endl;
      for (std::size_t ii=0; ii<residues.size(); ii++) {
         mmdb::Residue *residue_p = residues[ii];
         mmdb::Atom **residue_atoms = 0;
         int n_residue_atoms;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            int idx = -1;
            at->GetUDData(udd_atom_index_handle, idx);
            std::cout << "#### input residue atom " << coot::atom_spec_t(at) << " had udd index "
                      << idx << std::endl;
         }
      }
   }

   int n_flanker = 0; // a info/debugging counter

   if (residues.size() > 0) {

      // Also add the index of the reference residue (the one in molecules[imol].atom_selection.mol)
      // to the molecule that we are construction here. So that we can properly link
      // the residues in restraints_container (there we rather need to know the references indices,
      // not the indices from the fragment molecule)
      //

      std::pair<bool,std::string> use_alt_conf(false, "");
      if (! alt_conf.empty())
         use_alt_conf = std::pair<bool, std::string> (true, alt_conf);

      std::cout << "----------------- in create_mmdbmanager_from_res_vector() alt_conf is "
                << "\"" << alt_conf << "\"" << std::endl;
      std::cout << "----------------- in create_mmdbmanager_from_res_vector() use_alt_conf is "
                << use_alt_conf.first << "\"" << use_alt_conf.second << "\"" << std::endl;

      std::pair<bool, mmdb::Manager *> n_mol_1 =
         coot::util::create_mmdbmanager_from_residue_vector(residues, mol_in, use_alt_conf);

      // check that first is sane, so indent all this lot (when it works)

      if (n_mol_1.first) {

         int index_from_reference_residue_handle =
            n_mol_1.second->GetUDDHandle(mmdb::UDR_RESIDUE, "index from reference residue");

         if (false) { // debug
            int imod = 1;
            mmdb::Model *model_p = n_mol_1.second->GetModel(imod);
            if (model_p) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  int nres = chain_p->GetNumberOfResidues();
                  for (int ires=0; ires<nres; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     for (int iat=0; iat<n_atoms; iat++) {
                        mmdb::Atom *at = residue_p->GetAtom(iat);
                        int idx = -1;
                        at->GetUDData(index_from_reference_residue_handle, idx);
                        std::cout << "   create_mmdbmanager_from_residue_vector() returns this mol atom "
                                  << iat << " " << coot::atom_spec_t(at) << " with idx " << idx << std::endl;
                     }
                  }
               }
            }
         }

         new_mol = n_mol_1.second;
         mmdb::Model *model_p = new_mol->GetModel(1);

         // how many (free) residues were added to that model? (add them to rv)
         //
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               rv.push_back(residue_p);
            }
         }

         if (false) {
            for (std::size_t ir=0; ir<rv.size(); ir++) {
               mmdb::Residue *r = rv[ir];
               std::cout << "Moving Residue " << coot::residue_spec_t(r) << std::endl;
               mmdb::Atom **residue_atoms = 0;
               int n_residue_atoms;
               r->GetAtomTable(residue_atoms, n_residue_atoms);
               for (int iat=0; iat<n_residue_atoms; iat++) {
                  mmdb::Atom *at = residue_atoms[iat];
                  std::cout << "    " << coot::atom_spec_t(at) << std::endl;
               }
            }
         }

         short int whole_res_flag = 0;
         int atom_index_udd_handle = molecules[imol].atom_sel.UDDAtomIndexHandle;

         // Now the flanking residues:
         //
         std::vector<mmdb::Residue *> flankers_in_reference_mol;
         flankers_in_reference_mol.reserve(32); // say

         // find the residues that are close to the residues of
         // residues that are not part of residues
         //
         // We don't have quite the function that we need in coot-utils,
         // so we need to munge residues in to local_residues:
         std::vector<std::pair<bool, mmdb::Residue *> > local_residues;
         local_residues.resize(residues.size());
         for (std::size_t ires=0; ires<residues.size(); ires++)
            local_residues[ires] = std::pair<bool, mmdb::Residue *>(false, residues[ires]);
         std::map<mmdb::Residue *, std::set<mmdb::Residue *> > rnr =
            coot::residues_near_residues(local_residues, mol_in, dist_crit);
         // now fill @var{flankers_in_reference_mol} from rnr, avoiding residues
         // already in @var{residues}.
         std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it;
         for (it=rnr.begin(); it!=rnr.end(); ++it) {
            const std::set<mmdb::Residue *> &s = it->second;
            std::set<mmdb::Residue *>::const_iterator its;
            for (its=s.begin(); its!=s.end(); ++its) {
               mmdb::Residue *tres = *its;
               if (std::find(residues.begin(), residues.end(), tres) == residues.end())
                  if (std::find(flankers_in_reference_mol.begin(), flankers_in_reference_mol.end(), tres) == flankers_in_reference_mol.end())
                     flankers_in_reference_mol.push_back(tres);
            }
         }

         // So we have a vector of residues that were flankers in the
         // reference molecule, we need to add copies of those to
         // new_mol (making sure that they go into the correct chain).
         //
         if (false) { // debug
            std::cout << "debug:: ############ Found " << flankers_in_reference_mol.size()
                      << " flanking residues" << std::endl;

            for (unsigned int ires=0; ires<flankers_in_reference_mol.size(); ires++)
               std::cout << "     #### flankers_in_reference_mol: " << ires << " "
                         << coot::residue_spec_t(flankers_in_reference_mol[ires]) << std::endl;
         }


         for (unsigned int ires=0; ires<flankers_in_reference_mol.size(); ires++) {
            mmdb::Residue *r;

            std::string ref_res_chain_id = flankers_in_reference_mol[ires]->GetChainID();

            mmdb::Chain *chain_p = NULL;
            int n_new_mol_chains = model_p->GetNumberOfChains();
            for (int ich=0; ich<n_new_mol_chains; ich++) {
               if (ref_res_chain_id == model_p->GetChain(ich)->GetChainID()) {
                  chain_p = model_p->GetChain(ich);
                  break;
               }
            }

            if (! chain_p) {
               // Add a new one then.
               chain_p = new mmdb::Chain;
               chain_p->SetChainID(ref_res_chain_id.c_str());
               model_p->AddChain(chain_p);
            }

            if (false)
               std::cout << "debug:: flankers_in_reference_mol " << ires << " "
                         << coot::residue_spec_t(flankers_in_reference_mol[ires]) << " "
                         << "had index " << flankers_in_reference_mol[ires]->index
                         << std::endl;

            // get rid of this function at some stage
            bool embed_in_chain = false;
            r = coot::deep_copy_this_residue_old_style(flankers_in_reference_mol[ires],
                                                       alt_conf, whole_res_flag,
                                                       atom_index_udd_handle, embed_in_chain);

            if (r) {

               r->PutUDData(index_from_reference_residue_handle, flankers_in_reference_mol[ires]->index);

               // copy over the atom indices. UDDAtomIndexHandle in mol_n becomes UDDOldAtomIndexHandle
               // indices in the returned molecule

               int sni = find_serial_number_for_insert(r->GetSeqNum(), r->GetInsCode(), chain_p);

               if (false) { // debug
                  mmdb::Atom **residue_atoms = 0;
                  int n_residue_atoms;
                  std::cout << "Flanker Residue " << coot::residue_spec_t(r) << std::endl;
                  r->GetAtomTable(residue_atoms, n_residue_atoms);
                  for (int iat=0; iat<n_residue_atoms; iat++) {
                     mmdb::Atom *at = residue_atoms[iat];
                     std::cout << "    " << coot::atom_spec_t(at) << std::endl;
                  }
               }

               if (sni == -1)
                  chain_p->AddResidue(r); // at the end
               else
                  chain_p->InsResidue(r, sni);
               r->seqNum = flankers_in_reference_mol[ires]->GetSeqNum();
               r->SetResName(flankers_in_reference_mol[ires]->GetResName());
               n_flanker++;

               if (false)
                  std::cout << "debug:: create_mmdbmanager_from_residue_vector() inserted/added flanker "
                            << coot::residue_spec_t(r) << std::endl;

            }
         }

         // super-critical for correct peptide bonding in refinement!
         //
         coot::util::pdbcleanup_serial_residue_numbers(new_mol);

         if (debug) {
            int imod = 1;
            mmdb::Model *model_p = new_mol->GetModel(imod);
            if (model_p) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  int nres = chain_p->GetNumberOfResidues();
                  for (int ires=0; ires<nres; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     std::cout << "create_mmdb..  ^^^ " << coot::residue_spec_t(residue_p) << " "
                               << residue_p << " index " << residue_p->index
                               << std::endl;
                  }
               }
            }
         }

         if (debug)
            std::cout << "DEBUG:: in create_mmdbmanager_from_res_vector: " << rv.size()
                      << " free residues and " << n_flanker << " flankers" << std::endl;
      }
   }

   return std::pair <mmdb::Manager *, std::vector<mmdb::Residue *> > (new_mol, rv);
}



std::string
molecules_container_t::adjust_refinement_residue_name(const std::string &resname) const {

   std::string r = resname;
   if (resname == "UNK") r = "ALA"; // hack for KC/buccaneer.
   if (resname.length() > 2)
      if (resname[2] == ' ')
         r = resname.substr(0,2);
   return r;
}


// Return 0 (first) if any of the residues don't have a dictionary
// entry and a list of the residue type that don't have restraints.
//
std::pair<int, std::vector<std::string> >
molecules_container_t::check_dictionary_for_residue_restraints(int imol, mmdb::PResidue *SelResidues, int nSelResidues) {

   int status;
   bool status_OK = 1; // pass, by default
   std::vector<std::string> res_name_vec;

   for (int ires=0; ires<nSelResidues; ires++) {
      std::string resn(SelResidues[ires]->GetResName());
      std::string resname = adjust_refinement_residue_name(resn);
      status = geom.have_dictionary_for_residue_type(resname, imol, cif_dictionary_read_number);
      cif_dictionary_read_number++;
      if (! status) {
         status_OK = 0;
         res_name_vec.push_back(resname);
      }

      if (0)
         std::cout << "DEBUG:: have_dictionary_for_residues() on residue "
                   << ires << " of " << nSelResidues << ", "
                   << resname << " returns "
                   << status << std::endl;
      cif_dictionary_read_number++;
   }
   return std::pair<int, std::vector<std::string> > (status_OK, res_name_vec);
}

std::pair<int, std::vector<std::string> >
molecules_container_t::check_dictionary_for_residue_restraints(int imol, const std::vector<mmdb::Residue *> &residues) {

   std::vector<std::string> res_name_vec;
   std::pair<int, std::vector<std::string> > r(0, res_name_vec);
   for (unsigned int i=0; i<residues.size(); i++) {
      std::string resname = adjust_refinement_residue_name(residues[i]->GetResName());
      int status = geom.have_dictionary_for_residue_type(resname, imol, cif_dictionary_read_number);
      if (! status) {
         r.first = 0;
         r.second.push_back(resname);
      }
      cif_dictionary_read_number++; // not sure why this is needed.
   }
   return r;
}


std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > >
molecules_container_t::make_rotamer_torsions(const std::vector<std::pair<bool, mmdb::Residue *> > &local_residues) const {

   std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > > v;
   for (unsigned int i=0; i<local_residues.size(); i++) {
      if (! local_residues[i].first) {
         mmdb::Residue *residue_p = local_residues[i].second;
         std::string rn(residue_p->GetResName());
         if (coot::util::is_standard_amino_acid_name(rn)) {
            std::string alt_conf; // run through them all, ideally.
            coot::rotamer rot(residue_p, alt_conf, 1);
            coot::closest_rotamer_info_t cri = rot.get_closest_rotamer(rn);
            if (cri.residue_chi_angles.size() > 0) {
               std::vector<coot::dict_torsion_restraint_t> dictionary_vec;
               std::vector<std::vector<std::string> > rotamer_atom_names = rot.rotamer_atoms(rn);

               if (cri.residue_chi_angles.size() != rotamer_atom_names.size()) {

                  std::cout << "-------------- mismatch for " << coot::residue_spec_t(residue_p) << " "
                            << cri.residue_chi_angles.size() << " "  << rotamer_atom_names.size()
                            << " ---------------" << std::endl;
               } else {

                  for (unsigned int ichi=0; ichi<cri.residue_chi_angles.size(); ichi++) {
                     // we have to convert chi angles to atom names
                     double esd = 3.0; // 20210315-PE was 10.0. I want them tighter than that.
                     int per = 1;
                     std::string id = "chi " + coot::util::int_to_string(cri.residue_chi_angles[ichi].first);
                     const std::string &atom_name_1 = rotamer_atom_names[ichi][0];
                     const std::string &atom_name_2 = rotamer_atom_names[ichi][1];
                     const std::string &atom_name_3 = rotamer_atom_names[ichi][2];
                     const std::string &atom_name_4 = rotamer_atom_names[ichi][3];
                     double torsion = cri.residue_chi_angles[ichi].second;
                     coot::dict_torsion_restraint_t dr(id, atom_name_1, atom_name_2, atom_name_3, atom_name_4,
                                                       torsion, esd, per);
                     dictionary_vec.push_back(dr);
                  }

                  if (dictionary_vec.size() > 0) {
                     std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > p(residue_p, dictionary_vec);
                     v.push_back(p);
                  }
               }
            }
         }
      }
   }
   return v;
}



atom_selection_container_t
molecules_container_t::make_moving_atoms_asc(mmdb::Manager *residues_mol,
                                             const std::vector<mmdb::Residue *> &residues) const {

   // This also rebonds the imol_moving_atoms molecule

   atom_selection_container_t local_moving_atoms_asc;
   local_moving_atoms_asc.UDDAtomIndexHandle = -1;
   local_moving_atoms_asc.UDDOldAtomIndexHandle = residues_mol->GetUDDHandle(mmdb::UDR_ATOM, "old atom index");

   int SelHnd = residues_mol->NewSelection();

   for (unsigned int ir=0; ir<residues.size(); ir++) {
      const char *chain_id = residues[ir]->GetChainID();
      const char *inscode = residues[ir]->GetInsCode();
      int resno = residues[ir]->GetSeqNum();
      residues_mol->Select(SelHnd, mmdb::STYPE_ATOM,
                           0, chain_id,
                           resno, // starting resno, an int
                           inscode, // any insertion code
                           resno, // ending resno
                           inscode, // ending insertion code
                           "*", // any residue name
                           "*", // atom name
                           "*", // elements
                           "*",  // alt loc.
                           mmdb::SKEY_OR);
   }

   local_moving_atoms_asc.mol = residues_mol;
   local_moving_atoms_asc.SelectionHandle = SelHnd;
   residues_mol->GetSelIndex(local_moving_atoms_asc.SelectionHandle,
                             local_moving_atoms_asc.atom_selection,
                             local_moving_atoms_asc.n_selected_atoms);


   if (false) {
      std::cout << "returning an atom selection for all moving atoms "
                << local_moving_atoms_asc.n_selected_atoms << " atoms "
                << std::endl;
   }

   // This new block added so that we don't draw atoms in the "static" molecule when we have the
   // corresponding atoms in the moving atoms.
   //
#if 0 // 20221018-PE there is no drawing at the momment
   const atom_selection_container_t &imol_asc = molecules[imol_moving_atoms].atom_sel;
   std::set<int> atom_set = coot::atom_indices_in_other_molecule(imol_asc, local_moving_atoms_asc);

   if (false) { // debug atoms in other molecule
      std::set<int>::const_iterator it;
      for(it=atom_set.begin(); it!=atom_set.end(); it++) {
         int idx = *it;
         mmdb::Atom *at = imol_asc.atom_selection[idx];
         coot::atom_spec_t as(at);
         std::cout << " this is a moving atom: " << idx << " " << as << std::endl;
      }
   }

   if (false) { // debug old atom index
      for (int i=0; i<local_moving_atoms_asc.n_selected_atoms; i++) {
         mmdb::Atom *at = local_moving_atoms_asc.atom_selection[i];
         coot::atom_spec_t as(at);
         int idx = -1;
         at->GetUDData(local_moving_atoms_asc.UDDOldAtomIndexHandle, idx);
         std::cout << "DEBUG:: in make_moving_atoms_asc " << as << " idx " << idx << std::endl;
      }
   }
   // now rebond molecule imol without bonds to atoms in atom_set
   if (atom_set.size() > 0) {
      if (regenerate_bonds_needs_make_bonds_type_checked_flag) {
         molecules[imol_moving_atoms].make_bonds_type_checked(atom_set, __FUNCTION__);
      }
   }
#endif

   return local_moving_atoms_asc;
}

// static
void
molecules_container_t::all_atom_pulls_off() {
   for (std::size_t i=0; i<atom_pulls.size(); i++)
      atom_pulls[i].off();
   atom_pulls.clear();
}


// return the state of having found restraints.
bool
molecules_container_t::make_last_restraints(const std::vector<std::pair<bool,mmdb::Residue *> > &local_residues,
                                      const std::vector<mmdb::Link> &links,
                                      const coot::protein_geometry &geom,
                                      mmdb::Manager *mol_for_residue_selection,
                                      const std::vector<coot::atom_spec_t> &fixed_atom_specs,
                                      coot::restraint_usage_Flags flags,
                                      bool use_map_flag,
                                      const clipper::Xmap<float> *xmap_p) {

   bool do_torsion_restraints = true; // make this a data member
   double torsion_restraints_weight = 10.0;
   bool convert_dictionary_planes_to_improper_dihedrals_flag = false;
   double geometry_vs_map_weight = 25.5;
   bool do_trans_peptide_restraints = true;
   double rama_plot_restraints_weight = 20.0;
   bool do_rama_restraints = false;
   bool make_auto_h_bond_restraints_flag = false;
   coot::pseudo_restraint_bond_type pseudo_bonds_type = coot::NO_PSEUDO_BONDS;
   bool use_harmonic_approximation_for_NBCs = false;
   double pull_restraint_neighbour_displacement_max_radius = 1.0;
   double lennard_jones_epsilon = 1.0;
   int restraints_rama_type = 1;
   bool do_rotamer_restraints = false;
   double geman_mcclure_alpha = 0.1;
   bool do_numerical_gradients =  false;
   bool draw_gl_ramachandran_plot_flag = false;
   bool use_graphics_interface_flag = false;


   if (last_restraints) {
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << "    ERROR:: A: last_restraints not cleared up " << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
   }

   if (false) { // these are the passed residues, nothing more.
      std::cout << "debug:: on construction of restraints_container_t local_residues: "
                << std::endl;
      for (std::size_t jj=0; jj<local_residues.size(); jj++) {
         std::cout << "   " << coot::residue_spec_t(local_residues[jj].second)
                   << " is fixed: " << local_residues[jj].first << std::endl;
      }
   }

   // moving_atoms_extra_restraints_representation.clear();
   continue_threaded_refinement_loop = true; // no longer set in refinement_loop_threaded()

   // the refinment of torsion seems a bit confused? If it's in flags, why does it need an flag
   // of its own? I suspect that it doesn't. For now I will keep it (as it was).
   //
   bool do_residue_internal_torsions = false;
   if (do_torsion_restraints) {
      do_residue_internal_torsions = 1;
   }

   last_restraints = new
      coot::restraints_container_t(local_residues,
                                   links,
                                   geom,
                                   mol_for_residue_selection,
                                   fixed_atom_specs, xmap_p);

   std::cout << "debug:: on creation last_restraints is " << last_restraints << std::endl;

   last_restraints->set_torsion_restraints_weight(torsion_restraints_weight);

   if (convert_dictionary_planes_to_improper_dihedrals_flag) {
      last_restraints->set_convert_plane_restraints_to_improper_dihedral_restraints(true);
   }

   // This seems not to work yet.
   // last_restraints->set_dist_crit_for_bonded_pairs(9.0);

   if (use_map_flag)
      last_restraints->add_map(geometry_vs_map_weight);

   unsigned int n_threads = coot::get_max_number_of_threads();
   if (n_threads > 0)
      last_restraints->thread_pool(&thread_pool, n_threads);

   all_atom_pulls_off();
   particles_have_been_shown_already_for_this_round_flag = false;

   // elsewhere do this:
   // gtk_widget_remove_tick_callback(glareas[0], wait_for_hooray_refinement_tick_id);

   // moving_atoms_visited_residues.clear(); // this is used for HUD label colour

   int n_restraints = last_restraints->make_restraints(imol_moving_atoms,
                                                       geom, flags,
                                                       do_residue_internal_torsions,
                                                       do_trans_peptide_restraints,
                                                       rama_plot_restraints_weight,
                                                       do_rama_restraints,
                                                       true, true, make_auto_h_bond_restraints_flag,
                                                       pseudo_bonds_type);
                                                       // link and flank args default true

   if (use_harmonic_approximation_for_NBCs) {
      std::cout << "INFO:: using soft harmonic restraints for NBC" << std::endl;
      last_restraints->set_use_harmonic_approximations_for_nbcs(true);
   }

   if (pull_restraint_neighbour_displacement_max_radius > 1.99) {
      last_restraints->set_use_proportional_editing(true);
      last_restraints->pull_restraint_neighbour_displacement_max_radius =
         pull_restraint_neighbour_displacement_max_radius;
   }

   last_restraints->set_geman_mcclure_alpha(geman_mcclure_alpha);
   last_restraints->set_lennard_jones_epsilon(lennard_jones_epsilon);
   last_restraints->set_rama_type(restraints_rama_type);
   last_restraints->set_rama_plot_weight(rama_plot_restraints_weight); // >2? danger of non-convergence
                                                                       // if planar peptide restraints are used
   // Oh, I see... it's not just the non-Bonded contacts of the hydrogens.
   // It's the planes, chiral and angles too. Possibly bonds too.
   // How about marking non-H atoms in restraints that contain H atoms as
   // "invisible"? i.e. non-H atoms are not influenced by the positions of the
   // Hydrogen atoms (but Hydrogen atoms *are* influenced by the positions of the
   // non-Hydrogen atoms). This seems like a lot of work. Might be easier to turn
   // off angle restraints for H-X-X (but not H-X-H) in the first instance, that
   // should go most of the way to what "invisible" atoms would do, I imagine.
   // is_H_non_bonded_contact should be renamed to is_H_turn_offable_restraint
   // or something.
   //
   // last_restraints->set_apply_H_non_bonded_contacts(false);

   if (do_rotamer_restraints) {
      std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > > rotamer_torsions = make_rotamer_torsions(local_residues);
      std::cout << "debug:: calling add_or_replace_torsion_restraints_with_closest_rotamer_restraints() from make_last_restraints() " << std::endl;
      last_restraints->add_or_replace_torsion_restraints_with_closest_rotamer_restraints(rotamer_torsions);
   }

   if (molecules[imol_moving_atoms].extra_restraints.has_restraints()) {
      std::cout << "debug:: calling add_extra_restraints() from make_last_restraints() " << std::endl;
      last_restraints->add_extra_restraints(imol_moving_atoms, "user-defined restraints called from make_last_restraints()",
                                            molecules[imol_moving_atoms].extra_restraints, geom);
   }

   if (do_numerical_gradients)
      last_restraints->set_do_numerical_gradients();

   bool found_restraints_flag = false;

   if (last_restraints->size() > 0) {

      last_restraints->analyze_for_bad_restraints();
      thread_for_refinement_loop_threaded();
      found_restraints_flag = true;
      // rr.found_restraints_flag = true;
      draw_gl_ramachandran_plot_flag = true;

      // are you looking for conditionally_wait_for_refinement_to_finish() ?

      if (refinement_immediate_replacement_flag) {
         // wait until refinement finishes
         while (restraints_lock) {
            std::this_thread::sleep_for(std::chrono::milliseconds(7));
            std::cout << "INFO:: make_last_restraints() [immediate] restraints locked by "
                      << restraints_locking_function_name << std::endl;
         }
      }

   } else {
      continue_threaded_refinement_loop = false;
      if (use_graphics_interface_flag) {
         // GtkWidget *widget = create_no_restraints_info_dialog();
         // GtkWidget *widget = widget_from_builder("no_restraints_info_dialog");
         // gtk_widget_show(widget);
      }
   }

   return found_restraints_flag;
}


// simple mmdb::Residue * interface to refinement.  20081216
coot::refinement_results_t
molecules_container_t::generate_molecule_and_refine(int imol,  // needed for UDD Atom handle transfer
                                                    const std::vector<mmdb::Residue *> &residues_in,
                                                    const std::string &alt_conf,
                                                    mmdb::Manager *mol,
                                                    bool use_map_flag) {

   // 20221018-PE make a function in the class
   auto set_refinement_flags = [] () {
      return coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
   };
   int cif_dictionary_read_number = 44; // make this a class member

   bool do_torsion_restraints = true;
   bool do_rama_restraints = false; // or true?
   bool moving_atoms_have_hydrogens_displayed = false;


   coot::refinement_results_t rr(0, GSL_CONTINUE, "");

   if (is_valid_map_molecule(imol_refinement_map) || (! use_map_flag)) {
      // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
      coot::restraint_usage_Flags flags = set_refinement_flags();
      bool do_residue_internal_torsions = false;
      if (do_torsion_restraints) {
         do_residue_internal_torsions = 1;
         flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
      }

      if (do_rama_restraints)
         // flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
         flags = coot::ALL_RESTRAINTS;

      std::vector<coot::atom_spec_t> fixed_atom_specs = molecules[imol].get_fixed_atoms();

      // refinement goes a bit wonky if there are multiple occurrances of the same residue
      // in input residue vector, so let's filter out duplicates here
      //
      std::vector<mmdb::Residue *> residues;
      std::set<mmdb::Residue *> residues_set;
      std::set<mmdb::Residue *>::const_iterator it;
      for (std::size_t i=0; i<residues_in.size(); i++)
         residues_set.insert(residues_in[i]);
      residues.reserve(residues_set.size());
      for(it=residues_set.begin(); it!=residues_set.end(); ++it)
         residues.push_back(*it);

      // OK, so the passed residues are the residues in the graphics_info_t::molecules[imol]
      // molecule.  We need to do 2 things:
      //
      // convert the mmdb::Residue *s of the passed residues to the mmdb::Residue *s of residues mol
      //
      // and
      //
      // in create_mmdbmanager_from_res_vector() make sure that that contains the flanking atoms.
      // The create_mmdbmanager_from_res_vector() from this class is used, not coot::util
      //
      // The flanking atoms are fixed the passed residues are not fixed.
      // Keep a clear head.

      std::vector<std::string> residue_types = coot::util::residue_types_in_residue_vec(residues);
      // use try_dynamic_add()
      bool have_restraints = geom.have_restraints_dictionary_for_residue_types(residue_types, imol, cif_dictionary_read_number);
      cif_dictionary_read_number += residue_types.size();

      if (have_restraints) {

         std::string residues_alt_conf = alt_conf;
         imol_moving_atoms = imol;
         std::pair<mmdb::Manager *, std::vector<mmdb::Residue *> > residues_mol_and_res_vec =
            create_mmdbmanager_from_res_vector(residues, imol, mol, residues_alt_conf);

         if (true) { // debug
            mmdb::Manager *residues_mol = residues_mol_and_res_vec.first;
            int imod = 1;
            mmdb::Model *model_p = residues_mol->GetModel(imod);
            if (model_p) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  std::cout << "DEBUG:: in generate_molecule_and_refine() residues_mol_and_res_vec mol: chain: "
                            << chain_p->GetChainID() << std::endl;
                  int nres = chain_p->GetNumberOfResidues();
                  for (int ires=0; ires<nres; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     std::cout << "DEBUG:: in generate_molecule_and_refine() residues_mol_and_res_vec mol:   residue "
                               << coot::residue_spec_t(residue_p) << " residue "
                               << residue_p << " chain " << residue_p->chain << " index "
                               << residue_p->index << std::endl;
                  }
               }
            }
         }

         // We only want to act on these new residues and molecule, if
         // there is something there.
         //
         if (residues_mol_and_res_vec.first) {

            // Now we want to do an atom name check.  This stops exploding residues.
            //
            bool check_hydrogens_too_flag = false;
            std::pair<bool, std::vector<std::pair<mmdb::Residue *, std::vector<std::string> > > >
               icheck_atoms = geom.atoms_match_dictionary(imol, residues, check_hydrogens_too_flag, false);

            if (! icheck_atoms.first) {

               std::cout << "WARNING:: non-matching atoms! " << std::endl;

            } else {

               moving_atoms_have_hydrogens_displayed = true;
               if (! molecules[imol].hydrogen_atom_should_be_drawn())
                  moving_atoms_have_hydrogens_displayed = false;

               atom_selection_container_t local_moving_atoms_asc =
                  make_moving_atoms_asc(residues_mol_and_res_vec.first, residues);

               // 20221018-PE make_moving_atoms_graphics_object(imol, local_moving_atoms_asc); not today!

               int n_cis = coot::util::count_cis_peptides(local_moving_atoms_asc.mol);
               // moving_atoms_n_cis_peptides = n_cis; // 20221018-PE not today

               std::vector<std::pair<bool,mmdb::Residue *> > local_residues;  // not fixed.
               for (unsigned int i=0; i<residues_mol_and_res_vec.second.size(); i++)
                  local_residues.push_back(std::pair<bool, mmdb::Residue *>(0, residues_mol_and_res_vec.second[i]));

               moving_atoms_asc_type = NEW_COORDS_REPLACE;

               int imol_for_map = imol_refinement_map;
               clipper::Xmap<float> *xmap_p = dummy_xmap;

               if (is_valid_map_molecule(imol_for_map))
                  xmap_p = &molecules[imol_for_map].xmap;

               bool found_restraints_flag = make_last_restraints(local_residues,
                                                                 local_moving_atoms_asc.links,
                                                                 geom,
                                                                 residues_mol_and_res_vec.first,
                                                                 fixed_atom_specs,
                                                                 flags, use_map_flag, xmap_p);

               if (last_restraints) {
                  // 20220423-PE I can't do this here because setup_minimize() has not been called yet
                  // rr = last_restraints->get_refinement_results();
               }
               rr.found_restraints_flag = found_restraints_flag;

            }
         }
      } else {

         // we didn't have restraints for everything.
         //
         // If we are in this state, we need to make that apparent to the calling function
         rr.found_restraints_flag = false;
         rr.info_text = "Missing or incomplete dictionaries";

         std::pair<int, std::vector<std::string> > icheck =
            check_dictionary_for_residue_restraints(imol, residues);
         if (icheck.first == 0) {
            std::cout << "WARNING:: <some info here about missing residue types> " << std::endl;
         }
      }
   }
   return rr;
}


int
molecules_container_t::mutate(int imol, const std::string &cid, const std::string &new_residue_type) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, cid);
      coot::residue_spec_t residue_spec(atom_spec);
      status = molecules[imol].mutate(residue_spec, new_residue_type);
      set_updating_maps_need_an_update(imol);
      // qstd::cout << "mutate status " << status << std::endl;
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

#include "coot-utils/blob-line.hh"

std::pair<bool, clipper::Coord_orth>
molecules_container_t::go_to_blob(float x1, float y1, float z1, float x2, float y2, float z2, float contour_level) {

   std::pair<bool, clipper::Coord_orth> p;

   clipper::Coord_orth p1(x1,y1,z1);
   clipper::Coord_orth p2(x2,y2,z2);

   // iterate through all the maps (another day)

   if (is_valid_map_molecule(imol_refinement_map)) {
      const clipper::Xmap<float> &xmap = molecules[imol_refinement_map].xmap;
      std::pair<bool, clipper::Coord_orth> pp = coot::find_peak_along_line_favour_front(p1, p2, contour_level, xmap);
      p = pp;
   }
   return p;
}


int
molecules_container_t::side_chain_180(int imol, const std::string &atom_cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, atom_cid);
      coot::residue_spec_t residue_spec(atom_spec);
      status = molecules[imol].side_chain_180(residue_spec, atom_spec.alt_conf, &geom);
      set_updating_maps_need_an_update(imol); // won't change much
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;

}

std::string
molecules_container_t::jed_flip(int imol, const std::string &atom_cid, bool invert_selection) {

   std::string message;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, atom_cid);
      coot::residue_spec_t res_spec(atom_spec);
      std::string atom_name = atom_spec.atom_name;
      std::string alt_conf  = atom_spec.alt_conf;
      message = molecules[imol].jed_flip(res_spec, atom_name, alt_conf, invert_selection, &geom);
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return message;
}


#include "ligand/ligand.hh"

int
molecules_container_t::add_waters(int imol_model, int imol_map) {

   int n_waters_added = -1;
   int ligand_water_n_cycles = 3;

   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
         coot::ligand lig;
         int n_cycles = ligand_water_n_cycles; // 3 by default

         try {
            // n_cycles = 1; // for debugging.
            short int mask_waters_flag; // treat waters like other atoms?
            // mask_waters_flag = g.find_ligand_mask_waters_flag;
            mask_waters_flag = 1; // when looking for waters we should not
            // ignore the waters that already exist.
            // short int do_flood_flag = 0;    // don't flood fill the map with waters for now.

            lig.import_map_from(molecules[imol_map].xmap, molecules[imol_map].get_map_rmsd_approx());
            // lig.set_masked_map_value(-2.0); // sigma level of masked map gets distorted
            lig.set_map_atom_mask_radius(1.9); // Angstroms
            lig.set_water_to_protein_distance_limits(ligand_water_to_protein_distance_lim_max,
                                                     ligand_water_to_protein_distance_lim_min);
            lig.set_variance_limit(ligand_water_variance_limit);
            lig.mask_map(molecules[imol_model].atom_sel.mol, mask_waters_flag);
            // lig.output_map("masked-for-waters.map");
            // std::cout << "debug:: add_waters(): using n-sigma cut off " << ligand_water_sigma_cut_off << std::endl;
            logger.log(log_t::DEBUG, logging::function_name_t("add_waters"),
                       "using n-sigma cut-ff ", ligand_water_sigma_cut_off);

            lig.water_fit(ligand_water_sigma_cut_off, n_cycles);

            coot::minimol::molecule water_mol = lig.water_mol();
            molecules[imol_model].insert_waters_into_molecule(water_mol, "HOH");
            n_waters_added = water_mol.count_atoms();
            set_updating_maps_need_an_update(imol_model);
         }
         catch (const std::bad_alloc &e) {
            std::cout << "WARNING::" << e.what() << std::endl;
         }
      }
   }
   return n_waters_added;
}

//! Flood with dummy atoms
//! @param imol is the model molecule index
//! @param imol_map is the map molecule index
//!
//! @return the number of waters added on a success, -1 on failure.
int
molecules_container_t::flood(int imol_model, int imol_map, float n_rmsd) {

   int n_waters_added = -1;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
         mmdb::Manager *mol = get_mol(imol_model);
         if (mol) {
            float flood_atom_mask_radius = 1.0; // was 1.4
            coot::ligand lig;
            short int mask_waters_flag = true;
            lig.import_map_from(molecules[imol_map].xmap);
            lig.mask_map(mol, mask_waters_flag);
            lig.set_cluster_size_check_off();
            lig.set_chemically_sensible_check_off();
            lig.set_sphericity_test_off();

            lig.set_map_atom_mask_radius(flood_atom_mask_radius);
            // lig.set_water_to_protein_distance_limits(10.0, 1.5); // should not be
            // used in lig.
            float water_to_protein_max_dist = 99.9;
            float water_to_protein_min_dist = 1.5;
            lig.set_water_to_protein_distance_limits(water_to_protein_max_dist,
                                                     water_to_protein_min_dist);
            float map_rmsd = get_map_rmsd_approx(imol_map);
            lig.flood2(n_rmsd * map_rmsd);
            coot::minimol::molecule water_mol = lig.water_mol();
            molecules[imol_model].insert_waters_into_molecule(water_mol, "DUM");
            n_waters_added = water_mol.get_number_of_atoms();
         }
      }
   }
   return n_waters_added;
}



std::vector<coot::molecule_t::interesting_place_t>
molecules_container_t::unmodelled_blobs(int imol_model, int imol_map, float rmsd_cut_off) const {

   std::vector<coot::molecule_t::interesting_place_t> v;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {

         coot::ligand lig;

         short int mask_waters_flag = true;
         float sigma = molecules[imol_map].get_map_rmsd_approx();
         lig.import_map_from(molecules[imol_map].xmap, sigma);
         lig.set_map_atom_mask_radius(1.9); // Angstrom
         lig.mask_map(molecules[imol_model].atom_sel.mol, mask_waters_flag);
         float sigma_cut_off = rmsd_cut_off;
         std::cout << "Unmodelled blobs using sigma cut off " << sigma_cut_off << std::endl;
         int n_cycles = 1;
         lig.water_fit(sigma_cut_off, n_cycles);
         std::vector<std::pair<clipper::Coord_orth, double> > big_blobs = lig.big_blobs();
         int n_big_blobs = lig.big_blobs().size();
         if (n_big_blobs > 0) {
            for (unsigned int i=0; i<big_blobs.size(); i++) {
               std::string l = std::string("Blob ") + std::to_string(i+1);
               clipper::Coord_orth pt = big_blobs[i].first;
               coot::molecule_t::interesting_place_t ip("Unmodelled Blob", pt, l);
               ip.set_feature_value(big_blobs[i].second);
               v.push_back(ip);
            }
         }
      }
   }
   return v;
}





std::pair<int, unsigned int>
molecules_container_t::delete_side_chain(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) {

   int status = 0;
   // 20221025-PE Fill me later
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec(chain_id, res_no, ins_code);
      status = molecules[imol].delete_side_chain(res_spec);
      if (status) {
         set_updating_maps_need_an_update(imol);
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   unsigned int atom_count = get_number_of_atoms(imol);
   return std::make_pair(status, atom_count);
}

//! delete side chain
//! @return 1 on successful deletion, return 0 on failure to delete.
std::pair<int, unsigned int>
molecules_container_t::delete_side_chain_using_cid(int imol, const std::string &cid) {

   int status = 0;

   // 20221025-PE Fill me later
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec = residue_cid_to_residue_spec(imol, cid);
      if (! res_spec.unset_p()) {
         status = molecules[imol].delete_side_chain(res_spec);
         set_updating_maps_need_an_update(imol);
      } else {
         std::cout << "WARNING:: in delete_side_chain_using_cid didn't find residue from cid " << cid << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   unsigned int atom_count = get_number_of_atoms(imol);
   return std::make_pair(status, atom_count);
}



int
molecules_container_t::fill_partial_residue(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) {

   int status = 0;
   std::string alt_conf;

   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec(chain_id, res_no, ins_code);
      if (is_valid_map_molecule(imol_refinement_map)) {
         const clipper::Xmap<float> &xmap = molecules.at(imol_refinement_map).xmap;
         molecules[imol].fill_partial_residue(res_spec, alt_conf, xmap, geom);
         set_updating_maps_need_an_update(imol);
      } else {
         std::cout << "WARNING:: fill_partial_residue() incorrect imol_refinement_map " << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

//! fill the specified residue
//! @return 1 on a successful fill, 0 on failure.
int
molecules_container_t::fill_partial_residue_using_cid(int imol, const std::string &cid) {

   int status = 0;
   std::string alt_conf;

   if (is_valid_model_molecule(imol)) {
      std::pair<bool, coot::residue_spec_t> res_spec_pair = molecules[imol].cid_to_residue_spec(cid);
      if (res_spec_pair.first) {
         const auto &res_spec = res_spec_pair.second;
         if (is_valid_map_molecule(imol_refinement_map)) {
            const clipper::Xmap<float> &xmap = molecules.at(imol_refinement_map).xmap;
            molecules[imol].fill_partial_residue(res_spec, alt_conf, xmap, geom);
            set_updating_maps_need_an_update(imol);
         } else {
            std::cout << "WARNING:: fill_partial_residue_using_cid() incorrect imol_refinement_map " << std::endl;
         }
      } else {
         std::cout << "fill_partial_residue_using_cid() residue not found " << cid << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;

}



int
molecules_container_t::fill_partial_residues(int imol) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_refinement_map)) {
         const clipper::Xmap<float> &xmap = molecules.at(imol_refinement_map).xmap;
         status = molecules[imol].fill_partial_residues(xmap, &geom);
         set_updating_maps_need_an_update(imol);
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;

}

//! Add N-linked glycosylation
//!
//! @param imol_model is the model molecule index
//! @param imol_map is the map molecule index
//! @param glycosylation_name is the type of glycosylation, one of:
//!       "NAG-NAG-BMA" or "high-mannose" or "hybrid" or "mammalian-biantennary" or "plant-biantennary"
//! @param asn_chain_id is the chain-id of the ASN to which the carbohydrate is to be added
//! @param asn_res_no is the residue number of the ASN to which the carbohydrate is to be added
void
molecules_container_t::add_named_glyco_tree(int imol_model, int imol_map, const std::string &glycosylation_name,
                                            const std::string &asn_chain_id, int asn_res_no) {

   int status = 0;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
         const clipper::Xmap<float> &xmap = molecules.at(imol_map).xmap;
         molecules[imol_model].add_named_glyco_tree(glycosylation_name, asn_chain_id, asn_res_no, xmap, &geom);
      } else {
         std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_map << std::endl;
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_model << std::endl;
   }

}

std::vector<std::string>
molecules_container_t::get_chains_in_model(int imol) const {

   std::vector<std::string> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].chains_in_model();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}

std::vector<std::pair<coot::residue_spec_t, std::string> >
molecules_container_t::get_single_letter_codes_for_chain(int imol, const std::string &chain_id) const {

   std::vector<std::pair<coot::residue_spec_t, std::string> > v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_single_letter_codes_for_chain(chain_id);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}


std::vector<std::string>
molecules_container_t::get_residue_names_with_no_dictionary(int imol) const {

   std::vector<std::string> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_residue_names_with_no_dictionary(geom);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}



int
molecules_container_t::apply_transformation_to_atom_selection(int imol, const std::string &atoms_selection_cid,
                                                              int n_atoms, // for validation of the atom selection
                                                              float m00, float m01, float m02,
                                                              float m10, float m11, float m12,
                                                              float m20, float m21, float m22,
                                                              float c0, float c1, float c2, // the centre of the rotation
                                                              float t0, float t1, float t2) { // translation

   int n_atoms_moved = 0;
   if (is_valid_model_molecule(imol)) {
      clipper::Coord_orth rotation_centre(c0, c1, c2);
      clipper::Coord_orth t(t0, t1, t2);
      clipper::Mat33<double> m(m00, m01, m02, m10, m11, m12, m20, m21, m22);
      clipper::RTop_orth rtop_orth(m, t);
      n_atoms_moved = molecules[imol].apply_transformation_to_atom_selection(atoms_selection_cid, n_atoms, rotation_centre, rtop_orth);
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return n_atoms_moved;

}


int
molecules_container_t::new_positions_for_residue_atoms(int imol, const std::string &residue_cid, std::vector<coot::api::moved_atom_t> &moved_atoms) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].new_positions_for_residue_atoms(residue_cid, moved_atoms);
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

int
molecules_container_t::new_positions_for_atoms_in_residues(int imol, const std::vector<coot::api::moved_residue_t> &moved_residues) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].new_positions_for_atoms_in_residues(moved_residues);
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;

}

// put this in a new file molecules_container_validation.cc

#include "coot-utils/pepflip-using-difference-map.hh"


std::vector<coot::molecule_t::interesting_place_t>
molecules_container_t::pepflips_using_difference_map(int imol_coords, int imol_difference_map, float n_sigma) const {

   auto mmdb_to_clipper = [] (mmdb::Atom *at) {
      return clipper::Coord_orth(at->x, at->y, at->z);
   };

   std::vector<coot::molecule_t::interesting_place_t> v;

   if (is_valid_model_molecule(imol_coords)) {
      if (is_valid_map_molecule(imol_difference_map)) {
         if (molecules[imol_difference_map].is_difference_map_p()) {
            const clipper::Xmap<float> &diff_xmap = molecules[imol_difference_map].xmap;
            mmdb::Manager *mol = get_mol(imol_coords);
            coot::pepflip_using_difference_map pf(mol, diff_xmap);
            std::vector<coot::residue_spec_t> flips = pf.get_suggested_flips(n_sigma);
            for (std::size_t i=0; i<flips.size(); i++) {
               const auto &res_spec = flips[i];
               mmdb::Residue *residue_this_p = get_residue(imol_coords, res_spec);
               if (residue_this_p) {
                  coot::residue_spec_t res_spec_next =  res_spec.next();
                  mmdb::Residue *residue_next_p = get_residue(imol_coords, res_spec);
                  if (residue_next_p) {
                     std::string feature_type = "Difference Map Suggest Pepflip";
                     std::string label = "Flip: " + res_spec.format();
                     mmdb::Atom *at_1 = residue_this_p->GetAtom(" CA ");
                     mmdb::Atom *at_2 = residue_next_p->GetAtom(" CA ");
                     if (at_1 && at_2) {
                        clipper::Coord_orth pt_1 = mmdb_to_clipper(at_1);
                        clipper::Coord_orth pt_2 = mmdb_to_clipper(at_2);
                        clipper::Coord_orth pos = 0.5 * (pt_1 + pt_2);
                        float f = static_cast<float>(i)/static_cast<float>(flips.size());
                        float badness = 20.0 + 50.0 * (1.0 - f);
                        coot::molecule_t::interesting_place_t ip(feature_type, res_spec, pos, label);
                        ip.set_badness_value(badness);
                        v.push_back(ip);
                     }
                  }
               }
            }
         }
      }
   }
   std::cout << "DEBUG:: pepflips_using_difference_map() returns " << v.size() << " flips" << std::endl;
   return v;

}

//! @return a vector of residue specifiers for the ligand residues - the residue name is encoded
//! in the `string_user_data` data item of the residue specifier
std::vector<coot::residue_spec_t>
molecules_container_t::get_non_standard_residues_in_molecule(int imol) const {

   std::vector<coot::residue_spec_t> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_non_standard_residues_in_molecule();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}


//! Try to read the dictionaries for any residue type in imol that as yet does not have
//! a dictionary
//!
//! @param imol is the model molecule index
//! @return true if there were no dictionary for new types that couldn't be read.
bool
molecules_container_t::try_read_dictionaries_for_new_residue_types(int imol) {

   bool status = true;
   std::vector<std::string> v = get_residue_names_with_no_dictionary(imol);
   if (v.empty()) {
      return true;
   } else {
      int read_number = 50;
      int imol_enc_any = get_imol_enc_any();
      for (const auto &rn : v) {
         int ss = geom.check_and_try_dynamic_add(rn, imol_enc_any, read_number);
         read_number++;
      }
   }
   return status;
}

//! \brief set the residue properties
//!
//! a list of propperty maps such as `{"chain-id": "A", "res-no": 34, "ins-code": "", "worm-radius": 1.2}`
//!
//! @param json_string is the properties in JSON format
//! @return true
bool molecules_container_t::set_residue_properties(int imol, const std::string &json_string) {

   bool status = false;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].set_residue_properties(json_string);
   }
   return status;
}

// \brief clear the reisidue properties
//!
//! @param imol is the model molecule index
void molecules_container_t::clear_residue_properties(int imol) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].clear_residue_properties();
   }
}


coot::simple_mesh_t
molecules_container_t::get_molecular_representation_mesh(int imol, const std::string &cid, const std::string &colour_scheme,
                                                         const std::string &style, int secondary_structure_usage_flag) {

   coot::simple_mesh_t mesh;
   if (is_valid_model_molecule(imol)) {

#if 0 // testing the colour rules
      add_colour_rule(imol, "//N", "tomato");
      add_colour_rule(imol, "//A", "pink");
      add_colour_rule(imol, "//B", "skyblue");
      add_colour_rule(imol, "//G", "cyan");
      add_colour_rule(imol, "//P", "brown");
      add_colour_rule(imol, "//R", "yellow");
      add_colour_rule(imol, "//A/40-46",   "green");
      add_colour_rule(imol, "//A/207-214", "green");
      add_colour_rule(imol, "//A/217-223", "green");
      add_colour_rule(imol, "//A/243-249", "green");
      add_colour_rule(imol, "//A/276-283", "green");
#endif

      mesh = molecules[imol].get_molecular_representation_mesh(cid, colour_scheme, style, secondary_structure_usage_flag);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return mesh;
}


//! get interesting places (does not work yet)
//! @return a vector of `validation_information_t`
std::vector<coot::molecule_t::interesting_place_t>
molecules_container_t::get_interesting_places(int imol, const std::string &mode) const {

   std::vector<coot::molecule_t::interesting_place_t> v;
   std::cout << "Nothing here yet" << std::endl;

   return v;
}


//! get Gaussian surface representation
coot::simple_mesh_t
molecules_container_t::get_gaussian_surface(int imol, float sigma, float contour_level,
                                            float box_radius, float grid_scale, float fft_b_factor) const {

   coot::simple_mesh_t mesh;
   if (is_valid_model_molecule(imol)) {
      mesh = molecules[imol].get_gaussian_surface(sigma, contour_level, box_radius, grid_scale, fft_b_factor);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return mesh;

}

coot::simple_mesh_t
molecules_container_t::get_gaussian_surface_for_atom_selection(int imol, const std::string &cid,
                                                               float sigma, float contour_level,
                                                               float box_radius, float grid_scale, float b_factor) const {
   coot::simple_mesh_t mesh;
   if (is_valid_model_molecule(imol)) {
      mesh = molecules[imol].get_gaussian_surface_for_atom_selection(cid, sigma, contour_level, box_radius,
                                                                     grid_scale, b_factor);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return mesh;
}

#include "coot-utils/coot-map-heavy.hh"

int molecules_container_t::gaussian_surface_to_map_molecule(int imol_ref, int imol_model, const std::string &cid,
                                                            float sigma, float box_radius, float fft_b_factor) {

   int imol_new = -1;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_ref)) {
         const clipper::Xmap<float> &xmap_ref = molecules[imol_ref].xmap;
         mmdb::Manager *mol = molecules[imol_model].get_mol();
         if (mol) {
            int sel_hnd = mol->NewSelection(); // d
            mol->Select(sel_hnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_NEW);
            clipper::Xmap<float> xmap = coot::util::make_gaussian_atom_map_for_mask(xmap_ref, mol, sel_hnd, sigma, box_radius);
            imol_new = molecules.size();
            std::string name = "Gaussian Map";
            if (fft_b_factor != 0.0f) {
               clipper::Xmap<float> xmap_blur = coot::util::sharpen_blur_map(xmap, fft_b_factor);
               xmap = xmap_blur;
            }
            bool is_em_map = true; // not sure
            coot::molecule_t m = coot::molecule_t(name, imol_new, xmap, is_em_map);
            molecules.push_back(m);
            mol->DeleteSelection(sel_hnd);
         }
      }
   }
   return imol_new;
}

#include "density-contour/gaussian-surface.hh"

// Make a map from a gaussian surface
//!
//! Waters are not included in the surface calculation
//!
//! @param imol is the model molecule index
//! @param cid is the atom selection CID
//! @param sigma default 4.4
//! @param contour_level default 4.0
//! @param box_radius default 5.0
//! @param grid_scale default 0.7
//! @param b_factor default 100.0 (use 0.0 for no FFT-B-factor smoothing)
//!
//! @return a new molecule index for the map or -1 on failur
int molecules_container_t::gaussian_surface_to_map_molecule_v2(int imol, const std::string &cid,
                                                              float sigma, float box_radius,
                                                              float grid_scale, float fft_b_factor) {

   int imol_new = -1;
   if (is_valid_model_molecule(imol)) {
      int chain_cid_mode = 0; // use cid to make the atom selection
      mmdb::Manager *mol = molecules[imol].get_mol();
      if (mol) {
         coot::gaussian_surface_t gauss_surf(mol, cid, chain_cid_mode,
                                             sigma, 0.5,
                                             box_radius, grid_scale, fft_b_factor);
         clipper::Xmap<float> xmap = gauss_surf.get_xmap();
         imol_new = molecules.size();
         std::string name = "Gaussian Map";
         bool is_em_map = true; // not sure
         coot::molecule_t m = coot::molecule_t(name, imol, xmap, is_em_map);
         molecules.push_back(m);
      }
   }
   return imol_new;
}


//! get chemical feaatures for the given residue
coot::simple_mesh_t
molecules_container_t::get_chemical_features_mesh(int imol, const std::string &cid) const {

   coot::simple_mesh_t mesh;
   if (is_valid_model_molecule(imol)) {
      mesh = molecules[imol].get_chemical_features_mesh(cid, geom);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return mesh;
}


//! add an alternative conformation for the specified residue
int
molecules_container_t::add_alternative_conformation(int imol_model, const std::string &cid) {

   int status = 0;
   if (is_valid_model_molecule(imol_model)) {
      status = molecules[imol_model].add_alternative_conformation(cid);
      set_updating_maps_need_an_update(imol_model);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_model << std::endl;
   }
   return status;
}


//! return@ an object that has information about residues without dictionaries and residues with missing atom
//! in the the specified molecule
coot::util::missing_atom_info
molecules_container_t::missing_atoms_info_raw(int imol) {

   coot::util::missing_atom_info mai;

   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = molecules[imol].atom_sel.mol;
      bool do_missing_hydrogen_atoms_flag = false;
      mai = coot::util::missing_atoms(mol, do_missing_hydrogen_atoms_flag, &geom);
   }
   return mai;
}


//! @return an object that has information about residues without dictionaries and residues with missing atom
//! in the the specified molecule
std::vector<coot::residue_spec_t>
molecules_container_t::residues_with_missing_atoms(int imol) {

   std::vector<coot::residue_spec_t> v;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = molecules[imol].atom_sel.mol;
      bool do_missing_hydrogen_atoms_flag = false;
      coot::util::missing_atom_info mai = coot::util::missing_atoms(mol, do_missing_hydrogen_atoms_flag, &geom);
      for (unsigned int i=0; i<mai.residues_with_missing_atoms.size(); i++) {
         mmdb::Residue *r = mai.residues_with_missing_atoms[i];
         v.push_back(coot::residue_spec_t(r));
      }
   }
   return v;
}

//! @return the instanced mesh for the specified ligand
coot::instanced_mesh_t
molecules_container_t::contact_dots_for_ligand(int imol, const std::string &cid,
                                               unsigned int num_subdivisions) const {

   coot::instanced_mesh_t im;
   if (is_valid_model_molecule(imol)) {
      im = molecules[imol].contact_dots_for_ligand(cid, geom, num_subdivisions);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return im;
}


//! @return the instanced mesh for the specified molecule
coot::instanced_mesh_t
molecules_container_t::all_molecule_contact_dots(int imol, unsigned int num_subdivisions) const {

   coot::instanced_mesh_t im;
   if (is_valid_model_molecule(imol)) {
      im = molecules[imol].all_molecule_contact_dots(geom, num_subdivisions);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return im;
}

//! If any colour rule has been set for this molecule, then we will use those. Otherwise, colorChainsScheme() will be called
//! (and that his its internal colour-by-chain colouring scheme).
//!
void
molecules_container_t::add_colour_rule(int imol, const std::string &selection_cid, const std::string &colour) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].add_colour_rule(selection_cid, colour);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

}

//! add multiple colour rules, combined like the following "//A/1^#cc0000|//A/2^#cb0002|//A/3^#c00007"
//!
void
molecules_container_t::add_colour_rules_multi(int imol, const std::string &selections_and_colours_combo_string) {

   if (is_valid_model_molecule(imol)) {
      std::vector<std::string> sel_col_pairs = coot::util::split_string(selections_and_colours_combo_string, "|");
      for (const auto &pair_string : sel_col_pairs) {
         std::vector<std::string> parts = coot::util::split_string(pair_string, "^");
         if (parts.size() == 2) {
            const std::string &selection_cid = parts[0];
            const std::string &colour        = parts[1];
            molecules[imol].add_colour_rule(selection_cid, colour);
         }
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

}


//! delete the colour rules for the given molecule
void
molecules_container_t::delete_colour_rules(int imol) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].delete_colour_rules();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}


//! print the colour rules
void
molecules_container_t::print_colour_rules(int imol) const {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].print_colour_rules();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}


//! get the colour rules
std::vector<std::pair<std::string, std::string> >
molecules_container_t::get_colour_rules(int imol) const {

   std::vector<std::pair<std::string, std::string> > v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_colour_rules();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}


//! Update float parameter for MoleculesToTriangles molecular mesh
void
molecules_container_t::M2T_updateFloatParameter(int imol, const std::string &param_name, float value) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].M2T_updateFloatParameter(param_name, value);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}

//! Update int parameter for MoleculesToTriangles molecular mesh
void
molecules_container_t::M2T_updateIntParameter(int imol, const std::string &param_name, int value) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].M2T_updateIntParameter(param_name, value);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

}


//! add waters, updating imol_model (of course)
//! @return 1 on a successful move, 0 on failure.
int
molecules_container_t::add_hydrogen_atoms(int imol_model) {
   int status = 0;
   if (is_valid_model_molecule(imol_model)) {
      status = molecules[imol_model].add_hydrogen_atoms(&geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_model << std::endl;
   }
   return status;
}


//! delete hydrogen atoms, updating imol_model (of course)
//! @return 1 on a successful move, 0 on failure.
int
molecules_container_t::delete_hydrogen_atoms(int imol_model) {
   int status = 0;
   if (is_valid_model_molecule(imol_model)) {
      status = molecules[imol_model].delete_hydrogen_atoms();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_model << std::endl;
   }
   return status;
}

//! generate GM self restraints
int
molecules_container_t::generate_self_restraints(int imol, float local_dist_max) {

   int status = -1;
   if (is_valid_model_molecule(imol)) {
      molecules[imol].generate_self_restraints(local_dist_max, geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status; // nothing useful.
}


//! generate GM self restraints for the given chain
void
molecules_container_t::generate_chain_self_restraints(int imol,
                                                      float local_dist_max,
                                                      const std::string &chain_id) {
   if (is_valid_model_molecule(imol)) {
      molecules[imol].generate_chain_self_restraints(local_dist_max, chain_id, geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}

//! generate GM self restraints for the given residues.
//! `residue_cids" is a "||"-separated list of residues, e.g. "//A/12||//A/14||/B/56"
void
molecules_container_t::generate_local_self_restraints(int imol, float local_dist_max,
                                                      const std::string &multi_selection_cid) {

   std::string residue_cids = multi_selection_cid; // 20231220-PE old style, residue by residue
   bool do_old_style = false;
   if (is_valid_model_molecule(imol)) {
      if (do_old_style) {
         std::vector<coot::residue_spec_t> residue_specs;
         std::vector<std::string> parts = coot::util::split_string(residue_cids, "||");
         for (const auto &part : parts) {
            coot::residue_spec_t rs = residue_cid_to_residue_spec(imol, part);
            if (! rs.empty())
               residue_specs.push_back(rs);
         }
         molecules[imol].generate_local_self_restraints(local_dist_max, residue_specs, geom);
      } else {
         molecules[imol].generate_local_self_restraints(local_dist_max, multi_selection_cid, geom);
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}



//! generate parallel plane restraints (for RNA and DNA)
void
molecules_container_t::add_parallel_plane_restraint(int imol,
                                                    const std::string &residue_cid_1,
                                                    const std::string &residue_cid_2) {

   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t rs_1 = residue_cid_to_residue_spec(imol, residue_cid_1);
      coot::residue_spec_t rs_2 = residue_cid_to_residue_spec(imol, residue_cid_1);
      molecules[imol].add_parallel_plane_restraint(rs_1, rs_2);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}

//! clear the extra restraints

void
molecules_container_t::clear_extra_restraints(int imol) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].clear_extra_restraints();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

}



// ----------------------- map utils
int
molecules_container_t::sharpen_blur_map(int imol_map, float b_factor, bool in_place_flag) {

   int imol_new = -1;
   if (is_valid_map_molecule(imol_map)) {
      const clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
      clipper::Xmap<float> xmap_new = coot::util::sharpen_blur_map(xmap, b_factor);
      if (in_place_flag) {
         molecules[imol_map].xmap = xmap_new;
      } else {
         std::string name = molecules[imol_map].get_name();
         if (b_factor < 0.0)
            name += " Sharpen ";
         else
            name += " Blur ";
         name += std::to_string(b_factor);
         imol_new = molecules.size();
         bool is_em = molecules[imol_map].is_EM_map();
         coot::molecule_t cm(name, imol_new, is_em);
         cm.xmap = xmap_new;
         molecules.push_back(cm);
      }
   }
   return imol_new;
}

//! create a new map that is blurred/sharpened
//! @return the molecule index of the new map or -1 on failure or if `in_place_flag` was true.
int
molecules_container_t::sharpen_blur_map_with_resample(int imol_map, float b_factor, float resample_factor, bool in_place_flag) {

   int imol_new = -1;
   if (is_valid_map_molecule(imol_map)) {
      const clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
      clipper::Xmap<float> xmap_new = coot::util::sharpen_blur_map_with_resample(xmap, b_factor, resample_factor);
      if (in_place_flag) {
         molecules[imol_map].xmap = xmap_new;
      } else {
         std::string name = molecules[imol_map].get_name();
         if (b_factor < 0.0)
            name += " Sharpen ";
         else
            name += " Blur ";
         name += std::to_string(b_factor);
         if (resample_factor < 0.999 || resample_factor > 1.001) {
            name += " Resample ";
            name += coot::util::float_to_string_using_dec_pl(resample_factor, 2);
         }
         imol_new = molecules.size();
         coot::molecule_t cm(name, imol_new);
         cm.xmap = xmap_new;
         molecules.push_back(cm);
      }
   }
   return imol_new;
}




//! Make a vector of maps that are split by chain-id of the input imol
//! @return a vector of the map molecule indices.
std::vector<int>
molecules_container_t::make_masked_maps_split_by_chain(int imol, int imol_map) {

   std::vector<int> v;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         coot::ligand lig;
         mmdb::Manager *mol = molecules[imol].atom_sel.mol;
         lig.set_map_atom_mask_radius(3.3);
         lig.import_map_from(molecules[imol_map].xmap);
         // monster
         std::vector<std::pair<std::string, clipper::Xmap<float> > > maps = lig.make_masked_maps_split_by_chain(mol);
         std::cout << "INFO:: made " << maps.size() << " masked maps" << std::endl;
         std::string orig_map_name = molecules[imol_map].get_name();
         bool is_em_flag = molecules[imol_map].is_EM_map();
         for(unsigned int i=0; i<maps.size(); i++) {
            std::string map_name = std::string("Map for chain ") + maps[i].first;
            map_name += std::string(" of ") + orig_map_name;
            int idx = molecules.size();
            coot::molecule_t cm(map_name, idx, maps[i].second, is_em_flag);
            molecules.push_back(cm);
            v.push_back(idx);
         }
      } else {
         std::cout << "WARNING:: molecule " << imol_map << " is not a valid map molecule"
                   << std::endl;
      }
   } else {
      std::cout << "WARNING:: molecule " << imol_map << " is not a valid model molecule"
                << std::endl;
   }
   return v;
}


//! mask map by atom selection (note the argument order is reversed compared to the coot api).
//!
//! the ``invert_flag`` changes the parts of the map that are masked, so to highlight the density
//! for a ligand one would pass the ``cid`` for the ligand and invert_flag as true, so that the
//! parts of the map that are not the ligand are suppressed.
//!
//! @return the index of the new map - or -1 on failure
int
molecules_container_t::mask_map_by_atom_selection(int imol_coords, int imol_map, const std::string &multi_cids,
                                                  float radius, bool invert_flag) {

   int imol_map_new = -1;
   if (is_valid_model_molecule(imol_coords)) {
      if (is_valid_map_molecule(imol_map)) {
         coot::ligand lig;
         lig.import_map_from(molecules[imol_map].xmap);

         float map_mask_atom_radius = 1.5; // check
         lig.set_map_atom_mask_radius(map_mask_atom_radius);

         int selectionhandle = molecules[imol_coords].atom_sel.mol->NewSelection();
         // molecules[imol_coords].atom_sel.mol->Select(selectionhandle, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_NEW);

         std::vector<std::string> parts = coot::util::split_string(multi_cids, "||");
         for (const auto &part : parts) {
            std::cout << "-------------------------- selecting part: " << part << std::endl;
            molecules[imol_coords].atom_sel.mol->Select(selectionhandle, mmdb::STYPE_ATOM, part.c_str(), mmdb::SKEY_OR);
         }

         if (radius > 0.0) lig.set_map_atom_mask_radius(radius);
         lig.mask_map(molecules[imol_coords].atom_sel.mol, selectionhandle, invert_flag);
         imol_map_new = molecules.size();
         std::string name = get_molecule_name(imol_map);
         std::string new_name = name + " Masked Map";
         bool is_em_map_flag = molecules[imol_map].is_EM_map();
         coot::molecule_t cm(new_name, imol_map_new, lig.masked_map(), is_em_map_flag);
         molecules.push_back(cm);
         molecules[imol_coords].atom_sel.mol->DeleteSelection(selectionhandle);
      } else {
         std::cout << "WARNING:: molecule " << imol_map << " is not a valid map molecule"
                   << std::endl;
      }
   } else {
      std::cout << "WARNING:: molecule " << imol_map << " is not a valid model molecule"
                << std::endl;
   }
   return imol_map_new;
}

//! @return the "EM" status of this molecule. Return false on not-a-map.
bool
molecules_container_t::is_EM_map(int imol) const {

   bool status = false;
   if (is_valid_map_molecule(imol)) {
      status = molecules[imol].is_EM_map();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}


// --------------------- symmetry

//! symmetry

// std::vector<std::pair<symm_trans_t, Cell_Translation> >
coot::symmetry_info_t
molecules_container_t::get_symmetry(int imol, float symmetry_search_radius, float x, float y, float z) const {

   coot::symmetry_info_t si;
   if (is_valid_model_molecule(imol)) {
      coot::Cartesian symmetry_centre(x, y, z);
      si = molecules[imol].get_symmetry(symmetry_search_radius, symmetry_centre);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return si;
}


//! set the colour wheel rotation base for the specified molecule
void
molecules_container_t::set_colour_wheel_rotation_base(int imol, float r) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].set_colour_wheel_rotation_base(r);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}


//! set the base colour - to be used as a base for colour wheel rotation
void
molecules_container_t::set_base_colour_for_bonds(int imol, float r, float g, float b) {
   if (is_valid_model_molecule(imol)) {
      molecules[imol].set_base_colour_for_bonds(r,g,b);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}

//! @return the string of the contents of the given file-name.
std::string
molecules_container_t::file_name_to_string(const std::string &file_name) const {

   std::string s;
   std::ifstream f(file_name.c_str(), std::ios::binary);
   if (!f) {
      std::cout << "WARNING:: Failed to open " << file_name << std::endl;
   } else {
      std::ostringstream ostrm;
      ostrm << f.rdbuf();
      s = ostrm.str();
   }
   return s;
}

//! the stored data set file name
std::string
molecules_container_t::get_data_set_file_name(int imol) const {

   std::string r;
   if (is_valid_model_molecule(imol))
      r = molecules[imol].Refmac_mtz_filename();

   return r;
}


coot::simple::molecule_t
molecules_container_t::get_simple_molecule(int imol, const std::string &residue_cid, bool draw_hydrogen_atoms_flag) {

   coot::simple::molecule_t sm;
   if (is_valid_model_molecule(imol)) {
      sm = molecules[imol].get_simple_molecule(imol, residue_cid, draw_hydrogen_atoms_flag, &geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return sm;
}



//! @return a vector of lines for non-bonded contacts and hydrogen bonds
generic_3d_lines_bonds_box_t
molecules_container_t::make_exportable_environment_bond_box(int imol, coot::residue_spec_t &spec, float max_dist) {

   // this function is non-const because the Bonds_lines function needs a mutable protein_geometry

   generic_3d_lines_bonds_box_t bonds_box;
   if (is_valid_model_molecule(imol)) {
      bonds_box = molecules[imol].make_exportable_environment_bond_box(spec, max_dist, geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return bonds_box;

}


//! use bespoke carbon atom colour
void
molecules_container_t::set_use_bespoke_carbon_atom_colour(int imol, bool state) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].set_use_bespoke_carbon_atom_colour(state);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}

//! set bespoke carbon atom colour
void
molecules_container_t::set_bespoke_carbon_atom_colour(int imol, const coot::colour_t &col) {
   if (is_valid_model_molecule(imol)) {
      molecules[imol].set_bespoke_carbon_atom_colour(col);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}


void
molecules_container_t::add_target_position_restraint(int imol, const std::string &atom_cid, float pos_x, float pos_y, float pos_z) {
   if (is_valid_model_molecule(imol)) {
      molecules[imol].add_target_position_restraint(atom_cid, pos_x, pos_y, pos_z);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}

void
molecules_container_t::init_refinement_of_molecule_as_fragment_based_on_reference(int imol_frag, int imol_ref, int imol_map) {

   // make last_restraints
   if (is_valid_model_molecule(imol_frag)) {
      if (is_valid_model_molecule(imol_ref)) {
         if (is_valid_map_molecule(imol_map)) {
            mmdb::Manager *mol_ref = molecules[imol_ref].atom_sel.mol;
            // this is a fragment molecule - a few residues. mol_ref is used for the NBC an peptide links
            // a the end of the fragment
            const clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
            std::cout << "debug:: in init_refinement_of_molecule_as_fragment_based_on_reference() "
                      << " cell " << xmap.cell().descr().format() << std::endl;
            molecules[imol_frag].init_all_molecule_refinement(imol_ref, geom, xmap, map_weight, &thread_pool);
         } else {
            std::cout << "WARNING:: in init_refinement_of_molecule_as_fragment_based_on_reference()"
                      << " not a valid map" << std::endl;
         }
      } else {
         std::cout << "WARNING:: in init_refinement_of_molecule_as_fragment_based_on_reference()"
                   << " not a valid ref model" << std::endl;
      }
   } else {
         std::cout << "WARNING:: in init_refinement_of_molecule_as_fragment_based_on_reference()"
                   << " not a valid frag model" << std::endl;
   }
}


coot::instanced_mesh_t
molecules_container_t::add_target_position_restraint_and_refine(int imol, const std::string &atom_cid,
                                                                float pos_x, float pos_y, float pos_z,
                                                                int n_cycles) {

   coot::instanced_mesh_t m;
   if (is_valid_model_molecule(imol)) {
      m = molecules[imol].add_target_position_restraint_and_refine(atom_cid, pos_x, pos_y, pos_z, n_cycles, &geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return m;
}


//! clear any and all drag-atom target position restraints
void
molecules_container_t::clear_target_position_restraints(int imol) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].clear_target_position_restraints();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}

//! clear target_position restraint
void
molecules_container_t::clear_target_position_restraint(int imol, const std::string &atom_cid) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].clear_target_position_restraint(atom_cid);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}

//! clear target_position restraint if it is (or they are) close to their target position
void
molecules_container_t::turn_off_when_close_target_position_restraint(int imol) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].turn_off_when_close_target_position_restraint();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}



void
molecules_container_t::clear_refinement(int imol) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].clear_refinement();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}


//! fix atoms during refinement
void
molecules_container_t::fix_atom_selection_during_refinement(int imol, const std::string &atom_selection_cid) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].fix_atom_selection_during_refinement(atom_selection_cid);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

}

//! Run some cycles of refinement and return a mesh
//! That way we can see the molecule animate as it refines
std::pair<int, coot::instanced_mesh_t>
molecules_container_t::refine(int imol, int n_cycles) {

   coot::instanced_mesh_t im;
   int status = 0;
   if (is_valid_model_molecule(imol)) {

      std::cout << "debug:: in mc::refine() calling refine_using_last_restraints() using imol " << imol << std::endl;
      status = molecules[imol].refine_using_last_restraints(n_cycles);
      std::string mode = "COLOUR-BY-CHAIN-AND-DICTIONARY";
      bool draw_hydrogen_atoms_flag = true; // use data member as we do for draw_missing_residue_loops_flag?
      bool show_atoms_as_aniso_flag = false;
      bool show_aniso_atoms_as_ortep = false;
      bool show_aniso_atoms_as_empty = false;
      float aniso_probability = 0.5f;
      unsigned int smoothness_factor = 1;
      im = molecules[imol].get_bonds_mesh_instanced(mode, &geom, true, 0.12, 1.4,
                                                    show_atoms_as_aniso_flag,
                                                    aniso_probability,
                                                    show_aniso_atoms_as_ortep,
                                                    show_aniso_atoms_as_empty,
                                                    smoothness_factor,
                                                    draw_hydrogen_atoms_flag, draw_missing_residue_loops_flag);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return std::make_pair(status, im);
}

//! get the mesh for extra restraints (currently an empty object is returned)
coot::instanced_mesh_t
molecules_container_t::get_extra_restraints_mesh(int imol, int mode) {

   coot::instanced_mesh_t m;
   if (is_valid_model_molecule(imol)) {
      m = molecules[imol].get_extra_restraints_mesh(mode);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

   return m;
}

//! flip the hand of the map.
//! @return the molecule index of the new map, or -1 on failure.
int
molecules_container_t::flip_hand(int imol_map) {

   int imol_new = -1;
   if (is_valid_map_molecule(imol_map)) {
      clipper::Xmap<float> xmap = molecules[imol_map].xmap;
      coot::util::flip_hand(&xmap);
      std::string name = "Flipped Hand of " + molecules[imol_map].get_name();
      imol_new = molecules.size();
      molecules.push_back(coot::molecule_t(name, imol_new, xmap, true));
   }
   return imol_new;
}


//! @return the suggested initial contour level. Return -1 on not-a-map
float
molecules_container_t::get_suggested_initial_contour_level(int imol) const {

   float l = -1;
   if (is_valid_map_molecule(imol)) {
      l = molecules[imol].get_suggested_initial_contour_level();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return l;

}

//! get the mesh for ligand validation vs dictionary, coloured by badness.
//! greater then 3 standard deviations is fully red.
//! Less than 0.5 standard deviations is fully green.
coot::simple_mesh_t
molecules_container_t::get_mesh_for_ligand_validation_vs_dictionary(int imol, const std::string &ligand_cid) {

   coot::simple_mesh_t m;
   if (is_valid_model_molecule(imol)) {
      m = molecules[imol].get_mesh_for_ligand_validation_vs_dictionary(ligand_cid, geom, thread_pool);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return m;

}

//! ligand validation - basically we do the same as the above function, but the
//! return type is validation data, not a mesh
//!
//! @return a vector of `geometry_distortion_info_container_t`
std::vector<coot::geometry_distortion_info_pod_container_t>
molecules_container_t::get_ligand_validation_vs_dictionary(int imol,
                                                           const std::string &ligand_cid,
                                                           bool with_nbcs) {

   std::vector<coot::geometry_distortion_info_pod_container_t> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].geometric_distortions_for_one_residue_from_mol(ligand_cid, with_nbcs, geom, thread_pool);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}

//! General fragment distortion analysis
//!
//! @param imol is the model molecule index
//! @param selection_cid is the selection CID e.g "//A/15-23"
//! @param include_non_bonded_contacts is the flag to include non bonded contacts
//!
//! @return a vector/list of interesting geometry
std::vector<coot::geometry_distortion_info_pod_container_t>
molecules_container_t::get_validation_vs_dictionary_for_selection(int imol,
                                                                  const std::string &selection_cid,
                                                                  bool include_non_bonded_contacts) {

   std::vector<coot::geometry_distortion_info_pod_container_t> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].geometric_distortions_for_selection_from_mol(selection_cid,
                                                                       include_non_bonded_contacts,
                                                                       geom,
                                                                       thread_pool);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}


//! Get ligand distortion
//!
//! a more simple interface to the above
//!
//! @return a pair: the first is the status (1 for OK, 0 for fail)
//!
// should this be const?
std::pair<int, double>
molecules_container_t::get_ligand_distortion(int imol, const std::string &ligand_cid, bool with_nbcs) {

   int status = 0;
   double d = 0;
   if (is_valid_model_molecule(imol)) {
      std::pair<int, double> p =
         molecules[imol].simple_geometric_distortions_from_mol(ligand_cid, with_nbcs, geom, thread_pool);
      status = p.first;
      d = p.second;
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return std::make_pair(status, d);
}



//! set the map saturation
void
molecules_container_t::set_map_colour_saturation(int imol, float s) {

   if (is_valid_map_molecule(imol)) {
      molecules[imol].set_map_colour_saturation(s);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid map molecule " << imol << std::endl;
   }
}


//! @return the map histogram
coot::molecule_t::histogram_info_t
molecules_container_t::get_map_histogram(int imol, unsigned int n_bins, float zoom_factor) const {

   coot::molecule_t::histogram_info_t hi;
   if (is_valid_map_molecule(imol)) {
      hi = molecules[imol].get_map_histogram(n_bins, zoom_factor);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a map model molecule " << imol << std::endl;
   }
   return hi;

}


//! read extra restraints (e.g. from ProSMART)
int
molecules_container_t::read_extra_restraints(int imol, const std::string &file_name) {

   int n = -1;
   if (is_valid_model_molecule(imol)) {
      n = molecules[imol].read_extra_restraints(file_name);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return n;
}


#include "coot-utils/find-water-baddies.hh"

//! check waters, implicit OR
//! return a vector of atom specifiers
std::vector <coot::atom_spec_t>
molecules_container_t::find_water_baddies(int imol_model, int imol_map,
                                          float b_factor_lim,
                                          float outlier_sigma_level,
                                          float min_dist, float max_dist,
                                          bool ignore_part_occ_contact_flag,
                                          bool ignore_zero_occ_flag) {

   std::vector <coot::atom_spec_t> v;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {

         float map_sigma = molecules[imol_map].get_map_rmsd_approx();
         v = coot::find_water_baddies_OR(molecules[imol_model].atom_sel,
                                         b_factor_lim,
                                         molecules[imol_map].xmap,
                                         map_sigma,
                                         outlier_sigma_level,
                                         min_dist, max_dist,
                                         ignore_part_occ_contact_flag,
                                         ignore_zero_occ_flag);

         std::cout << "........... find_water_baddies_OR() returned " << v.size() << " water baddies " << std::endl;
      } else {
         std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_model << std::endl;
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_map << std::endl;
   }
   return v;

}

//! @return the dictionary read for the give residue type, return an empty string on failure
//! to lookup the residue type
std::string
molecules_container_t::get_cif_file_name(const std::string &comp_id, int imol_enc) const {

   std::string fn = geom.get_cif_file_name(comp_id, imol_enc);
   return fn;
}

//! @return a string that is the contents of a dictionary cif file
std::string
molecules_container_t::get_cif_restraints_as_string(const std::string &comp_id, int imol_enc) const {

   // make this a util function, or a class function at least
   auto file_to_string = [] (const std::string &file_name) {
      std::string s;
      std::string line;
      std::ifstream f(file_name.c_str());
      if (!f) {
         std::cout << "get_cif_restraints_as_string(): Failed to open " << file_name << std::endl;
      } else {
         while (std::getline(f, line)) {
            s += line;
            s += "\n";
         }
      }
      return s;
   };

   std::string r;
   std::pair<bool, coot::dictionary_residue_restraints_t> r_p =
      geom.get_monomer_restraints(comp_id, imol_enc);

   if (r_p.first) {
      const auto &dict = r_p.second;
      std::string fn("tmp.cif");
      dict.write_cif(fn);
      if (coot::file_exists(fn)) {
         r = file_to_string(fn);
      }
   }
   return r;
}


//! @return a list of residues specs that have atoms within dist of the atoms of the specified residue
std::vector<coot::residue_spec_t>
molecules_container_t::get_residues_near_residue(int imol, const std::string &residue_cid, float dist) const {

   std::vector<coot::residue_spec_t> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].residues_near_residue(residue_cid, dist);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;

}


//! Get the chains that are related by NCS:
std::vector<std::vector<std::string> >
molecules_container_t::get_ncs_related_chains(int imol) const {

   std::vector<std::vector<std::string> > v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_ncs_related_chains();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}


//! @return the moldel molecule imol as a string. Return emtpy string on error
std::string
molecules_container_t::molecule_to_PDB_string(int imol) const {

   std::string s;
   if (is_valid_model_molecule(imol)) {
      s = molecules[imol].molecule_to_PDB_string();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return s;
}

//! @return the moldel molecule imol as a string. Return emtpy string on error
std::string
molecules_container_t::molecule_to_mmCIF_string(int imol) const {

   std::string s;
   if (is_valid_model_molecule(imol)) {
      s = molecules[imol].molecule_to_mmCIF_string();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return s;
}


//! return the hb_tye for the given atom. On failure return an empty string
std::string
molecules_container_t::get_hb_type(const std::string &compound_id, int imol_enc, const std::string &atom_name) const {

   coot::hb_t hbt = geom.get_h_bond_type(atom_name, compound_id, imol_enc);
   std::string hb;
   if (hbt == coot::HB_UNASSIGNED) hb = "HB_UNASSIGNED";
   if (hbt == coot::HB_NEITHER)    hb = "HB_NEITHER";
   if (hbt == coot::HB_DONOR)      hb = "HB_DONOR";
   if (hbt == coot::HB_ACCEPTOR)   hb = "HB_ACCEPTOR";
   if (hbt == coot::HB_BOTH)       hb = "HB_BOTH";
   if (hbt == coot::HB_HYDROGEN)   hb = "HB_HYDROGEN";
   return hb;
}


#include "utils/coot-utils.hh"

//! set the maximum number of threads in a thread pool and vector of threads
void
molecules_container_t::set_max_number_of_threads(unsigned int n_threads) {
   coot::set_max_number_of_threads(n_threads);
   thread_pool.resize(n_threads);
}

// call the above function
void
molecules_container_t::set_max_number_of_threads_in_thread_pool(unsigned int n_threads) {
   set_max_number_of_threads(n_threads);
}


//! get the time to run test test function in miliseconds
double
molecules_container_t::test_the_threading(int n_threads) {

   auto reference_data = [] (const std::string &file) {
      char *env = getenv("MOORHEN_TEST_DATA_DIR");
      if (env) {
         std::string joined = coot::util::append_dir_file(env, file);
         return joined;
      } else {
         return file;
      }
   };

   int imol_map = read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);
   coot::set_max_number_of_threads(n_threads);
   float radius = 50;
   auto tp_0 = std::chrono::high_resolution_clock::now();
   coot::simple_mesh_t map_mesh = get_map_contours_mesh(imol_map, 55,10,10, radius, 0.12);
   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   close_molecule(imol_map);
   return d10;
}

double
molecules_container_t::test_launching_threads(unsigned int n_threads_per_batch, unsigned int n_batches) const {

   auto sum = [] (unsigned int i, unsigned int j) {
      return i+j;
   };

   if (n_threads_per_batch == 0) {
      return -1.0;
   } else {
      if (n_batches == 0) {
         return -2.0;
      } else {
         auto tp_0 = std::chrono::high_resolution_clock::now();
         for (unsigned int i=0; i<n_batches; i++) {
            std::vector<std::thread> threads;
            for (unsigned int j=0; j<n_threads_per_batch; j++)
               threads.push_back(std::thread(sum, i, j));
            for (unsigned int j=0; j<n_threads_per_batch; j++)
               threads[j].join();
         }
         auto tp_1 = std::chrono::high_resolution_clock::now();
         auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
         double time_per_patch = d10/static_cast<double>(n_batches);
         return time_per_patch;
      }
   }
}

//! @return time in microsections
double
molecules_container_t::test_thread_pool_threads(unsigned int n_threads) {

   auto sum = [] (unsigned int thread_index, unsigned int i, unsigned int j, std::atomic<unsigned int> &done_count_for_threads) {
      done_count_for_threads++;
      return i+j;
   };

   double t = 0;
   auto tp_0 = std::chrono::high_resolution_clock::now();
   std::atomic<unsigned int> done_count_for_threads(0);

   for (unsigned int i=0; i<n_threads; i++) {
      thread_pool.push(sum, i, i, std::ref(done_count_for_threads));
   }
   while (done_count_for_threads < n_threads)
      std::this_thread::sleep_for(std::chrono::nanoseconds(300));

   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
   t = d10;
   return t;

}

namespace mmcif_tests {
   int run_tests(bool last_test_only);
}

//! a test for mmdb/gemmi/mmcif functionality
int
molecules_container_t::mmcif_tests(bool last_test_only) {

   int status = mmcif_tests::run_tests(last_test_only);
   return status;

}

#include "coot-utils/cfc.hh"

void
molecules_container_t::test_function(const std::string &s) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   // test cfc here.

   // when extracting/reworking this for a real chapi function,
   // pass the output files names for the features clusters and
   // the water cluster (and maybe residue clusters later).
   // the input will be a vector of pair of molecule indices
   // and ligand residue types - and maybe a centre position.

   // --------------------- main line -------------------

   std::vector<std::pair<std::string, std::string> > mol_info =
      { std::pair("brd1/5PB8.pdb", "AC6"),
	std::pair("brd1/5PB9.pdb", "53C"),
	std::pair("brd1/5POW.pdb", "8UA"),
	std::pair("brd1/5PBF.pdb", "8HJ"),
	std::pair("brd1/5PBD.pdb", "TYZ"),
	std::pair("brd1/5PBE.pdb", "TYL"),
	std::pair("brd1/5PBA.pdb", "53B"),
	std::pair("brd1/5PB7.pdb", "8H4"),
	std::pair("brd1/5PBB.pdb", "ES1"),
	std::pair("brd1/5PBC.pdb", "ES3")};

   std::vector<cfc::input_info_t> mol_infos;

   for (const auto &m : mol_info) {
      std::string fn = m.first;
      int imol = read_coordinates(fn);
      mmdb::Manager *mol = get_mol(imol);
      std::string res_name = m.second;
      cfc::input_info_t mi(mol, imol, res_name);
      mol_infos.push_back(mi);
   }

   auto cfc = cfc::chemical_feature_clustering(mol_infos, geom);

#endif // MAKE_ENHANCED_LIGAND_TOOLS

}


//! @return a vector of string pairs that were part of a gphl_chem_comp_info.
//!  return an empty vector on failure to find any such info.
std::vector<std::pair<std::string, std::string> >
molecules_container_t::get_gphl_chem_comp_info(const std::string &compound_id, int imol_enc) {

   std::vector<std::pair<std::string, std::string> > v;
   std::pair<bool, coot::dictionary_residue_restraints_t> r_p =
      geom.get_monomer_restraints(compound_id, imol_enc);
   if (r_p.first) {
      v = r_p.second.gphl_chem_comp_info.info;
   }
   return v;
}

//! get a list of atom names and their associated atedrg atom types
//!
//! @return a list of atom names and their associated atedrg atom types, return an empty list
//! on failure (atoms types are not in the dictionary or atom failure to look up the compound id)l
std::vector<std::pair<std::string, std::string> >
molecules_container_t::get_acedrg_atom_types(const std::string &compound_id, int imol_enc) const {

   std::vector<std::pair<std::string, std::string> > v;
   std::pair<bool, coot::dictionary_residue_restraints_t> r_p =
      geom.get_monomer_restraints(compound_id, imol_enc);
   if (r_p.first) {
      const auto &restraints = r_p.second;
      const auto &atom_info = restraints.atom_info;
      for (unsigned int iat=0; iat<atom_info.size(); iat++) {
         const auto &atom = atom_info[iat];
         const auto &atom_id = atom.atom_id;
         const auto &acedrg_atom_type = atom.acedrg_atom_type;
         if (! acedrg_atom_type.empty()) {
            auto pair = std::make_pair(atom_id, acedrg_atom_type);
            v.push_back(pair);
         }
      }
   }
   return v;

}

//! get acedrg types for ligand bonds
//! @return a vector of `acedrg_types_for_residue_t`
coot::acedrg_types_for_residue_t
molecules_container_t::get_acedrg_atom_types_for_ligand(int imol, const std::string &residue_cid) const {

   coot::acedrg_types_for_residue_t types;

   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *residue_p = molecules[imol].get_residue(residue_cid);
      if (residue_p) {
         int imol_enc = imol;
         types = coot::get_acedrg_types_for_residue(residue_p, imol_enc, geom);
      }
   }
   return types;
}




//! export map molecule as glTF
void
molecules_container_t::export_map_molecule_as_gltf(int imol, float pos_x, float pos_y, float pos_z, float radius, float contour_level,
                                                   const std::string &file_name) {

   if (is_valid_map_molecule(imol)) {
      clipper::Coord_orth pos(pos_x, pos_y, pos_z);
      molecules[imol].export_map_molecule_as_gltf(pos, radius, contour_level, file_name);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid map molecule " << imol << std::endl;
   }


}

//! export model molecule as glTF - This API will change - we want to specify surfaces and ribbons too.
void
molecules_container_t::export_model_molecule_as_gltf(int imol,
                                                     const std::string &selection_cid,
                                                     const std::string &mode,
                                                     bool against_a_dark_background,
                                                     float bonds_width, float atom_radius_to_bond_width_ratio, int smoothness_factor,
                                                     bool draw_hydrogen_atoms_flag, bool draw_missing_residue_loops,
                                                     const std::string &file_name) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].export_model_molecule_as_gltf(mode, selection_cid, &geom, against_a_dark_background,
                                                    bonds_width, atom_radius_to_bond_width_ratio, smoothness_factor,
                                                    draw_hydrogen_atoms_flag, draw_missing_residue_loops, file_name);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}

void
molecules_container_t::export_molecular_representation_as_gltf(int imol, const std::string &atom_selection_cid,
                                                              const std::string &colour_scheme, const std::string &style,
                                                              int secondary_structure_usage_flag,
                                                              const std::string &file_name) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].export_molecular_representation_as_gltf(atom_selection_cid, colour_scheme, style,
                                                             secondary_structure_usage_flag, file_name);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}

//! export chemical features for the specified residue
//!
void molecules_container_t::export_chemical_features_as_gltf(int imol, const std::string &cid, const std::string &file_name) const {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      molecules[imol].export_chemical_features_as_gltf(cid, geom, file_name);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}


//! set the gltf PBR roughness factor
//!
//! @param imol is the model molecule index
//! @param roughness_factor is the factor for the roughness (0.0 to 1.0)
void
molecules_container_t::set_gltf_pbr_roughness_factor(int imol, float roughness_factor) {

   bool is_valid = false;
   if (is_valid_model_molecule(imol)) is_valid =  true;
   if (is_valid_map_molecule(imol)) is_valid =  true;
   if (is_valid) {
      molecules[imol].gltf_pbr_roughness = roughness_factor;
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}

//! set the gltf PBR metalicity factor
//!
//! @param imol is the model molecule index
//! @param metalicity is the factor for the roughness (0.0 to 1.0)
void
molecules_container_t::set_gltf_pbr_metalicity_factor(int imol, float metalicity) {

   bool is_valid = false;
   if (is_valid_model_molecule(imol)) is_valid =  true;
   if (is_valid_map_molecule(imol)) is_valid =  true;
   if (is_valid) {
      molecules[imol].gltf_pbr_metalicity = metalicity;
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}

//! get density at position
//! @return density value
float
molecules_container_t::get_density_at_position(int imol_map, float x, float y, float z) const {

   float f = -1;
   if (is_valid_map_molecule(imol_map)) {
      clipper::Coord_orth pt(x,y,z);
      f = molecules[imol_map].get_density_at_position(pt);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_map << std::endl;
   }
   return f;
}


   //! get residue name
std::string
molecules_container_t::get_residue_name(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) const {

   std::string n;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec(chain_id, res_no, ins_code);
      n = molecules[imol].get_residue_name(res_spec);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return n;


}


//! Get the SMILES string for the give residue type
//!
//! @param residue name the compound-id
//! @param is the molecule index for the residue type/compound_id
//! @return the SMILES string if the residue type can be foound in the dictionary store
//!         or the empty string on a failure.
std::string
molecules_container_t::get_SMILES_for_residue_type(const std::string &residue_name, int imol_enc) const {

   std::string s = geom.Get_SMILES_for_comp_id(residue_name, imol_enc);
   return s;
}


//! @return an estimate of the diameter of the model molecule (-1 on failure)
float
molecules_container_t::get_molecule_diameter(int imol) const {

   float r = -1;
   if (is_valid_model_molecule(imol)) {
      r = molecules[imol].get_molecule_diameter();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return r;


}



//! the caller has access to a compressed file that contains the rotamer probabilities.
//! libcootapi will fill the rotamer probabilities tables from this compressed data stream.
//! (placeholder only)
void
molecules_container_t::accept_rotamer_probability_tables_compressed_data(const std::string &data_stream) {

   // now do something with that data stread - it need not be written to disk

   // uncompress it first

}


//! Interactive B-factor refinement (fun).
//! "factor" might typically be say 0.9 or 1.1
void
molecules_container_t::multiply_residue_temperature_factors(int imol, const std::string &cid, float factor) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].multiply_residue_temperature_factors(cid, factor);
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

}


//! change the chain id
//! @return -1 on a conflict
//! 1 on good.
//! 0 on did nothing
//! return also an information/error message
std::pair<int, std::string>
molecules_container_t::change_chain_id(int imol,
                                       const std::string &from_chain_id,
                                       const std::string &to_chain_id,
                                       bool use_resno_range,
                                       int start_resno, int end_resno) {

   std::pair<int, std::string> status(0, "");
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].change_chain_id(from_chain_id, to_chain_id, use_resno_range, start_resno, end_resno);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

//! split an NMR model into multiple models - all in MODEL 1.
//! @return the vector of new molecule indices.
std::vector<int>
molecules_container_t::split_multi_model_molecule(int imol) {

   std::vector<int> v;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = get_mol(imol);
      if (mol) {
         std::vector<mmdb::Manager *> mv = coot::util::split_multi_model_molecule(mol);
         for (unsigned int i=0; i<mv.size(); i++) {
            auto asc = make_asc(mv[i]);
            std::string name = "split-molecule" + std::to_string(i+1);
            int idx = molecules.size();
            molecules.push_back(coot::molecule_t(asc, idx, name));
            v.push_back(idx);
         }
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}


//! make a multi-model molecule given the input molecules
//! ``model_molecules_list`` is a colon-separated list of molecules, *e.g.* "2:3:4"
//! @return the new molecule index - -1 if no models were found in the ``model_molecules_list``
int
molecules_container_t::make_ensemble(const std::string &model_molecule_list) {

   // make this be a member function.
   // Add a test for valid model molecule when you do so
   auto model_molecule_string_list_to_molecule_index_vec = [] (const std::string &model_molecule_list) {
      std::vector<int> mols;
      std::vector<std::string> number_strings = coot::util::split_string(model_molecule_list, ":");
      for (const auto &item : number_strings) {
         int idx = coot::util::string_to_int(item);
         mols.push_back(idx);
      }
      return mols;
   };

   int imol_new = -1;
   mmdb::Manager *mol_sumo = new mmdb::Manager;
   std::vector<int> mols = model_molecule_string_list_to_molecule_index_vec(model_molecule_list);
   unsigned int n_done = 0;
   for (unsigned int i=0; i<mols.size(); i++) {
      unsigned int idx = mols[i];
      if (is_valid_model_molecule(idx)) {
         mmdb::Manager *mol = molecules[idx].atom_sel.mol;
         if (mol) {
            for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
               mmdb::Model *model_p = mol->GetModel(imod);
               mmdb::Model *new_model = new mmdb::Model;
               new_model->Copy(model_p);
               mol_sumo->AddModel(new_model);
               n_done++;
            }
         }
      }
   }

   if (n_done > 0) {

      // now create a new molecule
      std::string name = "Ensemble " + model_molecule_list;
      imol_new = molecules.size();
      atom_selection_container_t asc = make_asc(mol_sumo);
      coot::molecule_t m(asc, imol_new, name);
      molecules.push_back(m);

   } else {
      // clean up
      delete mol_sumo;
   }
   return imol_new;
}


//! Fourier Shell Correlation (FSC) between maps
//! @return a vector or pairs of graph points (resolution, correlation)
std::vector<std::pair<double, double> >
molecules_container_t::fourier_shell_correlation(int imol_map_1, int imol_map_2) const {

   std::vector<std::pair<double, double> > v;

   if (is_valid_map_molecule(imol_map_1)) {
      if (is_valid_map_molecule(imol_map_2)) {
         const clipper::Xmap<float> &xmap_1 = molecules[imol_map_1].xmap;
         const clipper::Xmap<float> &xmap_2 = molecules[imol_map_2].xmap;
         auto fsc = coot::util::fsc(xmap_1, xmap_2);
         if (! fsc.empty()) {
            v.resize(fsc.size());
            for (unsigned int i=0; i<fsc.size(); i++) {
               v[i].first  = fsc[i].first.invresolsq_limit();
               v[i].second = fsc[i].second;
            }
         }
      }
   }
   return v;
}


//! not const because it can dynamically add dictionaries
coot::atom_overlaps_dots_container_t
molecules_container_t::get_overlap_dots(int imol) {

   coot::atom_overlaps_dots_container_t aodc;
   if (is_valid_model_molecule(imol)) {
      aodc = molecules[imol].get_overlap_dots(&geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return aodc;
}

//! not const because it can dynamically add dictionaries
coot::atom_overlaps_dots_container_t
molecules_container_t::get_overlap_dots_for_ligand(int imol, const std::string &cid_ligand) {

   coot::atom_overlaps_dots_container_t aodc;
   if (is_valid_model_molecule(imol)) {
      aodc = molecules[imol].get_overlap_dots_for_ligand(cid_ligand, &geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return aodc;
}



//! not const because it can dynamically add dictionaries
std::vector<coot::plain_atom_overlap_t>
molecules_container_t::get_atom_overlaps(int imol) {

   std::vector<coot::plain_atom_overlap_t> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_atom_overlaps(&geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}

//! Get the atom overlap score
//!
//! @param imol the model molecule index
//! @return the overlap score - a negative number indicates failure
float
molecules_container_t::get_atom_overlap_score(int imol) {

   float v = -1.0;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_atom_overlap_score(&geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

   return v;

}



//! not const because it can dynamically add dictionaries
std::vector<coot::plain_atom_overlap_t>
molecules_container_t::get_overlaps_for_ligand(int imol, const std::string &cid_ligand) {

   std::vector<coot::plain_atom_overlap_t> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_overlaps_for_ligand(cid_ligand, &geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}

void
molecules_container_t::print_secondary_structure_info(int imol) const {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].print_secondary_structure_info();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}


//! copy the dictionary that is specific for imol_current so that it can be used with imol_new
bool
molecules_container_t::copy_dictionary(const std::string &monomer_name, int imol_current, int imol_new) {

   std::cout << "--------------------------   debug:: calling copy_monomer_restraints() "
             << monomer_name << " " << imol_current << " " << imol_new << std::endl;
   bool status = geom.copy_monomer_restraints(monomer_name, imol_current, imol_new);

   std::pair<bool, coot::dictionary_residue_restraints_t> r =
      geom.get_monomer_restraints(monomer_name, imol_new);

   std::cout << "-------------- r " << r.first << std::endl;
   
   return status;

}


#include "coot-utils/pae.hh"
//! @return a string of a png
std::string
molecules_container_t::pae_png(const std::string &pae_file_name) const {

// This test acts the way we want it to, but it's not a good name
// something like "HAVE_CAIRO" would be prefered.
//
#if RDKIT_HAS_CAIRO_SUPPORT // Cairo is not allowed in Moorhen.
   int n_pixels = 600;
   pae_t pae(pae_file_name, n_pixels);
   return pae.get_image();
#else
   return "no-cairo";
#endif

}


//! get the median temperature factor for the model
//! @return a negative number on failure.
float
molecules_container_t::get_median_temperature_factor(int imol) const {

   float b_factor = -1.1;
   if (is_valid_model_molecule(imol)) {
      b_factor = molecules[imol].get_median_temperature_factor();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return b_factor;
}

//! Get the atom temperature factor
//!
//! @param imol is the model molecule index
//! @param atom_cid is the selection cid for the atom
//!
//! @return a negative number on failure, otherwise the temperature factor
float
molecules_container_t::get_temperature_factor_of_atom(int imol, const std::string &atom_cid) const {

   float b_factor = -1.1;
   if (is_valid_model_molecule(imol)) {
      b_factor = molecules[imol].get_temperature_factor_of_atom(atom_cid);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return b_factor;
}



// return the atom name match on superposing the atoms of the given dictionaries
std::map<std::string, std::string>
molecules_container_t::dictionary_atom_name_map(const std::string &comp_id_1, int imol_1, const std::string &comp_id_2, int imol_2) {

   std::map<std::string, std::string> m;

   std::pair<bool, coot::dictionary_residue_restraints_t> r_p_1 = geom.get_monomer_restraints(comp_id_1, imol_1);
   std::pair<bool, coot::dictionary_residue_restraints_t> r_p_2 = geom.get_monomer_restraints(comp_id_2, imol_2);
   if (r_p_1.first) {
      if (r_p_2.first) {
         const coot::dictionary_residue_restraints_t &dict_1 = r_p_1.second;
         const coot::dictionary_residue_restraints_t &dict_2 = r_p_2.second;
         coot::dictionary_match_info_t dm = dict_1.match_to_reference(dict_2, nullptr, comp_id_1, "dummy");
         if (false) {
            std::cout << "There are " << dm.same_names.size() << " atoms with the same name" << std::endl;
            std::cout << "There are " << dm.name_swaps.size() << " atoms with the different names" << std::endl;
         }
         for (const auto &name : dm.same_names)
            m[name] = name;
         for (const auto &name : dm.name_swaps) {
            m[name.first] = name.second;
         }
      }
   }
   return m;
}



//! Get the residue CA position
//!
//! @return a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values
std::vector<double>
molecules_container_t::get_residue_CA_position(int imol, const std::string &cid) const {

   std::vector<double> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_residue_CA_position(cid);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;

}

//! Get the average residue position
//!
//! @return a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values
std::vector<double>
molecules_container_t::get_residue_average_position(int imol, const std::string &cid) const {

   std::vector<double> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_residue_average_position(cid);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}

//! Get the average residue side-chain position
//!
//! @return a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values
std::vector<double>
molecules_container_t::get_residue_sidechain_average_position(int imol, const std::string &cid) const {

   std::vector<double> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_residue_sidechain_average_position(cid);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}

//! Get the torsion of the specified atom in the specified residue
//!
//! @param imol is the model molecule index
//! @param cid is the selection CID, e.g. //A/15 (residue 15 in chain A)
//! @param atom_names is a list of atom names, e.g. ["CA", "CB", "CG", "CD"]
//!
//! @return a pair, the first of which is a succes status (1 success, 0 failure), the second is the torsion in degrees
std::pair<int, double>
molecules_container_t::get_torsion(int imol, const std::string &cid, const std::vector<std::string> &atom_names) {
   std::pair<int, double> p(0,0);

   if (is_valid_model_molecule(imol)) {
      p = molecules[imol].get_torsion(cid, atom_names);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return p;
}

//! set the occupancy for the given atom selection
//!
//! @param imol is the model molecule index
//! @param cod is the atom selection CID
//! @param is the new occupancy
void
molecules_container_t::set_occupancy(int imol, const std::string &cid, float occ_new) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].set_occupancy(cid, occ_new);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

}

#include "utils/subprocess.hpp"


//! External refinement using servalcat
int
molecules_container_t::servalcat_refine_xray(int imol, int imol_map, const std::string &output_prefix) {

   std::map<std::string, std::string> kvm;
   return servalcat_refine_xray_internal(imol, imol_map, output_prefix, kvm);

}


int
molecules_container_t::servalcat_refine_xray_internal(int imol, int imol_map, const std::string &output_prefix,
                                                      const std::map<std::string, std::string> &key_value_pairs) {

   int imol_refined_model = -1;

   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {

         bool set_weight = false;
         std::string weight_str;
         if (! key_value_pairs.empty()) {
            for (const auto &kv : key_value_pairs) {
               if (kv.first == "weight") {
                  set_weight = true;
                  weight_str = kv.second;
               }
            }
         }

         bool clibd_mon_is_set = false;
         char *e = getenv("CLIBD_MON");
         if (e) {
            std::string env(e);
            if (std::filesystem::exists(env))
               clibd_mon_is_set = true;
         }
         if (clibd_mon_is_set) {
            std::string mtz_file   = molecules[imol_map].refmac_mtz_filename;
            std::string fobs_col   = molecules[imol_map].refmac_fobs_col;
            std::string sigfob_col = molecules[imol_map].refmac_sigfobs_col;
            std::string r_free_col = molecules[imol_map].refmac_r_free_col;
            if (! mtz_file.empty()) {

               bool read_pdb_output = false; // this gets set to true if the output pdb is sane

               std::string c(",");
               std::string labin = fobs_col + c + sigfob_col + c + r_free_col;

               std::string dir_1 = "coot-servalcat";
               coot::util::create_directory(dir_1);
               std::string prefix = coot::util::append_dir_file(dir_1, output_prefix);
               std::string  input_pdb_file_name = prefix + std::string("-in.pdb");
               std::string output_pdb_file_name = prefix + std::string(".pdb"); // named by servalcat
               int status = molecules[imol].write_coordinates(input_pdb_file_name);
               // see https://www.ebi.ac.uk/pdbe/docs/cldoc/object/cl_obj_rdwr.html#CMMDBManager::WritePDBASCII
               if (status == 0) {
                  bool output_pdb_file_name_exists = false;
                  std::filesystem::file_time_type output_pdb_file_name_time;
                  std::filesystem::path p(output_pdb_file_name);
                  if (std::filesystem::exists(p)) {
                     output_pdb_file_name_exists = true;
                     output_pdb_file_name_time = std::filesystem::last_write_time(p);
                  }
                  std::vector<std::string> cmd_list = {"servalcat", "refine_xtal_norefmac",
                                                       "-s", "xray", "--model", input_pdb_file_name,
                                                       "--hklin", mtz_file, "--labin", labin,
                                                       "-o", prefix};
                  if (set_weight) {
                     cmd_list.push_back("--weight");
                     cmd_list.push_back(weight_str);
                  }

                  if (true) {
                     std::cout << "commandline: ";
                     for (unsigned int i=0; i<cmd_list.size(); i++) std::cout << " " << cmd_list[i];
                     std::cout << "\n";
                  }
                  std::cout << "running servalcat..." << std::endl;
                  subprocess::OutBuffer obuf = subprocess::check_output(cmd_list);
                  if (std::filesystem::exists(p)) {
                     if (output_pdb_file_name_exists) {
                        std::filesystem::file_time_type new_output_pdb_file_name_time = std::filesystem::last_write_time(p);
                        auto t1 =     output_pdb_file_name_time.time_since_epoch();
                        auto t2 = new_output_pdb_file_name_time.time_since_epoch();
                        auto tt1 = std::chrono::duration_cast<std::chrono::seconds>(t1).count();
                        auto tt2 = std::chrono::duration_cast<std::chrono::seconds>(t2).count();
                        auto d = tt2 - tt1;
                        if (d > 0)
                           read_pdb_output = true;
                     } else {
                        read_pdb_output = true;
                     }
                     if (read_pdb_output) {
                        imol_refined_model = read_coordinates(output_pdb_file_name);
                     }
                  } else {
                     std::cout << "WARNING:: " << __FUNCTION__ << "(): path does not exist " << p << std::endl;
                  }
               } else {
                  std::cout << "WARNING::" << __FUNCTION__ << "(): bad status on writing servalcat input file" << std::endl;
               }
            } else {
               std::cout << "WARNING::" << __FUNCTION__ << "(): mtz file_name was empty" << std::endl;
            }
         } else {
            std::cout << "WARNING::" << __FUNCTION__ << "(): CLIBD_MON was not set correctly" << std::endl;
         }
      } else {
         std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_map << std::endl;
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return imol_refined_model;
}


// void
// molecules_container_t::servalcat_refine_cryoem(int imol,
//                                                const std::string &half_map_1,
//                                                const std::string &half_map_2,
//                                                const std::string &mask,
//                                                const std::string &output_prefix);


//! split a residue into alt-confs
//!
//! do nothing if the residue already has alt-confs.
//!
//! @param imol the modified model
//! @param residue_cid the modified residue
//! @param the difference map that is used to determine the residue split
//! @return split success status
int
molecules_container_t::split_residue_using_map(int imol, const std::string &residue_cid, int imol_diff_map) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_diff_map)) {
         mmdb::Residue *residue_p = molecules[imol].get_residue(residue_cid);
         if (residue_p) {
            mmdb::Manager *mol = get_mol(imol);
            if (mol) {
               const clipper::Xmap<float> &xmap = molecules[imol_diff_map].xmap;
               coot::util::split_residue_using_map(residue_p, mol, xmap);
            }
         }
      }
   }
   return status;
}


std::vector<coot::residue_range_t>
molecules_container_t::get_missing_residue_ranges(int imol) const {

   std::vector<coot::residue_range_t> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_missing_residue_ranges();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}


//! Get the sequence information
//!
//! @param imol is the molecule index
//! @return the sequence information
std::vector<std::pair<std::string, std::string> >
molecules_container_t::get_sequence_info(int imol) const {

   std::vector<std::pair<std::string, std::string> > v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_sequence_info();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}


//! return the mismatches/mutations:
coot::chain_mutation_info_container_t
molecules_container_t::get_mutation_info(int imol) const {

  coot::chain_mutation_info_container_t mci;
  if (is_valid_model_molecule(imol)) {
    mci = molecules[imol].get_mutation_info();
  } else {
    std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
  }
  return mci;
}

//! Change the B factors
//!
//! @param imol is the model molecule index
//! @param cid is the selection CID, e.g. //A/15 (residue 15 in chain A)
//! @param temp_fact is the isotropic ADP/temperature factor, e.g.,  22
void
molecules_container_t::set_temperature_factors_using_cid(int imol, const std::string &cid, float temp_fact) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].set_temperature_factors_using_cid(cid, temp_fact);
   }
}

//! Residue is nucleic acid?
//!
//! Every residue in the selection is checked
//!
//! @param imol is the model molecule index
//! @param cid is the selection CID e.g "//A/15" (residue 15 of chain A)
//!
//! @return a bool
bool
molecules_container_t::residue_is_nucleic_acid(int imol, const std::string &cid) const {

   bool status = false;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].residue_is_nucleic_acid(cid);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}


//! Get the residue type
//!
//! @param imol is the model molecule index
//! @param cid is the selection CID e.g "//A/16" (residue 16 of chain A)
//! @return a string. Return an empty string on failure
std::string
molecules_container_t::get_residue_type(int imol, const std::string &cid) const {

   std::string r;
   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *res_p = get_residue_using_cid(imol, cid);
      if (res_p) {
         std::string chain_id = res_p->GetChainID();
         std::string ins_code = res_p->GetInsCode();
         int res_no = res_p->GetSeqNum();
         coot::residue_spec_t rs(chain_id, res_no, ins_code);
         r = molecules[imol].get_residue_name(rs);
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return r;
}


//! get atom distances
//! other stuff here
std::vector<coot::atom_distance_t>
molecules_container_t::get_distances_between_atoms_of_residues(int imol,
							       const std::string &cid_res_1,
							       const std::string &cid_res_2,
							       float dist_max) const {
  std::vector<coot::atom_distance_t> v;
  if (is_valid_model_molecule(imol)) {
     v = molecules[imol].get_distances_between_atoms_of_residues(cid_res_1, cid_res_2, dist_max);
  } else {
     std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
  }

  return v;
}


std::vector<std::string>
molecules_container_t::get_types_in_molecule(int imol) const {

   std::vector<std::string> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_types_in_molecule();
   } else {
      // std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
      logger.log(log_t::WARNING, logging::function_name_t(__FUNCTION__),
		 "not a valid model molecule", imol);
   }
   return v;
}


//! Get Radius of Gyration
//!
//! @param imol is the model molecule index
//!
//! @return the molecule centre. If the number is less than zero, there
//! was a problem finding the molecule or atoms.
double
molecules_container_t::get_radius_of_gyration(int imol) const {

   double d = -1.0; // failure
   if (is_valid_model_molecule(imol)) {
      d = molecules[imol].get_radius_of_gyration();
   } else {
      // std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
      logger.log(log_t::WARNING, logging::function_name_t(__FUNCTION__),
		 "not a valid model molecule", imol);
   }
   return d;

}

//! Get atom selection as json
//!
//! @param imol is the model molecule index
//! @param cid is the atom selection CID e.g "//A/15/OH" (atom OH in residue 15 of chain A)
std::string molecules_container_t::get_molecule_selection_as_json(int imol, const std::string &cid) const {

   std::string s;
   if (is_valid_model_molecule(imol)) {
      s = molecules[imol].get_molecule_selection_as_json(cid);
   } else {
      logger.log(log_t::WARNING, logging::function_name_t(__FUNCTION__),
		 "not a valid model molecule", imol);
   }
   return s;
}

//! get pucker info
//!
//! @param imol2 is the model molecule index
//! @return a json string or an empty string on failure
std::string
molecules_container_t::get_pucker_analysis_info(int imol) const {

   std::string pai;
   if (is_valid_model_molecule(imol)) {
      pai = molecules[imol].get_pucker_analysis_info();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return pai;
}

