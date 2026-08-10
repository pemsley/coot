/*
 * api/molecules-container-refine.cc
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


// Refinement, restraints and updating-maps functions of molecules_container_t

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
#include "coot-utils/read-amber-trajectory.hh"
#include "coot-utils/json.hpp"

#include "coords/Bond_lines.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-extras.hh"

#include "mmdb2/mmdb_atom.h"
#include "utils/logging.hh"
extern logging logger;

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
      // std::cout << "INFO:: refine_residues_vec() with altconf \"" << alt_conf << "\"" << std::endl;
      logger.log(log_t::INFO, "refine_residues_vec() with altconf", alt_conf);

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
      // std::cout << "INFO:: using soft harmonic restraints for NBC" << std::endl;
      logger.log(log_t::INFO, "using soft harmonic restraints for NBC");
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
            // std::cout << "INFO:: make_last_restraints() [immediate] restraints locked by "
            //           << restraints_locking_function_name << std::endl;
            logger.log(log_t::INFO, "make_last_restraints() [immediate] restraints locked by", restraints_locking_function_name);
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


