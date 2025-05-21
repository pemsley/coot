/*
 * api/molecules-container-maps.cc
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

#include "compat/coot-sysdep.h"
#include "clipper/core/clipper_types.h"
#include "clipper/core/coords.h"
#include "mmdb2/mmdb_atom.h"
#include "molecules-container.hh"
#include "validation-information.hh"

#include "clipper-ccp4-map-file-wrapper.hh"
#include "coot-utils/slurp-map.hh"

int
molecules_container_t::read_ccp4_map(const std::string &file_name, bool is_a_difference_map) {

   int imol = -1; // currently unset
   int imol_in_hope = molecules.size();
   bool done = false;

   if (! coot::file_exists(file_name)) {
      std::cout << "WARNING:: file does not exist " << file_name << std::endl;
      return imol;
   }

   if (true) {
      if (coot::util::is_basic_em_map_file(file_name) == coot::util::slurp_map_result_t::IS_SLURPABLE_EM_MAP) {
         std::cout << "::::: read_ccp4_map() is_basic_em_map_file() finds a slurpable map " << std::endl;
      } else {
         std::cout << "::::: read_ccp4_map() is_basic_em_map_file() - not a slurpable map " << std::endl;
      }
   }


   if (coot::util::is_basic_em_map_file(file_name) == coot::util::slurp_map_result_t::IS_SLURPABLE_EM_MAP) {

      std::cout << "DEBUG:: mc::read_ccp4_map() returns true for is_basic_em_map_file() "
                << file_name << std::endl;

      // fill xmap
      bool check_only = false;
      short int is_em_map = 1; // this is the correct type - it can be -1.
      coot::molecule_t m(file_name, imol_in_hope, is_em_map);
      short int m_em_status = m.is_EM_map();
      std::cout << "m_em_status " << m_em_status << std::endl;
      clipper::Xmap<float> &xmap = m.xmap;
      coot::util::slurp_map_result_t smr = coot::util::slurp_fill_xmap_from_map_file(file_name, &xmap, check_only);
      if (smr == coot::util::slurp_map_result_t::OK) {
         molecules.push_back(m);
         imol = imol_in_hope;
      }
   }

   if (false) {
      if (is_valid_map_molecule(imol)) {
         short int em_status = molecules[imol].is_EM_map();
         std::cout << "here with imol " << imol << " molecules size " << molecules.size() << std::endl;
         std::cout << "here with imol " << imol << " done " << done << std::endl;
         std::cout << "here with imol " << imol << " is_em_map:  " << em_status << std::endl;
      }
   }

   if (! done) {
      std::cout << "INFO:: attempting to read CCP4 map: " << file_name << " via non-slurp method" << std::endl;
      // clipper::CCP4MAPfile file;
      clipper_map_file_wrapper w_file;
      try {
         w_file.open_read(file_name);

         // em = set_is_em_map(file_name);

         if (true) {
            clipper::Cell fcell = w_file.cell();
            double vol = fcell.volume();
            if (vol < 1.0) {
               std::cout << "WARNING:: read_ccp4_map(): non-sane unit cell volume " << vol << " - skip read"
                         << std::endl;
               // bad_read = true;
            } else {
               try {
                  clipper::CCP4MAPfile file;
                  file.open_read(file_name);
                  clipper::Xmap<float> xmap;
                  file.import_xmap(xmap);
                  if (xmap.is_null()) {
                     std::cout << "ERROR:: failed to read the map" << file_name << std::endl;
                  } else {
                     std::string name = file_name;
                     coot::molecule_t m(name, imol_in_hope);
                     m.xmap = xmap;
                     if (is_a_difference_map)
                        m.set_map_is_difference_map(true);
                     molecules.push_back(m); // oof.
                     imol = imol_in_hope;
                  }
               }
               catch (const clipper::Message_generic &exc) {
                  std::cout << "WARNING:: failed to read " << file_name
                            << " Bad ASU (inconsistant gridding?)." << std::endl;
                  // bad_read = true;
               }
            }
         }
      } catch (const clipper::Message_base &exc) {
         std::cout << "WARNING:: failed to open " << file_name << std::endl;
         // bad_read = true;
         imol = -3; // clipper error
      }
   }
   return imol;
}



coot::validation_information_t
molecules_container_t::density_fit_analysis(int imol_model, int imol_map) const {

   coot::validation_information_t r;
   r.name = "Density fit analysis";
#ifdef EMSCRIPTEN
   r.type = "DENSITY";
#else
   r.type = coot::DENSITY;
#endif
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
         // fill these
         mmdb::PResidue *SelResidues = 0;
         int nSelResidues = 0;

         auto atom_sel = molecules[imol_model].atom_sel;
         int selHnd = atom_sel.mol->NewSelection(); // yes, it's deleted.
         int imod = 1; // multiple models don't work on validation graphs

         atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, imod,
                              "*", // chain_id
                              mmdb::ANY_RES, "*",
                              mmdb::ANY_RES, "*",
                              "*",  // residue name
                              "*",  // Residue must contain this atom name?
                              "*",  // Residue must contain this Element?
                              "*",  // altLocs
                              mmdb::SKEY_NEW // selection key
                              );
         atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

         for (int ir=0; ir<nSelResidues; ir++) {
            mmdb::Residue *residue_p = SelResidues[ir];
            coot::residue_spec_t res_spec(residue_p);
            mmdb::PAtom *residue_atoms=0;
            int n_residue_atoms;
            residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
            double residue_density_score =
               coot::util::map_score(residue_atoms, n_residue_atoms, molecules[imol_map].xmap, 1);
            std::string l = res_spec.label();
            std::string atom_name = coot::util::intelligent_this_residue_mmdb_atom(residue_p)->GetAtomName();
            const std::string &chain_id = res_spec.chain_id;
            int this_resno = res_spec.res_no;
            coot::atom_spec_t atom_spec(chain_id, this_resno, res_spec.ins_code, atom_name, "");
            coot::residue_validation_information_t rvi(res_spec, atom_spec, residue_density_score, l);
            r.add_residue_validation_information(rvi, chain_id);
         }
         atom_sel.mol->DeleteSelection(selHnd);
      }
   }
   r.set_min_max();
   return r;
}

//! @return the sum of the density of the given atoms in the specified CID
//!  return -1001 on failure to find the residue or any atoms in the residue or if imol_map is not a map
double
molecules_container_t::get_sum_density_for_atoms_in_residue(int imol, const std::string &cid,
                                                            const std::vector<std::string> &atom_names,
                                                            int imol_map) {
   double v = 1001.0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         const clipper::Xmap<float> &xmap = molecules.at(imol_map).xmap;
         v = molecules[imol].sum_density_for_atoms_in_residue(cid, atom_names, xmap);
      } else {
         std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_map << std::endl;
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;

}


//! density correlation validation information
coot::validation_information_t
molecules_container_t::density_correlation_analysis(int imol_model, int imol_map) const {

   coot::validation_information_t r;
   r.name = "Density correlation analysis";
#ifdef EMSCRIPTEN
   r.type = "CORRELATION";
#else
   r.type = coot::CORRELATION;
#endif
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {

         mmdb::Manager *mol = molecules[imol_model].atom_sel.mol;
         const clipper::Xmap<float> &xmap = molecules.at(imol_map).xmap;

         unsigned short int atom_mask_mode = 0;
         float atom_radius = 2.0;

         std::vector<coot::residue_spec_t> residue_specs;
         std::vector<mmdb::Residue *> residues = coot::util::residues_in_molecule(mol);
         for (unsigned int i=0; i<residues.size(); i++)
            residue_specs.push_back(coot::residue_spec_t(residues[i]));

         std::vector<std::pair<coot::residue_spec_t, float> > correlations =
            coot::util::map_to_model_correlation_per_residue(mol,
                                                             residue_specs,
                                                             atom_mask_mode,
                                                             atom_radius, // for masking
                                                             xmap);

         std::vector<std::pair<coot::residue_spec_t, float> >::const_iterator it;
         for (it=correlations.begin(); it!=correlations.end(); ++it) {
            const auto &r_spec(it->first);
            const auto &correl(it->second);

            std::string atom_name = " CA ";
            coot::atom_spec_t atom_spec(r_spec.chain_id, r_spec.res_no, r_spec.ins_code, atom_name, "");
            std::string label = "Correl: ";
            coot::residue_validation_information_t rvi(r_spec, atom_spec, correl, label);
            r.add_residue_validation_information(rvi, r_spec.chain_id);
         }

      } else {
         std::cout << "debug:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_map << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_model << std::endl;
   }
   r.set_min_max();
   return r;
}

// This function does no normalisztion of the scales,
// presuming that they are pre-normalized.
// @return the index of the new map, or -1 on failure.
int
molecules_container_t::average_map(const std::string &imol_maps_string, std::vector<float> &scales) {

   int imol_new = -1;
   if (! scales.empty()) {
      std::vector<std::string> number_strings = coot::util::split_string(imol_maps_string, ":");
      std::vector<std::pair<int, float> > weights_and_maps;
      unsigned int scales_index = 0;
      for (const auto &item : number_strings) {
         int idx = coot::util::string_to_int(item);
         if (is_valid_map_molecule(idx)) {
            if (scales_index < scales.size()) {
               weights_and_maps.push_back(std::make_pair(idx, scales[scales_index]));
               scales_index++;
            }
         }
      }

      if (weights_and_maps.size() == scales.size()) {
         // now check that the grids are the same
         int idx = weights_and_maps[0].first;
         clipper::Grid_sampling gs_ref = molecules[idx].xmap.grid_sampling();
         unsigned int n_match = 0;
         for (unsigned int i=0; i<weights_and_maps.size(); i++) {
            idx = weights_and_maps[i].first;
            bool match = false;
            clipper::Grid_sampling gs_this = molecules[idx].xmap.grid_sampling();
            if (gs_ref.nu() == gs_this.nu())
               if (gs_ref.nv() == gs_this.nv())
                  if (gs_ref.nw() == gs_this.nw())
                     match = true;
            if (match) n_match++;
         }
         if (n_match == weights_and_maps.size()) {
            idx = weights_and_maps[0].first;
            std::string name = "Average map " + imol_maps_string;
            bool is_em_map = molecules[idx].is_EM_map();
            imol_new = molecules.size();
            std::vector<std::pair<clipper::Xmap<float>, float> > maps_and_weights(scales.size());
            // heavyweight stuff going on here
            for (unsigned int i=0; i<weights_and_maps.size(); i++) {
               idx = weights_and_maps[i].first;
               std::pair<clipper::Xmap<float>, float> p(molecules[idx].xmap, weights_and_maps[i].second);
               maps_and_weights[i] = p;
            }
            for (unsigned int i=0; i<weights_and_maps.size(); i++) {
	      std::cout << "debug map and weight " << weights_and_maps[i].first << " "
			<< weights_and_maps[i].second << "\n";
	    }
            clipper::Xmap<float> xmap_new = coot::util::average_map(maps_and_weights);
            molecules.push_back(coot::molecule_t(name, imol_new, xmap_new, is_em_map));
         }
      }
   }
   return imol_new;
}

// This function does no normalisztion of the scales, presuming that they are pre-normalized.
// @return success status
bool
molecules_container_t::regen_map(int imol_map, const std::string &imol_maps_string, const std::vector<float> &scales) {

   bool status = false;
   // molecules.push_back(coot::molecule_t(name, imol_new, xmap_ref, is_em_map));

   if (is_valid_map_molecule(imol_map)) {
      if (! scales.empty()) {
         std::vector<std::string> number_strings = coot::util::split_string(imol_maps_string, ":");
         std::vector<std::pair<clipper::Xmap<float> *, float> > maps_and_scales_vec;
         unsigned int scales_index = 0;
         for (const auto &item : number_strings) {
            int idx = coot::util::string_to_int(item);
            if (is_valid_map_molecule(idx)) {
               if (scales_index < scales.size()) {
                  maps_and_scales_vec.push_back(std::make_pair(&molecules[idx].xmap, scales[scales_index]));
                  scales_index++;
               }
            }
         }

         if (maps_and_scales_vec.size() == scales.size()) {
            coot::util::regen_weighted_map(&molecules[imol_map].xmap, maps_and_scales_vec);
            status = true;
         }
      }
   }

   return status;
}

texture_as_floats_t
molecules_container_t::get_map_section_texture(int imol, int section_index, int axis,
                                               float data_value_for_bottom, float data_value_for_top) const {

   texture_as_floats_t t;
   if (is_valid_map_molecule(imol)) {
       t = molecules[imol].get_map_section_texture(section_index, axis, data_value_for_bottom, data_value_for_top);
   }
   return t;
}

//! @return the number of section in the map along the give axis.
//! (0 for X-axis, 1 for y-axis, 2 for Z-axis).
//! return -1 on failure.
int
molecules_container_t::get_number_of_map_sections(int imol_map, int axis_id) const {

   int n = -1;
   if (is_valid_map_molecule(imol_map)) {
      n = molecules[imol_map].get_number_of_map_sections(axis_id);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_map << std::endl;
   }
   return n;
}

#include "coot-utils/xmap-stats.hh"
#include "coot-utils/q-score.hh"

coot::validation_information_t
molecules_container_t::get_q_score_validation_information(mmdb::Manager *mol, int udd_q_score, bool do_per_atom) const {

   coot::validation_information_t vi; // return this
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain = 0; ichain < n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         coot::chain_validation_information_t cvi(chain_p->GetChainID());
         int n_res = chain_p->GetNumberOfResidues();
         for (int ires = 0; ires < n_res; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               int n_atoms = residue_p->GetNumberOfAtoms();

               if (n_atoms > 0) {
                  coot::residue_spec_t residue_spec(residue_p);
                  coot::atom_spec_t spec_for_atom_in_residue(
                      residue_p->GetAtom(0)); // say
                  double qed_residue_sum = 0;
                  unsigned int n_qed_atoms = 0;

                  // accumulate QED for the residue
                  for (int iat = 0; iat < n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     if (!at->isTer()) {
                       mmdb::realtype q;
                       if (at->GetUDData(udd_q_score, q) == mmdb::UDDATA_Ok) {
                         if (q > -1000.0) { // test for not invalid QED on atom
                           qed_residue_sum += q;
                           n_qed_atoms++;
                         }
                       }
                     }
                  }

                  if (n_qed_atoms > 0) {
                     double qed_residue =
                         qed_residue_sum / static_cast<double>(n_qed_atoms);
                     std::string fs = coot::util::float_to_string_using_dec_pl(
                         qed_residue, 2);
                     std::string label = residue_spec.chain_id + " " +
                                         std::to_string(residue_spec.res_no) +
                                         " " + fs;
                     coot::residue_validation_information_t rvi(
                         residue_spec, spec_for_atom_in_residue, qed_residue,
                         label);
                     cvi.add(rvi);

                     if (do_per_atom) {
                       // now the per-atom score
                       for (int iat = 0; iat < n_atoms; iat++) {
                         mmdb::Atom *at = residue_p->GetAtom(iat);
                         if (!at->isTer()) {
                           mmdb::realtype q;
                           at->GetUDData(udd_q_score, q);
                           if (false)
                             std::cout << " " << coot::atom_spec_t(at) << " B "
                                       << at->tempFactor << "  Q-Score: " << q
                                       << std::endl;
                         }
                       }
                     }
                  }
               }
            }
         }
         vi.add(cvi);
      }
   }
   return vi;
}

//! Get the Pintile et al. Q Score
//!
//! @return a coot::validation_information_t object
coot::validation_information_t
molecules_container_t::get_q_score(int imol_model, int imol_map) const {


   // The Q score is stored as a UDD per atom.
   //
   // We need to extract the Q Score for each atom and average them
   // to make scores for each residue

   coot::validation_information_t vi;
   bool do_per_atom = false; // 20240629-PE we currently don't have a (good) container for this.
                             //  Per residue will do for now.

   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
         const clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
         mean_and_variance<float> mv = map_density_distribution(xmap, 1000, false, true);
         mmdb::Manager *mol = molecules[imol_model].atom_sel.mol;
         coot::q_score_t q_score(mol);
         q_score.calc(xmap, mv.mean, std::sqrt(mv.variance));
         int udd_q_score = mol->GetUDDHandle(mmdb::UDR_ATOM, "Q Score");
         vi = get_q_score_validation_information(mol, udd_q_score, do_per_atom);
         q_score.close(); // delete selection
      }
   }
   return vi;
}

//! Get the Pintile et al. Q Score for a particular residue (typically a ligand)
//!
//! @param cid If the `cid` matches more than one residue the score will be returned for all of the
//! residues covered in the `cid`. Typically, of course the `cid` will be something like
//! "//A/301".
//!
//! @return a coot::validation_information_t object
coot::validation_information_t
molecules_container_t::get_q_score_for_cid(int imol_model, const std::string &cid, int imol_map) const {

   coot::validation_information_t vi;
   bool do_per_atom = false; // 20240629-PE we currently don't have a (good) container for this.
                             //  Per residue will do for now.

   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
         const clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
         mean_and_variance<float> mv = map_density_distribution(xmap, 1000, false, true);
         mmdb::Manager *mol = molecules[imol_model].atom_sel.mol;
         std::vector<mmdb::Residue *> residues = molecules[imol_model].cid_to_residues(cid);
         coot::q_score_t q_score(mol, cid);
         q_score.calc(xmap, mv.mean, std::sqrt(mv.variance));
         int udd_q_score = mol->GetUDDHandle(mmdb::UDR_ATOM, "Q Score");
         vi = get_q_score_validation_information(mol, udd_q_score, do_per_atom);
         q_score.close(); // delete selection
      }
   }

   return vi;
}


//! Make a FSC-scaled map
//!
//! @return the molecule index of the new map
int
molecules_container_t::make_power_scaled_map(int imol_ref, int imol_map_for_scaling) {

  int idx = -1;
  if (is_valid_map_molecule(imol_ref)) {
    if (is_valid_map_molecule(imol_map_for_scaling)) {
      clipper::Xmap<float> &xmap_ref   = molecules[imol_ref].xmap;
      clipper::Xmap<float> &xmap_other = molecules[imol_map_for_scaling].xmap;
      clipper::Xmap<float> scaled = coot::util::power_scale(xmap_ref, xmap_other);
      bool is_em_map = molecules[idx].is_EM_map();
      int imol_new = molecules.size();
      std::string name = std::string("Copy of map ") + std::to_string(imol_map_for_scaling) + " scaled to " + std::to_string(imol_ref);
      molecules.push_back(coot::molecule_t(name, imol_new, scaled, is_em_map));
      idx = imol_new;
    }
  }
  return idx;
}


std::vector<int>
molecules_container_t::partition_map_by_chain(int imol_map, int imol_model) {

   std::vector<int> v;
   if (is_valid_map_molecule(imol_map)) {
      if (is_valid_model_molecule(imol_model)) {
         const clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
         mmdb::Manager *mol = molecules[imol_model].atom_sel.mol;
         std::string state_string;
         std::vector<std::pair<std::string, clipper::Xmap<float> > > maps_info =
            coot::util::partition_map_by_chain(xmap, mol, &state_string);
         if (! maps_info.empty()) {
            bool is_em_map = molecules[imol_map].is_EM_map();
            for (const auto &mi : maps_info) {
               std::string chain_id = mi.first;
               const clipper::Xmap<float> &xmap_new = mi.second;
               int imol_for_map = molecules.size();
               std::string label = "Partitioned map Chain " + chain_id;
               molecules.push_back(coot::molecule_t(label, imol_for_map, xmap_new, is_em_map));
               v.push_back(imol_for_map);
            }
         }
      }
   }
   return v;
}


//! make a masked map
//!
//! @return the index of the newly created mask. Return -1 on failure.
int
molecules_container_t::make_mask(int imol_map_ref, int imol_model, const std::string &atom_selection_cid, float radius) {

   int imol_map_new = -1;
   if (is_valid_map_molecule(imol_map_ref)) {
      if (is_valid_model_molecule(imol_model)) {
         const clipper::Xmap<float> &xmap = molecules[imol_map_ref].xmap;
         mmdb::Manager *mol = molecules[imol_model].atom_sel.mol;
         clipper::Cell cell = xmap.cell();
         clipper::Spacegroup spacegroup = xmap.spacegroup();
         clipper::Grid_sampling gs = xmap.grid_sampling();

         int selhandle = mol->NewSelection(); // d
         mol->Select(selhandle, mmdb::STYPE_ATOM, atom_selection_cid.c_str(), mmdb::SKEY_NEW);
         clipper::Xmap<float> xmap_new = coot::util::make_map_mask(spacegroup, cell, gs, mol, selhandle, radius, 1.0f);
         mol->DeleteSelection(selhandle);

         imol_map_new = molecules.size();
         std::string label = "Mask created by selection " + atom_selection_cid;
         bool is_em_map = molecules[imol_map_ref].is_EM_map();
         molecules.push_back(coot::molecule_t(label, imol_map_new, xmap_new, is_em_map));
      }
   }
   return imol_map_new;

}

//! transform a map and create a new map
//! @return the molecule index of the new map, -1 for failure
int
molecules_container_t::transform_map_using_lsq_matrix(int imol_map, lsq_results_t lsq_matrix,
                                                      float x, float y, float z, float radius) {

   auto make_rtop_from_lsq_results = [] (const lsq_results_t &lsq_mat) {
      clipper::Coord_orth t(lsq_mat.translation[0],lsq_mat.translation[1], lsq_mat.translation[2]);
      const std::vector<double> &m = lsq_mat.rotation_matrix;
      clipper::Mat33<double> rm(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]);
      return clipper::RTop_orth(rm,t);
   };

   int imol_map_new = -1;
   if (is_valid_map_molecule(imol_map)) {
      if (! lsq_matrix.empty()) {
         clipper::Coord_orth about_pt(x,y,z);
         clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
         clipper::RTop_orth rtop = make_rtop_from_lsq_results(lsq_matrix);
         clipper::Xmap<float> xmap_new = coot::util::transform_map(xmap, xmap.spacegroup(), xmap.cell(), rtop,
                                                                   about_pt, radius);
         imol_map_new = molecules.size();
         std::string name = "Transformed map from " + molecules[imol_map].get_name();
         bool is_em_map = molecules[imol_map].is_EM_map();
         molecules.push_back(coot::molecule_t(name, imol_map_new, xmap_new, is_em_map));
      }
   }
   return imol_map_new;
}

//! Scale map
//!
//! @param imol is the model molecule index
void
molecules_container_t::scale_map(int imol, float factor) {

   if (is_valid_map_molecule(imol)) {
      molecules[imol].scale_map(factor);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}


//! Get map vertices histogram
//!
//! @param imol is the map molecule index
//! @param n_bins is the number of bins - 200 is a reasonable default.
//! @param zoom_factor (reduces the range by the given factor)
//! centred around the median (typically 1.0 but usefully can vary until ~20.0).
//!
//! @return the map histogram
coot::molecule_t::histogram_info_t
molecules_container_t::get_map_vertices_histogram(int imol, int imol_map_for_sampling,
						  double position_x, double position_y, double position_z,
						  float radius, float contour_level,
						  unsigned int n_bins) {
   coot::molecule_t::histogram_info_t hi;
   if (is_valid_map_molecule(imol)) {
      if (is_valid_map_molecule(imol_map_for_sampling)) {
	 clipper::Coord_orth p(position_x, position_y, position_z);
	 const clipper::Xmap<float> &other_xmap = molecules[imol_map_for_sampling].xmap;
	 hi = molecules[imol].get_map_vertices_histogram(other_xmap, p, radius, contour_level,
							 map_is_contoured_using_thread_pool_flag, &thread_pool, n_bins);
      }
   }
   return hi;

}

