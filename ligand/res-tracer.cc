/*
 * ligand/res-tracer.cc
 *
 * Copyright 2021 by Medical Research Council
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

#include <clipper/ccp4/ccp4_map_io.h>

#include "utils/ctpl.h"
#include "utils/coot-fasta.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "rama-rsr-extend-fragments.hh"
#include "ligand.hh"
#include "libres-tracer.hh"

void
test_best_fit_phi_psi_func() {

   std::string hklin_file_name = "coot-download/1gwd_map.mtz";
   std::string    f_col_label   = "FWT";
   std::string phi_col_label = "PHWT";
   std::string pir_file_name = "1gwd.pir";
   std::string pdb_files = "1gwd-A-no-sidechains.pdb";

   clipper::Xmap<float> xmap;
   coot::util::map_fill_from_mtz(&xmap, hklin_file_name, f_col_label, phi_col_label, "", 0, 0);

   mmdb::Manager *mol = new mmdb::Manager;
   mmdb::ERROR_CODE err = mol->ReadCoorFile("1gwd-A-no-sidechains.pdb");
   if (err) {
      std::cout << "ERROR:: reading pdb file " << std::endl;
   } else {
      coot::residue_spec_t spec("A", 5, "");
      mmdb::Residue *residue_sel_1_p = coot::util::get_residue(spec, mol);
      mmdb::Residue *residue_sel_2_p = coot::util::get_residue(spec.next(), mol);

      if (residue_sel_1_p && residue_sel_2_p) {

         std::vector<mmdb::Residue *> residues = {residue_sel_1_p, residue_sel_2_p };
         auto mol_for_residue = coot::util::create_mmdbmanager_from_residue_vector(residues, mol);
         if (mol_for_residue.first) {
            unsigned int n_threads = coot::get_max_number_of_threads();
            ctpl::thread_pool thread_pool(n_threads);
            float weight = 60.0;
            unsigned int n_phi_psi_trials = 1000;
            coot::protein_geometry geom;
            geom.init_standard();
            unsigned int update_count = 0;
            std::pair<float, float> mv = coot::util::mean_and_variance(xmap);
            float xmap_rmsd = std::sqrt(mv.second);
            rama_rsr_extend_fragments(mol_for_residue.second, xmap, xmap_rmsd, &thread_pool, n_threads, weight, n_phi_psi_trials, geom, &update_count);
            mol_for_residue.second->WritePDBASCII("test_best_fit_phi_psi_result.pdb");
         } else {
            std::cout << "No mol_for_residue" << std::endl;
         }
      } else {
         std::cout << "No mol for residue " << spec << std::endl;
      }
   }
}

#include "coot-utils/atom-selection-container.hh"

void
learn_lyso() {

   std::string hklin_file_name = "1gwd_map.mtz";
   std::string f_col_label   = "FWT";
   std::string phi_col_label = "PHWT";
   std::string pir_file_name = "1gwd.pir";
   std::string pdb_file = "1gwd.cif";
   try {
      coot::fasta_multi fam;
      fam.read(pir_file_name);
      clipper::Xmap<float> xmap;
      std::cout << "Read mtz file " << hklin_file_name << " " << f_col_label << " " << phi_col_label << std::endl;
      bool use_weights = false;
      coot::util::map_fill_from_mtz(&xmap, hklin_file_name, f_col_label, phi_col_label, "", use_weights);
      atom_selection_container_t asc = get_atom_selection(pdb_file, false);
      double variation = 0.4;
      unsigned int n_top_spin_pairs = 70000; // deep into the mire (scores are for backwards pairs also)
      float flood_atom_mask_radius = 0.6; // was 1.0
      float weight = 8.0f; // calculate this (using rmsd)
      mmdb::Manager *working_mol = new mmdb::Manager;
      float rmsd_cuffoff = 2.3;
      std::pair<float, float> mv = coot::util::mean_and_variance(xmap);
      float xmap_rmsd = std::sqrt(mv.second);
      int imol = 0;
      watch_res_tracer_data_t tracer_data(working_mol, imol);
      res_tracer_learn(xmap, weight, xmap_rmsd, rmsd_cuffoff, fam, variation, n_top_spin_pairs, flood_atom_mask_radius, asc.mol);

   }
   catch (clipper::Message_fatal &rte) {
      std::cout << "ERROR::" << rte.text() << std::endl;
   }
}

void
learn() {

   learn_lyso();
}

int main(int argc, char **argv) {

   int status = 0;
   bool test_best_fit_phi_psi = false;
   bool with_ncs = false;

   std::string hklin_file_name = "rnasa-1.8-all_refmac1.mtz";
   std::string f_col_label = "FWT";
   std::string phi_col_label = "PHWT";
   std::string pir_file_name = "../src/rnase.pir";

   std::string map_file_name = "1gwd_map_mainchain-masked.map";

   if (false) {
      hklin_file_name = "../src/tm-A-1.0.mtz";
      f_col_label   = "FC";
      phi_col_label = "PHIC";
   }

   if (true) {
      hklin_file_name = "1gwd_map.mtz";
      f_col_label   = "FWT";
      phi_col_label = "PHWT";
      pir_file_name = "1gwd.pir";
   }

   if (false) { // this is hard
      hklin_file_name = "65_bucc3atest3_fphiout_acorn.mtz";
      f_col_label   = "F";
      phi_col_label = "PHI";
      // pir_file_name = "all-e-coli.fasta";
      pir_file_name = "1gwd.pir";
   }

   if (false) { // DPS map
      map_file_name = "emd_11210.map";
      pir_file_name = "6zgl.fasta";
   }

   bool do_learn = false;
   if (argc > 1) {
      std::string argv_1 = argv[1];
      if (argv_1 == "test-best-fit-phi-psi") test_best_fit_phi_psi = true;
      if (argv_1 == "learn") do_learn = true;
   }

   if (test_best_fit_phi_psi) {

      test_best_fit_phi_psi_func();

   } else {

      if (do_learn) {

         learn();

      } else {

         coot::fasta_multi fam;
         fam.read(pir_file_name);

         bool test_from_crystallography = true;
         bool test_from_map = false;

         // ===================================================================
         //                user input map and fam file
         // ===================================================================
         if (argc == 3) {
            map_file_name = argv[1];
            pir_file_name = argv[2];
            test_from_crystallography = false;
            test_from_map = true;
         }

         double variation = 0.8; // how much is a CA-CA distance allowed to vary from 3.81
         // before it is too long or short to be considered a potential
         // peptide?
         // make bigger at lower resolutions (maybe up to 1.0?)
         variation = 0.4; // speed
         // variation = 0.8; // 20240811-PE 1.4 crashes at the close of grow_trees_v4() atm, 1.0 is OK, 1.2 crashes, 1.1 crashes

         unsigned int n_top_spin_pairs = 500; // Use for tracing at most this many spin score pairs (which have been sorted).
         // This and variation affect the run-time (and results?)
         n_top_spin_pairs = 200000; // was 1000
         n_top_spin_pairs = 1000;

         unsigned int n_top_fragments = 1000; // was 4000 // The top 1000 fragments at least are all the same trace for no-side-chain lyso test
         n_top_fragments = 10000;
         n_top_fragments = 1000;

         // pass the atom radius from here too

         float flood_atom_mask_radius = 1.0; // was 0.6 for emdb

         unsigned int n_phi_psi_trials = 40000; // was 5000 // was 40000 // 40000 seems good for production

         float weight = 8.0f; // calculate this (using rmsd)

         if (false) {     // EMD-22898
            // map_file_name = "emd_22898_blur_20.map";
            weight = 6.0;
            flood_atom_mask_radius = 0.7;
            n_top_fragments = 2000;
            n_top_spin_pairs = 30000; // was 10000
            variation = 0.8;
            n_phi_psi_trials = 10000;// was 50000
         }

         // for RNAse test
         if (false) {
            n_phi_psi_trials = 10000;
            n_top_spin_pairs =  3000;
            n_top_fragments  =  5000;
            flood_atom_mask_radius = 0.8;
            variation = 0.3;
            with_ncs = true;
         }

         // for acorn
         if (false) {
            n_phi_psi_trials = 10000;
            n_top_spin_pairs =  3000;
            n_top_fragments  =  10000;
            flood_atom_mask_radius = 0.8;
            variation = 0.4;
            with_ncs = true;
         }

         if (test_from_crystallography) {

            clipper::Xmap<float> xmap;
            if (test_from_map) {
               if (coot::file_exists(map_file_name)) {
                  clipper::CCP4MAPfile file;
                  file.open_read(map_file_name);
                  file.import_xmap(xmap);
                  std::cout << "Imported map from file " << map_file_name << std::endl;
               } else {
                  std::cout << "Map file does not exist " << map_file_name << std::endl;
               }
            } else {

               if (! f_col_label.empty()) {
                  if (! phi_col_label.empty()) {
                     try {
                        std::cout << "Read mtz file " << hklin_file_name << " " << f_col_label << " " << phi_col_label << std::endl;
                        bool use_weights = false;
                        coot::util::map_fill_from_mtz(&xmap, hklin_file_name, f_col_label, phi_col_label, "", use_weights);
                     }
                     catch (clipper::Message_fatal &rte) {
                        std::cout << "ERROR::" << rte.text() << std::endl;
                     }
                  }
               }
            }
            if (! xmap.is_null()) {

               mmdb::Manager *working_mol = new mmdb::Manager;
               unsigned int update_count = 0;

               float rmsd_cuffoff = 2.3;
               std::pair<float, float> mv = coot::util::mean_and_variance(xmap);
               float xmap_rmsd = std::sqrt(mv.second);

               int imol = 0;
               watch_res_tracer_data_t tracer_data(working_mol, imol);
               res_tracer_proc(xmap, xmap_rmsd, fam, variation, n_top_spin_pairs, n_top_fragments, rmsd_cuffoff, flood_atom_mask_radius,
                               weight, n_phi_psi_trials, with_ncs, &tracer_data);

            } else {
               if (test_from_map)
                  std::cout << "Map from " << map_file_name << " is null" << std::endl;
            }

         } else {  // not crystallography

            try {

               if (coot::file_exists(map_file_name)) {

                  std::cout << "::::::::: cryo-em map path" << map_file_name << std::endl;
                  clipper::CCP4MAPfile file;
                  file.open_read(map_file_name);
                  clipper::Xmap<float> xmap;
                  fam.read(pir_file_name);
                  file.import_xmap(xmap);
                  mmdb::Manager *working_mol = new mmdb::Manager;
                  int imol = 0;
                  flood_atom_mask_radius = 0.6;
                  watch_res_tracer_data_t tracer_data(working_mol, imol);
                  std::pair<float, float> mv = coot::util::mean_and_variance(xmap);
                  float xmap_rmsd = std::sqrt(mv.second);

                  bool do_reinterp = true;
                  float reinterp_factor = 2.0; // was 3.0;

                  if (do_reinterp) {
                     // 1gwd has ~0.3A/grid - is that a factor?
                     // 20240809-PE doesn't seem to make a difference.
                     std::cout << "reinterpreting map " << reinterp_factor << std::endl;
                     clipper::Xmap<float> xmap_fine = coot::util::reinterp_map(xmap, reinterp_factor);
                     xmap = xmap_fine; // replace it
                  }

                  res_tracer_proc(xmap, xmap_rmsd, fam, variation, n_top_spin_pairs, n_top_fragments, 4.5, flood_atom_mask_radius, weight, n_phi_psi_trials,
                                  with_ncs, &tracer_data); // working_mol, &update_count
               } else {
                  std::cout << "No such file " << map_file_name << std::endl;
               }
            }
            catch (clipper::Message_fatal &rte) {
               std::cout << "ERROR::" << rte.text() << std::endl;
            }
         }
      }
   }

   return status;
}
//
