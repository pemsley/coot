/*
 * ligand/libres-tracer.hh
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
#ifndef LIBRES_TRACER_HH
#define LIBRES_TRACER_HH

#include <atomic>

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>
#include "utils/coot-fasta.hh"

class watch_res_tracer_data_t {
public:
   mmdb::Manager *working_mol;
   int imol_new;
   bool finished;
   bool update_flag;
   unsigned int update_count;
   std::atomic<bool> mol_edit_lock;
   watch_res_tracer_data_t(mmdb::Manager *working_mol_in, int imol_new_in) : working_mol(working_mol_in), imol_new(imol_new_in), mol_edit_lock(false) {
      update_flag = false;
      update_count = 0;
      finished = false;
   }
};

// working_mol is so that the graphics can show were we are in the building.
//
void res_tracer_proc(const clipper::Xmap<float> &xmap, float xmap_rmsd, const coot::fasta_multi &fam, double variation,
                     unsigned int n_top_spin_pairs, unsigned int n_top_fragments,
                     float rmsd_cut_off_for_flood, float flood_atom_mask_radius, float weight, unsigned int n_phi_psi_trials,
                     bool with_ncs, watch_res_tracer_data_t *watch_res_tracer_data_p);

void res_tracer_learn(const clipper::Xmap<float> &xmap, float weight, float xmap_rmsd, float rmsd_cut_off,
                      const coot::fasta_multi &fam, double variation,
                      unsigned int n_top_spin_pairs, float flood_atom_mask_radius, mmdb::Manager *reference_mol);

#endif // LIBRES_TRACER_HH
