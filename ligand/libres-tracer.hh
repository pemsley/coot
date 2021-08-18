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

// working_mol is so that the graphics can show were we are in the building. Change count_p when working_mol gets updated.
// If count_p is null, then don't try to update working_mol
//
// When the work is finished, set *finished_p to true;
void res_tracer_proc(const clipper::Xmap<float> &xmap, const coot::fasta_multi &fam, double variation,
                     unsigned int n_top_spin_pairs, unsigned int n_top_fragments,
                     float rmsd_cut_off_for_flood, float flood_atom_mask_radius, float weight, unsigned int n_phi_psi_trials,
                     bool with_ncs, watch_res_tracer_data_t *watch_res_tracer_data_p);

#endif // LIBRES_TRACER_HH
