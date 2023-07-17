#ifndef COOT_PHI_PSI_HH
#define COOT_PHI_PSI_HH

// yet another phi psi container.
// 20230611-PE good grief!

// from rama_plot_with_canvas.hh

#include <mmdb2/mmdb_manager.h>
#include "clipper/core/ramachandran.h"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-rama.hh"

namespace coot {
   class phi_psis_for_model_t {

      void process(mmdb::Manager *mol, int imodel);

   public:
      int model_number;
      std::map<residue_spec_t, util::phi_psi_with_residues_t> phi_psi;
      explicit phi_psis_for_model_t(int model_number_in) {
         model_number = model_number_in;
      }
      explicit phi_psis_for_model_t(mmdb::Manager *mol, int model_number_in) {
         model_number = model_number_in;
         process(mol, model_number_in);
      }
      void add_phi_psi(const residue_spec_t &spec, const util::phi_psi_with_residues_t &phi_psi_in) {
         phi_psi[spec] = phi_psi_in;
      }
      util::phi_psi_with_residues_t operator[](const residue_spec_t &spec) {
         return phi_psi[spec];
      }
      unsigned int size() { return phi_psi.size(); }
   };

   class diff_sq_t {

      double v_;

      util::phi_psi_t pp_1;
      util::phi_psi_t pp_2;
      residue_spec_t res_1;
      residue_spec_t res_2;

   public:
      diff_sq_t() {};
      diff_sq_t(const util::phi_psi_t &pp_1_in, const util::phi_psi_t &pp_2_in,
                const residue_spec_t &r1, const residue_spec_t r2,
                double d) : pp_1(pp_1_in), pp_2(pp_2_in), res_1(r1), res_2(r2) {
         v_ = d;
      }
      util::phi_psi_t phi_psi_1() const { return pp_1; }
      util::phi_psi_t phi_psi_2() const { return pp_2; }
      double v() const { return v_;}
   };

   // and this is pulled out from the rama_plot class:
   // std::pair<bool, util::phi_psi_t> get_phi_psi(mmdb::PResidue *SelResidue);

}

#endif // COOT_PHI_PSI_HH
