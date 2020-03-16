
#ifdef USE_PYTHON
#include <Python.h>
#endif

#include "graphics-info.h"

/*! \brief shiftfield B-factor refinement */
void 
graphics_info_t::shiftfield_b_factor_refinement(int imol) {

   int imol_map = Imol_Refinement_Map();
   if (is_valid_map_molecule(imol_map)) {
      molecules[imol_map].fill_fobs_sigfobs(); // caches
      const clipper::HKL_data<clipper::data32::F_sigF> &fobs_data = molecules[imol_map].get_original_fobs_sigfobs();
      const clipper::HKL_data<clipper::data32::Flag> &free_flag   = molecules[imol_map].get_original_rfree_flags();
      molecules[imol].shiftfield_b_factor_refinement(fobs_data, free_flag);
   }

}

/*! \brief shiftfield xyz refinement */
void
graphics_info_t::shiftfield_xyz_factor_refinement(int imol) {

}


