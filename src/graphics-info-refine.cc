
#ifdef USE_PYTHON
#include <Python.h>
#endif

#include "graphics-info.h"

/*! \brief shiftfield B-factor refinement */
void 
graphics_info_t::shiftfield_b_factor_refinement(int imol) {

   int imol_map = Imol_Refinement_Map();
   if (is_valid_map_molecule(imol_map)) {
      try {
         molecules[imol_map].fill_fobs_sigfobs(); // caches
         const clipper::HKL_data<clipper::data32::F_sigF> *fobs_data = molecules[imol_map].get_original_fobs_sigfobs();
         const clipper::HKL_data<clipper::data32::Flag> *free_flag   = molecules[imol_map].get_original_rfree_flags();
         if (fobs_data && free_flag) {
            molecules[imol].shiftfield_b_factor_refinement(*fobs_data, *free_flag);
         } else {
            std::cout << "ERROR:: null pointer in function " << __FUNCTION__ << std::endl;
         }
      }
      catch (const std::runtime_error &rte) {
         std::cout << "ERROR:: " << rte.what() << std::endl;
      }
   }

}

/*! \brief shiftfield xyz refinement */
void
graphics_info_t::shiftfield_xyz_factor_refinement(int imol) {

   std::cout << "Not implemented." << std::endl;

}


//static
bool
graphics_info_t::showing_intermediate_atoms_from_refinement() {

   if (use_graphics_interface_flag) {
      if (moving_atoms_asc) {
         if (last_restraints) {
            return true;
         }
      }
   }
   return false;
}
