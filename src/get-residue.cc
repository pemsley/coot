
#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "graphics-info.h"
#include "get-residue.hh"

mmdb::Residue *get_residue(int molecule_number, const coot::residue_spec_t &res_spec) {

   mmdb::Residue *r = NULL;
   graphics_info_t g;
   if (g.is_valid_model_molecule(molecule_number)) {
      r = g.molecules[molecule_number].get_residue(res_spec);
   }
   return r;
} 
