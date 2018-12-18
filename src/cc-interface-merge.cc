
// --------------------------------------------
#ifdef USE_PYTHON
#include <Python.h> // add first for _XOPEN_SOURCE order issues
#endif

#include <cstddef> // define std::ptrdiff_t

#ifdef USE_GUILE
#include <libguile.h>
#endif

#include "graphics-info.h"

// --------------------------------------------
#include "c-interface.h" // for is_valid_model_molecule()
#include "cc-interface.hh"

/*  ------------------------------------------------------------------------ */
/*                         merge fragments                                   */
/*  ------------------------------------------------------------------------ */
//! \name Merge Fragments
//! \{
//! \brief merge fragments
// 
//! each fragment is presumed to be in its own chain.
//
int merge_fragments(int imol) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = graphics_info_t::molecules[imol].merge_fragments();
      graphics_draw();
   }

   return status;

}
//! \}
