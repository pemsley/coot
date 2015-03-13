
#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "graphics-info.h"

#include "c-interface.h" // needed for is_valid_model_molecule()


#include "globularize.hh"

void 
globularize(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].globularize();
      graphics_draw();
   }
} 
