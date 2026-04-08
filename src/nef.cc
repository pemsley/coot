
#include "nef.hh"
#include "graphics-info.h"

void read_nef(int imol, const std::string &file_name) {


#ifdef USE_GEMMI
   if (graphics_info_t::is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].read_nef(file_name);
   }
#else
   std::cout << "WARNING:: this build not compiled with GEMMI" << std::endl;
#endif // USE_GEMMI

}

void show_nef_stuff(int imol) {

   if (graphics_info_t::is_valid_model_molecule(imol)) {
      std::cout << "show nef stuff " << imol << std::endl;
   }

}

