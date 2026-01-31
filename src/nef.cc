
#include "nef.hh"
#include "graphics-info.h"

void read_nef(int imol, const std::string &file_name) {

   if (graphics_info_t::is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].read_nef(file_name);
   }

}

void show_nef_stuff(int imol) {

   if (graphics_info_t::is_valid_model_molecule(imol)) {
      std::cout << "show nef stuff " << imol << std::endl;
   }

}

