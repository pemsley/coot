
#include "M2T-interface.hh"
#include "graphics-info.h"

//! Update float parameter for MoleculesToTriangles molecular mesh
void M2T_updateFloatParameter(int imol, const std::string &param_name, float value) {

   if (graphics_info_t::is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].M2T_updateFloatParameter(param_name, value);
   }
   graphics_info_t::graphics_draw();
}

//! Update int parameter for MoleculesToTriangles molecular mesh
void M2T_updateIntParameter(int imol, const std::string &param_name, int value) {

   if (graphics_info_t::is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].M2T_updateIntParameter(param_name, value);
   }
   graphics_info_t::graphics_draw();
}
