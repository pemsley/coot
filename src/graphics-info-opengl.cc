
#ifdef USE_PYTHON
#include <Python.h>
#endif


#include "graphics-info.h"

void
graphics_info_t::init_shaders() {

   std::cout << "--------------------- init_shaders() --------" << std::endl;
   graphics_info_t::shader_for_maps.init("map.shader", Shader::Entity_t::MAP);
   graphics_info_t::shader_for_models.init("model.shader", Shader::Entity_t::MODEL);
   graphics_info_t::shader_for_central_cube.init("central-cube.shader", Shader::Entity_t::INFRASTRUCTURE);
   graphics_info_t::shader_for_origin_cube.init("central-cube.shader", Shader::Entity_t::INFRASTRUCTURE);
   graphics_info_t::shader_for_hud_text.init("hud-text.shader", Shader::Entity_t::HUD_TEXT);
   // we use the above to make an image/texture in the framebuffer and use then
   // shader_for_screen to convert that framebuffer to the screen buffer.
   graphics_info_t::shader_for_screen.init("screen.shader", Shader::Entity_t::SCREEN);
   graphics_info_t::shader_for_blur.init("blur.shader", Shader::Entity_t::SCREEN);



}

