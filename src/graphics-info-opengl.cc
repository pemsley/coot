
#ifdef USE_PYTHON
#include <Python.h>
#endif

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>
#include "graphics-info.h"

void
graphics_info_t::init_shaders() {

   std::vector<std::reference_wrapper<Shader> > shaders = {shader_for_maps,
                                                           shader_for_models,
                                                           shader_for_central_cube,
                                                           shader_for_origin_cube,
                                                           shader_for_hud_text,
                                                           shader_for_atom_labels,
                                                           shader_for_moleculestotriangles,
                                                           shader_for_lines,
                                                           shader_for_screen,
                                                           shader_for_blur};
   std::string p = coot::package_data_dir();
   std::string d = coot::util::append_dir_dir(p, "shaders");
   std::cout << "INFO:: shader default dir: " << d << std::endl;
   std::vector<std::reference_wrapper<Shader> >::iterator it;
   for (it=shaders.begin(); it!=shaders.end(); it++)
      it->get().set_default_directory(d);

   shader_for_maps.init("map.shader", Shader::Entity_t::MAP);
   shader_for_map_caps.init("draw-map-cap.shader", Shader::Entity_t::MAP);
   shader_for_models.init("model.shader", Shader::Entity_t::MODEL);
   shader_for_central_cube.init("central-cube.shader", Shader::Entity_t::INFRASTRUCTURE);
   shader_for_origin_cube.init("central-cube.shader", Shader::Entity_t::INFRASTRUCTURE);
   shader_for_hud_text.init("hud-text.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_hud_geometry_bars.init("hud-bars.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_hud_geometry_labels.init("hud-labels.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_atom_labels.init("atom-label.shader", Shader::Entity_t::MODEL);
   shader_for_moleculestotriangles.init("moleculestotriangles.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);
   shader_for_lines.init("lines.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);
   shader_for_rama_balls.init("rama-balls.shader", Shader::Entity_t::MODEL);
   shader_for_particles.init("particles.shader", Shader::Entity_t::MODEL);
   shader_for_instanced_objects.init("instanced-objects.shader", Shader::Entity_t::INSTANCED_DISPLAY_OBJECT);

   // we use the above to make an image/texture in the framebuffer and use then
   // shader_for_screen to convert that framebuffer to the screen buffer.
   shader_for_screen.init("screen.shader", Shader::Entity_t::SCREEN);
   shader_for_blur.init("blur.shader", Shader::Entity_t::SCREEN);

}



glm::vec4
graphics_info_t::unproject(float x, float y, float z) {

   // z is 1 and -1 for front and back (or vice verse).

   if (! glareas[0]) return glm::vec4(0,0,0,0);

   GtkAllocation allocation;
   gtk_widget_get_allocation(glareas[0], &allocation);
   int w = allocation.width;
   int h = allocation.height;

   float mouseX = x / (w * 0.5f) - 1.0f;
   float mouseY = (h - y) / (h * 0.5f) - 1.0f;

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 vp_inv = glm::inverse(mvp);
   glm::vec4 screenPos_f = glm::vec4(mouseX, mouseY, z, 1.0f); // maybe +1
   glm::vec4 worldPos_f = vp_inv * screenPos_f;
   if (false) {
      std::cout << "unproject(" << x << "," << y << "," << z << ")   " << glm::to_string(mvp) << std::endl;
      std::cout << "unproject(" << x << "," << y << "," << z << ")   " << glm::to_string(vp_inv) << std::endl;
      std::cout << "unproject(" << x << "," << y << "," << z << ")   " << glm::to_string(screenPos_f) << std::endl;
      std::cout << "unproject(" << x << "," << y << "," << z << ")   " << glm::to_string(worldPos_f) << std::endl;
   }

   // to turn these into points in the world space, don't forget to divide by w.
   return worldPos_f;

}

// projected_coords are in clip space, so mouse position will have to be converted.
//
// static
glm::vec3
graphics_info_t::unproject_to_world_coordinates(const glm::vec3 &projected_coords) {

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 vp_inv = glm::inverse(mvp);
   glm::vec4 screenPos = glm::vec4(projected_coords, 1.0f);
   glm::vec4 c = vp_inv * screenPos;

   // std::cout << "debug:: unproject_to_world_coordinates() c " << glm::to_string(c) << std::endl;
   double oow = 1.0/c.w;
   return glm::vec3(c.x * oow, c.y * oow, c.z * oow);

}


int
graphics_info_t::blob_under_pointer_to_screen_centre() {

   // a quick slide would be better than a jump

   graphics_info_t g; // needed?
   int r = 0;
   if (use_graphics_interface_flag) {
      int imol_map = Imol_Refinement_Map();
      if (imol_map != -1) {
	 // OK we have a map to search.
	 // coot::Cartesian front = unproject(0.0);
	 // coot::Cartesian back  = unproject(1.0);
         // glm::vec4 glm_front = new_unproject(-0.3);
         // glm::vec4 glm_back  = new_unproject( 1.0);

         GtkAllocation allocation = graphics_info_t::get_glarea_allocation();
         int w = allocation.width;
         int h = allocation.height;

         glm::mat4 mvp = graphics_info_t::get_molecule_mvp(); // modeglml matrix includes orientation with the quaternion
         glm::mat4 vp_inv = glm::inverse(mvp);

         float mouseX_2 = g.mouse_current_x  / (w * 0.5f) - 1.0f;
         float mouseY_2 = g.mouse_current_y  / (h * 0.5f) - 1.0f;
         // I revered the sign here - it does the right thing now.
         glm::vec4 screenPos_1 = glm::vec4(mouseX_2, -mouseY_2, -1.0f, 1.0f);
         glm::vec4 screenPos_2 = glm::vec4(mouseX_2, -mouseY_2,  1.0f, 1.0f);
         glm::vec4 worldPos_1 = vp_inv * screenPos_1;
         glm::vec4 worldPos_2 = vp_inv * screenPos_2;

         double oowp_1 = 1.0/worldPos_1.w;
         double oowp_2 = 1.0/worldPos_2.w;
         coot::Cartesian front(worldPos_1.x * oowp_1, worldPos_1.y * oowp_1, worldPos_1.z * oowp_1);
	 coot::Cartesian  back(worldPos_2.x * oowp_2, worldPos_2.y * oowp_2, worldPos_2.z * oowp_2);

	 clipper::Coord_orth p1(front.x(), front.y(), front.z());
	 clipper::Coord_orth p2( back.x(),  back.y(),  back.z());
         if (false) {
            std::cout << "blob_under_pointer_to_screen_centre() " << glm::to_string(screenPos_1) << " "
                      << glm::to_string(screenPos_2) << std::endl;
            std::cout << "blob_under_pointer_to_screen_centre() " << front << " " << back << std::endl;
            // std::cout << "blob_under_pointer_to_screen_centre() " << p1.format() << " "
            // << p2.format() << std::endl;
         }
         coot::Cartesian rc = g.RotationCentre();

	 try {
	    clipper::Coord_orth blob =
	       molecules[imol_refinement_map].find_peak_along_line_favour_front(p1, p2);
	    coot::Cartesian cc(blob.x(), blob.y(), blob.z());
            // coot::Cartesian cc = front.mid_point(back);
            coot::Cartesian delta = rc - cc;
            // std::cout << "Delta: " << delta << std::endl;
	    g.setRotationCentre(cc);
	    for(int ii=0; ii<n_molecules(); ii++) {
	       molecules[ii].update_map();
	       molecules[ii].update_symmetry();
	    }
	    g.make_pointer_distance_objects();
	    graphics_draw();
	 }
	 catch (const std::runtime_error &mess) {
            std::cout << "debug:: given front " << front << " and back " << back << std::endl;
	    std::cout << mess.what() << std::endl;
	 }
      } else {
	 std::string s = "WARNING:: Refinement map not selected - no action";
	 std::cout << s << std::endl;
	 // add_status_bar_text(s.c_str());
	 info_dialog(s.c_str());
      }
   }
   return r;
}


void
graphics_info_t::set_clipping_front(float v) {

   if (perspective_projection_flag) {
      double l = eye_position.z;
      float screen_z_near_perspective_limit = l * 0.99;
      if (v < screen_z_near_perspective_limit)
         if (v > 2.0)
            screen_z_near_perspective = v;
   } else {
      clipping_front = v;
   }
   graphics_draw();
}


void
graphics_info_t::set_clipping_back(float v) {

   if (perspective_projection_flag) {
      double l = eye_position.z;
      float screen_z_far_perspective_limit = l * 1.01;
      if (v > screen_z_far_perspective_limit)
         if (v < 1000.0)
            screen_z_far_perspective = v;
   } else {
      clipping_back = v;
   }
   graphics_draw();
}


void
graphics_info_t::adjust_clipping(float d) {

   if (! graphics_info_t::perspective_projection_flag) {

      clipping_front = clipping_front * (1.0 + d);
      clipping_back  = clipping_back  * (1.0 + d);

   } else {

      // --- perspective ---

      double l = eye_position.z;
      double zf = screen_z_far_perspective;
      double zn = screen_z_near_perspective;

      // we should (and now do) concern ourselves with the distance to
      // the rotation centre so that the clipping planes are not
      // changed so that rotation centre is clipped.

      if (d < 0) {

         // close down (narrow)

         screen_z_near_perspective = l - (l-zn) * 0.97;
         screen_z_far_perspective  = l + (zf-l) * 0.95;

      } else {

         // expand

         screen_z_far_perspective  = l + (zf-l) * 1.05;
         screen_z_near_perspective = l - (l-zn) * 1.03;

      }

      float screen_z_near_perspective_limit = l * 0.99;
      float screen_z_far_perspective_limit  = l * 1.01;
      if (screen_z_near_perspective > screen_z_near_perspective_limit)
         screen_z_near_perspective = screen_z_near_perspective_limit;
      if (screen_z_far_perspective < screen_z_far_perspective_limit)
         screen_z_far_perspective = screen_z_far_perspective_limit;

      if (screen_z_near_perspective <    2.0) screen_z_near_perspective =    2.0;
      if (screen_z_far_perspective  > 1000.0) screen_z_far_perspective  = 1000.0;

      std::cout << "adjust_clipping(): debug l " << l
                << "    post-manip: " << screen_z_near_perspective << " "
                << screen_z_far_perspective << std::endl;
   }
}

//static
void
graphics_info_t::update_view_quaternion(int area_width, int area_height) {

   graphics_info_t g;
   float tbs = g.get_trackball_size();

   glm::quat tb_quat =
      g.trackball_to_quaternion((2.0*g.GetMouseBeginX() - area_width)/area_width,
                                (area_height - 2.0*g.GetMouseBeginY())/area_height,
                                (2.0*g.mouse_current_x - area_width)/area_width,
                                (area_height - 2.0*g.mouse_current_y)/area_height,
                                tbs);

   tb_quat = glm::conjugate(tb_quat); // hooray, no more "backwards" mouse motion
   glm::quat product = tb_quat * glm_quat;
   glm_quat = glm::normalize(product);

}


