
#include "graphics-info.h"

// 20250725-PE these functions were remove from  graphics-info.h


//! \brief enable/disable the model animation (on or off)
void
graphics_info_t::set_model_animation_state(unsigned int model_index, bool state) {
   if (model_index < models.size()) {
      auto &model = models[model_index];
      model.set_do_animation(state);
   }
}

// static
std::string
graphics_info_t::stringify_error_code(GLenum err) {

      std::string r = std::to_string(err);
      if (err == GL_INVALID_ENUM)      r = "GL_INVALID_ENUM";
      if (err == GL_INVALID_VALUE)     r = "GL_INVALID_VALUE";
      if (err == GL_INVALID_OPERATION) r = "GL_INVALID_OPERATION";
      return r;
};

// static
void
graphics_info_t::attach_buffers(const char *calling_function) {

   // During high-res screendump, redirect to the screendump FBO instead of the GTK default
   if (screendump_target_framebuffer != 0) {
      glBindFramebuffer(GL_FRAMEBUFFER, screendump_target_framebuffer);
      return;
   }

   bool print_errors = false;
   if (use_graphics_interface_flag) {
      if (print_errors) {
         GLenum err = glGetError();
         if (err) {
            std::cout << "GL ERROR:: attach_buffers --- start --- "
                      << stringify_error_code(err) <<  " \n";
#ifdef USE_BACKWARD
            backward::StackTrace st;
            backward::Printer p;
            st.load_here(32);
            p.print(st);
#endif
         }
         auto gl_area = glareas[0];
         gtk_gl_area_attach_buffers(GTK_GL_AREA(gl_area));
         err = glGetError();
         if (err) {
            if (calling_function)
               std::cout << "GL ERROR:: attach_buffers() --- post gtk_gl_area_attach_buffers() "
                         << stringify_error_code(err) << " with gl_area " << gl_area
                         << " calling function: " << calling_function << "()\n";
            else
               std::cout << "GL ERROR:: attach_buffers() --- post gtk_gl_area_attach_buffers() "
                         << stringify_error_code(err) << " with gl_area " << gl_area << "\n";
#ifdef USE_BACKWARD
            backward::StackTrace st;
            backward::Printer p;
            st.load_here(32);
            p.print(st);
#endif
         }
      } else {

         // cleaner output
         auto gl_area = glareas[0];
         gtk_gl_area_attach_buffers(GTK_GL_AREA(gl_area));

      }
   }
}


// static
void
graphics_info_t::graphics_grab_focus() {

      if (use_graphics_interface_flag) {
         if (! glareas.empty()) {
            GtkWidget *glarea = glareas[0];
            gtk_widget_grab_focus(glarea);
         }
      }
   }

//static
void
graphics_info_t::graphics_draw() {

      // Don't put timing things here - it's not called when tick function is used (somehow). Put it in render()
      if (use_graphics_interface_flag) {
         if (! glareas.empty()) {
            for (unsigned int i=0; i<glareas.size(); i++) {
               GtkWidget *glarea = glareas[i];
               gtk_widget_queue_draw(glarea);
               if (make_movie_flag)
                  dump_a_movie_image();
            }
         }
      }
      if (! smooth_scroll_on_going) // 20230423-PE exclude other animations too?
         draw_rama_plots(); // the widgeted rama plots, not the in-window one.
   }

// static
gl_context_info_t
graphics_info_t::get_gl_context_info() {

   gl_context_info_t glc; // null default
   if (glareas.size() > 0) glc.widget_1 = glareas[0];
   if (glareas.size() > 1) glc.widget_2 = glareas[1];
   return glc;
}

// static
void
graphics_info_t::add_a_tick() {
   // needs glarea-tick-func.hh
   if (! tick_function_is_active())
      tick_function_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
}

void
graphics_info_t::set_shadow_texture_resolution_multiplier(unsigned int m) {
   if (m != 0) {
      if (m < 8) {
         if (shadow_texture_multiplier != m) {
            shadow_texture_multiplier = m;
            shadow_texture_width  = 1024 * m;
            shadow_texture_height = 1024 * m;
            // rengerate the framebuffer texture
            glBindTexture(GL_TEXTURE_2D, shadow_depthMap_texture);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, shadow_texture_width, shadow_texture_height,
                         0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
         }
      }
   }
}


// static
void
graphics_info_t::add_model(const Model &model) {
   models.push_back(model);
}

void
graphics_info_t::scale_model(unsigned int idx, float scale_factor) {
   attach_buffers(); // because the glbufferdata is changed
   if (idx < models.size())
      models[idx].scale(scale_factor);
}




// new file for these generic object functions?


// static
void
graphics_info_t::on_generic_objects_dialog_object_check_button_toggled(GtkButton       *button,
                                                                       gpointer         user_data) {

   // std::cout << "in on_generic_objects_dialog_object_check_button_toggled() " << std::endl;
   int generic_object_number = GPOINTER_TO_INT(user_data);
   int state = 0;
   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(button)))
      state = 1;
   set_display_generic_object_simple(generic_object_number, state);
   graphics_draw();
}


void
graphics_info_t::generic_objects_dialog_grid_add_object_internal(const meshed_generic_display_object &gdo,
                                                                 GtkWidget *dialog,
                                                                 GtkWidget *grid,
                                                                 int io) {

   // std::cout << "generic_objects_dialog_grid_add_object_internal() --- start --- " << std::endl;

   if (! gdo.mesh.is_closed()) {
      // std::cout << "generic_objects_dialog_grid_add_object_internal() no-closed " << io << std::endl;
      GtkWidget *checkbutton = gtk_check_button_new_with_mnemonic (("Display"));
      std::string label_str = gdo.mesh.name;
      GtkWidget *label = gtk_label_new(label_str.c_str());

      std::string stub = "generic_object_" + std::to_string(io);
      std::string toggle_button_name = stub + "_toggle_button";
      std::string label_name = stub + "_label";

      // set the names of these widgets so that they can be
      // looked up and toggled/hidden dynamically.

      if (dialog) {
         g_object_set_data(G_OBJECT(dialog), toggle_button_name.c_str(), checkbutton);
         g_object_set_data(G_OBJECT(dialog), label_name.c_str(), label);
      } else {
         std::cout << "WARNING:: null dialog in generic_objects_dialog_grid_add_object_internal()" << std::endl;
      }

      // grid child left top width height
      gtk_grid_attach (GTK_GRID (grid), label,       0, io, 1, 1);
      gtk_grid_attach (GTK_GRID (grid), checkbutton, 1, io, 1, 1);

      if (gdo.mesh.get_draw_this_mesh())
         gtk_check_button_set_active(GTK_CHECK_BUTTON(checkbutton), TRUE);

      g_signal_connect(G_OBJECT(checkbutton), "toggled",
                       G_CALLBACK(on_generic_objects_dialog_object_check_button_toggled),
                       GINT_TO_POINTER(io));

      gtk_widget_set_visible (label, TRUE);
      gtk_widget_set_visible (checkbutton, TRUE);

   }

}

int
graphics_info_t::new_generic_object_number(const std::string &name) {
   Mesh mesh(name);
   meshed_generic_display_object meshed(mesh);
   generic_display_objects.push_back(meshed);
   int n_new = generic_display_objects.size() - 1;
   if (use_graphics_interface_flag) {
      GtkWidget *grid = widget_from_builder("generic_objects_dialog_grid"); // changed 20211020-PE
      if (grid) {
         // 20240420-PE the class variable generic_objects_dialog needs to be removed
         GtkWidget *generic_objects_dialog = widget_from_builder("generic_objects_dialog");
         const meshed_generic_display_object &gdo = generic_display_objects[n_new];
         generic_objects_dialog_grid_add_object_internal(gdo,
                                                         generic_objects_dialog,
                                                         grid,
                                                         n_new);
      }
   }
   return n_new;
}

int
graphics_info_t::new_generic_object_number_for_molecule(const std::string &name, int imol) {
   int idx = new_generic_object_number(name);
   generic_display_objects.at(idx).imol = imol;
   return idx;
}

// static
int
graphics_info_t::generic_object_index(const std::string &name) {
   int index = -1;
   int nobjs = generic_display_objects.size();
   for (int iobj=0; iobj<nobjs; iobj++) {
      if (generic_display_objects[iobj].mesh.name == name) {
         if (!generic_display_objects[iobj].mesh.this_mesh_is_closed) {
            index = iobj;
            break;
         }
      }
   }
   return index;
}

// static
void
graphics_info_t::set_mouse_previous_position(double x, double y) {
   mouse_previous_position.first = x; mouse_previous_position.second = y; }

// static
double
graphics_info_t::get_mouse_previous_position_x() { return mouse_previous_position.first; }

// static
double
graphics_info_t::get_mouse_previous_position_y() { return mouse_previous_position.second; }
