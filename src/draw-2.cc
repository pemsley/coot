#ifdef USE_PYTHON
#include <Python.h>
#endif // USE_PYTHON

#define GLM_ENABLE_EXPERIMENTAL // # for norm things
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>  // to_string()

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include "c-interface.h" // for update_go_to_atom_from_current_position()
#include "globjects.h"
#include "graphics-info.h"

#include "draw.hh"
#include "draw-2.hh"
#include "framebuffer.hh"

#include "text-rendering-utils.hh"
#include "cc-interface-scripting.hh"
#include "cylinder-with-rotation-translation.hh"

gint idle_contour_function(gpointer data);


glm::vec3
get_camera_up_direction(const glm::mat4 &mouse_quat_mat) {

   glm::vec4 z_p(0.0f, 1.0f, 0.0f, 1.0f);
   glm::vec4 r = z_p * mouse_quat_mat;
   glm::vec3 r3(r);
   return r3;
}

// static
gboolean
graphics_info_t::tick_function_is_active() {

   if (false)
      std::cout << "tick_function_is_active() " << do_tick_particles << " " << do_tick_spin << " " << do_tick_boids << " "
                << do_tick_hydrogen_bonds_mesh << " " << do_tick_happy_face_residue_markers << " "
                << do_tick_constant_draw << std::endl;

   if (do_tick_particles ||
       do_tick_spin      ||
       do_tick_rock      ||
       do_tick_boids     ||
       do_tick_constant_draw       ||
       do_tick_hydrogen_bonds_mesh ||
       do_tick_outline_for_active_residue ||
       do_tick_happy_face_residue_markers)
      return gboolean(TRUE);
   else
      return gboolean(FALSE);
}

// Put this and the above into graphics_info_t. And in it's own file.

gboolean
glarea_tick_func(GtkWidget *widget,
                 GdkFrameClock *frame_clock,
                 gpointer data) {

   graphics_info_t::tick_function_is_active();

   if (graphics_info_t::do_tick_particles) {
      if (graphics_info_t::particles.empty()) {
         graphics_info_t::do_tick_particles = false;
      } else {
         gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0])); // needed?
         // std::cout << "glarea_tick_func() calls update_particles() " << std::endl;
         graphics_info_t::particles.update_particles();
         graphics_info_t::mesh_for_particles.update_instancing_buffer_data_for_particles(graphics_info_t::particles);
      }
   }

   if (graphics_info_t::do_tick_spin) {
      float delta = 0.004 * graphics_info_t::idle_function_rotate_angle;
      // delta *= 10.0;
      glm::vec3 EulerAngles(0, delta, 0);
      glm::quat quat_delta(EulerAngles);
      glm::quat normalized_quat_delta(glm::normalize(quat_delta));
      glm::quat product = normalized_quat_delta * graphics_info_t::glm_quat;
      graphics_info_t::glm_quat = glm::normalize(product);
   }

   if (graphics_info_t::do_tick_rock) {
      std::chrono::time_point<std::chrono::high_resolution_clock> tp_now = std::chrono::high_resolution_clock::now();
      auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(tp_now - graphics_info_t::time_holder_for_rocking);
      double angle = delta.count() * 0.0007 * graphics_info_t::idle_function_rock_freq_scale_factor;
      double theta = 0.008 * graphics_info_t::idle_function_rock_amplitude_scale_factor * sin(angle);
      glm::vec3 EulerAngles(0, theta, 0);
      glm::quat quat_delta(EulerAngles);
      glm::quat normalized_quat_delta(glm::normalize(quat_delta));
      glm::quat product = normalized_quat_delta * graphics_info_t::glm_quat;
      graphics_info_t::glm_quat = glm::normalize(product);
   }

   if (graphics_info_t::do_tick_outline_for_active_residue > 0) {
      graphics_info_t::outline_for_active_residue_frame_count--;
      if (graphics_info_t::outline_for_active_residue_frame_count == 0) {
         graphics_info_t::do_tick_outline_for_active_residue = false;
      }
   }

   if (graphics_info_t::do_tick_constant_draw) {
      // don't change anything - I just want to remind you (well myself, I suppose) that it's here
   }

   if (graphics_info_t::do_tick_hydrogen_bonds_mesh) {
      // 20211210-PE  the rotation is done in instanced-object.shader now
   }

   if (graphics_info_t::do_tick_boids) {
      graphics_info_t::boids.update();
      std::vector<glm::mat4> mats(graphics_info_t::boids.size());
      for (unsigned int ii=0; ii<graphics_info_t::boids.size(); ii++)
         mats[ii] = graphics_info_t::boids[ii].make_mat();
      graphics_info_t::mesh_for_boids.update_instancing_buffer_data_standard(mats);
   }

   if (false)
      std::cout << "### in the glarea_tick_func() "
                << graphics_info_t::do_tick_happy_face_residue_markers << "  "
                << graphics_info_t::draw_count_for_happy_face_residue_markers << std::endl;

   if (graphics_info_t::do_tick_happy_face_residue_markers) {
      // this is a texture mesh, currently direct access to the draw flag.
      if (graphics_info_t::tmesh_for_happy_face_residues_markers.draw_this_mesh) {
         graphics_info_t::draw_count_for_happy_face_residue_markers += 1;
         graphics_info_t g;
         if (g.draw_count_for_happy_face_residue_markers >= g.draw_count_max_for_happy_face_residue_markers) {
            graphics_info_t::do_tick_happy_face_residue_markers = false;
            graphics_info_t::draw_count_for_happy_face_residue_markers = 0;

            graphics_info_t::tmesh_for_happy_face_residues_markers.draw_this_mesh = false;
         }

         // repeating code in setup_draw_for_happy_face_residue_markers()

         const std::vector<glm::vec3> &positions = graphics_info_t::happy_face_residue_marker_starting_positions;
         glm::vec3 up_uv = g.get_screen_y_uv();
         unsigned int draw_count = g.draw_count_for_happy_face_residue_markers;
         unsigned int draw_count_max = g.draw_count_max_for_happy_face_residue_markers;
         g.tmesh_for_happy_face_residues_markers.update_instancing_buffer_data_for_happy_faces(positions,
                                                                                               draw_count,
                                                                                               draw_count_max,
                                                                                               up_uv);
      }
   }

   gtk_widget_queue_draw(widget); // needed? 20210904-PE yeah... I  think so

   return graphics_info_t::tick_function_is_active();
}


#include <gdk/gdk.h>
#include "screendump-tga.hh"

void
on_glarea_realize(GtkGLArea *glarea) {

   auto setup_test_texture = [] () {
                                graphics_info_t g;
                                // g.texture_for_camera_facing_quad.init("some-test-label.png");
                                g.texture_for_camera_facing_quad.init("hud-label-rama.png");
                                // camera facing quad test
                                float image_apect_ratio = static_cast<float>(395)/static_cast<float>(93); // testt-label.png pixels
                                g.tmesh_for_camera_facing_quad.setup_camera_facing_quad(image_apect_ratio, 1.0, 0.0, 0.0);
                                GLenum err = glGetError(); if (err) std::cout << "realize() D err " << err << std::endl;
                                g.tmesh_for_hud_image_testing.setup_quad();
                                err = glGetError(); if (err) std::cout << "realize() D err " << err << std::endl;
                             };

   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(glarea), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   // 1 and 1 here!
   // std::cout << "in on_glarea_realize() " << w << " " << h << std::endl;

   // std::cout << "debug:: on_glarea_realize() about to make_current()" << std::endl;
   gtk_gl_area_make_current(glarea);
   GLenum err = glGetError();
   err = glGetError(); if (err) std::cout << "on_glarea_realize() A err " << err << std::endl;
   if (gtk_gl_area_get_error(glarea) != NULL) {
      std::cout << "OOPS:: on_glarea_realize() error on gtk_gl_area_make_current()" << std::endl;
      return;
   }
   if (gtk_gl_area_get_error(GTK_GL_AREA(glarea)) != NULL) {
      std::cout << "ERROR:: GLArea in an error state - goodbye " << std::endl;
      return;
   }

   gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(glarea), TRUE);

   // GLX_SAMPLE_BUFFERS_ARB
   // https://www.khronos.org/registry/OpenGL/extensions/ARB/ARB_multisample.txt
   // gdk/x11/gdkglcontext-x11.c
   // glXGetConfig(dpy, &visual_list[0], GLX_SAMPLE_BUFFERS_ARB, &gl_info[i].num_multisample);

   // glEnable(GL_MULTISAMPLE); // seems not to work at the moment. Needs work on the GTK->OpenGL interface

   const char *s1 = reinterpret_cast<const char *>(glGetString(GL_VERSION));
   const char *s2 = reinterpret_cast<const char *>(glGetString(GL_SHADING_LANGUAGE_VERSION));
   const char *s3 = reinterpret_cast<const char *>(glGetString(GL_RENDERER));
   const char *s4 = reinterpret_cast<const char *>(glGetString(GL_VENDOR));
   if (s1 && s2 && s3 && s4) {
      std::string ss1(s1);
      std::string ss2(s2);
      std::string ss3(s3);
      std::string ss4(s4);
      std::cout << "INFO:: GL Version:                  " << ss1 << std::endl;
      std::cout << "INFO:: GL Shading Language Version: " << ss2 << std::endl;
      std::cout << "INFO:: GL Renderer:                 " << ss3 << std::endl;
      std::cout << "INFO:: GL Vendor:                   " << ss4 << std::endl;
   } else {
      std::cout << "error:: on_glarea_realize() null from glGetString()" << std::endl;
   }

   graphics_info_t g;
   bool status = true; // was g.init_shaders();

   if (status == true) {
      // happy path

      g.init_buffers();

      err = glGetError();
      if (err) std::cout << "error:: on_glarea_realize() post init_shaders() err is " << err << std::endl;

      // graphics_info_t::shader_for_screen.Use(); // needed?

      err = glGetError();
      if (err) std::cout << "error:: start on_glarea_realize() err is " << err << std::endl;

      // g.init_shaders(); // here are errors

      if (graphics_info_t::use_framebuffers) {
         unsigned int index_offset = 0;
         graphics_info_t::screen_framebuffer.init(w, h, index_offset, "screen/occlusion");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post screen_framebuffer init() err is "
                                                << err << std::endl;
         index_offset = 1;
         graphics_info_t::blur_y_framebuffer.init(w, h, index_offset, "blur-y");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_y_framebuffer init() err is "
                                                << err << std::endl;
         index_offset = 2;
         graphics_info_t::blur_x_framebuffer.init(w, h, index_offset, "blur-x");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_x_framebuffer init() err is "
                                                << err << std::endl;
         index_offset = 3;
         graphics_info_t::combine_textures_using_depth_framebuffer.init(w, h, index_offset, "new-blur");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_combine framebuffer init() err is "
                                                << err << std::endl;
         index_offset = 4;
         graphics_info_t::blur_framebuffer.init(w, h, index_offset, "blur");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_framebuffer init() err is "
                                                << err << std::endl;

	 g.init_shaders();

      }

      // g.init_shaders();  - no validation failure here

      // std::cout << "DEBUG:: calling setup_hud_text for shader " << g.shader_for_hud_text.name << std::endl;
      setup_hud_text(w, h, graphics_info_t::shader_for_hud_text, false);
      // std::cout << "DEBUG:: calling setup_hud_text for shader " << g.shader_for_atom_labels.name << std::endl;
      setup_hud_text(w, h, graphics_info_t::shader_for_atom_labels, true);

      g.tmesh_for_hud_refinement_dialog_arrow = HUDTextureMesh("HUD tmesh for refinement dialog arrow");
      g.tmesh_for_hud_refinement_dialog_arrow.setup_quad();
      g.texture_for_hud_refinement_dialog_arrow             = Texture("refinement-dialog-arrrow.png", Texture::DIFFUSE);
      g.texture_for_hud_refinement_dialog_arrow_highlighted = Texture("refinement-dialog-arrrow-highlighted.png", Texture::DIFFUSE);

      g.tmesh_for_shadow_map.setup_quad();

      gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(glarea), TRUE);

      glEnable(GL_DEPTH_TEST);

      // At some stage, enable this.  Currently (I think) the winding on the atoms is the wrong
      // way around - replace with octaspheres and octahemispheres.
      // I have just now changed the winding on the "solid" map triangles - and now it looks
      // fine.
      // 
      // glEnable(GL_CULL_FACE); // if I enable this, then I get to see the back side
      // of the atoms. It's a weird look.

      // Make antialised lines - not in this
      if (false) {
         glEnable(GL_BLEND);
         glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
         glEnable(GL_LINE_SMOOTH);
      }

      g.setup_lights();

      float x_scale = 4.4;  // what are these numbers!?
      float y_scale = 1.2;
      x_scale = 1.002;
      y_scale = 1.002;
      g.tmesh_for_labels.setup_camera_facing_quad(x_scale, y_scale, 0.0, 0.0);

      g.setup_hud_geometry_bars();

      g.setup_hud_buttons();

      g.setup_rama_balls();

      g.setup_key_bindings();

      float double_rama_size = 0.8; // scaled by 0.5 in the gl-rama draw call.
      g.gl_rama_plot.setup_buffers(double_rama_size); // rama relative size, put it into graphics_info_t
                                                      // and allow it to be set in the API

      g.setup_draw_for_happy_face_residue_markers_init();

      g.setup_draw_for_bad_nbc_atom_pair_markers();

      g.setup_draw_for_anchored_atom_markers_init();

      g.lines_mesh_for_hud_lines.set_name("lines mesh for fps graph");

      if (false) { // testing how textures work
         setup_test_texture();
      }

      g.init_framebuffers();
      g.init_joey_ssao_stuff();

      err = glGetError();
      if (err) std::cout << "#### GL ERROR on_glarea_realize() --end-- with err " << err << std::endl;

      std::chrono::time_point<std::chrono::high_resolution_clock> tp_now = std::chrono::high_resolution_clock::now();
      graphics_info_t::previous_frame_time_for_per_second_counter = tp_now;

      unsigned int frame_time_history_list_max_n_elements = 500;
      std::vector<s_generic_vertex> empty_vertices(frame_time_history_list_max_n_elements + 40); // +40 for base and grid lines
      std::vector<unsigned int> empty_indices(1500, 0); // or some number
      g.lines_mesh_for_hud_lines.setup_vertices_and_indices(empty_vertices, empty_indices);

      // GdkGLContext *context = gtk_gl_area_get_context(GTK_GL_AREA(glarea));
      // gboolean legacy_flag = gdk_gl_context_is_legacy(context);
      // std::cout << "INFO:: gdk_gl_context_is_legacy() returns " << legacy_flag << std::endl;

      Material dummy_material;
      std::vector<s_generic_vertex> outline_empty_vertices(1000);
      std::vector<g_triangle> outline_empty_triangles(1000);
      g.mesh_for_outline_of_active_residue.import(outline_empty_vertices, outline_empty_triangles);
      g.mesh_for_outline_of_active_residue.setup(dummy_material);

      // Is this the place to set the window as unresizable?
      // No, it isn't. As far as I can see gtk_window_set_resizable() expands
      // the window in Y fully.  That's not what I want, of course.
      if (false) {
         GtkWidget *window = widget_from_builder("main_window");
	 gtk_window_set_resizable(GTK_WINDOW(window), FALSE);
      }

   } else {
      std::cout << "ERROR:: Shader compilation (init_shaders()) failed " << std::endl;
      exit(1);
   }

}

gboolean
on_glarea_render(GtkGLArea *glarea) {

   return graphics_info_t::render(false); // not to screendump framebuffer.

}


void
on_glarea_resize(GtkGLArea *glarea, gint width, gint height) {

   graphics_info_t g;

   // for the GL widget, not the window.
   g.graphics_x_size = width;
   g.graphics_y_size = height;

   std::cout << "INFO:: GtkGLArea widget dimensions " << width << " " << height << std::endl;

   // why do I need to do this?
   // setup_hud_text(width, height, g.shader_for_hud_text, false);
   // setup_hud_text(width, height, g.shader_for_atom_labels, true); // change the function name

   // g.setup_hud_geometry_bars(); // because they depend on the aspect ratio - but can't that be
                                   // passed as a uniform?

   // std::cout << "INFO:: Reset frame buffers " << width << "x" << height << std::endl;
   g.reset_frame_buffers(width, height);

   g.resize_framebuffers_textures_renderbuffers(width, height); // 20220131-PE added from crows merge

   g.reset_hud_buttons_size_and_position();
}



void
do_drag_pan_gtk3(GtkWidget *widget) {

   // This should be a graphics_info_t function

   GtkAllocation allocation;
   gtk_widget_get_allocation(widget, &allocation);
   int w = allocation.width;
   int h = allocation.height;

   graphics_info_t g;
   glm::mat4 mvp = g.get_molecule_mvp(); // modeglml matrix includes orientation with the quaternion

   float mouseX_1 = g.GetMouseBeginX() / (w * 0.5f) - 1.0f;
   float mouseY_1 = g.GetMouseBeginY() / (h * 0.5f) - 1.0f;
   float mouseX_2 = g.mouse_current_x  / (w * 0.5f) - 1.0f;
   float mouseY_2 = g.mouse_current_y  / (h * 0.5f) - 1.0f;

   glm::mat4 vp_inv = glm::inverse(mvp);

   glm::vec4 screenPos_1 = glm::vec4(mouseX_1, -mouseY_1, 1.0f, 1.0f);
   glm::vec4 screenPos_2 = glm::vec4(mouseX_2, -mouseY_2, 1.0f, 1.0f);
   glm::vec4 worldPos_1 = vp_inv * screenPos_1;
   glm::vec4 worldPos_2 = vp_inv * screenPos_2;

   glm::vec4 delta(worldPos_1 / worldPos_1.w - worldPos_2 / worldPos_2.w);
   glm::vec3 delta_v3(delta);

   g.add_to_rotation_centre(delta_v3);

   // g.update_maps();
   // if (graphics_info_t::glareas.size() > 0)
   // int contour_idle_token = g_idle_add(idle_contour_function, g.glareas[0]);

   g.update_things_on_move(); // 20211013-PE do I need the _and_redraw() version of this function?
}


void on_glarea_drag_begin(GtkGestureDrag *gesture,
                          double          x,
                          double          y,
                          GtkWidget      *area) {

   std::cout << "drag begin" << std::endl;
}

void on_glarea_drag_update(GtkGestureDrag *gesture,
                          double          x,
                          double          y,
                          GtkWidget      *area) {

   std::cout << "drag update" << std::endl;
}

void
on_glarea_swipe(GtkGestureSwipe *gesture,
          gdouble          velocity_x,
          gdouble          velocity_y,
          GtkWidget        *widget) {
   std::cout << "swipe" << std::endl;
}

void
my_glarea_add_signals_and_events(GtkWidget *glarea) {

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
   // 20220528-PE FIXME events
#else

   gtk_widget_add_events(glarea, GDK_SCROLL_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON_PRESS_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON_RELEASE_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON1_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON2_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON3_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_POINTER_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON1_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON2_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON3_MASK);
   gtk_widget_add_events(glarea, GDK_KEY_PRESS_MASK);

   // key presses for the glarea:

   gtk_widget_set_can_focus(glarea, TRUE);
   gtk_widget_grab_focus(glarea);

   g_signal_connect(glarea, "realize", G_CALLBACK(on_glarea_realize), NULL);
   g_signal_connect(glarea, "render",  G_CALLBACK(on_glarea_render),  NULL);
   g_signal_connect(glarea, "resize",  G_CALLBACK(on_glarea_resize),  NULL);
   g_signal_connect(glarea, "scroll-event",          G_CALLBACK(on_glarea_scroll),             NULL);
   g_signal_connect(glarea, "button-press-event",    G_CALLBACK(on_glarea_button_press),       NULL);
   g_signal_connect(glarea, "button-release-event",  G_CALLBACK(on_glarea_button_release),     NULL);
   g_signal_connect(glarea, "motion-notify-event",   G_CALLBACK(on_glarea_motion_notify),      NULL);
   g_signal_connect(glarea, "key-press-event",       G_CALLBACK(on_glarea_key_press_notify),   NULL);
   g_signal_connect(glarea, "key-release-event",     G_CALLBACK(on_glarea_key_release_notify), NULL);

   // testing
   GtkGesture *swipe = gtk_gesture_swipe_new(glarea);
   g_signal_connect(swipe, "swipe", G_CALLBACK(on_glarea_swipe), NULL);

#if 0
   // 20220415-PE new event controllers - this means turning off the motion and button press event callbacks above.
   //             Another time.

#if (GTK_MAJOR_VERSION >= 4)
   GtkGesture *drag = gtk_gesture_drag_new();
#else
   GtkGesture *drag = gtk_gesture_drag_new(glarea);
#endif

   g_signal_connect(drag, "drag-begin",  G_CALLBACK(on_glarea_drag_begin),  glarea);
   g_signal_connect(drag, "drag-update", G_CALLBACK(on_glarea_drag_update), glarea);

#endif // commented code

#endif // GTK VERSION

}
