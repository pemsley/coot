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
      glm::quat product = normalized_quat_delta * graphics_info_t::view_quaternion;
      graphics_info_t::view_quaternion = glm::normalize(product);
   }

   if (graphics_info_t::do_tick_rock) {
      std::chrono::time_point<std::chrono::high_resolution_clock> tp_now = std::chrono::high_resolution_clock::now();
      auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(tp_now - graphics_info_t::time_holder_for_rocking);
      double angle = delta.count() * 0.0007 * graphics_info_t::idle_function_rock_freq_scale_factor;
      double theta = 0.008 * graphics_info_t::idle_function_rock_amplitude_scale_factor * sin(angle);
      glm::vec3 EulerAngles(0, theta, 0);
      glm::quat quat_delta(EulerAngles);
      glm::quat normalized_quat_delta(glm::normalize(quat_delta));
      glm::quat product = normalized_quat_delta * graphics_info_t::view_quaternion;
      graphics_info_t::view_quaternion = glm::normalize(product);
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


// void on_glarea_drag_begin(GtkGestureDrag *gesture,
//                           double          x,
//                           double          y,
//                           GtkWidget      *area) {

//    std::cout << "drag begin" << std::endl;
// }

// void on_glarea_drag_update(GtkGestureDrag *gesture,
//                           double          x,
//                           double          y,
//                           GtkWidget      *area) {

//    std::cout << "drag update" << std::endl;
// }

void
on_glarea_swipe(GtkGestureSwipe *gesture,
                gdouble          velocity_x,
                gdouble          velocity_y,
                GtkWidget        *widget) {
   std::cout << "swipe" << std::endl;
}

void
my_glarea_add_signals_and_events(GtkWidget *glarea) {

// #if (GTK_MAJOR_VERSION == 4)

//    // 20220528-PE FIXME events
//    std::cout << "-------------------- signal-connect realize" << std::endl;
//    std::cout << "-------------------- signal-connect render" << std::endl;
//    std::cout << "-------------------- signal-connect resize" << std::endl;
//    g_signal_connect(glarea, "realize", G_CALLBACK(on_glarea_realize), NULL);
//    g_signal_connect(glarea, "render",  G_CALLBACK(on_glarea_render),  NULL);
//    g_signal_connect(glarea, "resize",  G_CALLBACK(on_glarea_resize),  NULL);

// #else

//    gtk_widget_add_events(glarea, GDK_SCROLL_MASK);
//    gtk_widget_add_events(glarea, GDK_BUTTON_PRESS_MASK);
//    gtk_widget_add_events(glarea, GDK_BUTTON_RELEASE_MASK);
//    gtk_widget_add_events(glarea, GDK_BUTTON1_MOTION_MASK);
//    gtk_widget_add_events(glarea, GDK_BUTTON2_MOTION_MASK);
//    gtk_widget_add_events(glarea, GDK_BUTTON3_MOTION_MASK);
//    gtk_widget_add_events(glarea, GDK_POINTER_MOTION_MASK);
//    gtk_widget_add_events(glarea, GDK_BUTTON1_MASK);
//    gtk_widget_add_events(glarea, GDK_BUTTON2_MASK);
//    gtk_widget_add_events(glarea, GDK_BUTTON3_MASK);
//    gtk_widget_add_events(glarea, GDK_KEY_PRESS_MASK);

//    // key presses for the glarea:

//    gtk_widget_set_can_focus(glarea, TRUE);
//    gtk_widget_grab_focus(glarea);

//    g_signal_connect(glarea, "realize", G_CALLBACK(on_glarea_realize), NULL);
//    g_signal_connect(glarea, "render",  G_CALLBACK(on_glarea_render),  NULL);
//    g_signal_connect(glarea, "resize",  G_CALLBACK(on_glarea_resize),  NULL);
//    g_signal_connect(glarea, "scroll-event",          G_CALLBACK(on_glarea_scroll),             NULL);
//    g_signal_connect(glarea, "button-press-event",    G_CALLBACK(on_glarea_button_press),       NULL);
//    g_signal_connect(glarea, "button-release-event",  G_CALLBACK(on_glarea_button_release),     NULL);
//    g_signal_connect(glarea, "motion-notify-event",   G_CALLBACK(on_glarea_motion_notify),      NULL);
//    g_signal_connect(glarea, "key-press-event",       G_CALLBACK(on_glarea_key_press_notify),   NULL);
//    g_signal_connect(glarea, "key-release-event",     G_CALLBACK(on_glarea_key_release_notify), NULL);

//    // testing
//    GtkGesture *swipe = gtk_gesture_swipe_new(glarea);
//    g_signal_connect(swipe, "swipe", G_CALLBACK(on_glarea_swipe), NULL);

// #if 0
//    // 20220415-PE new event controllers - this means turning off the motion and button press event callbacks above.
//    //             Another time.

// #if (GTK_MAJOR_VERSION >= 4)
//    GtkGesture *drag = gtk_gesture_drag_new();
// #else
//    GtkGesture *drag = gtk_gesture_drag_new(glarea);
// #endif

//    g_signal_connect(drag, "drag-begin",  G_CALLBACK(on_glarea_drag_begin),  glarea);
//    g_signal_connect(drag, "drag-update", G_CALLBACK(on_glarea_drag_update), glarea);

// #endif // commented code

// #endif // GTK VERSION

}
