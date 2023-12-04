
#include "cc-interface.hh"

#include "widget-from-builder.hh"

#include "graphics-info.h"

// 20211019-PE these have moved into graphics_info_t now because I want to add a
// key-binding "Esc" to do an unfullscreen()

void fullscreen() {
   graphics_info_t::fullscreen();
}

void unfullscreen() {
   graphics_info_t::unfullscreen();

}

void set_use_trackpad(short int state) {
   graphics_info_t::using_trackpad = state;
}

// maybe this function should have its own file?
//
//! \brief display the SMILES entry. This is the simple version - no dictionary
//! is generated.
void do_smiles_to_simple_3d_overlay_frame() {

   GtkWidget *frame = widget_from_builder("smiles_to_simple_3d_frame");
   if (frame)
      gtk_widget_set_visible(frame, TRUE);

}

void show_coot_points_frame() {

   auto coot_points_frame_callback = +[] (gpointer user_data) {
      GtkWidget *frame = widget_from_builder("coot-points-frame");
      if (frame) {
         gtk_widget_set_visible(frame, FALSE);
      }
      return FALSE;
   };

   GtkWidget *frame = widget_from_builder("coot-points-frame");
   if (frame) {
      gtk_widget_set_visible(frame, TRUE);
      GSourceFunc cb = G_SOURCE_FUNC(coot_points_frame_callback);
      g_timeout_add(4000, cb, nullptr);
   }

}
