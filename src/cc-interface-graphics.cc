
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

// maybe this function should have its own file?
//
//! \brief display the SMILES entry. This is the simple version - no dictionary
//! is generated.
void do_smiles_to_simple_3d_overlay_frame() {

   GtkWidget *frame = widget_from_builder("smiles_to_simple_3d_frame");
   if (frame)
      gtk_widget_set_visible(frame, TRUE);

}
