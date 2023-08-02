#ifndef LIGAND_BUILDER_UI_HPP
#define LIGAND_BUILDER_UI_HPP
#include <gtk/gtk.h>
#include "ligand_editor_canvas.hpp"
#include "ligand_builder_state.hpp"

namespace coot::ligand_editor {

// This doesn't need to be public
// void setup_actions(LigandBuilderState* state, GtkApplicationWindow* win, GtkBuilder* builder);

/// Setups the main window using GtkBuilder (created from 'layla.ui').
///
/// Instantiates CootLigandEditor canvas and connects its' signals.
/// Setups window actions.
/// Instantiates the global instance of LigandBuilderState
GtkApplicationWindow* setup_main_window(GtkApplication* app, GtkBuilder* builder);

}

#endif //  LIGAND_BUILDER_UI_HPP