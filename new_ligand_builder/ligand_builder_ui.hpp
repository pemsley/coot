#ifndef LIGAND_BUILDER_UI_HPP
#define LIGAND_BUILDER_UI_HPP
#include <gtk/gtk.h>
#include "ligand_editor_canvas.hpp"
#include "ligand_builder_state.hpp"

namespace coot::ligand_editor {

void setup_actions(LigandBuilderState* state, GtkApplicationWindow* win, GtkBuilder* builder);

}

#endif //  LIGAND_BUILDER_UI_HPP