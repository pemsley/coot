#ifndef LIGAND_BUILDER_UI_HPP
#define LIGAND_BUILDER_UI_HPP
#include <gtk/gtk.h>
#include "ligand_editor_canvas.hpp"

namespace coot::ligand_editor {

void build_main_window(GtkWindow* win, CootLigandEditorCanvas* canvas, GtkLabel* status_label);
GMenu *build_menu(GtkApplication* app, CootLigandEditorCanvas* canvas, GtkWindow* win);

}

#endif //  LIGAND_BUILDER_UI_HPP