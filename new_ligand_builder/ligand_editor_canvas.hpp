#ifndef COOT_LIGAND_EDITOR_CANVAS_HPP
#define COOT_LIGAND_EDITOR_CANVAS_HPP

#include "ligand_editor_canvas/model.hpp"
#include <cstdarg>
#include <cstdint>
#include <cstdlib>
#include <ostream>
#include <new>
#include "gtk/gtk.h"
// GObject declaration 
G_BEGIN_DECLS   

#define COOT_LIGAND_EDITOR_CANVAS_TYPE (coot_ligand_editor_canvas_get_type ())
G_DECLARE_FINAL_TYPE  (CootLigandEditorCanvas, coot_ligand_editor_canvas, COOT, COOT_LIGAND_EDITOR_CANVAS, GtkWidget)

G_END_DECLS


extern "C" {

CootLigandEditorCanvas *coot_ligand_editor_canvas_new();

} // extern "C"

#endif // COOT_LIGAND_EDITOR_CANVAS_HPP
