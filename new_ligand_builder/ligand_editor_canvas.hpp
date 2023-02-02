#ifndef COOT_LIGAND_EDITOR_CANVAS_HPP
#define COOT_LIGAND_EDITOR_CANVAS_HPP

#include "ligand_editor_canvas/model.hpp"
#include <memory>
#include "gtk/gtk.h"
// GObject declaration 
G_BEGIN_DECLS   

#define COOT_LIGAND_EDITOR_CANVAS_TYPE (coot_ligand_editor_canvas_get_type ())
G_DECLARE_FINAL_TYPE  (CootLigandEditorCanvas, coot_ligand_editor_canvas, COOT, COOT_LIGAND_EDITOR_CANVAS, GtkWidget)

G_END_DECLS


extern "C" {

CootLigandEditorCanvas *coot_ligand_editor_canvas_new();


} // extern "C"

void coot_ligand_editor_set_active_tool(CootLigandEditorCanvas* self, std::unique_ptr<coot::ligand_editor_canvas::ActiveTool>&& active_tool);

#endif // COOT_LIGAND_EDITOR_CANVAS_HPP
