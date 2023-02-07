#ifndef COOT_LIGAND_EDITOR_CANVAS_HPP
#define COOT_LIGAND_EDITOR_CANVAS_HPP

#include "ligand_editor_canvas/model.hpp"
#include "ligand_editor_canvas/tools.hpp"
#include <rdkit/GraphMol/RWMol.h>
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
void coot_ligand_editor_append_molecule(CootLigandEditorCanvas* self, std::shared_ptr<RDKit::RWMol> rdkit_mol) noexcept;

#endif // COOT_LIGAND_EDITOR_CANVAS_HPP
