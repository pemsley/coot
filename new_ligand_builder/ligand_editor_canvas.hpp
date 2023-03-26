#ifndef COOT_LIGAND_EDITOR_CANVAS_HPP
#define COOT_LIGAND_EDITOR_CANVAS_HPP

#include "ligand_editor_canvas/core.hpp"
#include "ligand_editor_canvas/model.hpp"
#include "ligand_editor_canvas/tools.hpp"

/// GObject declaration had to be moved to "core.hpp"

extern "C" {

CootLigandEditorCanvas *coot_ligand_editor_canvas_new();


} // extern "C"


void coot_ligand_editor_set_active_tool(CootLigandEditorCanvas* self, std::unique_ptr<coot::ligand_editor_canvas::ActiveTool>&& active_tool);
void coot_ligand_editor_append_molecule(CootLigandEditorCanvas* self, std::shared_ptr<RDKit::RWMol> rdkit_mol) noexcept;
void coot_ligand_editor_set_scale(CootLigandEditorCanvas* self, float scale) noexcept;
float coot_ligand_editor_get_scale(CootLigandEditorCanvas* self) noexcept;

void coot_ligand_editor_undo_edition(CootLigandEditorCanvas* self) noexcept;
void coot_ligand_editor_redo_edition(CootLigandEditorCanvas* self) noexcept;

/// Returns the pointer to the molecule with the given index or nullptr when index is out of range.
///
/// Canvas owns the returned pointer
const RDKit::ROMol* coot_ligand_editor_get_rdkit_molecule(CootLigandEditorCanvas* self, unsigned int index) noexcept;
unsigned int coot_ligand_editor_get_molecule_count(CootLigandEditorCanvas* self) noexcept;

#endif // COOT_LIGAND_EDITOR_CANVAS_HPP
