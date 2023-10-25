#include "embind.hpp"
#include "lhasa.hpp"
#include "../layla/ligand_editor_canvas.hpp"
#include "../layla/utils.hpp"

EMSCRIPTEN_BINDINGS(lhasa) {
  function("remove_non_polar_hydrogens", &coot::layla::remove_non_polar_hydrogens);
  // function("rdkit_mol_from_smiles", &lhasa::rdkit_mol_from_smiles);
  // function("rdkit_mol_to_smiles", &lhasa::rdkit_mol_to_smiles);
  class_<CootLigandEditorCanvas>("LhasaCanvas")
    .constructor<>()
    .function("set_active_tool", &CootLigandEditorCanvas::set_active_tool)
    .function("append_molecule", &CootLigandEditorCanvas::append_molecule)
    .function("set_scale", &CootLigandEditorCanvas::set_scale)
    .function("get_scale", &CootLigandEditorCanvas::get_scale)
    .function("undo_edition", &CootLigandEditorCanvas::undo_edition)
    .function("redo_edition", &CootLigandEditorCanvas::redo_edition)
    .function("get_molecule_count", &CootLigandEditorCanvas::get_molecule_count)
    .function("set_allow_invalid_molecules", &CootLigandEditorCanvas::set_allow_invalid_molecules)
    .function("get_allow_invalid_molecules", &CootLigandEditorCanvas::get_allow_invalid_molecules)
    .function("get_display_mode", &CootLigandEditorCanvas::get_display_mode)
    .function("set_display_mode", &CootLigandEditorCanvas::set_display_mode)
    .function("get_smiles", &CootLigandEditorCanvas::get_smiles)
    .function("get_smiles_for_molecule", &CootLigandEditorCanvas::get_smiles_for_molecule)
    .function("clear_molecules", &CootLigandEditorCanvas::clear_molecules);
}