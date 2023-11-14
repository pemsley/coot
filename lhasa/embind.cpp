#ifndef __EMSCRIPTEN__
// For LSP 
#define HAD_EMSCRIPTEN_ENV_DEFINED 0
#define __EMSCRIPTEN__
#else
#define HAD_EMSCRIPTEN_ENV_DEFINED 1
#endif

#include "embind.hpp"
#include "lhasa.hpp"
#include "../layla/ligand_editor_canvas.hpp"
#include "../layla/utils.hpp"

#if HAD_EMSCRIPTEN_ENV_DEFINED == 0
#undef __EMSCRIPTEN__
#endif


using namespace coot::ligand_editor_canvas;

EMSCRIPTEN_BINDINGS(lhasa) {
  function("remove_non_polar_hydrogens", &coot::layla::remove_non_polar_hydrogens);
  function("append_from_smiles", &lhasa::append_from_smiles);
  // function("rdkit_mol_from_smiles", &lhasa::rdkit_mol_from_smiles);
  // function("rdkit_mol_to_smiles", &lhasa::rdkit_mol_to_smiles);
  enum_<DisplayMode>("LhasaDisplayMode")
    .value("Standard", DisplayMode::Standard)
    .value("AtomIndices", DisplayMode::AtomIndices)
    .value("AtomNames", DisplayMode::AtomNames);
  class_<impl::Renderer>("LhasaRenderer")
    .constructor<std::string>();
  class_<impl::Renderer::DrawingCommand>("LhasaDrawingCommand");
    //.constructor()
  class_<DeleteTool>("LhasaDeleteTool")
    .constructor<>();
  class_<ChargeModifier>("LhasaChargeModifier")
    .constructor<>();
  class_<GeometryModifier>("LhasaGeometryModifier")
    .constructor<>();
  class_<FormatTool>("LhasaFormatTool")
    .constructor<>();
  class_<RemoveHydrogensTool>("LhasaRemoveHydrogensTool")
    .constructor<>();
  class_<ActiveTool>("LhasaActiveTool")
    // ActiveTool(ElementInsertion insertion) noexcept;
    // ActiveTool(BondModifier modifier) noexcept;
    // ActiveTool(TransformTool) noexcept;
    // ActiveTool(StructureInsertion insertion) noexcept;
    // ActiveTool(FlipTool) noexcept;
    // .constructor<ChargeModifier>()
    // .constructor<GeometryModifier>()
    // .constructor<FormatTool>()
    // .constructor<RemoveHydrogensTool>()
    .constructor<DeleteTool>()
    .constructor<>();
    // .smart_ptr<std::unique_ptr<ActiveTool>>("UniquePtrLhasaActiveTool");
  class_<CootLigandEditorCanvas>("LhasaCanvas")
    .constructor<>()
    .function("set_active_tool", &CootLigandEditorCanvas::set_active_tool)
    .function("append_molecule", &CootLigandEditorCanvas::append_molecule)
    .function("set_scale", &CootLigandEditorCanvas::set_scale)
    .function("get_scale", &CootLigandEditorCanvas::get_scale)
    .function("undo_edition", &CootLigandEditorCanvas::undo)
    .function("redo_edition", &CootLigandEditorCanvas::redo)
    .function("get_molecule_count", &CootLigandEditorCanvas::get_molecule_count)
    .function("set_allow_invalid_molecules", &CootLigandEditorCanvas::set_allow_invalid_molecules)
    .function("get_allow_invalid_molecules", &CootLigandEditorCanvas::get_allow_invalid_molecules)
    .function("get_display_mode", &CootLigandEditorCanvas::get_display_mode)
    .function("set_display_mode", &CootLigandEditorCanvas::set_display_mode)
    .function("get_smiles", &CootLigandEditorCanvas::get_smiles)
    .function("get_smiles_for_molecule", &CootLigandEditorCanvas::get_smiles_for_molecule)
    .function("clear_molecules", &CootLigandEditorCanvas::clear_molecules)
    .function("on_hover", &CootLigandEditorCanvas::on_hover)
    .function("on_left_click", &CootLigandEditorCanvas::on_left_click)
    .function("on_left_click_released", &CootLigandEditorCanvas::on_left_click_released)
    .function("on_right_click", &CootLigandEditorCanvas::on_right_click)
    .function("on_right_click_released", &CootLigandEditorCanvas::on_right_click_released)
    .function("render", &CootLigandEditorCanvas::render);
}