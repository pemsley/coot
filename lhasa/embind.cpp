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
  register_vector<impl::Renderer::DrawingCommand>("LhasaDrawingCommandVector");
  class_<impl::Renderer>("LhasaRenderer")
    .constructor<std::string>()
    .function("get_commands", &impl::Renderer::get_commands);
  value_object<impl::Renderer::BrushStyle>("LhasaBrushStyle")
    .field("r", &impl::Renderer::BrushStyle::r)
    .field("g", &impl::Renderer::BrushStyle::g)
    .field("b", &impl::Renderer::BrushStyle::b)
    .field("a", &impl::Renderer::BrushStyle::a)
    .field("line_width", &impl::Renderer::BrushStyle::line_width);
  value_object<graphene_point_t>("GraphenePoint")
    .field("x", &graphene_point_t::x)
    .field("y", &graphene_point_t::y);
  value_object<impl::Renderer::Line>("LhasaLine")
    .field("start", &impl::Renderer::Line::start)
    .field("end", &impl::Renderer::Line::end)
    .field("style", &impl::Renderer::Line::style);
  value_object<impl::Renderer::Arc>("LhasaArc")
    .field("origin", &impl::Renderer::Arc::origin)
    .field("radius", &impl::Renderer::Arc::radius)
    .field("angle_one", &impl::Renderer::Arc::angle_one)
    // todo: fill
    .field("angle_two", &impl::Renderer::Arc::angle_two);
  value_object<impl::Renderer::Path>("LhasaPath")
    // todo: fill
    .field("commands", &impl::Renderer::Path::commands);
  class_<impl::Renderer::TextSpan>("LhasaTextSpan")
    // todo: everything
    .property("caption", &impl::Renderer::TextSpan::caption);    
  register_vector<impl::Renderer::TextSpan>("LhasaTextSpanVector");
  class_<impl::Renderer::Text>("LhasaText")
    // todo: everything
    .property("spans", &impl::Renderer::Text::spans);
  class_<impl::Renderer::DrawingCommand>("LhasaDrawingCommand")
    .function("is_path", &impl::Renderer::DrawingCommand::is_path)
    .function("is_arc", &impl::Renderer::DrawingCommand::is_arc)
    .function("is_line", &impl::Renderer::DrawingCommand::is_line)
    .function("is_text", &impl::Renderer::DrawingCommand::is_text)
    .function("as_path", &impl::Renderer::DrawingCommand::as_path)
    .function("as_arc", &impl::Renderer::DrawingCommand::as_arc)
    .function("as_line", &impl::Renderer::DrawingCommand::as_line)
    .function("as_text", &impl::Renderer::DrawingCommand::as_text);
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
  value_object<CootLigandEditorCanvas::SizingInfo>("LhasaSizingInfo")
    .field("requested_size", &CootLigandEditorCanvas::SizingInfo::requested_size);
  enum_<CootLigandEditorCanvas::MeasurementDirection>("LhasaMeasurementDirection")
    .value("HORIZONTAL", CootLigandEditorCanvas::MeasurementDirection::HORIZONTAL)
    .value("VERTICAL", CootLigandEditorCanvas::MeasurementDirection::VERTICAL);
  class_<impl::WidgetCoreData>("LhasaImplWidgetCoreData");
  class_<CootLigandEditorCanvas, base<impl::WidgetCoreData>>("LhasaCanvas")
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
    .function("render", &CootLigandEditorCanvas::render)
    .function("measure", &CootLigandEditorCanvas::measure);
}