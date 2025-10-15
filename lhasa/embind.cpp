
/* lhasa/embind.cpp
 * 
 * Copyright 2023 by Global Phasing Ltd.
 * Author: Jakub Smulski
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */
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
  // TODO: RDKit typedefinitions
  // function("remove_non_polar_hydrogens", &coot::layla::remove_non_polar_hydrogens);
  function("append_from_smiles", &lhasa::append_from_smiles);
  function("append_from_pickle_base64", &lhasa::append_from_pickle_base64);
  // TODO: RDKit typedefinitions
  // function("rdkit_mol_from_smiles", &lhasa::rdkit_mol_from_smiles);
  // TODO: RDKit typedefinitions
  // function("rdkit_mol_to_smiles", &lhasa::rdkit_mol_to_smiles);
  // TODO: RDKit typedefinitions
  // function("rdkit_mol_from_pickle", &lhasa::rdkit_mol_from_pickle);
  // TODO: RDKit typedefinitions
  // function("rdkit_mol_to_pickle", &lhasa::rdkit_mol_to_pickle);
  enum_<DisplayMode>("DisplayMode")
    .value("Standard", DisplayMode::Standard)
    .value("AtomIndices", DisplayMode::AtomIndices)
    .value("AtomNames", DisplayMode::AtomNames);
  register_vector<impl::Renderer::DrawingCommand>("DrawingCommandVector");
  register_vector<impl::Renderer::PathElement>("PathElementVector");
  class_<impl::Renderer::TextMeasurementCache>("TextMeasurementCache")
    .constructor<>();
  class_<impl::Renderer>("Renderer")
    .constructor<emscripten::val>()
    .constructor<emscripten::val, impl::Renderer::TextMeasurementCache&>()
    .function("get_commands", &impl::Renderer::get_commands);
  value_object<impl::Renderer::Color>("Color")
    .field("r", &impl::Renderer::Color::r)
    .field("g", &impl::Renderer::Color::g)
    .field("b", &impl::Renderer::Color::b)
    .field("a", &impl::Renderer::Color::a);
  value_object<impl::Renderer::BrushStyle>("BrushStyle")
    .field("color", &impl::Renderer::BrushStyle::color)
    .field("line_width", &impl::Renderer::BrushStyle::line_width);
  value_object<graphene_point_t>("GraphenePoint")
    .field("x", &graphene_point_t::x)
    .field("y", &graphene_point_t::y);
  value_object<impl::Renderer::Line>("Line")
    .field("start", &impl::Renderer::Line::start)
    .field("end", &impl::Renderer::Line::end);
  value_object<impl::Renderer::Arc>("Arc")
    .field("origin", &impl::Renderer::Arc::origin)
    .field("radius", &impl::Renderer::Arc::radius)
    .field("angle_one", &impl::Renderer::Arc::angle_one)
    .field("angle_two", &impl::Renderer::Arc::angle_two);
  class_<impl::Renderer::PathElement>("PathElement")
    .function("is_arc", &impl::Renderer::PathElement::is_arc)
    .function("as_arc", &impl::Renderer::PathElement::as_arc)
    .function("as_line", &impl::Renderer::PathElement::as_line)
    .function("is_line", &impl::Renderer::PathElement::is_line);
  class_<impl::Renderer::Path>("Path")
    .property("fill_color", &impl::Renderer::Path::fill_color)
    .property("has_fill", &impl::Renderer::Path::has_fill)
    .property("stroke_style", &impl::Renderer::Path::stroke_style)
    .property("has_stroke", &impl::Renderer::Path::has_stroke)
    .function("get_elements", &impl::Renderer::Path::get_elements);
  enum_<impl::Renderer::TextPositioning>("TextPositioning")
    .value("Normal", impl::Renderer::TextPositioning::Normal)
    .value("Sub", impl::Renderer::TextPositioning::Sub)
    .value("Super", impl::Renderer::TextPositioning::Super);
  class_<impl::Renderer::TextStyle>("TextStyle")
    .property("positioning", &impl::Renderer::TextStyle::positioning)
    .property("weight", &impl::Renderer::TextStyle::weight)
    .property("size", &impl::Renderer::TextStyle::size)
    .property("color", &impl::Renderer::TextStyle::color)
    .property("specifies_color", &impl::Renderer::TextStyle::specifies_color)
    .constructor();
  class_<impl::Renderer::TextSpan>("TextSpan")
    .property("style", &impl::Renderer::TextSpan::style)
    .property("specifies_style", &impl::Renderer::TextSpan::specifies_style)
    .function("is_newline", &impl::Renderer::TextSpan::is_newline)
    .function("has_subspans", &impl::Renderer::TextSpan::has_subspans) 
    .function("as_caption", select_const(&impl::Renderer::TextSpan::as_caption))
    .function("as_subspans", select_const(&impl::Renderer::TextSpan::as_subspans))
    .constructor<std::vector<impl::Renderer::TextSpan>>()
    .constructor();
  register_vector<impl::Renderer::TextSpan>("TextSpanVector");
  value_object<impl::Renderer::TextSize>("TextSize")
    .field("width", &impl::Renderer::TextSize::width)
    .field("height", &impl::Renderer::TextSize::height);
  class_<impl::Renderer::Text>("Text")
    .property("origin", &impl::Renderer::Text::origin)
    .property("style", &impl::Renderer::Text::style)
    .property("spans", &impl::Renderer::Text::spans)
    .constructor();
  class_<impl::Renderer::DrawingCommand>("DrawingCommand")
    .function("is_path", &impl::Renderer::DrawingCommand::is_path)
    .function("as_path", &impl::Renderer::DrawingCommand::as_path)
    .function("is_text", &impl::Renderer::DrawingCommand::is_text)
    .function("as_text", &impl::Renderer::DrawingCommand::as_text);
  class_<DeleteTool>("DeleteTool")
    .constructor<>();
  class_<ChargeModifier>("ChargeModifier")
    .constructor<>();
  class_<GeometryModifier>("GeometryModifier")
    .constructor<>();
  class_<FormatTool>("FormatTool")
    .constructor<>();
  class_<RemoveHydrogensTool>("RemoveHydrogensTool")
    .constructor<>();
  // I had to rename it to "LhasaElement" due to name clash in Moorhen
  enum_<ElementInsertion::Element>("LhasaElement")
    .value("C", ElementInsertion::Element::C)
    .value("N", ElementInsertion::Element::N)
    .value("O", ElementInsertion::Element::O)
    .value("S", ElementInsertion::Element::S)
    .value("P", ElementInsertion::Element::P)
    .value("H", ElementInsertion::Element::H)
    .value("F", ElementInsertion::Element::F)
    .value("Cl", ElementInsertion::Element::Cl)
    .value("Br", ElementInsertion::Element::Br)
    .value("I", ElementInsertion::Element::I);
  class_<ElementInsertion>("ElementInsertion")
    // Doesn't compile: https://github.com/emscripten-core/emscripten/issues/11274
    // .constructor([](std::string element){
    //   return ElementInsertion(element.c_str());
    // });
    .constructor<ElementInsertion::Element>();
  // Workaround for constructor overload not being available:
  function("element_insertion_from_symbol", &lhasa::element_insertion_from_symbol);
  // I had to rename it to "LhasaStructure" due to name clash in Moorhen
  enum_<StructureInsertion::Structure>("LhasaStructure")
    .value("CycloPropaneRing", StructureInsertion::Structure::CycloPropaneRing)
    .value("CycloButaneRing", StructureInsertion::Structure::CycloButaneRing)
    .value("CycloPentaneRing", StructureInsertion::Structure::CycloPentaneRing)
    .value("CycloHexaneRing", StructureInsertion::Structure::CycloHexaneRing)
    .value("BenzeneRing", StructureInsertion::Structure::BenzeneRing)
    .value("CycloHeptaneRing", StructureInsertion::Structure::CycloHeptaneRing)
    .value("CycloOctaneRing", StructureInsertion::Structure::CycloOctaneRing);
  class_<StructureInsertion>("StructureInsertion")
    .constructor<StructureInsertion::Structure>();
  enum_<BondModifier::BondModifierMode>("BondModifierMode")
    .value("Single", BondModifier::BondModifierMode::Single)
    .value("Double", BondModifier::BondModifierMode::Double)
    .value("Triple", BondModifier::BondModifierMode::Triple);
  class_<BondModifier>("BondModifier")
    .constructor<BondModifier::BondModifierMode>();
  enum_<TransformManager::Mode>("TransformMode")
    .value("Rotation", TransformManager::Mode::Rotation)
    .value("Translation", TransformManager::Mode::Translation);
  class_<TransformTool>("TransformTool")
    .constructor<TransformManager::Mode>();
  enum_<FlipMode>("FlipMode")
    .value("Horizontal", FlipMode::Horizontal)
    .value("Vertical", FlipMode::Vertical);
  class_<FlipTool>("FlipTool")
    .constructor<FlipMode>();
  class_<ActiveTool>("ActiveTool")
    // Embind doesn't currently support type-based overloads
    // so I needed a generic constructor from emscripten::val.
    // The code below does not compile.
    // I've no idea what's wrong with that
    // but that's why 'make_active_tool' exists.
    // EDIT: Likely doesn't compile because of this: https://github.com/emscripten-core/emscripten/issues/11274
    // .constructor<emscripten::val>(select_overload<ActiveTool(emscripten::val)>([](emscripten::val t) -> ActiveTool {
    //   std::string type_name = t.typeOf().as<std::string>();

    //   return ActiveTool();
    // }))
    .constructor<>();
  function("make_active_tool", &lhasa::make_active_tool);
  value_object<CootLigandEditorCanvas::SizingInfo>("SizingInfo")
    .field("requested_size", &CootLigandEditorCanvas::SizingInfo::requested_size);
  enum_<CootLigandEditorCanvas::MeasurementDirection>("MeasurementDirection")
    .value("HORIZONTAL", CootLigandEditorCanvas::MeasurementDirection::HORIZONTAL)
    .value("VERTICAL", CootLigandEditorCanvas::MeasurementDirection::VERTICAL);
  value_object<CanvasMolecule::QEDInfo>("QEDInfo")
    .field("number_of_hydrogen_bond_acceptors", &CanvasMolecule::QEDInfo::number_of_hydrogen_bond_acceptors)
    .field("number_of_hydrogen_bond_donors",&CanvasMolecule::QEDInfo:: number_of_hydrogen_bond_donors)
    .field("number_of_rotatable_bonds", &CanvasMolecule::QEDInfo::number_of_rotatable_bonds)
    .field("number_of_aromatic_rings", &CanvasMolecule::QEDInfo::number_of_aromatic_rings)
    .field("number_of_alerts", &CanvasMolecule::QEDInfo::number_of_alerts)
    .field("molecular_weight", &CanvasMolecule::QEDInfo::molecular_weight)
    .field("alogp", &CanvasMolecule::QEDInfo::alogp)
    .field("molecular_polar_surface_area", &CanvasMolecule::QEDInfo::molecular_polar_surface_area)
    .field("ads_mw", &CanvasMolecule::QEDInfo::ads_mw)
    .field("ads_alogp", &CanvasMolecule::QEDInfo::ads_alogp)
    .field("ads_hba", &CanvasMolecule::QEDInfo::ads_hba)
    .field("ads_hbd", &CanvasMolecule::QEDInfo::ads_hbd)
    .field("ads_psa", &CanvasMolecule::QEDInfo::ads_psa)
    .field("ads_rotb", &CanvasMolecule::QEDInfo::ads_rotb)
    .field("ads_arom", &CanvasMolecule::QEDInfo::ads_arom)
    .field("ads_alert", &CanvasMolecule::QEDInfo::ads_alert)
    .field("qed_score", &CanvasMolecule::QEDInfo::qed_score);
  class_<impl::WidgetCoreData>("ImplWidgetCoreData");
  register_map<unsigned int, std::string>("SmilesMap");
  // type already registered
  // register_map<unsigned int, std::string>("InchiKeyMap");
  // Without this, Emscripten errors out
  register_vector<unsigned int>("MoleculeIdVector");
  class_<CootLigandEditorCanvas, base<impl::WidgetCoreData>>("Canvas")
    .constructor<>()
    .function("set_active_tool", &CootLigandEditorCanvas::set_active_tool)
    // TODO: RDKit typedefinitions
    // .function("append_molecule", &CootLigandEditorCanvas::append_molecule)
    .function("update_molecule_from_smiles", &CootLigandEditorCanvas::update_molecule_from_smiles)
    .function("set_scale", &CootLigandEditorCanvas::set_scale)
    .function("get_scale", &CootLigandEditorCanvas::get_scale)
    .function("undo_edition", &CootLigandEditorCanvas::undo)
    .function("redo_edition", &CootLigandEditorCanvas::redo)
    .function("get_molecule_count", &CootLigandEditorCanvas::get_molecule_count)
    .function("get_idx_of_first_molecule", &CootLigandEditorCanvas::get_idx_of_first_molecule)
    .function("get_max_molecule_idx", &CootLigandEditorCanvas::get_max_molecule_idx)
    .function("set_allow_invalid_molecules", &CootLigandEditorCanvas::set_allow_invalid_molecules)
    .function("get_allow_invalid_molecules", &CootLigandEditorCanvas::get_allow_invalid_molecules)
    .function("get_display_mode", &CootLigandEditorCanvas::get_display_mode)
    .function("set_display_mode", &CootLigandEditorCanvas::set_display_mode)
    .function("get_smiles", &CootLigandEditorCanvas::get_smiles)
    .function("get_smiles_for_molecule", &CootLigandEditorCanvas::get_smiles_for_molecule)
    .function("get_inchi_keys", &CootLigandEditorCanvas::get_inchi_keys)
    .function("get_inchi_key_for_molecule", &CootLigandEditorCanvas::get_inchi_key_for_molecule)
    .function("get_pickled_molecule", &CootLigandEditorCanvas::get_pickled_molecule)
    .function("get_pickled_molecule_base64", &CootLigandEditorCanvas::get_pickled_molecule_base64)
    .function("clear_molecules", &CootLigandEditorCanvas::clear_molecules)
    .function("on_hover", &CootLigandEditorCanvas::on_hover)
    .function("on_scroll", &CootLigandEditorCanvas::on_scroll)
    .function("on_left_click", &CootLigandEditorCanvas::on_left_click)
    .function("on_left_click_released", &CootLigandEditorCanvas::on_left_click_released)
    .function("on_right_click", &CootLigandEditorCanvas::on_right_click)
    .function("on_right_click_released", &CootLigandEditorCanvas::on_right_click_released)
    .function("render", &CootLigandEditorCanvas::render)
    .function("measure", &CootLigandEditorCanvas::measure)
    .function("connect", &CootLigandEditorCanvas::connect);
}