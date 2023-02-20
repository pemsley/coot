#include "tools.hpp"
#include "core.hpp"
#include "model.hpp"
#include <stdexcept>
#include <variant>

using namespace coot::ligand_editor_canvas;

ActiveTool::ActiveTool() noexcept {
    this->variant = ActiveTool::Variant::None;
}

ActiveTool::ActiveTool(ElementInsertion insertion) noexcept {
    this->variant = ActiveTool::Variant::ElementInsertion;
    this->element_insertion = insertion;
}

ActiveTool::ActiveTool(BondModifier modifier) noexcept {
    this->variant = ActiveTool::Variant::BondModifier;
    this->bond_modifier = modifier;
}

ActiveTool::Variant ActiveTool::get_variant() const noexcept {
    return this->variant;
}

void ActiveTool::set_core_widget_data(impl::CootLigandEditorCanvasPriv* owning_widget) noexcept {
    this->widget_data = static_cast<impl::WidgetCoreData*>(owning_widget);
}

void ActiveTool::check_variant(ActiveTool::Variant expected) {
    if (expected != this->variant) {
        throw std::runtime_error("Unexpected ActiveTool variant.");
    }
}



void ActiveTool::insert_atom(int x, int y) {
    check_variant(Variant::ElementInsertion);
    auto& element_insertion = this->element_insertion;
    const char* el_name = element_insertion.get_element_symbol();
    g_debug("Inserting element '%s' at %i %i.",el_name,x,y);
    //1. Find what we've clicked at
    auto click_result = this->widget_data->resolve_click(x, y);
    try{
        auto [bond_or_atom,molecule_idx] = click_result.value();
        if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
            auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
            g_debug("Resolved insertion destination atom: idx=%i, symbol=%s",atom.idx,atom.symbol.c_str());
        } else { // a bond
            auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
            g_warning("TODO: Implement handling insertion at bonds");
        }
    } catch(std::exception& e) {
        // Nothing has been clicked on.
        g_debug("The click could not be resolved to any atom or bond.");
    }
}

ElementInsertion::ElementInsertion(ElementInsertion::Element el) noexcept {
    this->element = el;
}

ElementInsertion::Element ElementInsertion::get_element() const noexcept {
    return this->element;
}

const char* ElementInsertion::get_element_symbol() const noexcept {
    const char* el_name;
    switch (this->get_element()) {
        case ElementInsertion::Element::C:{
            el_name = "C";
            break;
        }
        case ElementInsertion::Element::N:{
            el_name = "N";
            break;
        }
        case ElementInsertion::Element::O:{
            el_name = "O";
            break;
        }
        case ElementInsertion::Element::S:{
            el_name = "S";
            break;
        }
        case ElementInsertion::Element::P:{
            el_name = "P";
            break;
        }
        case ElementInsertion::Element::H:{
            el_name = "H";
            break;
        }
        case ElementInsertion::Element::F:{
            el_name = "F";
            break;
        }
        case ElementInsertion::Element::Cl:{
            el_name = "Cl";
            break;
        }
        case ElementInsertion::Element::Br:{
            el_name = "Br";
            break;
        }
        case ElementInsertion::Element::I:{
            el_name = "I";
            break;
        }
        case ElementInsertion::Element::X:{
            //todo: fixme
            el_name = "X";
            break;
        }
        default: {
            break;
        }
    }
    return el_name;
}