#include "tools.hpp"

using namespace coot::ligand_editor_canvas;

ActiveTool::ActiveTool() noexcept {
    this->variant = ActiveTool::Variant::None;
}

ActiveTool::ActiveTool(ElementInsertion insertion) noexcept {
    this->variant = ActiveTool::Variant::ElementInsertion;
    this->element_insertion = insertion;
}

ActiveTool::Variant ActiveTool::get_variant() const noexcept {
    return this->variant;
}

void ActiveTool::insert_atom(int x, int y) noexcept {
    // optional safety check (with exceptions)
    auto& element_insertion = this->element_insertion;
    const char* el_name = element_insertion.get_element_symbol();
    g_debug("Inserting element '%s' at %i %i.",el_name,x,y);
    //1. Find what we've clicked at
    //2. Decide what to do based on that
    //3. Do it
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