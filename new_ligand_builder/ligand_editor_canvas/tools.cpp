#include "tools.hpp"
#include "core.hpp"
#include "model.hpp"
#include <exception>
#include <stdexcept>
#include <utility>
#include <variant>
#include <rdkit/GraphMol/MolOps.h>

using namespace coot::ligand_editor_canvas;

ActiveTool::ActiveTool() noexcept {
    this->variant = ActiveTool::Variant::None;
}

BondModifier::BondModifier(BondModifierMode mode) noexcept {
    this->mode = mode;
}

MoveTool::MoveTool() noexcept {
    this->in_move = false;
}

void MoveTool::begin_move(int x, int y) noexcept {
    this->prev_move_pos = std::make_pair(x,y);
    this->in_move = true;
    this->current_move_pos = std::make_pair(x, y);
}

std::pair<int,int> MoveTool::end_move() {
    this->in_move = false;
    auto ret = this->get_current_offset();
    this->current_move_pos = std::nullopt;
    this->prev_move_pos = std::nullopt;
    this->canvas_mol_idx = std::nullopt;
    return ret.value();
}

void MoveTool::update_current_move_pos(int x, int y) noexcept {
    this->prev_move_pos = this->current_move_pos;
    this->current_move_pos = std::make_pair(x, y);
}

std::optional<std::pair<int,int>> MoveTool::get_current_offset() const {
    if(!this->current_move_pos.has_value() || !this->prev_move_pos.has_value()) {
        return std::nullopt;
    }
    auto [x1,y1] = this->prev_move_pos.value();
    auto [x2,y2] = this->current_move_pos.value();
    return std::make_pair(x2 - x1, y2 - y1);
}

bool MoveTool::is_in_move() const noexcept {
    return this->in_move;
}

void MoveTool::set_canvas_molecule_index(unsigned int idx) noexcept {
    this->canvas_mol_idx = idx;
}

std::optional<unsigned int> MoveTool::get_canvas_molecule_index() const noexcept {
    return this->canvas_mol_idx;
}


CanvasMolecule::BondType BondModifier::get_target_bond_type() const noexcept {
    switch (this->mode) {
        case BondModifierMode::Single:{
            return CanvasMolecule::BondType::Single;
        }
        case BondModifierMode::Double:{
            return CanvasMolecule::BondType::Double;
        }
        case BondModifierMode::Triple:{
            return CanvasMolecule::BondType::Triple;
        }
    }
}

ActiveTool::ActiveTool(ElementInsertion insertion) noexcept {
    this->variant = ActiveTool::Variant::ElementInsertion;
    this->element_insertion = insertion;
}

ActiveTool::ActiveTool(BondModifier modifier) noexcept {
    this->variant = ActiveTool::Variant::BondModifier;
    this->bond_modifier = modifier;
}

ActiveTool::ActiveTool(DeleteTool deltool) noexcept {
    this->variant = ActiveTool::Variant::Delete;
    this->delete_tool = deltool;
}

ActiveTool::ActiveTool(ChargeModifier chargemod) noexcept {
    this->variant = ActiveTool::Variant::ChargeModifier;
    this->charge_modifier = chargemod;
}

ActiveTool::ActiveTool(MoveTool mov) noexcept {
    this->variant = ActiveTool::Variant::MoveTool;
    this->move_tool = mov;
}

ActiveTool::Variant ActiveTool::get_variant() const noexcept {
    return this->variant;
}

void ActiveTool::set_core_widget_data(impl::CootLigandEditorCanvasPriv* owning_widget) noexcept {
    this->widget_data = static_cast<impl::WidgetCoreData*>(owning_widget);
}

void ActiveTool::check_variant(ActiveTool::Variant expected) const {
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
    if(click_result.has_value()) {
        try{
            auto [bond_or_atom,molecule_idx] = click_result.value();
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
                g_debug("Resolved insertion destination atom: idx=%i, symbol=%s",atom.idx,atom.symbol.c_str());
                auto& rdkit_mol = this->widget_data->rdkit_molecules->at(molecule_idx);
                auto* new_atom = new RDKit::Atom(std::string(el_name));
                rdkit_mol->replaceAtom(atom.idx, new_atom);
                auto& canvas_mol = this->widget_data->molecules->at(molecule_idx);
                canvas_mol.lower_from_rdkit();
            } else { // a bond
                auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
                g_warning("TODO: Implement handling insertion at bonds (if any should happen)");
            }
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
        }
    } else {
        g_debug("The click could not be resolved to any atom or bond.");
    }
}

void ActiveTool::alter_bond(int x, int y) {
    check_variant(Variant::BondModifier);
    BondModifier& mod = this->bond_modifier;
    g_warning("TODO: Implement ActiveTool::alter_bond");
    
    auto click_result = this->widget_data->resolve_click(x, y);
    if(click_result.has_value()) {
        try{
            auto [bond_or_atom,molecule_idx] = click_result.value();
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                //g_warning("The BondModifier tool does not operate on atoms. Nothing to do.");
                auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
                //g_debug("Resolved insertion destination atom: idx=%i, symbol=%s",atom.idx,atom.symbol.c_str());
                auto& rdkit_mol = this->widget_data->rdkit_molecules->at(molecule_idx);
                auto* new_atom = new RDKit::Atom(6);
                RDKit::MolOps::Kekulize(*rdkit_mol.get());
                auto new_atom_idx = rdkit_mol->addAtom(new_atom,false,true);
                rdkit_mol->addBond(new_atom_idx,atom.idx,CanvasMolecule::bond_type_to_rdkit(mod.get_target_bond_type()));
                g_info("New atom added: idx=%i",new_atom_idx);
                RDKit::MolOps::sanitizeMol(*rdkit_mol.get());
                auto& canvas_mol = this->widget_data->molecules->at(molecule_idx);
                canvas_mol.lower_from_rdkit();
            } else {
                auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
                auto& rdkit_mol = this->widget_data->rdkit_molecules->at(molecule_idx);
                auto* rdkit_bond = rdkit_mol->getBondBetweenAtoms(bond.first_atom_idx,bond.second_atom_idx);
                RDKit::MolOps::Kekulize(*rdkit_mol.get());
                auto old_type = rdkit_bond->getBondType();
                rdkit_bond->setBondType(CanvasMolecule::bond_type_to_rdkit(mod.get_target_bond_type()));
                try {
                    RDKit::MolOps::sanitizeMol(*rdkit_mol.get());
                }catch(std::exception& e) {
                    // rollback
                    g_warning("Rolling back invalid molecule change that makes it unable to sanitize with the following error: %s",e.what());
                    rdkit_bond->setBondType(old_type);
                    RDKit::MolOps::sanitizeMol(*rdkit_mol.get());
                    // rethrow
                    throw std::runtime_error(std::string("Invalid bond modification: ") + e.what());
                }
                auto& canvas_mol = this->widget_data->molecules->at(molecule_idx);
                canvas_mol.lower_from_rdkit();
            }
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
        }
    } else {
        // Nothing has been clicked on.
        g_debug("The click could not be resolved to any atom or bond.");
    }
}

void ActiveTool::alter_charge(int x, int y) {
    check_variant(Variant::ChargeModifier);
    g_warning("TODO: Implement ActiveTool::alter_charge");
    auto click_result = this->widget_data->resolve_click(x, y);
    if(click_result.has_value()) {
        try{
            auto [bond_or_atom,molecule_idx] = click_result.value();
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
            } else { // a bond
                auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
            }
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
        }
    } else {
        // Nothing has been clicked on.
        g_debug("The click could not be resolved to any atom or bond.");
    }
}

void ActiveTool::delete_at(int x, int y) {
    check_variant(Variant::Delete);
    g_warning("TODO: Implement ActiveTool::delete_at");
    auto click_result = this->widget_data->resolve_click(x, y);
    if(click_result.has_value()) {
        try{
            auto [bond_or_atom,molecule_idx] = click_result.value();
            auto& rdkit_mol = this->widget_data->rdkit_molecules->at(molecule_idx);
            RDKit::MolOps::Kekulize(*rdkit_mol.get());
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
                rdkit_mol->removeAtom(atom.idx);
            } else { // a bond
                auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
                rdkit_mol->removeBond(bond.first_atom_idx, bond.second_atom_idx);
            }
            RDKit::MolOps::sanitizeMol(*rdkit_mol.get());
            auto& canvas_mol = this->widget_data->molecules->at(molecule_idx);
            canvas_mol.lower_from_rdkit();
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
        }
    } else {
        // Nothing has been clicked on.
        g_debug("The click could not be resolved to any atom or bond.");
    }
}

void ActiveTool::insert_structure(int x, int y) {
    check_variant(Variant::StructureInsertion);
    g_warning("TODO: Implement ActiveTool::insert_structure");
    auto click_result = this->widget_data->resolve_click(x, y);
    if(click_result.has_value()) {
        try{
            auto [bond_or_atom,molecule_idx] = click_result.value();
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
            } else { // a bond
                auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
            }
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
        }
    } else {
        // Nothing has been clicked on.
        g_debug("The click could not be resolved to any atom or bond.");
    }
}

void ActiveTool::update_move_cursor_pos(int x, int y) {
    check_variant(Variant::MoveTool);
    auto& move_tool = this->move_tool;
    if(move_tool.is_in_move()) {
        move_tool.update_current_move_pos(x, y);
        auto [offset_x,offset_y] = move_tool.get_current_offset().value();
        auto mol_idx_opt = move_tool.get_canvas_molecule_index();
        this->widget_data->molecules->at(mol_idx_opt.value()).apply_canvas_translation(offset_x, offset_y);
    }
}

void ActiveTool::end_move() {
    check_variant(Variant::MoveTool);
    auto& move_tool = this->move_tool;
    if(move_tool.is_in_move()) {
        auto mol_idx_opt = move_tool.get_canvas_molecule_index();
        auto [offset_x,offset_y] = move_tool.end_move();
        this->widget_data->molecules->at(mol_idx_opt.value()).apply_canvas_translation(offset_x, offset_y);
    }
}

void ActiveTool::begin_move(int x, int y) {
    check_variant(Variant::MoveTool);
    auto& move_tool = this->move_tool;
    auto mol_opt = this->widget_data->resolve_click(x, y);
    if(mol_opt.has_value()) {
        auto [atom_or_bond,mol_id] = mol_opt.value();
        move_tool.begin_move(x, y);
        move_tool.set_canvas_molecule_index(mol_id);
    }
}

bool ActiveTool::is_in_move() const {
    check_variant(Variant::MoveTool);
    auto& move_tool = this->move_tool;
    return move_tool.is_in_move();
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