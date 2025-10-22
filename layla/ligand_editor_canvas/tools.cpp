/* layla/ligand_editor_canvas/tools.cpp
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

#include "tools.hpp"
#include "core.hpp"
#include "model.hpp"
#include <exception>
#include <algorithm>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <boost/range/iterator_range.hpp>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/PeriodicTable.h>
#include <rdkit/GraphMol/QueryOps.h>
#include "../ligand_editor_canvas.hpp"
#include "../utils.hpp"

using namespace coot::ligand_editor_canvas;

void Tool::on_load(impl::WidgetCoreData& widget_data) {
    // nothing by default
}

void Tool::on_click(ClickContext& ctx, int x, int y) {
    // nothing by default
}

void Tool::on_blank_space_click(ClickContext& ctx, int x, int y) {
    g_debug("The click could not be resolved to any atom or bond.");
}

void Tool::on_release(ClickContext& ctx, int x, int y) {
    // nothing by default
}

void Tool::after_molecule_click(MoleculeClickContext& ctx) {
    // nothing by default
}

bool Tool::on_molecule_click(MoleculeClickContext& ctx) {
    // nothing by default
    return true;
}

void Tool::on_bond_click(MoleculeClickContext& ctx, CanvasMolecule::Bond&) {
    // nothing by default
    g_debug("The tool does not operate on bonds.");
}

void Tool::on_atom_click(MoleculeClickContext& ctx, CanvasMolecule::Atom&) {
    // nothing by default
    g_debug("The tool does not operate on atoms.");
}

void Tool::on_right_click(ClickContext& ctx, int x, int y) {
    // nothing by default
}

void Tool::on_blank_space_right_click(ClickContext& ctx, int x, int y) {
    g_debug("The click could not be resolved to any atom or bond.");
}

bool Tool::on_molecule_right_click(MoleculeClickContext& ctx) {
    // nothing by default
    return true;
}

void Tool::on_bond_right_click(MoleculeClickContext& ctx, CanvasMolecule::Bond&) {
    // nothing by default
    g_debug("The tool does not handle right-click on bonds.");
}

void Tool::on_atom_right_click(MoleculeClickContext& ctx, CanvasMolecule::Atom&) {
    // nothing by default
    g_debug("The tool does not handle right-click on atoms.");
}

void Tool::after_molecule_right_click(MoleculeClickContext& ctx) {
    // nothing by default
}

bool Tool::on_hover(ClickContext& ctx, int x, int y) {
    // not processing hover events by default
    return false;
}

void Tool::on_blank_space_hover(ClickContext& ctx, int x, int y) {
    // nothing by default
}

bool Tool::on_molecule_hover(MoleculeClickContext& ctx) {
    // nothing by default
    return true;
}

void Tool::on_bond_hover(MoleculeClickContext& ctx, CanvasMolecule::Bond&) {
    // nothing by default
}

void Tool::on_atom_hover(MoleculeClickContext& ctx, CanvasMolecule::Atom&) {
    // nothing by default
}

void Tool::after_molecule_hover(MoleculeClickContext& ctx) {
    // nothing by default
}

std::string Tool::get_exception_message_prefix() const noexcept {
    return "An error occured: ";
}

Tool::~Tool() {

}

Tool::ClickContext::ClickContext(impl::WidgetCoreData& widget_data) :widget_data(widget_data) {

}

Tool::MoleculeClickContext::MoleculeClickContext(ClickContext super, unsigned int mol_idx, std::shared_ptr<RDKit::RWMol>& rdkit_mol, CanvasMolecule& canvas_mol) 
:ClickContext(super), rdkit_mol(rdkit_mol), canvas_mol(canvas_mol) {
    this->mol_idx = mol_idx;
}

void ActiveTool::on_load() {
    if(this->tool) {
        this->tool->on_load(*this->widget_data);
    }
}

void ActiveTool::on_click(bool alt_pressed, bool ctrl_pressed, int x, int y, bool right_click) {
    if(!this->tool) {
        return;
    }
    Tool::ClickContext ctx(*this->widget_data);
    ctx.control_pressed = ctrl_pressed;
    ctx.alt_pressed = alt_pressed;


    if(!right_click) {
        this->tool->on_click(ctx, x, y);
    } else {
        this->tool->on_right_click(ctx, x, y);
    }
    auto click_result = this->widget_data->resolve_click(x, y);
    if(click_result.has_value()) {
        auto [bond_or_atom,molecule_idx] = *click_result;
        auto& rdkit_mol = *this->widget_data->rdkit_molecules->at(molecule_idx);
        auto& canvas_mol = *this->widget_data->molecules->at(molecule_idx);
        try{
            Tool::MoleculeClickContext mctx(ctx, molecule_idx, rdkit_mol, canvas_mol);
            if(!right_click) {
                if(!this->tool->on_molecule_click(mctx)) {
                    return;
                }
            } else {
                if(!this->tool->on_molecule_right_click(mctx)) {
                    return;
                }
            }
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
                if(!right_click) {
                    this->tool->on_atom_click(mctx, atom);
                } else {
                    this->tool->on_atom_right_click(mctx, atom);
                }
            } else { // a bond
                auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
                if(!right_click) {
                    this->tool->on_bond_click(mctx, bond);
                } else {
                    this->tool->on_bond_right_click(mctx, bond);
                }
            }
            if(!right_click) {
                this->tool->after_molecule_click(mctx);
            } else {
                this->tool->after_molecule_right_click(mctx);
            }
        // } catch(RDKit::AtomSanitizeException& e) {
        //     g_warning("Invalid change: %s",e.what());
        //     std::string msg = "Invalid change: "; msg += e.what();
        //     this->widget_data->update_status(msg.c_str());
        //     this->widget_data->rollback_current_edition();
        //     canvas_mol.highlight_atom(e.getAtomIdx(), CanvasMolecule::HighlightType::Error);
            
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
            std::string msg = this->tool->get_exception_message_prefix() + e.what();
            this->widget_data->update_status(msg.c_str());
            this->widget_data->rollback_current_edition();
        }
    } else {
        if(!right_click) {
            this->tool->on_blank_space_click(ctx, x, y);
        } else {
            this->tool->on_blank_space_right_click(ctx, x, y);
        }
    }
}

void ActiveTool::on_release(bool alt_pressed, bool ctrl_pressed, int x, int y, bool right_click) {
    if(!this->tool) {
        return;
    }
    Tool::ClickContext ctx(*this->widget_data);
    ctx.control_pressed = ctrl_pressed;
    ctx.alt_pressed = alt_pressed;

    if(right_click) {
        g_warning("todo: Add support for releasing right-clicks in the tools API (when needed).");
        return;
    }

    this->tool->on_release(ctx, x, y);
}

void ActiveTool::on_hover(bool alt_pressed, bool ctrl_pressed, int x, int y) {
    if(!this->tool) {
        return;
    }
    Tool::ClickContext ctx(*this->widget_data);
    ctx.control_pressed = ctrl_pressed;
    ctx.alt_pressed = alt_pressed;

    if(!this->tool->on_hover(ctx, x, y)) {
        return;
    }
    auto hover_result = this->widget_data->resolve_click(x, y);
    if(hover_result.has_value()) {
        auto [bond_or_atom,molecule_idx] = *hover_result;
        auto& rdkit_mol = *this->widget_data->rdkit_molecules->at(molecule_idx);
        auto& canvas_mol = *this->widget_data->molecules->at(molecule_idx);
        try{
            Tool::MoleculeClickContext mctx(ctx, molecule_idx, rdkit_mol, canvas_mol);
            if(!this->tool->on_molecule_hover(mctx)) {
                return;
            }
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
                this->tool->on_atom_hover(mctx, atom);
            } else { // a bond
                auto bond = std::get<CanvasMolecule::Bond>(std::move(bond_or_atom));
                this->tool->on_bond_hover(mctx, bond);
            }
            this->tool->after_molecule_hover(mctx);     
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
            std::string msg = this->tool->get_exception_message_prefix() + e.what();
            this->widget_data->update_status(msg.c_str());
            this->widget_data->rollback_current_edition();
        }
    } else {
        this->tool->on_blank_space_hover(ctx, x, y);
    }
}


void TransformTool::on_click(ClickContext& ctx, int x, int y) {
    auto mol_opt = ctx.widget_data.resolve_click(x, y);
    if(mol_opt.has_value()) {
        auto [atom_or_bond,mol_id] = mol_opt.value();
        this->transform_manager->begin_transform(x, y, this->mode);
        this->transform_manager->set_canvas_molecule_index(mol_id);
        ctx.widget_data.begin_edition();
    }
}

bool TransformTool::on_molecule_click(MoleculeClickContext& ctx) {
    return false;
}

TransformTool::TransformTool(TransformManager::Mode mode) noexcept {
    this->mode = mode;
    this->transform_manager = nullptr;
}

void TransformTool::set_transform_manager(TransformManager* mgr) noexcept {
    this->transform_manager = mgr;
}


BondModifier::BondModifier(BondModifierMode mode) noexcept {
    this->mode = mode;
    this->molecule_idx_and_first_atom_of_new_bond = std::nullopt;
    this->is_in_drag = false;
}

bool BondModifier::is_creating_bond() const noexcept {
    return this->is_in_drag;
}

void BondModifier::begin_creating_bond(unsigned int molecule_idx, unsigned int atom_idx) noexcept {
    this->is_in_drag = true;
    this->molecule_idx_and_first_atom_of_new_bond = std::make_pair(molecule_idx,atom_idx);
}

void BondModifier::finish_creating_bond() noexcept {
    this->is_in_drag = false;
    this->molecule_idx_and_first_atom_of_new_bond = std::nullopt;
}

std::optional<std::pair<unsigned int,unsigned int>> BondModifier::get_molecule_idx_and_first_atom_of_new_bond() const noexcept {
    return this->molecule_idx_and_first_atom_of_new_bond;
}


TransformManager::TransformManager() noexcept {
    this->state = IdleState();
}

bool TransformManager::is_active() const noexcept {
    return !std::holds_alternative<IdleState>(this->state);
}

void TransformManager::begin_transform(int x, int y, Mode mode) noexcept {
    switch(mode) {
        case Mode::Rotation: {
            auto rot = RotationState();
            rot.original_rotation_pos = std::make_pair(x,y);
            rot.current_rotation_pos = std::make_pair(x, y);
            rot.last_absolute_angle = 0.f;
            this->state = rot;
            break;
        }
        case Mode::Translation: {
            auto tr = TranslationState();
            tr.prev_move_pos = std::make_pair(x,y);
            tr.current_move_pos = std::make_pair(x, y);
            this->state = tr;
            break;
        }
    }
}

double TransformManager::RotationState::get_current_absolute_angle(bool snap_to_angle) const {
    auto [x1,y1] = this->original_rotation_pos;
    auto [x2,y2] = this->current_rotation_pos;
    auto diff_x = x2 - x1;
    auto diff_y = y1 - y2;
    auto diff_original = ((double)(diff_x + diff_y)) / 125.0;
    if(snap_to_angle) {
        const double snap = 15.f / 180 * M_PI;
        int divided = diff_original / snap;
        if(divided != 0)  {
            diff_original = divided * snap;
        } else {
            diff_original = 0;
        }
    }
    return diff_original;
}

void TransformManager::update_current_cursor_pos(int x, int y, bool snap) noexcept {
    RotationState* rot = std::get_if<RotationState>(&this->state);
    if(rot) {
        rot->last_absolute_angle = rot->get_current_absolute_angle(snap);
        rot->current_rotation_pos = std::make_pair(x, y);
        return;
    }
    TranslationState* tr = std::get_if<TranslationState>(&this->state);
    if(tr) {
        tr->prev_move_pos = tr->current_move_pos;
        tr->current_move_pos = std::make_pair(x, y);
        return;
    }
}

void TransformManager::end_transform() noexcept {
    this->state = IdleState();
    this->canvas_mol_idx = std::nullopt;
}

void TransformManager::set_canvas_molecule_index(unsigned int idx) noexcept {
    this->canvas_mol_idx = idx;
}

void TransformManager::apply_current_transform_state(impl::WidgetCoreData* widget_data, bool snap_to_angle, bool about_to_end) const {
    if(!this->canvas_mol_idx.has_value()) {
        return;
    }
    auto& mol = *widget_data->molecules->at(this->canvas_mol_idx.value());
    const RotationState* rot = std::get_if<RotationState>(&this->state);
    if(rot) {
        auto angle = rot->get_current_angle_diff(snap_to_angle);
        auto abs_angle = rot->get_current_absolute_angle(snap_to_angle) / M_PI * 180;
        mol.rotate_by_angle(angle);
        mol.lower_from_rdkit(!widget_data->allow_invalid_molecules, false);
        std::string msg;
        if (about_to_end) {
            msg = "Molecule rotated by: " + std::to_string(abs_angle) + " degrees.";
        } else {
            msg = "Current molecule rotation: " + std::to_string(abs_angle) + " degrees.";
        }
        widget_data->update_status(msg.c_str());
        return;
    }
    const TranslationState* tr = std::get_if<TranslationState>(&this->state);
    if(tr) {
        auto [offset_x,offset_y] = tr->get_current_offset();
        mol.apply_canvas_translation(offset_x, offset_y, widget_data->scale);
        return;
    }
}

std::pair<int,int> TransformManager::TranslationState::get_current_offset() const {
    auto [x1,y1] = this->prev_move_pos;
    auto [x2,y2] = this->current_move_pos;
    return std::make_pair(x2 - x1, y2 - y1);
}

double TransformManager::RotationState::get_current_angle_diff(bool snap_to_angle) const {
    auto original = this->get_current_absolute_angle(snap_to_angle);
    auto diff = original - this->last_absolute_angle;
    return diff;
}

StructureInsertion::StructureInsertion(StructureInsertion::Structure st) noexcept 
:structure(st) {

}

StructureInsertion::Structure StructureInsertion::get_structure() const noexcept {
    return this->structure;
}


CanvasMolecule::BondType BondModifier::get_target_bond_type() const noexcept {
    switch (this->mode) {
        default:
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

ActiveTool::ActiveTool() noexcept {
    // nothing
}

ActiveTool::ActiveTool(ElementInsertion insertion) noexcept {
    this->tool = std::make_unique<ElementInsertion>(std::move(insertion));
}

ActiveTool::ActiveTool(BondModifier modifier) noexcept {
    this->tool = std::make_unique<BondModifier>(std::move(modifier));
}

ActiveTool::ActiveTool(DeleteTool deltool) noexcept {
    this->tool = std::make_unique<DeleteTool>(std::move(deltool));
}

ActiveTool::ActiveTool(ChargeModifier chargemod) noexcept {
    this->tool = std::make_unique<ChargeModifier>(std::move(chargemod));
}

ActiveTool::ActiveTool(TransformTool mov) noexcept {
    mov.set_transform_manager(&this->transform_manager);
    this->tool = std::make_unique<TransformTool>(std::move(mov));
}

ActiveTool::ActiveTool(StructureInsertion insertion) noexcept {
    this->tool = std::make_unique<StructureInsertion>(std::move(insertion));
}

ActiveTool::ActiveTool(GeometryModifier modifier) noexcept {
    this->tool = std::make_unique<GeometryModifier>(std::move(modifier));
}

ActiveTool::ActiveTool(FormatTool fmt) noexcept {
    this->tool = std::make_unique<FormatTool>(std::move(fmt));
}

ActiveTool::ActiveTool(FlipTool flip) noexcept {
    this->tool = std::make_unique<FlipTool>(std::move(flip));
}

ActiveTool::ActiveTool(RemoveHydrogensTool rh) noexcept {
    this->tool = std::make_unique<RemoveHydrogensTool>(std::move(rh));
}

void ActiveTool::set_core_widget_data(impl::CootLigandEditorCanvasPriv* owning_widget) noexcept {
    this->widget_data = static_cast<impl::WidgetCoreData*>(owning_widget);
}

bool ElementInsertion::on_molecule_click(MoleculeClickContext& ctx) {
    ctx.widget_data.begin_edition();
    return true;
}

void ElementInsertion::on_bond_click(MoleculeClickContext& ctx, CanvasMolecule::Bond& bond) {
    unsigned int atomic_number = this->get_atomic_number();
    // auto el_name = RDKit::PeriodicTable::getTable()->getElementSymbol(atomic_number);
    auto* new_atom = new RDKit::Atom(atomic_number);
    ctx.rdkit_mol->removeBond(bond.first_atom_idx, bond.second_atom_idx);
    unsigned int new_atom_idx = ctx.rdkit_mol->addAtom(new_atom, true, true);
    ctx.rdkit_mol->addBond(bond.first_atom_idx, new_atom_idx, RDKit::Bond::SINGLE);
    ctx.rdkit_mol->addBond(new_atom_idx, bond.second_atom_idx, RDKit::Bond::SINGLE);
    ctx.widget_data.update_status("Atom has been inserted.");
}

void ElementInsertion::on_atom_click(MoleculeClickContext& ctx, CanvasMolecule::Atom& atom) {
    unsigned int atomic_number = this->get_atomic_number();
    auto el_name = RDKit::PeriodicTable::getTable()->getElementSymbol(atomic_number);
    g_debug("Appending element '%u' (%s) to destination atom: idx=%i, symbol=%s.",atomic_number,el_name.c_str(),atom.idx,atom.symbol.c_str());
    auto* new_atom = new RDKit::Atom(std::string(el_name));
    ctx.rdkit_mol->replaceAtom(atom.idx, new_atom);
    ctx.widget_data.update_status("Atom has been replaced.");
}

void ElementInsertion::after_molecule_click(MoleculeClickContext& ctx) {
    ctx.canvas_mol.lower_from_rdkit(!ctx.widget_data.allow_invalid_molecules);
    ctx.widget_data.finalize_edition();
}

std::string ElementInsertion::get_exception_message_prefix() const noexcept {
    return "Could not insert atom: ";
}

void ElementInsertion::on_blank_space_click(ClickContext& ctx, int x, int y) {
    g_debug("The click could not be resolved to any atom or bond.");
    // This 'if' should be removed once we implement merging molecules
    if(ctx.widget_data.get_molecule_count_impl() == 0) {
        g_debug("There are no molecules. Element insertion will therefore create a new one.");
        auto rdkit_mol = std::make_shared<RDKit::RWMol>();
        rdkit_mol->addAtom(new RDKit::Atom(this->get_atomic_number()),false,true);
        RDKit::MolOps::sanitizeMol(*rdkit_mol);
        // append_molecule() already calls begin_edition() and finalize_edition()
        // ctx.widget_data.begin_edition();
        #ifndef __EMSCRIPTEN__
        auto* widget_ptr = static_cast<impl::CootLigandEditorCanvasPriv*>(&ctx.widget_data);
        coot_ligand_editor_canvas_append_molecule(COOT_COOT_LIGAND_EDITOR_CANVAS(widget_ptr), rdkit_mol);
        #else // __EMSCRIPTEN__ defined
        // Lhasa-specific includes/definitions
        auto* widget_ptr = static_cast<::CootLigandEditorCanvas*>(&ctx.widget_data);
        coot_ligand_editor_canvas_append_molecule(widget_ptr, rdkit_mol);
        #endif
        ctx.widget_data.update_status("New molecule created from an atom.");
    } else {
        g_debug("There are already molecules. Element insertion on blank space is a no-op.");
    }
}

bool BondModifier::on_molecule_click(MoleculeClickContext& ctx) {
    ctx.widget_data.begin_edition();
    return true;
}

void BondModifier::on_bond_click(MoleculeClickContext& ctx, CanvasMolecule::Bond& bond) {
    auto* rdkit_bond = ctx.rdkit_mol->getBondBetweenAtoms(bond.first_atom_idx,bond.second_atom_idx);
    RDKit::MolOps::Kekulize(*ctx.rdkit_mol.get());
    auto old_type = rdkit_bond->getBondType();
    rdkit_bond->setBondType(CanvasMolecule::bond_type_to_rdkit(this->get_target_bond_type()));
    try {
        this->sanitize_molecule(ctx.widget_data, *ctx.rdkit_mol.get());
    }catch(std::exception& e) {
        // rollback
        g_warning("Rolling back invalid molecule change that makes it unable to sanitize with the following error: %s",e.what());
        rdkit_bond->setBondType(old_type);
        this->sanitize_molecule(ctx.widget_data, *ctx.rdkit_mol.get());
        // rethrow
        throw std::runtime_error(std::string("Invalid bond modification: ") + e.what());
    }
    ctx.widget_data.update_status("Bond has been altered.");
    ctx.canvas_mol.lower_from_rdkit(!ctx.widget_data.allow_invalid_molecules);
    ctx.widget_data.finalize_edition();
}

void BondModifier::on_atom_click(MoleculeClickContext& ctx, CanvasMolecule::Atom& atom) {
    this->begin_creating_bond(ctx.mol_idx, atom.idx);
}

std::string BondModifier::get_exception_message_prefix() const noexcept {
    return "Could not alter/create bond: ";
}

void BondModifier::on_blank_space_click(ClickContext& ctx, int x, int y) {
    g_debug("The click could not be resolved to any atom or bond.");
    // This 'if' should be removed once we implement merging molecules
    if(ctx.widget_data.get_molecule_count_impl() == 0) {
        g_debug("There are no molecules. Element insertion will therefore create a new one.");
        auto rdkit_mol = std::make_shared<RDKit::RWMol>();
        auto first_carbon_idx = rdkit_mol->addAtom(new RDKit::Atom(6),false,true);
        auto second_carbon_idx = rdkit_mol->addAtom(new RDKit::Atom(6),false,true);
        rdkit_mol->addBond(first_carbon_idx,second_carbon_idx,CanvasMolecule::bond_type_to_rdkit(this->get_target_bond_type()));
        // This theoretically may throw but it has no reason to as we just created a valid molecule
        RDKit::MolOps::sanitizeMol(*rdkit_mol);
        // append_molecule() already calls begin_edition() and finalize_edition()
        // ctx.widget_data.begin_edition();
        #ifndef __EMSCRIPTEN__
        auto* widget_ptr = static_cast<impl::CootLigandEditorCanvasPriv*>(&ctx.widget_data);
        coot_ligand_editor_canvas_append_molecule(COOT_COOT_LIGAND_EDITOR_CANVAS(widget_ptr), rdkit_mol);
        #else // __EMSCRIPTEN__ defined
        // Lhasa-specific includes/definitions
        auto* widget_ptr = static_cast<::CootLigandEditorCanvas*>(&ctx.widget_data);
        coot_ligand_editor_canvas_append_molecule(widget_ptr, rdkit_mol);
        #endif
        ctx.widget_data.update_status("New molecule created from an atom.");
    } else {
        g_debug("There are already molecules. Element insertion on blank space is a no-op.");
    }
}

bool ActiveTool::is_creating_bond() const noexcept {
    const BondModifier* mod = dynamic_cast<const BondModifier*>(this->tool.get());
    if (mod) {
        return mod->is_creating_bond();
    }
    return false;
}

void BondModifier::on_release(ClickContext& ctx, int x, int y) {
    if(!this->is_creating_bond()) {
        return;
    }
    impl::WidgetCoreData& widget_data = ctx.widget_data;

    auto click_result = widget_data.resolve_click(x, y);
    auto [original_molecule_idx, first_atom_idx] = this->get_molecule_idx_and_first_atom_of_new_bond().value();
    this->finish_creating_bond();
    widget_data.currently_created_bond = std::nullopt;

    if(click_result.has_value()) {
        try{
            auto [bond_or_atom,molecule_idx] = *click_result;
            if(std::holds_alternative<CanvasMolecule::Atom>(bond_or_atom)) {
                auto second_atom = std::get<CanvasMolecule::Atom>(std::move(bond_or_atom));
                if(original_molecule_idx != molecule_idx) {
                    widget_data.update_status("Cannot create bond between different molecules!");
                    widget_data.rollback_current_edition();
                    return;
                }
                auto& rdkit_mol = *widget_data.rdkit_molecules->at(molecule_idx);
                RDKit::MolOps::Kekulize(*rdkit_mol.get());
                if(first_atom_idx == second_atom.idx) {
                    auto* new_atom = new RDKit::Atom(6);
                    auto new_atom_idx = rdkit_mol->addAtom(new_atom,false,true);
                    rdkit_mol->addBond(new_atom_idx,second_atom.idx,CanvasMolecule::bond_type_to_rdkit(this->get_target_bond_type()));
                    g_info("New atom added: idx=%i",new_atom_idx);
                    widget_data.update_status("New carbon atom added.");
                } else {
                    rdkit_mol->addBond(first_atom_idx,second_atom.idx,CanvasMolecule::bond_type_to_rdkit(this->get_target_bond_type()));
                    widget_data.update_status("Created new bond between atoms.");
                }
                this->sanitize_molecule(widget_data, *rdkit_mol.get());
                auto& canvas_mol = *widget_data.molecules->at(molecule_idx);
                canvas_mol.lower_from_rdkit(!widget_data.allow_invalid_molecules);
                widget_data.finalize_edition();
            } else {
                widget_data.update_status("Can't link bond to a bond!");
                widget_data.rollback_current_edition();
            }
        } catch(std::exception& e) {
            g_warning("An error occured: %s",e.what());
            std::string msg = std::string("Could not alter/create bond: ") + e.what();
            widget_data.update_status(msg.c_str());
            widget_data.rollback_current_edition();
        }
    } else {
        std::string msg = "The new bond goes nowhere.";
        widget_data.update_status(msg.c_str());
        widget_data.rollback_current_edition();
    }
}

void GeometryModifier::on_bond_click(MoleculeClickContext& ctx, CanvasMolecule::Bond& bond) {
    ctx.widget_data.begin_edition();
    auto* rdkit_bond = ctx.rdkit_mol->getBondBetweenAtoms(bond.first_atom_idx,bond.second_atom_idx);

    auto bond_geometry = CanvasMolecule::bond_geometry_from_rdkit(rdkit_bond->getBondDir());
    auto new_bond_geometry = CanvasMolecule::cycle_bond_geometry(bond_geometry);
    g_debug("Target bond geometry: %u",static_cast<unsigned int>(new_bond_geometry));
    rdkit_bond->setBondDir(CanvasMolecule::bond_geometry_to_rdkit(new_bond_geometry));

    ctx.widget_data.update_status("Geometry of bond has been altered.");
    ctx.canvas_mol.lower_from_rdkit(!ctx.widget_data.allow_invalid_molecules);
    g_debug("Final bond geometry: %u",static_cast<unsigned int>(CanvasMolecule::bond_geometry_from_rdkit(rdkit_bond->getBondDir())));
    ctx.widget_data.finalize_edition();
}

std::string GeometryModifier::get_exception_message_prefix() const noexcept {
    return "Could not alter bond geometry: ";
}

void ChargeModifier::on_atom_right_click(MoleculeClickContext& ctx, CanvasMolecule::Atom& atom) {
    ctx.widget_data.begin_edition();
    // Do we need this here?
    RDKit::MolOps::Kekulize(*ctx.rdkit_mol.get());
    auto* rdkit_atom = ctx.rdkit_mol->getAtomWithIdx(atom.idx);
    int old_charge = rdkit_atom->getFormalCharge();
    
    rdkit_atom->setFormalCharge(old_charge + 1);

    ctx.widget_data.update_status("Charge of atom has been increased.");

    Tool::sanitize_molecule(ctx.widget_data, *ctx.rdkit_mol.get());
    ctx.canvas_mol.lower_from_rdkit(!ctx.widget_data.allow_invalid_molecules);

    ctx.widget_data.finalize_edition();
}

void ChargeModifier::on_atom_click(MoleculeClickContext& ctx, CanvasMolecule::Atom& atom) {
    ctx.widget_data.begin_edition();
    // Do we need this here?
    RDKit::MolOps::Kekulize(*ctx.rdkit_mol.get());
    auto* rdkit_atom = ctx.rdkit_mol->getAtomWithIdx(atom.idx);
    int old_charge = rdkit_atom->getFormalCharge();

    // const RDKit::PeriodicTable *table = RDKit::PeriodicTable::getTable();
    // auto valence_list = table->getValenceList(rdkit_atom->getAtomicNum());
    // auto valence_list_to_string = [&](){
    //     std::string valence_list_string = "[";
    //     auto vl_it = valence_list.begin();
    //     while(vl_it != valence_list.end()) {
    //         valence_list_string += std::to_string(*vl_it);
    //         vl_it++;
    //         if(vl_it != valence_list.end()) {
    //             valence_list_string += ", ";
    //         }
    //     }
    //     valence_list_string += "]";
    //     return valence_list_string;
    // };

    // auto valence_list_string = valence_list_to_string();
    // g_debug(
    //     "Valence list from RDKit: %s Num of outer shell electrons: %i", 
    //     valence_list_string.c_str(),
    //     table->getNouterElecs(rdkit_atom->getAtomicNum())
    // );

    // g_warning("todo: Fix-up computing plausible charges for atoms in the charge tool.");
    // // We won't have to clear the list when it'll have the right contents (something to be fixed)
    // valence_list.clear();
    // // valence_list.push_back(0);
    // for(int i = -4; i != 5; i++) {
    //     valence_list.push_back(i);
    // }

    // auto it = std::find(valence_list.begin(),valence_list.end(),old_charge);
    // if(it != valence_list.end()) {
    //     it++;
    // }
    // if(it == valence_list.end()) {
    //     it = valence_list.begin();
    // }

    // int new_charge = *it;

    // valence_list_string = valence_list_to_string();
    // g_info("Old formal charge: %i New formal charge: %i List: %s",old_charge,new_charge,valence_list_string.c_str());
    // rdkit_atom->setFormalCharge(new_charge);

    rdkit_atom->setFormalCharge(old_charge - 1);

    ctx.widget_data.update_status("Charge of atom has been decreased.");

    Tool::sanitize_molecule(ctx.widget_data, *ctx.rdkit_mol.get());
    ctx.canvas_mol.lower_from_rdkit(!ctx.widget_data.allow_invalid_molecules);

    ctx.widget_data.finalize_edition();
}

std::string ChargeModifier::get_exception_message_prefix() const noexcept {
    return "Could not alter charge: ";
}

DeleteTool::ListOfAtomsOrBonds DeleteTool::trace_chain_impl(const RDKit::ROMol* mol, std::set<unsigned int>& processed_atoms, RDKit::Atom const* rdatom) {
    DeleteTool::ListOfAtomsOrBonds ret;
    ret.push_back(rdatom->getIdx());
    processed_atoms.emplace(rdatom->getIdx());
    for(const auto& i: boost::make_iterator_range(mol->getAtomNeighbors(rdatom))) {
        ret.push_back(std::make_tuple(rdatom->getIdx(), (unsigned int) i));
        if(processed_atoms.find(i) != processed_atoms.end()) {
            continue;
        }
        processed_atoms.emplace(i);
        // First, push the bond between this atom and the neighbor atom
        auto list = trace_chain_impl(mol, processed_atoms, mol->getAtomWithIdx(i));
        // Append list to ret
        ret.insert(ret.end(), list.begin(), list.end());
    }
    return ret;
}

bool DeleteTool::chain_contains_majority_of_atoms(const  ListOfAtomsOrBonds& chain, const RDKit::ROMol* mol) {
    unsigned int no_atoms = std::count_if(chain.cbegin(), chain.cend(), [](const AtomOrBond& el){
        return std::holds_alternative<unsigned int>(el);
    });
    return no_atoms >= mol->getNumAtoms() / 2;
}

DeleteTool::ListOfAtomsOrBonds DeleteTool::trace_rchain(const MoleculeClickContext& ctx, const CanvasMolecule::Atom& atom) {
    RDKit::Atom const* rdatom = ctx.rdkit_mol->getAtomWithIdx(atom.idx);
    const RDKit::ROMol* mol = ctx.rdkit_mol.get();

    std::set<unsigned int> processed_atoms;
    DeleteTool::ListOfAtomsOrBonds ret;
    ret.push_back(rdatom->getIdx());
    processed_atoms.emplace(rdatom->getIdx());
    unsigned int neighbors_count = 0;
    // Is there no better idea?
    for(const auto& _i: boost::make_iterator_range(mol->getAtomNeighbors(rdatom))) {
        neighbors_count += 1;
    }
    // Skip tracing in-ring atoms
    if(neighbors_count <= 1 || (neighbors_count == 2 && RDKit::queryIsAtomInRing(rdatom))) {
        return ret;
    }

    std::vector<DeleteTool::ListOfAtomsOrBonds> branches;
    for(const auto& i: boost::make_iterator_range(mol->getAtomNeighbors(rdatom))) {
        if(processed_atoms.find(i) != processed_atoms.end()) {
            continue;
        }
        DeleteTool::ListOfAtomsOrBonds list;
        // First, push the bond between this atom and the branch atom
        list.push_back(std::make_tuple(rdatom->getIdx(), (unsigned int) i));
        auto list2 = trace_chain_impl(mol, processed_atoms, mol->getAtomWithIdx(i));
        // Append list2 to list
        list.insert(list.end(), list2.begin(), list2.end());
        branches.push_back(list);
        processed_atoms.emplace(i);
    }
    // Iterator to the shortest branch
    auto it = std::min_element(branches.begin(), branches.end(), [](const std::vector<DeleteTool::AtomOrBond>& a, const std::vector<DeleteTool::AtomOrBond>& b){
        return a.size() < b.size();
    });
    if(it != branches.end()) {
        if(!chain_contains_majority_of_atoms(*it, mol)) {
            // Append the smallest branch to our result
            ret.insert(ret.end(), it->begin(), it->end());
        }
    }
    return ret;
}

DeleteTool::ListOfAtomsOrBonds DeleteTool::trace_rchain(const MoleculeClickContext& ctx, const CanvasMolecule::Bond& bond) {
    const RDKit::ROMol* mol = ctx.rdkit_mol.get();

    DeleteTool::ListOfAtomsOrBonds ret;
    ret.push_back(std::make_tuple(bond.first_atom_idx, bond.second_atom_idx));

    // Skip in-ring bonds
    if(RDKit::queryIsBondInRing(mol->getBondBetweenAtoms(bond.first_atom_idx, bond.second_atom_idx))) {
        return ret;
    }

    std::set<unsigned int> processed_atoms1;
    processed_atoms1.insert(bond.first_atom_idx);
    processed_atoms1.insert(bond.second_atom_idx);
    std::set<unsigned int> processed_atoms2 = processed_atoms1;

    auto case1 = trace_chain_impl(mol, processed_atoms1, mol->getAtomWithIdx(bond.first_atom_idx));
    auto case2 = trace_chain_impl(mol, processed_atoms2, mol->getAtomWithIdx(bond.second_atom_idx));
    
    if(case1.size() <= case2.size() && ! chain_contains_majority_of_atoms(case1, mol)) {
        ret.insert(ret.end(), case1.begin(), case1.end());
    } else if(! chain_contains_majority_of_atoms(case2, mol)) {
        ret.insert(ret.end(), case2.begin(), case2.end());
    }

    return ret;
}

bool DeleteTool::on_hover(ClickContext& ctx, int x, int y) {
    return true;
}

bool DeleteTool::on_molecule_hover(MoleculeClickContext& ctx) {
    if(ctx.control_pressed && !ctx.alt_pressed) {
        // Highlight whole molecule for deletion
        for(unsigned int i = 0; i < ctx.rdkit_mol->getNumAtoms(); i++) {
            ctx.canvas_mol.add_atom_highlight(i, CanvasMolecule::HighlightType::Hover);
        }
        ctx.canvas_mol.add_highlight_to_all_bonds(CanvasMolecule::HighlightType::Hover);
        return false;
    }
    return true;
}

void DeleteTool::on_bond_hover(MoleculeClickContext& ctx, CanvasMolecule::Bond& bond) {
    if(ctx.control_pressed && ctx.alt_pressed) {
        // No need to do anything for single-atom mode
        return;
    }
    auto rchain = trace_rchain(ctx, bond);
    this->highlight_rchain(ctx, rchain);
}

void DeleteTool::on_atom_hover(MoleculeClickContext& ctx, CanvasMolecule::Atom& atom) {
    if(ctx.control_pressed && ctx.alt_pressed) {
        // No need to do anything for single-atom mode
        return;
    }
    auto rchain = trace_rchain(ctx, atom);
    this->highlight_rchain(ctx, rchain);
}

bool DeleteTool::on_molecule_click(MoleculeClickContext& ctx) {
    if(ctx.control_pressed) {
        // Single-atom deletion mode
        if(ctx.alt_pressed) {
            // Cannot do it here
            // widget_data.begin_edition();
            RDKit::MolOps::Kekulize(*ctx.rdkit_mol.get());
            return true;
        }
        // Whole molecule deletion mode
        ctx.widget_data.delete_molecule_with_idx(ctx.mol_idx);
        return false;
    } else {
        // Regular R-chain deletion mode
        return true;
    }
}

void DeleteTool::remove_rchain(const MoleculeClickContext& ctx, const ListOfAtomsOrBonds& chain) {
    // First, remove bonds
    for(const auto& i: chain) {
        if(std::holds_alternative<std::tuple<unsigned int, unsigned int>>(i)) { // bond
            auto [idx_a, idx_b] = std::get<std::tuple<unsigned int, unsigned int>>(i);
            //g_debug("RCHAIN ITEM: BOND %u->%u", idx_a, idx_b);
            ctx.rdkit_mol->removeBond(idx_a, idx_b);
        }
    }
    // Then atoms
    std::vector<unsigned int> atoms;
    for(const auto& i: chain) {
        if(std::holds_alternative<unsigned int>(i)) { // atom
            unsigned int atom_idx = std::get<unsigned int>(i);
            atoms.push_back(atom_idx);
        }
    }
    std::sort(atoms.begin(), atoms.end(), [](unsigned int a, unsigned int b){
        return a >= b;
    });
    
    for(const auto& atom_idx: atoms) {
        //g_debug("RCHAIN ITEM: ATOM %u", atom_idx);
        ctx.rdkit_mol->removeAtom(atom_idx);
        ctx.canvas_mol.update_cached_atom_coordinate_map_after_atom_removal(atom_idx);
    }
}

void DeleteTool::highlight_rchain(const MoleculeClickContext& ctx, const ListOfAtomsOrBonds& chain) {
    for(const auto& i: chain) {
        if(std::holds_alternative<unsigned int>(i)) { // atom
            unsigned int atom_idx = std::get<unsigned int>(i);
            ctx.canvas_mol.add_atom_highlight(atom_idx, CanvasMolecule::HighlightType::Hover);
        } else { // bond
            auto [idx_a, idx_b] = std::get<std::tuple<unsigned int, unsigned int>>(i);
            ctx.canvas_mol.add_bond_highlight(idx_a, idx_b, CanvasMolecule::HighlightType::Hover);
        }
    }
}

void DeleteTool::on_bond_click(MoleculeClickContext& ctx, CanvasMolecule::Bond& bond) {
    ctx.widget_data.begin_edition();
    if(ctx.control_pressed && ctx.alt_pressed) { // Single-atom mode
        ctx.rdkit_mol->removeBond(bond.first_atom_idx, bond.second_atom_idx);
        ctx.widget_data.update_status("Bond has been deleted.");
    } else { // R-chain mode
        auto rchain = trace_rchain(ctx, bond);
        this->remove_rchain(ctx, rchain);
    }
}

void DeleteTool::on_atom_click(MoleculeClickContext& ctx, CanvasMolecule::Atom& atom) {
    if(ctx.rdkit_mol->getNumAtoms() <= 1) {
        // For 1-atom molecules, we switch to the whole-moleculue deletion mode.
        // This is triggered by not beginning edition
        // Look at after_molecule_click for reference
        return;
    }
    ctx.widget_data.begin_edition();
    if(ctx.control_pressed && ctx.alt_pressed) { // Single-atom mode
        ctx.rdkit_mol->removeAtom(atom.idx);
        ctx.canvas_mol.update_cached_atom_coordinate_map_after_atom_removal(atom.idx);
        ctx.widget_data.update_status("Atom has been deleted.");
    } else { // R-chain mode
        auto rchain = trace_rchain(ctx, atom);
        this->remove_rchain(ctx, rchain);
    }
}

void DeleteTool::after_molecule_click(MoleculeClickContext& ctx) {
    if(ctx.widget_data.is_in_edition()) {
        Tool::sanitize_molecule(ctx.widget_data, *ctx.rdkit_mol.get());
        ctx.canvas_mol.lower_from_rdkit(!ctx.widget_data.allow_invalid_molecules);
        ctx.widget_data.finalize_edition();
    } else {
        ctx.widget_data.delete_molecule_with_idx(ctx.mol_idx);
    }
}

std::string DeleteTool::get_exception_message_prefix() const noexcept {
    return "Could not delete atom/bond: ";
}


unsigned int StructureInsertion::append_carbon(RDKit::RWMol* rdkit_mol_ptr, unsigned int target_idx, RDKit::Bond::BondType bond_type) {
    RDKit::Atom* new_carbon = new RDKit::Atom(6);
    auto new_carbon_idx = rdkit_mol_ptr->addAtom(new_carbon,false,true);
    rdkit_mol_ptr->addBond(target_idx,new_carbon_idx,bond_type);
    return new_carbon_idx;  
}

unsigned int StructureInsertion::append_carbon_chain(RDKit::RWMol* rdkit_mol_ptr, unsigned int chain_start_idx, std::size_t atom_count) {
    unsigned int current_atom = chain_start_idx;
    for(std::size_t i = 0; i < atom_count; i++) {
        current_atom = append_carbon(rdkit_mol_ptr, current_atom);
    }
    return current_atom;
}

void StructureInsertion::append_structure_to_atom(RDKit::RWMol* rdkit_mol_ptr, unsigned int atom_idx, bool spiro) const {
    unsigned int first_carbon;
    if(spiro) {
        first_carbon = atom_idx;
    } else {
        first_carbon = append_carbon(rdkit_mol_ptr, atom_idx);
    }
    switch (this->get_structure()) {
        case Structure::CycloPropaneRing: {
            auto third_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,2);
            rdkit_mol_ptr->addBond(first_carbon,third_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloButaneRing: {
            auto last_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,3);
            rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloPentaneRing: {
            auto last_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,4);
            rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloHexaneRing: {
            auto last_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,5);
            rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::BenzeneRing: {
            auto second_carbon = append_carbon(rdkit_mol_ptr, first_carbon);
            auto third_carbon = append_carbon(rdkit_mol_ptr, second_carbon, RDKit::Bond::DOUBLE);
            auto fourth_carbon = append_carbon(rdkit_mol_ptr, third_carbon);
            auto fifth_carbon = append_carbon(rdkit_mol_ptr, fourth_carbon, RDKit::Bond::DOUBLE);
            auto sixth_carbon = append_carbon(rdkit_mol_ptr, fifth_carbon);
            rdkit_mol_ptr->addBond(first_carbon,sixth_carbon,RDKit::Bond::DOUBLE);
            break;
        }
        case Structure::CycloHeptaneRing: {
            auto last_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,6);
            rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloOctaneRing: {
            auto last_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,7);
            rdkit_mol_ptr->addBond(first_carbon,last_carbon,RDKit::Bond::SINGLE);
            break;
        }
    }
}


bool StructureInsertion::on_molecule_click(MoleculeClickContext& ctx) {
    ctx.widget_data.begin_edition();
    // This causes problems with appending benzene to a benzene
    // RDKit::MolOps::Kekulize(*ctx.rdkit_mol.get());
    return true;
}

void StructureInsertion::on_bond_click(MoleculeClickContext& ctx, CanvasMolecule::Bond& bond) {
    auto last_carbon = bond.second_atom_idx;
    auto first_carbon = bond.first_atom_idx;
    auto structure_kind = this->get_structure();
    auto& rdkit_mol = ctx.rdkit_mol;
    auto* rdkit_mol_ptr = rdkit_mol.get();
    switch (structure_kind) {
        case Structure::CycloPropaneRing: {
            auto second_carbon = append_carbon(rdkit_mol_ptr, first_carbon);
            rdkit_mol->addBond(second_carbon,last_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloButaneRing: {
            auto fourth_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,2);
            rdkit_mol->addBond(last_carbon,fourth_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloPentaneRing: {
            auto fifth_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,3);
            rdkit_mol->addBond(last_carbon,fifth_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloHexaneRing: {
            auto sixth_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,4);
            rdkit_mol->addBond(last_carbon,sixth_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::BenzeneRing: {
            auto second_carbon = append_carbon(rdkit_mol_ptr, first_carbon,RDKit::Bond::DOUBLE);
            auto third_carbon = append_carbon(rdkit_mol_ptr, second_carbon);
            auto fourth_carbon = append_carbon(rdkit_mol_ptr, third_carbon,RDKit::Bond::DOUBLE);
            auto fifth_carbon = append_carbon(rdkit_mol_ptr, fourth_carbon);
            rdkit_mol->addBond(last_carbon,fifth_carbon,RDKit::Bond::DOUBLE);
            break;
        }
        case Structure::CycloHeptaneRing: {
            auto seventh_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,5);
            rdkit_mol->addBond(last_carbon,seventh_carbon,RDKit::Bond::SINGLE);
            break;
        }
        case Structure::CycloOctaneRing: {
            auto eigth_carbon = append_carbon_chain(rdkit_mol_ptr, first_carbon,6);
            rdkit_mol->addBond(last_carbon,eigth_carbon,RDKit::Bond::SINGLE);
            break;
        }
    }
    ctx.widget_data.update_status("Carbon ring has been added, adjacent to the bond.");
}

void StructureInsertion::on_atom_click(MoleculeClickContext& ctx, CanvasMolecule::Atom& atom) {
    bool spiro = ctx.control_pressed;
    this->append_structure_to_atom(ctx.rdkit_mol.get(),atom.idx, spiro);
    ctx.widget_data.update_status("Carbon ring has been appended to the atom.");
}

void StructureInsertion::after_molecule_click(MoleculeClickContext& ctx) {
    this->sanitize_molecule(ctx.widget_data, *ctx.rdkit_mol);
    ctx.canvas_mol.lower_from_rdkit(!ctx.widget_data.allow_invalid_molecules);
    ctx.widget_data.finalize_edition();
}

std::string StructureInsertion::get_exception_message_prefix() const noexcept {
    return "Could not insert structure: ";
}

void StructureInsertion::on_blank_space_click(ClickContext& ctx, int x, int y) {
    g_debug("The click could not be resolved to any atom or bond.");
    // This 'if' should be removed once we implement merging molecules
    if(ctx.widget_data.get_molecule_count_impl() == 0) {
        g_debug("There are no molecules. Structure insertion will therefore create a new one.");
        auto rdkit_mol = std::make_shared<RDKit::RWMol>();
        rdkit_mol->addAtom(new RDKit::Atom(6),false,true);
        append_structure_to_atom(rdkit_mol.get(),0, false);
        rdkit_mol->removeAtom((unsigned int) 0);
        // This function calls "begin_edition" and "finalize_edition", 
        // so we can't call "begin_edition" here above.
        RDKit::MolOps::sanitizeMol(*rdkit_mol);

        #ifndef __EMSCRIPTEN__
        auto* widget_ptr = static_cast<impl::CootLigandEditorCanvasPriv*>(&ctx.widget_data);
        coot_ligand_editor_canvas_append_molecule(COOT_COOT_LIGAND_EDITOR_CANVAS(widget_ptr), rdkit_mol);
        #else // __EMSCRIPTEN__ defined
        // Lhasa-specific includes/definitions
        auto* widget_ptr = static_cast<::CootLigandEditorCanvas*>(&ctx.widget_data);
        coot_ligand_editor_canvas_append_molecule(widget_ptr, rdkit_mol);
        #endif
        
        ctx.widget_data.update_status("New molecule created from carbon ring.");
        // todo: make sure that this is crash-safe vs edit/undo
    } else {
        g_debug("There are already molecules. Structure insertion on blank space is a no-op.");
    }
}

bool RemoveHydrogensTool::on_molecule_click(MoleculeClickContext& ctx) {
    ctx.widget_data.begin_edition();
    layla::remove_non_polar_hydrogens(*ctx.rdkit_mol);
    this->sanitize_molecule(ctx.widget_data, *ctx.rdkit_mol);
    ctx.canvas_mol.lower_from_rdkit(!ctx.widget_data.allow_invalid_molecules);
    ctx.widget_data.finalize_edition();
    ctx.widget_data.update_status("Non-polar hydrogens have been removed.");
    return false;
}

void RemoveHydrogensTool::on_load(impl::WidgetCoreData& widget_data) {
    if(widget_data.get_molecule_count_impl() == 1) {
        auto idx_of_first = widget_data.get_first_molecule_idx();
        auto& canvas_mol = *widget_data.molecules->at(idx_of_first);
        auto& rdkit_mol = *widget_data.rdkit_molecules->at(idx_of_first);
        ClickContext ctx(widget_data);
        ctx.control_pressed = false;
        MoleculeClickContext mctx(ctx, 0, rdkit_mol, canvas_mol);
        this->on_molecule_click(mctx);
    }
}

std::string RemoveHydrogensTool::get_exception_message_prefix() const noexcept {
    return "Could not remove hydrogens from molecule: ";
}

void FlipTool::on_load(impl::WidgetCoreData& widget_data) {
    if(widget_data.get_molecule_count_impl() == 1) {
        auto idx_of_first = widget_data.get_first_molecule_idx();
        auto& canvas_mol = *widget_data.molecules->at(idx_of_first);
        auto& rdkit_mol = *widget_data.rdkit_molecules->at(idx_of_first);
        ClickContext ctx(widget_data);
        ctx.control_pressed = false;
        MoleculeClickContext mctx(ctx, 0, rdkit_mol, canvas_mol);
        this->on_molecule_click(mctx);
    }
}

bool FlipTool::on_molecule_click(MoleculeClickContext& ctx) {
    ctx.widget_data.begin_edition();
    ctx.canvas_mol.perform_flip(this->mode);
    ctx.canvas_mol.lower_from_rdkit(!ctx.widget_data.allow_invalid_molecules);
    ctx.widget_data.finalize_edition();
    ctx.widget_data.update_status("Molecule has been flipped.");
    return false;
}

std::string FlipTool::get_exception_message_prefix() const noexcept {
    return "Could not flip molecule: ";
}

void FormatTool::on_load(impl::WidgetCoreData& widget_data) {
    if(widget_data.get_molecule_count_impl() == 1) {
        auto idx_of_first = widget_data.get_first_molecule_idx();
        auto& canvas_mol = *widget_data.molecules->at(idx_of_first);
        auto& rdkit_mol = *widget_data.rdkit_molecules->at(idx_of_first);
        ClickContext ctx(widget_data);
        ctx.control_pressed = false;
        MoleculeClickContext mctx(ctx, 0, rdkit_mol, canvas_mol);
        this->on_molecule_click(mctx);
    }
}

bool FormatTool::on_molecule_click(MoleculeClickContext& ctx) {
    ctx.widget_data.begin_edition();
    ctx.canvas_mol.clear_cached_atom_coordinate_map();
    ctx.canvas_mol.lower_from_rdkit(!ctx.widget_data.allow_invalid_molecules, false);
    ctx.widget_data.finalize_edition();
    ctx.widget_data.update_status("Molecule has been formatted.");
    return false;
}

std::string FormatTool::get_exception_message_prefix() const noexcept {
    return "Could not format molecule: ";
}

void ActiveTool::begin_transform(int x, int y, TransformManager::Mode mode) {
    auto mol_opt = this->widget_data->resolve_click(x, y);
    if(mol_opt.has_value()) {
        auto [atom_or_bond,mol_id] = mol_opt.value();
        transform_manager.begin_transform(x, y, mode);
        transform_manager.set_canvas_molecule_index(mol_id);
        this->widget_data->begin_edition();
    }
}

bool ActiveTool::is_in_transform() const noexcept {
    return this->transform_manager.is_active();
}

void ActiveTool::update_transform_cursor_pos(int x, int y, bool snap_to_angle) noexcept {
    if(!transform_manager.is_active()) {
        return;
    }
    transform_manager.update_current_cursor_pos(x, y, snap_to_angle);
    transform_manager.apply_current_transform_state(this->widget_data, snap_to_angle, false);
}

void ActiveTool::end_transform(bool snap_to_angle) {
    if(!transform_manager.is_active()) {
        return;
    }
    transform_manager.apply_current_transform_state(this->widget_data, snap_to_angle, true);
    transform_manager.end_transform();
    this->widget_data->finalize_edition();
}

void Tool::sanitize_molecule(impl::WidgetCoreData& widget_data, RDKit::RWMol& mol) {
    if (!widget_data.allow_invalid_molecules) {
        RDKit::MolOps::sanitizeMol(mol);
    }
}

std::optional<std::pair<unsigned int,unsigned int>> ActiveTool::get_molecule_idx_and_first_atom_of_new_bond() const noexcept {
    const BondModifier* mod = dynamic_cast<const BondModifier*>(this->tool.get());
    if(mod) {
        return mod->get_molecule_idx_and_first_atom_of_new_bond();
    } else {
        return std::nullopt;
    }
}

ElementInsertion::ElementInsertion(ElementInsertion::Element el) noexcept {
    this->element = el;
}

ElementInsertion::ElementInsertion(const char* symbol) {
    RDKit::PeriodicTable* t = RDKit::PeriodicTable::getTable();
    unsigned int atomic_number = t->getAtomicNumber(symbol);
    this->element = atomic_number;
}

unsigned int ElementInsertion::get_atomic_number() const noexcept {
    if(std::holds_alternative<Element>(this->element)) {
        switch (std::get<Element>(this->element)) {
            default:
            case ElementInsertion::Element::C:{
                return 6;
            }
            case ElementInsertion::Element::N:{
                return 7;
            }
            case ElementInsertion::Element::O:{
                return 8;
            }
            case ElementInsertion::Element::S:{
                return 16;
            }
            case ElementInsertion::Element::P:{
                return 15;
            }
            case ElementInsertion::Element::H:{
                return 1;
            }
            case ElementInsertion::Element::F:{
                return 9;
            }
            case ElementInsertion::Element::Cl:{
                return 17;
            }
            case ElementInsertion::Element::Br:{
                return 35;
            }
            case ElementInsertion::Element::I:{
                return 53;
            }
        }
    } else { // raw atomic number
        return std::get<unsigned int>(this->element);
    }
}

FlipTool::FlipTool(FlipMode mode) noexcept {
    this->mode = mode;
}
