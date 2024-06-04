/* layla/ligand_editor_canvas/core.cpp
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

#include "core.hpp"
#include "../ligand_editor_canvas.hpp"
#include <iterator>
#ifdef __EMSCRIPTEN__
#include <emscripten/val.h>
#endif

using namespace coot::ligand_editor_canvas;
using namespace coot::ligand_editor_canvas::impl;

StateSnapshot::StateSnapshot(const WidgetCoreData& core_data) {
    this->molecules = std::make_unique<std::vector<CanvasMolecule>>(*core_data.molecules);
    const std::vector<std::shared_ptr<RDKit::RWMol>>& original_rdkit_molecules = *core_data.rdkit_molecules;
    std::vector<std::shared_ptr<RDKit::RWMol>> copied_rdkit_molecules;
    for(const auto& rwmol_shptr: original_rdkit_molecules) {
        const auto& rwmol = *rwmol_shptr.get();
        copied_rdkit_molecules.push_back(std::make_shared<RDKit::RWMol>(RDKit::RWMol(rwmol)));
    }
    for(unsigned int i = 0;i < this->molecules->size();i++) {
        auto& mol = this->molecules->at(i);
        mol.update_source_molecule(copied_rdkit_molecules.at(i));
    }
    this->rdkit_molecules = std::make_unique<std::vector<std::shared_ptr<RDKit::RWMol>>>(std::move(copied_rdkit_molecules));
}

const unsigned int WidgetCoreData::MAX_STATE_STACK_LENGTH = 100;
const unsigned int WidgetCoreData::STATE_STACK_TRIM_BATCH_SIZE = 30;

WidgetCoreData::MaybeAtomOrBondWithMolIdx WidgetCoreData::resolve_click(int x, int y) const noexcept {
    const auto* molecules_vec = this->molecules.get();
    unsigned int idx = 0;
    for(const auto& mol: *molecules_vec) {
        auto result = mol.resolve_click(x, y);
        if(result.has_value()) {
            return std::pair(result.value(),idx);
        }
        idx++;
    }
    return std::nullopt;
}

void WidgetCoreData::update_status(const char* status_text) const noexcept {
    #ifndef __EMSCRIPTEN__
    auto* widget_ptr = static_cast<const CootLigandEditorCanvasPriv*>(this);
    #else
    auto* widget_ptr = static_cast<const ::CootLigandEditorCanvas*>(this);
    #endif
    _LIGAND_EDITOR_SIGNAL_EMIT_ARG(widget_ptr, status_updated_signal, status_text);
}


void WidgetCoreData::undo_edition() {
    if(this->state_stack->size() > this->state_stack_pos + 1) {
        if(this->state_stack_pos == -1) {
            // Special case. 
            // We need to backup the current state
            // so that we can later go back to it via "Redo"
            this->state_stack->push_back(std::make_unique<StateSnapshot>(*this));
            this->state_stack_pos++;
        }
        this->state_stack_pos++;
        auto iterator = this->state_stack->rbegin();
        std::advance(iterator,this->state_stack_pos);
        auto& target_state = *iterator;
        this->molecules = std::make_unique<std::vector<CanvasMolecule>>(*target_state->molecules);
        this->rdkit_molecules = std::make_unique<std::vector<std::shared_ptr<RDKit::RWMol>>>(*target_state->rdkit_molecules);
        update_status("");
        // this is done elsewhere, namely in the caller of this function
        // auto* widget_ptr = static_cast<const CootLigandEditorCanvasPriv*>(this);
        // gtk_widget_queue_draw(GTK_WIDGET(widget_ptr));
    } else {
        //g_debug("Nothing to be undone. Stack size: %zu Stack pos: %i",this->state_stack->size(),this->state_stack_pos);
        update_status("Nothing to be undone.");
    }
}


void WidgetCoreData::redo_edition() {
    if(this->state_stack_pos == 0) { 
        g_error("Internal error: Undo/Redo stack position should never stay at 0.");
    }
    if(this->state_stack_pos != -1) {
        //g_debug("Initial stack pos: %u",this->state_stack_pos);
        this->state_stack_pos--;
        //g_debug("Target stack pos: %u",this->state_stack_pos);
        auto iterator = this->state_stack->rbegin();
        std::advance(iterator,this->state_stack_pos);
        
        auto& target_state = *iterator;
        this->molecules = std::make_unique<std::vector<CanvasMolecule>>(*target_state->molecules);
        this->rdkit_molecules = std::make_unique<std::vector<std::shared_ptr<RDKit::RWMol>>>(*target_state->rdkit_molecules);
        if(this->state_stack_pos == 0) { 
            // Special case. We're now at pos 0 which means 
            // that we've reached the newest state.
            // We should remove it from the state history stack
            // and set the position to -1 to denote that we're "fresh", at the newest change.

            // g_debug(
            //     "We're now at position 0. The newest state has been reached. "
            //     "Resetting position to -1 and removing the last element from the stack to prevent duplication."
            // );
            this->state_stack->pop_back();
            this->state_stack_pos = -1;
        }
        update_status("");
    } else {
        //g_debug("Position in stack is at -1 (fresh change). Nothing to redo. Stack size: %zu",this->state_stack->size());
        update_status("Nothing to be redone.");
    }
}

bool WidgetCoreData::is_in_edition() {
    return this->state_before_edition.get() != nullptr;
}

void WidgetCoreData::rollback_current_edition() {
    if(this->state_before_edition) {
        this->molecules = std::move(this->state_before_edition->molecules);
        this->rdkit_molecules = std::move(this->state_before_edition->rdkit_molecules);
        this->state_before_edition.reset(nullptr);
    }
}


void WidgetCoreData::begin_edition() {
    this->state_before_edition = std::make_unique<StateSnapshot>(*this);
}

void WidgetCoreData::finalize_edition() {
    if(this->state_before_edition) {
        if (this->state_stack_pos != -1) {
            auto& state_stack = *this->state_stack;
            auto it1 = state_stack.begin();
            //g_debug("Finalizing edition with pos in stack: %i. The stack has to be trimmed.",this->state_stack_pos);
            //g_debug("Current stack size: %zu Iterator advanced by: %zu",this->state_stack->size(),this->state_stack->size() - this->state_stack_pos - 1);
            std::advance(it1,state_stack.size() - this->state_stack_pos - 1);
            state_stack.erase(it1);
            this->state_stack_pos = -1;
            //g_debug("Stack size after trim: %zu",this->state_stack->size());
        }
        this->state_stack->push_back(std::move(this->state_before_edition));

        if(this->state_stack->size() > MAX_STATE_STACK_LENGTH) {
            auto last_iter = state_stack->begin();
            std::advance(last_iter, STATE_STACK_TRIM_BATCH_SIZE);
            this->state_stack->erase(this->state_stack->begin(), last_iter);
        }

        this->queue_resize();
        this->queue_redraw();
        #ifndef __EMSCRIPTEN__
        auto* widget_ptr = static_cast<const CootLigandEditorCanvasPriv*>(this);
        #else
        auto* widget_ptr = static_cast<const ::CootLigandEditorCanvas*>(this);
        #endif
        _LIGAND_EDITOR_SIGNAL_EMIT(widget_ptr, smiles_changed_signal);
    }
}

void WidgetCoreData::delete_molecule_with_idx(unsigned int idx) noexcept {
    if(idx < this->rdkit_molecules->size()) {
        this->begin_edition();
        auto iter = this->molecules->begin();
        std::advance(iter, idx);
        this->molecules->erase(iter);

        auto iter2 = this->rdkit_molecules->begin();
        std::advance(iter2, idx);
        this->rdkit_molecules->erase(iter2);

        this->finalize_edition();
        this->update_status("Molecule deleted.");

        this->queue_redraw();

        #ifndef __EMSCRIPTEN__
        auto* widget_ptr = static_cast<const CootLigandEditorCanvasPriv*>(this);
        #else
        auto* widget_ptr = static_cast<const ::CootLigandEditorCanvas*>(this);
        #endif
        _LIGAND_EDITOR_SIGNAL_EMIT_ARG(widget_ptr, molecule_deleted_signal, idx);
    }
}

std::string WidgetCoreData::build_smiles_string() const {
    std::string ret;
    auto it = this->rdkit_molecules->cbegin();
    auto append_smiles = [&](){
        RDKit::RWMol* mol_ptr = it->get();
        ret += RDKit::MolToSmiles(*mol_ptr);
        
    };
    if(it != this->rdkit_molecules->cend()) {
        append_smiles();
        it++;
    }
    while(it != this->rdkit_molecules->cend()) {
        ret.push_back('\n');
        append_smiles();
        it++;
    }
    return ret;
}

void WidgetCoreData::render(Renderer& ren) {
    if (this->molecules) {
        for(auto& drawn_molecule: *this->molecules) {
            drawn_molecule.set_canvas_scale(this->scale);
            drawn_molecule.draw(ren,this->display_mode);
        }
    } else {
        g_error("Molecules vector not initialized!");
    }
    if(this->currently_created_bond.has_value()) {
        auto& bond = this->currently_created_bond.value();
        ren.set_line_width(4.0);
        ren.set_source_rgb(1.0, 0.5, 1.0);
        ren.move_to(bond.first_atom_x, bond.first_atom_y);
        ren.line_to(bond.second_atom_x, bond.second_atom_y);
        ren.stroke();
    }
}

void WidgetCoreData::queue_redraw() const noexcept {
    #ifndef __EMSCRIPTEN__
    auto* widget_ptr = static_cast<const CootLigandEditorCanvasPriv*>(this);
    gtk_widget_queue_draw(GTK_WIDGET(widget_ptr));
    #else
    auto* widget_ptr = static_cast<const::CootLigandEditorCanvas*>(this);
    _LIGAND_EDITOR_SIGNAL_EMIT(widget_ptr, queue_redraw_signal);
    #endif
}

void WidgetCoreData::queue_resize() const noexcept {
    #ifndef __EMSCRIPTEN__
    auto* widget_ptr = static_cast<const CootLigandEditorCanvasPriv*>(this);
    gtk_widget_queue_resize(GTK_WIDGET(widget_ptr));
    #else
    auto* widget_ptr = static_cast<const::CootLigandEditorCanvas*>(this);
    _LIGAND_EDITOR_SIGNAL_EMIT(widget_ptr, queue_resize_signal);
    #endif
}

#ifdef __EMSCRIPTEN__

void CootLigandEditorCanvas::set_active_tool(std::unique_ptr<coot::ligand_editor_canvas::ActiveTool> active_tool) {
    coot_ligand_editor_canvas_set_active_tool(this, std::move(active_tool));
}

void CootLigandEditorCanvas::append_molecule(std::shared_ptr<RDKit::RWMol> rdkit_mol) noexcept {
    coot_ligand_editor_canvas_append_molecule(this, std::move(rdkit_mol));
}

void CootLigandEditorCanvas::set_scale(float scale) noexcept {
    coot_ligand_editor_canvas_set_scale(this, scale);
}

float CootLigandEditorCanvas::get_scale() noexcept {
    return coot_ligand_editor_canvas_get_scale(this);
}

void CootLigandEditorCanvas::undo() noexcept {
    coot_ligand_editor_canvas_undo_edition(this);
}

void CootLigandEditorCanvas::redo() noexcept {
    coot_ligand_editor_canvas_redo_edition(this);
}

// const RDKit::ROMol& CanvasMolecule::get_rdkit_molecule(unsigned int index) noexcept {
//     RDKit::ROMol* coot_ligand_editor_canvas_get_rdkit_molecule(CootLigandEditorCanvas* self, unsigned int index) noexcept;
// }

unsigned int CootLigandEditorCanvas::get_molecule_count() noexcept {
    return coot_ligand_editor_canvas_get_molecule_count(this);
}

void CootLigandEditorCanvas::set_allow_invalid_molecules(bool value) noexcept {
    coot_ligand_editor_canvas_set_allow_invalid_molecules(this, value);
}

bool CootLigandEditorCanvas::get_allow_invalid_molecules() noexcept {
    return coot_ligand_editor_canvas_get_allow_invalid_molecules(this);
}

coot::ligand_editor_canvas::DisplayMode CootLigandEditorCanvas::get_display_mode() noexcept {
    return coot_ligand_editor_canvas_get_display_mode(this);
}

void CootLigandEditorCanvas::set_display_mode(coot::ligand_editor_canvas::DisplayMode value) noexcept {
    coot_ligand_editor_canvas_set_display_mode(this, value);
}

std::string CootLigandEditorCanvas::get_smiles() noexcept {
    return coot_ligand_editor_canvas_get_smiles(this);
}

std::string CootLigandEditorCanvas::get_smiles_for_molecule(unsigned int molecule_idx) noexcept {
    return coot_ligand_editor_canvas_get_smiles_for_molecule(this, molecule_idx);
}

std::string CootLigandEditorCanvas::get_pickled_molecule(unsigned int molecule_idx) noexcept {
    return coot_ligand_editor_canvas_get_pickled_molecule(this, molecule_idx);
}

void CootLigandEditorCanvas::clear_molecules() noexcept {
    coot_ligand_editor_canvas_clear_molecules(this);
}

void CootLigandEditorCanvas::connect(std::string signal_name, emscripten::val callback) {

    if(signal_name == "status_updated") {
        status_updated_signal.connect([=](const char* new_status){
            callback(std::string(new_status));
        });
    } else if(signal_name == "scale_changed") {
        scale_changed_signal.connect([=](float new_scale){
            callback(new_scale);
        });
    } else if(signal_name == "smiles_changed") {
        smiles_changed_signal.connect([=](){
            callback();
        });
    } else if(signal_name == "molecule_deleted") {
        molecule_deleted_signal.connect([=](int molecule_idx){
            callback(molecule_idx);
        });
    } else if(signal_name == "queue_redraw") {
        queue_redraw_signal.connect([=](){
            callback();
        });
    } else if(signal_name == "queue_resize") {
        queue_resize_signal.connect([=](){
            callback();
        });
    } else {
        g_critical("No such signal!");
    }
}
#endif // EMSCRIPTEN