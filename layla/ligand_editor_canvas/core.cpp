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
    this->molecules = std::make_unique<std::vector<std::optional<CanvasMolecule>>>(*core_data.molecules);
    const std::vector<std::optional<std::shared_ptr<RDKit::RWMol>>>& original_rdkit_molecules = *core_data.rdkit_molecules;
    std::vector<std::optional<std::shared_ptr<RDKit::RWMol>>> copied_rdkit_molecules;
    for(const auto& rwmol_shptr_opt: original_rdkit_molecules) {
        if(rwmol_shptr_opt.has_value()) {
            const auto& rwmol = *rwmol_shptr_opt->get();
            copied_rdkit_molecules.push_back(std::make_shared<RDKit::RWMol>(RDKit::RWMol(rwmol)));
        } else {
            copied_rdkit_molecules.push_back(std::nullopt);
        }
    }
    for(unsigned int i = 0;i < this->molecules->size();i++) {
        auto& mol_opt = this->molecules->at(i);
        if(mol_opt.has_value()) {
            mol_opt->update_source_molecule(*copied_rdkit_molecules.at(i));
        }
    }
    this->rdkit_molecules = std::make_unique<std::vector<std::optional<std::shared_ptr<RDKit::RWMol>>>>(std::move(copied_rdkit_molecules));
}

#ifdef __EMSCRIPTEN__
const unsigned int WidgetCoreData::MAX_STATE_STACK_LENGTH = 20;
const unsigned int WidgetCoreData::STATE_STACK_TRIM_BATCH_SIZE = 5;
#else
const unsigned int WidgetCoreData::MAX_STATE_STACK_LENGTH = 100;
const unsigned int WidgetCoreData::STATE_STACK_TRIM_BATCH_SIZE = 30;
#endif
const unsigned int WidgetCoreData::VIEWPORT_OFFSET_PIXEL_MARGIN = 30;

WidgetCoreData::MaybeAtomOrBondWithMolIdx WidgetCoreData::resolve_click(int raw_x, int raw_y) const noexcept {
    const int x = raw_x + this->viewport_origin_offset.first;
    const int y = raw_y + this->viewport_origin_offset.second;
    const auto* molecules_vec = this->molecules.get();
    unsigned int idx = 0;
    for(const auto& mol_opt: *molecules_vec) {
        if(mol_opt.has_value()) {
            auto result = mol_opt->resolve_click(x, y, this->scale);
            if(result.has_value()) {
                return std::pair(result.value(),idx);
            }
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

void WidgetCoreData::emit_mutation_signals() const noexcept {
    #ifndef __EMSCRIPTEN__
    auto* widget_ptr = static_cast<const CootLigandEditorCanvasPriv*>(this);
    #else
    auto* widget_ptr = static_cast<const ::CootLigandEditorCanvas*>(this);
    #endif
    _LIGAND_EDITOR_SIGNAL_EMIT(widget_ptr, smiles_changed_signal);
    const auto* molecules_vec = this->molecules.get();
    for(unsigned int id = 0; id < molecules_vec->size(); id++) {
        const auto& mol_opt = molecules_vec->at(id);
        if(mol_opt.has_value()) {
            auto qed_info = mol_opt->get_qed_info();
            if(qed_info.has_value()) {
                const auto& qed_info_value = *qed_info;
                _LIGAND_EDITOR_SIGNAL_EMIT_ARG(widget_ptr, qed_info_updated_signal, id, &qed_info_value);
            }
        }
    }
}


void WidgetCoreData::undo_edition() {
    int state_stack_size = this->state_stack->size();
    if(state_stack_size > this->state_stack_pos + 1) {
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
        this->molecules = std::make_unique<std::vector<std::optional<CanvasMolecule>>>(*target_state->molecules);
        this->rdkit_molecules = std::make_unique<std::vector<std::optional<std::shared_ptr<RDKit::RWMol>>>>(*target_state->rdkit_molecules);
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
        this->molecules = std::make_unique<std::vector<std::optional<CanvasMolecule>>>(*target_state->molecules);
        this->rdkit_molecules = std::make_unique<std::vector<std::optional<std::shared_ptr<RDKit::RWMol>>>>(*target_state->rdkit_molecules);
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
        this->emit_mutation_signals();
    }
}

void WidgetCoreData::delete_molecule_with_idx(unsigned int idx, bool integrate_with_edit_undo) noexcept {
    if(idx >= this->rdkit_molecules->size()) {
        return;
    }
    auto& mol = this->molecules->at(idx);
    if(mol.has_value()) {
        if(integrate_with_edit_undo) {
            this->begin_edition();
        }
        mol = std::nullopt;
        this->rdkit_molecules->at(idx) = std::nullopt;

        if(integrate_with_edit_undo) {
            this->finalize_edition();
            this->update_status("Molecule deleted.");
        }

        // Already called by finalize_edition()
        // this->queue_redraw();

        #ifndef __EMSCRIPTEN__
        auto* widget_ptr = static_cast<const CootLigandEditorCanvasPriv*>(this);
        #else
        auto* widget_ptr = static_cast<const ::CootLigandEditorCanvas*>(this);
        #endif
        _LIGAND_EDITOR_SIGNAL_EMIT_ARG(widget_ptr, molecule_deleted_signal, idx);
    }
}

coot::ligand_editor_canvas::SmilesMap WidgetCoreData::build_smiles() const {
    coot::ligand_editor_canvas::SmilesMap ret;
    auto it = this->rdkit_molecules->cbegin();
    unsigned int idx = 0;

    auto append_smiles = [&](){
        const auto& mol_ptr_opt = *it;
        if(mol_ptr_opt.has_value()) {
            RDKit::RWMol* mol_ptr = mol_ptr_opt->get();
            ret.emplace(idx, RDKit::MolToSmiles(*mol_ptr));
        }
    };

    while(it != this->rdkit_molecules->cend()) {
        append_smiles();
        it++;
        idx++;
    }
    return ret;
}

coot::ligand_editor_canvas::InchiKeyMap WidgetCoreData::build_inchi_keys() const {
    coot::ligand_editor_canvas::InchiKeyMap ret;
    auto it = this->rdkit_molecules->cbegin();
    unsigned int idx = 0;

    auto append_inchi_key = [&](){
        const auto& mol_ptr_opt = *it;
        if(mol_ptr_opt.has_value()) {
            RDKit::RWMol* mol_ptr = mol_ptr_opt->get();
            #ifdef RDK_BUILD_INCHI_SUPPORT
            ret.emplace(idx, RDKit::MolToInchiKey(*mol_ptr));
            #else
            ret.emplace(idx, "");
            #warning Your version of RDKit was built without InChI support. Molecule InChI key lookup will not be available.
            #endif
        }
    };

    while(it != this->rdkit_molecules->cend()) {
        append_inchi_key();
        it++;
        idx++;
    }
    return ret;
}

unsigned int WidgetCoreData::get_molecule_count_impl() const noexcept {
    unsigned int ret = 0;
    for(const auto& mol_opt: *this->rdkit_molecules) {
        if(mol_opt.has_value()) {
            ret += 1;
        }
    }
    return ret;
}

int WidgetCoreData::get_first_molecule_idx() const noexcept {
    int idx = 0;
    for(const auto& mol_opt: *this->rdkit_molecules) {
        if(mol_opt.has_value()) {
            return idx;
        }
        idx += 1;
    }
    return -1;
}

void WidgetCoreData::render(Renderer& ren) const {
    if (this->molecules) {
        for(auto& drawn_molecule_opt: *this->molecules) {
            if(drawn_molecule_opt.has_value()) {
                drawn_molecule_opt->draw(ren, this->display_mode, this->viewport_origin_offset, this->scale);
            }
        }
    } else {
        g_error("Molecules vector not initialized!");
    }
    if(this->currently_created_bond.has_value()) {
        auto& bond = *this->currently_created_bond;
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

graphene_rect_t WidgetCoreData::get_on_screen_bounding_rect() const noexcept {
    graphene_rect_t bounding_rect_for_all;
    graphene_rect_init(&bounding_rect_for_all, 0, 0, 0, 0);

    for(const auto& a: *this->molecules) {
        if(a.has_value()) {
            auto bounding_rect = a->get_on_screen_bounding_rect(this->viewport_origin_offset, this->scale);
            graphene_rect_union(&bounding_rect_for_all, &bounding_rect, &bounding_rect_for_all);
        }
    }
    return bounding_rect_for_all;
}

void WidgetCoreData::queue_resize() noexcept {
    const auto bounding_rect_for_all = this->get_on_screen_bounding_rect();
    std::pair<int, int> viewport_offset = {0, 0};
    if(bounding_rect_for_all.origin.x < 0) {
        viewport_offset.first = bounding_rect_for_all.origin.x - impl::WidgetCoreData::VIEWPORT_OFFSET_PIXEL_MARGIN;
        g_debug("Viewport x-offset set to %i", viewport_offset.first);
    }
    if(bounding_rect_for_all.origin.y < 0) {
        viewport_offset.second = bounding_rect_for_all.origin.y - impl::WidgetCoreData::VIEWPORT_OFFSET_PIXEL_MARGIN;
        g_debug("Viewport y-offset set to %i", viewport_offset.second);
    }
    this->viewport_origin_offset = viewport_offset;
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

int CootLigandEditorCanvas::append_molecule(std::shared_ptr<RDKit::RWMol> rdkit_mol) noexcept {
    return coot_ligand_editor_canvas_append_molecule(this, std::move(rdkit_mol));
}

 void CootLigandEditorCanvas::update_molecule_from_smiles(unsigned int molecule_idx, const std::string& smiles) {
    coot_ligand_editor_canvas_update_molecule_from_smiles(this, molecule_idx, smiles.c_str());
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

unsigned int CootLigandEditorCanvas::get_idx_of_first_molecule() noexcept {
    return coot_ligand_editor_canvas_get_idx_of_first_molecule(this);
}

unsigned int CootLigandEditorCanvas::get_max_molecule_idx() noexcept {
    return coot_ligand_editor_canvas_get_max_molecule_idx(this);
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

coot::ligand_editor_canvas::InchiKeyMap CootLigandEditorCanvas::get_inchi_keys() noexcept {
    return coot_ligand_editor_canvas_get_inchi_keys(this);
}

coot::ligand_editor_canvas::SmilesMap CootLigandEditorCanvas::get_smiles() noexcept {
    return coot_ligand_editor_canvas_get_smiles(this);
}

std::string CootLigandEditorCanvas::get_smiles_for_molecule(unsigned int molecule_idx) noexcept {
    return coot_ligand_editor_canvas_get_smiles_for_molecule(this, molecule_idx);
}

std::string CootLigandEditorCanvas::get_inchi_key_for_molecule(unsigned int molecule_idx) noexcept {
    return coot_ligand_editor_canvas_get_inchi_key_for_molecule(this, molecule_idx);
}

std::string CootLigandEditorCanvas::get_pickled_molecule(unsigned int molecule_idx) noexcept {
    return coot_ligand_editor_canvas_get_pickled_molecule(this, molecule_idx);
}

std::string CootLigandEditorCanvas::get_pickled_molecule_base64(unsigned int molecule_idx) noexcept {
    return coot_ligand_editor_canvas_get_pickled_molecule_base64(this, molecule_idx);
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
    } else if(signal_name == "qed_info_updated") {
        qed_info_updated_signal.connect([=](int molecule_idx, const coot::ligand_editor_canvas::CanvasMolecule::QEDInfo* qed_info){
            callback(molecule_idx, *qed_info);
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
