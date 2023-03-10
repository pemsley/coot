#include "core.hpp"
#include <iterator>

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

void WidgetCoreData::update_status(const gchar* status_text) const noexcept {
    auto* widget_ptr = static_cast<const CootLigandEditorCanvasPriv*>(this);
    g_signal_emit((gpointer) widget_ptr, status_updated_signal, 0, status_text);
}


void WidgetCoreData::undo_edition() {
    g_warning("EDIT->UNDO/REDO IS EXPERIMENTAL UNTESTED CODE. BUGS ARE EXPECTED");
    if(this->state_stack->size() > this->state_stack_pos) {
        auto iterator = this->state_stack->rbegin();
        if(this->state_stack_pos == 0) { // special case. the current state has to be snapshotted
            this->state_stack->push_back(std::make_unique<StateSnapshot>(*this));
            g_debug("Special case: stack pos is at 0. Current state has been snapshotted and appended to the stack.");
            iterator = this->state_stack->rbegin();
            iterator++;
        } else {
            std::advance(iterator,this->state_stack_pos);
            g_debug("Reverse iterator advanced by %u, thus expected to point at [last - %u] with stack size: %zu",this->state_stack_pos, this->state_stack_pos,this->state_stack->size());
        }
        auto& target_state = *iterator;
        this->molecules = std::make_unique<std::vector<CanvasMolecule>>(*target_state->molecules);
        this->rdkit_molecules = std::make_unique<std::vector<std::shared_ptr<RDKit::RWMol>>>(*target_state->rdkit_molecules);
        this->state_stack_pos++;
    } else {
        g_debug("Nothing to be undone.");
    }
}


void WidgetCoreData::redo_edition() {
    g_warning("EDIT->UNDO/REDO IS EXPERIMENTAL UNTESTED CODE. BUGS ARE EXPECTED");
    if(this->state_stack_pos > 0) {
        g_debug("Initial stack pos: %u",this->state_stack_pos);
        auto old_state_stack_pos = this->state_stack_pos;
        this->state_stack_pos--;
        g_debug("Target stack pos: %u",this->state_stack_pos);
        auto iterator = this->state_stack->rbegin();
        std::advance(iterator,this->state_stack_pos);
        g_debug("Reverse iterator advanced by %u, thus expected to point at [last - %u] with stack size: %zu",this->state_stack_pos, this->state_stack_pos,this->state_stack->size());
        auto& target_state = *iterator;
        this->molecules = std::make_unique<std::vector<CanvasMolecule>>(*target_state->molecules);
        this->rdkit_molecules = std::make_unique<std::vector<std::shared_ptr<RDKit::RWMol>>>(*target_state->rdkit_molecules);
        if(old_state_stack_pos == 1) { 
            // special case. If we're now at pos 0
            // the last thing from the stack needs to be removed
            this->state_stack->pop_back();
            g_debug("We're not at position 0. The newest element has been removed from stack to prevent duplication.");
        }
    } else {
        g_debug("Position in stack is 0 (the newest element). Nothing to redo. Stack size: %zu",this->state_stack->size());
    }
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
        if (this->state_stack_pos > 0) {
            g_warning("EDIT->UNDO/REDO IS EXPERIMENTAL UNTESTED CODE. BUGS ARE EXPECTED");
            auto& state_stack = *this->state_stack;
            auto it1 = state_stack.begin();
            g_debug("Finalizing edition with pos in stack: %u. The stack has to be trimmed.",this->state_stack_pos);
            g_debug("Current stack size: %zu Iterator advanced by: %zu",this->state_stack->size(),this->state_stack->size() - this->state_stack_pos);
            std::advance(it1,state_stack.size() - this->state_stack_pos);
            state_stack.erase(it1);
            this->state_stack_pos = 0;
            g_debug("Stack size after trim: %zu",this->state_stack->size());
        }
        this->state_stack->push_back(std::move(this->state_before_edition));
    }
}