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
        // this is done elsewhere
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

        auto* widget_ptr = static_cast<const CootLigandEditorCanvasPriv*>(this);
        gtk_widget_queue_draw(GTK_WIDGET(widget_ptr));
    }
}