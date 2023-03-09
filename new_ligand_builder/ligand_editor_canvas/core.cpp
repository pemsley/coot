#include "core.hpp"
#include <iterator>

using namespace coot::ligand_editor_canvas;
using namespace coot::ligand_editor_canvas::impl;

StateSnapshot::StateSnapshot(const WidgetCoreData& core_data) {
    this->molecules = std::make_unique<std::vector<CanvasMolecule>>(*core_data.molecules);
    this->rdkit_molecules = std::make_unique<std::vector<std::shared_ptr<RDKit::RWMol>>>(*core_data.rdkit_molecules);
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
    auto iterator = this->state_stack->rbegin();
    std::advance(iterator,this->state_stack_pos);
    auto& target_state = *iterator;
    this->molecules = std::make_unique<std::vector<CanvasMolecule>>(*target_state->molecules);
    this->rdkit_molecules = std::make_unique<std::vector<std::shared_ptr<RDKit::RWMol>>>(*target_state->rdkit_molecules);
    this->state_stack_pos++;
}


void WidgetCoreData::redo_edition() {
    g_warning("EDIT->UNDO/REDO IS EXPERIMENTAL UNTESTED CODE. BUGS ARE EXPECTED");
    if(this->state_stack_pos > 0) {
        this->state_stack_pos--;
        auto iterator = this->state_stack->rbegin();
        std::advance(iterator,this->state_stack_pos);
        auto& target_state = *iterator;
        this->molecules = std::make_unique<std::vector<CanvasMolecule>>(*target_state->molecules);
        this->rdkit_molecules = std::make_unique<std::vector<std::shared_ptr<RDKit::RWMol>>>(*target_state->rdkit_molecules);
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
            std::advance(it1,state_stack.size() - this->state_stack_pos);
            state_stack.erase(it1);
            this->state_stack_pos = 0;
        }
        this->state_stack->push_back(std::move(this->state_before_edition));
    }
}