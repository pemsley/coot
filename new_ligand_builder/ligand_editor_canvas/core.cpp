#include "core.hpp"

using namespace coot::ligand_editor_canvas;
using namespace coot::ligand_editor_canvas::impl;

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
    
}


void WidgetCoreData::redo_edition() {
    
}


void WidgetCoreData::rollback_current_edition() {
    
}


void WidgetCoreData::begin_edition() {
    
}

void WidgetCoreData::finalize_edition() {
    
}