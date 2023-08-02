#include <gtk/gtk.h>
#include "ligand_builder_state.hpp"

using namespace coot::ligand_editor;
using namespace coot::ligand_editor_canvas;
using BondModifierMode = coot::ligand_editor_canvas::BondModifier::BondModifierMode;
using Element = coot::ligand_editor_canvas::ElementInsertion::Element;
using Structure = coot::ligand_editor_canvas::StructureInsertion::Structure;
using TransformMode = coot::ligand_editor_canvas::TransformManager::Mode;

extern "C" G_MODULE_EXPORT
void
layla_on_close(GtkButton* button, gpointer user_data) {
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");

    state->file_exit();
}

extern "C" G_MODULE_EXPORT
void
layla_on_apply(GtkButton* button, gpointer user_data) {
    g_warning("TODO: Implement 'Apply'");
}

extern "C" G_MODULE_EXPORT
void
layla_on_invalid_molecule_toggled(GtkCheckButton* check_button, gpointer user_data) {
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    coot_ligand_editor_set_allow_invalid_molecules(canvas, gtk_check_button_get_active(check_button));
}

extern "C" G_MODULE_EXPORT
void
layla_on_show_alerts_toggled(GtkCheckButton* check_button, gpointer user_data) {
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    g_warning("TODO: Implement 'Show Alerts'");
}

extern "C" G_MODULE_EXPORT
void
layla_on_scale_spinbutton_value_changed(GtkSpinButton* self,gpointer user_data){
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    double new_scale = gtk_spin_button_get_value(self);
    // This should prevent infinite cascade of signals being emited
    if (coot_ligand_editor_get_scale(canvas) != new_scale) {
        coot_ligand_editor_set_scale(canvas, new_scale);
    }
}

extern "C" G_MODULE_EXPORT
void
layla_on_X_button_clicked(GtkButton* _btn, gpointer user_data){
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    state->run_choose_element_dialog();
}

extern "C" G_MODULE_EXPORT
void
layla_on_C_button_clicked(GtkButton* _btn, gpointer user_data){
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::C)));
}

extern "C" G_MODULE_EXPORT
void
layla_on_N_button_clicked(GtkButton* _btn, gpointer user_data){
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::N)));
}

extern "C" G_MODULE_EXPORT
void
layla_on_O_button_clicked(GtkButton* _btn, gpointer user_data){
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::O)));
}

extern "C" G_MODULE_EXPORT
void
layla_on_S_button_clicked(GtkButton* _btn, gpointer user_data){
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::S)));
}

extern "C" G_MODULE_EXPORT
void
layla_on_P_button_clicked(GtkButton* _btn, gpointer user_data){
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::P)));
}

extern "C" G_MODULE_EXPORT
void
layla_on_H_button_clicked(GtkButton* _btn, gpointer user_data){
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::H)));
}

extern "C" G_MODULE_EXPORT
void
layla_on_F_button_clicked(GtkButton* _btn, gpointer user_data){
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::F)));
}

extern "C" G_MODULE_EXPORT
void
layla_on_Cl_button_clicked(GtkButton* _btn, gpointer user_data){
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::Cl)));
}

extern "C" G_MODULE_EXPORT
void
layla_on_Br_button_clicked(GtkButton* _btn, gpointer user_data){
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::Br)));
}

extern "C" G_MODULE_EXPORT
void
layla_on_I_button_clicked(GtkButton* _btn, gpointer user_data){
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
 
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::I)));
}

extern "C" G_MODULE_EXPORT
void
layla_on_button_3C_clicked(GtkButton* _btn, gpointer user_data){
    CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(StructureInsertion(Structure::CycloPropaneRing)));

}

extern "C" G_MODULE_EXPORT
void
layla_on_button_4C_clicked(GtkButton* _btn, gpointer user_data){
    CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(StructureInsertion(Structure::CycloButaneRing)));

}

extern "C" G_MODULE_EXPORT
void
layla_on_button_5C_clicked(GtkButton* _btn, gpointer user_data){
    CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(StructureInsertion(Structure::CycloPentaneRing)));

}

extern "C" G_MODULE_EXPORT
void
layla_on_button_6C_clicked(GtkButton* _btn, gpointer user_data){
    CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(StructureInsertion(Structure::CycloHexaneRing)));

}

extern "C" G_MODULE_EXPORT
void
layla_on_button_6Arom_clicked(GtkButton* _btn, gpointer user_data){
    CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(StructureInsertion(Structure::BenzeneRing)));

}

extern "C" G_MODULE_EXPORT
void
layla_on_button_7C_clicked(GtkButton* _btn, gpointer user_data){
    CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(StructureInsertion(Structure::CycloHeptaneRing)));

}

extern "C" G_MODULE_EXPORT
void
layla_on_button_8C_clicked(GtkButton* _btn, gpointer user_data){
    CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
    coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(StructureInsertion(Structure::CycloOctaneRing)));


}