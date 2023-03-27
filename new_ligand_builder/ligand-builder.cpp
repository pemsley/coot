#include "ligand-builder.hpp"
#include "gtk/gtktypebuiltins.h"
#include "ligand_editor_canvas.hpp"
#include <exception>
#include <memory>
#include <stdexcept>
#include <gtk/gtk.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/FileParsers/FileParsers.h>

using namespace coot::ligand_editor;
/// Structure holding the state of the editor

LigandBuilderState::LigandBuilderState(CootLigandEditorCanvas* canvas_widget, GtkWindow* win, GtkLabel* status_label) noexcept {
    this->canvas = canvas_widget;
    this->main_window = win;
    this->status_label = status_label;
}

void LigandBuilderState::update_status(const char* new_status) noexcept {
    if(this->status_label) {
        gtk_label_set_text(this->status_label, new_status);
    }
}

void LigandBuilderState::append_molecule(RDKit::RWMol* molecule_ptr) {
    coot_ligand_editor_append_molecule(this->canvas, std::shared_ptr<RDKit::RWMol>(molecule_ptr));
}

void LigandBuilderState::load_from_smiles() {
    auto* load_dialog = gtk_dialog_new();
    gtk_window_set_transient_for(GTK_WINDOW(load_dialog), this->main_window);
    // This isn't the best practice but it tremendously simplifies things
    // by saving us from unnecessary boilerplate.
    g_object_set_data(G_OBJECT(load_dialog), "ligand_builder_instance", this);
    gtk_window_set_title(GTK_WINDOW(load_dialog), "Load from SMILES");
    
    auto* dialog_body = gtk_box_new(GTK_ORIENTATION_VERTICAL,10);
    gtk_widget_set_margin_bottom(dialog_body, 10);
    gtk_widget_set_margin_end(dialog_body, 10);
    gtk_widget_set_margin_start(dialog_body, 10);
    gtk_widget_set_margin_top(dialog_body, 10);

    auto* label = gtk_label_new("Insert SMILES code");
    gtk_box_append(GTK_BOX(dialog_body),label);

    auto* entry_buf = gtk_entry_buffer_new("", 0);
    auto* entry = gtk_entry_new_with_buffer(entry_buf);

    gtk_box_append(GTK_BOX(dialog_body),entry);

    auto* submit_button = gtk_button_new_with_label("Submit");
    gtk_box_append(GTK_BOX(dialog_body),submit_button);
    g_signal_connect(submit_button, "clicked", G_CALLBACK(+[](GtkButton* btn, gpointer user_data){
        gtk_dialog_response(GTK_DIALOG(user_data), GTK_RESPONSE_ACCEPT);
    }), load_dialog);

    g_signal_connect(load_dialog, "response", G_CALLBACK(+[](GtkDialog* dialog, gint response_id, gpointer user_data){
        if(response_id != GTK_RESPONSE_ACCEPT) {
            g_debug("Ignoring unhandled response type: %s",g_enum_to_string(gtk_response_type_get_type(), response_id));
            return;
        }
        g_info("Importing SMILES...");
        auto* text_buf = GTK_ENTRY_BUFFER(user_data);
        try {

            RDKit::RWMol* molecule = RDKit::SmilesToMol(gtk_entry_buffer_get_text(text_buf));
            if(!molecule) {
                throw std::runtime_error("RDKit::RWMol* is a nullptr. The SMILES code is probably invalid.");
            }
            g_info("SMILES Import: Molecule constructed.");
            LigandBuilderState* state = (LigandBuilderState*) g_object_get_data(G_OBJECT(dialog), "ligand_builder_instance");
            state->append_molecule(molecule);
            gtk_window_destroy(GTK_WINDOW(dialog));
        } catch (std::exception& e) {
            g_warning("SMILES Import error: %s",e.what());
            auto* message = gtk_message_dialog_new(
                GTK_WINDOW(dialog), 
                GTK_DIALOG_DESTROY_WITH_PARENT, 
                GTK_MESSAGE_ERROR, 
                GTK_BUTTONS_CLOSE, 
                "Error: Molecule could not be constructed.\n%s", 
                e.what()
            );
            g_signal_connect(message,"response",G_CALLBACK(+[](GtkDialog* message_dialog, gint response_id, gpointer user_data){
                gtk_window_close(GTK_WINDOW(message_dialog));
            }),nullptr);
            gtk_widget_show(message);
        }
    }), entry_buf);

    gtk_window_set_child(GTK_WINDOW(load_dialog),dialog_body);
    gtk_window_present(GTK_WINDOW(load_dialog));
}

void LigandBuilderState::file_import_molecule() {
    g_warning("TODO: Implement void LigandBuilderState::file_import_molecule()");
}

void LigandBuilderState::file_fetch_molecule() {
    g_warning("TODO: Implement void LigandBuilderState::file_fetch_molecule()");
}

void LigandBuilderState::file_new() {
    g_warning("TODO: Implement void LigandBuilderState::file_new()");
}

void LigandBuilderState::file_save() {
    g_warning("TODO: Finish implementing void LigandBuilderState::file_save()");
    file_save_as();

}

void LigandBuilderState::file_save_as() {
    g_warning("TODO: Finish mplementing void LigandBuilderState::file_save_as()");
    auto* save_dialog = gtk_file_dialog_new();
    gtk_file_dialog_save(save_dialog, this->main_window, NULL, +[](GObject* source_object, GAsyncResult* res, gpointer user_data){
        GError** e = NULL;
        GFile* file = gtk_file_dialog_save_finish(GTK_FILE_DIALOG(source_object), res, e);
        LigandBuilderState* self = (LigandBuilderState*) user_data;
        if(file) {
            //g_info("I have a file");
            const char* path = g_file_get_path(file);
            try {
                /// todo: Add support for multiple molecules
                const auto* mol = coot_ligand_editor_get_rdkit_molecule(self->canvas, 0);
                RDKit::MolToMolFile(*mol,std::string(path));
                g_info("MolFile Save: Molecule file saved.");
                self->update_status("File saved.");
            } catch(std::exception& e) {
                g_warning("MolFile Save error: %s",e.what());
                auto* message = gtk_message_dialog_new(
                    GTK_WINDOW(source_object), 
                    GTK_DIALOG_DESTROY_WITH_PARENT, 
                    GTK_MESSAGE_ERROR, 
                    GTK_BUTTONS_CLOSE, 
                    "Error: Molecule could not be saved to file.\n%s", 
                    e.what()
                );
            }
            g_object_unref(file);
        }
        if(e) {
            g_info("Save File: No file was given.");
            g_object_unref(*e);
        }
    }, this);
}

void LigandBuilderState::file_open() {
    auto* open_dialog = gtk_file_dialog_new();
    gtk_file_dialog_open(open_dialog, this->main_window, NULL, +[](GObject* source_object, GAsyncResult* res, gpointer user_data){
        GError** e = NULL;
        GFile* file = gtk_file_dialog_open_finish(GTK_FILE_DIALOG(source_object), res, e);
        LigandBuilderState* self = (LigandBuilderState*) user_data;
        if(file) {
            //g_info("I have a file");
            const char* path = g_file_get_path(file);
            try {
                RDKit::RWMol* mol = RDKit::MolFileToMol(std::string(path),true,false,false);
                if(!mol) {
                    throw std::runtime_error("RDKit::RWMol* is a nullptr. The MolFile could not be loaded.");
                }
                g_info("MolFile Import: Molecule constructed.");
                self->append_molecule(mol);
            } catch(std::exception& e) {
                g_warning("MolFile Import error: %s",e.what());
                auto* message = gtk_message_dialog_new(
                    GTK_WINDOW(source_object), 
                    GTK_DIALOG_DESTROY_WITH_PARENT, 
                    GTK_MESSAGE_ERROR, 
                    GTK_BUTTONS_CLOSE, 
                    "Error: Molecule could not be loaded.\n%s", 
                    e.what()
                );
            }
            g_object_unref(file);
        }
        if(e) {
            g_info("Open File: No file was given.");
            g_object_unref(*e);
        }
    }, this);

}

void LigandBuilderState::file_export(ExportMode mode) {
    g_warning("TODO: Implement exports.");
    switch (mode) {
        case ExportMode::PDF: {
            break;
        }
        case ExportMode::PNG: {
            break;
        }
        case ExportMode::SVG: {
            break;
        }
        default: {
            break;
        }
    }
}

void LigandBuilderState::edit_undo() {
    coot_ligand_editor_undo_edition(this->canvas);
}

void LigandBuilderState::edit_redo() {
    coot_ligand_editor_redo_edition(this->canvas);
}

void coot::ligand_editor::initialize_global_instance(CootLigandEditorCanvas* canvas, GtkWindow* win, GtkLabel* status_label) {
    global_instance = new LigandBuilderState(canvas,win,status_label);
    g_info("Global instance of LigandBuilderState has been initialized at: %p",global_instance);
}