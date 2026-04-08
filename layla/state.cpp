/* layla/state.cpp
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

#include "state.hpp"
#include "geometry/protein-geometry.hh"
#include "ligand_editor_canvas.hpp"
#include <exception>
#include <memory>
#include <optional>
#include <stdexcept>
#include <cairo.h>
#include <cairo/cairo-pdf.h>
#include <cairo/cairo-svg.h>
#include <gtk/gtk.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/FileParsers/FileParsers.h>
#include "utils.hpp"
// todo: remove dependency on lidia-core
#include "lidia-core/rdkit-interface.hh"
#include "utils/coot-utils.hh"
#include "inchi_key_database.hpp"
#include <string>
#include "python_utils.hpp"

using namespace coot::layla;

LaylaState::LaylaState(CootLigandEditorCanvas* canvas_widget, GtkWindow* win, GtkLabel* status_label) noexcept {
    this->canvas = canvas_widget;
    this->main_window = win;
    this->status_label = status_label;
    this->monomer_library_info_store = std::make_unique<protein_geometry>();
    this->unsaved_changes = false;
    this->notifier = coot_layla_notifier_new();
    g_signal_connect(canvas_widget, "molecule-deleted", G_CALLBACK(+[](CootLigandEditorCanvas* self, unsigned int deleted_mol_idx, gpointer user_data){
        LaylaState* state = (LaylaState*) user_data;
        if (state->current_filesave_molecule.has_value()) {
            unsigned int idx = state->current_filesave_molecule.value();
            if(idx == deleted_mol_idx) {
                // The "current" molecule no longer exists
                state->current_filesave_molecule = std::nullopt;
                state->current_filesave_filename = std::nullopt;
            }
        }
    }), this);


    std::string package_dir = coot::package_data_dir();
    typedef std::tuple<std::string, std::unique_ptr<InchiKeyDatabase>*> HelperTuple1;
    HelperTuple1* helper_tuple = new std::tuple(coot::util::append_dir_dir(package_dir, "data"), &this->inchi_key_database);
    g_idle_add(+[](gpointer user_data){
        HelperTuple1* helper_tuple = (HelperTuple1*) user_data;
        std::string inchi_key_database_path = coot::util::append_dir_file(std::get<0>(*helper_tuple), "Components-inchikey.ich");
        g_info("Loading InChIKey database from: %s", inchi_key_database_path.c_str());
        GFile* file = g_file_new_for_path(inchi_key_database_path.c_str());
        std::unique_ptr<InchiKeyDatabase>* load_target = std::get<1>(*helper_tuple);

        typedef std::tuple<std::unique_ptr<InchiKeyDatabase>*> HelperTuple2;
        HelperTuple2* ht2 = new std::tuple(load_target);
        g_file_load_contents_async(file, nullptr, +[](GObject* source_object, GAsyncResult* res, gpointer data) {
            HelperTuple2* ht2 = (HelperTuple2*) data;
            GError* error = nullptr;
            gchar* contents = nullptr;
            gsize length = 0;
            gboolean success = g_file_load_contents_finish(G_FILE(source_object), res, &contents, &length, nullptr, &error);
            if (success && contents) {
                try {
                    g_info("InchiKey database file loaded successfully. Size=%zu bytes.", length);
                    std::istringstream iss(std::string(contents, length));
                    std::unique_ptr<InchiKeyDatabase>* tgt = std::get<0>(*ht2);
                    tgt->reset(new InchiKeyDatabase(parseInchikeyDatabase(iss)));
                    g_info("InchiKey database parsed successfully. Entry count=%zu", (*tgt)->size());
                } catch (const std::exception& e) {
                    g_warning("Failed to load InChIKey database: %s", e.what());
                }
            } else {
                g_warning("Failed to read InChIKey database file: %s", error ? error->message : "unknown error");
                if (error) { 
                    g_error_free(error);
                }
            }
            if (contents) { 
                g_free(contents);
            }
            delete ht2;
            g_object_unref(G_FILE(source_object));
        }, ht2);
        delete helper_tuple;
        return gboolean(FALSE);
    }, helper_tuple);


    g_signal_connect(canvas_widget, "smiles-changed", G_CALLBACK(+[](CootLigandEditorCanvas* self, gpointer user_data){
        LaylaState* state = (LaylaState*) user_data;
        state->unsaved_changes = true;
    }), this);
    g_signal_connect(win, "close-request", G_CALLBACK(+[](GtkWindow* win, gpointer user_data){
        auto state = (LaylaState*) user_data;
        if(state->has_unsaved_changes()) {
            state->unsaved_changes_dialog_purpose = UnsavedChangesDialogPurpose::CloseEditor;
            auto* win = gtk_builder_get_object(global_layla_gtk_builder, "layla_unsaved_changes_dialog");
            gtk_window_present(GTK_WINDOW(win));
            return true;
        }
        // returning false closes the window
        return false;
    }), this);
}

LaylaState::~LaylaState() noexcept {
    if(this->notifier) {
        g_object_unref(this->notifier);
    }
}

void LaylaState::reset() noexcept {
    // for now this is sufficient.
    // Consider removing edit history too.
    this->current_filesave_filename = std::nullopt;
    this->current_filesave_filename = std::nullopt;
    this->unsaved_changes = false;
    coot_ligand_editor_canvas_clear_molecules(this->canvas);
    this->update_status("");
}

void LaylaState::unsaved_changes_dialog_accepted() {
    if(!this->unsaved_changes_dialog_purpose.has_value()) {
        throw std::runtime_error("unsaved_changes_dialog_accepted() called with no 'unsaved_changes_dialog_purpose' set.");
    }
    this->unsaved_changes = false;
    switch(this->unsaved_changes_dialog_purpose.value()) {
        case UnsavedChangesDialogPurpose::CloseEditor: {
            this->file_exit();
            break;
        }
        case UnsavedChangesDialogPurpose::FileNew: {
            this->file_new();
            break;
        }
        default: {
            throw std::runtime_error("Unknown value of UnsavedChangesDialogPurpose.");
            break;
        }
    }
}

bool LaylaState::has_unsaved_changes() const noexcept {
    return this->unsaved_changes;
}

CootLigandEditorCanvas* LaylaState::get_canvas() const noexcept {
    return this->canvas;
}

CootLaylaNotifier* LaylaState::get_notifier() const noexcept {
    return g_object_ref(this->notifier);
}

void LaylaState::update_status(const char* new_status) noexcept {
    if(this->status_label) {
        gtk_label_set_text(this->status_label, new_status);
    }
}

int LaylaState::append_molecule(RDKit::RWMol* molecule_ptr) {
    if(!coot_ligand_editor_canvas_get_allow_invalid_molecules(this->canvas)) {
        RDKit::MolOps::sanitizeMol(*molecule_ptr);
    }
    return coot_ligand_editor_canvas_append_molecule(this->canvas, std::shared_ptr<RDKit::RWMol>(molecule_ptr));
}

void LaylaState::load_from_smiles() {

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

    g_signal_connect(entry, "activate", G_CALLBACK(+[](GtkEntry* entry, gpointer user_data){
        gtk_dialog_response(GTK_DIALOG(user_data), GTK_RESPONSE_ACCEPT);
    }), load_dialog);

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
            RDKit::RWMol* molecule = RDKit::SmilesToMol(gtk_entry_buffer_get_text(text_buf), 0, false);
            if(!molecule) {
                throw std::runtime_error("RDKit::RWMol* is a nullptr. The SMILES code is probably invalid.");
            }
            // We don't need that here, do we?
            // RDKit::MolOps::sanitizeMol(*molecule);
            g_info("SMILES Import: Molecule constructed.");
            LaylaState* state = (LaylaState*) g_object_get_data(G_OBJECT(dialog), "ligand_builder_instance");
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

void LaylaState::file_import_molecule() {

    GtkWidget *load_dialog = gtk_dialog_new();
    gtk_window_set_transient_for(GTK_WINDOW(load_dialog), this->main_window);
    g_object_set_data(G_OBJECT(load_dialog), "ligand_builder_instance", this);
    gtk_window_set_title(GTK_WINDOW(load_dialog), "Monomer Import");
    GtkWidget *dialog_body = gtk_box_new(GTK_ORIENTATION_VERTICAL,10);
    gtk_widget_set_margin_bottom(dialog_body, 10);
    gtk_widget_set_margin_end(dialog_body, 10);
    gtk_widget_set_margin_start(dialog_body, 10);
    gtk_widget_set_margin_top(dialog_body, 10);

    GtkWidget *label = gtk_label_new("Insert Monomer Code");
    gtk_box_append(GTK_BOX(dialog_body),label);

    GtkEntryBuffer *entry_buf = gtk_entry_buffer_new("", 0);
    GtkWidget *entry = gtk_entry_new_with_buffer(entry_buf);

    gtk_box_append(GTK_BOX(dialog_body),entry);

    GtkWidget* remove_hydrogens_checkbutton = gtk_check_button_new_with_label("Remove hydrogens");
    gtk_box_append(GTK_BOX(dialog_body), remove_hydrogens_checkbutton);

    GtkWidget *submit_button = gtk_button_new_with_label("Submit");
    gtk_box_append(GTK_BOX(dialog_body), submit_button);

    auto submit_callback = +[] (GtkWidget* widget, gpointer user_data) {
        gtk_dialog_response(GTK_DIALOG(user_data), GTK_RESPONSE_ACCEPT);
    };
    g_signal_connect(submit_button, "clicked", G_CALLBACK(submit_callback), load_dialog);
    g_signal_connect(entry, "activate", G_CALLBACK(submit_callback), load_dialog);

    gtk_window_set_child(GTK_WINDOW(load_dialog), dialog_body);
    gtk_window_present(GTK_WINDOW(load_dialog));

    struct ImportDialogWidgets {
        GtkEntryBuffer* entry_buf;
        GtkCheckButton* remove_hydrogens_checkbutton;
    };

    ImportDialogWidgets* dialog_widgets = g_slice_new0(ImportDialogWidgets);
    dialog_widgets->entry_buf = entry_buf;
    dialog_widgets->remove_hydrogens_checkbutton = GTK_CHECK_BUTTON(remove_hydrogens_checkbutton);

    auto dialog_response = +[](GtkDialog* dialog, gint response_id, gpointer user_data) {
        if(response_id != GTK_RESPONSE_ACCEPT) {
            g_debug("Ignoring unhandled response type: %s", g_enum_to_string(gtk_response_type_get_type(), response_id));
            return;
        } else {
            ImportDialogWidgets* dialog_widgets = (ImportDialogWidgets*) user_data;
            const char *text_buf = gtk_entry_buffer_get_text(dialog_widgets->entry_buf);
            std::string monomer_type(text_buf);
            int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
            LaylaState* self = static_cast<LaylaState*>(g_object_get_data(G_OBJECT(dialog),
                                                                                        "ligand_builder_instance"));
            // what is 42???
            self->monomer_library_info_store->try_dynamic_add(monomer_type, 42);
            std::pair<bool, dictionary_residue_restraints_t> p =
                self->monomer_library_info_store->get_monomer_restraints(monomer_type, imol_enc);
            if (p.first) {
                bool should_remove_hydrogens = gtk_check_button_get_active(dialog_widgets->remove_hydrogens_checkbutton);
                // todo: it'd be best to rewrite this function.
                // It's a mess.
                auto mol = std::make_unique<RDKit::RWMol>(coot::rdkit_mol(p.second));
                if (should_remove_hydrogens) {
                    remove_non_polar_hydrogens(*mol);
                }
                RDKit::MolOps::sanitizeMol(*mol);
                int new_mol_id = self->append_molecule(mol.release());
                if(new_mol_id >= 0) {
                    self->current_filesave_molecule = static_cast<unsigned int>(new_mol_id);
                }
                gtk_window_destroy(GTK_WINDOW(dialog));
                g_slice_free(ImportDialogWidgets, user_data);
            } else {
                g_warning("Failed to find monomer \"%s\"", monomer_type.c_str());
                auto* message = gtk_message_dialog_new(
                    GTK_WINDOW(dialog), 
                    GTK_DIALOG_DESTROY_WITH_PARENT, 
                    GTK_MESSAGE_ERROR, 
                    GTK_BUTTONS_CLOSE, 
                    "Error: Monomer \"%s\" could not be found.\n", 
                    monomer_type.c_str()
                );
                g_signal_connect(message,"response",G_CALLBACK(+[](GtkDialog* message_dialog, gint response_id, gpointer user_data){
                    gtk_window_close(GTK_WINDOW(message_dialog));
                }),nullptr);
                gtk_widget_show(message);
            }
        }
   };

   g_signal_connect(load_dialog, "response", G_CALLBACK(dialog_response), dialog_widgets);
}

void LaylaState::run_choose_element_dialog() {
    auto* choose_element_dialog = gtk_dialog_new();
    gtk_window_set_transient_for(GTK_WINDOW(choose_element_dialog), this->main_window);
    // This isn't the best practice but it tremendously simplifies things
    // by saving us from unnecessary boilerplate.
    g_object_set_data(G_OBJECT(choose_element_dialog), "ligand_builder_instance", this);
    gtk_window_set_title(GTK_WINDOW(choose_element_dialog), "Pick chemical element");
    
    auto* dialog_body = gtk_box_new(GTK_ORIENTATION_VERTICAL,10);
    gtk_widget_set_margin_bottom(dialog_body, 10);
    gtk_widget_set_margin_end(dialog_body, 10);
    gtk_widget_set_margin_start(dialog_body, 10);
    gtk_widget_set_margin_top(dialog_body, 10);

    auto* label = gtk_label_new("Element symbol");
    gtk_box_append(GTK_BOX(dialog_body),label);

    auto* entry_buf = gtk_entry_buffer_new("", 0);
    auto* entry = gtk_entry_new_with_buffer(entry_buf);

    gtk_box_append(GTK_BOX(dialog_body),entry);

    auto* submit_button = gtk_button_new_with_label("Submit");
    gtk_box_append(GTK_BOX(dialog_body),submit_button);
    g_signal_connect(submit_button, "clicked", G_CALLBACK(+[](GtkButton* btn, gpointer user_data){
        gtk_dialog_response(GTK_DIALOG(user_data), GTK_RESPONSE_ACCEPT);
    }), choose_element_dialog);

    g_signal_connect(choose_element_dialog, "response", G_CALLBACK(+[](GtkDialog* dialog, gint response_id, gpointer user_data){
        if(response_id != GTK_RESPONSE_ACCEPT) {
            g_debug("Ignoring unhandled response type: %s",g_enum_to_string(gtk_response_type_get_type(), response_id));
            return;
        }
        auto* text_buf = GTK_ENTRY_BUFFER(user_data);
        try {
            auto insertion_tool = std::make_unique<ligand_editor_canvas::ActiveTool>(ligand_editor_canvas::ElementInsertion(gtk_entry_buffer_get_text(text_buf)));
            LaylaState* state = (LaylaState*) g_object_get_data(G_OBJECT(dialog), "ligand_builder_instance");
            coot_ligand_editor_canvas_set_active_tool(state->canvas, std::move(insertion_tool));
            gtk_window_destroy(GTK_WINDOW(dialog));
        } catch (std::exception& e) {
            g_warning("Could not pick element: %s",e.what());
            auto* message = gtk_message_dialog_new(
                GTK_WINDOW(dialog), 
                GTK_DIALOG_DESTROY_WITH_PARENT, 
                GTK_MESSAGE_ERROR, 
                GTK_BUTTONS_CLOSE, 
                "Error: Invalid symbol:\n%s", 
                e.what()
            );
            g_signal_connect(message,"response",G_CALLBACK(+[](GtkDialog* message_dialog, gint response_id, gpointer user_data){
                gtk_window_close(GTK_WINDOW(message_dialog));
            }),nullptr);
            gtk_widget_show(message);
        }
    }), entry_buf);

    gtk_window_set_child(GTK_WINDOW(choose_element_dialog),dialog_body);
    gtk_window_present(GTK_WINDOW(choose_element_dialog));
}

void LaylaState::file_fetch_molecule() {

    GtkWidget *load_dialog = gtk_dialog_new();
    gtk_window_set_transient_for(GTK_WINDOW(load_dialog), this->main_window);
    g_object_set_data(G_OBJECT(load_dialog), "ligand_builder_instance", this);
    gtk_window_set_title(GTK_WINDOW(load_dialog), "Fetch Molecule");
    GtkWidget *dialog_body = gtk_box_new(GTK_ORIENTATION_VERTICAL,10);
    gtk_widget_set_margin_bottom(dialog_body, 10);
    gtk_widget_set_margin_end(dialog_body, 10);
    gtk_widget_set_margin_start(dialog_body, 10);
    gtk_widget_set_margin_top(dialog_body, 10);

    GtkWidget *label = gtk_label_new("Type in drug name");
    gtk_box_append(GTK_BOX(dialog_body),label);

    GtkEntryBuffer *entry_buf = gtk_entry_buffer_new("", 0);
    GtkWidget *entry = gtk_entry_new_with_buffer(entry_buf);

    gtk_box_append(GTK_BOX(dialog_body),entry);

    GtkWidget *submit_button = gtk_button_new_with_label("Fetch");
    gtk_box_append(GTK_BOX(dialog_body), submit_button);

    auto submit_callback = +[] (GtkWidget* widget, gpointer user_data) {
        gtk_dialog_response(GTK_DIALOG(user_data), GTK_RESPONSE_ACCEPT);
    };
    g_signal_connect(submit_button, "clicked", G_CALLBACK(submit_callback), load_dialog);
    g_signal_connect(entry, "activate", G_CALLBACK(submit_callback), load_dialog);

    gtk_window_set_child(GTK_WINDOW(load_dialog), dialog_body);
    gtk_window_present(GTK_WINDOW(load_dialog));

    auto dialog_response = +[](GtkDialog* dialog, gint response_id, gpointer user_data) {
        if(response_id != GTK_RESPONSE_ACCEPT) {
            g_debug("Ignoring unhandled response type: %s", g_enum_to_string(gtk_response_type_get_type(), response_id));
            return;
        } else {
            const char *text_buf = gtk_entry_buffer_get_text(GTK_ENTRY_BUFFER(user_data));
            auto res = coot::layla::get_drug_via_wikipedia_and_drugbank_curl(std::string(text_buf));
            LaylaState* self = static_cast<LaylaState*>(g_object_get_data(G_OBJECT(dialog),
                                                                                        "ligand_builder_instance"));
            try {
                if(res.empty()) {
                    throw std::runtime_error("Could not fetch MolFile from the internet.");
                }
                RDKit::RWMol* mol = RDKit::MolFileToMol(res,true,false,false);
                if(!mol) {
                    throw std::runtime_error("RDKit::RWMol* is a nullptr. The MolFile could not be loaded.");
                }
                g_info("Molecule Fetch: Molecule constructed.");
                int new_mol_id = self->append_molecule(mol);
                if(new_mol_id >= 0) {
                    self->current_filesave_molecule = static_cast<unsigned int>(new_mol_id);
                    self->current_filesave_filename = res;
                }
                gtk_window_destroy(GTK_WINDOW(dialog));
                // todo: optionally delete the file
            } catch(std::exception& e) {
                g_warning("MolFile Fetch error: %s",e.what());
                auto* message = gtk_message_dialog_new(
                    GTK_WINDOW(dialog), 
                    GTK_DIALOG_DESTROY_WITH_PARENT, 
                    GTK_MESSAGE_ERROR, 
                    GTK_BUTTONS_CLOSE, 
                    "Error: Molecule could not be fetched.\n%s", 
                    e.what()
                );
                g_signal_connect(message,"response",G_CALLBACK(+[](GtkDialog* message_dialog, gint response_id, gpointer user_data){
                    gtk_window_close(GTK_WINDOW(message_dialog));
                }),nullptr);
                gtk_widget_show(message);
            }
        }
   };

   g_signal_connect(load_dialog, "response", G_CALLBACK(dialog_response), entry_buf);
}

void LaylaState::file_new() {
    if(this->has_unsaved_changes()) {
        this->unsaved_changes_dialog_purpose = UnsavedChangesDialogPurpose::FileNew;
        auto* win = gtk_builder_get_object(global_layla_gtk_builder, "layla_unsaved_changes_dialog");
        gtk_window_present(GTK_WINDOW(win));
        return;
    }
    this->reset();
    
}

void LaylaState::file_save() {
    if(this->current_filesave_filename.has_value() && this->current_filesave_molecule.has_value()) {
        save_file(this->current_filesave_molecule.value(), this->current_filesave_filename->c_str());
    } else {
        file_save_as();
    }

}

void LaylaState::save_file(unsigned int idx, const char* filename, GtkWindow* parent) noexcept {
    try {
        const auto* mol = coot_ligand_editor_canvas_get_rdkit_molecule(this->canvas, idx);
        RDKit::MolToMolFile(*mol,std::string(filename));
        g_info("MolFile Save: Molecule file saved.");
        this->update_status("File saved.");
        this->current_filesave_filename = std::string(filename);
        this->current_filesave_molecule = idx;
        this->unsaved_changes = false;
    } catch(std::exception& e) {
        g_warning("MolFile Save error: %s",e.what());
        auto* message = gtk_message_dialog_new(
            parent, 
            GTK_DIALOG_DESTROY_WITH_PARENT, 
            GTK_MESSAGE_ERROR, 
            GTK_BUTTONS_CLOSE, 
            "Error: Molecule could not be saved to file.\n%s", 
            e.what()
        );
        gtk_widget_show(message);
    }
}

void LaylaState::run_file_save_dialog(unsigned int molecule_idx) noexcept {

#if (GTK_MAJOR_VERSION == 4 && GTK_MINOR_VERSION >= 10) || (GTK_MAJOR_VERSION == 5)

    auto* save_dialog = gtk_file_dialog_new();
    // This isn't the best practice but it tremendously simplifies things
    // by saving us from unnecessary boilerplate.
    g_object_set_data(G_OBJECT(save_dialog), "ligand_builder_instance", this);
    gtk_file_dialog_save(save_dialog, this->main_window, NULL, +[](GObject* source_object, GAsyncResult* res, gpointer user_data){
        GError** e = NULL;
        GFile* file = gtk_file_dialog_save_finish(GTK_FILE_DIALOG(source_object), res, e);
        unsigned int molecule_idx = GPOINTER_TO_UINT(user_data);
        LaylaState* self = (LaylaState*) g_object_get_data(G_OBJECT(source_object), "ligand_builder_instance");
        if(file) {
            //g_info("I have a file");
            const char* path = g_file_get_path(file);
            self->save_file(molecule_idx, path,GTK_WINDOW(source_object));
            g_object_unref(file);
        }
        if(e) {
            g_info("Save File: No file was given.");
            g_error_free(*e);
        }
    }, GUINT_TO_POINTER(molecule_idx));
#else
#warning "You're compiling Layla with an unsupported version of GTK. Some functionality will be broken."
    g_warning("Layla has been compiled with an unsupported version of GTK. Some functionality is broken.");
#endif
}

void LaylaState::file_save_as() {
    auto mol_count = coot_ligand_editor_canvas_get_molecule_count(this->canvas);
    if(mol_count == 1) {
        auto m_idx = coot_ligand_editor_canvas_get_idx_of_first_molecule(this->canvas);
        run_file_save_dialog(m_idx);
    } else if(mol_count == 0) {
        update_status("Nothing to be saved!");
    } else {
        auto* mol_chooser_window = gtk_window_new();
        gtk_window_set_title(GTK_WINDOW(mol_chooser_window), "Molecule chooser");
        g_object_set_data(G_OBJECT(mol_chooser_window), "ligand_builder_instance", this);
        g_object_set_data(G_OBJECT(mol_chooser_window),"chosen_molecule",GINT_TO_POINTER(-1));
        gtk_window_set_transient_for(GTK_WINDOW(mol_chooser_window), this->main_window);
        auto* mol_chooser_box = gtk_box_new(GTK_ORIENTATION_VERTICAL,10);
        gtk_window_set_child(GTK_WINDOW(mol_chooser_window), mol_chooser_box);
        auto* mol_chooser_label = gtk_label_new("Choose molecule to be written to a file.");
        gtk_box_append(GTK_BOX(mol_chooser_box), mol_chooser_label);
        auto* mol_chooser_list_box = gtk_list_box_new();
        gtk_box_append(GTK_BOX(mol_chooser_box), mol_chooser_list_box);

        auto max_mol_id = coot_ligand_editor_canvas_get_max_molecule_idx(this->canvas);
        for(unsigned int i = 0; i <= max_mol_id; i++) {
            auto label_str = coot_ligand_editor_canvas_get_smiles_for_molecule(this->canvas,i);
            auto* label = gtk_label_new(label_str.c_str());
            gtk_list_box_append(GTK_LIST_BOX(mol_chooser_list_box),label);
        }

        g_signal_connect(mol_chooser_list_box, "row-activated", G_CALLBACK(+[](GtkListBox* self, GtkListBoxRow* row, gpointer user_data){
            auto idx = gtk_list_box_row_get_index(row);
            GtkWindow* window = GTK_WINDOW(user_data);
            g_object_set_data(G_OBJECT(window),"chosen_molecule",GINT_TO_POINTER(idx));
        }), mol_chooser_window);

        auto* mol_chooser_button_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL,10);
        gtk_box_append(GTK_BOX(mol_chooser_box), mol_chooser_button_box);
        auto* cancel_button =  gtk_button_new_with_label("Cancel");
        auto* ok_button =  gtk_button_new_with_label("Ok");
        gtk_box_append(GTK_BOX(mol_chooser_button_box), cancel_button);
        gtk_box_append(GTK_BOX(mol_chooser_button_box), ok_button);

        g_signal_connect(G_OBJECT(cancel_button), "clicked", G_CALLBACK(+[](GtkButton* button, gpointer userdata){
            GtkWindow* window = GTK_WINDOW(userdata);
            gtk_window_destroy(window);
        }), mol_chooser_window);

        g_signal_connect(G_OBJECT(ok_button), "clicked", G_CALLBACK(+[](GtkButton* button, gpointer userdata){
            GtkWindow* window = GTK_WINDOW(userdata);
            int chosen_molecule = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(window), "chosen_molecule"));
            LaylaState* self = (LaylaState*) g_object_get_data(G_OBJECT(window), "ligand_builder_instance");
            if(chosen_molecule == -1) {
               // gtk_message_dialog_new() is deprecated now
               GtkWidget* message = gtk_message_dialog_new(window,
                                                           GTK_DIALOG_DESTROY_WITH_PARENT, // flags
                                                           GTK_MESSAGE_ERROR, // type
                                                           GTK_BUTTONS_CLOSE, // buttons
                                                           "Nothing was chosen!" // message_format
                                                           ); // G_GNUC_PRINTF (5, 6);

                gtk_widget_set_visible(message, TRUE);
                g_info("Nothing was chosen.");
            } else {
                self->run_file_save_dialog(chosen_molecule);
                gtk_window_destroy(window);
            }
        }), mol_chooser_window);

        gtk_window_present(GTK_WINDOW(mol_chooser_window));
    }
}

void LaylaState::file_open() {

#if (GTK_MAJOR_VERSION == 4 && GTK_MINOR_VERSION >= 10) || (GTK_MAJOR_VERSION == 5)

    auto* open_dialog = gtk_file_dialog_new();
    gtk_file_dialog_open(open_dialog, this->main_window, NULL, +[](GObject* source_object, GAsyncResult* res, gpointer user_data){
        GError** e = NULL;
        GFile* file = gtk_file_dialog_open_finish(GTK_FILE_DIALOG(source_object), res, e);
        LaylaState* self = (LaylaState*) user_data;
        if(file) {
            //g_info("I have a file");
            const char* path = g_file_get_path(file);
            try {
                RDKit::RWMol* mol = RDKit::MolFileToMol(std::string(path),true,false,false);
                if(!mol) {
                    throw std::runtime_error("RDKit::RWMol* is a nullptr. The MolFile could not be loaded.");
                }
                g_info("MolFile Import: Molecule constructed.");
                int new_mol_id = self->append_molecule(mol);
                if(new_mol_id >= 0) {
                    self->current_filesave_molecule = static_cast<unsigned int>(new_mol_id);
                    self->current_filesave_filename = std::string(path);
                }
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
            g_error_free(*e);
        }
    }, this);
#else
#warning "You're compiling Layla with an unsupported version of GTK. Some functionality will be broken."
    g_warning("Layla has been compiled with an unsupported version of GTK. Some functionality is broken.");
#endif
}

void LaylaState::file_export(ExportMode mode) {

#if (GTK_MAJOR_VERSION == 4 && GTK_MINOR_VERSION >= 10) || (GTK_MAJOR_VERSION == 5)

    auto* export_dialog = gtk_file_dialog_new();
    auto* mode_ptr = new ExportMode(mode);
    // This isn't the best practice but it tremendously simplifies things
    // by saving us from unnecessary boilerplate.
    g_object_set_data(G_OBJECT(export_dialog), "ligand_builder_instance", this);
    gtk_file_dialog_save(export_dialog, this->main_window, NULL, +[](GObject* source_object, GAsyncResult* res, gpointer user_data){
        GError** e = NULL;
        GFile* file = gtk_file_dialog_save_finish(GTK_FILE_DIALOG(source_object), res, e);
        ExportMode* mode_ptr = (ExportMode*) user_data;
        LaylaState* self = (LaylaState*) g_object_get_data(G_OBJECT(source_object), "ligand_builder_instance");
        if(file) {
            //g_info("I have a file");
            auto path = std::string(g_file_get_path(file));
            cairo_surface_t* target = nullptr;
            auto draw = [&](){
                if(target) {
                    cairo_t* cr = cairo_create(target);
                    coot_ligand_editor_canvas_draw_on_cairo_surface(self->canvas, cr);
                }
            };
            auto ends_with = [](std::string const & value, std::string const & ending){
                if (ending.size() > value.size()) 
                    return false;
                return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
            };
            int width = gtk_widget_get_size(GTK_WIDGET(self->canvas),GTK_ORIENTATION_HORIZONTAL);
            int height = gtk_widget_get_size(GTK_WIDGET(self->canvas),GTK_ORIENTATION_VERTICAL);
            switch (*mode_ptr) {
                case ExportMode::PDF: {
                    if(!ends_with(path, ".pdf")) {
                        path += ".pdf";
                    }
                    target = cairo_pdf_surface_create(path.c_str(), width, height);
                    draw();
                    break;
                }
                case ExportMode::PNG: {
                    // target = cairo_image_surface_create(CAIRO_FORMAT_RGBA128F, width, height);
                    target = cairo_image_surface_create(CAIRO_FORMAT_RGB24, width, height);
                    draw();
                    if(!ends_with(path, ".png")) {
                        path += ".png";
                    }
                    cairo_surface_write_to_png(target, path.c_str());
                    break;
                }
                case ExportMode::SVG: {
                    if(!ends_with(path, ".svg")) {
                        path += ".svg";
                    }
                    target = cairo_svg_surface_create(path.c_str(), width, height);
                    draw();
                    break;
                }
                default: {
                    break;
                }
            }
            if(target) {
                cairo_surface_destroy(target);
            }
            g_object_unref(file);
        }
        if(e) {
            g_info("Export: No file was given.");
            g_error_free(*e);
        }
        delete mode_ptr;
    }, mode_ptr);
#else
#warning "You're compiling Layla with an unsupported version of GTK. Some functionality will be broken."
    g_warning("Layla has been compiled with an unsupported version of GTK. Some functionality is broken.");
#endif
}

void LaylaState::file_exit() {
    gtk_window_close(GTK_WINDOW(this->main_window));
}

void LaylaState::edit_undo() {
    coot_ligand_editor_canvas_undo_edition(this->canvas);
}

void LaylaState::edit_redo() {
    coot_ligand_editor_canvas_redo_edition(this->canvas);
}

void LaylaState::switch_display_mode(ligand_editor_canvas::DisplayMode mode) {
    coot_ligand_editor_canvas_set_display_mode(this->canvas, mode);
}

void coot::layla::initialize_global_instance(CootLigandEditorCanvas* canvas, GtkWindow* win, GtkLabel* status_label) {
    global_instance = new LaylaState(canvas,win,status_label);
    g_info("Global instance of LaylaState has been initialized at: %p",global_instance);
}

std::optional<LaylaState::InchiKeyLookupResult> LaylaState::lookup_inchi_key(const std::string& inchi_key) const noexcept {
    if(!this->inchi_key_database) {
        return std::nullopt;
    }
    const auto& db = *this->inchi_key_database;
    auto it = db.find(inchi_key);
    if (it == db.end()) {
        return std::nullopt;
    }
    return InchiKeyLookupResult{it->second.first, it->second.second};
}
