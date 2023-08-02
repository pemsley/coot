#include "ligand_builder_state.hpp"
#include "ligand_editor_canvas.hpp"
#include <exception>
#include <memory>
#include <stdexcept>
#include <gtk/gtk.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/FileParsers/FileParsers.h>
#include "lidia-core/rdkit-interface.hh"
#include <string>

using namespace coot::ligand_editor;

LigandBuilderState::LigandBuilderState(CootLigandEditorCanvas* canvas_widget, GtkWindow* win, GtkLabel* status_label) noexcept {
    this->canvas = canvas_widget;
    this->main_window = win;
    this->status_label = status_label;
    //g_object_set_data(G_OBJECT(win), "ligand_builder_instance", this);
}

CootLigandEditorCanvas* LigandBuilderState::get_canvas() const noexcept {
    return this->canvas;
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
    // g_warning("TODO: Implement void LigandBuilderState::file_import_molecule()");

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

    GtkWidget *submit_button = gtk_button_new_with_label("Submit");
    gtk_box_append(GTK_BOX(dialog_body), submit_button);

    auto submit_callback = +[] (GtkButton *button, gpointer user_data) {
        gtk_dialog_response(GTK_DIALOG(user_data), GTK_RESPONSE_ACCEPT);
    };
    g_signal_connect(submit_button, "clicked", G_CALLBACK(submit_callback), load_dialog);

    gtk_window_set_child(GTK_WINDOW(load_dialog), dialog_body);
    gtk_window_present(GTK_WINDOW(load_dialog));

    auto dialog_response = +[](GtkDialog* dialog, gint response_id, gpointer user_data) {
        if(response_id != GTK_RESPONSE_ACCEPT) {
            g_debug("Ignoring unhandled response type: %s", g_enum_to_string(gtk_response_type_get_type(), response_id));
            return;
        } else {
            const char *text_buf = gtk_entry_buffer_get_text(GTK_ENTRY_BUFFER(user_data));
            std::string monomer_type(text_buf);
            int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
            LigandBuilderState* self = static_cast<LigandBuilderState*>(g_object_get_data(G_OBJECT(dialog),
                                                                                        "ligand_builder_instance"));
            // what is 42???
            self->monomer_library_info_store.try_dynamic_add(monomer_type, 42);
            std::pair<bool, dictionary_residue_restraints_t> p =
                self->monomer_library_info_store.get_monomer_restraints(monomer_type, imol_enc);
            if (p.first) {
                bool show_hydrogens_status = false;
                auto mol = std::make_unique<RDKit::RWMol>(coot::rdkit_mol(p.second));
                if (! show_hydrogens_status) {
                    coot::remove_non_polar_Hs(mol.get());
                }
                self->append_molecule(mol.release());
                self->current_filesave_molecule = coot_ligand_editor_get_molecule_count(self->canvas) - 1;
                gtk_window_destroy(GTK_WINDOW(dialog));
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

   g_signal_connect(load_dialog, "response", G_CALLBACK(dialog_response), entry_buf);
}

void LigandBuilderState::run_choose_element_dialog() {
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
            RDKit::RWMol* molecule = RDKit::SmilesToMol(gtk_entry_buffer_get_text(text_buf));
            LigandBuilderState* state = (LigandBuilderState*) g_object_get_data(G_OBJECT(dialog), "ligand_builder_instance");
            coot_ligand_editor_set_active_tool(state->canvas, std::move(insertion_tool));
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

void LigandBuilderState::file_fetch_molecule() {
    g_warning("TODO: Implement void LigandBuilderState::file_fetch_molecule()");
}

void LigandBuilderState::file_new() {
    g_warning("TODO: Implement void LigandBuilderState::file_new()");
}

void LigandBuilderState::file_save() {
    if(this->current_filesave_filename.has_value() && this->current_filesave_molecule.has_value()) {
        save_file(this->current_filesave_molecule.value(), this->current_filesave_filename->c_str());
    } else {
        file_save_as();
    }

}

void LigandBuilderState::save_file(unsigned int idx, const char* filename, GtkWindow* parent) noexcept {
    try {
        const auto* mol = coot_ligand_editor_get_rdkit_molecule(this->canvas, idx);
        RDKit::MolToMolFile(*mol,std::string(filename));
        g_info("MolFile Save: Molecule file saved.");
        this->update_status("File saved.");
        this->current_filesave_filename = std::string(filename);
        this->current_filesave_molecule = idx;
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

void LigandBuilderState::run_file_save_dialog(unsigned int molecule_idx) noexcept {
    auto* save_dialog = gtk_file_dialog_new();
    // This isn't the best practice but it tremendously simplifies things
    // by saving us from unnecessary boilerplate.
    g_object_set_data(G_OBJECT(save_dialog), "ligand_builder_instance", this);
    gtk_file_dialog_save(save_dialog, this->main_window, NULL, +[](GObject* source_object, GAsyncResult* res, gpointer user_data){
        GError** e = NULL;
        GFile* file = gtk_file_dialog_save_finish(GTK_FILE_DIALOG(source_object), res, e);
        unsigned int molecule_idx = GPOINTER_TO_UINT(user_data);
        LigandBuilderState* self = (LigandBuilderState*) g_object_get_data(G_OBJECT(source_object), "ligand_builder_instance");
        if(file) {
            //g_info("I have a file");
            const char* path = g_file_get_path(file);
            self->save_file(molecule_idx, path,GTK_WINDOW(source_object));
            g_object_unref(file);
        }
        if(e) {
            g_info("Save File: No file was given.");
            g_object_unref(*e);
        }
    }, GUINT_TO_POINTER(molecule_idx));
}

void LigandBuilderState::file_save_as() {
    auto mol_count = coot_ligand_editor_get_molecule_count(this->canvas);
    if(mol_count == 1) {
        run_file_save_dialog(0);
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

        for(unsigned int i = 0; i < mol_count; i++) {
            auto label_str = coot_ligand_editor_get_smiles_for_molecule(this->canvas,i);
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
            LigandBuilderState* self = (LigandBuilderState*) g_object_get_data(G_OBJECT(window), "ligand_builder_instance");
            if(chosen_molecule == -1) {
                auto* message = gtk_message_dialog_new(
                    window, 
                    GTK_DIALOG_DESTROY_WITH_PARENT, 
                    GTK_MESSAGE_ERROR, 
                    GTK_BUTTONS_CLOSE,
                    "Nothing was chosen!",
                    NULL
                );
                gtk_widget_show(message);
                g_info("Nothing was chosen.");
            } else {
                self->run_file_save_dialog(chosen_molecule);
                gtk_window_destroy(window);
            }
        }), mol_chooser_window);

        gtk_window_present(GTK_WINDOW(mol_chooser_window));
    }
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
                self->current_filesave_molecule = coot_ligand_editor_get_molecule_count(self->canvas) - 1;
                self->current_filesave_filename = std::string(path);
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

void LigandBuilderState::file_exit() {
    // todo: this should probably do some checks before just closing
    gtk_window_close(GTK_WINDOW(this->main_window));
}

void LigandBuilderState::edit_undo() {
    coot_ligand_editor_undo_edition(this->canvas);
}

void LigandBuilderState::edit_redo() {
    coot_ligand_editor_redo_edition(this->canvas);
}

void LigandBuilderState::switch_display_mode(ligand_editor_canvas::DisplayMode mode) {
    coot_ligand_editor_set_display_mode(this->canvas, mode);
}

void LigandBuilderState::run_apply() {
    g_warning("TODO: Implement 'Apply'");
    auto* apply_dialog = gtk_dialog_new();
    gtk_window_set_transient_for(GTK_WINDOW(apply_dialog), this->main_window);
    // This isn't the best practice but it tremendously simplifies things
    // by saving us from unnecessary boilerplate.
    g_object_set_data(G_OBJECT(apply_dialog), "ligand_builder_instance", this);
    gtk_window_set_title(GTK_WINDOW(apply_dialog), "Apply");
    
    auto* dialog_body = gtk_box_new(GTK_ORIENTATION_VERTICAL,10);
    gtk_widget_set_margin_bottom(dialog_body, 10);
    gtk_widget_set_margin_end(dialog_body, 10);
    gtk_widget_set_margin_start(dialog_body, 10);
    gtk_widget_set_margin_top(dialog_body, 10);

    // On pressing the "Apply" button the user needs to be given 3 choices/questions

    // 1: What is the monomer id? (default "LIG", "INH" or allow the user to specify)
    // 2: What input format (SMILES or MolFile)
    // 3: Which generator to use ("Acedrg", "Grade2", "Pyrogen")

    auto* submit_button = gtk_button_new_with_label("Accept");
    gtk_box_append(GTK_BOX(dialog_body),submit_button);
    g_signal_connect(submit_button, "clicked", G_CALLBACK(+[](GtkButton* btn, gpointer user_data){
        gtk_dialog_response(GTK_DIALOG(user_data), GTK_RESPONSE_ACCEPT);
    }), apply_dialog);

    g_signal_connect(apply_dialog, "response", G_CALLBACK(+[](GtkDialog* dialog, gint response_id, gpointer user_data){
        if(response_id != GTK_RESPONSE_ACCEPT) {
            g_debug("Ignoring unhandled response type: %s",g_enum_to_string(gtk_response_type_get_type(), response_id));
            return;
        }
        
        gtk_window_destroy(GTK_WINDOW(dialog));

    }),nullptr);

    gtk_window_set_child(GTK_WINDOW(apply_dialog),dialog_body);
    gtk_window_present(GTK_WINDOW(apply_dialog));
}

void coot::ligand_editor::initialize_global_instance(CootLigandEditorCanvas* canvas, GtkWindow* win, GtkLabel* status_label) {
    global_instance = new LigandBuilderState(canvas,win,status_label);
    g_info("Global instance of LigandBuilderState has been initialized at: %p",global_instance);
}
