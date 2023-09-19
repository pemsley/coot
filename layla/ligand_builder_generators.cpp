/* layla/ligand_builder_generators.cpp
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

#include "ligand_builder_generators.hpp"
#include "ligand_builder_state.hpp"
#include <cstddef>
#include <cstring>
#include <glib.h>
#include <memory>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/FileParsers/FileParsers.h>

struct GeneratorTaskData {

    std::unique_ptr<coot::ligand_editor::GeneratorRequest> request;
    std::unique_ptr<std::string> file_contents;
    GtkProgressBar* progress_bar;
    GtkWindow* progress_dialog;
    GtkButton* progress_dialog_close_button;
    GtkTextBuffer* stdout_ui_textbuffer;
    GtkLabel* dialog_status_label;
    GtkSpinner* spinner;

    GSubprocess* subprocess;
    bool subprocess_running;

    GInputStream* input_stream;
    std::unique_ptr<std::string> stdout_read;

    void initialize(coot::ligand_editor::GeneratorRequest&& request) {
        g_warning("void GeneratorTaskData::initialize() called.");
        
        using namespace coot::ligand_editor;

        this->progress_bar = (GtkProgressBar*) gtk_builder_get_object(global_layla_gtk_builder, "layla_generator_progress_dialog_progress_bar");
        this->progress_dialog = (GtkWindow*) gtk_builder_get_object(global_layla_gtk_builder, "layla_generator_progress_dialog");
        this->progress_dialog_close_button = (GtkButton*) gtk_builder_get_object(global_layla_gtk_builder, "layla_generator_progress_dialog_close_button");
        this->stdout_ui_textbuffer = gtk_text_view_get_buffer(
            (GtkTextView*) gtk_builder_get_object(global_layla_gtk_builder, "layla_generator_progress_dialog_stdout_textview")
        );
        this->dialog_status_label = (GtkLabel*) gtk_builder_get_object(global_layla_gtk_builder, "layla_generator_progress_dialog_status_label");
        this->spinner = (GtkSpinner*) gtk_builder_get_object(global_layla_gtk_builder, "layla_generator_progress_dialog_spinner");
        this->request = std::make_unique<GeneratorRequest>(request);
        this->file_contents = nullptr;
        this->subprocess = nullptr;
        this->stdout_read = std::make_unique<std::string>();
        this->input_stream = nullptr;
        this->subprocess_running = false;
    }

    void cleanup() {
        g_warning("void GeneratorTaskData::cleanup() called.");
        if (this->subprocess) {
            g_object_unref(subprocess);
        }
        if (this->input_stream) {
            g_warning("todo: Make sure that `input_stream` does not leak and there's no crash.");
           // g_object_unref(input_stream);
        }
        this->request.reset();
        this->file_contents.reset();
        this->stdout_read.reset();
    }
};

std::string coot::ligand_editor::GeneratorRequest::get_filename() const {

    std::string file_name;
    switch (generator) {
        case Generator::Grade2: {
            file_name = "grade2-";
            break;
        }
        default:
        case Generator::Acedrg: {
            file_name = "acedrg-";
            break;
        }
    }
    file_name += monomer_id;
    switch(input_format) {
        case InputFormat::MolFile: {
            file_name += ".mol";
            break;
        }
        default:
        case InputFormat::SMILES: {
            file_name += ".smi";
            break;
        }
    }
    return file_name;
}

std::vector<std::string> coot::ligand_editor::GeneratorRequest::build_commandline() const {

    std::vector<std::string> ret;
    ret.push_back(this->executable_path.value());
    auto input_filename = this->get_filename();
    switch(generator) {
        case Generator::Grade2: {
            switch(input_format) {
                case InputFormat::MolFile: {
                    g_error("Todo: implement molfile for grade2");
                    break;
                }
                default:
                case InputFormat::SMILES: {
                    g_error("Todo: implement smiles for grade2");
                    break;
                }
            }
            break;
        }
        default:
        case Generator::Acedrg: {
            switch(input_format) {
                case InputFormat::MolFile: {
                    ret.push_back("-m");
                    ret.push_back(input_filename);
                    break;
                }
                default:
                case InputFormat::SMILES: {
                    // 20230916-PE "-p" and "-z" should be conditionally added, but
                    // let's just shove them in for now, while we are testing
                    ret.push_back("-p");
                    ret.push_back("-z");
                    ret.push_back("-i");
                    ret.push_back(input_filename);
                    break;
                }
            }
            ret.push_back("-r");
            ret.push_back(this->monomer_id);
            ret.push_back("-o");
            ret.push_back(std::string("acedrg-") + this->monomer_id);
            break;
        }
    }
    return ret;
}


// Forward declaration
void write_input_file_finish(GObject* file_object, GAsyncResult* res, gpointer user_data);

void write_input_file_async(GTask* task) {
    GCancellable* cancellable = g_task_get_cancellable(task);
    GeneratorTaskData* task_data = (GeneratorTaskData*) g_task_get_task_data(G_TASK(task));
    std::string file_contents;
    std::string file_name = task_data->request->get_filename();

    using InputFormat = coot::ligand_editor::GeneratorRequest::InputFormat;
    switch(task_data->request->input_format) {
        case InputFormat::MolFile: {
            RDKit::RWMol* mol = RDKit::SmilesToMol(task_data->request->molecule_smiles);
            file_contents = RDKit::MolToMolBlock(*mol);
            break;
        }
        default:
        case InputFormat::SMILES: {
            file_contents = task_data->request->molecule_smiles;
            break;
        }
    }
    GFile* file = g_file_new_for_path(file_name.c_str());
    task_data->file_contents = std::make_unique<std::string>(std::move(file_contents));
    gtk_label_set_text(task_data->dialog_status_label, "Writing input file...");

    g_file_replace_contents_async(
        file, 
        task_data->file_contents->c_str(), 
        task_data->file_contents->size(), 
        NULL, 
        FALSE, 
        G_FILE_CREATE_REPLACE_DESTINATION, 
        cancellable, 
        write_input_file_finish, 
        task
    );
}

// Forward declaration
void resolve_target_generator_executable(GTask* task);

void write_input_file_finish(GObject* file_object, GAsyncResult* res, gpointer user_data) {
    GTask* task = G_TASK(user_data);
    GeneratorTaskData* task_data = (GeneratorTaskData*) g_task_get_task_data(G_TASK(task));
    GFile* file = G_FILE(file_object);
    GError* err = NULL;
    bool file_io_res = g_file_replace_contents_finish(file, res, NULL, &err);
    g_object_unref(file);
    if(!file_io_res) {
        g_warning("Write failed");
        g_task_return_error(task, err);
        return;
    }
    g_warning("Write ok.");
    gtk_label_set_text(task_data->dialog_status_label, "Input file has been written.");
    resolve_target_generator_executable(task);
}

// Forward declaration
void launch_generator_async(GTask* task);

void resolve_target_generator_executable(GTask* task) {
    GCancellable* cancellable = g_task_get_cancellable(task);
    GeneratorTaskData* task_data = (GeneratorTaskData*) g_task_get_task_data(G_TASK(task));
    using Generator = coot::ligand_editor::GeneratorRequest::Generator;

    switch(task_data->request->generator) {
        default:
        case Generator::Acedrg: {
            g_warning("todo: Implement resolving acedrg executable");
            task_data->request->executable_path = "acedrg";
            break;
        }
        case Generator::Grade2: {
            g_error("todo: Implement resolving Grade2 executable");
            break;
        }
    }
    launch_generator_async(task);
}

// Forward declaration
void launch_generator_finish(GObject* subprocess_object, GAsyncResult* res, gpointer user_data);
void pipe_reader(gpointer user_data);

void launch_generator_async(GTask* task) {

    GCancellable* cancellable = g_task_get_cancellable(task);
    GeneratorTaskData* task_data = (GeneratorTaskData*) g_task_get_task_data(G_TASK(task));
    GSubprocessLauncher* launcher = g_subprocess_launcher_new(G_SUBPROCESS_FLAGS_STDOUT_PIPE);
    std::vector<std::string> argv = task_data->request->build_commandline();
    gsize slice_size = sizeof(gchar*) * (argv.size() + 1);
    const gchar **argv_raw = (const gchar **)g_slice_alloc0(slice_size);
    for(unsigned int i = 0; i < argv.size(); i++) {
        argv_raw[i] = argv[i].data();
    }
    GError* err = NULL;
    GSubprocess* subprocess = g_subprocess_launcher_spawnv(launcher, argv_raw, &err);
    g_object_unref(launcher);
    g_slice_free1(slice_size, argv_raw);
    if(!subprocess) {
        g_warning("The subprocess could not be spawned.");
        if(err) {
            g_task_return_error(task, err);
        } else {
            g_task_return_boolean(task, false);
        }
        return;
    }
    g_warning("Subprocess spawned!");
    task_data->subprocess = g_object_ref(subprocess);
    GInputStream* input_stream = g_subprocess_get_stdout_pipe(subprocess);
    task_data->input_stream = input_stream;
    task_data->subprocess_running = true;
    g_subprocess_wait_check_async(subprocess, cancellable, launch_generator_finish, task);
    gtk_label_set_text(task_data->dialog_status_label, "Child process has been launched.");
    g_timeout_add(50, [](gpointer user_data){
        GTask* task = G_TASK(user_data);
        GeneratorTaskData* task_data = (GeneratorTaskData*) g_task_get_task_data(G_TASK(task));
        int should_run = task_data->subprocess_running;
        if(!should_run) {
            g_object_unref(task);
            g_warning("ProgressBar loop exits.");
        } else {
            gtk_progress_bar_pulse(task_data->progress_bar);
        }
        return should_run;
    }, g_object_ref(task));

    g_idle_add_once(pipe_reader, g_object_ref(task));
}

void pipe_reader(gpointer user_data) {
    GTask* task = G_TASK(user_data);
    GCancellable* cancellable = g_task_get_cancellable(task);
    GeneratorTaskData* task_data = (GeneratorTaskData*) g_task_get_task_data(G_TASK(task));
    auto callback = [](GObject* obj, GAsyncResult* res, gpointer user_data){
        GTask* task = G_TASK(user_data);
        GeneratorTaskData* task_data = (GeneratorTaskData*) g_task_get_task_data(G_TASK(task));
        GInputStream* stream = G_INPUT_STREAM(obj);
        GError* err = NULL;
        GBytes* bytes = g_input_stream_read_bytes_finish(stream, res, &err);
        bool should_go_on = true;
        if(err) {
            g_warning("Stream reading operation ended due to error: %s", err->message);
            g_error_free(err);
            should_go_on = false;
        }
        if(bytes) {
            auto size = g_bytes_get_size(bytes);
            if(size == 0) {
                g_warning("Stream reading operation ended due to EOF");
                should_go_on = false;
            } else {
                std::string buf;
                buf.reserve(size + 1);
                memcpy(buf.data(), g_bytes_get_data(bytes, NULL), size);
                *(buf.data() + size) = (char) 0;
                g_debug("Read this: %s", buf.c_str());
                (*task_data->stdout_read) += buf;
                GtkTextIter iter;
                gtk_text_buffer_get_end_iter(task_data->stdout_ui_textbuffer, &iter);
                gtk_text_buffer_insert(task_data->stdout_ui_textbuffer, &iter, buf.c_str(), -1);
            }
            g_bytes_unref(bytes);
        }
        if(should_go_on) {
            g_idle_add_once(pipe_reader, g_object_ref(task));
        }
        g_object_unref(task);
    };
    g_input_stream_read_bytes_async(task_data->input_stream, 4, G_PRIORITY_DEFAULT, cancellable, callback, g_object_ref(task));
    g_object_unref(task);
};

void launch_generator_finish(GObject* subprocess_object, GAsyncResult* res, gpointer user_data) {
    GTask* task = G_TASK(user_data);
    GeneratorTaskData* task_data = (GeneratorTaskData*) g_task_get_task_data(G_TASK(task));
    GSubprocess* subprocess = G_SUBPROCESS(subprocess_object);
    GError* err = NULL;
    bool io_res = g_subprocess_wait_check_finish(subprocess, res, &err);
    task_data->subprocess_running = false;
    gint exit_status = g_subprocess_get_exit_status(subprocess);
    g_warning("Exit status: %i", exit_status);
    g_object_unref(subprocess);
    if(!io_res || exit_status != 0) {
        g_warning("Subprocess failed.");
        if(err) {
            g_task_return_error(task, err);
        } else {
            g_task_return_boolean(task, false);
        }
        return;
    } 
    g_warning("Subprocess exited ok.");
    g_task_return_boolean(task, true);
}

GCancellable* coot::ligand_editor::run_generator_request(GeneratorRequest request, CootLaylaNotifier* notifier) {
    // Useless dummy
    GObject* dummy = (GObject*)g_object_new(G_TYPE_OBJECT, NULL);
    GCancellable* cancellable = g_cancellable_new();
    GeneratorTaskData* task_data = g_slice_new0(GeneratorTaskData);
    task_data->initialize(std::move(request));

    GtkTextIter start, end;
    gtk_text_buffer_get_end_iter(task_data->stdout_ui_textbuffer, &end);
    gtk_text_buffer_get_start_iter(task_data->stdout_ui_textbuffer, &start);
    gtk_text_buffer_delete(task_data->stdout_ui_textbuffer, &start, &end); 

    gtk_label_set_text(task_data->dialog_status_label, "");
    gtk_spinner_set_spinning(task_data->spinner, true);

    auto task_completed_callback = [](GObject* obj, GAsyncResult* res, gpointer user_data){
        g_warning("Task completed callback!");
        GTask* task = G_TASK(res);
        GeneratorTaskData* task_data = (GeneratorTaskData*) g_task_get_task_data(task);
        CootLaylaNotifier* notifier = COOT_COOT_LAYLA_NOTIFIER(user_data);
        // todo: cleanup after child process (if any) and report results
        GError* err = NULL;
        if(!g_task_propagate_boolean(task, &err)) {
            if(err) {
                std::string label_text = "Operation failed: ";
                label_text += err->message;
                gtk_label_set_text(task_data->dialog_status_label, label_text.c_str());
                g_warning("Task failed. Error: %s", err->message);
                g_error_free(err);
            }
        } else {
            gtk_label_set_text(task_data->dialog_status_label, "Operation completed successfully!");
            g_warning("Task finished successfully!");
        }

        // Delete the useless object
        g_object_unref(obj);

        g_object_unref(task);
        g_object_unref(notifier);

        // We need to manually remove our own cancellable from here.
        g_object_unref(global_generator_request_task_cancellable);
        global_generator_request_task_cancellable = nullptr;
        // Apply necessary UI changes after the task is done
        auto* cancel_button = gtk_builder_get_object(global_layla_gtk_builder, "layla_generator_progress_dialog_cancel_button");
        gtk_widget_set_sensitive(GTK_WIDGET(cancel_button), FALSE);
        auto* accept_button = gtk_builder_get_object(global_layla_gtk_builder, "layla_apply_dialog_accept_button");
        gtk_widget_set_sensitive(GTK_WIDGET(accept_button), TRUE);
        gtk_widget_set_sensitive(GTK_WIDGET(task_data->progress_dialog_close_button), TRUE);

        gtk_spinner_set_spinning(task_data->spinner, false);
    };

    GTask* task = g_task_new(dummy,cancellable,task_completed_callback, notifier);
    g_task_set_task_data(task, task_data, [](gpointer task_data_ptr){
        GeneratorTaskData* task_data = (GeneratorTaskData*) task_data_ptr;
        task_data->cleanup();
        g_slice_free(GeneratorTaskData, task_data_ptr);
    });

    g_warning("Implement 'Apply'");

    g_idle_add_once([](gpointer user_data){
        write_input_file_async(G_TASK(user_data));
    }, task);
    
    // this segfaults:
    // g_object_unref(dummy);

    return cancellable;
    
}
