#include "ligand_builder_generators.hpp"
#include "ligand_builder_state.hpp"
#include <glib.h>
#include <memory>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/FileParsers/FileParsers.h>

struct GeneratorTaskData {
    std::unique_ptr<coot::ligand_editor::GeneratorRequest> request;
    GtkProgressBar* progress_bar;
    GtkWindow* progress_dialog;

    void initialize(coot::ligand_editor::GeneratorRequest&& request) {
        g_warning("void GeneratorTaskData::initialize() called.");
        
        using namespace coot::ligand_editor;

        this->progress_bar = (GtkProgressBar*) gtk_builder_get_object(global_layla_gtk_builder, "layla_generator_progress_dialog_progress_bar");
        this->progress_dialog = (GtkWindow*) gtk_builder_get_object(global_layla_gtk_builder, "layla_generator_progress_dialog");
        this->request = std::make_unique<GeneratorRequest>(request);

    }

    void cleanup() {
        g_warning("void GeneratorTaskData::cleanup() called.");

        this->request.reset();
    }
};

void launch_generator_async(GTask* task) {
    GCancellable* cancellable = g_task_get_cancellable(task);
    GeneratorTaskData* task_data = (GeneratorTaskData*) g_task_get_task_data(G_TASK(task));
    //todo: implement me
    g_task_return_boolean(task, true);
}

void launch_generator_finish(GObject* something, GAsyncResult* res, gpointer user_data) {

}

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

void write_input_file_finish(GObject* file_object, GAsyncResult* res, gpointer user_data) {
    GTask* task = G_TASK(user_data);
    GFile* file = G_FILE(file_object);
    GError* err;
    bool file_io_res = g_file_replace_contents_finish(file, res, NULL, &err);
    g_object_unref(file);
    if(!file_io_res) {
        g_warning("Write failed");
        g_task_return_boolean(task, false);
        if(err) {
            g_object_unref(err);
        }
        return;
    } 
    g_warning("Write ok.");
    resolve_target_generator_executable(task);
}

std::string coot::ligand_editor::GeneratorRequest::get_filename() const {
    std::string file_name;
    switch(generator) {
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
    g_file_replace_contents_async(
        file, 
        file_contents.c_str(), 
        file_contents.size(), 
        NULL, 
        FALSE, 
        G_FILE_CREATE_REPLACE_DESTINATION, 
        cancellable, 
        write_input_file_finish, 
        task
    );
}

GCancellable* coot::ligand_editor::run_generator_request(GeneratorRequest request) {
    // Useless dummy
    GObject* dummy = (GObject*)g_object_new(G_TYPE_OBJECT, NULL);
    GCancellable* cancellable = g_cancellable_new();
    GeneratorTaskData* task_data = g_slice_new0(GeneratorTaskData);
    task_data->initialize(std::move(request));

    auto task_completed_callback = [](GObject* obj, GAsyncResult* res, gpointer user_data){
        g_warning("Task completed callback!");
        GeneratorTaskData* task_data = (GeneratorTaskData*) g_task_get_task_data(G_TASK(res));
        gtk_window_close(task_data->progress_dialog);
        // todo: cleanup after child process (if any) and report results

        // Delete the useless object
        g_object_unref(obj);
        // does this deallocate the task?        
        g_object_unref(res);

        // We need to manually remove our own cancellable from here.
        g_object_unref(global_generator_request_task_cancellable);
        global_generator_request_task_cancellable = nullptr;
        // Apply necessary UI changes after the task is done
        auto* accept_button = gtk_builder_get_object(global_layla_gtk_builder, "layla_apply_dialog_accept_button");
        gtk_widget_set_sensitive(GTK_WIDGET(accept_button), TRUE);
    };

    GTask* task = g_task_new(dummy,cancellable,task_completed_callback, nullptr);
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
    
    // Cannnot deallocate task here because it's still runnning.
    // Where should I do this?

    return cancellable;
    
}