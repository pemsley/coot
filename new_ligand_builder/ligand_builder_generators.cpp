#include "ligand_builder_generators.hpp"
#include "ligand_builder_state.hpp"
#include <glib.h>
#include <memory>

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

GCancellable* coot::ligand_editor::run_generator_request(GeneratorRequest request) {
    GObject* dummy = (GObject*)g_object_new(G_TYPE_OBJECT, NULL);
    GCancellable* cancellable = g_cancellable_new();
    GeneratorTaskData* task_data = g_slice_new0(GeneratorTaskData);
    task_data->initialize(std::move(request));

    auto task_completed_callback = [](GObject* obj, GAsyncResult* res, gpointer user_data){
        g_warning("Task completed callback!");
        GeneratorTaskData* task_data = (GeneratorTaskData*) g_task_get_task_data(G_TASK(res));
        gtk_window_close(task_data->progress_dialog);
        // todo: cleanup after child process (if any) and report results

        g_object_unref(obj);
        // does this deallocate the task?        
        g_object_unref(res);
    };

    GTask* task = g_task_new(dummy,cancellable,task_completed_callback, nullptr);
    g_task_set_task_data(task, task_data, [](gpointer task_data_ptr){
        GeneratorTaskData* task_data = (GeneratorTaskData*) task_data_ptr;
        task_data->cleanup();
        g_slice_free(GeneratorTaskData, task_data_ptr);
    });

    g_timeout_add_once(5000, [](gpointer user_data){
        g_warning("5 seconds elapsed.");
        g_task_return_boolean(G_TASK(user_data), true);
    }, task);
    
    // this segfaults:
    // g_object_unref(dummy);
    
    // Cannnot deallocate task here because it's still runnning.
    // Where should I do this?

    return cancellable;
    
}