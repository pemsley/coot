#include "ligand_builder_generators.hpp"
#include "ligand_builder_state.hpp"
#include <glib.h>

struct GeneratorTaskData {

    void initialize() {
        g_warning("void initialize() called.");
    }

    void cleanup() {
        g_warning("void cleanup() called.");
    }
};

GCancellable* coot::ligand_editor::run_generator_request(GeneratorRequest request) {
    GObject* dummy = (GObject*)g_object_new(G_TYPE_OBJECT, NULL);
    GCancellable* cancellable = g_cancellable_new();
    GeneratorTaskData* task_data = g_slice_new0(GeneratorTaskData);
    task_data->initialize();

    auto task_completed_callback = +[](GObject* obj, GAsyncResult* res, gpointer user_data){
        g_warning("Task callback!");
        auto* progress_dialog = gtk_builder_get_object(global_layla_gtk_builder, "layla_generator_progress_dialog");
        gtk_window_close(GTK_WINDOW(progress_dialog));

        g_object_unref(obj);
        // does this deallocate the task?        
        g_object_unref(res);
    };

    GTask* task = g_task_new(dummy,cancellable,task_completed_callback, nullptr);
    g_task_set_task_data(task, task_data, +[](gpointer task_data_ptr){
        GeneratorTaskData* task_data = (GeneratorTaskData*) task_data_ptr;
        task_data->cleanup();
        g_slice_free(GeneratorTaskData, task_data_ptr);
    });

    g_timeout_add_once(5000, +[](gpointer user_data){
        g_warning("5 seconds elapsed.");
        g_task_return_boolean(G_TASK(user_data), true);
    }, task);
    
    //this segfaults:
    //g_object_unref(dummy);
    
    // Cannnot deallocate task here because it's still runnning.
    // Where should I do this?

    return cancellable;
    
}