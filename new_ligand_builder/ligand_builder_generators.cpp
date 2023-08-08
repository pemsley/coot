#include "ligand_builder_generators.hpp"
#include "ligand_builder_state.hpp"
#include <glib.h>

GCancellable* coot::ligand_editor::run_generator_request(GeneratorRequest request) {
    GObject* dummy = (GObject*)g_object_new(G_TYPE_OBJECT, NULL);
    GCancellable* cancellable = g_cancellable_new();

    auto task_callback = +[](GObject* obj, GAsyncResult* res, gpointer user_data){
        g_warning("Task callback!");
        auto* progress_dialog = gtk_builder_get_object(layla_gtk_builder, "layla_generator_progress_dialog");
        gtk_window_close(GTK_WINDOW(progress_dialog));

        g_object_unref(obj);        
    };


    GTask* task = g_task_new(dummy,cancellable,task_callback, nullptr);

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