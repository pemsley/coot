#include "ligand_builder_generators.hpp"
#include "ligand_builder_state.hpp"
#include <glib.h>
#include <gio/gio.h>
#include <thread>
#include <chrono>

void coot::ligand_editor::run_generator_request(GeneratorRequest request) {
    GObject* dummy = (GObject*)g_object_new(G_TYPE_OBJECT, NULL);
    GCancellable* cancellable = g_cancellable_new();

    auto task_callback = +[](GObject* obj, GAsyncResult* res, gpointer user_data){
        g_warning("Task callback!");
        auto* progress_dialog = gtk_builder_get_object(layla_gtk_builder,"layla_generator_progress_dialog");
        gtk_window_close(GTK_WINDOW(progress_dialog));

        g_object_unref(obj);
    };

    GTask* task = g_task_new(dummy,cancellable,task_callback,nullptr);

    auto task_func = +[](GTask* task, GObject* source_object, gpointer task_data, GCancellable* cancellable){
        g_warning("Task func!");
        std::this_thread::sleep_for(std::chrono::seconds(5));
        g_warning("Task func ends!");
    };

    g_task_run_in_thread(task, (GTaskThreadFunc) task_func);
    g_object_unref(task);
    
    //this segfaults:
    //g_object_unref(dummy);

    // todo: cancellable
    
}