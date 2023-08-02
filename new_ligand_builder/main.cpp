#include <gtk/gtk.h>
#include "ligand_builder_state.hpp"
#include "ligand_builder_ui.hpp"

int main() {
    using namespace coot::ligand_editor;

    gtk_init();
    
    GtkApplication* app = gtk_application_new("org.pemsley.Layla",G_APPLICATION_DEFAULT_FLAGS);
    GError *error = NULL;
    g_application_register(G_APPLICATION(app), NULL, &error);

    g_signal_connect(app,"activate",G_CALLBACK(+[](GtkApplication* app, gpointer user_data){
        // todo: Make this not use a relative path:
        GtkBuilder* builder = gtk_builder_new_from_file("layla.ui");
        GtkApplicationWindow* win = (GtkApplicationWindow*) gtk_builder_get_object(builder, "layla_window");
        gtk_window_set_application(GTK_WINDOW(win),app);
        GtkWidget* status_label = (GtkWidget*) gtk_builder_get_object(builder, "layla_status_label");
        GtkScrolledWindow* viewport = (GtkScrolledWindow*) gtk_builder_get_object(builder, "layla_canvas_viewport");
        auto* canvas = coot_ligand_editor_canvas_new();
        gtk_scrolled_window_set_child(viewport, GTK_WIDGET(canvas));
        coot::ligand_editor::initialize_global_instance(canvas,GTK_WINDOW(win),GTK_LABEL(status_label));
        coot::ligand_editor::setup_actions(win, canvas, builder);
        gtk_window_present(GTK_WINDOW(win));
        gtk_application_add_window(app,GTK_WINDOW(win));
        g_object_unref(builder);
    }),NULL);


    auto ret = g_application_run(G_APPLICATION(app),0,0);
    g_info("Exiting...");
    delete coot::ligand_editor::global_instance;
    return ret;
}