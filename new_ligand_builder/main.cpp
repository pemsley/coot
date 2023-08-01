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
        // GtkBuilder* builder = gtk_builder_new_from_file("layla.ui");
        // g_object_unref(builder);
        //GtkWindow* win = GTK_WINDOW(user_data);
        GtkWidget* win = gtk_application_window_new(app);
        gtk_window_set_title(GTK_WINDOW(win),"New Ligand Editor");
        gtk_application_window_set_show_menubar(GTK_APPLICATION_WINDOW(win), TRUE);
        gtk_window_set_application(GTK_WINDOW(win),app);
        auto* canvas = coot_ligand_editor_canvas_new();
        GtkWidget* status_label = gtk_label_new("");
        coot::ligand_editor::initialize_global_instance(canvas,GTK_WINDOW(win),GTK_LABEL(status_label));
        gtk_application_set_menubar(app, G_MENU_MODEL(build_menu(app,canvas,GTK_WINDOW(win))));
        gtk_application_add_window(app,GTK_WINDOW(win));
        build_main_window(GTK_WINDOW(win),canvas,GTK_LABEL(status_label));
        gtk_window_present(GTK_WINDOW(win));

    }),NULL);


    auto ret = g_application_run(G_APPLICATION(app),0,0);
    g_info("Exiting...");
    delete coot::ligand_editor::global_instance;
    return ret;
}