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
        // gtk_icon_theme_add_search_path(GtkIconTheme *self, const char *path)
        // todo: Make this not use a relative path:
        GtkBuilder* builder = gtk_builder_new_from_file("layla.ui");
        
        auto* win = coot::ligand_editor::setup_main_window(app,builder);
        coot::ligand_editor::layla_gtk_builder = builder;
        
        gtk_window_present(GTK_WINDOW(win));
        gtk_application_add_window(app,GTK_WINDOW(win));
    }),NULL);


    auto ret = g_application_run(G_APPLICATION(app),0,0);
    g_info("Exiting...");
    delete coot::ligand_editor::global_instance;
    return ret;
}