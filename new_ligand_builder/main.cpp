#include <gtk/gtk.h>
#include <thread>
#include "ligand_builder_state.hpp"
#include "ligand_builder_generators.hpp"
#include "ligand_builder_ui.hpp"
#include "../src/coot-setup-python.hh"

int main(int argc, char** argv) {
    using namespace coot::ligand_editor;

    std::thread python_init_thread([argc,argv](){
        setup_python_basic(argc, argv);
        setup_python_coot_module();
        PyRun_SimpleString("import coot_utils");
    });

    python_init_thread.detach();

    gtk_init();
    
    GtkApplication* app = gtk_application_new("org.pemsley.Layla",G_APPLICATION_DEFAULT_FLAGS);
    GError *error = NULL;
    g_application_register(G_APPLICATION(app), NULL, &error);

    g_signal_connect(app,"activate",G_CALLBACK(+[](GtkApplication* app, gpointer user_data){
        // todo: Make this not use a relative path:
        GtkBuilder* builder = gtk_builder_new_from_file("layla.ui");
        
        auto* win = coot::ligand_editor::setup_main_window(app,builder);

        // for now
        auto* icon_theme = gtk_icon_theme_get_for_display(gtk_widget_get_display(GTK_WIDGET(win)));
        gtk_icon_theme_add_search_path(icon_theme, "icons");
        
        coot::ligand_editor::global_layla_gtk_builder = builder;
        coot::ligand_editor::global_generator_request_task_cancellable = nullptr;
        
        gtk_window_present(GTK_WINDOW(win));
        gtk_application_add_window(app,GTK_WINDOW(win));
    }),NULL);


    auto ret = g_application_run(G_APPLICATION(app),0,0);
    g_info("Exiting...");
    delete coot::ligand_editor::global_instance;
    return ret;
}