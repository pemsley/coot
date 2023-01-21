#include "ligand-builder.hpp"
#include <gtk/gtk.h>



void build_main_window(GtkWindow* win) {

}

int main() {
    gtk_init();
      
      GtkApplication* app = gtk_application_new("org.pemsley.NewLigandEditor",G_APPLICATION_FLAGS_NONE);
      GError *error = NULL;
      g_application_register(G_APPLICATION(app), NULL, &error);

      g_signal_connect(app,"activate",G_CALLBACK(+[](GtkApplication* app, gpointer user_data){
         //GtkWindow* win = GTK_WINDOW(user_data);
         GtkWidget* win = gtk_application_window_new(app);
         gtk_application_add_window(app,GTK_WINDOW(win));
         gtk_window_set_application(GTK_WINDOW(win),app);
         build_main_window(GTK_WINDOW(win));
         gtk_widget_show(win);

      }),NULL);


      return g_application_run(G_APPLICATION(app),0,0);
}