
#include <iostream>
#include <gtk/gtk.h>

void on_draw(GtkWidget *w, gpointer user_data) {
   std::cout << "draw the data here " << std::endl;
}

void quit(GtkWidget *w, gpointer user_data) {
   std::cout << "bye" << std::endl;
}

int main(int argc, char *argv[]) {

   gtk_init();

   GtkApplication *app = gtk_application_new ("org.emsley.coot",
      (GApplicationFlags) (G_APPLICATION_NON_UNIQUE));

   GtkWidget *window = gtk_window_new();
   gtk_window_set_title(GTK_WINDOW(window), "2D Array Plot");
   gtk_window_set_default_size(GTK_WINDOW(window), 800, 600);

   GtkWidget *drawing_area = gtk_drawing_area_new();
   gtk_drawing_area_set_content_width(GTK_DRAWING_AREA(drawing_area), 512);
   gtk_drawing_area_set_content_height(GTK_DRAWING_AREA(drawing_area), 512);
   g_signal_connect(drawing_area, "draw", G_CALLBACK(on_draw), NULL);

   gtk_window_set_child(GTK_WINDOW(window), drawing_area);

   g_signal_connect(window, "destroy", G_CALLBACK(quit), NULL);

   gtk_widget_set_visible(window, TRUE);

   int status = g_application_run(G_APPLICATION(app), 1, argv);

   return status;
}
