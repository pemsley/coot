
#include <string>
#include "coot-utils/atom-selection-container.hh"
#include "sequence-view-widget.hh"

void fill(GtkWidget *window, mmdb::Manager *mol) {

   GtkWidget *scrolled_window = gtk_scrolled_window_new();
   GtkWidget *box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 2);
   GtkWidget *frame = gtk_frame_new("");
   gtk_widget_set_hexpand(scrolled_window, TRUE);
   gtk_widget_set_vexpand(scrolled_window, TRUE);
   gtk_widget_set_size_request(scrolled_window, 200, 240);
   gtk_widget_set_size_request(frame, 200, 250); // h size be bigger than the h-size for the scrolled window
   gtk_window_set_child(GTK_WINDOW(window), scrolled_window);
   CootSequenceView *sv = coot_sequence_view_new();
   int imol = 0; // for now
   coot_sequence_view_set_structure(sv, imol, mol);
   g_object_set_data(G_OBJECT(sv), "sv3-frame", frame);
   gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scrolled_window), GTK_WIDGET(frame));

   // gtk_box_append(GTK_BOX(box), GTK_WIDGET(sv));
   gtk_frame_set_child(GTK_FRAME(frame), GTK_WIDGET(sv));

   auto callback = +[] (CootSequenceView* self,
                        const box_info_t *residue_vip,
                        gpointer userdata) {

      std::cout << "residue-clicked handler" << std::endl;
   };

   g_signal_connect(sv, "residue-clicked", G_CALLBACK(callback), nullptr);

}

int main(int argc, char **argv) {

   int status = 0;
   if (argc > 1) {
      std::string pdb_file_name = argv[1];
      auto atom_sel = get_atom_selection(pdb_file_name, true, false, false);
      if (atom_sel.read_success) {
         gtk_init();
         GtkApplication* app = gtk_application_new("org.pemsley.test-sequence-view", G_APPLICATION_DEFAULT_FLAGS);
         GError *error = NULL;
         g_application_register(G_APPLICATION(app), NULL, &error);
         if (error) {
            std::cout << "ERROR" << error->message << std::endl;
         } else {
            // let's go
            auto activate_callback = +[] (GtkApplication *app, gpointer user_data) {
               mmdb::Manager *mol = static_cast<mmdb::Manager *>(user_data);
               GtkWidget *win = gtk_application_window_new(app);
               gtk_application_add_window(app, GTK_WINDOW(win));
               gtk_window_set_application(GTK_WINDOW(win), app);
               gtk_window_set_title(GTK_WINDOW(win), "Coot: Sequence View");
               int width  = 500;
               int height = 100;
               gtk_window_set_default_size(GTK_WINDOW(win), width, height);
               fill(win, mol);
               gtk_widget_set_visible(win, TRUE);
            };
            g_signal_connect(app, "activate", G_CALLBACK(activate_callback), atom_sel.mol);

            // Ideally the arguments should be argc, argv and
            // call above g_application_new with G_APPLICATION_HANDLES_OPEN
            // and connect a signal to "open":
            // g_signal_connect(app, "open", G_CALLBACK(open_callback), activate_data).
            return g_application_run(G_APPLICATION(app), 0, 0);
         }
      
      } else {
         status = 1; // bad input
      }
   } else {
      status = 1; // bad input
   }

   return status;
}
