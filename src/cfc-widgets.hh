

// header-here

#ifndef CFC_WIDGETS_HH
#define CFC_WIDGETS_HH

#include <gtk/gtk.h>

#include "cfc.hh"

namespace cfc {
   
   GtkWidget *wrapped_create_cfc_dialog(const cfc::extracted_cluster_info_from_python &extracted_cluster_info);
   void on_cfc_water_cluster_button_clicked(GtkButton *button, gpointer user_data);
   
}

#endif // CFC_WIDGETS_HH

