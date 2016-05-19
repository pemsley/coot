

// header-here

#ifndef CFC_WIDGETS_HH
#define CFC_WIDGETS_HH

#include <gtk/gtk.h>

#include "cfc.hh"

namespace cfc {
   
   GtkWidget *wrapped_create_cfc_dialog(cfc::extracted_cluster_info_from_python &extracted_cluster_info);
   void on_cfc_water_cluster_button_clicked(GtkButton *button, gpointer user_data);
   void on_cfc_water_cluster_structure_button_clicked(GtkButton *button, gpointer user_data);

   void on_cfc_pharmacophore_cluster_button_clicked(GtkButton *button, gpointer user_data);
   void on_cfc_pharmacophore_cluster_structure_button_clicked(GtkButton *button, gpointer user_data);

   void wrapped_create_cfc_dialog_add_waters(extracted_cluster_info_from_python &extracted_cluster_info,
					     GtkWidget *cfc_dialog);
   void wrapped_create_cfc_dialog_add_pharmacophores(extracted_cluster_info_from_python &extracted_cluster_info,
					     GtkWidget *cfc_dialog);
   
}

#endif // CFC_WIDGETS_HH

