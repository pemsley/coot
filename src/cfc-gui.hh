
#ifndef COOT_SRC_CFC_GUI_HH
#define COOT_SRC_CFC_GUI_HH

#include <gtk/gtk.h>
#include "coot-utils/cfc.hh"

class cfc_gui_t {
public:
   std::vector<cfc::typed_cluster_t> cluster_infos;
   std::vector<std::vector<cfc::water_info_t> > water_infos;
   std::string style_css;
   GtkBuilder *builder;
   GtkWidget *widget_from_builder(const std::string &s);
   cfc_gui_t() : builder(nullptr) {}
   void setup();
   GtkWidget *get_dialog() {
      if (!builder)
	 setup();
      return widget_from_builder("cfc-dialog");
   }
   // pre-sort indices
   void fill_ligands_grid(const std::vector<int> &generic_object_indices_for_features);
   void fill_waters_grid(const std::vector<int> &generic_object_indices_for_waters);
};

#endif // COOT_SRC_CFC_GUI_HH
