
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
   std::vector<int> imols_in_cluster_vec; // this is used for vertical sorting
   cfc_gui_t() : builder(nullptr) {}
   void setup();
   GtkWidget *get_dialog() {
      if (!builder)
	 setup();
      return widget_from_builder("cfc-dialog");
   }
   // store these so that they can be used when the ligand grid
   // gets sorted differently (and thus regenerated)
   std::vector<int> generic_object_indices_for_features;
   std::vector<int> generic_object_indices_for_waters;
   std::vector<int> generic_object_indices_for_contributors;
   // pre-sort indices
   void fill_ligands_grid();
   void fill_waters_grid();
   static void toggle_molecule_buttons(GtkToggleButton *toggle_button, int imol, int state);

   void sort_cluster_info_internals_by_mol_no();
   void sort_cluster_info_internals_by_number_of_chemical_features();

   void set_generic_object_indices_for_features(const std::vector<int> &generic_object_indices_for_features_in) {
     generic_object_indices_for_features = generic_object_indices_for_features_in;
   }
   void set_generic_object_indices_for_waters(const std::vector<int> &generic_object_indices_for_waters_in) {
     generic_object_indices_for_waters = generic_object_indices_for_waters_in;
   }

   // utility function
   void clear_grid(GtkWidget *grid);
   // when we rebuild the grid, it would be good if the displayed imols toggled buttons
   // are reproduced in the new grid - so let's save their state using this function
   // beefore the grid is deleted.
   std::set<int> get_imols_that_are_displayed(GtkWidget *grid);
   // note that fill_ligands_grid() doesn't take an argument, so we need to store the
   // above so that it can be used in fill_ligands_grid();
   std::set<int> imols_that_are_displayed;
   // likewise the chemical features (the vertical list on the left)
   std::set<std::string> get_chemical_features_that_are_displayed(GtkWidget *grid);
   std::set<std::string> chemical_features_that_are_displayed;

};

#endif // COOT_SRC_CFC_GUI_HH
