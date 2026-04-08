
#include <filesystem>

#include "cc-interface.hh"
#include "coot-utils/cfc.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "glib-object.h"
#include "glib.h"
#include "glibconfig.h"
#include "gtk/gtk.h"
#include "gtk/gtkshortcut.h"
#include "utils/coot-utils.hh"
#include "cfc-gui.hh"
#include "graphics-info.h"
#include "c-interface.h"
#include "c-interface-generic-objects.h"

void
cfc_gui_t::setup() {

   builder = gtk_builder_new();
   if (GTK_IS_BUILDER(builder)) {
      GError* error = NULL;
      std::filesystem::path pkg_data_dir = coot::package_data_dir();
      std::filesystem::path ui_dir = pkg_data_dir / "ui";
      std::filesystem::path cfc_ui_path = ui_dir / "cfc.ui";
      std::string ui_file_full = cfc_ui_path.string();
      gboolean status = gtk_builder_add_from_file(builder, ui_file_full.c_str(), &error);
      if (status == FALSE) {
         std::cout << "ERROR:: Failure to read or parse " << ui_file_full << std::endl;
         std::cout << error->message << std::endl;
      }
   } else {
      std::cout << "ERROR:: cfc bad builder" << std::endl;
   }

   style_css = ".custom-cfc-button {         \
                               background: #577040;        /* Green */ \
                               color: #aaaaaa;        \
                               padding: 0;         \
                               margin: 0;          \
                               min-width: 12px;   /* adjust as needed */ \
                               border-radius: 0px; \
                             } \
                .custom-cfc-button:checked{         \
                               background: #274020;        /* Green */ \
                               color: white;       \
                               padding: 0;         \
                               margin: 0;          \
                               min-width: 12px;   /* adjust as needed */ \
                               border-radius: 0px; \
                             } \
                .custom-cfc-blank-button {   \
                               background: #404040;        /* Green */ \
                               padding: 0;         \
                               color: #404040;     \
                               margin: 0;          \
                               min-width: 12px;     \
                               border-radius: 0px; \
                            } \
                .custom-cfc-blank-button:checked {   \
                               background: #303030;        /* Green */ \
                               padding: 0;         \
                               color: #505050;       \
                               margin: 0;          \
                               min-width: 12px;     \
                               border-radius: 0px; \
                            } \
                            ";
#if 0 // 20260122-PE remove compiler warnings
   // style_css = "button.custom-cfc-button { \
   //                             background: #576F50;        /* Green */ \
   //                             color: white;       \
   //                             padding: 0;         \
   //                             margin: 0;          \
   //                             min-width: 3px;   /* adjust as needed */ \
   //                             border-radius: 3px; \
   //                          } \
   //                         button.custom-cfc-blank-button { \
   //                             padding: 0;         \
   //                             color: #202020;     \
   //                             margin: 0;          \
   //                             min-width: 3px;     \
   //                             border-radius: 1px; \
   //                          } \
   //                          ";
#endif
}

GtkWidget *
cfc_gui_t::widget_from_builder(const std::string &widget_name) {

   GtkWidget *w = nullptr;
   if (! builder)
      setup();
   if (builder)
      w = GTK_WIDGET(gtk_builder_get_object(builder, widget_name.c_str()));
   return w;
}

void
cfc_gui_t::sort_cluster_info_internals_by_mol_no() {

   // we sort imols_in_cluster_vec, because that is what is used
   // to fill the "grid" of buttons

   std::sort(imols_in_cluster_vec.begin(), imols_in_cluster_vec.end());

}

void
cfc_gui_t::sort_cluster_info_internals_by_number_of_chemical_features() {

   // a mapping between molecule number and number of chemical features
   std::map<int, unsigned int> number_of_chemical_features;
   // fill it...
   for (const auto &ci : cluster_infos) {
      for (const auto &p : ci.imols_with_specs) {
         int imol = p.first;
         if (number_of_chemical_features.find(imol) == number_of_chemical_features.end())
            number_of_chemical_features[imol] = 1;
         else
            number_of_chemical_features[imol] += 1;
      }
   }

   auto sorter = [number_of_chemical_features] (int imol_1, int imol_2) {
      return number_of_chemical_features.at(imol_2) < number_of_chemical_features.at(imol_1);
   };

   std::sort(imols_in_cluster_vec.begin(), imols_in_cluster_vec.end(), sorter);

}


extern "C" G_MODULE_EXPORT
void
on_cfc_ligands_all_on_button_clicked(GtkButton       *button,
                                     gpointer         user_data) {

   graphics_info_t g;
   int n = g.n_molecules();
   for (int i=0; i<n; i++) {
      if (is_valid_model_molecule(i)) {
         if (g.molecules[i].is_displayed_p()) {
         } else {
            g.molecules[i].set_mol_is_displayed(true);
            set_display_control_button_state(i, "Displayed", true);
         }
      }
   }

   // now toggle the buttons of the chemical features

   GtkWidget *grid = g.cfc_gui.widget_from_builder("cfc-ligands-grid");
   if (grid) {
      int row = 0;
      GtkWidget *feature_toggle_button = gtk_grid_get_child_at(GTK_GRID(grid), 0, row);
      do {
         if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(feature_toggle_button))) {
         } else {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(feature_toggle_button), TRUE);
         }
         row++;
         feature_toggle_button = gtk_grid_get_child_at(GTK_GRID(grid), 0, row);
      } while(feature_toggle_button);
   }

   graphics_draw();
}

extern "C" G_MODULE_EXPORT
void
on_cfc_ligands_all_off_button_clicked(GtkButton       *button,
                                     gpointer         user_data) {

   graphics_info_t g;
   int n = g.n_molecules();
   for (int i=0; i<n; i++) {
      if (is_valid_model_molecule(i)) {
         if (g.molecules[i].is_displayed_p()) {
            g.molecules[i].set_mol_is_displayed(false);
            set_display_control_button_state(i, "Displayed", false);
         }
      }
   }

   // now toggle the buttons of the chemical features

   GtkWidget *grid = g.cfc_gui.widget_from_builder("cfc-ligands-grid");
   if (grid) {
      int row = 0;
      GtkWidget *feature_toggle_button = gtk_grid_get_child_at(GTK_GRID(grid), 0, row);
      do {
         if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(feature_toggle_button))) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(feature_toggle_button), FALSE);
         }
         row++;
         feature_toggle_button = gtk_grid_get_child_at(GTK_GRID(grid), 0, row);
      } while(feature_toggle_button);
   }

   graphics_draw();
}

extern "C" G_MODULE_EXPORT
void
on_cfc_ligands_show_chemical_features_button_clicked(GtkButton       *button,
                                                     gpointer         user_data) {

   graphics_info_t g;
   std::cout << "debug:: on_cfc_ligands_show_chemical_features_button_clicked() showing "
             << g.cfc_gui.generic_object_indices_for_features.size()
             << " generic objects" << std::endl;
   for (unsigned int i=0; i<g.cfc_gui.generic_object_indices_for_features.size(); i++) {
      int idx = g.cfc_gui.generic_object_indices_for_features[i];
      set_display_generic_object(idx, 1);
   }
}

extern "C" G_MODULE_EXPORT
void
on_cfc_ligands_hide_chemical_features_button_clicked(GtkButton       *button,
                                                     gpointer         user_data) {

   graphics_info_t g;
   std::cout << "debug:: on_cfc_ligands_hide_chemical_features_button_clicked() hiding "
             << g.cfc_gui.generic_object_indices_for_features.size()
             << " generic objects" << std::endl;
   for (unsigned int i=0; i<g.cfc_gui.generic_object_indices_for_features.size(); i++) {
      int idx = g.cfc_gui.generic_object_indices_for_features[i];
      set_display_generic_object(idx, 0);
   }
}

extern "C" G_MODULE_EXPORT
void on_cfc_hide_contributor_objects_button_clicked(GtkButton *button,
                                                    gpointer   user_data) {

   graphics_info_t g;

   std::cout << "debug:: in on_cfc_waters_hide_contributor_objects_button_clicked() "
             << g.cfc_gui.generic_object_indices_for_contributors.size()
             << " objects for contributors" << std::endl;

   for (unsigned int i=0; i<g.cfc_gui.generic_object_indices_for_contributors.size(); i++) {
      int idx = g.cfc_gui.generic_object_indices_for_contributors[i];
      set_display_generic_object(idx, 0);
   }
}


extern "C" G_MODULE_EXPORT
void
on_cfc_waters_all_on_button_clicked(GtkButton       *button,
                                    gpointer         user_data) {

   graphics_info_t g;
   int n = g.n_molecules();
   for (int i=0; i<n; i++) {
      if (is_valid_model_molecule(i)) {
         if (g.molecules[i].is_displayed_p()) {
         } else {
            g.molecules[i].set_mol_is_displayed(true);
            set_display_control_button_state(i, "Displayed", true);
         }
      }
   }

   // now toggle the buttons of the water clusters

   GtkWidget *grid = g.cfc_gui.widget_from_builder("cfc-waters-grid");
   if (grid) {
      int row = 0;
      GtkWidget *water_toggle_button = gtk_grid_get_child_at(GTK_GRID(grid), 0, row);
      do {
         if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(water_toggle_button))) {
         } else {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(water_toggle_button), TRUE);
         }
         row++;
         water_toggle_button = gtk_grid_get_child_at(GTK_GRID(grid), 0, row);
      } while(water_toggle_button);
   }

   graphics_draw();

}

extern "C" G_MODULE_EXPORT
void
on_cfc_waters_all_off_button_clicked(GtkButton       *button,
                                     gpointer         user_data) {

   graphics_info_t g;
   int n = g.n_molecules();
   for (int i=0; i<n; i++) {
      if (is_valid_model_molecule(i)) {
         if (g.molecules[i].is_displayed_p()) {
            g.molecules[i].set_mol_is_displayed(false);
            set_display_control_button_state(i, "Displayed", false);
         }
      }
   }

   // now toggle the buttons of the water clusters

   GtkWidget *grid = g.cfc_gui.widget_from_builder("cfc-waters-grid");
   if (grid) {
      int row = 0;
      GtkWidget *water_toggle_button = gtk_grid_get_child_at(GTK_GRID(grid), 0, row);
      do {
         if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(water_toggle_button))) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(water_toggle_button), FALSE);
         }
         row++;
         water_toggle_button = gtk_grid_get_child_at(GTK_GRID(grid), 0, row);
      } while(water_toggle_button);
   }

   graphics_draw();
}

extern "C" G_MODULE_EXPORT
void
on_cfc_waters_show_generic_objects_button_clicked(GtkButton       *button,
                                                  gpointer         user_data) {

   graphics_info_t g;
   for (unsigned int i=0; i<g.cfc_gui.generic_object_indices_for_waters.size(); i++) {
      int idx = g.cfc_gui.generic_object_indices_for_waters[i];
      set_display_generic_object(idx, 1);
   }
}

extern "C" G_MODULE_EXPORT
void
on_cfc_waters_hide_generic_objects_button_clicked(GtkButton       *button,
                                                  gpointer         user_data) {

   graphics_info_t g;
   for (unsigned int i=0; i<g.cfc_gui.generic_object_indices_for_waters.size(); i++) {
      int idx = g.cfc_gui.generic_object_indices_for_waters[i];
      set_display_generic_object(idx, 0);
   }
}

extern "C" G_MODULE_EXPORT
void
on_cfc_waters_hide_contributor_objects_button_clicked(GtkButton       *button,
                                                      gpointer         user_data) {
   // this is no signal for this atm
}

std::set<int>
cfc_gui_t::get_imols_that_are_displayed(GtkWidget *grid) {

   std::set<int> s;
   if (grid) {
      int row = 0;
      GtkWidget *feature_hbox = gtk_grid_get_child_at(GTK_GRID(grid), 1, row);
      do {
         if (feature_hbox) {
            GtkWidget *tb = gtk_widget_get_first_child(feature_hbox);
            while (tb) {
               if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(tb))) {
                  int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(tb), "imol"));
                  s.insert(imol);
               }
               tb = gtk_widget_get_next_sibling(tb);
            }
         }
         row++;
         feature_hbox = gtk_grid_get_child_at(GTK_GRID(grid), 1, row);
      } while(feature_hbox);

   }
   return s;
}

// likewise the chemical features (the vertical list on the left)
std::set<std::string>
cfc_gui_t::get_chemical_features_that_are_displayed(GtkWidget *grid) {

   std::set<std::string> s;
   if (grid) {
      int row = 0;
      GtkWidget *feature_button = gtk_grid_get_child_at(GTK_GRID(grid), 0, row);
      do {
         if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(feature_button))) {
            std::string l = gtk_button_get_label(GTK_BUTTON(feature_button));
            s.insert(l);
         }
         row++;
         feature_button = gtk_grid_get_child_at(GTK_GRID(grid), 0, row);
      } while(feature_button);
   }
   return s;
}

void
cfc_gui_t::clear_grid(GtkWidget *grid) {

   // for future reference, this is not a "proper" grid, it is just
   // two column - the second of which is an hbox.
   // I forget now why I didn't want to use a grid.

   if (grid) {

      // save state before deleting
      imols_that_are_displayed = get_imols_that_are_displayed(grid);
      chemical_features_that_are_displayed = get_chemical_features_that_are_displayed(grid);

      int row = 0;
      GtkWidget *feature_toggle_button = gtk_grid_get_child_at(GTK_GRID(grid), 0, row);
      GtkWidget *feature_hbox          = gtk_grid_get_child_at(GTK_GRID(grid), 1, row);
      do {
         if (feature_hbox) {
            GtkWidget *tb = gtk_widget_get_first_child(feature_hbox);
            while (tb) {
               // gtk_box_remove(GTK_BOX(feature_hbox), tb);
               tb = gtk_widget_get_next_sibling(tb);
            }
         }

         gtk_grid_remove(GTK_GRID(grid), feature_toggle_button);
         gtk_grid_remove(GTK_GRID(grid), feature_hbox);
         row++;
         feature_toggle_button = gtk_grid_get_child_at(GTK_GRID(grid), 0, row);
         feature_hbox          = gtk_grid_get_child_at(GTK_GRID(grid), 1, row);
      } while(feature_toggle_button);
   }
}

extern "C" G_MODULE_EXPORT
void
on_cfc_ligands_sort_by_molecule_number(GtkButton *button,
                                       gpointer   user_data) {

   graphics_info_t g;
   g.cfc_gui.sort_cluster_info_internals_by_mol_no();
   // regen the ligand grid
   GtkWidget *grid = g.cfc_gui.widget_from_builder("cfc-ligands-grid");
   g.cfc_gui.clear_grid(grid); // stash imols_that_are_displayed
   g.cfc_gui.fill_ligands_grid();

}

extern "C" G_MODULE_EXPORT
void
on_cfc_ligands_sort_by_number_of_chemical_features(GtkButton *button,
                                                   gpointer   user_data) {

   std::cout << "sort by number of chemical features" << std::endl;
   graphics_info_t g;
   g.cfc_gui.sort_cluster_info_internals_by_number_of_chemical_features();
   // regen the ligand grid
   GtkWidget *grid = g.cfc_gui.widget_from_builder("cfc-ligands-grid");
   g.cfc_gui.clear_grid(grid);  // stash imols_that_are_displayed
   std::cout << "calling fill_ligands_grid()" << std::endl;
   g.cfc_gui.fill_ligands_grid();

}

void
cfc_gui_t::fill_waters_grid() {

   auto sorter = +[] (const std::vector<cfc::water_info_t> &v1,
                      const std::vector<cfc::water_info_t> &v2) {
      return v1.size() > v2.size();
   };

   auto imol_is_part_of_cluster = [] (int imol, const std::set<int> &imols) {
      return imols.find(imol) != imols.end();
   };

   class water_cluster_t {
   public:
      water_cluster_t() : imol(-1), pos(0.0, 0.0, 0.0) {}
      water_cluster_t(int i, const RDGeom::Point3D &p) : imol(i), pos(p) {}
      int imol;
      RDGeom::Point3D pos;
   };

   GtkWidget *grid = widget_from_builder("cfc-waters-grid");
   if (grid) {

      // std::sort(water_infos.begin(), water_infos.end(), sorter);

      // Load CSS
      GtkCssProvider *provider = gtk_css_provider_new();
      gtk_css_provider_load_from_string(provider, style_css.c_str());
      gtk_style_context_add_provider_for_display(gdk_display_get_default(),
                                                 GTK_STYLE_PROVIDER(provider),
                                                 GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

      std::set<int> all_the_imols;
      for (unsigned int i=0; i<water_infos.size(); i++) {
         const auto &wi = water_infos[i];
         for (unsigned int jj=0; jj<wi.size(); jj++) {
            int imol = wi[jj].imol;
            all_the_imols.insert(imol);
         }
      }

      auto callback = +[] (GtkToggleButton* togglebutton, gpointer data) {
         std::cout << "water cluster toggle-button toggled" << std::endl;
         water_cluster_t *wc = static_cast<water_cluster_t *>(data);
         int state = 0;
         if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
            state = 1;
         } else {
            state = 0;
         }
         clipper::Coord_orth pt(wc->pos.x, wc->pos.y, wc->pos.z);
         // don't recentre if we are merely turnning off the representation
         if (state == 1)
            graphics_info_t::set_rotation_centre(pt);
         int godi = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(togglebutton),
                                                      "generic-display-object-index"));
         if (state == 1)
            set_display_generic_object_simple(godi, 1);
         else
            set_display_generic_object_simple(godi, 0);
      };

      auto tb_callback = +[] (GtkToggleButton *tb, gpointer data) {
         int imol = GPOINTER_TO_INT(data);
         if (gtk_toggle_button_get_active(tb))
            set_mol_displayed(imol, 1);
         else
            set_mol_displayed(imol, 0);
      };

      for (unsigned int i=0; i<water_infos.size(); i++) {
         const auto &wi = water_infos[i];

         std::set<int> imols_in_cluster;
         for (unsigned int jj=0; jj<wi.size(); jj++) {
            int imol = wi[jj].imol;
            imols_in_cluster.insert(imol);
         }

         GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
         gtk_grid_attach(GTK_GRID(grid), hbox, 1, i, 1,1);

         std::string label = "Water Cluster " + std::to_string(i);
         GtkWidget *tb = gtk_toggle_button_new_with_label(label.c_str());
         gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tb), TRUE);
         int generic_dislay_object_index = -1;
         if (i < generic_object_indices_for_waters.size())
            generic_dislay_object_index = generic_object_indices_for_waters[i];
         g_object_set_data(G_OBJECT(tb),
                           "generic-display-object-index",
                           GINT_TO_POINTER(generic_dislay_object_index));
         water_cluster_t *wc = new water_cluster_t(i, wi[0].pos);
         g_signal_connect(G_OBJECT(tb), "toggled", G_CALLBACK(callback), wc);
         gtk_grid_attach(GTK_GRID(grid), tb, 0, i, 1, 1);

         std::set<int>::const_iterator it;
         for (it=all_the_imols.begin(); it!=all_the_imols.end(); ++it) {
            int imol = *it;
            GtkWidget *b = gtk_toggle_button_new_with_label("~");
            g_signal_connect(G_OBJECT(b), "toggled", G_CALLBACK(tb_callback),
                             GINT_TO_POINTER(imol));
            GtkStyleContext *context = gtk_widget_get_style_context(b);
            if (imol_is_part_of_cluster(imol, imols_in_cluster))
               gtk_style_context_add_class(context, "custom-cfc-button");
            else
               gtk_style_context_add_class(context, "custom-cfc-blank-button");
            gtk_box_append(GTK_BOX(hbox), b);
         }
      }
   }
}

// static
void
toggle_molecule_buttons(GtkToggleButton *toggle_button, int imol, int state) {

   std::cout << "toggle_molecule_buttons called" << std::endl;

   GtkWidget *grid = GTK_WIDGET(g_object_get_data(G_OBJECT(toggle_button), "grid"));
   if (grid) {

   }
};

void
cfc_gui_t::fill_ligands_grid() {

   auto sorter = +[] (const cfc::typed_cluster_t &t1, const cfc::typed_cluster_t &t2) {
      return t1.imols_with_specs.size() > t2.imols_with_specs.size();
   };


   GtkWidget *grid = widget_from_builder("cfc-ligands-grid");
   if (grid) {

      // a set of imols, implicitly sorted when used.
      std::set<int> imols_in_cluster_set;
      for (unsigned int i=0; i<cluster_infos.size(); i++) {
         const auto &ci = cluster_infos[i];
         for (unsigned int jj=0; jj<ci.imols_with_specs.size(); jj++)
            imols_in_cluster_set.insert(ci.imols_with_specs[jj].first);
      }

      // We want to know the state of the buttons. Do we know?
      // The first time round, every button should be "on"/toggled.
      bool first_time = false;
      // if imols_in_cluster_vec is empty, this is the first time round
      if (imols_in_cluster_vec.empty()) {
         first_time = true;
         std::set<int>::const_iterator its;
         for (its=imols_in_cluster_set.begin(); its!=imols_in_cluster_set.end(); its++) {
            imols_in_cluster_vec.push_back(*its);
         }
      }

      // this is the "vertical" sort.
      std::vector<cfc::typed_cluster_t> sorted_cluster_infos = cluster_infos;
      std::sort(sorted_cluster_infos.begin(), sorted_cluster_infos.end(), sorter);

      // Load CSS
      GtkCssProvider *provider = gtk_css_provider_new();
      // gtk_css_provider_load_from_path(provider, "style.css");
      gtk_css_provider_load_from_string(provider, style_css.c_str());

      gtk_style_context_add_provider_for_display(gdk_display_get_default(),
                                                 GTK_STYLE_PROVIDER(provider),
                                                 GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);

      for (unsigned int i=0; i<sorted_cluster_infos.size(); i++) {
         const auto &ci = sorted_cluster_infos[i];

         GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
         gtk_grid_attach(GTK_GRID(grid), hbox, 1, i, 1,1);

         std::string label = ci.family + " " + ci.type + " " + std::to_string(ci.idx);
         GtkWidget *tb = gtk_toggle_button_new_with_label(label.c_str());
         if (chemical_features_that_are_displayed.find(label) != chemical_features_that_are_displayed.end())
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tb), TRUE);
         else
            if (first_time)
               gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tb), TRUE);
         int generic_dislay_object_index = -1;
         if (i < generic_object_indices_for_features.size())
            generic_dislay_object_index = generic_object_indices_for_features[i];
         g_object_set_data(G_OBJECT(tb),
                           "generic-display-object-index",
                           GINT_TO_POINTER(generic_dislay_object_index));
         gtk_grid_attach(GTK_GRID(grid), tb, 0, i, 1, 1);

         auto callback = +[] (GtkToggleButton* togglebutton, gpointer data) {

            cfc::typed_cluster_t *tc = static_cast<cfc::typed_cluster_t *>(data);
            int state = 0;
            if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
               state = 1;
            } else {
               state = 0;
            }
            clipper::Coord_orth pt(tc->pos.x, tc->pos.y, tc->pos.z);
            // don't recentre if we are merely turning off the representation
            if (state == 1)
               graphics_info_t::set_rotation_centre(pt);
            for (unsigned int ii=0; ii<tc->imols_with_specs.size(); ii++) {
               int imol = tc->imols_with_specs[ii].first;
               set_mol_displayed(imol, state);
            }
            int godi = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(togglebutton),
                                                         "generic-display-object-index"));
            if (state == 1)
               set_display_generic_object_simple(godi, 1);
            else
               set_display_generic_object_simple(godi, 0);

            // now we need to turn off (or on) each of the individual buttons
            GtkWidget *hbox = GTK_WIDGET(g_object_get_data(G_OBJECT(togglebutton), "hbox"));
            // std::cout << "fill_lignds_grid(): callback testing hbox " << hbox << std::endl;
            if (hbox) {
               GtkWidget *item_widget = gtk_widget_get_first_child(hbox);
               while (item_widget) {
                  GtkToggleButton *tb = GTK_TOGGLE_BUTTON(item_widget);
                  int imol_is_in_this_cluster =
                     GPOINTER_TO_INT(g_object_get_data(G_OBJECT(tb), "imol-is-in-this-cluster"));
                  if (imol_is_in_this_cluster == 1) {
                     if (state == 0)
                        if (gtk_toggle_button_get_active(tb))
                           gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tb), 0);
                     if (state == 1)
                        if (gtk_toggle_button_get_active(tb) == FALSE)
                           gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tb), 1);
                  }
                  item_widget = gtk_widget_get_next_sibling(item_widget);
               }
            } else {
               std::cout << "ERROR:: missing hbox" << std::endl;
            }
         };

         // the molecule toggle buttons (not the feature toggle buttons)
         auto toggle_button_toggled_callback = +[] (GtkToggleButton *b, gpointer data) {

            int imol = GPOINTER_TO_INT(data);
            int state = 2;
            if (gtk_toggle_button_get_active(b)) {
               set_mol_displayed(imol, 1);
               state = 1;
            } else {
               set_mol_displayed(imol, 0);
               state = 0;
            }
            GtkWidget *grid = GTK_WIDGET(g_object_get_data(G_OBJECT(b), "grid"));
            if (grid) {

               int row = 0;
               while (true) {
                  GtkWidget *box = gtk_grid_get_child_at(GTK_GRID(grid), 1, row);
                  // std::cout << "do loop with box " << box << std::endl;
                  if (box) {
                     box = gtk_grid_get_child_at(GTK_GRID(grid), 1, row);
                     row++; // for next time

                     // the children of this box are a row of toggle buttons
                     GtkWidget *item_widget = gtk_widget_get_first_child(box);
                     while (item_widget) {
                        GtkToggleButton *box_toggle_button = GTK_TOGGLE_BUTTON(item_widget);
                        if (box_toggle_button != b) {
                           int imol_box_toggle_button =
                              GPOINTER_TO_INT(g_object_get_data(G_OBJECT(box_toggle_button), "imol"));
                           if (imol_box_toggle_button == imol) {
                              // don't toggle unless needed
                              if (state == 1) {
                                 if (gtk_toggle_button_get_active(box_toggle_button) == FALSE)
                                    gtk_toggle_button_set_active(box_toggle_button, 1);
                              }
                              if (state == 0) {
                                 if (gtk_toggle_button_get_active(box_toggle_button))
                                    gtk_toggle_button_set_active(box_toggle_button, 0);
                              }
                           }
                        }
                        item_widget = gtk_widget_get_next_sibling(item_widget);
                     }
                  } else {
                     break;
                  }
               }

            } else {
               std::cout << "ERROR:: no grid" << std::endl;
            }

         };

         cfc::typed_cluster_t *tc = new cfc::typed_cluster_t(ci);
         g_signal_connect(G_OBJECT(tb), "toggled", G_CALLBACK(callback), tc);
         generic_dislay_object_index = -1;
         if (i < generic_object_indices_for_features.size())
            generic_dislay_object_index = generic_object_indices_for_features[i];
         g_object_set_data(G_OBJECT(tb), "hbox", hbox);
         g_object_set_data(G_OBJECT(tb),
                           "generic-display-object-index",
                           GINT_TO_POINTER(generic_dislay_object_index));

         std::vector<int>::const_iterator its;
         for (its=imols_in_cluster_vec.begin(); its!=imols_in_cluster_vec.end(); ++its) {
            int imol = *its;
            GtkWidget *b = gtk_toggle_button_new_with_label(".");
            if (first_time)
               gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), TRUE);
            else
               if (imols_that_are_displayed.find(imol) != imols_that_are_displayed.end())
                  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), TRUE);
            int imol_is_in_this_cluster = 0;
            if (ci.imol_is_part_of_cluster(imol)) imol_is_in_this_cluster = 1;
            g_object_set_data(G_OBJECT(b), "imol-is-in-this-cluster", GINT_TO_POINTER(imol_is_in_this_cluster));
            g_object_set_data(G_OBJECT(b), "imol", GINT_TO_POINTER(imol));
            g_object_set_data(G_OBJECT(b), "grid", grid);

            g_signal_connect(G_OBJECT(b), "toggled",
                             G_CALLBACK(toggle_button_toggled_callback),
                             GINT_TO_POINTER(imol));

            GtkStyleContext *context = gtk_widget_get_style_context(b);
            if (ci.imol_is_part_of_cluster(imol)) {
               // Add CSS class
               gtk_style_context_add_class(context, "custom-cfc-button");
            } else {
               // Add CSS class
               gtk_style_context_add_class(context, "custom-cfc-blank-button");
            }

            // gtk_grid_attach(GTK_GRID(grid), b, 1+jj, i, 1, 1);

            gtk_box_append(GTK_BOX(hbox), b);
         }
      }
   } else {
      std::cout << "ERROR:: no grid widget" << std::endl;
   }
}


extern "C" G_MODULE_EXPORT
void
on_cfc_dialog_close(GtkDialog       *dialog,
                    gpointer         user_data) {

   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_cfc_dialog_response(GtkDialog       *dialog,
                       gint             response_id,
                       gpointer         user_data) {

   std::cout << "on_cfc_dialog_response with response_id " << response_id << std::endl;

}

