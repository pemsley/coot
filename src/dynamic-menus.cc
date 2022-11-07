
#include "dynamic-menus.hh"

#include "graphics-info.h"
#include "c-interface.h"
#include "gtk-manual.h"

#include "widget-from-builder.hh"


void
create_initial_validation_graph_submenu_generic(GtkWidget *widget, // window
						const std::string &menu_name,
						const std::string &sub_menu_name) {

#if (GTK_MAJOR_VERSION >= 4)
   std::cout << "in create_initial_validation_graph_submenu_generic() FIXME" << std::endl;
#else

   // GtkWidget *b_factor_menu_item = lookup_widget(widget, menu_name.c_str());
   GtkWidget *b_factor_menu_item = widget_from_builder(menu_name);
   GtkWidget *b_factor_sub_menu = gtk_menu_new();

   // gtk_widget_ref(b_factor_sub_menu);
   g_object_set_data_full(G_OBJECT(widget), sub_menu_name.c_str(), b_factor_sub_menu, NULL);
   gtk_menu_item_set_submenu(GTK_MENU_ITEM(b_factor_menu_item), b_factor_sub_menu);

   if (false)
      std::cout << "create_initial_validation_graph_submenu_generic: for menu " << menu_name
                << " setting data for " << sub_menu_name << std::endl;

   g_object_set_data(G_OBJECT(b_factor_menu_item), sub_menu_name.c_str(), b_factor_sub_menu);
#endif

}

// #include "widget-headers.hh"

void create_dynamic_menus(GtkWidget *window1) {

   // 20211002-PE is window1 actually needed now?

	 // We need to somehow connect the submenu to the menu's (which are
	 // accessible via window1)
	 //
         // No longer put map colours in the menus
         // create_initial_map_color_submenu(window1);

	 create_initial_ramachandran_mol_submenu(window1);
	 create_initial_sequence_view_mol_submenu(window1);

	 // old style non-generic functions
	 //      create_initial_validation_graph_b_factor_submenu(window1);
	 //      create_initial_validation_graph_geometry_submenu(window1);
	 //      create_initial_validation_graph_omega_submenu(window1);

	 // OK, these things work thusly:
	 //
	 // probe_clashes1 is the name of the menu_item set/created in
	 // by glade and is in mapview.glade.
	 //
	 // probe_submenu is something I make up. It must be the same
	 // here and in c-interface-validate.cc's
	 // add_on_validation_graph_mol_options()
	 //
	 // attach a function to the menu item activate function
	 // created by glade in callbacks.c
	 // (e.g. on_probe_clashes1_activate).  The name that is used
	 // there to look up the menu is as above (e.g. probe_clashes1).
	 //
	 // The type defined there is that checked in
	 // c-interface-validate.cc's 
	 // add_on_validation_graph_mol_options()


	 create_initial_validation_graph_submenu_generic(window1, "peptide_omega_analysis1", "omega_submenu");
	 create_initial_validation_graph_submenu_generic(window1, "geometry_analysis1",      "geometry_submenu");
         create_initial_validation_graph_submenu_generic(window1, "temp_fact_analysis1",     "temp_factor_submenu");
	 create_initial_validation_graph_submenu_generic(window1, "rotamer_analysis1",       "rotamer_submenu");
	 create_initial_validation_graph_submenu_generic(window1, "density_fit_analysis1",   "density_fit_submenu");
	 create_initial_validation_graph_submenu_generic(window1, "ncs_differences1",        "ncs_diffs_submenu");
	 create_initial_validation_graph_submenu_generic(window1, "gln_and_asn_b_factor_outliers1",
							 "gln_and_asn_b_factor_outliers_submenu");
	 create_initial_validation_graph_submenu_generic(window1, "temp_fact_variance_analysis1",
							 "temp_factor_variance_submenu");
         create_initial_validation_graph_submenu_generic(window1, "pukka_puckers_1", "pucker_submenu");
}


void validation_graph_b_factor_mol_selector_activate (GMenuItem     *menuitem,
						      gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
      graphics_info_t g;
      g.b_factor_graphs(imol);

}

////B B GRAPH
void validation_graph_calc_b_factor_mol_selector_activate (GMenuItem     *menuitem,
                                                           gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);

      graphics_info_t g;
      g.calc_b_factor_graphs(imol);

}
////E B GRAPH

void validation_graph_geometry_mol_selector_activate (GMenuItem     *menuitem,
						      gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);

      graphics_info_t g;
      g.geometric_distortion(imol);
}

void validation_graph_omega_mol_selector_activate (GMenuItem     *menuitem,
						   gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);

   graphics_info_t g;
   g.omega_graphs(imol);

}

#include "cc-interface-scripting.hh"

void pukka_puckers_mol_selector_activate (GMenuItem     *menuitem,
                                          gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
   // 20211211-PE bleugh
   std::string imol_str = std::to_string(imol);
   std::string ss = "coot_gui.pukka_puckers_qm(" + imol_str + ")";
   safe_python_command("import coot_gui");
   safe_python_command(ss);

}

void validation_graph_rotamer_mol_selector_activate (GMenuItem     *menuitem,
						     gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
   graphics_info_t g;
   g.rotamer_graphs(imol);

}

void validation_graph_density_fit_mol_selector_activate (GMenuItem     *menuitem,
							 gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
   graphics_info_t g;
   g.density_fit_graphs(imol);
}

void probe_mol_selector_activate (GMenuItem     *menuitem,
 				  gpointer         user_data) {

    int imol = GPOINTER_TO_INT(user_data);
    // 20211002-PE goodbye probe stuff
}

void gln_and_asn_b_factor_outlier_mol_selector_activate (GMenuItem     *menuitem,
							 gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
   gln_asn_b_factor_outliers(imol);
}


GtkWidget *
add_validation_mol_menu_item(int imol,
			     const std::string &name,
			     GtkWidget *menu,
 			     GCallback callback) {

#if (GTK_MAJOR_VERSION >= 4)
   return nullptr;
   // 20220602-PE FIXME menus
   std::cout << "in add_validation_mol_menu_item() FIXME" << std::endl;
#else
    GtkWidget *menu_item = gtk_menu_item_new_with_label(name.c_str());
    gtk_container_add(GTK_CONTAINER(menu), menu_item);
    g_signal_connect(G_OBJECT(menu_item), "activate", callback, GINT_TO_POINTER(imol));
    gtk_widget_show(menu_item);
    return menu_item;
#endif
}

// #include "c-interface-gtk-widgets.h" // for validation_graph_ncs_diffs_mol_selector_activate

// -----------------------------------------------------
// The geometry graphs have a home on the range:
// -----------------------------------------------------


void validation_graph_ncs_diffs_mol_selector_activate (GMenuItem     *menuitem,
						       gpointer         user_data);

// the menu here is the one set in the glade file - and extracted by name using widget_from_builder()
// e.g. geometry_analysis1
//
void add_on_validation_graph_mol_options(GtkWidget *menu, const char *type_in) {

#if (GTK_MAJOR_VERSION >= 4)
   // 20220602-PE FIXME validation graphs - canvas
   std::cout << "in add_on_validation_graph_mol_options() FIXME big job" << std::endl;
#else

   graphics_info_t g;
   std::string validation_type(type_in);
   std::string sub_menu_name;
   GCallback callback = 0; // depends on type
   bool found_validation_type = false;

   if (validation_type == "ramachandran") {
      callback = G_CALLBACK(rama_plot_mol_selector_activate);
      found_validation_type = true;
      sub_menu_name = "rama_plot_submenu";
   }
   if (validation_type == "sequence_view") {
      callback = G_CALLBACK(sequence_view_mol_selector_activate);
      found_validation_type = true;
      sub_menu_name = "sequence_view_submenu";
   }
   if (validation_type == "b factor") {
      callback = G_CALLBACK(validation_graph_b_factor_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "temp_factor_variance_submenu";
   }
   if (validation_type == "calc b factor") {
      callback = G_CALLBACK(validation_graph_calc_b_factor_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "temp_factor_submenu";
   }
   if (validation_type == "geometry") {
      callback = G_CALLBACK(validation_graph_geometry_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "geometry_submenu";
   }
   if (validation_type == "omega") {
      callback = G_CALLBACK(validation_graph_omega_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "omega_submenu";
   }
   if (validation_type == "puckers") {
      callback = G_CALLBACK(pukka_puckers_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "pucker_submenu";
   }
   if (validation_type == "rotamer") {
      callback = G_CALLBACK(validation_graph_rotamer_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "rotamer_submenu";
   }
   if (validation_type == "density-fit") {
      callback = G_CALLBACK(validation_graph_density_fit_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "density_fit_submenu";
   }
#if 0 // 20211002-PE we don't need this any more.
   if (validation_type == "probe") {
      callback = G_CALLBACK(probe_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "probe_submenu";
   }
#endif
   if (validation_type == "gln_and_asn_b_factor_outliers") {
      callback = G_CALLBACK(gln_and_asn_b_factor_outlier_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "gln_and_asn_b_factor_outliers_submenu";
   }
   if (validation_type == "ncs-diffs") {
      callback = G_CALLBACK(validation_graph_ncs_diffs_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "ncs_diffs_submenu";
   }

   auto my_delete_validaton_graph_mol_option = [] (GtkWidget *w, void *data) {
                                                  gtk_container_remove(GTK_CONTAINER(data), w);
                                               };

   // 20211002-PE this is badly designed. Don't use this as a template.
   // I have added a hack to get it it to work with least amount of rewriting.

   GtkWidget *sub_menu = GTK_WIDGET(g_object_get_data(G_OBJECT(menu), sub_menu_name.c_str()));
   if (sub_menu) {

      gtk_container_foreach(GTK_CONTAINER(sub_menu),
			    my_delete_validaton_graph_mol_option,
			    (gpointer) sub_menu);

      for(int i=0; i<g.n_molecules(); i++) {
         if (g.molecules[i].has_model()) {
            std::string name = graphics_info_t::molecules[i].dotted_chopped_name();
                  if (false)
                     std::cout << "debug:: in add_on_validation_graph_mol_options sub_menu_name:"
                              << sub_menu_name << " " << sub_menu << std::endl;
            add_validation_mol_menu_item(i, name, sub_menu, callback);
         }
      }

      
   } else {
      std::cout << "ERROR:: in add_on_validation_graph_mol_options() sub menu not found: "
                << sub_menu_name << " for menu " << menu << std::endl;
   }
#endif

}
