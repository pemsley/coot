
#include <filesystem>

#include "utils/coot-utils.hh"
#include "cfc-gui.hh"
#include "graphics-info.h"
#include "c-interface.h"

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
cfc_gui_t::fill_ligands_grid() {

   auto sorter = +[] (const cfc::typed_cluster_t &t1, const cfc::typed_cluster_t &t2) {
      return t1.imols_with_specs.size() > t2.imols_with_specs.size();
   };

   std::string style_css = "button.custom-cfc-button { \
                               background: #576F50;        /* Green */ \
                               color: white;       \
                               padding: 0;         \
                               margin: 0;          \
                               min-width: 3px;   /* adjust as needed */ \
                               border-radius: 3px; \
                             } \
                               \
                            button.custom-cfc-blank-button { \
                               padding: 0;         \
                               color: #202020;     \
                               margin: 0;          \
                               min-width: 3px;     \
                               border-radius: 1px; \
                            } \
                            ";

GtkWidget *grid = widget_from_builder("cfc-ligands-grid");
   if (grid) {
      std::cout << "-------- grid " << grid
		<< " n-clusters: " << cluster_infos.size() << std::endl;

      std::set<int> imols_in_cluster;
      for (unsigned int i=0; i<cluster_infos.size(); i++) {
	 const auto &ci = cluster_infos[i];
	 for (unsigned int jj=0; jj<ci.imols_with_specs.size(); jj++)
	    imols_in_cluster.insert(ci.imols_with_specs[jj].first);
      }

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
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tb), TRUE);
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
	    // don't recentre if we are merely turnning off the representation
	    if (state == 1)
	       graphics_info_t::set_rotation_centre(pt);
	    for (unsigned int ii=0; ii<tc->imols_with_specs.size(); ii++) {
	       int imol = tc->imols_with_specs[ii].first;
	       set_mol_displayed(imol, state);
	    }
	 };

	 // the molecule toggle buttons (not the feature toggle buttons)
	 auto toggle_button_toggled_callback = +[] (GtkToggleButton *b, gpointer data) {
	    int imol = GPOINTER_TO_INT(data);
	    if (gtk_toggle_button_get_active(b))
	       set_mol_displayed(imol, 1);
	    else
	       set_mol_displayed(imol, 0);
	 };

	 cfc::typed_cluster_t *tc = new cfc::typed_cluster_t(ci);
	 g_signal_connect(G_OBJECT(tb), "toggled", G_CALLBACK(callback), tc);

	 std::set<int>::const_iterator it;
	 for (it=imols_in_cluster.begin(); it!=imols_in_cluster.end(); ++it) {
	    int imol = *it;
	    GtkWidget *b = gtk_toggle_button_new_with_label(".");

	    g_signal_connect(G_OBJECT(b), "toggled", G_CALLBACK(toggle_button_toggled_callback),
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

