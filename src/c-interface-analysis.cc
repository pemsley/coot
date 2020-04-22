/* src/c-interface.h
 * 
 * Copyright 2011, 2012 by The University of Oxford
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <iomanip> 
#include "compat/coot-sysdep.h"

// We don't use gtkgraph these days
// #include "libgtkgraph/gtkgraph.h"

#ifdef HAVE_GOOCANVAS
#include "goograph/goograph.hh"
#endif

#include "coot-utils/coot-map-utils.hh" // for variance map

#include "graphics-info.h"
#include "cc-interface.hh"
#include "c-interface.h"
#include "c-interface-generic-objects.h"

#include "coot-utils/coot-hole.hh"

/* ------------------------------------------------------------------------- */
/*                      HOLE                                                 */
/* ------------------------------------------------------------------------- */
/*! \name Coot's Hole implementation */
// if export_dots_file_name string length is zero, then don't try to
// export the surface dots.
void hole(int imol, float start_x, float start_y, float start_z,
	  float end_x, float end_y, float end_z,
	  float colour_map_multiplier, float colour_map_offset,
	  int n_runs, bool show_probe_radius_graph_flag,
	  std::string export_dots_file_name ) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      clipper::Coord_orth p_1(start_x, start_y, start_z);
      clipper::Coord_orth p_2(  end_x,   end_y,   end_z);

      coot::hole hole(mol, p_1, p_2, *g.Geom_p());
      hole.set_colour_shift(colour_map_multiplier, colour_map_offset);
      std::pair<std::vector<std::pair<clipper::Coord_orth, double> >,
		std::vector<coot::hole_surface_point_t> >
	 hole_path_and_surface = hole.generate();

      const std::vector<std::pair<clipper::Coord_orth, double> > &probe_path =
	 hole_path_and_surface.first;

      int obj_path    = new_generic_object_number("Probe path");
      int obj_surface = new_generic_object_number("Probe surface");
   
      for (unsigned int i=0; i<probe_path.size(); i++) {
	 to_generic_object_add_point(obj_path, "red", 3,
				     probe_path[i].first.x(),
				     probe_path[i].first.y(),
				     probe_path[i].first.z());
      }

      for (unsigned int i=0; i<hole_path_and_surface.second.size(); i++) { 
	 to_generic_object_add_point(obj_surface,
				     hole_path_and_surface.second[i].colour.hex().c_str(),
				     1, // pixel
				     hole_path_and_surface.second[i].position.x(),
				     hole_path_and_surface.second[i].position.y(),
				     hole_path_and_surface.second[i].position.z());
      }

      set_display_generic_object(obj_path,    1);
      set_display_generic_object(obj_surface, 1);

      std::string text;
      double path_length = sqrt((p_1-p_2).lengthsq());
      unsigned int n = probe_path.size();
      for (unsigned int i=0; i<n; i++) {
	 double f = path_length * double(i)/double(n);
	 std::string line;
	 line += coot::util::float_to_string_using_dec_pl(f, 4);
	 line += "      ";
	 line += coot::util::float_to_string_using_dec_pl(probe_path[i].second, 4);
	 line += "\n";
	 text += line;
      }

      std::string file_name("probe-radius.tab");
      std::ofstream radius_stream(file_name.c_str());
      if (! radius_stream) {
	 std::cout << "WARNING:: failure to open " << file_name<< std::endl;
      } else { 
	 for (unsigned int i=0; i<n; i++) {
	    double f = path_length * double(i)/double(n);
	    radius_stream << f << "    "  << probe_path[i].second
			  << "\n";
	 }
      }

      if (! export_dots_file_name.empty()) {
	 std::ofstream f(export_dots_file_name.c_str());
	 if (!f) {
	    std::cout << "WARNING:: Failed to open " << export_dots_file_name << std::endl;
	 } else {
	    unsigned int n_path = hole_path_and_surface.second.size();
	    for (unsigned int i=0; i<n_path; i++) {
	       f << "    "
		 << std::fixed << std::setprecision(5) << std::setw(10)
		 << hole_path_and_surface.second[i].position.x() << " " 
		 << std::fixed << std::setprecision(5) << std::setw(10)
		 << hole_path_and_surface.second[i].position.y() << " " 
		 << std::fixed << std::setprecision(5) << std::setw(10)
		 << hole_path_and_surface.second[i].position.z() << "    "
		 << std::fixed << std::setprecision(5) << std::setw(10)
		 << hole_path_and_surface.second[i].colour.red   << " " 
		 << std::fixed << std::setprecision(5) << std::setw(10)
		 << hole_path_and_surface.second[i].colour.green << " " 
		 << std::fixed << std::setprecision(5) << std::setw(10)
		 << hole_path_and_surface.second[i].colour.blue  << "  "
		 << hole_path_and_surface.second[i].colour.hex()
		 << "\n";
	    }
	 } 
      } 

      std::pair<int, int> geom(160, 400);
      simple_text_dialog("Probe radius data", text, geom.first, geom.second);

      if (show_probe_radius_graph_flag)
	 show_hole_probe_radius_graph(probe_path, path_length);


      bool export_map = true;
      if (export_map) {
	 int imol_map = imol_refinement_map();
	 if (is_valid_map_molecule(imol_map)) { 
	    const clipper::Xmap<float> &ref_map = g.molecules[imol_map].xmap;
	    hole.carve_a_map(probe_path, ref_map, "hole.map");
	 }
      }

      export_map = false; // for commit while debugging.
      if (export_map) {
	 clipper::NXmap<float> nxmap = hole.carve_a_map(probe_path, "hole_nx.map");
      }
   }
}

void show_hole_probe_radius_graph(const std::vector<std::pair<clipper::Coord_orth, double> > &hole_path, double path_length) {

#ifdef HAVE_GOOCANVAS
   show_hole_probe_radius_graph_goocanvas(hole_path, path_length);
#else
   show_hole_probe_radius_graph_basic(hole_path, path_length);
#endif    
}

void show_hole_probe_radius_graph_basic(const std::vector<std::pair<clipper::Coord_orth, double> > &hole_path, double path_length) {

   
   GtkWidget *d = gtk_dialog_new();
   g_object_set_data(G_OBJECT(d), "probe_radius_graph", d);
   gtk_window_set_title (GTK_WINDOW(d), "Probe Radius Graph");
   gtk_window_set_default_size(GTK_WINDOW(d), 600, 500);
   // GtkWidget *vbox = GTK_DIALOG(d)->vbox;
   GtkWidget *vbox = gtk_dialog_get_content_area(GTK_DIALOG(d));
   GtkWidget *vbox_inner = gtk_box_new(GTK_ORIENTATION_VERTICAL, 2);
   GtkWidget *scrolled_window = gtk_scrolled_window_new (NULL, NULL);
   // this is deprecated
   // gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrolled_window),
   // GTK_WIDGET(vbox_inner));
   // replace with:
   gtk_container_add(GTK_CONTAINER(scrolled_window), vbox_inner);

   gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(scrolled_window), TRUE, TRUE, 2);
   gtk_widget_show(scrolled_window);
   gtk_widget_show(vbox_inner);

/* Replace by goograph GTK-FIXME

   GtkWidget *graph = gtk_graph_new(XY);
   gtk_graph_set_title(GTK_GRAPH(graph), "Probe Radius along Path", NULL);
   gtk_container_add(GTK_CONTAINER(vbox_inner), graph);
   gtk_widget_show(graph);
   gtk_widget_show(d);
   GtkWidget *close_button = gtk_dialog_add_button(GTK_DIALOG(d), "Close", 2);
   gtk_widget_show(close_button);
   g_signal_connect(G_OBJECT(close_button), "clicked",
		    G_CALLBACK(probe_radius_graph_close_callback),
		    (gpointer) d);

   // graph needs to be realized before gtk_graph_trace_new()
   // returns a useful result.
   int trace = gtk_graph_trace_new(GTK_GRAPH(graph));
   unsigned int n = hole_path.size();
   gfloat *x = (gfloat *) malloc (n * sizeof(gfloat));
   gfloat *y = (gfloat *) malloc (n * sizeof(gfloat));
   for (unsigned int i=0; i<n; i++) {
      x[i] = path_length * double(i)/double(n);
      y[i] = hole_path[i].second;
   }
   gtk_graph_trace_set_data(GTK_GRAPH(graph), trace, x, y, n);
   gtk_graph_axis_format(GTK_GRAPH(graph),
			 GTK_GRAPH_AXIS_INDEPENDANT,
			 FLOATING_POINT,
			 1, "Position along Path");
   gtk_graph_axis_format(GTK_GRAPH(graph),
			 GTK_GRAPH_AXIS_DEPENDANT,
			 FLOATING_POINT,
			 2, "Radius (A)");
   gtk_graph_axis_set_limits(GTK_GRAPH(graph),
			     GTK_GRAPH_AXIS_INDEPENDANT,
			     path_length+1.0, 0.0);
   gtk_graph_axis_set_limits(GTK_GRAPH(graph),
			     GTK_GRAPH_AXIS_DEPENDANT,
			     3.5, 0.0);
   gtk_graph_axis_set_tick(GTK_GRAPH(graph), GTK_GRAPH_AXIS_INDEPENDANT, 5.0, 1.0);
   gtk_graph_trace_format_title(GTK_GRAPH(graph), trace, "Hole Radius");
*/
}

void show_hole_probe_radius_graph_goocanvas(const std::vector<std::pair<clipper::Coord_orth, double> > &hole_path, double path_length) {

#ifdef HAVE_GOOCANVAS
   if (graphics_info_t::use_graphics_interface_flag) { 
      coot::goograph* g = new coot::goograph;
      int trace = g->trace_new();
      std::vector<std::pair<double, double> > data;
      unsigned int n = hole_path.size();
      double y_max = 0.1;
      for (unsigned int i=0; i<n; i++) {
	 double x, y;
	 x = path_length * double(i)/double(n);
	 y = hole_path[i].second;
	 std::pair<double, double> p(x,y);
	 data.push_back(p);
	 if (y > y_max)
	    y_max = y;
      }
      g->set_plot_title("Probe Radius Graph");
      g->set_data(trace, data);
      g->set_extents(coot::goograph::Y_AXIS, 0, y_max);
      g->set_ticks(coot::goograph::X_AXIS, 5, 1);
      g->set_ticks(coot::goograph::Y_AXIS, 1, 0.2);
      g->set_axis_label(coot::goograph::X_AXIS, "Position along Path");
      g->set_axis_label(coot::goograph::Y_AXIS, "Radius");
      g->set_trace_type(trace, coot::graph_trace_info_t::PLOT_TYPE_LINE);
      g->show_dialog();
   }
#endif // HAVE_GOOCANVAS   
}



void probe_radius_graph_close_callback( GtkWidget *button,
 					GtkWidget *dialog) {

   gtk_widget_destroy(dialog);
} 


#ifdef USE_GUILE
SCM model_composition_statistics_scm(int imol) {

   SCM r = SCM_EOL;

   /* 
Model composition
   Non-hydrogen atoms
   Protein residues
   Zinc ions (or other ligands) 

RMS deviations
   Bonds
   Angles

Validation
   Number of C-beta deviation outliers
   Atom overlap score [Clash score]
   Good rotamers (%)
   Ramachandran plot
       Favoured
       Outliers
   */

   if (is_valid_model_molecule(imol)) {
      coot::model_composition_stats_t s = graphics_info_t::molecules[imol].get_model_composition_statistics();
   }

   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *model_composition_statistics_py(int imol) {

   PyObject *r = Py_False;

   if (PyBool_Check(r))
     Py_INCREF(r);
   return r;

}
#endif
