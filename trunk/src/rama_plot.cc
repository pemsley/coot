/* src/main.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2005 by Bernhard Lohkamp
 * Copyright 2009 by The University of Oxford
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#ifdef _MSC_VER
#define snprintf _snprintf
#endif

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
 
#ifndef RAMA_PLOT
#define RAMA_PLOT

#include <iostream>
#include <algorithm>

#include <gdk/gdkkeysyms.h> // for keyboarding.

#include "coot-utils.hh" // int to string

#include "rama_plot.hh" // has gtk/gtk.h which interface.h needs

#ifdef HAVE_GNOME_CANVAS
  typedef GnomeCanvas GtkCanvas;
  typedef GnomeCanvasItem GtkCanvasItem;
  typedef GnomeCanvasPoints GtkCanvasPoints;
  #define GTK_CANVAS GNOME_CANVAS
  #define GTK_CANVAS_TYPE_CANVAS_RECT GNOME_TYPE_CANVAS_RECT
  #define GTK_CANVAS_TYPE_CANVAS_LINE GNOME_TYPE_CANVAS_LINE
  #define GTK_CANVAS_TYPE_CANVAS_TEXT GNOME_TYPE_CANVAS_TEXT
  #define gtk_canvas_init gnome_canvas_init
  #define gtk_canvas_new  gnome_canvas_new
  #define gtk_canvas_root gnome_canvas_root
  #define gtk_canvas_item_new gnome_canvas_item_new
  #define gtk_canvas_points_new gnome_canvas_points_new
  #define gtk_canvas_points_free gnome_canvas_points_free
  #define gtk_canvas_item_w2i gnome_canvas_item_w2i
  #define gtk_canvas_item_grab gnome_canvas_item_grab
  #define gtk_canvas_item_lower_to_bottom gnome_canvas_item_lower_to_bottom
  #define gtk_canvas_item_lower gnome_canvas_item_lower
  #define gtk_canvas_set_scroll_region gnome_canvas_set_scroll_region
  #define gtk_canvas_item_raise_to_top gnome_canvas_item_raise_to_top
  #define gtk_canvas_item_raise gnome_canvas_item_raise
  #define gtk_canvas_item_move gnome_canvas_item_move
  #define gtk_canvas_item_ungrab gnome_canvas_item_ungrab
  #define gtk_canvas_rect_get_type gnome_canvas_rect_get_type
  #define gtk_canvas_set_pixels_per_unit gnome_canvas_set_pixels_per_unit
  #define gtk_canvas_window_to_world gnome_canvas_window_to_world
#endif

#include "interface.h"
#ifndef HAVE_SUPPORT_H
#define HAVE_SUPPORT_H
#include "support.h"
#endif /* HAVE_SUPPORT_H */

#include "rama_mousey.hh"
#include "clipper/core/coords.h"

#include "c-interface.h" // for the mapview callback in button_press()


using namespace std;


void
coot::rama_plot::init(int imol_in, const std::string &mol_name_in, float level_prefered, float level_allowed, float block_size, short int is_kleywegt_plot_flag_in) {

   imol = imol_in; 
   phipsi_edit_flag = 0;
   backbone_edit_flag = 0;
   init_internal(mol_name_in, level_prefered, level_allowed, block_size, 0,
		 is_kleywegt_plot_flag_in);
   // is_kleywegt_plot_flag = is_kleywegt_plot_flag_in;

}

// We could pass to this init the level_prefered and level_allowed
// here, if we wanted the phi/psi edit and backbone edit to have
// "non-standard" contour levels.
void
coot::rama_plot::init(const std::string &type) { 
 
   if (type == "phi/psi-edit") { 
      phipsi_edit_flag = 1;
      backbone_edit_flag = 0;
      imol = -9999; // magic number used in OK button callback.
      init_internal("Ramachandran Plot", 0.02, 0.002, 10);
      hide_stats_frame();
   }
   if (type == "backbone-edit") { 
      phipsi_edit_flag = 0;
      backbone_edit_flag = 1;
      imol = -9999; // magic number used in OK button callback.
      short int hide_buttons = 1;
      init_internal("Ramachandran Plot", 0.02, 0.002, 10, hide_buttons); 
      hide_stats_frame();
   }
   big_box_item = 0;
}


//  The mapview entry point
//
// hide_buttons is optional arg
void
coot::rama_plot::init_internal(const std::string &mol_name,
			       float level_prefered, float level_allowed,
			       float step_in, 
			       short int hide_buttons,
			       short int is_kleywegt_plot_flag_local) {

   fixed_font_str = "fixed";
#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
   fixed_font_str = "monospace";
#endif

   GtkWidget *app1 = create_dynarama_window ();
   if (hide_buttons == 1) {
      GtkWidget *w;
      w = lookup_widget(app1, "dynarama_ok_button");
      gtk_widget_destroy(w);
      w = lookup_widget(app1, "dynarama_cancel_button");
      gtk_widget_destroy(w);
   }
   dialog = app1;
#ifdef HAVE_GTK_CANVAS
   gtk_canvas_init(); 
#endif

   // set the title of of widget
   GtkWidget *label = lookup_widget(app1, "dynarama_label");
   if (label) 
      gtk_label_set_text(GTK_LABEL(label), mol_name.c_str());
   

   allow_seqnum_offset_flag = 0;

   canvas = GTK_CANVAS(gtk_canvas_new());
   gtk_widget_set_usize(GTK_WIDGET(canvas), 400, 400);

   int ysize = 500;
   if (! is_kleywegt_plot_flag_local) // extra space needed
      ysize = 535;

   gtk_widget_set_usize(app1, 400, ysize);

   gtk_widget_ref(GTK_WIDGET(canvas));
   gtk_object_set_data_full(GTK_OBJECT(app1), "canvas", canvas,
			    (GtkDestroyNotify) gtk_widget_unref);

   GtkWidget *scrolled_window = lookup_widget(GTK_WIDGET(app1),
					      "dynarama_scrolledwindow");
   
   gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrolled_window),
					 GTK_WIDGET(canvas));
   
   gtk_object_set_user_data(GTK_OBJECT(canvas), (char *) this); 
   
   gtk_widget_set_events(GTK_WIDGET(canvas),
			 GDK_EXPOSURE_MASK      |
			 GDK_BUTTON_PRESS_MASK  |
			 GDK_BUTTON_RELEASE_MASK|
			 GDK_POINTER_MOTION_MASK|
			 GDK_KEY_RELEASE_MASK   |
			 GDK_POINTER_MOTION_HINT_MASK);

   if (dialog_position_x > -1)
      gtk_widget_set_uposition(app1, dialog_position_x, dialog_position_y);

   gtk_widget_show (GTK_WIDGET(canvas));
   
   gtk_widget_show (app1);

   // Normally we have a plot from a molecule (and we communicate back
   // to graphics_info_t that we now have one), but occassionally we
   // want to edit the phipsi angle.
   // 
   if (! phipsi_edit_flag && ! backbone_edit_flag)
      // a c-interface function
      set_dynarama_is_displayed(GTK_WIDGET(canvas), imol); 

   setup_internal(level_prefered, level_allowed); 
   step = step_in;
   kleywegt_plot_uses_chain_ids = 0;

   green_box = coot::util::phi_psi_t(-999, -999, "", "", 0, "", ""); // set unsensible phi phi initially.
   setup_canvas(); 
}

// draw a big square after everything else for residue i:
// 
void
coot::rama_plot::big_square(int i,
			    const std::string &ins_code,
			    const std::string &chain_id) {


   for (unsigned int ich=0; ich<phi_psi_sets.size(); ich++) {
      if (phi_psi_sets[ich].chain_id == chain_id) {

	 if (phi_psi_sets[ich].phi_psi.size() > 0) {  // not a water/ligand chain
	 
	    int ires = i-ifirst_res[ich];
	    if (ires >= 0 && ires < int(phi_psi_sets[ich].phi_psi.size())) { 
	       if ((phi_psi_sets[ich].phi_psi[ires].residue_number == i) &&
		   (phi_psi_sets[ich].phi_psi[ires].ins_code == ins_code)) {
		  // normal case (hopefully)
		  draw_phi_psi_point_internal(ires, &phi_psi_sets[ich].phi_psi, 0, 4); // not white, size 4
	       } else {
		  for (unsigned int j=0; j< phi_psi_sets[ich].phi_psi.size(); j++) {
		     if ((phi_psi_sets[ich].phi_psi[j].residue_number == i) &&
			 (phi_psi_sets[ich].phi_psi[j].ins_code == ins_code)) { 
			ires = j;
			draw_phi_psi_point_internal(ires, &phi_psi_sets[ich].phi_psi, 0, 4); // not white, size 4
			break; 
		     }
		  }
	       }
	    } else {
	       std::cout << "trapped out of range phi/psi access: " << i << " " << chain_id
			 << " ich: " << ich << " ires: " << ires
			 << " phi_psi_sets[ich].phi_psi.size(): "
			 << phi_psi_sets[ich].phi_psi.size() << std::endl;
	    }
	 }
      }
   }
}


void
coot::rama_plot::clear_canvas_items() { 

   GtkCanvasItem *item;
   for (unsigned int i=0; i<canvas_item_vec.size(); i++) { 
      item = canvas_item_vec[i];
      if (canvas_item_vec[i] == NULL) {
	 std::cout << "oops - null canvas item" << std::endl;
      }
      gtk_object_destroy(GTK_OBJECT(item));
   }
   canvas_item_vec.resize(0);
}


void
coot::rama_plot::destroy_yourself() {

   if (dialog)
      gtk_widget_destroy(dialog);
   else
      std::cout << "ERROR:: could not get dialog from canvas in rama_plot::destroy_yourself"
		<< std::endl;
} 

void 
coot::rama_plot::clear_last_canvas_item() { 
   
   GtkCanvasItem *item;
   int n = canvas_item_vec.size() -1;

   item = canvas_item_vec[n];
   if (item == NULL) {
      std::cout << "oops - null canvas item" << std::endl;
   }
   gtk_object_destroy(GTK_OBJECT(item));
   
   canvas_item_vec.resize(n);
}

void 
coot::rama_plot::clear_last_canvas_items(int np) { 
   
   GtkCanvasItem *item;
   int n = canvas_item_vec.size() - np;

   for (int i=0; i<np; i++) { 
      item = canvas_item_vec[n+i];
      if (item == NULL) {
 	 std::cout << "ERROR:: - null canvas item in clear_last_canvas_items() " << i
 		   << " of " << n << std::endl;
      } else { 
// 	 std::cout << "clearing canvas item " <<  i << " of " << n << " "
// 		   << item << std::endl;
	 gtk_object_destroy(GTK_OBJECT(item));
      }
   }
   canvas_item_vec.resize(n);
}


void
coot::rama_plot::setup_internal(float level_prefered, float level_allowed) {

   zoom = 0.8;
   have_sticky_labels = 0; 

   n_diffs = 50; // default value.
   drawing_differences = 0; 

   tooltip_item = NULL;
   tooltip_item_text = NULL;

   rama.init(clipper::Ramachandran::All);
   displayed_rama_type = clipper::Ramachandran::All;

   // cliper defaults: 
   rama_threshold_preferred = 0.01; 
   rama_threshold_allowed = 0.0005; 

   // values that make the plot look similar to a procheck plot:
   rama_threshold_preferred = 0.06; // 0.05 
   rama_threshold_allowed = 0.0012; // 0.002

   // Procheck is old and too liberal
   rama_threshold_preferred = 0.1;
   rama_threshold_allowed = 0.005;

   // cliper defaults: 
   rama_threshold_preferred = 0.01; 
   rama_threshold_allowed = 0.0005;

   // nipped clipper values:
   rama_threshold_preferred = 0.02; 
   rama_threshold_allowed = 0.0012;
   
   // Lovell et al. 2003, 50, 437 Protein Structure, Function and Genetics values:
   rama_threshold_preferred = 0.02; 
   rama_threshold_allowed = 0.002;
   
   //clipper defaults: 0.01 0.0005

   rama.set_thresholds(level_prefered, level_allowed);
   //
   r_gly.init(clipper::Ramachandran::Gly);
   r_gly.set_thresholds(level_prefered, level_allowed);
   //
   r_pro.init(clipper::Ramachandran::Pro);
   r_pro.set_thresholds(level_prefered, level_allowed);
   // 
   r_non_gly_pro.init(clipper::Ramachandran::NonGlyPro);
   r_non_gly_pro.set_thresholds(level_prefered, level_allowed);
}

void
coot::rama_plot::set_n_diffs(int nd) {
   n_diffs = nd;
}

void 
coot::rama_plot::display_background() {


   // Surely we don't need item and itm? Fix when you're bored.
   
   GtkCanvasItem *item;

   // delete the old background, if there is one
   //
   GtkCanvasItem *itm;
   
   // std::cout << "deleting background" << std::endl;
   for (unsigned int i=0; i<canvas_item_vec.size(); i++) { 
      itm = canvas_item_vec[i];
      if (canvas_item_vec[i] == NULL) {
	 std::cout << "oops - null canvas item" << std::endl;
      }
      // debugging.  FIXME
      //
//       int n_match = 0;
//       for (int j=0; j<canvas_item_vec.size(); j++) {
// 	 if (canvas_item_vec[j] == itm) {
// 	    n_match++;
// 	 }
//       }
//       if (n_match != 1) {
// 	 std::cout << "Caught an error, n_match is " << n_match << std::endl;
//       } 

      gtk_object_destroy(GTK_OBJECT(itm));
   }
   canvas_item_vec.resize(0);

   // std::cout << "new underlay background" << std::endl;
   //
   basic_white_underlay(); 

   float x;
   float y;
   short int doit;
   std::string colour; 

   for (float i= -180.0; i<180.0; i += step) { 
      for (float j= -180.0; j<180.0; j += step) {

	 x =  clipper::Util::d2rad(i+((float) step)/2.0); 
	 y =  clipper::Util::d2rad(-(j+((float) step)/2.0));
	 doit = 0; 

	 if ( rama.allowed(x,y) ) {
	    colour = "grey50";
	    colour= "yellow"; // Procheck colour
             	              // colour = "PaleGoldenrod";
	    colour = "khaki";

	    cell_border(int(i), int(j), int(step));
	    item = gtk_canvas_item_new(gtk_canvas_root(canvas),
					 GTK_CANVAS_TYPE_CANVAS_RECT,
					 "x1", i+0.0,
					 "y1", j+0.0,
					 "x2", i+step+0.0,
					 "y2", j+step+0.0,
					 "fill_color", colour.c_str(),
					 NULL);
	    canvas_item_vec.push_back(item); 

	 }
      }
   }

   // We have 2 loops so that the favoured regions are drawn last and
   // that their borders are not stamped on by the allowed regions.
   // 
   for (float i= -180; i<180; i += step) { 
      for (float j= -180; j<180; j += step) {

	 x = clipper::Util::d2rad(i+((float) step)/2.0);
	 y = clipper::Util::d2rad(-(j+((float) step)/2.0)); 
	 doit = 0; 

	 if ( rama.favored(x,y) ) {
	    colour = "grey";
	    colour = "red";
	    colour = "pink";
	    colour = "HotPink";

	    cell_border(int(i), int(j), int(step));
	    item = gtk_canvas_item_new(gtk_canvas_root(canvas),
					 GTK_CANVAS_TYPE_CANVAS_RECT,
					 "x1", i+0.0,
					 "y1", j+0.0,
					 "x2", i+step+0.0,
					 "y2", j+step+0.0,
					 "fill_color", colour.c_str(),
					 NULL);
	    canvas_item_vec.push_back(item); 

	 }
      }
   }

   black_border();

}


// not const because we modify canvas_item_vec
// 
void
coot::rama_plot::basic_white_underlay() {

   GtkCanvasItem *item; 

   item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				GTK_CANVAS_TYPE_CANVAS_RECT,
				"x1", -180.0,
				"y1", -180.0,
				"x2", 180.0,
				"y2", 180.0,
				"fill_color", "grey80",
				"outline_color", "black",
				NULL);
   
   canvas_item_vec.push_back(item); 

   

} 

void
coot::rama_plot::setup_canvas() {

  // gtk_widget_show(GTK_WIDGET(canvas)); 

   gtk_canvas_set_pixels_per_unit(canvas,zoom);
   // gtk_canvas_set_scroll_region moves about the canvas in the
   // widget, it doesn't seem to change the amount that it is
   // scrollable.
   gtk_canvas_set_scroll_region(canvas,-220.0,-220.0,210.0,210.0);
   gtk_signal_connect (GTK_OBJECT(canvas), "button_press_event",
		       GTK_SIGNAL_FUNC(rama_button_press), NULL);
   
   gtk_signal_connect (GTK_OBJECT(canvas), "motion_notify_event",
		       GTK_SIGNAL_FUNC(rama_motion_notify), NULL);

   // also need expose_event, configure_event, realise

//    gtk_signal_connect(GTK_OBJECT(canvas), "key_press_event",
// 		      GTK_SIGNAL_FUNC(key_press_event), NULL);
   gtk_signal_connect(GTK_OBJECT(canvas), "key_release_event",
		      GTK_SIGNAL_FUNC(rama_key_release_event), NULL);

   /* set focus to canvas - we need this to get key presses. */
  GTK_WIDGET_SET_FLAGS(canvas, GTK_CAN_FOCUS);
  gtk_widget_grab_focus(GTK_WIDGET(canvas));


}

gint
coot::rama_plot::key_release_event(GtkWidget *widget, GdkEventKey *event) {

   switch (event->keyval) {
      
   case GDK_plus:
   case GDK_equal:  // unshifted plus, usually.

      zoom_in();
      break;

   case GDK_minus:
      zoom_out();
      break; 
   }

   /* prevent the default handler from being run */
   gtk_signal_emit_stop_by_name(GTK_OBJECT(canvas),"key_release_event");

   return 0; 
}


void
coot::rama_plot::black_border() {

   
   GtkCanvasItem *item;

   GtkCanvasPoints *points  = gtk_canvas_points_new(5);
   
   points->coords[0] = -180.0; 
   points->coords[1] = -180.0; 

   points->coords[2] = -180.0; 
   points->coords[3] =  180.0; 

   points->coords[4] = 180.0; 
   points->coords[5] = 180.0; 

   points->coords[6] =  180.0; 
   points->coords[7] = -180.0; 

   points->coords[8] = -180.0; 
   points->coords[9] = -180.0; 

   item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				GTK_CANVAS_TYPE_CANVAS_LINE,
				"width_pixels", 2, 
				"points", points,
				"fill_color", "black",
				NULL);
   canvas_item_vec.push_back(item); 

   gtk_canvas_points_free(points); 
}


void
coot::rama_plot::cell_border(int i, int j, int step) {

   // put a border round the canvas one pixel shifted right and up
   //
   GtkCanvasItem *item;
   GtkCanvasPoints *points  = gtk_canvas_points_new(5);
   
   points->coords[0] = i+1;
   points->coords[1] = j+1;

   points->coords[2] = i+step+1;
   points->coords[3] = j+1;

   points->coords[4] = i+step+1;
   points->coords[5] = j+step+1;

   points->coords[6] = i+1;
   points->coords[7] = j+step+1;

   points->coords[8] = i+1;
   points->coords[9] = j+1;

   
   item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				GTK_CANVAS_TYPE_CANVAS_LINE,
				"width_pixels", 2, 
				"points", points,
				"fill_color", "grey50",
				NULL);

   canvas_item_vec.push_back(item); 

   gtk_canvas_points_free(points); 

   
} 

// return the region of the point
int
coot::rama_plot::draw_phi_psi_point(int i,
				    const vector<coot::util::phi_psi_t> *phi_psi_vec,
				    short int as_white_flag) {

   return draw_phi_psi_point_internal(i,phi_psi_vec,as_white_flag, 2); 
}

int 
coot::rama_plot::draw_phi_psi_point_internal(int i,
					     const vector<coot::util::phi_psi_t> *phi_psi_vec,
					     short int as_white_flag,
					     int box_size) {

   int region = coot::rama_plot::RAMA_UNKNOWN; // initially unset
   
   std::string outline_color("black");

   if (box_size == 4) {
      draw_green_box((*phi_psi_vec)[i].phi(), (*phi_psi_vec)[i].psi());
   } else {

      if ((*phi_psi_vec)[i].residue_name() == "GLY") {
	 region = draw_phi_psi_as_gly(i,phi_psi_vec);
      } else {
	 std::string colour;
	 double phi = (*phi_psi_vec)[i].phi();
	 double psi = (*phi_psi_vec)[i].psi();

	 if (rama.allowed(clipper::Util::d2rad(phi),
			  clipper::Util::d2rad(psi))) {
	    colour = "blue";
	    region = coot::rama_plot::RAMA_ALLOWED;
	    if (rama.favored(clipper::Util::d2rad(phi),
			     clipper::Util::d2rad(psi))) {
	       region = coot::rama_plot::RAMA_PREFERRED;
	    }
	 } else {
	    colour = "red";
	    region = coot::rama_plot::RAMA_OUTLIER;
	 }

	 if ( as_white_flag == 1 ) {
	    colour = "white";
	 } else { 
	    if ((*phi_psi_vec)[i].residue_name() == "PRO") {
	       outline_color = "grey";
	 
	       if (r_pro.allowed(clipper::Util::d2rad(phi),
				 clipper::Util::d2rad(psi))) {
		  colour = "blue";
		  region = coot::rama_plot::RAMA_ALLOWED;
		  if (r_pro.favored(clipper::Util::d2rad(phi),
				    clipper::Util::d2rad(psi))) {
		     region = coot::rama_plot::RAMA_PREFERRED;
		  }
	       } else {
		  colour = "red3";
		  region = coot::rama_plot::RAMA_OUTLIER;
	       }
	    } else {

	       // conventional residue
	       if (r_non_gly_pro.allowed(clipper::Util::d2rad(phi),
					 clipper::Util::d2rad(psi))) {
		  region = coot::rama_plot::RAMA_ALLOWED;
		  colour = "blue"; 
		  if (r_non_gly_pro.favored(clipper::Util::d2rad(phi),
					    clipper::Util::d2rad(psi))) {
		     region = coot::rama_plot::RAMA_PREFERRED;
		  }
	       } else {
		  colour = "red3";
		  region = coot::rama_plot::RAMA_OUTLIER;
	       }
	    }
	 }
	 GtkCanvasItem *item = gtk_canvas_item_new(gtk_canvas_root(canvas),
						   GTK_CANVAS_TYPE_CANVAS_RECT,
						   "x1", phi-box_size,
						   "y1",-psi-box_size,
						   "x2", phi+box_size,
						   "y2",-psi+box_size,
						   "fill_color", colour.c_str(),
						   "outline_color", outline_color.c_str(),
						   NULL);
	 canvas_item_vec.push_back(item);
      }
   }

   // a GtkCanvasItem is not the right widget:
   // gtk_tooltips_set_tip(tooltips, GTK_WIDGET(item), "This is a
   // residue", NULL);

   // maybe this works - dont know how though.  It seems to be making
   // a tooltip window in the root window.
   // 
   // meta_fixed_tip_show( int(phi), int(psi), "testing tooltip");

   return region;
}

// store the green box position too
void
coot::rama_plot::draw_green_box(double phi, double psi) {


   // std::cout << "============ adding green box at " << phi << " " << psi << std::endl;
   std::string colour = "green";
   std::string outline_color = "black";
   int box_size = 4;
   
   if (big_box_item)
      gtk_object_destroy(GTK_OBJECT(big_box_item));
   GtkCanvasItem *item = gtk_canvas_item_new(gtk_canvas_root(canvas),
					     GTK_CANVAS_TYPE_CANVAS_RECT,
					     "x1", phi-box_size,
					     "y1",-psi-box_size,
					     "x2", phi+box_size,
					     "y2",-psi+box_size,
					     "fill_color", colour.c_str(),
					     "outline_color", outline_color.c_str(),
					     NULL);
   big_box_item = item;
   green_box = coot::util::phi_psi_t(phi, psi, "", "",  0, "", ""); 

}


void
coot::rama_plot::draw_green_box() {

   if (green_box_is_sensible(green_box))
      draw_green_box(green_box.phi(), green_box.psi());
}


// could/should be a member function of coot::phi_psi_t
bool
coot::rama_plot::green_box_is_sensible(coot::util::phi_psi_t gb) const { // have the phi and psi been set to something sensible?

   bool r = 0;
   if (gb.phi()>= -180.0)
      if (gb.phi()<= 180.0)
	 r = 1;
   return r;
}


int // return region
coot::rama_plot::draw_phi_psi_as_gly(int i, const vector<coot::util::phi_psi_t> *phi_psi_vec) {

   GtkCanvasItem *item;

   // clipper::Ramachandran r_tmp;
//    r_gly.init(clipper::Ramachandran::Gly);
//    r_gly.set_thresholds(0.05,0.002); //defaults 0.01 0.0005

   std::string colour;
   int region;

   double phi = (*phi_psi_vec)[i].phi();
   double psi = (*phi_psi_vec)[i].psi();
   
   if (r_gly.allowed(clipper::Util::d2rad(phi), clipper::Util::d2rad(psi))) {
      colour = "blue";
      region = coot::rama_plot::RAMA_ALLOWED;
      if (r_gly.favored(clipper::Util::d2rad(phi), clipper::Util::d2rad(psi))) {
	 region = coot::rama_plot::RAMA_PREFERRED;
      }
   } else {
      colour = "red3";
      region = coot::rama_plot::RAMA_OUTLIER;
   }

   GtkCanvasPoints *points  = gtk_canvas_points_new(4);

   points->coords[0] =  phi;
   points->coords[1] = -psi-3;

   points->coords[2] =  phi-3;
   points->coords[3] = -psi+3;

   points->coords[4] =  phi+3;
   points->coords[5] = -psi+3;
   
   points->coords[6] =  phi;
   points->coords[7] = -psi-3;

   item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				GTK_CANVAS_TYPE_CANVAS_LINE,
				"width_pixels", 2, 
				"points", points,
				"fill_color", colour.c_str(),
				NULL); 
   canvas_item_vec.push_back(item);
   gtk_canvas_points_free(points);

   return region;
}


// Uses a class data member
coot::rama_stats_container_t 
coot::rama_plot::draw_phi_psi_points(int ich) {


   coot::rama_stats_container_t counts;
   short int as_white_flag = 0; 
   int n_sets = phi_psi_sets[ich].phi_psi.size();
//    std::cout << "draw_phi_psi_points " << size
// 	     << " residues for chain number " << ich << std::endl;
   for (int i=0; i<n_sets; i++) {
      int type = draw_phi_psi_point(i, &phi_psi_sets[ich].phi_psi, as_white_flag);
      if (type != coot::rama_plot::RAMA_UNKNOWN) {
	 counts.n_ramas++; 
	 if (type == coot::rama_plot::RAMA_ALLOWED)
	    counts.n_allowed++;
	 if (type == coot::rama_plot::RAMA_PREFERRED)
	    counts.n_preferred++;
      }
   }
   return counts;
}

// fill phi_psi vector
//
// We want another version of this routine, one in which
// we get passed the chain_id.
// 
void
coot::rama_plot::generate_phi_psis(std::vector <phi_psi_set_container> *phi_psi_set_vec,
				   CMMDBManager *mol_in) {

   // clear out whatever data we have in the phi_psi's currently
   // for(int i=0; i<phi_psi_set_vec->size();)
   phi_psi_set_vec->resize(0);
   
   // InitMatType();

   // PCAtom atom1;
   CMMDBManager *mol = mol_in;

   // int selHnd = mol->NewSelection();
   // int nSelAtoms;
   // PPCAtom SelAtom;

   // int resno;
   std::string chain_id;
   // int nres;

   int nmodels = mol->GetNumberOfModels();
   // std::cout << "models: " << nmodels << std::endl;
   for (int imodel=1; imodel<=nmodels; imodel++) {
      int nchains = mol->GetNumberOfChains(imodel);

      ifirst_res.resize(nchains);
      PCChain chn;
      int nres;
      PCResidue res;
      int seqno;
      for (int ichain=0; ichain<nchains; ichain++) {

	 // int nres = mol->GetNumberOfResidues(imodel,ichain); old, delete me.

	 // or PCModel model = mol->GetModel(imodel)
	 // 
	 // or PCChain chn = model->GetChain(ichain);
	 chn = mol->GetChain(imodel, ichain);
	 nres = chn->GetNumberOfResidues();
	 std::string chain_name = chn->GetChainID();
	 std::cout << "    chain: " << chain_name << " "
		   << ichain << " has " << nres
		   << " residues" << std::endl;

	 // make a new chain for phi/psi.  We need phi_psi_set_vec to
	 // be consistent with ichain, so we sometimes add chains that
	 // have no ramachandran values (e.g.) ligands
	 // 
	 phi_psi_set_container a(chain_name);
	 phi_psi_set_vec->push_back(a);
	    
	 if (nres > 2) {
	    
	    // no shadowing of the loop:
	    res = mol->GetResidue(imodel,ichain,1);
	    seqno  = res->GetSeqNum();
	    ifirst_res[ichain] = seqno;
	    
	    // iresidues in the chain go from 0 to nres-1, but for getting
	    // ramachandran angles, we want limts: 1 to nres-2
	    //
	    CResidue *prev, *next; 
	    for (int ires=1; ires<(nres-1); ires++) {
	       res = mol->GetResidue(imodel,ichain,ires);
	       prev = mol->GetResidue(imodel,ichain,ires-1);
	       if (res && prev) { 
		  seqno  = res->GetSeqNum();
		  next = mol->GetResidue(imodel,ichain,ires+1);
		  // std::cout << "adding phi/psi for " << seqno << " "
		  // << chain_name << std::endl;
		  add_phi_psi(& (*phi_psi_set_vec)[ichain].phi_psi,
			      mol, chain_name.c_str(), prev, res, next);
	       }
	    }
	 }
      }
   }

// debugging   
//    for (int ich=0; ich<phi_psi_set_vec->size(); ich++) { 
//       cout << "we got "<< (*phi_psi_set_vec)[ich].phi_psi.size()
// 	   << " ramachandran angles for chain " << (*phi_psi_set_vec)[ich].chain_id
// 	   << endl;
//    }
}




void
coot::rama_plot::generate_phi_psis_debug()
{
   int ich = 0;
   std::string resname("GLY");
   std::string label("debug_label"); 
   for(double phi=-180; phi<180; phi += 50) { 
      for(double psi=-180; psi<180; psi += 50) { 
         phi_psi_sets[ich].phi_psi.push_back(coot::util::phi_psi_t(phi, psi,resname,
								   label, 1, "", "A"));
      }
   }
}

void
coot::rama_plot::map_mouse_pos(double x, double y) {

   // std::cout << "DEBUG:: canvas: " << canvas << std::endl;

   mouse_util_t t = mouse_point_check(x,y); 

   // if we are close to (over) a residue box/triangle, then display a
   // tooltip_like_box, if not, then if we don't have sticky labels,
   // then delete the one that's was being displayed.
   // 
   if (t.set_smallest_diff_flag) {

      if (drawing_differences) { 
	 residue_type_background_as(phi_psi_sets[t.ichain].phi_psi[t.smallest_diff_i].residue_name());
	 tooltip_like_box(t);
      } else { 
	 residue_type_background_as(phi_psi_sets[t.ichain].phi_psi[t.smallest_diff_i].residue_name());
	 tooltip_like_box(t);
      } 
   } else {

      if (! have_sticky_labels) {

	 if (tooltip_item != NULL) { 
	    gtk_object_destroy(GTK_OBJECT(tooltip_item));
	    tooltip_item = NULL;
	 } 
	 if (tooltip_item_text != NULL) { 
	    gtk_object_destroy(GTK_OBJECT(tooltip_item_text));
	    tooltip_item_text = NULL;
	 }
      }
   }
}

void
coot::rama_plot::mouse_motion_notify(GdkEventMotion *event, double x, double y) {

   if (phipsi_edit_flag)
      mouse_motion_notify_editphipsi(event, x, y);
   else 
      map_mouse_pos(x,y);
}

void 
coot::rama_plot::mouse_motion_notify_editphipsi(GdkEventMotion *event, double x, double y) { 

   GdkModifierType state;
   int x_as_int, y_as_int;
   // int x, y;

   if (event->is_hint) {
      gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
   }
    

   if ( state & GDK_BUTTON1_MASK ) {
      double worldx, worldy;
      gtk_canvas_window_to_world(canvas, x,y, &worldx, &worldy);
      // std::cout << "drag! " << worldx << " " << worldy << std::endl;
      
      // so worldx is phi and worldy is -psi.
      // 
      // clamp the vaules to -180 to + 180:
      worldx = worldx >  180.0 ?  180.0 : worldx;
      worldx = worldx < -180.0 ? -180.0 : worldx;
      worldy = worldy >  180.0 ?  180.0 : worldy;
      worldy = worldy < -180.0 ? -180.0 : worldy;

      // set the point to new values:
      phi_psi_sets[0].phi_psi[0] = coot::util::phi_psi_t(worldx, -worldy, "dum",
							 "moving", 0, "", "");
      clear_last_canvas_item();
      draw_phi_psi_points(0);

      set_moving_atoms(worldx, -worldy);
   
   }
} 

gint
coot::rama_plot::button_press (GtkWidget *widget, GdkEventButton *event) {

   if (phipsi_edit_flag) 
      return button_press_editphipsi(widget, event);
   else
      if (backbone_edit_flag) 
	 return button_press_backbone_edit(widget, event);
      else
	 return button_press_conventional(widget, event);
      
}

gint
coot::rama_plot::button_press_backbone_edit (GtkWidget *widget, GdkEventButton *event) {

   return 0;
}

gint
coot::rama_plot::button_press_conventional (GtkWidget *widget, GdkEventButton *event) {

   //
   
   double x,y;
   // int x_as_int, y_as_int;
   GdkModifierType state;

   x = event->x;
   y = event->y;
   state = (GdkModifierType) event->state;

   if (event->button == 1) { 

      mouse_util_t t = mouse_point_check(x,y);
      int ich = t.ichain;
      if (t.set_smallest_diff_flag) {
	 if (t.mouse_over_secondary_set) {
	    cout << "secondary: ";
	    cout << secondary_phi_psi_sets[ich].phi_psi[t.smallest_diff_i].label()
		 << ", phi=" << secondary_phi_psi_sets[ich].phi_psi[t.smallest_diff_i].phi() << " "
		 << " psi=" << secondary_phi_psi_sets[ich].phi_psi[t.smallest_diff_i].psi() << " "
		 << endl; 
	    set_go_to_atom_molecule(molecule_numbers_.second);
	    int ires = t.smallest_diff_i; 
	    set_go_to_atom_chain_residue_atom_name((gchar *)secondary_phi_psi_sets[t.ichain].phi_psi[ires].chain_id.c_str(),
						   secondary_phi_psi_sets[t.ichain].phi_psi[ires].residue_number,
						   "CA");
	 } else {
	    cout << "primary: ";
	    cout << phi_psi_sets[t.ichain].phi_psi[t.smallest_diff_i].label()
		 << ", phi=" << phi_psi_sets[t.ichain].phi_psi[t.smallest_diff_i].phi() << " "
		 <<  " psi=" << phi_psi_sets[t.ichain].phi_psi[t.smallest_diff_i].psi() << " "
		 << endl;
	    // the all important mapview callback
	    if (drawing_differences) 
	       set_go_to_atom_molecule(molecule_numbers_.first);
	    else 
	       set_go_to_atom_molecule(imol);
	    int ires = t.smallest_diff_i; 
	    set_go_to_atom_chain_residue_atom_name((gchar *)phi_psi_sets[t.ichain].phi_psi[ires].chain_id.c_str(),
						   phi_psi_sets[t.ichain].phi_psi[ires].residue_number,
						   "CA");
	 }
      }
   } else {

      if (event->button == 3) {
	 zoom += 0.2;
	 gtk_canvas_set_pixels_per_unit(canvas,zoom);
      }
   }
   return 0; 
}

gint
coot::rama_plot::button_press_editphipsi (GtkWidget *widget, GdkEventButton *event) {

   double x,y;
   // int x_as_int, y_as_int;
   GdkModifierType state;

   x = event->x;
   y = event->y;
   state = (GdkModifierType) event->state;

   if (event->button == 1) { 

      mouse_util_t t = mouse_point_check(x,y);
      if (t.set_smallest_diff_flag) {
	    cout << "primary: ";
	    cout << phi_psi_sets[t.ichain].phi_psi[t.smallest_diff_i].label()
		 << ", phi=" << phi_psi_sets[t.ichain].phi_psi[t.smallest_diff_i].phi() << " "
		 <<  " psi=" << phi_psi_sets[t.ichain].phi_psi[t.smallest_diff_i].psi() << " "
		 << endl;
      }
   }
   return 0; 
}

// provides a object that contains the residue number and the chain
// number (for phi_psi_set indexing) of the mouse position
// 
coot::mouse_util_t
coot::rama_plot::mouse_point_check(double x, double y) const {

   coot::mouse_util_t t;

   double worldx, worldy;
   double diff1x, diff1y;
   double diff2x, diff2y;
   double smallest_diff = 999999; 

//    int n = phi_psi_sets[t.ichain].phi_psi.size(); // FIXME?
//    int ich = t.ichain;

   t.set_smallest_diff_flag = 0;
   t.mouse_over_secondary_set = 0; 
   
   gtk_canvas_window_to_world(canvas, x,y, &worldx, &worldy);

   if (drawing_differences) {

      // need only to check over the top n_diffs differences.  (but we
      // do need to check the primary and the secondary phi_psi sets).
      //
      unsigned int i; 
      if (worldx <= 180 && worldx >= -180) { 
	 if (worldy <= 180 && worldy >= -180) {
	    unsigned int ich, ich_2;
	    
	    for (int j=0; j<int(diff_sq.size()) && j<n_diffs; j++) {
	       
	       i = diff_sq[j].i(); 
	       ich = diff_sq[j].ich();
	       ich_2 = diff_sq[j].ich2();

		     
	       if ((ich < phi_psi_sets.size()) && (ich < secondary_phi_psi_sets.size())) {
		  if (i < phi_psi_sets[ich].phi_psi.size()){ 
	       
		     diff1x = fabs(phi_psi_sets[ich].phi_psi[i].phi() - worldx); 
		     diff1y = fabs(phi_psi_sets[ich].phi_psi[i].psi() + worldy);
		     
		     diff2x = fabs(secondary_phi_psi_sets[ich_2].phi_psi[i].phi() - worldx); 
		     diff2y = fabs(secondary_phi_psi_sets[ich_2].phi_psi[i].psi() + worldy);
		     
		     if (diff1x < 3 && diff1y < 3){
			t.set_smallest_diff_flag = 1; 
			if (diff1x+diff1y< smallest_diff) {
			   t.mouse_over_secondary_set = 0; 
			   smallest_diff = diff1x + diff1y;
			   t.smallest_diff_i = i;
			   t.ichain = ich;
			}
		     }
		     if (diff2x < 3 && diff2y < 3){
			t.set_smallest_diff_flag = 1; 
			if (diff2x+diff2y< smallest_diff) {
			   t.mouse_over_secondary_set = 1; 
			   smallest_diff = diff2x + diff2y;
			   t.smallest_diff_i = i;
			   t.ichain = ich; // correct? CHECKME
			}
		     }
		  }
	       }
	    }
	 }
      }
   } else {

      // single molecule, check over all ramachandran angles in the
      // first list:
      // 
      if (worldx <= 180 && worldx >= -180) { 
	 if (worldy <= 180 && worldy >= -180) {
	    
	    for (unsigned int ich=0; ich<phi_psi_sets.size(); ich++) { 
	       
	       for (unsigned int i=0; i<phi_psi_sets[ich].phi_psi.size(); i++) {
		  
		  diff1x = fabs(phi_psi_sets[ich].phi_psi[i].phi() - worldx); 
		  diff1y = fabs(phi_psi_sets[ich].phi_psi[i].psi() + worldy);
		  
		  if (diff1x < 3 && diff1y < 3) {
		     
		     t.set_smallest_diff_flag = 1; // we have at least one
		     if (diff1x+diff1y< smallest_diff) {
			smallest_diff = diff1x + diff1y;
			t.smallest_diff_i = i;
			t.ichain = ich;
		     }
		  }
	       } 
	    }
	 }
      }
   }
   return t; 
} 

// redraw everything with the given background
// 
void
coot::rama_plot::residue_type_background_as(std::string res) {

   if (res == "GLY" &&
       (displayed_rama_type != clipper::Ramachandran::Gly)) {
      all_plot(clipper::Ramachandran::Gly);
   }
   if (res == "PRO" &&
       (displayed_rama_type != clipper::Ramachandran::Pro)) {
      all_plot(clipper::Ramachandran::Pro);
   }
   if ( res != "GLY"  ) {
      if (res != "PRO") {
	 if (displayed_rama_type != clipper::Ramachandran::NonGlyPro) {
	    all_plot(clipper::Ramachandran::NonGlyPro);
	 }
      }
   }
}

// and bits means zero lines and axes labels and ticks.
// 
void
coot::rama_plot::all_plot_background_and_bits(clipper::Ramachandran::TYPE type) {

   rama.init(type);
   rama.set_thresholds(rama_threshold_preferred,rama_threshold_allowed);
   display_background();
   draw_axes();
   draw_zero_lines();

   displayed_rama_type = type;

} 

void
coot::rama_plot::all_plot(clipper::Ramachandran::TYPE type) {

   //
   all_plot_background_and_bits(type);
   if (drawing_differences) {
      // draws top n_diffs differences.
      draw_phi_psi_differences();
   } else {
      for (unsigned int ich=0; ich<phi_psi_sets.size(); ich++) 
	 draw_phi_psi_points(ich);
   }

   draw_green_box();
}

void
coot::rama_plot::tooltip_like_box(const mouse_util_t &t) {

   int i = t.smallest_diff_i;
   short int is_secondary = t.mouse_over_secondary_set;

   // make a box slightly below the mouse position
   // (i.e. at highter psi)
   int ich = t.ichain;

   // GtkCanvasItem *item;
   double phi =  phi_psi_sets[ich].phi_psi[i].phi(); 
   double psi = -phi_psi_sets[ich].phi_psi[i].psi();
   const char *lab = phi_psi_sets[ich].phi_psi[i].label().c_str();

   if (is_secondary) {
      phi =  secondary_phi_psi_sets[ich].phi_psi[i].phi(); 
      psi = -secondary_phi_psi_sets[ich].phi_psi[i].psi();
      lab =  secondary_phi_psi_sets[ich].phi_psi[i].label().c_str(); 
   }

   if (tooltip_item)
      gtk_object_destroy(GTK_OBJECT(tooltip_item)); 
   if (tooltip_item_text)
      gtk_object_destroy(GTK_OBJECT(tooltip_item_text));

   // cout << phi_psi[i] << endl;

   // Now we deal with some non-straightforward coordinates.  We want
   // a box with text label overlaying.  However, the coordinates of
   // the box (in world coordintes (i.e. phi/psi space)) depend on the
   // zooming.  i.e.we don't want the box 15 degrees below a residue
   // when the zoom is 4.  Similarly, the size of the box should not
   // be 30-10. It need to be less with more zooming.
   //
   // tw is in pixels not in world coordinates. 
   
   double tw=75; 
   tooltip_item_text = gtk_canvas_item_new(gtk_canvas_root(canvas),
					   GTK_CANVAS_TYPE_CANVAS_TEXT,
					   "text", lab,
					   "x", phi-tw/2,
					   "y", psi+15,
					   "anchor",GTK_ANCHOR_WEST,
					   "font", fixed_font_str.c_str(),
					   "fill_color", "black",
					   NULL);
   gtk_object_get(GTK_OBJECT(tooltip_item_text),
		  "text_width", &tw,
		  NULL);
   gtk_object_destroy(GTK_OBJECT(tooltip_item_text)); 

   tooltip_item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				GTK_CANVAS_TYPE_CANVAS_RECT,
				"x1", phi-(tw/(2))-2/zoom,
				"y1", psi+12/zoom,
				"x2", phi+(tw/(2))+2/zoom,
				"y2", psi+30/zoom,
				"fill_color", "PaleGreen",
				"outline_color", "black",
					NULL);

   tooltip_item_text = gtk_canvas_item_new(gtk_canvas_root(canvas),
					     GTK_CANVAS_TYPE_CANVAS_TEXT,
					     "text", lab,
					     "x", phi-tw/2.0,
					     "y", psi+20.0/zoom,
					     "anchor",GTK_ANCHOR_WEST,
					     "font", fixed_font_str.c_str(),
					     "fill_color", "black",
					     NULL);
}



// The helper type contains a flag the signals validity (this residue
// existed in the both the first and seconde molecules).
//
coot::util::phi_psi_pair_helper_t
coot::rama_plot::make_phi_psi_pair(CMMDBManager *mol1,
				   CMMDBManager *mol2,
				   const std::string &chain_id1,
				   const std::string &chain_id2,
				   int i_seq_num,
				   const std::string &ins_code) const {
   PCResidue *SelResidue1;
   PCResidue *SelResidue2;
   int nSelResidues1, nSelResidues2;
   coot::util::phi_psi_pair_helper_t r;
   r.is_valid_pair_flag = 0;
   
   int selHnd1 = mol1->NewSelection();
   mol1->Select ( selHnd1, STYPE_RESIDUE, 1, // .. TYPE, iModel
		  (char *) chain_id1.c_str(), // Chain id
		  i_seq_num-1,"*",  // starting res
		  i_seq_num+1,"*",  // ending res
		  "*",  // residue name
		  "*",  // Residue must contain this atom name?
		  "*",  // Residue must contain this Element?
		  "*",  // altLocs
		  SKEY_NEW // selection key
		  );
   mol1->GetSelIndex (selHnd1, SelResidue1, nSelResidues1);
   int i_seq_num2 = i_seq_num;
   if (allow_seqnum_offset_flag)
      i_seq_num2 = get_seqnum_2(i_seq_num);
   if (nSelResidues1 == 3) {
      int selHnd2 = mol2->NewSelection();
      mol2->Select ( selHnd2, STYPE_RESIDUE, 1, // .. TYPE, iModel
		     (char *) chain_id2.c_str(), // Chain id
		     i_seq_num2-1,"*",  // starting res
		     i_seq_num2+1,"*",  // ending res
		     "*",  // residue name
		     "*",  // Residue must contain this atom name?
		     "*",  // Residue must contain this Element?
		     "*",  // altLocs
		     SKEY_NEW // selection key
		     );
      mol2->GetSelIndex (selHnd2, SelResidue2, nSelResidues2);
      if (nSelResidues2 == 3) {
	 std::pair<bool, coot::util::phi_psi_t> phi_psi_info_1 = get_phi_psi(SelResidue1);
	 if (phi_psi_info_1.first) {
	    std::pair<bool, coot::util::phi_psi_t> phi_psi_info_2 = get_phi_psi(SelResidue2);
	    if (phi_psi_info_2.first) {
	       bool is_valid = 1;
	       r = coot::util::phi_psi_pair_helper_t(phi_psi_info_1.second,
						     phi_psi_info_2.second,
						     is_valid);
	    } else {
	       // std::cout << "didn't get good phi_psi_info_2 " << std::endl;
	    } 
	 } else {
	    // std::cout << "didn't get good phi_psi_info_1 " << std::endl;
	 } 
      } 
      mol2->DeleteSelection(selHnd2);
   }
   mol1->DeleteSelection(selHnd1);
   return r;
}


// SelResidue is guaranteed to have 3 residues (there is no protection
// for that in this function).
std::pair<bool, coot::util::phi_psi_t> coot::rama_plot::get_phi_psi(PCResidue *SelResidue) const {
   return coot::util::get_phi_psi(SelResidue);
}




// add to the phi_psi vector
void
coot::rama_plot::add_phi_psi(vector <coot::util::phi_psi_t> *phi_psi_vec,
			     CMMDBManager *mol_in, const char *segid,
			     CResidue *prev, CResidue *this_res, CResidue *next_res) {


   // we need the C (for phi) of the previous residue and the N of the
   // next one (for psi).
   // 

   //    CMMDBManager *mol = (CMMDBManager *)mol_in;
   PCAtom *res_selection;
   int i_no_res_atoms;
   
   int natom = 0;
   clipper::Coord_orth c_prev, n_this, ca_this, c_this, n_next;

   if (prev && this_res && next_res) { 

      // prev:
      
      prev->GetAtomTable(res_selection, i_no_res_atoms);
      if (i_no_res_atoms > 0) {
	 for (int j=0; j<i_no_res_atoms; j++) {
	    if ( std::string(res_selection[j]->name) == " C  ") {
	       c_prev = clipper::Coord_orth(res_selection[j]->x,
					    res_selection[j]->y,
					    res_selection[j]->z);
	       natom++;
	    }
	 }
      }

      // this:
      this_res->GetAtomTable(res_selection, i_no_res_atoms);
      if (i_no_res_atoms > 0) {
	 for (int j=0; j<i_no_res_atoms; j++) {
	    if ( std::string(res_selection[j]->name) == " C  ") {
	       c_this = clipper::Coord_orth(res_selection[j]->x,
					    res_selection[j]->y,
					    res_selection[j]->z);
	       natom++;
	    }
	    if ( std::string(res_selection[j]->name) == " CA ") {
	       ca_this = clipper::Coord_orth(res_selection[j]->x,
					     res_selection[j]->y,
					     res_selection[j]->z);
	       natom++;
	    }
	    if ( std::string(res_selection[j]->name) == " N  ") {
	       n_this = clipper::Coord_orth(res_selection[j]->x,
					    res_selection[j]->y,
					    res_selection[j]->z);
	       natom++;
	    }
	 }
      }

      // next
      next_res->GetAtomTable(res_selection, i_no_res_atoms);
      if (i_no_res_atoms > 0) {
	 for (int j=0; j<i_no_res_atoms; j++) {
	    if (std::string(res_selection[j]->name) == " N  ") {
	       n_next = clipper::Coord_orth(res_selection[j]->x,
					    res_selection[j]->y,
					    res_selection[j]->z);
	       natom++;
	    }
	 }
      }

      // So that we don't accidently join fragments that should not be
      // linked, we do a distance check too.
      
      if (natom == 5) {
	 float max_dist = 2.0; // A
	 if ((clipper::Coord_orth::length(c_prev, n_this) <= max_dist) &&
	     (clipper::Coord_orth::length(c_this, n_next) <= max_dist)) {


	    std::string label = coot::util::int_to_string(this_res->GetSeqNum());
	    label += this_res->GetInsCode();
	    label += " ";
	    label += this_res->name;
	    label += " ";      
	    label += segid;
	    std::string inscode = this_res->GetInsCode();
	    
	    // std::cout << " constructed label " << label << std::endl;
	    double phi = clipper::Util::rad2d(ca_this.torsion(c_prev, n_this,
							      ca_this, c_this));
	    double psi = clipper::Util::rad2d(ca_this.torsion(n_this, ca_this,
							      c_this, n_next));
	    phi_psi_vec->push_back(coot::util::phi_psi_t(phi, psi,
							 this_res->name,
							 label.c_str(),
							 this_res->GetSeqNum(),
							 inscode,
							 segid));
	 }
      }
   }
}



void
coot::rama_plot::draw_phi_psis_on_canvas(char *filename) {

   
   // Do the background first, of course.
   //
   display_background();

   draw_axes();
   draw_zero_lines(); 

   // or however you want to get your mmdbmanager
   // 
   CMMDBManager *mol = rama_get_mmdb_manager(filename); 

   // put the results in the *primary* list
   generate_phi_psis(&phi_psi_sets, mol);

   // generate_phi_psis_debug(); 
   for (unsigned int ich=0; ich<phi_psi_sets.size(); ich++) 
      draw_phi_psi_points(ich);
   
   delete mol; 
}


void
coot::rama_plot::draw_it(CMMDBManager *mol) {

   display_background();
   draw_axes();
   draw_zero_lines(); 
   generate_phi_psis(&phi_psi_sets, mol);
   coot::rama_stats_container_t counts;
   for (unsigned int ich=0; ich<phi_psi_sets.size(); ich++) 
      counts += draw_phi_psi_points(ich);
   counts_to_stats_frame(counts);
}

void
coot::rama_plot::draw_it(int imol1, int imol2,
			 CMMDBManager *mol1, CMMDBManager *mol2) {

   molecule_numbers_ = std::pair<int, int> (imol1, imol2); // save for later
   display_background();
   draw_axes();
   draw_zero_lines(); 
   draw_2_phi_psi_sets_on_canvas(mol1, mol2);
   if (is_kleywegt_plot())
      hide_stats_frame();
}

// the CMMDBManager could have gone out of date when we come to redraw
// the widget after refinement, so we pass the imol1 and imol2 so that
// we can ask globjects if the molecule numbers are still valid.
// 
void
coot::rama_plot::draw_it(int imol1, int imol2,
			 CMMDBManager *mol1, CMMDBManager *mol2,
			 const std::string &chain_id_1, const std::string &chain_id_2) {

   molecule_numbers_ = std::pair<int, int> (imol1, imol2); // save for later
   chain_ids_ = std::pair<std::string, std::string> (chain_id_1, chain_id_2);
   display_background();
   draw_axes();
   draw_zero_lines();
   if (allow_seqnum_offset_flag)
      set_seqnum_offset(imol1, imol2, mol1, mol2, chain_id_1, chain_id_2);
   draw_2_phi_psi_sets_on_canvas(mol1, mol2, chain_id_1, chain_id_2);
   if (is_kleywegt_plot())
      hide_stats_frame();
}

// Was from a shelx molecule with A 1->100 and B 201->300.
// For shelx molecule as above, what do we need to add to seqnum_1 to get the
// corresponding residue in the B chain (in the above example it is 100).
// 
void
coot::rama_plot::set_seqnum_offset(int imol1, int imol2,
				   CMMDBManager *mol1,
				   CMMDBManager *mol2,
				   const std::string &chain_id_1,
				   const std::string &chain_id_2) {

   seqnum_offset = MinInt4;
   if (is_valid_model_molecule(imol1)) { 
      if (is_valid_model_molecule(imol2)) {

	 int imod = 1;
      
	 CModel *model_p_1 = mol1->GetModel(imod);
	 CChain *chain_p_1;
	 // run over chains of the existing mol
	 int nchains_1 = model_p_1->GetNumberOfChains();
	 for (int ichain_1=0; ichain_1<nchains_1; ichain_1++) {
	    chain_p_1 = model_p_1->GetChain(ichain_1);
	    if (chain_id_1 == chain_p_1->GetChainID()) {
	       int nres_1 = chain_p_1->GetNumberOfResidues();
	       PCResidue residue_p_1;
	       if (nres_1 > 0) {
		  residue_p_1 = chain_p_1->GetResidue(0);
		  CModel *model_p_2 = mol2->GetModel(imod);
		  CChain *chain_p_2;
		  int nchains_2 = model_p_2->GetNumberOfChains();
		  for (int ichain_2=0; ichain_2<nchains_2; ichain_2++) {
		     chain_p_2 = model_p_2->GetChain(ichain_2);
		     if (chain_id_2 == chain_p_2->GetChainID()) {
			int nres_2 = chain_p_2->GetNumberOfResidues();
			PCResidue residue_p_2;
			if (nres_2 > 0) {
			   residue_p_2 = chain_p_2->GetResidue(0);

			   seqnum_offset = residue_p_2->GetSeqNum() - residue_p_1->GetSeqNum();
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   std::cout << "DEBUG:: seqnum_offset is: " << seqnum_offset << std::endl;
   
   if (seqnum_offset == MinInt4) {
      std::cout << "WARNING:: Ooops! Failed to set the Chain Residue numbering different\n"
		<< "WARNING::        offset correctly." << std::endl;
      std::cout << "WARNING:: Ooops! Bad Kleywegts will result!" << std::endl;
      seqnum_offset = 0;
   }
} 

int
coot::rama_plot::get_seqnum_2(int seqnum_1) const {

   return seqnum_1 + seqnum_offset;

}

void
coot::rama_plot::hide_stats_frame() {

   if (canvas) {
      GtkWidget *frame = lookup_widget(GTK_WIDGET(canvas), "rama_stats_frame");
      gtk_widget_hide(frame);
   } else {
      std::cout << "ERROR:: null widget in hide_stats_frame\n";
   } 
}

void
coot::rama_plot::counts_to_stats_frame(const coot::rama_stats_container_t &sc) {

   if (sc.n_ramas > 0) { 
      float pref_frac = float(sc.n_preferred)/float(sc.n_ramas);
      float allow_frac = float(sc.n_allowed)/float(sc.n_ramas);
      int n_outliers = sc.n_ramas - sc.n_preferred - sc.n_allowed;
      float outlr_frac = float(n_outliers)/float(sc.n_ramas);

      std::string pref_str = "In Preferred Regions:  ";
      pref_str += coot::util::int_to_string(sc.n_preferred);
      pref_str += "  (";
      pref_str += coot::util::float_to_string(100.0*pref_frac);
      pref_str += "%)";
	 
      std::string allow_str = "In Allowed Regions:  ";
      allow_str += coot::util::int_to_string(sc.n_allowed);
      allow_str += "  (";
      allow_str += coot::util::float_to_string(100.0*allow_frac);
      allow_str += "%)";
	 
      std::string outlr_str = "Outliers:  ";
      outlr_str += coot::util::int_to_string(n_outliers);
      outlr_str += "  (";
      outlr_str += coot::util::float_to_string(100.0*outlr_frac);
      outlr_str += "%)";

      GtkWidget *label1 = lookup_widget(GTK_WIDGET(canvas), "rama_stats_label_1");
      GtkWidget *label2 = lookup_widget(GTK_WIDGET(canvas), "rama_stats_label_2");
      GtkWidget *label3 = lookup_widget(GTK_WIDGET(canvas), "rama_stats_label_3");

      gtk_label_set_text(GTK_LABEL(label1),  pref_str.c_str());
      gtk_label_set_text(GTK_LABEL(label2), allow_str.c_str());
      gtk_label_set_text(GTK_LABEL(label3), outlr_str.c_str());
	 
   } else {
      hide_stats_frame();
   } 
} 


void
coot::rama_plot::draw_it(const coot::util::phi_psi_t &phipsi) {

  display_background();
  draw_axes();
  draw_zero_lines();
  phi_psi_sets.resize(1);
  phi_psi_sets[0].phi_psi.push_back(phipsi);
  for (unsigned int ich=0; ich<phi_psi_sets.size(); ich++) 
      draw_phi_psi_points(ich);
 
}

void
coot::rama_plot::draw_it(const std::vector<coot::util::phi_psi_t> &phipsi) {

//   display_background();
//   draw_axes();
//   draw_zero_lines();
  phi_psi_sets.resize(1);
  phi_psi_sets[0].phi_psi.clear();
  clear_last_canvas_items(phipsi.size());
  for (unsigned int ipp=0; ipp<phipsi.size(); ipp++) { 
     phi_psi_sets[0].phi_psi.push_back(phipsi[ipp]);
  }

  // std::cout << "DEBUG:: phi_psi_sets.size(): " << phi_psi_sets.size() << std::endl;
  for (unsigned int ich=0; ich<phi_psi_sets.size(); ich++)   
     // std::cout << "DEBUG:: " << ich << " " << phi_psi_sets[ich].phi_psi.size() << std::endl;

  for (unsigned int ich=0; ich<phi_psi_sets.size(); ich++) 
      draw_phi_psi_points(ich);
}


void
coot::rama_plot::draw_2_phi_psi_sets_on_canvas(char *file1, char *file2) {

   CMMDBManager *mol1 = rama_get_mmdb_manager(file1); 
   CMMDBManager *mol2 = rama_get_mmdb_manager(file2); 

   draw_2_phi_psi_sets_on_canvas(mol1, mol2);
   delete mol1; 
   delete mol2;
}

void 
coot::rama_plot::draw_2_phi_psi_sets_on_canvas(CMMDBManager *mol1, 
					       CMMDBManager *mol2) { 

//    generate_phi_psis(&phi_psi_sets, mol1);
//    generate_phi_psis(&secondary_phi_psi_sets, mol2);


   // Are you sure that you want to edit this function? (not the one below?)

   int imod = 1;
      
   CModel *model_p = mol1->GetModel(imod);
   CChain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   int set_index = 0;
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      PCResidue residue_p;
      phi_psi_sets.push_back(coot::phi_psi_set_container());
      secondary_phi_psi_sets.push_back(coot::phi_psi_set_container());
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 std::string chain_id = residue_p->GetChainID();
	 int res_num = residue_p->GetSeqNum();
	 std::string inscode = residue_p->GetInsCode();
	 coot::util::phi_psi_pair_helper_t p_i =
	    coot::rama_plot::make_phi_psi_pair(mol1, mol2, chain_id, chain_id, res_num, inscode);

	 if (p_i.is_valid_pair_flag) {
	    phi_psi_sets[set_index].phi_psi.push_back(p_i.first);
	    secondary_phi_psi_sets[set_index].phi_psi.push_back(p_i.second);
	 }
      }
      set_index++;
   }
   
   // all_plot_background_and_bits(clipper::Ramachandran::All);
   // std::cout << "finding differences..." << std::endl;
   find_phi_psi_differences_whole_molecule(); 
   // std::cout << "drawing differences..." << std::endl;
   draw_phi_psi_differences();

   // do we need to do something like this?
//    for (int ich=0; ich<phi_psi_sets.size(); ich++) 
//       draw_phi_psi_points(ich);

   drawing_differences = 1; // set flag for later drawing
}

void 
coot::rama_plot::draw_2_phi_psi_sets_on_canvas(CMMDBManager *mol1, 
					       CMMDBManager *mol2,
					       std::string chainid1,
					       std::string chainid2) {

   // std::cout << "--------------- draw 2 on canvas -------------" << std::endl;

   int imod = 1;

   phi_psi_sets.clear();
   secondary_phi_psi_sets.clear();
      
   CModel *model_p = mol1->GetModel(imod);
   CChain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   int set_index = 0;
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      std::string this_chain = chain_p->GetChainID();
      PCResidue residue_p;
      phi_psi_sets.push_back(coot::phi_psi_set_container());
      secondary_phi_psi_sets.push_back(coot::phi_psi_set_container());
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 std::string chain_id = residue_p->GetChainID();
	 int res_num = residue_p->GetSeqNum();
	 std::string inscode = residue_p->GetInsCode();
	 coot::util::phi_psi_pair_helper_t p_i =
	   coot::rama_plot::make_phi_psi_pair(mol1, mol2, chainid1, chainid2, res_num, inscode);

	 if (p_i.is_valid_pair_flag) {
// 	    std::cout << "Adding pair for " << p_i.first.label() << " "
// 		      << p_i.second.label() << std::endl;
	    phi_psi_sets[set_index].phi_psi.push_back(p_i.first);
	    secondary_phi_psi_sets[set_index].phi_psi.push_back(p_i.second);
	 } else {
	    // std::cout << "pair for " <<  p_i.first.label() << " was not valid!" << std::endl;
	 }
      }
      set_index++;
   }
   
   // all_plot_background_and_bits(clipper::Ramachandran::All);
   std::cout << "finding differences..." << std::endl;
   find_phi_psi_differences_by_chain(); 
   std::cout << "drawing differences..." << std::endl;
   draw_phi_psi_differences();

   drawing_differences = 1; // set flag for later drawing
}


void
coot::rama_plot::find_phi_psi_differences_whole_molecule() {

   double d1, d2;

   diff_sq.clear();
   for (unsigned int ich=0; ich<phi_psi_sets.size(); ich++) { 
      for (unsigned int ires=0; ires<phi_psi_sets[ich].phi_psi.size(); ires++) {

	 // std::cout << "ich: " << ich << " i: " << i << std::endl;
	 
	 d1 = fabs(phi_psi_sets[ich].phi_psi[ires].phi() -
		  secondary_phi_psi_sets[ich].phi_psi[ires].phi()); 
	 d2 = fabs(phi_psi_sets[ich].phi_psi[ires].psi() -
		  secondary_phi_psi_sets[ich].phi_psi[ires].psi());
	 
	 if (d1 > 180.0) d1 = 360.0 - d1;
	 if (d2 > 180.0) d2 = 360.0 - d2;
	 
	 diff_sq.push_back(diff_sq_t(ires, sqrt(d1*d1+d2*d2), ich, ich));
      }
   
   }
   // sort diff_sq
   std::sort(diff_sq.begin(), diff_sq.end(),
	     compare_phi_psi_diffs(diff_sq));

}

void
coot::rama_plot::find_phi_psi_differences_by_chain() {

   double d1, d2;

   diff_sq.clear();
   int ich = 0;
   for (unsigned int ires=0; ires<phi_psi_sets[ich].phi_psi.size(); ires++) {
      
      // std::cout << "ich: " << ich << " i: " << i << std::endl;
      
      d1 = fabs(phi_psi_sets[ich].phi_psi[ires].phi() -
		secondary_phi_psi_sets[ich].phi_psi[ires].phi()); 
      d2 = fabs(phi_psi_sets[ich].phi_psi[ires].psi() -
		secondary_phi_psi_sets[ich].phi_psi[ires].psi());
      
      if (d1 > 180.0) d1 = 360.0 - d1;
      if (d2 > 180.0) d2 = 360.0 - d2;
      
      diff_sq.push_back(diff_sq_t(ires, sqrt(d1*d1+d2*d2), ich, ich));
   }
   // sort diff_sq
   std::sort(diff_sq.begin(), diff_sq.end(),
	     compare_phi_psi_diffs(diff_sq));

//    // debug
//    for (int id=0; id<diff_sq.size(); id++) {
//       std::cout << "diffs: " << id << " "
// 		<< phi_psi_sets[ich].phi_psi[diff_sq[id].i()].residue_number << " "
// 		<< phi_psi_sets[ich].phi_psi[diff_sq[id].i()].chain_id << "   "
// 		<< diff_sq[id].i() << " " << diff_sq[id].v() 
// 		<< std::endl;
//    }
}

// n_diffs is the max number of differences to be displayed.
// 
void
coot::rama_plot::draw_phi_psi_differences() {


   // First calculate the differences
   //
   // then sort the differences
   //
   // then plot the top n differences
   // 
   // The first point should be hollowed out (white centre), then an
   // arrow should point to a conventionally drawn point from the
   // second set.
   //
   // So we will need to change things so that the drawing routine
   // will know whether to put a white interior or not.
   //
   // Note that glycines will not be affected since they are made from
   // (canvas) lines, not rectangles.
   // 
   // Maybe I could make the starting
   // triangle be drawn in grey in that case.

   int ich = 0;

   if (phi_psi_sets.size() == 0 ) {

      std::cout << "oops! phi psi data data!" << std::endl;

   } else { 

      if (phi_psi_sets[ich].phi_psi.size() != secondary_phi_psi_sets[ich].phi_psi.size()) {
	 std::cout << "wrong number of residues in draw_phi_psi_differences "
		   << phi_psi_sets[ich].phi_psi.size() << " "
		   << secondary_phi_psi_sets[ich].phi_psi.size() << std::endl;

      } else { 

// 	 std::cout << "--- DEBUG:: running over " << n_diffs
// 		   << " diffs with diffs.size(): " << diff_sq.size()
// 		   << std::endl;

	 GtkCanvasPoints *points = gtk_canvas_points_new(2);
	 int i; 
      
	 for (int j=0; j<int(diff_sq.size()) && j<n_diffs; j++) {

	 
	    i = diff_sq[j].i();
	    ich = diff_sq[j].ich();
	    // std::cout << i << " " << diff_sq[j].v() << "   ";

	    draw_phi_psi_point(i, &phi_psi_sets[ich].phi_psi, 1);
	    draw_phi_psi_point(i, &secondary_phi_psi_sets[ich].phi_psi, 0);

	    draw_kleywegt_arrow(phi_psi_sets[ich].phi_psi[i],
				secondary_phi_psi_sets[ich].phi_psi[i],
				points);
	    
	 }
      }
   }
}


void
coot::rama_plot::draw_kleywegt_arrow(const coot::util::phi_psi_t &phi_psi_primary,
				     const coot::util::phi_psi_t &phi_psi_secondary,
				     GtkCanvasPoints *points) {

   
   coot::rama_kleywegt_wrap_info wi = test_kleywegt_wrap(phi_psi_primary,
							 phi_psi_secondary);
   
   GtkCanvasItem *item;

   if (wi.is_wrapped == 0) {
      // the normal case
      points->coords[0] =  phi_psi_primary.phi(); 
      points->coords[1] = -phi_psi_primary.psi();
      
      points->coords[2] =  phi_psi_secondary.phi(); 
      points->coords[3] = -phi_psi_secondary.psi();
      item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				 GTK_CANVAS_TYPE_CANVAS_LINE,
				 "width_pixels", 1,
				 "points", points,
				 "last_arrowhead", 1,
				 "arrow_shape_a", 5.0, 
				 "arrow_shape_b", 8.0, 
				 "arrow_shape_c", 4.0, 
				 "fill_color", "black",
				 NULL);
      canvas_item_vec.push_back(item);
   } else {
      // border crosser:
      
      // line to the border
      points->coords[0] =  phi_psi_primary.phi(); 
      points->coords[1] = -phi_psi_primary.psi();
      
      points->coords[2] =  wi.primary_border_point.first;
      points->coords[3] = -wi.primary_border_point.second;

      std::cout << "Borderline 1 : "
		<< "(" << phi_psi_primary.phi() << "," << phi_psi_primary.psi() << ")"
		<< " to "
		<< "("
		<< wi.primary_border_point.first
		<< ", " << wi.primary_border_point.second
		<< ")" << std::endl;
      
      item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				 GTK_CANVAS_TYPE_CANVAS_LINE,
				 "width_pixels", 1,
				 "points", points,
				 // "last_arrowhead", 1,
// 				 "arrow_shape_a", 5.0, 
// 				 "arrow_shape_b", 8.0, 
// 				 "arrow_shape_c", 4.0, 
				 "fill_color", "black",
				 NULL);
      canvas_item_vec.push_back(item);
      
      // line from the border
      
      points->coords[0] =  wi.secondary_border_point.first;
      points->coords[1] = -wi.secondary_border_point.second;
      
      points->coords[2] =  phi_psi_secondary.phi(); 
      points->coords[3] = -phi_psi_secondary.psi();

      std::cout << "Borderline 2: "
		<< "(" << wi.secondary_border_point.first << ","
		<< wi.secondary_border_point.second << ")"
		<< " to " 
		<< "("  << phi_psi_secondary.phi()
		<< ", " << phi_psi_secondary.psi()
		<< ")" << std::endl;
      
      item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				 GTK_CANVAS_TYPE_CANVAS_LINE,
				 "width_pixels", 1,
				 "points", points,
				 "last_arrowhead", 1,
				 "arrow_shape_a", 5.0, 
				 "arrow_shape_b", 8.0, 
				 "arrow_shape_c", 4.0, 
				 "fill_color", "black",
				 NULL);
      canvas_item_vec.push_back(item);
      
   } 
}

coot::rama_kleywegt_wrap_info
coot::rama_plot::test_kleywegt_wrap(const coot::util::phi_psi_t &phi_psi_primary,
				    const coot::util::phi_psi_t &phi_psi_secondary) const {

   coot::rama_kleywegt_wrap_info wi;


   if (fabs(phi_psi_primary.phi() - phi_psi_secondary.phi()) > 200.0) {
       std::cout << "DEBUG:: PHI Outlier detected: " << phi_psi_primary.chain_id
		 << " " << phi_psi_primary.residue_number << std::endl;
      wi.is_wrapped = 1;
      
      float phi_1 = phi_psi_primary.phi();
      float psi_1 = phi_psi_primary.psi();
      float phi_2 = phi_psi_secondary.phi();
      float psi_2 = phi_psi_secondary.psi();

      float psi_diff = psi_2 - psi_1; 
      float psi_gradient = 999999999.9;
      if (fabs(psi_diff) > 0.000000001)
	 psi_gradient = (180.0 - phi_1)/(phi_2 + 360.0 - phi_1);

      float psi_critical = psi_1 + psi_gradient * (psi_2 - psi_1);
      wi.primary_border_point.first = 180.0;
      wi.primary_border_point.second = psi_critical;
      wi.secondary_border_point.first = -180.0;
      wi.secondary_border_point.second = psi_critical;
   } 
   if (fabs(phi_psi_primary.psi() - phi_psi_secondary.psi()) > 200.0) {
      std::cout << "DEBUG:: PSI Outlier detected: " << phi_psi_primary.chain_id
		<< " " << phi_psi_primary.residue_number << std::endl;
      wi.is_wrapped = 1;

      float phi_1 = phi_psi_primary.phi();
      float psi_1 = phi_psi_primary.psi();
      float phi_2 = phi_psi_secondary.phi();
      float psi_2 = phi_psi_secondary.psi();

      float psi_diff = psi_2 - psi_1; 
      float psi_gradient = 999999999.9;
      if (fabs(psi_diff) > 0.000000001)
	 psi_gradient = (-180.0 - (psi_2 - 360.0))/(psi_1 - (psi_2 - 360.0));

      float phi_critical = phi_2 + psi_gradient * (phi_1 - phi_2);
      wi.primary_border_point.first = phi_critical;
      wi.primary_border_point.second = -180.0;
      wi.secondary_border_point.first = phi_critical;
      wi.secondary_border_point.second = 180.0;
   } 
   return wi;
} 


void
coot::rama_plot::draw_zero_lines() {

   GtkCanvasPoints *points = gtk_canvas_points_new(2);

   points->coords[0] = 0.0;
   points->coords[1] = -180.0; 

   points->coords[2] = 0.0; 
   points->coords[3] = 180.0; 

   // consider also:
   // fill_stipple		GdkBitmap*
   // 				"line_style", 2,
   GtkCanvasItem *item =
      gtk_canvas_item_new(gtk_canvas_root(canvas),
			    GTK_CANVAS_TYPE_CANVAS_LINE,
			    "width_pixels", 1,
			    // 			    "line_style", 2, 
			    "points", points,
			    "fill_color", "grey",
			    NULL);

   canvas_item_vec.push_back(item);

   points->coords[0] = -180.0;
   points->coords[1] = 0.0;

   points->coords[2] = 180.0;
   points->coords[3] = 0.0;

   item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				GTK_CANVAS_TYPE_CANVAS_LINE,
				"width_pixels", 1,
				"points", points,
				"fill_color", "grey",
				NULL);

   canvas_item_vec.push_back(item);
   gtk_canvas_points_free(points); 
}

// Tick marks and text labels for phi and phi axes.
// 
void
coot::rama_plot::draw_axes() {

   // First do the text for the axes labels.
   //

   GtkCanvasItem *item;
   //  line_style		GdkLineStyle
   item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				GTK_CANVAS_TYPE_CANVAS_TEXT,
				"text", "Phi", 
				"x", -10.0,
				"y", 220.0,
				"anchor",GTK_ANCHOR_WEST,
				"font", fixed_font_str.c_str(),
				"fill_color", "black",
				NULL);
   canvas_item_vec.push_back(item);

   item = gtk_canvas_item_new(gtk_canvas_root(canvas),
			      GTK_CANVAS_TYPE_CANVAS_TEXT,
			      "text", "Psi",
			      "x", -230.0,
			      "y", 10.0,
			      "anchor",GTK_ANCHOR_WEST,
			      "font", fixed_font_str.c_str(),
			      "fill_color", "black",
			      NULL);
   canvas_item_vec.push_back(item);

   // Ticks
   vector<canvas_tick_t> pnts;

   // x axis
   pnts.push_back(canvas_tick_t(0,-180.0,180.0));
   pnts.push_back(canvas_tick_t(0,-120.0,180.0));
   pnts.push_back(canvas_tick_t(0,-60.0,180.0));
   pnts.push_back(canvas_tick_t(0,0.0,180.0));
   pnts.push_back(canvas_tick_t(0,60.0,180.0));
   pnts.push_back(canvas_tick_t(0,120.0,180.0));
   pnts.push_back(canvas_tick_t(0,180.0,180.0));
   

   // y axis
    pnts.push_back(canvas_tick_t(1,-180.0,-180.0));
    pnts.push_back(canvas_tick_t(1,-180.0,-120.0));
    pnts.push_back(canvas_tick_t(1,-180.0,-60.0));
    pnts.push_back(canvas_tick_t(1,-180.0,0.0));
    pnts.push_back(canvas_tick_t(1,-180.0,60.0));
    pnts.push_back(canvas_tick_t(1,-180.0,120.0));
    pnts.push_back(canvas_tick_t(1,-180.0,180.0));
   
   GtkCanvasPoints *points = gtk_canvas_points_new(2);

   for (unsigned int i=0; i<pnts.size(); i++) { 
      points->coords[0] = pnts[i].start_x();
      points->coords[1] = pnts[i].start_y();

      points->coords[2] = pnts[i].end_x();
      points->coords[3] = pnts[i].end_y();

      item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				   GTK_CANVAS_TYPE_CANVAS_LINE,
				   "width_pixels", 1,
				   "points", points,
				   "fill_color", "black",
				   NULL);
      canvas_item_vec.push_back(item);

   }

   // Ticks text

   // x axis

   vector<int> tick_text;
   tick_text.push_back(-180); 
   tick_text.push_back(-120); 
   tick_text.push_back(-60); 
   tick_text.push_back(0); 
   tick_text.push_back(60); 
   tick_text.push_back(120); 
   tick_text.push_back(180);
   char text[20];

   for (unsigned int i=0; i<tick_text.size(); i++) {

      snprintf(text,19,"%d",tick_text[i]);
   
      item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				   GTK_CANVAS_TYPE_CANVAS_TEXT,
				   "text", text,
				   "x", -220.0,
				   "y", -tick_text[i] +0.0,
				   "anchor",GTK_ANCHOR_WEST,
				 "font", fixed_font_str.c_str(),
				   "fill_color", "black",
				   NULL);

      canvas_item_vec.push_back(item);

      
//    // y axis

      item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				   GTK_CANVAS_TYPE_CANVAS_TEXT,
				   "text", text,
				   "x", tick_text[i] - 10.0,
				   "y", 200.0,
				   "anchor", GTK_ANCHOR_WEST,
				   "font", fixed_font_str.c_str(),
				   "fill_color", "black",
				   NULL);
      canvas_item_vec.push_back(item);

   }
   gtk_canvas_points_free(points); 
} 




CMMDBManager *
coot::rama_plot::rama_get_mmdb_manager(std::string pdb_name) {

   int err;
   CMMDBManager* MMDBManager;

   // Needed for the error message printing: 
   // MMDBManager->GetInputBuffer(S, lcount);
   // Used by reference and as a pointer.  Grimness indeed.
   int  error_count;
   char error_buf[500];

   //   Make routine initializations
   //
   InitMatType();

   MMDBManager = new CMMDBManager;

   cout << "Reading coordinate file: " << pdb_name.c_str() << "\n";
   err = MMDBManager->ReadCoorFile((char *)pdb_name.c_str());
   
   if (err) {
      // does_file_exist(pdb_name.c_str());
      cout << "There was an error reading " << pdb_name.c_str() << ". \n";
      cout << "ERROR " << err << " READ: "
	   << GetErrorDescription(err) << endl;
      //
      // This makes my stomach churn too. Sorry.
      // 
      MMDBManager->GetInputBuffer(error_buf, error_count);
      if (error_count >= 0) { 
	 cout << "         LINE #" << error_count << "\n     "
	      << error_buf << endl << endl;
      } else {
	 if (error_count == -1) { 
	    cout << "       CIF ITEM: " << error_buf << endl << endl;
	 }
      }

      //
   } else {
      // we read the coordinate file OK.
      //
      switch (MMDBManager->GetFileType())  {
      case MMDB_FILE_PDB    :  cout << " PDB"         ;
	 break;
      case MMDB_FILE_CIF    :  cout << " mmCIF"       ; 
	 break;
      case MMDB_FILE_Binary :  cout << " MMDB binary" ;
	 break;
      default:
	 cout << " Unknown (report as a bug!)\n";
      }
      cout << " file " << pdb_name.c_str() << " has been read.\n";
   }
   
    return MMDBManager;
}

void
coot::rama_plot::zoom_out() {

   float canvas_scale = zoom/(zoom + 0.2);
   zoom -= 0.2;
   gtk_canvas_set_pixels_per_unit(canvas,zoom);
   float wf = GTK_WIDGET(canvas)->allocation.width;
   float hf = GTK_WIDGET(canvas)->allocation.height;
   gtk_widget_set_usize(GTK_WIDGET(canvas),
			int(wf*canvas_scale), int(hf*canvas_scale));


}
void
coot::rama_plot::zoom_in() {

   float canvas_scale = (zoom + 0.2)/zoom;
   zoom += 0.2;
   gtk_canvas_set_pixels_per_unit(canvas,zoom);
   float wf = GTK_WIDGET(canvas)->allocation.width;
   float hf = GTK_WIDGET(canvas)->allocation.height;
   gtk_widget_set_usize(GTK_WIDGET(canvas),
			int(wf*canvas_scale), int(hf*canvas_scale));
}

void
coot::rama_plot::allow_seqnum_offset() {

   allow_seqnum_offset_flag = 1;
   // debug();
} 


void
coot::rama_plot::debug() const { 

   std::cout << std::endl;
   std::cout << "ramadebug: imol is " << imol << std::endl;
   std::cout << "ramadebug: canvas is " << canvas << std::endl;
   std::cout << "ramadebug: big_box_item is " << big_box_item << std::endl;
   std::cout << "ramadebug: step is " << step << std::endl;
   std::cout << "ramadebug: is ifirst_res size " << ifirst_res.size() << std::endl;
   std::cout << "ramadebug: is phi_psi_sets.size " << phi_psi_sets.size() << std::endl;
   std::cout << "ramadebug: allow_seqnum_offset_flag " << allow_seqnum_offset_flag << std::endl;
   //    std::cout << "ramadebug: is " << << std::endl;
}


// Notes for the workings of editphipsi
// 
// So, the motion callback function calls map_mouse_pos for this object
// 
// We need also need button_press member function and button_release
// function (we don't have that currently) for use when we "drop" the
// drag.  Hmmm... maybe not, actually, we just leave it where it was
// because the position will be updated on mouse motion.  We need to
// check on mouse motion if the left-button is being pressed (as I
// guess we do in globject.cc's glarea_motion_notify().
//
// So, we *do* need a mechanism to update the moving_atoms in the
// graphics window to the current phi/psi position.  
// 
// We need to know what phi/psi is given the mouse position - I guess
// that there is a function for that already.
// 
// We need to spot when we are over the green box (perhaps with a
// larger "margin of error" than we have with conventional points)
// when the mouse has been clicked.
//
// I think communication between the rama_plot and the graphics can be
// via a phi_psi_t
// 

#endif // RAMA_PLOT
#endif // HAVE_GTK_CANVAS

