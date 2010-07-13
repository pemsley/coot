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

   imol = imol_in; // is this used? (yes, sort of, but not handling
		   // the click on a residue)
   molecule_numbers_.first = imol;
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
coot::rama_plot::big_square(const std::string &chain_id,
			    int resno,
			    const std::string &ins_code) {

   big_square(1, chain_id, resno, ins_code);
   
}

void
coot::rama_plot::big_square(int model_number,
			    const std::string &chain_id,
			    int resno,
			    const std::string &ins_code) {

   coot::residue_spec_t res_spec(chain_id, resno, ins_code);
   if (model_number >= 1) { 
      if (model_number < phi_psi_model_sets.size()) {
	 coot::util::phi_psi_t pp = phi_psi_model_sets[model_number][res_spec];
	 if (pp.is_filled()) {
	    draw_phi_psi_point_internal(pp, 0, 4);
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
coot::rama_plot::draw_phi_psi_point(const coot::util::phi_psi_t &phi_psi,
				    bool as_white_flag) {

   return draw_phi_psi_point_internal(phi_psi, as_white_flag, 2); 
}

int 
coot::rama_plot::draw_phi_psi_point_internal(const coot::util::phi_psi_t &phi_psi,
					     bool as_white_flag,
					     int box_size) {

   int region = coot::rama_plot::RAMA_UNKNOWN; // initially unset
   
   std::string outline_color("black");

   if (box_size == 4) {
      draw_green_box(phi_psi.phi(), phi_psi.psi());
   } else {

      if (phi_psi.residue_name() == "GLY") {
	 region = draw_phi_psi_as_gly(phi_psi);
      } else {
	 std::string colour;
	 double phi = phi_psi.phi();
	 double psi = phi_psi.psi();

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
	    if (phi_psi.residue_name() == "PRO") {
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
coot::rama_plot::draw_phi_psi_as_gly(const coot::util::phi_psi_t &phi_psi) {

   GtkCanvasItem *item;

   // clipper::Ramachandran r_tmp;
//    r_gly.init(clipper::Ramachandran::Gly);
//    r_gly.set_thresholds(0.05,0.002); //defaults 0.01 0.0005

   std::string colour;
   int region;

   double phi = phi_psi.phi();
   double psi = phi_psi.psi();
   
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
coot::rama_plot::draw_phi_psi_points() {

   coot::rama_stats_container_t counts;
   short int as_white_flag = 0;

   for (unsigned int imod=1; imod<phi_psi_model_sets.size(); imod++) {
      counts += draw_phi_psi_points_for_model(phi_psi_model_sets[imod]);
   }
   return counts; 
}


coot::rama_stats_container_t
coot::rama_plot::draw_phi_psi_points_for_model(const coot::phi_psis_for_model_t &pp_set) {

   coot::rama_stats_container_t counts;
   bool as_white_flag = 0;

   std::map<coot::residue_spec_t, coot::util::phi_psi_t>::const_iterator it;
   
   for (it=pp_set.phi_psi.begin(); it!=pp_set.phi_psi.end(); it++) {
      int type = draw_phi_psi_point(it->second, as_white_flag);
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
void
coot::rama_plot::generate_phi_psis(CMMDBManager *mol_in) {
   bool is_primary = 1;
   generate_phi_psis(mol_in, is_primary);
}

// fill phi_psi vector
//
void
coot::rama_plot::generate_secondary_phi_psis(CMMDBManager *mol_in) {
   bool is_primary = 0;
   generate_phi_psis(mol_in, is_primary);
}


// fill phi_psi vector
//
void
coot::rama_plot::generate_phi_psis(CMMDBManager *mol_in, bool is_primary) {

   int n_models = mol_in->GetNumberOfModels();
   phi_psi_model_sets.clear();
   // add a place-holder for the "0-th" model
   coot::phi_psis_for_model_t empty(0);
   phi_psi_model_sets.push_back(empty);
   for (int imod=1; imod<=n_models; imod++) {
      coot::phi_psis_for_model_t model_phi_psis(imod);
      CModel *model_p = mol_in->GetModel(imod);
      if (model_p) { 
	 CChain *chain_p;
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    CResidue *residue_p;
	    if (nres > 2) { 
	       for (int ires=1; ires<(nres-1); ires++) { 
		  residue_p = chain_p->GetResidue(ires);

		  // this could be improved
		  CResidue *res_prev = chain_p->GetResidue(ires-1);
		  CResidue *res_next = chain_p->GetResidue(ires+1);

		  if (res_prev && residue_p && res_next) {
		     try {
			// coot::phi_psi_t constructor can throw an error
			// (e.g. bonding atoms too far apart).
			coot::residue_spec_t spec(residue_p);
			coot::util::phi_psi_t pp(res_prev, residue_p, res_next);
			model_phi_psis.add_phi_psi(spec, pp);
		     }
		     catch (std::runtime_error rte) {
			// nothing too bad, just don't add that residue
			// to the plot
		     }
		  }
	       }
	    }
	 }
      }
      if (is_primary) 
	 phi_psi_model_sets.push_back(model_phi_psis);
      else 
	 secondary_phi_psi_model_sets.push_back(model_phi_psis);
   }
}

void
coot::rama_plot::generate_phi_psis_by_selection(CMMDBManager *mol,
						bool is_primary,
						int SelectionHandle) {

   if (is_primary) { 
      phi_psi_model_sets.clear();
      coot::phi_psis_for_model_t empty(0);
      phi_psi_model_sets.push_back(empty);
   } else {
      secondary_phi_psi_model_sets.clear();
      coot::phi_psis_for_model_t empty(0);
      secondary_phi_psi_model_sets.push_back(empty);
   }
   PCResidue *residues = NULL;
   int n_residues;
   mol->GetSelIndex(SelectionHandle, residues, n_residues);
   coot::phi_psis_for_model_t model_phi_psis(1); // model 1.
   for (int ires=1; ires<(n_residues-1); ires++) {
      CResidue *res_prev = residues[ires-1];
      CResidue *res_this = residues[ires];
      CResidue *res_next = residues[ires+1];
      std::string chain_id_1 = res_prev->GetChainID();
      std::string chain_id_2 = res_this->GetChainID();
      std::string chain_id_3 = res_next->GetChainID();
      if (chain_id_1 == chain_id_2) {
	 if (chain_id_2 == chain_id_3) {
	    try { 
	       coot::util::phi_psi_t pp(res_prev, res_this, res_next);
	       coot::residue_spec_t spec(res_this);
	       model_phi_psis.add_phi_psi(spec, pp);
	    }
	    catch (std::runtime_error rte) {
	       // nothing bad.
	    }
	 }
      }
   }
   if (is_primary) 
      phi_psi_model_sets.push_back(model_phi_psis);
   else 
      secondary_phi_psi_model_sets.push_back(model_phi_psis);
}

void
coot::rama_plot::generate_phi_psis_debug()
{
   // 
}



void
coot::rama_plot::map_mouse_pos(double x, double y) {


   mouse_util_t t = mouse_point_check(x,y); 

   // if we are close to (over) a residue box/triangle, then display a
   // tooltip_like_box, if not, then if we don't have sticky labels,
   // then delete the one that's was being displayed.
   // 
   if (!t.spec.unset_p()) {


      std::string residue_name;
      if (drawing_differences) {
	 if (! t.mouse_over_secondary_set) { 
	    residue_name = phi_psi_model_sets[t.model_number].phi_psi[t.spec].residue_name();
	    residue_type_background_as(residue_name);
	    tooltip_like_box(t);
	 } else {
	    residue_name = secondary_phi_psi_model_sets[t.model_number].phi_psi[t.spec].residue_name();
	    residue_type_background_as(residue_name);
	    tooltip_like_box(t);
	 } 
      } else {

	 // OK, so is this a pre-PRO? (i.e. is the next residue a PRO?
	 // (first, check that the residue exists before asking what
	 // it's residue type is - not that we can't use [], because
	 // that might add a map item (in that case that the index is
	 // not found in the map - and that should happed at the end
	 // of every fragment/chain)).
	 //
	 coot::residue_spec_t next_res_spec(t.spec.chain,
					    t.spec.resno+1,
					    t.spec.insertion_code);

	 std::map<residue_spec_t, util::phi_psi_t>::const_iterator it = 
	    phi_psi_model_sets[t.model_number].phi_psi.find(next_res_spec);

	 if (it != phi_psi_model_sets[t.model_number].phi_psi.end()) { 

	    std::string next_residue_name = it->second.residue_name();

	    if (next_residue_name == "PRO")
	       residue_name = "PRE-PRO";
	 }

	 residue_name = phi_psi_model_sets[t.model_number].phi_psi[t.spec].residue_name();
	 residue_type_background_as(residue_name);
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
      phi_psi_model_sets[1].phi_psi[0] = coot::util::phi_psi_t(worldx, -worldy, "dum",
							 "moving", 0, "", "");
      clear_last_canvas_item();
      draw_phi_psi_points();

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

      if (! t.spec.unset_p()) {
	 if (! t.mouse_over_secondary_set) 
	    std::cout << "   " << t.spec << "  "
		      << phi_psi_model_sets[t.model_number].phi_psi[t.spec] << std::endl;
      } 

      recentre_graphics_maybe(t);

   } else {

      if (event->button == 3) {
	 zoom += 0.2;
	 gtk_canvas_set_pixels_per_unit(canvas,zoom);
      }
   }
   return 1; // Handled this, right?
}

void
coot::rama_plot::recentre_graphics_maybe(mouse_util_t t) {

   if (t.model_number != coot::mouse_util_t::MODEL_NUMBER_UNSET) {
      if (t.mouse_over_secondary_set) { 
	 set_go_to_atom_molecule(molecule_numbers_.second);
      } else {
	 set_go_to_atom_molecule(molecule_numbers_.first);
      }
      set_go_to_atom_chain_residue_atom_name(t.spec.chain.c_str(), t.spec.resno, " CA ");
   }

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
      if (! t.spec.unset_p()) {
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

//    int n = phi_psi_sets[t.ichain].phi_psi.size(); // FIXME?
//    int ich = t.ichain;

   gtk_canvas_window_to_world(canvas, x,y, &worldx, &worldy);

   if (drawing_differences) {

      // need only to check over the top n_diffs differences.  (but we
      // do need to check the primary and the secondary phi_psi sets).
      //
      t = mouse_point_check_differences(worldx, worldy);
      
   } else {
      
      // single molecule, check over all ramachandran angles in the
      // first list:
      // 
      if (worldx <= 180 && worldx >= -180) { 
	 if (worldy <= 180 && worldy >= -180) {
	    for (unsigned int imod=1; imod<phi_psi_model_sets.size(); imod++) {
	       bool is_secondary = 0;
	       t = mouse_point_check_internal(phi_psi_model_sets[imod], imod, worldx, worldy, is_secondary);
	    }
	 }
      }
   }
   return t; 
}

coot::mouse_util_t
coot::rama_plot::mouse_point_check_internal(const coot::phi_psis_for_model_t &phi_psi_set,
					    int imod, 
					    double worldx, double worldy,
					    bool is_secondary) const {

   coot::mouse_util_t t;
   std::map<coot::residue_spec_t, coot::util::phi_psi_t>::const_iterator it;
   double diff1x, diff1y;
   double smallest_diff = 999999; 
   for (it=phi_psi_set.phi_psi.begin(); it!=phi_psi_set.phi_psi.end(); it++) {
      diff1x = fabs(it->second.phi() - worldx); 
      diff1y = fabs(it->second.psi() + worldy);
      if ((diff1x < 3) && (diff1y < 3)) {
	 if ((diff1x+diff1y) < smallest_diff) {
	    t.spec = it->first;
	    t.model_number = imod;
	 }
      }
   }
   t.mouse_over_secondary_set = is_secondary;
   return t;
} 


coot::mouse_util_t
coot::rama_plot::mouse_point_check_differences(double worldx, double worldy) const {

   coot::mouse_util_t t;
   
   for (unsigned int imod=1; imod<phi_psi_model_sets.size(); imod++) {
      bool is_secondary = 0;
      t = mouse_point_check_internal(phi_psi_model_sets[imod], imod, worldx, worldy, is_secondary);
   }

   if (t.spec.unset_p()) {
      for (unsigned int imod=1; imod<secondary_phi_psi_model_sets.size(); imod++) {
	 bool is_secondary = 1;
	 t = mouse_point_check_internal(secondary_phi_psi_model_sets[imod],
					imod, worldx, worldy, is_secondary);
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

	 // PRE-PRO fix up needed here.
	 
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
      draw_phi_psi_points();
   }

   draw_green_box();
}

void
coot::rama_plot::tooltip_like_box(const mouse_util_t &t) {

   if (! t.spec.unset_p()) { 
      bool is_secondary = t.mouse_over_secondary_set;

      // make a box slightly below the mouse position
      // (i.e. at highter psi)
      double phi =  phi_psi_model_sets[t.model_number][t.spec].phi();
      double psi = -phi_psi_model_sets[t.model_number][t.spec].psi();
      std::string label = phi_psi_model_sets[t.model_number][t.spec].label();

//       std::cout << "      debug:: tooltip_like_box for "
// 		<< t.spec << " phi: " << phi << " psi: " << psi
// 		<< label << std::endl;
	 
      if (is_secondary) {
	  phi   =  secondary_phi_psi_model_sets[t.model_number][t.spec].phi();
	  psi   = -secondary_phi_psi_model_sets[t.model_number][t.spec].psi();
	  label = secondary_phi_psi_model_sets[t.model_number][t.spec].label();
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
					      "text", label.c_str(),
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
					      "text", label.c_str(),
					      "x", phi-tw/2.0,
					      "y", psi+20.0/zoom,
					      "anchor",GTK_ANCHOR_WEST,
					      "font", fixed_font_str.c_str(),
					      "fill_color", "black",
					      NULL);
   }
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
   generate_phi_psis(mol);

   draw_phi_psi_points();
   
   delete mol; 
}


void
coot::rama_plot::draw_it(CMMDBManager *mol) {

   display_background();
   draw_axes();
   draw_zero_lines(); 
   generate_phi_psis(mol);
   coot::rama_stats_container_t counts = draw_phi_psi_points();
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
   coot::phi_psis_for_model_t phi_psi_set(1);
   coot::residue_spec_t spec("", 0, "");
   phi_psi_set.add_phi_psi(spec, phipsi);
   phi_psi_model_sets.push_back(phi_psi_set);
   draw_phi_psi_points();
}


void
coot::rama_plot::draw_it(const std::vector<coot::util::phi_psi_t> &phi_psi_s) {

   phi_psi_model_sets.clear();
   clear_last_canvas_items(phi_psi_s.size());
   
   coot::phi_psis_for_model_t phi_psi_set_m0(0);
   coot::phi_psis_for_model_t phi_psi_set(1);
   phi_psi_model_sets.push_back(phi_psi_set_m0); // dummy/unused for model 0
   for (unsigned int i=0; i<phi_psi_s.size(); i++) {
      coot::residue_spec_t spec("", i, "");
      phi_psi_set.add_phi_psi(spec, phi_psi_s[i]);
   }
   phi_psi_model_sets.push_back(phi_psi_set);
   draw_phi_psi_points();
   
}

void 
coot::rama_plot::draw_2_phi_psi_sets_on_canvas(CMMDBManager *mol1, 
					       CMMDBManager *mol2) { 

//    generate_phi_psis(&phi_psi_sets, mol1);
//    generate_phi_psis(&secondary_phi_psi_sets, mol2);


   // Are you sure that you want to edit this function? (not the one below?)

   bool primary = 1;
   generate_phi_psis(mol1,  primary);
   generate_phi_psis(mol2, !primary);
   
   // all_plot_background_and_bits(clipper::Ramachandran::All);
   // std::cout << "finding differences..." << std::endl;
   find_phi_psi_differences(); 
   // std::cout << "drawing differences..." << std::endl;
   draw_phi_psi_differences();

   // do we need to do something like this?
//    for (int ich=0; ich<phi_psi_sets.size(); ich++) 
//       draw_phi_psi_points(ich);

   drawing_differences = 1; // set flag for later drawing
}

// Kleywegt plots call this function
void 
coot::rama_plot::draw_2_phi_psi_sets_on_canvas(CMMDBManager *mol1, 
					       CMMDBManager *mol2,
					       std::string chainid1,
					       std::string chainid2) {

   int imod = 1;

   phi_psi_model_sets.clear();
   secondary_phi_psi_model_sets.clear();

   int selhnd_1 = mol1->NewSelection();
   int selhnd_2 = mol2->NewSelection();

   mol1->Select(selhnd_1, STYPE_RESIDUE, 0,
		chainid1.c_str(),
		ANY_RES, "*",
		ANY_RES, "*",
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"*",   // altLocs
		SKEY_NEW // selection key
		);

   mol2->Select(selhnd_2, STYPE_RESIDUE, 0,
		chainid2.c_str(),
		ANY_RES, "*",
		ANY_RES, "*",
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"*",   // altLocs
		SKEY_NEW // selection key
		);
   
   generate_phi_psis_by_selection(mol1, 1, selhnd_1);
   generate_phi_psis_by_selection(mol2, 0, selhnd_2);

   mol1->DeleteSelection(selhnd_1);
   mol2->DeleteSelection(selhnd_2);
   
   // all_plot_background_and_bits(clipper::Ramachandran::All);
   // std::cout << "finding differences..." << std::endl;
   find_phi_psi_differences(chainid1, chainid2); 
   // std::cout << "drawing differences..." << std::endl;
   draw_phi_psi_differences();

   drawing_differences = 1; // set flag for later drawing
}

void
coot::rama_plot::find_phi_psi_differences() {
   std::string chain_id1;
   std::string chain_id2;
   find_phi_psi_differences_internal(chain_id1, chain_id2, 0);
}

void
coot::rama_plot::find_phi_psi_differences(const std::string &chain_id1,
					  const std::string &chain_id2) { 
   find_phi_psi_differences_internal(chain_id1, chain_id2, 1);
}


void
coot::rama_plot::find_phi_psi_differences_internal(const std::string &chain_id1,
						   const std::string &chain_id2,
						   bool use_chain_ids) {

   double d1, d2;

   std::map<coot::residue_spec_t, coot::util::phi_psi_t>::const_iterator it;

   diff_sq.clear();
   for (unsigned int imod=1; imod<phi_psi_model_sets.size(); imod++) {
      std::cout << "in find_phi_psi_differences_internal with model number " << imod
		<< " primary size " << phi_psi_model_sets[imod].size() << std::endl;
      for (it=phi_psi_model_sets[imod].phi_psi.begin();
	   it!=phi_psi_model_sets[imod].phi_psi.end(); it++) {
	 coot::residue_spec_t spec_1 = it->first;
	 coot::residue_spec_t spec_2 = it->first;
	 // now (re)set chain_id of spec_2
	 // FIXME
	 if (use_chain_ids) {
	    spec_2.chain = chain_id2;
	 }
	 coot::util::phi_psi_t pp_1 = it->second;
	 coot::util::phi_psi_t pp_2 = secondary_phi_psi_model_sets[imod][spec_2];
	 
	 if (pp_2.is_filled()) {
	    d1 = fabs(pp_1.phi() - pp_2.phi());
	    d2 = fabs(pp_1.psi() - pp_2.psi());
	    if (d1 > 180) d1 -= 360; 
	    if (d2 > 180) d2 -= 360;
	    double v = sqrt(d1*d1 + d2*d2);
	    diff_sq_t ds(pp_1, pp_2, spec_1, spec_2, v); // something
	    diff_sq.push_back(ds);
	 } else {
	    std::cout << "Not found " << spec_2 << " from ref spec " << spec_1
		      << " in " << secondary_phi_psi_model_sets[imod].size()
		      << " residue in secondary set" << std::endl;
	 } 
      }
   }
      
   // sort diff_sq

   std::cout << " debug:: -- generated " << diff_sq.size() << " rama differences " << std::endl;
   std::sort(diff_sq.begin(), diff_sq.end(), compare_phi_psi_diffs);

}


bool coot::compare_phi_psi_diffs(const diff_sq_t &d1, const diff_sq_t d2) {
      return (d1.v() > d2.v());
}


// n_diffs is the max number of differences to be displayed, 50 typically.
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

   if (phi_psi_model_sets.size() == 0 ) {

      std::cout << "oops! No phi psi data!" << std::endl;

   } else {

      
      int n_vsize = diff_sq.size();

      // std::cout << "debug:: draw_phi_psi_differences() " << n_vsize << " " << n_diffs << std::endl;

      for (int j=0; j<n_diffs && j<n_vsize; j++) {

	 coot::util::phi_psi_t pp1 = diff_sq[j].phi_psi_1();
	 coot::util::phi_psi_t pp2 = diff_sq[j].phi_psi_2();
	 
	 draw_phi_psi_point(pp1, 1);
	 draw_phi_psi_point(pp2, 0);
	 GtkCanvasPoints *points = gtk_canvas_points_new(2);

	 draw_kleywegt_arrow(pp1, pp2, points);
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

      if (0)
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
   std::cout << "ramadebug: is phi_psi_sets.size " << phi_psi_model_sets.size() << std::endl;
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

