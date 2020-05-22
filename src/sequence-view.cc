/* src/sequence-view.cc
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006 The University of York
 * Author: Paul Emsley
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */
 
#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

#if defined _MSC_VER
#define snprintf _snprintf
#include <windows.h>
#endif

#include <string>
#include "sequence-view.hh"
#include <mmdb2/mmdb_tables.h>
#include "seq-view-interface.h"
#include "graphics-info.h" // for the callback



#ifdef SEQ_VIEW_STANDALONE
coot::sequence_view *coot::sequence_view_object_t::seq_view = NULL; 
#endif


#ifdef HAVE_GNOME_CANVAS
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
  #define gtk_canvas_window_to_world gnome_canvas_window_to_world
#endif

coot::sequence_view::sequence_view(mmdb::Manager *mol_in, std::string name, int coot_mol_no_in) {

   GtkWidget *top_lev = create_sequence_view_dialog();
   gtk_widget_set_size_request(GTK_WIDGET(top_lev), 500, 160);
   molecule_names.push_back(name);
   setup_internal(mol_in);
   mol.push_back(mol_in);
   coot_mol_no = coot_mol_no_in;
   gtk_object_set_user_data(GTK_OBJECT(canvas), (void *) this);
   // this user datum is used in the dialog destroy method (so that
   // the graphics_info_t sequence_view_is_displayed vector knows that
   // this object no longer is displayed)
   // 
   // (on_sequence_view_dialog_destroy).
   gtk_object_set_user_data(GTK_OBJECT(top_lev), GINT_TO_POINTER(coot_mol_no_in));

   // connect canvas (which was created in setup_internal) to top_lev:
   //
   GtkWidget *scrolled_window = seq_lookup_widget(GTK_WIDGET(top_lev),
					      "sequence_view_scrolledwindow");
   gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrolled_window),
					 GTK_WIDGET(canvas));
   gtk_widget_show(top_lev);

}

void
coot::sequence_view::setup_internal(mmdb::Manager *mol_in) {

//    GtkWidget *sequence_view_window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
//    gtk_container_set_border_width (GTK_CONTAINER(sequence_view_window), 2);
//    gtk_window_set_title (GTK_WINDOW(sequence_view_window), ("Sequence View"));

   res_offset = 30;
   res_scale = 6;
   row_offset = 10;
   row_scale = 14;
   fixed_font = "fixed";
#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
   fixed_font = "monospace";
   res_scale = 8;
#endif   
   GdkFont *font;
   font = gdk_font_load(fixed_font.c_str());
   gint res_width = gdk_string_width(font, "m");
   //std::cout <<"BL DEBUG:: font width calc "<<res_width <<" and set " << res_scale<< std::endl;
   res_scale = res_width + 2;

   tooltip_item = NULL;
   tooltip_item_text = NULL;

   int mnr = max_number_of_residues_in_a_chain(mol_in);
   mmdb::Model *model_p = mol_in->GetModel(1);
   int n_chains = model_p->GetNumberOfChains();
   
   setup_canvas(mnr, n_chains);
   // draw_debugging_box();
   mol_to_canvas(mol_in);


}

void
coot::sequence_view::draw_debugging_box() const { 

   GtkCanvasItem *item;

   item = gtk_canvas_item_new(gtk_canvas_root(canvas),
			      GTK_CANVAS_TYPE_CANVAS_RECT,
			      "x1", 0.0,
			      "y1", 0.0,
			      "x2", 100.0,
			      "y2",  50.0,
			      "fill_color", "grey80",
			      "outline_color", "black",
			      NULL);

   item = gtk_canvas_item_new(gtk_canvas_root(canvas),
			      GTK_CANVAS_TYPE_CANVAS_RECT,
			      "x1", -400.0,
			      "y1",    0.0,
			      "x2", -300.0,
			      "y2",   50.0,
			      "fill_color", "grey90",
			      "outline_color", "black",
			      NULL);

} 

void
coot::sequence_view::setup_canvas(int max_n_res, int n_chains) {

   // scrolling notes: I don't know the right values for the
   // set_scroll_region, but we have got a good x scroll with top
   // top_level widget being smaller usize than the canvas widget
   // (300 160) vs (800 60).
   //
   // Now we have the problem that when the canvas is started it's
   // putting (0,0) in canvas coordinates at the centre of the screen.
   // We would rather have (0,0) top left somewhere. That is fixed by
   // setting the first 2 parameters to gtk_canvas_set_scroll_region.
   // I don't know what the other parameters do.

#ifdef HAVE_GTK_CANVAS
   gtk_canvas_init(); 
#endif
   canvas = GTK_CANVAS(gtk_canvas_new());
//    std::cout << "DEBUG::  max_n_res = " << max_n_res 
// 	     << "   res_offset = " << res_offset << std::endl;

   int usize_x = 330 + int (res_offset) + int (res_scale*float(max_n_res));
   int usize_y = n_chains * 15 + 5;
   double scroll_width;   // right 
   double scroll_height;  // lower
   
   scroll_width = usize_x - 330; // This offset moves the left hand
				 // edge leftwards, so that the
				 // molecule label can be seen.
				 // Ideally, it would be related to
				 // the length of the longest chain
				 // label, but for now it isn't.
   scroll_height = usize_y;

   gtk_widget_set_size_request(GTK_WIDGET(canvas), usize_x, usize_y);
   gtk_widget_show(GTK_WIDGET(canvas));

   gtk_widget_set_events(GTK_WIDGET(canvas),
			 GDK_EXPOSURE_MASK      |
			 GDK_BUTTON_PRESS_MASK  |
			 GDK_BUTTON_RELEASE_MASK|
			 GDK_POINTER_MOTION_MASK|
			 GDK_KEY_RELEASE_MASK   |
			 GDK_POINTER_MOTION_HINT_MASK);


   // the scroll_height only affects at what widget height we get
   // recentreing - i.e. with a low scroll_height, the text is
   // recentered at small widget height, bigger scroll_height means
   // the text is recentred only when the widget is higher (longer).
   // These are not the 'droids you're looking for.

//    std::cout << "usize: " << usize_x
// 	     << " scroll_width: " << scroll_width << std::endl;
//  * @x1: Leftmost limit of the scrolling region.
//  * @y1: Upper limit of the scrolling region.
//  * @x2: Rightmost limit of the scrolling region.
//  * @y2: Lower limit of the scrolling region.

//    gtk_canvas_set_scroll_region(canvas, scroll_width, 18.0, 30.0, 40.0);
   double left_limit =  00.0;
   double upper_limit = 00.0;
//    std::cout << "DEBUG: canvas scroll_region: " 
// 	     << left_limit << " " << upper_limit << " " 
// 	     << scroll_width << " " << scroll_height 
// 	     << " (left upper right lower) "<< std::endl;
//    std::cout << "DEBUG: canvas size:          " 
// 	     << usize_x << " " << usize_y << std::endl;

   gtk_canvas_set_scroll_region(canvas, left_limit, upper_limit, scroll_width, scroll_height);

   gtk_signal_connect (GTK_OBJECT(canvas), "button_press_event",
		       GTK_SIGNAL_FUNC(seq_view_button_press), NULL);
   gtk_signal_connect (GTK_OBJECT(canvas), "motion_notify_event",
		       GTK_SIGNAL_FUNC(seq_view_motion_notify), NULL);

}

// size  scroll
//  800    600   size/1.4 + 
// 1400   1000

// static
// coot::sequence_view_res_info_t
// coot::sequence_view::get_res_info_from_event(GtkWidget *widget,
// 					     GdkEventButton *event) {

   
// }


// static
gint
coot::sequence_view::seq_view_button_press (GtkWidget *widget,
					    GdkEventButton *event) {
   
   // we need to convert x, y to sequence number space
   // Now, x and y are in widget space, 
   // 
   double x,y;
   // int x_as_int, y_as_int;
   // GdkModifierType state;
   x = event->x;
   y = event->y;

   // so now we need to get to a canvas to convert from widget space
   // to canvas world coordinates.
   double worldx, worldy;
   gtk_canvas_window_to_world(GTK_CANVAS(widget), x,y, &worldx, &worldy);
   coot::sequence_view *seq_view;

#ifdef SEQ_VIEW_STANDALONE
   seq_view = coot::sequence_view_object_t::seq_view;
#else
//    int *imol = (int *) gtk_object_get_user_data(GTK_OBJECT(widget));
//    graphics_info_t g;
//    seq_view = g.get_sequence_view(*imol);
   seq_view = (coot::sequence_view *) gtk_object_get_user_data(GTK_OBJECT(widget));
#endif
   
   if (seq_view) {
       coot::sequence_view_res_info_t si =
 	 seq_view->get_sequence_view_res_info(worldx, worldy);
       if (si.residue_serial_number > -1 ) {

#ifndef SEQ_VIEW_STANDALONE

	  graphics_info_t g;
	  g.set_go_to_atom_molecule(seq_view->Coot_Mol_No());
	  si = seq_view->chain_and_resno(si);
	  if (si.residue_serial_number > -1) {
// 	     g.set_go_to_atom_chain_residue_atom_name(si.chain.c_str(),
// 						      si.residue_seq_num,
// 						      " CA ");
	     g.set_go_to_residue_intelligent(si.chain, si.residue_seq_num,
					     si.ins_code);
	     int success = g.try_centre_from_new_go_to_atom(); 
	     if (success)
		g.update_things_on_move_and_redraw(); 

	  } else {
	     std::cout << "ERROR: Oops can't find residue" << std::endl;
	  }
#endif
      } else {
	 // clear the tooltip box then
	 std::cout << "missed\n";
	 seq_view->clear_tooltip_box();
      }
   } else {
      std::cout << "No seq view" << std::endl;
   }
   return 0;
}

coot::sequence_view_res_info_t
coot::sequence_view::chain_and_resno(const coot::sequence_view_res_info_t &in) const {

   coot::sequence_view_res_info_t out = in;
   // out.residue_serial_number = -1; // signal error

   if (in.residue_serial_number >= 0) {

      // Problem here?  We presume that mol is valid.  What happens if
      // it was closed behind our back?
      //
      // We must make sure that we don't come here then - close the widget.
      // 
      mmdb::Model *model_p = mol[in.molecule_number]->GetModel(1);
      mmdb::Chain *chain_p = model_p->GetChain(in.chain_number);

      if (! chain_p) {
	 std::cout << "ERROR:: missing (NULL) chain! " << std::endl;
      } else { 
	 mmdb::Residue *residue_p = chain_p->GetResidue(in.residue_serial_number);
	 if (!residue_p) {
	    out.residue_serial_number = -1; // signal an error in finding residue
	 } else {
	    //       std::cout << "Real chain " << chain_p->GetChainID()
	    // 		<< " real resno " << residue_p->GetResName()
	    // 		<< " " << residue_p->GetSeqNum() << std::endl;
	    out.chain = chain_p->GetChainID();
	    out.residue_seq_num = residue_p->GetSeqNum();
	    out.resname = residue_p->GetResName();
	    out.ins_code = residue_p->GetInsCode();
	 }

      }
   } else {
      out.residue_serial_number = -1; // signal an error
   }
   return out;
}


// static
gint
coot::sequence_view::seq_view_motion_notify(GtkWidget *widget, GdkEventMotion *event) {

   // std::cout << "mouse motion" << std::endl;
   // we need to convert x, y to sequence number space
   // Now, x and y are in widget space, 
   //
   // 
   double x,y;
   int x_as_int, y_as_int;
   GdkModifierType state;

   // This is very important if we want events to keep happening as we
   // move the mouse
   // 
   if (event->is_hint) {
      // gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
      GdkModifierType mask;
      GdkSeat *seat = gdk_display_get_default_seat(gdk_display_get_default());
      GdkDevice *mouse = gdk_seat_get_pointer(seat);
      gdk_window_get_device_position(event->window, mouse, &x_as_int, &y_as_int, &mask);
   }
   x = event->x;
   y = event->y;

   // so now we need to get to a canvas to convert from widget space
   // to canvas world coordinates.
   double worldx, worldy;
   gtk_canvas_window_to_world(GTK_CANVAS(widget), x,y, &worldx, &worldy);

   coot::sequence_view *seq_view;

#ifdef SEQ_VIEW_STANDALONE
   seq_view = coot::sequence_view_object_t::seq_view;
#else
   // widget is the canvas
//    int *imol = (int *) gtk_object_get_user_data(GTK_OBJECT(widget));
//    graphics_info_t g;
//    seq_view = g.get_sequence_view(*imol);
   seq_view = (coot::sequence_view *) gtk_object_get_user_data(GTK_OBJECT(widget));
#endif
   
   if (seq_view) {
       coot::sequence_view_res_info_t si =
 	 seq_view->get_sequence_view_res_info(worldx, worldy);
      if (si.residue_serial_number > -1 ) {
	 seq_view->tooltip_like_box(si);
      } else {
	 // clear the tooltip box then
	 seq_view->clear_tooltip_box();
      }
   } else {
      std::cout << "No seq view" << std::endl;
   }
   return 0;
}

coot::sequence_view_res_info_t
coot::sequence_view::get_sequence_view_res_info(double worldx, double worldy) const {

   coot::sequence_view_res_info_t si;
   int iserial = int ((worldx - res_offset)/res_scale);
   si.residue_serial_number = iserial;

   // and now the chain and molecule info:
   int row_no = int ((worldy - row_offset + 8)/row_scale);

   if ((row_no < int(sequence_row.size())) && (row_no >= 0)) {
      std::pair<int, int> p = sequence_row[row_no];
      si.row = row_no;
      si.molecule_number = p.first;
      si.chain_number = p.second;
      si.ins_code = "";
   } else {
      // std::cout << "no row\n";
      si.residue_serial_number = -1; // an impossible serial number
   } 
   return si;
}

void 
coot::sequence_view::tooltip_like_box(const coot::sequence_view_res_info_t &si) {

   // So what is the real chain id and seqnum?

   // It is quite resonable for si to have have a
   // residue_serial_number that is off the end of the chain, because
   // we have not yet made a check for the real residue in the
   // molecule, we only are passing the index, chain_and_resno() does
   // the real checking.
   // 
   coot::sequence_view_res_info_t si_new = chain_and_resno(si);
   
   if (si_new.residue_serial_number < 0) {
      return;
   }

   std::string label = molecule_names[si.molecule_number];
   label += " ";
   label += si_new.chain;
   label += ": ";
   label += seq_int_to_string(si_new.residue_seq_num);
   label += si_new.ins_code;
   label += " ";
   label += si_new.resname;

   float x1 = float (res_offset + res_scale*si.residue_serial_number) - 60;
   float y1 = float (row_offset + si.row*row_scale) + 10 ;

#ifdef WINDOWS_MINGW
// BL says: again we want to make it bigger to fit
  float tw = label.size() * 8.0 + 10.0;
#else
  float tw = label.size() * 6.0 + 10.0;
#endif // MINGW

   clear_tooltip_box();

   // std::cout << "square at " << x1 << " " << y1 << std::endl;

   tooltip_item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				      GTK_CANVAS_TYPE_CANVAS_RECT,
				      "x1", x1,
				      "y1", y1,
				      "x2", x1 + tw,
				      "y2", y1 + 16,
				      "fill_color", "PaleGreen",
				      "outline_color", "black",
				      NULL);

   tooltip_item_text = gtk_canvas_item_new(gtk_canvas_root(canvas),
					   GTK_CANVAS_TYPE_CANVAS_TEXT,
					   "text", label.c_str(),
					   "x", x1 + 5,
					   "y", y1 + 7,
					   "anchor",GTK_ANCHOR_WEST,
					   "font", fixed_font.c_str(),
					   "fill_color", "black",
					   NULL);
}

void
coot::sequence_view::clear_tooltip_box() {
   if (tooltip_item)
      gtk_object_destroy(GTK_OBJECT(tooltip_item));
   if (tooltip_item_text)
      gtk_object_destroy(GTK_OBJECT(tooltip_item_text));
   tooltip_item = NULL;
   tooltip_item_text = NULL;
} 

std::string
coot::sequence_view::seq_int_to_string(int i) const {
   char s[100];
   snprintf(s,99,"%d",i);
   return std::string(s);
}


// The canvas is all correct by the the time this has been called.
//
void
coot::sequence_view::mol_to_canvas(mmdb::Manager *mol_in) {
   
   // Insertion codes are ignored.

   GtkCanvasItem *item;
   std::string res_code;
   mmdb::Model *model_p = mol_in->GetModel(1);

   std::cout << "calculating secondary structure...";

   int status = model_p->CalcSecStructure(1); // Hmm. Used to have an atomselhnd arg.
   std::cout << "done.\n";

   if (status == mmdb::SSERC_Ok) {
      std::cout << "INFO:: SSE status was OK\n";
   } else {
      std::cout << "INFO:: SSE status was bad\n" << status << "\n";
   }
   

   coot::util::print_secondary_structure_info(model_p);

   // Sometimes we get a chain with no residues in it.  That means
   // that we need to keep a count of actual rows, not chains (we act
   // as if the chain with no residues was not there).  The first row
   // should be positions at row 0, so there is an offset of -1 in
   // canvas item position because for the first one we have already
   // done a row++ by the time we get there.
   // 
   int row = 0;

   std::string colour;
   mmdb::Chain *chain_p;
   int n_chains = model_p->GetNumberOfChains(); 
   for (int i_chain=0; i_chain<n_chains; i_chain++) {
      chain_p = model_p->GetChain(i_chain); 
      std::string mol_chain(chain_p->GetChainID());

//       std::cout << " mol_to_canvas molecule_names has size "
// 		<< molecule_names.size() << std::endl;

      int nres = chain_p->GetNumberOfResidues();
      // std::cout << "chain " << i_chain << " has " << nres << " residues:\n";
      if (nres > 0) {
	 draw_mol_chain_label(mol_chain, row, 0);
	 std::pair<int, int> p(mol.size(), i_chain);
	 sequence_row.push_back(p);
	 row++;
      } 
      mmdb::Residue *residue_p;
      for (int ires=0; ires<nres; ires++) { // ires is a serial number
	 residue_p = chain_p->GetResidue(ires); 
// 	 std::cout << "DEBUG:: GetResName: " << residue_p->GetResName() 
// 		   << " " << strlen(residue_p->GetResName()) << std::endl;
	 res_code = coot::util::three_letter_to_one_letter_with_specials(residue_p->GetResName());
	 // std::cout << res_code;

	 colour = colour_by_secstr(residue_p, model_p);

	 item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				    GTK_CANVAS_TYPE_CANVAS_TEXT,
				    "text", res_code.c_str(), 
				    "x", res_offset + ires*res_scale,
				    "y", row_offset + (row - 1)*row_scale,
				    "anchor",GTK_ANCHOR_WEST,
				    "font", fixed_font.c_str(),
				    "fill_color", colour.c_str(),
				    NULL);
	 canvas_item_vec.push_back(item);
      }
      // std::cout << "\n";
   }
}

std::string
coot::sequence_view::colour_by_secstr(mmdb::Residue *residue_p, mmdb::Model *model_p) const {

   std::string s("black");

   switch (residue_p->SSE)  {

   case mmdb::SSE_Strand : s = "firebrick3";  break;
   case mmdb::SSE_Bulge  : s = "firebrick1";  break;
   case mmdb::SSE_3Turn  : s = "MediumBlue";  break;
   case mmdb::SSE_4Turn  : s = "SteelBlue4";  break;
   case mmdb::SSE_5Turn  : s = "DodgerBlue4"; break;
   case mmdb::SSE_Helix  : s = "navy";        break;
   case mmdb::SSE_None   : s = "black";       break;
   } 

   return s;
}


coot::sequence_view::sequence_view(mmdb::Manager *mol_in, GtkWidget *container_widget_in) {

   // we need to create a canvas in container_widget
   //
   setup_internal(mol_in); 
}

void
coot::sequence_view::draw_mol_chain_label(std::string mol_chain, int i_row, int mol_no) {

   GtkCanvasItem *item;
   // std::cout << "molecule_names has size " << molecule_names.size() << std::endl;
   std::string label = molecule_names[mol_no];
   label += " ";
   if (mol_chain == "")
      label += " ";
   else
      label += mol_chain;
   label += ":";
   
   item = gtk_canvas_item_new(gtk_canvas_root(canvas),
			      GTK_CANVAS_TYPE_CANVAS_TEXT,
			      "text", label.c_str(), 
			      "x", -res_scale*(label.length() + 1),
			      "y", row_offset + i_row*row_scale,
			      "anchor",GTK_ANCHOR_WEST,
			      "font", fixed_font.c_str(),
			      "fill_color", "black",
			      NULL);
   canvas_item_vec.push_back(item);


}


void 
coot::sequence_view::generate_from(mmdb::Manager *mol_in) {
   mol.push_back(mol_in);
   mol_to_canvas(mol_in);
}


void 
coot::sequence_view::regnerate() { // use mol pointer (it has had its atoms changed)
   mol_to_canvas(mol[0]);  // Fixme
}



void
coot::sequence_view::undisplay(int coot_molecule_no) {

   GtkCanvasItem *item;
   for (unsigned int i=0; i<canvas_item_vec.size(); i++) { 
      item = canvas_item_vec[i];
      if (canvas_item_vec[i] == NULL) {
	 std::cout << "oops - null canvas item" << std::endl;
      }
      gtk_object_destroy(GTK_OBJECT(item));
   }
   canvas_item_vec.clear();
} 


int
coot::sequence_view::max_number_of_residues_in_a_chain(mmdb::Manager *mol_in) const {

   int r = 0;
   int nres;
   mmdb::Model *model_p = mol_in->GetModel(1);
   if (model_p) {
      mmdb::Chain *chain_p;
      int n_chains = model_p->GetNumberOfChains(); 
      for (int i_chain=0; i_chain<n_chains; i_chain++) {
	 chain_p = model_p->GetChain(i_chain);
	 nres = chain_p->GetNumberOfResidues();
	 if (nres > r) {
	    r = nres;
	 }
      }
   }
   return r;
}

 
#endif //  HAVE_GTK_CANVAS
