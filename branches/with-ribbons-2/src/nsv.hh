/* src/nsv.hh
 * 
 * Copyright 2009 by The University of Oxford
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

#ifndef HAVE_STRING
#include <string>
#define HAVE_STRING
#endif

#include <mmdb/mmdb_manager.h>

#ifdef HAVE_GTK_CANVAS

#include <gtk/gtk.h>
#include <gdk_imlib.h>
#include <gtk-canvas.h>

#else 

   #ifdef HAVE_GNOME_CANVAS
      #include <gtk/gtk.h>
      #include <libgnomecanvas/libgnomecanvas.h>
      typedef GnomeCanvas     GtkCanvas; 
      typedef GnomeCanvasItem GtkCanvasItem; 
   #endif
#endif

#include "coot-coord-utils.hh"

namespace exptl {

   class nsv {

      class chain_length_residue_units_t {
      public:
	 std::string chain_id;
	 int n_residue_count;
	 int first_res_no;
	 int max_resno;
	 chain_length_residue_units_t(std::string chain_id_in,
				      int n_res_count_in,
				      int first_in,
				      int max_resno_in) {
	    chain_id = chain_id_in;
	    n_residue_count = n_res_count_in;
	    first_res_no = first_in;
	    max_resno = max_resno_in;
	 }
      };

      class spec_and_object {
      public:
	 int mol_no;
	 GtkCanvasItem *obj;
	 coot::atom_spec_t atom_spec;
	 int position_number;
  	 spec_and_object(int molecule_number_in, coot::atom_spec_t &spec_in,
			 int position_number_in) { 
	    mol_no = molecule_number_in;
	    atom_spec = spec_in;
	    position_number = position_number_in;
	 } 
	 void add_rect_attribs(GtkCanvasItem *rect_item_in) { 
	    obj = rect_item_in;
	 }
      };


     

      int molecule_number;
      GtkCanvas *canvas;
      std::vector<GtkCanvasItem *> canvas_item_vec;
      void setup_canvas(CMMDBManager *mol, GtkWidget *scrolled_window);
      std::vector<chain_length_residue_units_t> get_residue_counts(CMMDBManager *mol) const;
      bool use_graphics_interface_flag;
      static void on_nsv_close_button_clicked (GtkButton *button,
					       gpointer         user_data);
      static void on_nsv_dialog_destroy (GtkObject *obj,
					 gpointer user_data);
      static gint letter_event (GtkObject *obj, GdkEvent *event, gpointer data);
      static gint rect_event   (GtkObject *obj, GdkEvent *event, gpointer data);
      void draw_axes(std::vector<chain_length_residue_units_t>, int l, int b, double x_offset);
      std::string fixed_font_str;
      int pixels_per_letter;
      void mol_to_canvas(CMMDBManager *mol, int lowest_resno, double x_offset);
      void chain_to_canvas(CChain *chain_p, int position_number, int lowest_resno, double x_offset);
      void origin_marker();
      int tick_start_number(int l) const;
      int pixels_per_chain;
      bool add_text_and_rect(CResidue *residue_p, int pos_number, int lowest_resno, double x_offset);
      std::string colour_by_secstr(CResidue *residue_p) const;
      
   public:
      nsv(CMMDBManager *mol,
	  const std::string &molecule_name,
	  int molecule_number_in,
	  bool use_graphics_interface);
      void regenerate(CMMDBManager *mol);
      GtkWidget *Canvas() const { return GTK_WIDGET(canvas); }
   };

}


#endif // a CANVAS
