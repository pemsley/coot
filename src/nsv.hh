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

#if defined(HAVE_GNOME_CANVAS)

#ifndef HAVE_STRING
#include <string>
#define HAVE_STRING
#endif

#include <mmdb2/mmdb_manager.h>

#ifdef HAVE_GNOME_CANVAS
#include <gtk/gtk.h>
#include <libgnomecanvas/libgnomecanvas.h>
#endif
#ifdef HAVE_GOOCANVAS
#include <goocanvas.h>
#endif

#include "coot-utils/coot-coord-utils.hh"

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
	 GnomeCanvasItem *obj;
	 coot::atom_spec_t atom_spec;
         coot::residue_spec_t residue_spec;
	 int position_number;
  	 spec_and_object(int molecule_number_in, coot::atom_spec_t &spec_in,
			 int position_number_in) { 
	    mol_no = molecule_number_in;
	    atom_spec = spec_in;
            residue_spec = coot::residue_spec_t(atom_spec);
	    position_number = position_number_in;
	 } 
	 void add_rect_attribs(GnomeCanvasItem *rect_item_in) { 
	    obj = rect_item_in;
	 }
      };


      int molecule_number;
      GtkWidget *scrolled_window; // we use this for regenerate().
      std::vector<GnomeCanvasItem *> canvas_item_vec;
#ifdef HAVE_GOOCANVAS
      GooCanvasItem *canvas_group;
      GtkWidget *canvas;
      std::map<mmdb::Residue *, GooCanvasItem *> rect_residue_map;
      mmdb::Residue *current_highlight_residue;
#else
      GnomeCanvas *canvas;
#endif
      
      // return the canvas y size
      int setup_canvas(mmdb::Manager *mol);
      std::vector<chain_length_residue_units_t> get_residue_counts(mmdb::Manager *mol) const;
      bool use_graphics_interface_flag;
      static void on_nsv_close_button_clicked (GtkButton *button,
					       gpointer         user_data);
      static void on_nsv_dialog_destroy (GtkObject *obj,
					 gpointer user_data);
      static gint letter_event (GtkObject *obj, GdkEvent *event, gpointer data);
#ifdef HAVE_GOOCANVAS
      static gboolean rect_notify_event (GooCanvasItem *item,
                                         GooCanvasItem *target,
                                         GdkEvent *event,
                                         gpointer data);
      static gboolean rect_button_event (GooCanvasItem *item,
                                        GooCanvasItem *target,
                                        GdkEventButton *event,
                                        gpointer data);
#else
      static gint rect_event   (GtkObject *obj, GdkEvent *event, gpointer data);
#endif
      void draw_axes(std::vector<chain_length_residue_units_t>, int l, int b, double x_offset);
      std::string fixed_font_str;
      int pixels_per_letter;
      void mol_to_canvas(mmdb::Manager *mol, int lowest_resno, double x_offset);
      void chain_to_canvas(mmdb::Chain *chain_p, int position_number, int lowest_resno, double x_offset);
      void origin_marker();
      int tick_start_number(int l) const;
      int pixels_per_chain;
      bool add_text_and_rect(mmdb::Residue *residue_p, int pos_number, int lowest_resno, double x_offset);
      std::string colour_by_secstr(mmdb::Residue *residue_p) const;
      int points_max;
      void clear_canvas();
      void helix(mmdb::Chain *chain_p, int resno_low, int resno_high, double x_offet);
      void helix_single_inner(int i_turn_number, double x_start, double y_start, double hexix_scale);
      void strand(mmdb::Chain *, int resno_low, int resno_high, double x_offset, double y_offset, double scale);

   public:
      nsv(mmdb::Manager *mol,
	  const std::string &molecule_name,
	  int molecule_number_in,
          GtkWidget *main_window_vbox,
	  bool use_graphics_interface);
      nsv(mmdb::Manager *mol,
	  const std::string &molecule_name,
	  int molecule_number_in,
          GtkWidget *main_window_vbox,
	  bool use_graphics_interface,
	  int canvas_pixel_limit);
      void regenerate(mmdb::Manager *mol);
      GtkWidget *Canvas() const { return GTK_WIDGET(canvas); }
      // default is 22500 
      void set_points_max(int v) { points_max = v; }
      std::string sequence_letter_background_colour;

      // is this used?
      static gboolean on_canvas_button_press(GtkWidget      *canvas,
                                             GdkEventButton *event,
                                             gpointer        data);
      // used for the docked sequence view.
      static gint close_docked_sequence_view(GtkWidget *menu_item, GdkEventButton *event);

      // JED feature request 
      void highlight_residue(mmdb::Residue *residue_p);
   };
}


#endif // a CANVAS
