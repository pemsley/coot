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

#ifndef HAVE_SEQUENCE_VIEW
#define HAVE_SEQUENCE_VIEW

#include <gtk/gtk.h>
#include <goocanvas.h>

#include <mmdb2/mmdb_manager.h>
#include <string>
#include <vector>

namespace coot {

   class sequence_view_res_info_t {
   public:
      int residue_serial_number;
      int molecule_number;
      int chain_number;
      int row;
      std::string chain;
      std::string resname;
      int residue_seq_num;
      std::string ins_code;
   };

   class sequence_view {

      std::vector<std::string> molecule_names;
      // a place to convert from canvas row to molecule number and
      // chain number thereof.
      std::vector<std::pair<int, int> > sequence_row;
      GtkWidget *container_widget;
      GooCanvas *canvas;
      int coot_mol_no; // for the coot callback
      float res_offset;
      float res_scale;
      float row_scale;
      float row_offset;
      std::vector<mmdb::Manager *> mol;
      void setup_internal(mmdb::Manager *mol);
      void mol_to_canvas(mmdb::Manager *mol); 
      std::string three_letter_to_one_letter(const std::string &resname) const;
      void setup_canvas(int max_n_res, int n_chains);
      void draw_mol_chain_label(std::string mol_chain, int i_chain, int mol_no);
      void tooltip_like_box(const coot::sequence_view_res_info_t &si);
      void clear_tooltip_box();
      GooCanvasItem *tooltip_item;
      GooCanvasItem *tooltip_item_text;
      std::string seq_int_to_string(int i) const;
      std::vector<GooCanvasItem *> canvas_item_vec;
      int max_number_of_residues_in_a_chain(mmdb::Manager *mol_in) const;
      std::string colour_by_secstr(mmdb::Residue *res, mmdb::Model *model) const;
      std::string fixed_font;
      void draw_debugging_box() const;


   public:
      sequence_view(mmdb::Manager *mol_in,
		    std::string mol_name,
		    int coot_mol_no_in);  // create a toplevel
      //! create the sequence-view inside container_widget:
      sequence_view(mmdb::Manager *mol_in, GtkWidget *container_widget);
      void generate_from(mmdb::Manager *mol_in);
      void regnerate(); // use mol pointer (it has had its atoms changed)
      sequence_view_res_info_t
      get_sequence_view_res_info(double worldx, double worldy) const;

      // button callbacks
      static gint seq_view_button_press (GtkWidget *widget, GdkEventButton *event);
      static gint seq_view_motion_notify(GtkWidget *widget, GdkEventMotion *event);
      int Coot_Mol_No() const { return coot_mol_no; }
      sequence_view_res_info_t
      chain_and_resno(const coot::sequence_view_res_info_t &in) const;

      // When the coordinates are closed, we want to undisplay all
      // chains that are associated with that molecule, so that we
      // don't try to look up the chain and residue information in a
      // molecule that had been deleted.
      // 
      void undisplay(int coot_molecule_number);
      GtkWidget *Canvas() { return GTK_WIDGET(canvas); }

   };

#ifdef SEQ_VIEW_STANDALONE   
   // This should be in graphics_info_t at some future stage
   class sequence_view_object_t {
   public:
      static sequence_view *seq_view;
   };
#endif
   
}

#endif //  HAVE_SEQUENCE_VIEW
