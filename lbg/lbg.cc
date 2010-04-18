/* lbg/lbg.cc
 * 
 * Author: Paul Emsley
 * Copyright 2010 by The University of Oxford
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

#include <sys/types.h>  // for stating
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <algorithm>


#include <cairo.h>
#if CAIRO_HAS_PDF_SURFACE
#include <cairo-pdf.h>
#endif
#include "lbg.hh"

GtkWidget *get_canvas_from_scrolled_win(GtkWidget *canvas) {

   return canvas;
}

void
lbg_info_t::untoggle_others_except(GtkToggleToolButton *button_toggled_on) {
   
   std::map<std::string, GtkToggleToolButton *>::const_iterator it;
   for (it=widget_names.begin(); it!=widget_names.end(); it++) {
      if (it->second != button_toggled_on) {
	 if (gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(it->second))) {
	    gtk_toggle_tool_button_set_active(GTK_TOGGLE_TOOL_BUTTON(it->second), FALSE);
	 }
      }
   }
}

/* This handles button presses on the canvas */
static gboolean
on_canvas_button_press (GtkWidget *widget, GdkEventButton *event)
{
   int x_as_int, y_as_int;
   GdkModifierType state;
   gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
   GtkObject *obj = GTK_OBJECT(widget);
   if (obj) {
      gpointer gp = gtk_object_get_user_data(obj);
      if (gp) {
	 lbg_info_t *l = static_cast<lbg_info_t *> (gp);
	 l->set_mouse_pos_at_click(x_as_int, y_as_int); // save for dragging
	 // std::cout << "mouse_at_click: " << x_as_int << " " << y_as_int << std::endl;
	 
	 if (0) 
	    std::cout << "   on click: scale_correction  " << l->mol.scale_correction.first << " "
		      << l->mol.scale_correction.second
		      << " centre_correction " << l->mol.centre_correction << std::endl;

	    
	 if (l->in_delete_mode_p()) { 
	    l->handle_item_delete(event);
	 } else {
	    l->handle_item_add(event);
	 }
      }
   }
   return TRUE;
}

static gboolean
on_canvas_motion (GtkWidget *widget, GdkEventMotion *event) {
   int x_as_int=0, y_as_int=0;
   if (event->is_hint) {
      GdkModifierType state;
      gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
      GtkObject *obj = GTK_OBJECT(widget);
      if (obj) {
	 gpointer gp = gtk_object_get_user_data(obj);
	 // std::cout << "got gp: " << gp << std::endl;
	 if (gp) { 
	    lbg_info_t *l = static_cast<lbg_info_t *> (gp);
	    bool highlight_status = l->item_highlight_maybe(x_as_int, y_as_int);
	    if (! highlight_status) {
	       if (state & GDK_BUTTON1_MASK) {
// 		  std::cout << "dragging canvas given mouse point "
// 			    << x_as_int << " " << y_as_int <grep< std::endl;
		  l->drag_canvas(x_as_int, y_as_int);
	       }
	    }
	 } else {
	    std::cout << "Failed to get gpointer from object on canvas motion"
		      << std::endl;
	 } 
      } else {
	 std::cout << "Failed to get GtkObject from widget on canvas motion"
		   << std::endl;
      }
   }
   return TRUE;
}

void
lbg_info_t::drag_canvas(int mouse_x, int mouse_y) {

   double delta_x = (double(mouse_x) - mouse_at_click.x) * 1;
   double delta_y = (double(mouse_y) - mouse_at_click.y) * 1;

   mouse_at_click.x = mouse_x;
   mouse_at_click.y = mouse_y;

   // GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   // double top_left_left = GOO_CANVAS(canvas)->hadjustment->lower - delta_x;
   // double top_left_top  = GOO_CANVAS(canvas)->vadjustment->lower - delta_y;
   // std::cout << "scrolling to " << top_left_left << " " << top_left_top << std::endl;
   // goo_canvas_scroll_to(GOO_CANVAS(canvas), top_left_left, top_left_top);

   lig_build::pos_t delta(delta_x, delta_y);
   canvas_drag_offset += delta;

   // should santize delta *here* (and check for not 0,0) before
   // passing to these functions.

   if (is_sane_drag(delta)) { 
   
      clear_canvas();
      translate_residue_circles(delta);
      widgeted_molecule_t new_mol = translate_molecule(delta); // and do a canvas update

      render_from_molecule(new_mol);
      draw_bonds_to_ligand();
      draw_residue_circles(residue_circles);
   }
}


widgeted_molecule_t
lbg_info_t::translate_molecule(const lig_build::pos_t &delta) {  // and do a canvas update

   // we can't translate mol, because that gets wiped in render_from_molecule.
   widgeted_molecule_t new_mol = mol;
   new_mol.translate(delta);
   return new_mol;
}

void
lbg_info_t::translate_residue_circles(const lig_build::pos_t &delta) {
   
   for (unsigned int i=0; i<residue_circles.size(); i++) 
      residue_circles[i].pos += delta;

}

bool
lbg_info_t::is_sane_drag(const lig_build::pos_t &delta) const {

   bool sane = 0;
   if (delta.length() > 0.5) { // don't move with a delta of 0,0.
      if (delta.length() < 50 ) { // don't move with an absurd delta.
	 sane = 1;
      }
   }
   return sane;
} 



bool
lbg_info_t::item_highlight_maybe(int x, int y) {

   bool highlight_status = 0;
   remove_bond_and_atom_highlighting();
   std::pair<bool, widgeted_bond_t> bond_info = mol.highlighted_bond_p(x, y);
   if (bond_info.first) {
      highlight_status = 1;
      highlight_bond(bond_info.second, in_delete_mode_p());
   } else {
      std::pair<int, widgeted_atom_t> atom_info = mol.highlighted_atom_p(x, y);
      if (atom_info.first != UNASSIGNED_INDEX) {
	 highlight_status = 1;
	 highlight_atom(atom_info.second, atom_info.first, in_delete_mode_p());
      } 
   }
   return highlight_status;
}

void
lbg_info_t::remove_bond_and_atom_highlighting() {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   if (highlight_data.has_contents())
      highlight_data.clear(root);
}

// set highlight_data
void
lbg_info_t::highlight_bond(const lig_build::bond_t &bond, bool delete_mode) {

   int i1 = bond.get_atom_1_index();
   int i2 = bond.get_atom_2_index();
   std::pair<int, int> bond_indices(i1, i2);
   lig_build::pos_t A = mol.atoms[i1].atom_position;
   lig_build::pos_t B = mol.atoms[i2].atom_position;

   std::string col = "#20cc20";
   if (delete_mode)
      col = "#D03030";
   
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));
   GooCanvasItem *h_line =
      goo_canvas_polyline_new_line(root,
				   A.x, A.y,
				   B.x, B.y,
				   "line-width", 7.0,
				   "stroke-color", col.c_str(),
				   NULL);
   highlight_data = highlight_data_t(h_line, bond_indices, A, B);
}

// set highlight_data
void
lbg_info_t::highlight_atom(const lig_build::atom_t &atom, int atom_index, bool delete_mode) {

   lig_build::pos_t A = atom.atom_position;
   std::string col = "#20cc20";
   if (delete_mode)
      col = "#D03030";
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));
   double width = 10; // atom.name.length() * 3;
   double height = 14;
   double x1 = A.x - width/2;
   double y1 = A.y - height/2;

   // std::cout << x1 << " " << x2 << " "<< width << " " << height << std::endl;

   GooCanvasItem *rect_item =
      rect_item = goo_canvas_rect_new (root, x1, y1, width, height,
				       "line-width", 1.0,
				       "stroke-color", col.c_str(),
				       NULL);
   highlight_data = highlight_data_t(rect_item, A, atom_index);
   // std::cout << "Highlight atom with index " << atom_index << std::endl;

}


void
lbg_info_t::handle_item_add(GdkEventButton *event) {

   bool changed_status = 0;

   int x_as_int, y_as_int;
   GdkModifierType state;
   bool shift_is_pressed = 0;
   if (event->state & GDK_SHIFT_MASK)
      shift_is_pressed = 1;

   gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
   if (canvas_addition_mode == lbg_info_t::PENTAGON)
      changed_status = try_stamp_polygon(5, x_as_int, y_as_int, shift_is_pressed, 0);
   if (canvas_addition_mode == lbg_info_t::HEXAGON)
      changed_status = try_stamp_polygon(6, x_as_int, y_as_int, shift_is_pressed, 0);
   if (canvas_addition_mode == lbg_info_t::HEXAGON_AROMATIC)
      changed_status = try_stamp_hexagon_aromatic(x_as_int, y_as_int, shift_is_pressed);
   if (canvas_addition_mode == lbg_info_t::TRIANLE)
      changed_status = try_stamp_polygon(3, x_as_int, y_as_int, shift_is_pressed, 0);
   if (canvas_addition_mode == lbg_info_t::SQUARE)
      changed_status = try_stamp_polygon(4, x_as_int, y_as_int, shift_is_pressed, 0);
   if (canvas_addition_mode == lbg_info_t::HEPTAGON)
      changed_status = try_stamp_polygon(7, x_as_int, y_as_int, shift_is_pressed, 0);
   if (canvas_addition_mode == lbg_info_t::OCTAGON)
      changed_status = try_stamp_polygon(8, x_as_int, y_as_int, shift_is_pressed, 0);

   if (is_atom_element(canvas_addition_mode)) {
      changed_status = try_change_to_element(canvas_addition_mode);
   }

   if (is_bond(canvas_addition_mode)) {
      changed_status = try_add_or_modify_bond(canvas_addition_mode, x_as_int, y_as_int);
   }

   if (changed_status)
      save_molecule();
}

void
lbg_info_t::handle_item_delete(GdkEventButton *event) {

   int x_as_int, y_as_int;
   GdkModifierType state;
   gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));

   if (highlight_data.has_contents()) {
      if (highlight_data.single_atom()) {
	 mol.close_atom(highlight_data.get_atom_index(), root);
      } else {
	 int ind_1 = highlight_data.get_bond_indices().first;
	 int ind_2 = highlight_data.get_bond_indices().second;
	 int bond_index = mol.get_bond_index(ind_1, ind_2);
	 mol.close_bond(bond_index, root, 1);
      }
      save_molecule();
   }
}


// that's canvas_addition_mode (from the button press).
bool
lbg_info_t::is_bond(int addition_mode) const {
   bool r = 0; 
   if (addition_mode == lbg_info_t::ADD_SINGLE_BOND)
      r = 1;
   if (addition_mode == lbg_info_t::ADD_DOUBLE_BOND)
      r = 1;
   if (addition_mode == lbg_info_t::ADD_TRIPLE_BOND)
      r = 1;
   if (addition_mode == lbg_info_t::ADD_STEREO_OUT_BOND)
      r = 1;
   return r;
}

bool
lbg_info_t::is_atom_element(int addition_mode) const {

   bool r = 0;
   if (addition_mode == lbg_info_t::ATOM_C)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_N)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_O)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_S)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_P)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_F)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_CL)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_I)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_BR)
      r = 1;
   if (addition_mode == lbg_info_t::ATOM_X)
      r = 1;
   return r;
}

std::string
lbg_info_t::to_element(int addition_mode) const {

   std::string r = "";
   if (addition_mode == lbg_info_t::ATOM_C)
      r = "C";
   if (addition_mode == lbg_info_t::ATOM_N)
      r = "N";
   if (addition_mode == lbg_info_t::ATOM_O)
      r = "O";
   if (addition_mode == lbg_info_t::ATOM_S)
      r = "S";
   if (addition_mode == lbg_info_t::ATOM_P)
      r = "P";
   if (addition_mode == lbg_info_t::ATOM_F)
      r = "F";
   if (addition_mode == lbg_info_t::ATOM_CL)
      r = "Cl";
   if (addition_mode == lbg_info_t::ATOM_I)
      r = "I";
   if (addition_mode == lbg_info_t::ATOM_BR)
      r = "Br";
   if (addition_mode == lbg_info_t::ATOM_X)
      r = "X";  // wrong.
   return r;
}


bool
lbg_info_t::try_change_to_element(int addition_element_mode) {

   bool changed_status = 0;
   if (highlight_data.has_contents()) {
      if (highlight_data.single_atom()) {
	 int atom_index = highlight_data.get_atom_index();
	 if (atom_index != UNASSIGNED_INDEX) { 
	    std::string new_ele = to_element(addition_element_mode);
	    changed_status = mol.atoms[atom_index].change_element(new_ele);
	    std::string fc = font_colour(addition_element_mode);
	    if (changed_status) {
	       change_atom_element(atom_index, new_ele, fc);
	    }
	 }
      }
   }
   return changed_status;
}

void
lbg_info_t::change_atom_id_maybe(int atom_index) {

   std::string ele = mol.atoms[atom_index].element;
   std::string fc = font_colour(ele);
   change_atom_element(atom_index, ele, fc);

} 

bool
lbg_info_t::change_atom_element(int atom_index, std::string new_ele, std::string fc) {

   
   bool changed_status = 0;
   std::vector<int> local_bonds = mol.bonds_having_atom_with_atom_index(atom_index);
   lig_build::pos_t pos = mol.atoms[atom_index].atom_position;
	    
   std::string atom_id = mol.make_atom_id_by_using_bonds(new_ele, local_bonds);
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));

   // std::cout << "   calling update_name_maybe(" << atom_name << ") " << std::endl;
   changed_status = mol.atoms[atom_index].update_atom_id_maybe(atom_id, fc, root);
   // std::cout << "update_name_maybe return changed_status: " << changed_status << std::endl;

   if (changed_status) { 
      for (unsigned int ib=0; ib<local_bonds.size(); ib++) {
	 // make a new line for the bond 
	 //
	 int index_1 = mol.bonds[local_bonds[ib]].get_atom_1_index();
	 int index_2 = mol.bonds[local_bonds[ib]].get_atom_2_index();
	 lig_build::atom_t atom_1 = mol.atoms[index_1];
	 lig_build::atom_t atom_2 = mol.atoms[index_2];

	 mol.bonds[local_bonds[ib]].make_new_canvas_item(atom_1, atom_2, root);
      }
   }
   return changed_status;
}


bool
lbg_info_t::try_add_or_modify_bond(int canvas_addition_mode, int x_mouse, int y_mouse) {

   bool changed_status = 0;
   if (! highlight_data.has_contents()) {
      if (mol.is_empty()) {
	 try_stamp_bond_anywhere(canvas_addition_mode, x_mouse, y_mouse);
	 changed_status = 1; // try_stamp_bond_anywhere always modifies
      }
   } else {
      if (highlight_data.single_atom()) {
	 int atom_index = highlight_data.get_atom_index();
	 if (atom_index != UNASSIGNED_INDEX) {
	    changed_status = add_bond_to_atom(atom_index, canvas_addition_mode);
	 }
      } else {
	 // highlighted item was a bond then.
	 int ind_1 = highlight_data.get_bond_indices().first;
	 int ind_2 = highlight_data.get_bond_indices().second;
	 int bond_index = mol.get_bond_index(ind_1, ind_2);
	 lig_build::bond_t::bond_type_t bt =
	    addition_mode_to_bond_type(canvas_addition_mode);
	 if (bond_index != UNASSIGNED_INDEX) {
	    // we need to pass the atoms so that we know if and how to
	    // shorten the bonds (canvas items) to account for atom names.
	    lig_build::atom_t at_1 = mol.atoms[ind_1];
	    lig_build::atom_t at_2 = mol.atoms[ind_2];
	    GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));
	    if (canvas_addition_mode == lbg_info_t::ADD_TRIPLE_BOND)
	       mol.bonds[bond_index].change_bond_order(at_1, at_2, 1, root);
	    else 
	       mol.bonds[bond_index].change_bond_order(at_1, at_2, root); // single to double
						                          // or visa versa

	    // Now that we have modified this bond, the atoms
	    // compromising the bond may need to have their atom name
	    // changed (e.g.) N -> NH,
	    // 
	    // we do this for both of the atoms of the bond.
	    //
	    change_atom_id_maybe(ind_1);
	    change_atom_id_maybe(ind_2);
	    changed_status = 1;
	 }
      }
   }
   return changed_status;
}

lig_build::bond_t::bond_type_t
lbg_info_t::addition_mode_to_bond_type(int canvas_addition_mode) const {
   
   lig_build::bond_t::bond_type_t bt = lig_build::bond_t::SINGLE_BOND;
   if (canvas_addition_mode == lbg_info_t::ADD_DOUBLE_BOND)
      bt = lig_build::bond_t::DOUBLE_BOND;
   if (canvas_addition_mode == lbg_info_t::ADD_TRIPLE_BOND)
      bt = lig_build::bond_t::TRIPLE_BOND;
   if (canvas_addition_mode == lbg_info_t::ADD_STEREO_OUT_BOND)
      bt = lig_build::bond_t::OUT_BOND;
   return bt;
}

bool
lbg_info_t::add_bond_to_atom(int atom_index, int canvas_addition_mode) {

   bool changed_status = 0;
   std::vector<int> bonds = mol.bonds_having_atom_with_atom_index(atom_index);

   // std::cout << "here in add_bond_to_atom there are " << bonds.size() << " neighbours " << std::endl;

   switch (bonds.size()) {

   case 0:
      add_bond_to_atom_with_0_neighbours(atom_index, canvas_addition_mode);
      changed_status = 1;
      break;

   case 1:
      add_bond_to_atom_with_1_neighbour(atom_index, canvas_addition_mode, bonds[0]);
      changed_status = 1;
      break;

   case 2:
      add_bond_to_atom_with_2_neighbours(atom_index, canvas_addition_mode, bonds);
      changed_status = 1;
      break;

   case 3:
      add_bond_to_atom_with_3_neighbours(atom_index, canvas_addition_mode, bonds);
      changed_status = 1;
      break;

   default:
      std::cout << "not handled yet: " << bonds.size() << " neighbouring bonds" << std::endl;
   }
   return changed_status;
}


std::string
lbg_info_t::font_colour(int addition_element_mode) const {
   std::string font_colour = "black";
   if (addition_element_mode == lbg_info_t::ATOM_O)
      font_colour = "red";
   if (addition_element_mode == lbg_info_t::ATOM_N)
      font_colour = "blue";
   if (addition_element_mode == lbg_info_t::ATOM_S)
      font_colour = "#999900";
   if (addition_element_mode == lbg_info_t::ATOM_F)
      font_colour = "#00aa00";
   if (addition_element_mode == lbg_info_t::ATOM_CL)
      font_colour = "#229900";
   if (addition_element_mode == lbg_info_t::ATOM_I)
      font_colour = "#220077";
   return font_colour;
}


std::string
lbg_info_t::font_colour(const std::string &ele) const {
   std::string font_colour = dark;

   if (ele == "N")
      font_colour = "blue";
   if (ele == "O") 
      font_colour = "red";
   if (ele == "S") 
      font_colour = "#999900";
   if (ele == "F") 
      font_colour = "#00aa00";
   if (ele == "CL") 
      font_colour = "#229900";
   if (ele == "I") 
      font_colour = "#220077";
   
   return font_colour;
}

void
lbg_info_t::add_bond_to_atom_with_0_neighbours(int atom_index, int canvas_addition_mode) {

   // certain change
   
   widgeted_atom_t atom = mol.atoms[atom_index];
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));

   lig_build::pos_t current_atom_pos = atom.atom_position;
   lig_build::pos_t a_bond(SINGLE_BOND_CANVAS_LENGTH, 0.0);
   a_bond.rotate(60);
   lig_build::pos_t new_atom_pos = current_atom_pos + a_bond;
      
   widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
   int new_index = mol.add_atom(new_atom).second;
   lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
   widgeted_bond_t b(atom_index, new_index, atom, new_atom, bt, root);
   mol.add_bond(b);
}


void
lbg_info_t::add_bond_to_atom_with_1_neighbour(int atom_index, int canvas_addition_mode,
					      int bond_index) {
   
   int index_1 = mol.bonds[bond_index].get_atom_1_index();
   int index_2 = mol.bonds[bond_index].get_atom_2_index();

   int other_atom_index = index_1;
   if (index_1 == atom_index)
      other_atom_index = index_2;
   
   widgeted_atom_t atom = mol.atoms[atom_index];

   lig_build::pos_t pos_1 = mol.atoms[atom_index].atom_position;
   lig_build::pos_t pos_2 = mol.atoms[other_atom_index].atom_position;

   lig_build::pos_t diff = pos_1 - pos_2;
   lig_build::pos_t du = diff.unit_vector();
   lig_build::pos_t candidate_new_vec_1 = du.rotate( 60);
   lig_build::pos_t candidate_new_vec_2 = du.rotate(-60);

   lig_build::pos_t candidate_pos_1 = pos_1 + candidate_new_vec_1 * SINGLE_BOND_CANVAS_LENGTH;
   lig_build::pos_t candidate_pos_2 = pos_1 + candidate_new_vec_2 * SINGLE_BOND_CANVAS_LENGTH;

   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));

   std::vector<int> avoid_atoms(2);
   avoid_atoms[0] = atom_index;
   avoid_atoms[1] = other_atom_index;

   std::pair<bool,double> d1 = mol.dist_to_other_atoms_except(avoid_atoms, candidate_pos_1);
   std::pair<bool,double> d2 = mol.dist_to_other_atoms_except(avoid_atoms, candidate_pos_2);

   lig_build::pos_t new_atom_pos = candidate_pos_1;
   if (d1.first && d2.first) {
      if (d2.second > d1.second)
	 new_atom_pos = candidate_pos_2;
   }

   // Make an atom, add it, make a bond (using that atom and atom_index) and add it.


   // GooCanvasItem *ci = canvas_line_bond(pos_1, new_atom_pos, root, canvas_addition_mode);
   
   widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
   int new_index = mol.add_atom(new_atom).second;
   lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
   widgeted_bond_t b(atom_index, new_index, atom, new_atom, bt, root);
   mol.add_bond(b);
}

void
lbg_info_t::add_bond_to_atom_with_2_neighbours(int atom_index, int canvas_addition_mode,
					       const std::vector<int> &bond_indices) {

   widgeted_atom_t atom = mol.atoms[atom_index];
   int atom_index_1 = mol.bonds[bond_indices[0]].get_atom_1_index();
   int atom_index_2 = mol.bonds[bond_indices[0]].get_atom_2_index();
   int atom_index_3 = mol.bonds[bond_indices[1]].get_atom_1_index();
   int atom_index_4 = mol.bonds[bond_indices[1]].get_atom_2_index();

   int A_index = atom_index_1;
   if (atom_index_1 == atom_index)
      A_index = atom_index_2;
   int B_index = atom_index_3;
   if (atom_index_3 == atom_index)
      B_index = atom_index_4;

   lig_build::pos_t A_neighb = mol.atoms[A_index].atom_position;
   lig_build::pos_t B_neighb = mol.atoms[B_index].atom_position;
   lig_build::pos_t atom_pos = mol.atoms[atom_index].atom_position;

   lig_build::pos_t mp =
      lig_build::pos_t::mid_point(A_neighb, B_neighb);
   lig_build::pos_t mpa = atom_pos - mp;
   lig_build::pos_t mpa_unit = mpa.unit_vector();

   lig_build::pos_t new_atom_pos = atom_pos + mpa_unit * SINGLE_BOND_CANVAS_LENGTH;
   
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));
   widgeted_atom_t at(new_atom_pos, "C", 0, NULL);
   int new_index = mol.add_atom(at).second;
   
   lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
   widgeted_bond_t b(atom_index, new_index, atom, at, bt, root);
   mol.add_bond(b);


   // Now, what about the atom that has been added to?  Now that we
   // have a extra bond, that atoms name could go from NH -> N (or C
   // -> C+ for 5 bonds).
   //
   // std::cout << "calling change_atom_element(" << atom_index << ") " << std::endl;
   change_atom_id_maybe(atom_index); 
}


void
lbg_info_t::add_bond_to_atom_with_3_neighbours(int atom_index, int canvas_addition_mode,
					       const std::vector<int> &bond_indices) {

   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));

   // are (at least) 2 of the bonds attached to atom_index terminated
   // at their end?
   //
   // If so, we shall remove and replace those bonds before we add
   // this new one.
   //
   std::pair<bool, std::vector<int> > pr = have_2_stubs_attached_to_atom(atom_index, bond_indices);
   std::vector<int> attached_bonds = pr.second;
   
   if (pr.first) {
      int l = attached_bonds.size();
      orthogonalise_2_bonds(atom_index, attached_bonds, bond_indices);
      // now add a third

      lig_build::bond_t existing_bond = mol.bonds[attached_bonds[l-3]];
      lig_build::pos_t p1 =
	 mol.atoms[mol.bonds[bond_indices[attached_bonds[l-3]]].get_other_index(atom_index)].atom_position;
      lig_build::pos_t central_atom_pos = mol.atoms[atom_index].atom_position;
      lig_build::pos_t existing_bond_dir = central_atom_pos - p1;
      lig_build::pos_t ebd_uv = existing_bond_dir.unit_vector();

      lig_build::pos_t new_atom_pos = central_atom_pos + ebd_uv * SINGLE_BOND_CANVAS_LENGTH;
      widgeted_atom_t atom(new_atom_pos, "C", 0, NULL);
      int new_atom_index = mol.add_atom(atom).second;
      lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
      widgeted_bond_t b(atom_index, new_atom_index, mol.atoms[atom_index], atom, bt, root);
      mol.add_bond(b);

   } else {
      // bond_indices are bonds have an atom that is atom_index.
      squeeze_in_a_4th_bond(atom_index, canvas_addition_mode, bond_indices);
   } 
   
}

void
lbg_info_t::squeeze_in_a_4th_bond(int atom_index, int canvas_addition_mode,
				  const std::vector<int> &bond_indices) {

   
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));
   std::vector<double> angles = get_angles(atom_index, bond_indices);
   std::vector<double> sorted_angles = angles;
   std::sort(sorted_angles.begin(), sorted_angles.end());
   if (sorted_angles[0] > 115.0) { // smallest angle, i.e. each perfect 120.
      if (all_closed_rings(atom_index, bond_indices)) { 
	 lig_build::pos_t centre = mol.atoms[atom_index].atom_position;
	 bool found_from_pos = 0;
	 lig_build::pos_t from_pos(0,0); // unfilled
	 for (unsigned int i=0; i<bond_indices.size(); i++) { 
	    int idx = mol.bonds[bond_indices[i]].get_other_index(atom_index);
	    lig_build::pos_t pos = mol.atoms[idx].atom_position;
	    lig_build::pos_t diff = pos - centre;
	    double ori = diff.axis_orientation();
	    if (ori > -20) { 
	       if (ori < 90) { 
		  from_pos = pos;
		  found_from_pos = 1;
	       }
	    }
	 }
	 if (! found_from_pos) {
	    // use the first one.
	    int idx = mol.bonds[bond_indices[0]].get_other_index(atom_index);
	    from_pos = mol.atoms[idx].atom_position;
	 }

	 // now rotate from_pos by 60 degrees
	 lig_build::pos_t diff = from_pos - centre;
	 lig_build::pos_t d_uv = diff.unit_vector();
	 lig_build::pos_t d_uv_60 = d_uv.rotate(-60); // upside down canvas
	 lig_build::pos_t new_atom_pos = centre + d_uv_60 * SINGLE_BOND_CANVAS_LENGTH * 0.7;
	 widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
	 int new_atom_index = mol.add_atom(new_atom).second;
	 lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
	 widgeted_bond_t b(atom_index, new_atom_index, mol.atoms[atom_index], new_atom, bt, root);
	 mol.add_bond(b);

      } else {
	 
	 // OK, there were not centres on all side of all the bonds.
	 // Go place the new atom not towards a ring centre

	 std::cout << "Go place a bond not towards a ring centre" << std::endl;
	 lig_build::pos_t new_atom_pos = get_new_pos_not_towards_ring_centres(atom_index, bond_indices);
	 widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
	 int new_atom_index = mol.add_atom(new_atom).second;
	 lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
	 widgeted_bond_t b(atom_index, new_atom_index, mol.atoms[atom_index], new_atom, bt, root);
	 mol.add_bond(b);
	 
      }
      

   } else {
      lig_build::pos_t new_atom_pos = new_pos_by_bisection(atom_index, bond_indices, angles, root);
      widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
      int new_atom_index = mol.add_atom(new_atom).second;
      lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
      widgeted_bond_t b(atom_index, new_atom_index, mol.atoms[atom_index], new_atom, bt, root);
      mol.add_bond(b);
   } 
}


// bond_indices is a vector of bond indices that have an atom with
// index atom_index.
//
bool
lbg_info_t::all_closed_rings(int atom_index, const std::vector<int> &bond_indices) const {

   bool status = 0;

   if (bond_indices.size() > 3) {
      std::vector<lig_build::pos_t> centres = get_centres_from_bond_indices(bond_indices);
      if (centres.size() > 2)
	 status = 1;
   }
   return status;
}

std::vector<lig_build::pos_t>
lbg_info_t::get_centres_from_bond_indices(const std::vector<int> &bond_indices) const {

   std::vector<lig_build::pos_t> centres;
   for (unsigned int ib=0; ib<bond_indices.size(); ib++) {
      if (mol.bonds[bond_indices[ib]].have_centre_pos()) {
	 lig_build::pos_t test_centre = mol.bonds[bond_indices[ib]].centre_pos();

	 // add test_centre to centres only if it wasnt there already.
	 bool found_centre = 0;
	 for (unsigned int j=0; j<centres.size(); j++) {
	    if (test_centre.close_point(centres[j])) {
	       found_centre = 1; // it was already there
	       break;
	    }
	 }
	 if (! found_centre)
	    centres.push_back(test_centre);
      }
   }
   return centres;
} 

lig_build::pos_t
lbg_info_t::get_new_pos_not_towards_ring_centres(int atom_index,
						 const std::vector<int> &bond_indices) const { 

   lig_build::pos_t centre = mol.atoms[atom_index].atom_position;
   lig_build::pos_t p;
   std::vector<lig_build::pos_t> centres = get_centres_from_bond_indices(bond_indices);
   if (centres.size() == 2) { 
      lig_build::pos_t diff = centres[1] - centres[0];
      lig_build::pos_t diff_uv_90 = diff.unit_vector().rotate(90);
      lig_build::pos_t extra = diff_uv_90 * SINGLE_BOND_CANVAS_LENGTH;
      lig_build::pos_t candidate_1 = centre + extra;
      lig_build::pos_t candidate_2 = centre - extra;
      lig_build::pos_t best_candidate;
      double d_1 = 10000;
      double d_2 = 10000;
      for (unsigned int i=0; i<mol.atoms.size(); i++) { 
	 double dt_1 = lig_build::pos_t::length(candidate_1, mol.atoms[i].atom_position);
	 double dt_2 = lig_build::pos_t::length(candidate_2, mol.atoms[i].atom_position);
	 if (dt_1 < d_1)
	    d_1 = dt_1;
	 if (dt_2 < d_2)
	    d_2 = dt_2;
      }
      if (d_1 < 9999) { 
	 p = candidate_1;
	 if (d_2 > d_1)
	    p = candidate_2;
      }
   } else {
      // build away from ring_centre
      for (unsigned int i=0; i<bond_indices.size(); i++) { 
	 if (mol.bonds[bond_indices[i]].have_centre_pos()) {
	    lig_build::pos_t ring_centre = mol.bonds[bond_indices[i]].centre_pos();
	    int other_index = mol.bonds[bond_indices[i]].get_other_index(atom_index);
	    lig_build::pos_t p1 = mol.atoms[other_index].atom_position;
	    lig_build::pos_t p2 = mol.atoms[atom_index ].atom_position;
	    lig_build::pos_t bond_dir = p2 - p1;
	    lig_build::pos_t bond_dir_uv = bond_dir.unit_vector();
	    p = centre + bond_dir_uv * SINGLE_BOND_CANVAS_LENGTH * 0.8;
	    break;
	 }
      }
   }

   return p;
}


lig_build::pos_t
lbg_info_t::new_pos_by_bisection(int atom_index,
				 const std::vector<int> &bond_indices,
				 const std::vector<double> &angles,
				 GooCanvasItem *root) const {

   lig_build::pos_t p(0,0);
   std::vector<double> sorted_angles = angles;
   std::sort(sorted_angles.begin(), sorted_angles.end());
   lig_build::pos_t centre = mol.atoms[atom_index].atom_position;
   for (unsigned int i=0; i<bond_indices.size(); i++) { 
      int j = i+1;
      if (j == bond_indices.size())
	 j = 0;
      int idx_1 = mol.bonds[bond_indices[i]].get_other_index(atom_index);
      int idx_2 = mol.bonds[bond_indices[j]].get_other_index(atom_index);
      lig_build::pos_t pos_1 = mol.atoms[idx_1].atom_position;
      lig_build::pos_t pos_2 = mol.atoms[idx_2].atom_position;
      lig_build::pos_t d_1 = pos_1 - centre;
      lig_build::pos_t d_2 = pos_2 - centre;
      double cross = lig_build::pos_t::cross(d_1, d_2);
      double cos_theta = cross/(d_1.length() * d_2.length());
      double angle = acos(cos_theta)/DEG_TO_RAD;
      if (angle > (sorted_angles.back()-0.01)) {
	 lig_build::pos_t mp_diff = lig_build::pos_t::mid_point(pos_1, pos_2) - centre;
	 lig_build::pos_t mp_diff_uv = mp_diff.unit_vector();
	 p = centre + mp_diff_uv * SINGLE_BOND_CANVAS_LENGTH * 0.9;
	 break;
      }
   }

   return p;
}

std::vector<double>
lbg_info_t::get_angles(int atom_index, const std::vector<int> &bond_indices) const {

   std::vector<double> v(bond_indices.size());

   lig_build::pos_t centre = mol.atoms[atom_index].atom_position;

   for (unsigned int i=0; i<bond_indices.size(); i++) {
      int j = i+1;
      if (j == bond_indices.size())
	 j = 0;
      int idx_1 = mol.bonds[bond_indices[i]].get_other_index(atom_index);
      int idx_2 = mol.bonds[bond_indices[j]].get_other_index(atom_index);
      lig_build::pos_t pos_1 = mol.atoms[idx_1].atom_position;
      lig_build::pos_t pos_2 = mol.atoms[idx_2].atom_position;
      lig_build::pos_t d_1 = pos_1 - centre;
      lig_build::pos_t d_2 = pos_2 - centre;
      double cross = lig_build::pos_t::cross(d_1, d_2);
      double cos_theta = cross/(d_1.length() * d_2.length());
      double angle = acos(cos_theta)/DEG_TO_RAD;
      v[i] = angle;
   }
   return v;
}


// delete the last 2 bonds of attached_bonds. attached_bonds must be
// of size 3 at least.
//
void
lbg_info_t::orthogonalise_2_bonds(int atom_index,
				  const std::vector<int> &attached_bonds,
				  const std::vector<int> &bond_indices) { 

   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS (canvas));

   int l = attached_bonds.size();
   int ind_1 = attached_bonds[l-1];
   int ind_2 = attached_bonds[l-2];
   std::string ele_1 = mol.atoms[ind_1].element;
   std::string ele_2 = mol.atoms[ind_2].element;
   mol.close_atom(ind_1, root);
   mol.close_atom(ind_2, root);

   lig_build::bond_t existing_bond = mol.bonds[attached_bonds[l-3]];
   lig_build::pos_t p1 =
      mol.atoms[mol.bonds[bond_indices[attached_bonds[l-3]]].get_other_index(atom_index)].atom_position;
   lig_build::pos_t central_atom_pos = mol.atoms[atom_index].atom_position;
   lig_build::pos_t existing_bond_dir = central_atom_pos - p1;
   lig_build::pos_t ebd_uv = existing_bond_dir.unit_vector();
   lig_build::pos_t ebd_uv_90 = ebd_uv.rotate(90);

   std::cout << "existing_bond_dir uv: " << ebd_uv << " and rotated: " << ebd_uv_90
	     << std::endl;
      
   lig_build::pos_t new_atom_pos_1 = central_atom_pos + ebd_uv_90 * SINGLE_BOND_CANVAS_LENGTH;
   lig_build::pos_t new_atom_pos_2 = central_atom_pos - ebd_uv_90 * SINGLE_BOND_CANVAS_LENGTH;
   widgeted_atom_t atom_1(new_atom_pos_1, ele_1, 0, NULL);
   widgeted_atom_t atom_2(new_atom_pos_2, ele_2, 0, NULL);
   int new_atom_index_1 = mol.add_atom(atom_1).second;
   int new_atom_index_2 = mol.add_atom(atom_2).second;
   lig_build::bond_t::bond_type_t bt = addition_mode_to_bond_type(canvas_addition_mode);
   widgeted_bond_t b_1(atom_index, new_atom_index_1, mol.atoms[atom_index], atom_1, bt, root);
   widgeted_bond_t b_2(atom_index, new_atom_index_2, mol.atoms[atom_index], atom_2, bt, root);
   mol.add_bond(b_1);
   mol.add_bond(b_2);

}



// return a status and a vector of atoms (bonded to atom_index) having
// only one bond.
// 
std::pair<bool, std::vector<int> > 
lbg_info_t::have_2_stubs_attached_to_atom(int atom_index, const std::vector<int> &bond_indices) const {

   std::vector<int> v;
   for (unsigned int i=0; i<bond_indices.size(); i++) { 
      int other_index = mol.bonds[bond_indices[i]].get_other_index(atom_index);
      // now, does other_index have only one bond?
      std::vector<int> local_bonds = mol.bonds_having_atom_with_atom_index(other_index);
      if (local_bonds.size() == 1) {
	 v.push_back(other_index);
      }
   }
   bool status = 0;
   if (v.size() > 1)
      status = 1;
   return std::pair<bool, std::vector<int> > (status, v);
}



// Set the addition mode and turn off other toggle buttons (if
// needed).
void
lbg_handle_toggle_button(GtkToggleToolButton *tb, GtkWidget *canvas, int mode) {

   gpointer gp = gtk_object_get_user_data(GTK_OBJECT(canvas));
   lbg_info_t *l = static_cast<lbg_info_t *> (gp);

   if (l) { 
      if (gtk_toggle_tool_button_get_active(tb)) {

	 l->untoggle_others_except(tb);
	 l->canvas_addition_mode = mode;
      } else {
	 l->canvas_addition_mode = lbg_info_t::NONE;
      }
   } else {
      std::cout << "ERROR:: in lbg_handle_toggle_button() failed to get lbg pointer!"
		<< std::endl;
   } 
}


bool
lbg_info_t::try_stamp_hexagon_aromatic(int x_cen, int y_cen, bool shift_is_pressed) {

   bool is_aromatic = 1;
   return try_stamp_polygon(6, x_cen, y_cen, shift_is_pressed, is_aromatic);
}

// shift_is_pressed implies spiro.  Return changed status.
bool
lbg_info_t::try_stamp_polygon(int n_edges, int x_cen, int y_cen, bool is_spiro, bool is_aromatic) {

   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));
   bool changed_status = 0;

   if (mol.is_empty()) {
      stamp_polygon_anywhere(n_edges, x_cen, y_cen, is_aromatic, root);
      changed_status = 1; // stamp_polygon_anywhere() always modifies
   } else { 
      if (highlight_data.has_contents()) {
	 std::vector<int> new_atoms =
	    try_stamp_polygon_using_highlighted_data(n_edges, is_spiro, is_aromatic, root);
	 if (new_atoms.size() > 0)
	    changed_status = 1;
	 if (highlight_data.single_atom()) {
	    if (new_atoms.size() > 0) {
	       // we need to add the spur bond, if this is not a spiro compound

	       if (! is_spiro) {

		  // what are the atom indices of the spur bond?

		  int highlighted_atom_index  = highlight_data.get_atom_index();
		  if (highlighted_atom_index != new_atoms[0]) { 
		     std::cout << "adding spur bond...!" << std::endl;
		     
		     if (highlighted_atom_index != UNASSIGNED_INDEX) {
			lig_build::bond_t::bond_type_t bt = lig_build::bond_t::SINGLE_BOND;
			widgeted_bond_t bond(highlighted_atom_index, new_atoms[0],
					     mol.atoms[highlighted_atom_index],
					     mol.atoms[new_atoms[0]],
					     bt, root);
			mol.add_bond(bond);
		     }
		  }
	       }
	    }
	 }
      }
   }
   return changed_status;
}

// always modifies
// 
void
lbg_info_t::try_stamp_bond_anywhere(int canvas_addition_mode, int x_mouse, int y_mouse) {

   bool changed_status = 0;

   double dx = double(x_mouse);
   double dy = double(y_mouse);
   lig_build::pos_t new_atom_pos(dx, dy);

    widgeted_atom_t new_atom(new_atom_pos, "C", 0, NULL);
    int new_atom_index = mol.add_atom(new_atom).second;
    add_bond_to_atom(new_atom_index, canvas_addition_mode);
}


double
lbg_info_t::radius(int n_edges) const {

   double angle_step = 360.0/double(n_edges);
   double r = SINGLE_BOND_CANVAS_LENGTH/(2*sin(angle_step*DEG_TO_RAD*0.5));
   return r;
}

// always modifies
void
lbg_info_t::stamp_polygon_anywhere(int n_edges, int x_cen, int y_cen,
				   bool is_aromatic,
				   GooCanvasItem *root) {

   double angle_off = 0;
   lig_build::polygon_position_info_t ppi(x_cen, y_cen, angle_off);
   stamp_polygon(n_edges, ppi, is_aromatic, root);
}


// is_aromatic means (for now) that double and single bonds alternate.
// 
// always modifies
// 
std::vector<int>
lbg_info_t::stamp_polygon(int n_edges, lig_build::polygon_position_info_t ppi,
			  bool is_aromatic, GooCanvasItem *root) {

   // use the next two to construct the returned value.
   std::vector<int> atom_indices;
      
   double angle_step = 360.0/double(n_edges);
   double dx_cen = double(ppi.pos.x);
   double dy_cen = double(ppi.pos.y);
   lig_build::pos_t centre(dx_cen, dy_cen);
   std::vector<std::pair<lig_build::pos_t, std::pair<bool,int> > > atom_pos_index;

   // This makes the bond length the same, no matter how many edges.
   //
   double r = radius(n_edges);

   // offset so that the first bond lies along the X axis:
   // 
   // double internal_angle_offset = -angle_step/2.0;
   // double orientation_correction = 180;
   
   double internal_angle_offset = 0;
   double orientation_correction = 180;
   if (! ppi.apply_internal_angle_offset_flag) { 
      internal_angle_offset = 0;
      orientation_correction = 0;
   }
   
   for (int i=0; i<n_edges; i++) {
      double theta_deg =
	 (i * angle_step + orientation_correction + internal_angle_offset + ppi.angle_offset);
      double theta = theta_deg * DEG_TO_RAD;
      double pt_1_x = r * sin(theta) + dx_cen;
      double pt_1_y = r * cos(theta) + dy_cen;
      lig_build::pos_t pt(pt_1_x, pt_1_y);
      GooCanvasItem *aci = 0; // initially carbon, no canvas item
      widgeted_atom_t at(pt, "C", 0, aci);
      std::pair<bool,int> atom_index = mol.add_atom(at);
      std::pair<lig_build::pos_t, std::pair<bool,int> > pr(pt, atom_index);
      atom_pos_index.push_back(pr);
      if (0) { // debug polygon vertex placement
	 std::string stroke_colour = get_stroke_colour(i,atom_pos_index.size());
	 GooCanvasItem *rect_item  =
	    goo_canvas_rect_new (root,
				 pt.x - 5.0, 
				 pt.y - 5.0,
				 10.0, 10.0,
				 "line-width", 7.0,
				 "stroke-color", stroke_colour.c_str(),
				 NULL);
      }
   }
   
   for (unsigned int i=0; i<atom_pos_index.size(); i++) {
      int j = i + 1;
      if (j == atom_pos_index.size())
	 j = 0;
      int i_index = atom_pos_index[i].second.second;
      int j_index = atom_pos_index[j].second.second;


      // if either the first or the second atoms were new atoms then
      // draw the bond (and of course, of the first and second atom
      // indexes were both from pre-existing atoms then *don't* draw
      // the bond).
      // 
      if (atom_pos_index[i].second.first || atom_pos_index[j].second.first) { 
	 lig_build::bond_t::bond_type_t bt = lig_build::bond_t::SINGLE_BOND;
	 if (is_aromatic) {
	    if (((i/2)*2) == i)
	       bt = lig_build::bond_t::DOUBLE_BOND;
	 }
	 widgeted_bond_t bond(i_index, j_index,
			      mol.atoms[atom_pos_index[i].second.second],
			      mol.atoms[atom_pos_index[j].second.second], centre, bt, root);
	 mol.add_bond(bond);

	 // is this (the bond vector) needed?
	 //
	 // No it isn't. Remove std::vector<lig_build::bond_t> > from the
	 // return type of this function.
	 // 
      }
   }

   // transfer from atom_pos_index to atom_indices
   atom_indices.resize(atom_pos_index.size());
   for (unsigned int i=0; i<atom_pos_index.size(); i++)
      atom_indices[i] = atom_pos_index[i].second.second;
   
   return atom_indices;
}


std::vector<int>
lbg_info_t::try_stamp_polygon_using_highlighted_data(int n_edges,
						     bool spiro_flag,
						     bool is_aromatic,
						     GooCanvasItem *root) {

   std::vector<int> new_atoms;
   
   double alpha = 2 * M_PI/double(n_edges);
   double corr_r = radius(n_edges) * cos(alpha/2.0);
   double rad_std = radius(n_edges);
   
   lig_build::polygon_position_info_t ppi =
      highlight_data.get_new_polygon_centre(n_edges, spiro_flag, rad_std, corr_r, mol);
   if (ppi.can_stamp) 
      new_atoms = stamp_polygon(n_edges, ppi, is_aromatic, root);
   return new_atoms;
}

lig_build::polygon_position_info_t
lbg_info_t::highlight_data_t::get_new_polygon_centre(int n_edges, bool spiro_flag,
						     const double &rad_std,
						     const double &cor_radius,
						     const widgeted_molecule_t &mol) const {

   lig_build::polygon_position_info_t ppi;

   if (n_atoms_ == 2)
      ppi = get_new_polygon_centre_using_2_atoms(n_edges, cor_radius);

   if (n_atoms_ == 1)
      ppi = get_new_polygon_centre_using_1_atom(n_edges, spiro_flag, rad_std, cor_radius, mol);

   return ppi;
}


lig_build::polygon_position_info_t
lbg_info_t::highlight_data_t::get_new_polygon_centre_using_2_atoms(int n_edges,
								   const double &cor_radius) const {

   // we pass the cor_radius, how much to move the centre of the new
   // polygon from the mid point of "this" bond (which is not the
   // radius because the bond mid-point is "inside" the circle for the
   // new polygon points.

   double x_cen = 0;
   double y_cen = 0;

   // Make a bisector (vector) of the bond.  The new polygon centre
   // should lie on that bisector.
   //
   // The distance along the bisector depends on the number of edges.
   //
   // Bond is A->B, M is the midpoint.  M->P is the bisector.  C, the
   // centre of the new polygon lies along M->P, at r from M. Let's
   // say that M->P is a unit vector, d is |A->M|.
   //
   //
   //              A
   //              |\                           . 
   //              | \                               .
   //              |  \ .
   //              |   \   .
   //            d |    \           .
   //              |     \                  .
   //              |      \              .
   //            M |-------\ P ................... C
   //              |   1 
   //              |
   //              |
   //              |
   //              |
   //              B
   //
   // So the question here is: where is P?
   // 
   //    Take a unit vector along A->B.  Rotate it 90 degress about M.
   //
   //

   lig_build::pos_t a = pos_1_;
   lig_build::pos_t b = pos_2_;
   lig_build::pos_t m = lig_build::pos_t::mid_point(a,b);
   lig_build::pos_t ab_unit = (b - a).unit_vector();
   lig_build::pos_t mp_unit = ab_unit.rotate(90);
   lig_build::pos_t c_1 = m + mp_unit * cor_radius;
   lig_build::pos_t c_2 = m + mp_unit * cor_radius;

   lig_build::pos_t c = c_1;
   if (has_ring_centre_flag) {
      double l1 = lig_build::pos_t::length(c_1, ring_centre);
      double l2 = lig_build::pos_t::length(c_2, ring_centre);
      if (l1 < l2)
	 c = c_2;
   }

   //  0    for a triangle : 120
   //  45   for a square   : 90
   //  0    for pentagon   : 72
   // 30    for hexagon    : 60
   // 51.43 for pentagon   : 51.43
   // 22.5 for octagon : 45

   // orientation correction
   // 
   // double oc = 0;

   lig_build::pos_t ab = b - a;

   // canvas is upside down:
   lig_build::pos_t ab_cor(ab.x, -ab.y);
   
   double ao =  ab_cor.axis_orientation();

   // std::cout << "axis orientation of " << ab_cor << " is " << ao << std::endl;

   double angle_step = 360.0/double(n_edges);
   double internal_angle_offset = -angle_step/2.0;

   double angle_off = ao + internal_angle_offset;
   
   lig_build::polygon_position_info_t ppi(c, angle_off);
   return ppi;
   
}

// This creates a polygon on a spur (unless spiro_flag is set)
// 
lig_build::polygon_position_info_t
lbg_info_t::highlight_data_t::get_new_polygon_centre_using_1_atom(int n_edges,
								  bool spiro_flag,
								  const double &radius,
								  const double &cor_radius,
								  const widgeted_molecule_t &mol) const {

   lig_build::polygon_position_info_t ppi(pos_1_, 0);
   lig_build::pos_t A = pos_1_;

   int atom_index = get_atom_index();
   std::vector<int> bv = mol.bonds_having_atom_with_atom_index(atom_index);

   std::cout << "in get_new_polygon_centre_using_1_atom() bv size is " << bv.size() << std::endl;

   if (bv.size() == 2) {
      std::vector<lig_build::pos_t> neighbours;
      lig_build::pos_t test_pos = mol.atoms[mol.bonds[bv[0]].get_atom_1_index()].atom_position;
      if (! test_pos.near_point(A, 2))
	 neighbours.push_back(test_pos);
      test_pos = mol.atoms[mol.bonds[bv[0]].get_atom_2_index()].atom_position;
      if (! test_pos.near_point(A, 2))
	 neighbours.push_back(test_pos);
      test_pos = mol.atoms[mol.bonds[bv[1]].get_atom_1_index()].atom_position;
      if (! test_pos.near_point(A, 2))
	 neighbours.push_back(test_pos);
      test_pos = mol.atoms[mol.bonds[bv[1]].get_atom_2_index()].atom_position;
      if (! test_pos.near_point(A, 2))
	 neighbours.push_back(test_pos);

      if (neighbours.size() == 2) {
	 lig_build::pos_t mp =
	    lig_build::pos_t::mid_point(neighbours[0], neighbours[1]);
	 lig_build::pos_t mpa = A - mp;
	 lig_build::pos_t mpa_unit = mpa.unit_vector();

	 lig_build::pos_t new_centre = A + mpa_unit * radius;
	 if (! spiro_flag) {
	    new_centre += mpa_unit * SINGLE_BOND_CANVAS_LENGTH;
	 } 

	 ppi.pos = new_centre;
	 lig_build::pos_t bond_vector = neighbours[1] - neighbours[0];

	 // upside-down canvas
	 double angle = -mpa.axis_orientation();
	 // std::cout << "spur orientation " << angle << std::endl;
	 
	 ppi.angle_offset = angle - 90;
	 ppi.apply_internal_angle_offset_flag = 0; // don't internally correct orientation
      }
   }

   // now, say we are adding a ring to the end of a chain of single bonds
   // 
   if (bv.size() == 1) {
      // We need the coordinates of the bond atom, We have one of
      // those (A), what is the other?
      lig_build::pos_t pos_other = mol.atoms[mol.bonds[bv[0]].get_atom_1_index()].atom_position;
      if (A.close_point(pos_other))
	 pos_other = mol.atoms[mol.bonds[bv[0]].get_atom_2_index()].atom_position;
      lig_build::pos_t other_to_A = (A - pos_other);
      lig_build::pos_t other_to_A_unit = other_to_A.unit_vector();
      lig_build::pos_t new_centre = A + other_to_A_unit * radius;
      double angle = -other_to_A_unit.axis_orientation();
      ppi.pos = new_centre;
      ppi.angle_offset = angle - 90;
      ppi.apply_internal_angle_offset_flag = 0; // don't internally correct orientation
   }
   return ppi;
}


bool
lbg_info_t::save_togglebutton_widgets(GtkBuilder *builder) {

   std::vector<std::string> w_names;
   // w_names.push_back("charge_toggle_toolbutton");

   w_names.push_back("single_toggle_toolbutton");
   w_names.push_back("double_toggle_toolbutton");
   w_names.push_back("triple_toggle_toolbutton");
   w_names.push_back("stereo_out_toggle_toolbutton");
   w_names.push_back("c3_toggle_toolbutton");
   w_names.push_back("c4_toggle_toolbutton");
   w_names.push_back("c5_toggle_toolbutton");
   w_names.push_back("c6_toggle_toolbutton");
   w_names.push_back("c6_arom_toggle_toolbutton");
   w_names.push_back("c7_toggle_toolbutton");
   w_names.push_back("c8_toggle_toolbutton");
   w_names.push_back("carbon_toggle_toolbutton");
   w_names.push_back("nitrogen_toggle_toolbutton");
   w_names.push_back("oxygen_toggle_toolbutton");
   w_names.push_back("sulfur_toggle_toolbutton");
   w_names.push_back("phos_toggle_toolbutton");
   w_names.push_back("fluorine_toggle_toolbutton");
   w_names.push_back("chlorine_toggle_toolbutton");
   w_names.push_back("bromine_toggle_toolbutton");
   w_names.push_back("iodine_toggle_toolbutton");
   w_names.push_back("other_element_toggle_toolbutton");
   w_names.push_back("delete_item_toggle_toolbutton");

   for (unsigned int i=0; i<w_names.size(); i++) {
      GtkToggleToolButton *tb =
	 GTK_TOGGLE_TOOL_BUTTON(gtk_builder_get_object (builder, w_names[i].c_str()));
      widget_names[w_names[i]] = tb;
   }
   return TRUE;
} 



void
lbg_info_t::clear() {

   clear_canvas();

   // clear the molecule
   mol.clear();

}

void
lbg_info_t::clear_canvas() {

   // clear the canvas
   GooCanvasItem *root = goo_canvas_get_root_item(GOO_CANVAS(canvas));

   gint n_children = goo_canvas_item_get_n_children (root);
   for (int i=0; i<n_children; i++)
      goo_canvas_item_remove_child(root, 0);
}




void
lbg_info_t::init(GtkBuilder *builder) {

   GtkWidget *window = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_window"));
   gtk_widget_show (window);

   about_dialog =    GTK_WIDGET (gtk_builder_get_object (builder, "lbg_aboutdialog"));
   search_combobox = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_search_combobox"));
   open_dialog     = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_open_filechooserdialog"));
   save_as_dialog  = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_save_as_filechooserdialog"));
   lbg_sbase_search_results_dialog = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_sbase_search_results_dialog"));
   lbg_sbase_search_results_vbox = GTK_WIDGET (gtk_builder_get_object (builder, "lbg_sbase_search_results_vbox"));
   lbg_export_as_pdf_dialog =  GTK_WIDGET (gtk_builder_get_object (builder, "lbg_export_as_pdf_filechooserdialog"));
   lbg_export_as_png_dialog =  GTK_WIDGET (gtk_builder_get_object (builder, "lbg_export_as_png_filechooserdialog"));


   GtkWidget *lbg_scrolled_win =
      GTK_WIDGET(gtk_builder_get_object (builder, "lbg_scrolledwindow"));

   canvas = goo_canvas_new();

   // gtk_widget_set(GTK_WIDGET(canvas), "automatic-bounds",  1, NULL);
   gtk_widget_set(GTK_WIDGET(canvas), "bounds-padding", 50.0, NULL);
   gtk_object_set_user_data(GTK_OBJECT(canvas), (gpointer) this);
   // std::cout << "attached this lbg_info_t pointer to canvas: " << this << std::endl;
   
   save_togglebutton_widgets(builder);
   gtk_container_add(GTK_CONTAINER(lbg_scrolled_win), canvas);
   goo_canvas_set_bounds (GOO_CANVAS (canvas), 0, 0, 1200, 700);

   GooCanvas *gc = GOO_CANVAS(canvas);

   if (0) 
      std::cout << "   after setting bounds: canvas adjustments: h-lower " 
		<< gc->hadjustment->lower << " h-upper-page "
		<< gc->hadjustment->upper - gc->hadjustment->page_size  << " v-lower "
		<< gc->vadjustment->lower  << " v-upper-page"
		<< gc->vadjustment->upper - gc->vadjustment->page_size  << " h-page "
		<< gc->hadjustment->page_size  << " "
		<< gc->vadjustment->page_size << std::endl;

   g_signal_connect(canvas, "button_press_event",
		    (GtkSignalFunc) on_canvas_button_press, NULL);

   g_signal_connect(canvas, "motion_notify_event",
		    (GtkSignalFunc) on_canvas_motion, NULL);
   gtk_widget_show(canvas);

   // search combobox
   add_search_combobox_text();


   // ------------ sbase ---------------------------
   // 
   init_sbase(".");

   
   // ------------ watch for new files from coot ---------------------------
   // 

   int timeout_handle = gtk_timeout_add(500, watch_for_mdl_from_coot, this);

}

// static
gboolean
lbg_info_t::watch_for_mdl_from_coot(gpointer user_data) {

   lbg_info_t *l = static_cast<lbg_info_t *> (user_data);
   
   // the pdb file has the atom names.  When we read the mol file, we
   // will stuff the atom names into the widgeted_molecule (by
   // matching the coordinates).
   //

   std::string coot_dir = "../../build-coot-ubuntu-64bit/src"; // make this user-settable.
   std::string coot_ccp4_dir = coot_dir + "/coot-ccp4";
   // in use (non-testing), coot_dir will typically be ".";

   
   std::string coot_mdl_file = coot_ccp4_dir + "/.coot-to-lbg-mol";
   std::string coot_pdb_file = coot_ccp4_dir + "/.coot-to-lbg-pdb";
   std::string ready_file    = coot_ccp4_dir + "/.coot-to-lbg-mol-ready";
   std::string sa_file       = coot_dir      + "/coot-tmp-fle-view-solvent-accessibilites.txt";

   struct stat buf;
   int err = stat(ready_file.c_str(), &buf);
   if (! err) {
      time_t m = buf.st_mtime;
      // std::cout << "here 1 in watch_for_mdl_from_coot" << std::endl;
      if (m > l->coot_mdl_ready_time) {
	 if (l->coot_mdl_ready_time != 0) {

	    // std::cout << "here 2 in watch_for_mdl_from_coot" << std::endl;
	    CMMDBManager *flat_pdb_mol = l->get_cmmdbmanager(coot_pdb_file);
	    lig_build::molfile_molecule_t mm;
	    mm.read(coot_mdl_file);
	    widgeted_molecule_t wmol = l->import(mm, coot_mdl_file, flat_pdb_mol);
	    std::vector<solvent_accessible_atom_t> solvent_accessible_atoms =
	     l->read_solvent_accessibilities(sa_file);
	    wmol.map_solvent_accessibilities_to_atoms(solvent_accessible_atoms);
	    l->render_from_molecule(wmol);
	 }
	 l->coot_mdl_ready_time = m;
      }
   } else {
      // std::cout << "failed to stat " << ready_file << std::endl;
   }
   return 1; // keep running
}

CMMDBManager *
lbg_info_t::get_cmmdbmanager(const std::string &file_name) const {

   CMMDBManager *mol = new CMMDBManager;

   int err = mol->ReadCoorFile(file_name.c_str());
   if (err) {
      std::cout << "WARNING:: Problem reading coordinates file " << file_name << std::endl;
      delete mol;
      mol = NULL;
   } 

   return mol;

}

int
main(int argc, char *argv[]) {

   InitMatType(); // mmdb program. 
   
   gtk_init (&argc, &argv);
        
   std::string glade_file = "lbg.glade";

   bool glade_file_exists = 0;
   struct stat buf;
   int err = stat(glade_file.c_str(), &buf);
   if (! err)
      glade_file_exists = 1;

   if (glade_file_exists) { 
      
      GtkBuilder *builder = gtk_builder_new ();
      gtk_builder_add_from_file (builder, glade_file.c_str(), NULL);
      lbg_info_t *lbg = new lbg_info_t;
      lbg->init(builder);
      gtk_builder_connect_signals (builder, lbg->canvas);
      g_object_unref (G_OBJECT (builder));

      if (argc > 1) {
	 std::string file_name(argv[1]);
	 lig_build::molfile_molecule_t mm;
	 CMMDBManager *mol = NULL; // no atom names to transfer
	 mm.read(file_name);
	 widgeted_molecule_t wmol = lbg->import(mm, file_name, mol);
	 lbg->render_from_molecule(wmol);
      }
	 
      gtk_main ();

   } else {
      std::cout << "ERROR:: glade file " << glade_file << " not found" << std::endl;
      return 1; 
   } 
   return 0;
}


void
lbg_info_t::add_search_combobox_text() const {

   GtkTreeIter   iter;
   GtkListStore *list_store_similarities =
      gtk_list_store_new (1, G_TYPE_STRING);
   gtk_combo_box_set_model(GTK_COMBO_BOX(search_combobox), GTK_TREE_MODEL(list_store_similarities));
   for (unsigned int i=0; i<6; i++) {
      double f = 0.75 + i*0.05;
      std::string s = coot::util::float_to_string(f);
      gtk_list_store_append(GTK_LIST_STORE(list_store_similarities), &iter);
      gtk_list_store_set(GTK_LIST_STORE(list_store_similarities), &iter,
			 0, s.c_str(), -1);
   }
   gtk_combo_box_set_active(GTK_COMBO_BOX(search_combobox), 0);

   GtkCellRenderer *renderer = gtk_cell_renderer_text_new ();
   gtk_cell_layout_pack_start (GTK_CELL_LAYOUT (search_combobox), renderer, TRUE);
   gtk_cell_layout_add_attribute (GTK_CELL_LAYOUT (search_combobox), renderer, "text", 0);
   

}

std::string
lbg_info_t::get_stroke_colour(int i, int n) const {

   std::string s("fedcba9876543210");

   double f =  1 + 12 * double(i)/double(n);
   int f_int = int(f);

   char c = s[f];
   std::string r = "#";
   for (i=0; i<6; i++)
      r += c;
   return r;
}

void
lbg_info_t::render_from_molecule(const widgeted_molecule_t &mol_in) {


   make_saves_mutex = 0; // stop saving changes (restored at end)
   clear();
   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   
   int re_index[mol_in.atoms.size()]; // map form mol_in atom indexing
				      // the this molecule atom
				      // indexing.
   for (unsigned int i=0; i<mol_in.atoms.size(); i++)
      re_index[i] = UNASSIGNED_INDEX;

   // add in atoms
   for (unsigned int iat=0; iat<mol_in.atoms.size(); iat++) {
      if (!mol_in.atoms[iat].is_closed()) { 
	 GooCanvasItem *ci = NULL;
	 lig_build::pos_t pos = mol_in.atoms[iat].atom_position;
	 widgeted_atom_t new_atom = widgeted_atom_t(pos,
						    mol_in.atoms[iat].element,
						    mol_in.atoms[iat].charge,
						    ci);
	 double sa = mol_in.atoms[iat].get_solvent_accessibility();
	 new_atom.set_atom_name(mol_in.atoms[iat].get_atom_name());
	 if (0)
	    std::cout << "in render_from_molecule() atom with name :"
		      << mol_in.atoms[iat].get_atom_name()
		      << ": has solvent_accessibility " << sa << std::endl;
	 new_atom.add_solvent_accessibility(sa);
	 re_index[iat] = mol.add_atom(new_atom).second;
	 
	 if (0)
	    // clever c++
	    std::cout << " in render_from_molecule: old sa: "
		      << mol_in.atoms[iat].get_solvent_accessibility() << " new: "
		      << mol.atoms[re_index[iat]].get_solvent_accessibility() << std::endl;
	 if (sa > 0) {
	    // std::cout << "draw solvent accessibility " << sa << " at " << pos << std::endl;
	    draw_solvent_accessibility_of_atom(pos, sa, root);
	 }
      }
   }


   // add in bonds
   for (unsigned int ib=0; ib<mol_in.bonds.size(); ib++) {
      int idx_1 = re_index[mol_in.bonds[ib].get_atom_1_index()];
      int idx_2 = re_index[mol_in.bonds[ib].get_atom_2_index()];
      if ((idx_1 != UNASSIGNED_INDEX) && (idx_2 != UNASSIGNED_INDEX)) { 
	 lig_build::bond_t::bond_type_t bt = mol_in.bonds[ib].get_bond_type();
	 if (mol_in.bonds[ib].have_centre_pos()) {
	    lig_build::pos_t centre_pos = mol_in.bonds[ib].centre_pos();
	    widgeted_bond_t bond(idx_1, idx_2, mol.atoms[idx_1], mol.atoms[idx_2], centre_pos, bt, root);
	    mol.add_bond(bond);
	 } else {
	    widgeted_bond_t bond(idx_1, idx_2, mol.atoms[idx_1], mol.atoms[idx_2], bt, root);
	    mol.add_bond(bond);
	 }
      }
   }

   // redo the atoms, this time with widgets.
   for (unsigned int iat=0; iat<mol_in.atoms.size(); iat++) {

      std::vector<int> local_bonds = mol.bonds_having_atom_with_atom_index(iat);
      std::string ele = mol.atoms[iat].element;
      std::string atom_id = mol.make_atom_id_by_using_bonds(ele, local_bonds);
      std::string fc = font_colour(ele);
      if (ele != "C") 
	 mol.atoms[iat].update_atom_id_forced(atom_id, fc, root);
   }

   // for input_coords_to_canvas_coords() to work:
   // 
   mol.centre_correction = mol_in.centre_correction;
   mol.scale_correction  = mol_in.scale_correction;
   mol.mol_in_min_y = mol_in.mol_in_min_y;
   mol.mol_in_max_y = mol_in.mol_in_max_y;
   
   // 
   make_saves_mutex = 1; // allow saves again.
}

void
lbg_info_t::undo() {

   save_molecule_index--;
   if (save_molecule_index >= 0) { 
      widgeted_molecule_t saved_mol = previous_molecules[save_molecule_index];
      std::cout << "undo... reverting to save molecule number " << save_molecule_index << std::endl;
      render_from_molecule(saved_mol);
   } else {
      clear();
   } 
} 



void
lbg_info_t::write_pdf(const std::string &file_name) const { 

#if CAIRO_HAS_PDF_SURFACE

   cairo_surface_t *surface;
   cairo_t *cr;

   surface = cairo_pdf_surface_create(file_name.c_str(), 9 * 72, 10 * 72);
   cr = cairo_create (surface);

   /* Place it in the middle of our 9x10 page. */
   cairo_translate (cr, 20, 130);

   goo_canvas_render (GOO_CANVAS(canvas), cr, NULL, 1.0);
   cairo_show_page (cr);
   cairo_surface_destroy (surface);
   cairo_destroy (cr);

#else
   std::cout << "No PDF (no PDF Surface in Cairo)" << std::endl;
#endif    

}

void
lbg_info_t::write_png(const std::string &file_name) const {

   cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 700, 700);
   cairo_t *cr = cairo_create (surface);

   double scale_factor = 2.5;
   goo_canvas_render (GOO_CANVAS(canvas), cr, NULL, scale_factor); // sc doesn't do anything?
   cairo_surface_write_to_png(surface, file_name.c_str());
   cairo_surface_destroy (surface);
   cairo_destroy (cr);
}



void
lbg_info_t::save_molecule() {

   if (make_saves_mutex) { 
      previous_molecules.push_back(mol);
      save_molecule_index = previous_molecules.size() - 1;
      // std::cout << "saved molecule to index: " << save_molecule_index << std::endl;
   } else {
      std::cout << "debug:: save_molecule() excluded" << std::endl;
   }
}

// pdb_mol is the pdb representation of the (flat) ligand - and it has
// the atom names.  We will add the atom names into mol by matching
// coordinates.
// 
widgeted_molecule_t
lbg_info_t::import(const lig_build::molfile_molecule_t &mol_in, const std::string &file_name,
		   CMMDBManager *pdb_mol) {

   widgeted_molecule_t new_mol(mol_in, pdb_mol);
   mdl_file_name = file_name;
   save_molecule();
   return new_mol;
}


void
lbg_info_t::read_draw_residues(const std::string &file_name) {

   // read residues need to happen after the ligand has been placed on the canvas
   // i.e. mol.input_coords_to_canvas_coords() is called in read_residues();
   // 
   residue_circles = read_residues(file_name);
   std::vector<residue_circle_t> current_circles = residue_circles;

   for (int iround=0; iround<120; iround++) {
      std::cout << ":::::::::::::::::::: minimization round " << iround
		<< " :::::::::::::::::::::::::::::::::::::::::::\n";
      std::pair<int, std::vector<lbg_info_t::residue_circle_t> > new_c =
	 optimise_residue_circle_positions(residue_circles, current_circles);
      current_circles = new_c.second;
      // draw_residue_circles(current_circles);
      if (new_c.first == GSL_ENOPROG)
	 break;
      if (new_c.first == GSL_SUCCESS) { 
	 break;
      }
   }

   // save the class member data
   residue_circles = current_circles;
   draw_bonds_to_ligand();
   draw_residue_circles(current_circles);
}

// minimise layout energy
std::pair<int, std::vector<lbg_info_t::residue_circle_t> >
lbg_info_t::optimise_residue_circle_positions(const std::vector<lbg_info_t::residue_circle_t> &r,
					      const std::vector<lbg_info_t::residue_circle_t> &c) const { 

   if (r.size() > 0) {
      if (c.size() == r.size()) { 
	 optimise_residue_circles orc(r, c, mol);
	 int status = orc.get_status();
	 return orc.solution();
      }
   }
   std::vector<lbg_info_t::residue_circle_t> dv; // dummy 
   return std::pair<int, std::vector<lbg_info_t::residue_circle_t> > (0, dv);
}


std::vector<lbg_info_t::residue_circle_t>
lbg_info_t::read_residues(const std::string &file_name) const {
      
   std::vector<residue_circle_t> v;
   std::ifstream f(file_name.c_str());
   if (!f) {
      std::cout << "Failed to open " << file_name << std::endl;
   } else {

      std::cout << "opened residue_circle file: " << file_name << std::endl;
      
      std::vector<std::string> lines;
      std::string line;
      while (std::getline(f, line)) { 
	 lines.push_back(line);
      }

      for (unsigned int i=0; i<lines.size(); i++) {
	 std::vector<std::string> words = coot::util::split_string_no_blanks(lines[i], " ");
	 
	 // debug input
	 if (0) { 
	    std::cout << i << " " << lines[i] << "\n";
	    for (unsigned int j=0; j<words.size(); j++) { 
	       std::cout << "  " << j << " " << words[j] << " ";
	    }
	    std::cout << "\n";
	 }

	 if (words.size() > 5) {
	    if (words[0] == "RES") {
	       try { 
		  double pos_x = lig_build::string_to_float(words[1]);
		  double pos_y = lig_build::string_to_float(words[2]);
		  double pos_z = lig_build::string_to_float(words[3]);
		  std::string res_type = words[4];
		  std::string label = words[5];
		  residue_circle_t rc(pos_x, pos_y, pos_z, res_type, label);
		  clipper::Coord_orth cp(pos_x, pos_y, pos_z);
		  lig_build::pos_t pos = mol.input_coords_to_canvas_coords(cp);
		  rc.set_canvas_pos(pos);
		  v.push_back(residue_circle_t(rc));
	       }
	       catch (std::runtime_error rte) {
		  std::cout << "failed to parse :" << lines[i] << ":" << std::endl;
	       }
	    }
	 }
	 if (words.size() == 6) {
	    if (words[0] == "BOND") { // written with space
	       try {
		  double bond_l = lig_build::string_to_float(words[3]);
		  std::string atom_name = lines[i].substr(5,4);
		  int bond_type = lig_build::string_to_int(words[5]);
		  lbg_info_t::bond_to_ligand_t btl(atom_name, bond_l);
		  btl.bond_type = bond_type;
		  std::cout << "adding bond " << v.size() << " to :"
			    << atom_name << ": " << bond_l << std::endl;
		  if (v.size())
		     v.back().add_bond_to_ligand(btl);
	       }
	       catch (std::runtime_error rte) {
		  std::cout << "failed to parse :" << lines[i] << ":" << std::endl;
	       }
	    }
	 } 
      }
   }
   std::cout << "found " << v.size() << " residue centres" << std::endl;
   return v;
} 

// "must take exactly one argument" problem
// 
// std::ostream &
// lbg_info_t::operator<<(std::ostream &s, residue_circle_t rc) {

//    s << "res-circ{" << rc.pos << " " << rc.label << " with bond_to_ligand length "
//      << bond_to_ligand.bond_length << "}";

// }


void
lbg_info_t::draw_residue_circles(const std::vector<residue_circle_t> &residue_circles) {

   for (unsigned int i=0; i<residue_circles.size(); i++) {
      lig_build::pos_t pos = residue_circles[i].pos;
      add_residue_circle(residue_circles[i], pos);
   }
}

void
lbg_info_t::add_residue_circle(const residue_circle_t &residue_circle,
			       const lig_build::pos_t &pos) {

   if (0)
      std::cout << "   adding cirles " << residue_circle.residue_type
		<< " at init pos " << pos << " and canvas_drag_offset "
		<< canvas_drag_offset << std::endl;

   lig_build::pos_t circle_pos = pos;
      
   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   GooCanvasItem *group = goo_canvas_group_new (root, "stroke-color", dark,
						NULL);

   // Capitalise the residue type (takes less space than upper case).
   std::string rt = residue_circle.residue_type.substr(0,1);
   rt += coot::util::downcase(residue_circle.residue_type.substr(1));

   // fill colour and stroke colour
   std::pair<std::string, std::string> col = get_residue_circle_colour(residue_circle.residue_type);
   double line_width = 1.0;
   if (col.second != dark)
      line_width = 3.0;

   if (col.first != "") {
      GooCanvasItem *cirle = goo_canvas_ellipse_new(root,
						    circle_pos.x, circle_pos.y,
						    16.0, 16.0,
						    "line_width", line_width,
						    "fill-color",   col.first.c_str(),
						    "stroke-color", col.second.c_str(),
						    NULL);
   } else {
      GooCanvasItem *cirle = goo_canvas_ellipse_new(root,
						    circle_pos.x, circle_pos.y,
						    16.0, 16.0,
						    "line_width", line_width,
						    "stroke-color", col.second.c_str(),
						    NULL);
   }
   
   GooCanvasItem *text_1 = goo_canvas_text_new(root, rt.c_str(),
					       circle_pos.x, circle_pos.y-5, -1,
					       GTK_ANCHOR_CENTER,
					       "font", "Sans 9",
					       "fill_color", dark,
					       NULL);

   GooCanvasItem *text_2 = goo_canvas_text_new(root, residue_circle.residue_label.c_str(),
					       circle_pos.x, circle_pos.y+5, -1,
					       GTK_ANCHOR_CENTER,
					       "font", "Sans 7",
					       "fill_color", dark,
					       NULL);
}

// Return the fill colour and the stroke colour.
// 
std::pair<std::string, std::string>
lbg_info_t::get_residue_circle_colour(const std::string &residue_type) const {

   std::string fill_colour = "";
   std::string stroke_colour = dark;

   std::string grease = "#bbffbb";
   std::string purple = "#eeccee";
   std::string red    = "#cc0000";
   std::string blue   = "#0000cc";

   if (residue_type == "ALA")
      fill_colour = grease;
   if (residue_type == "TRP")
      fill_colour = grease;
   if (residue_type == "PHE")
      fill_colour = grease;
   if (residue_type == "LEU")
      fill_colour = grease;
   if (residue_type == "PRO")
      fill_colour = grease;
   if (residue_type == "ILE")
      fill_colour = grease;
   if (residue_type == "VAL")
      fill_colour = grease;

   if (residue_type == "GLY")
      fill_colour = purple;
   if (residue_type == "ASP")
      fill_colour = purple;
   if (residue_type == "ASN")
      fill_colour = purple;
   if (residue_type == "CYS")
      fill_colour = purple;
   if (residue_type == "GLN")
      fill_colour = purple;
   if (residue_type == "GLU")
      fill_colour = purple;
   if (residue_type == "HIS")
      fill_colour = purple;
   if (residue_type == "LYS")
      fill_colour = purple;
   if (residue_type == "LYS")
      fill_colour = purple;
   if (residue_type == "MET")
      fill_colour = purple;
   if (residue_type == "MSE")
      fill_colour = purple;
   if (residue_type == "ARG")
      fill_colour = purple;
   if (residue_type == "SER")
      fill_colour = purple;
   if (residue_type == "THR")
      fill_colour = purple;
   if (residue_type == "TYR")
      fill_colour = purple;
   if (residue_type == "HOH")
      fill_colour = "white";

   if (residue_type == "ASP")
      stroke_colour = red;
   if (residue_type == "GLU")
      stroke_colour = red;
   if (residue_type == "LYS")
      stroke_colour = blue;
   if (residue_type == "ARG")
      stroke_colour = blue;
   if (residue_type == "HIS")
      stroke_colour = blue;
	 
   return std::pair<std::string, std::string> (fill_colour, stroke_colour);

} 


std::vector<solvent_accessible_atom_t>
lbg_info_t::read_solvent_accessibilities(const std::string &file_name) const {

   // return this
   std::vector<solvent_accessible_atom_t> solvent_accessible_atoms;
   
   std::ifstream f(file_name.c_str());
   if (!f) {
      std::cout << "Failed to open " << file_name << std::endl;
   } else {

      std::cout << "reading solvent accessibilites file: " << file_name << std::endl;
      
      std::vector<std::string> lines;
      std::string line;
      while (std::getline(f, line)) { 
	 lines.push_back(line);
      }

      for (unsigned int i=0; i<lines.size(); i++) {
	 std::vector<std::string> words = coot::util::split_string_no_blanks(lines[i], " ");
	 if (words.size() > 5) {
	    
	    if (words[0] == "ATOM:") {
	       std::string atom_name = lines[i].substr(5,4);
	       try {
		  double pos_x = lig_build::string_to_float(words[2]);
		  double pos_y = lig_build::string_to_float(words[3]);
		  double pos_z = lig_build::string_to_float(words[4]);
		  double sa    = lig_build::string_to_float(words[5]);
		  clipper::Coord_orth pt(pos_x, pos_y, pos_z);
		  
		  if (0) 
		     std::cout << "got atom name :" << atom_name << ": and pos "
			       << pt.format() << " and accessibility: " << sa
			       << std::endl;
		  solvent_accessible_atom_t saa(atom_name, pt, sa);
		  solvent_accessible_atoms.push_back(saa);
	       }
	       catch (std::runtime_error rte) {
		  std::cout << "failed to parse :" << lines[i] << ":" << std::endl;
	       }
	    }
	 }
      }
   }
   return solvent_accessible_atoms;
}

void
lbg_info_t::draw_solvent_accessibility_of_atom(const lig_build::pos_t &pos, double sa,
					       GooCanvasItem *root) {

   int n_circles = int(sa*40) + 1; // needs fiddling?
   if (n_circles> 10) n_circles = 10; // needs fiddling?

   for (unsigned int i=0; i<n_circles; i++) { 
      double rad = 3.0 * double(i+1); // needs fiddling?
      GooCanvasItem *cirle = goo_canvas_ellipse_new(root,
						    pos.x, pos.y,
						    rad, rad,
						    "line_width", 0.0,
						    "fill-color-rgba", 0x7755bb30,
						    NULL);
   }
}

void
lbg_info_t::draw_bonds_to_ligand() {
   
   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   for (unsigned int ic=0; ic<residue_circles.size(); ic++) { 
      if (residue_circles[ic].bonds_to_ligand.size()) {
	 for (unsigned int ib=0; ib<residue_circles[ic].bonds_to_ligand.size(); ib++) { 

	    if (0) 
	       std::cout << "====================== yay bond to ligand from "
			 << residue_circles[ic].residue_label << " "
			 << residue_circles[ic].residue_type <<  " to ligand atom " 
			 << residue_circles[ic].bonds_to_ligand[ib].ligand_atom_name
			 << " ====== " << std::endl;
	    
	    try {

	       lig_build::pos_t pos = residue_circles[ic].pos;
	       std::string at_name = residue_circles[ic].bonds_to_ligand[ib].ligand_atom_name;
	       lig_build::pos_t lig_at_pos = mol.get_atom_canvas_position(at_name);
	       lig_build::pos_t rc_to_lig_at = lig_at_pos - pos;
	       lig_build::pos_t rc_to_lig_at_uv = rc_to_lig_at.unit_vector();
	       lig_build::pos_t B = lig_at_pos - rc_to_lig_at_uv * 8;
	       lig_build::pos_t A = pos + rc_to_lig_at_uv * 20;

	       // some colours
	       std::string blue = "blue";
	       std::string green = "darkgreen";
	       std::string olive = "#666600";
	       std::string stroke_colour = dark; // unset

	       // arrows (acceptor/donor) and stroke colour (depending
	       // on mainchain or sidechain interaction)
	       // 
	       gboolean start_arrow = 0;
	       gboolean   end_arrow = 0;
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_DONOR_SIDECHAIN) {
		  end_arrow = 1;
		  stroke_colour = green;
	       }
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_DONOR_MAINCHAIN) {
		  end_arrow = 1;
		  stroke_colour = blue;
	       }
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_ACCEPTOR_SIDECHAIN) {
		  start_arrow = 1;
		  stroke_colour = green;
	       }
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_ACCEPTOR_MAINCHAIN) {
		  start_arrow = 1;
		  stroke_colour = blue;
	       }
	       if (residue_circles[ic].residue_type == "HOH") { 
		  stroke_colour = olive;
		  start_arrow = 0;
		  end_arrow = 0;
	       }
	    
	       GooCanvasLineDash *dash = goo_canvas_line_dash_new (2, 2.5, 2.5);
	    
	       GooCanvasItem *item = goo_canvas_polyline_new_line(root,
								  A.x, A.y,
								  B.x, B.y,
								  "line-width", 2.5,
								  "line-dash", dash,
 								  "start_arrow", start_arrow,
 								  "end_arrow",   end_arrow,
								  "stroke-color", stroke_colour.c_str(),
								  NULL);
	    }
	    catch (std::runtime_error rte) {
	       std::cout << "WARNING:: " << rte.what() << std::endl;
	    }
	 }
	 
      } else {
	 if (0) 
	    std::cout << "... no bond to ligand from residue circle "
		      << residue_circles[ic].residue_label << " "
		      << residue_circles[ic].residue_type << std::endl;
      }
   }
}
