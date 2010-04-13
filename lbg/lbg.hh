/* lbg/lbg.hh
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


#ifndef LBG_HH
#define LBG_HH

#include <iostream>
#include <map>
#include <queue>

#include <gtk/gtk.h>
#include <goocanvas.h>

#include "mmdb_manager.h"
#include "mmdb_sbase.h"
#define MONOMER_DIR_STR "COOT_SBASE_DIR"

#include "lig-build.hh"
#include "lbg-molfile.hh"

#include "some-coot-utils.hh"

#define dark "#111111"

static double SINGLE_BOND_CANVAS_LENGTH=35.27;

bool save_togglebutton_widgets(GtkBuilder *builder);

void lbg_handle_toggle_button(GtkToggleToolButton *tb, GtkWidget *canvas, int mode);
GtkWidget *get_canvas_from_scrolled_win(GtkWidget *scrolled_window);

// Unfortunately, the bonds and atoms must also have canvas items as
// part of their make-up.  Why?
//
// Consider changing a carbon to a nitrogen (or vica vera).  We want
// the molecular description of the *bond* to stay the same, we want
// the atom to change its element (and the widget representing it to
// be changed or generated or deleted (when going back to Carbon)) and
// the representation of the bond should be changed from a line to the
// atom point to a line that approaches (but does not touch) the atom
// point.  Which means that the bond and its representation are at
// different places.
//
// Now try to delete that bond. Which is the widget that needs to be
// removed from the canvas?
//
// We can only know that if the canvas item is part of the bond
// description.
// 

// ====================================================================
//                     widgeted_atom_t
// ====================================================================

class widgeted_atom_t : public lig_build::atom_t {
   std::string font_colour;
   GooCanvasItem *ci;
   void clear(GooCanvasItem *root) {
      gint child_index = goo_canvas_item_find_child(root, ci);
      if (child_index != -1) {
	 goo_canvas_item_remove_child(root, child_index);
      }
      ci = NULL;
   }

   // then general form, use this not the below 2 (which are used by
   // this function)
   // 
   GooCanvasItem *make_canvas_text_item(const std::string atom_name_in,
					const std::string &fc,
					GooCanvasItem *root) {
      GooCanvasItem *text_item = NULL;
      if (atom_name_in != "C") {
	 if (atom_name_in.length() > 2) {
	    text_item = make_subscripted_canvas_item(atom_name_in, fc, root);
	 } else { 
	    text_item = make_convention_canvas_item(atom_name_in, fc, root);
	 }
      }
      return  text_item;
   }
   
   GooCanvasItem *make_convention_canvas_item(const std::string atom_name_in,
					       const std::string &fc,
					       GooCanvasItem *root) const {

      GooCanvasItem *text_item = goo_canvas_text_new(root, atom_name_in.c_str(),
						     atom_position.x, atom_position.y, -1,
						     GTK_ANCHOR_CENTER,
						     "font", "Sans 10",
						     "fill_color", fc.c_str(),
						     NULL);
      return text_item;
   }
   
   GooCanvasItem *make_subscripted_canvas_item(const std::string atom_name_in,
					       const std::string &fc,
					       GooCanvasItem *root) const {

      GooCanvasItem *group = goo_canvas_group_new (root,
						   "fill_color", fc.c_str(),
						   NULL);
      std::string p1 = atom_name_in.substr(0,2);
      std::string p2 = atom_name_in.substr(2);
      
      GooCanvasItem *item_1 = goo_canvas_text_new(group, p1.c_str(),
						  atom_position.x, atom_position.y, -1,
						  GTK_ANCHOR_CENTER,
						  "font", "Sans 10",
						  "fill_color", fc.c_str(),
						  NULL);

      double pos_p2_x = atom_position.x + 13;
      double pos_p2_y = atom_position.y + 2;
      
      GooCanvasItem *item_2 = goo_canvas_text_new(group, p2.c_str(),
						  pos_p2_x, pos_p2_y, -1,
						  GTK_ANCHOR_CENTER,
						  "font", "Sans 8",
						  "fill_color", fc.c_str(),
						  NULL);
      return group;
   }

public:
   widgeted_atom_t(lig_build::atom_t &at_in, GooCanvasItem *ci_in) : lig_build::atom_t(at_in) {
      ci = ci_in;
      font_colour = "yellow";
   }
   widgeted_atom_t(lig_build::pos_t pos_in,
		   std::string ele_in,
		   int charge_in,
		   GooCanvasItem *ci_in) :
      lig_build::atom_t(pos_in, ele_in, charge_in) {
      ci = ci_in;
      font_colour = "hotpink";
   }
   GooCanvasItem *get_canvas_item() const { return ci; }
   void update_canvas_item(GooCanvasItem *new_item, GooCanvasItem *root) {
      clear(root);
      ci = new_item;
   }
   bool update_name_maybe(const std::string &atom_name_in,
			  GooCanvasItem *root) {
      return update_name_maybe(atom_name_in, font_colour, root);
   }
   bool update_name_maybe(const std::string &atom_name_in,
			  const std::string &fc,
			  GooCanvasItem *root) {
      bool changed_status = 0;
      font_colour = fc;
      GooCanvasItem *text_item = NULL;
      if (! is_closed()) { 
	 if (atom_name_in != get_atom_name()) {
	    // std::cout << "in update_name_maybe: from " << get_atom_name() << " to "
	    // << atom_name_in << std::endl;
	    changed_status = set_atom_name(atom_name_in);
	    // std::cout << "in update_name_maybe: changed_status " << changed_status
	    // << std::endl;
	    if (changed_status) {
	       text_item = make_canvas_text_item(atom_name_in, fc, root);
	       update_canvas_item(text_item, root);
	    }
	 }
      } else {
	 update_canvas_item(text_item, root); // close atom, replace with null.
      } 
      return changed_status;
   }
   void update_name_forced(const std::string &atom_name_in,
			   const std::string &fc, 
			   GooCanvasItem *root) {
      set_atom_name(atom_name_in);
      GooCanvasItem *text_item = make_canvas_text_item(atom_name_in, fc, root);
      update_canvas_item(text_item, root);
   }
   
   void update_name_forced(const std::string &atom_name_in,
			   GooCanvasItem *root) {
      update_name_forced(atom_name_in, font_colour, root);
   }
   void close(GooCanvasItem *root) {
      // std::cout << " closing subclass atom" << std::endl;
      lig_build::atom_t::close();
      update_canvas_item(NULL, root);
   }
};

// ====================================================================
//                     widgeted_bond_t
// ====================================================================

class widgeted_bond_t : public lig_build::bond_t {
   GooCanvasItem *ci;
   void clear(GooCanvasItem *root) {
      gint child_index = goo_canvas_item_find_child(root, ci);
      if (child_index != -1) {
	 goo_canvas_item_remove_child(root, child_index);
      }
      ci = NULL;
   }

   void construct_internal(const lig_build::atom_t &atom_first,
			   const lig_build::atom_t &atom_second,
			   bond_type_t bt, GooCanvasItem *root) {
      bool shorten_first = 0;
      bool shorten_second = 0;
      if (atom_first.name != "C") { 
	 shorten_first = 1;
      } 
      if (atom_second.name != "C") { 
	 shorten_second = 1;
      }
      lig_build::pos_t pos_1 =  atom_first.atom_position;
      lig_build::pos_t pos_2 = atom_second.atom_position;
      ci = canvas_item_for_bond(pos_1, pos_2, shorten_first, shorten_second, bt, root);
   }

   // all bonds are made this way...
   // 
   GooCanvasItem *canvas_item_for_bond(const lig_build::pos_t &pos_1,
				       const lig_build::pos_t &pos_2,
				       bool shorten_first,
				       bool shorten_second,
				       bond_type_t bt,
				       GooCanvasItem *root) const;

   // We need to make a shorter bond canvas line because we have (say)
   // changed a carbon to a N (bond canvas line now does not
   // completely extend to the atom position).
   // 
   void make_new_canvas_item_given_type(const lig_build::atom_t &atom_changed,
					const lig_build::atom_t &atom_other,
					lig_build::bond_t::bond_type_t bt,
					GooCanvasItem *root) {

      lig_build::pos_t A = atom_changed.atom_position;
      lig_build::pos_t B =   atom_other.atom_position;

      bool shorten_first = 0;
      bool shorten_second = 0;
      if (atom_changed.name != "C")
	 shorten_first = 1;
      if (atom_other.name != "C")
	 shorten_second = 1;
      GooCanvasItem *new_line = canvas_item_for_bond(A, B, shorten_first, shorten_second,
						     bt, root);
      update_canvas_item(new_line, root);
   }

   GooCanvasItem * canvas_item_double_bond(const lig_build::pos_t &pos_1,
					   const lig_build::pos_t &pos_2,
					   GooCanvasItem *root) const;
					
   GooCanvasItem * canvas_item_double_aromatic_bond(const lig_build::pos_t &pos_1,
						    const lig_build::pos_t &pos_2,
						    GooCanvasItem *root) const;
					
   GooCanvasItem * make_wedge_bond_item(const lig_build::pos_t &pos_1,
					const lig_build::pos_t &pos_2,
					const lig_build::bond_t::bond_type_t &bt,
					GooCanvasItem *root) const;
   GooCanvasItem * make_wedge_out_bond_item(const lig_build::pos_t &pos_1,
					    const lig_build::pos_t &pos_2,
					    GooCanvasItem *root) const;
   GooCanvasItem * make_wedge_in_bond_item(const lig_build::pos_t &pos_1,
					   const lig_build::pos_t &pos_2,
					   GooCanvasItem *root) const;


   // -------------------- widgeted_bond_t public --------------------------

public:

   // this is for widgeted_bond_t that are invalid (to be assigned later).
   widgeted_bond_t() : lig_build::bond_t() {
      ci = NULL;
   } 
   
   // Now we use a constructor that does the creation of the canvas item too
   //
   widgeted_bond_t(int first, int second, 
		   const lig_build::atom_t &atom_first, const lig_build::atom_t &atom_second,
		   bond_type_t bt, GooCanvasItem *root) :
      lig_build::bond_t(first, second, bt) {
      construct_internal(atom_first, atom_second, bt, root);
   }
   // as above, but we give the centre of the ring too.
   widgeted_bond_t(int first, int second, 
		   const lig_build::atom_t &atom_first, const lig_build::atom_t &atom_second,
		   lig_build::pos_t centre_pos_in,
		   bond_type_t bt, GooCanvasItem *root) :
      bond_t(first, second, centre_pos_in, bt) {
      construct_internal(atom_first, atom_second, bt, root);
   }
   
   // old constructor.  Consider deleting.
   // 
   // widgeted_bond_t(int first, int second, bond_type_t bt, GooCanvasItem *ci_in) :
   // lig_build::bond_t(first, second, bt) {
   // ci = ci_in;
   // is_closed_ = 0;
   // }

   // old constructor.  Consider deleting.
   // 
   // widgeted_bond_t(int first, int second,
   // lig_build::pos_t centre_pos_in,
   // lig_build::bond_t::bond_type_t bt,
   // GooCanvasItem *ci_in) :
   // lig_build::bond_t(first, second, centre_pos_in, bt) {
   // ci = ci_in;
   // is_closed_ = 0;
   // }


   void update_canvas_item(GooCanvasItem *new_item, GooCanvasItem *root) {
      clear(root);
      ci = new_item;
   }

   // We need to make a shorter bond canvas line because we have (say)
   // changed a carbon to a N (bond canvas line now does not
   // completely extend to the atom position).
   void make_new_canvas_item(const lig_build::atom_t &atom_changed,
			     const lig_build::atom_t &atom_other,
			     GooCanvasItem *root) {

      lig_build::bond_t::bond_type_t bt = get_bond_type();
      make_new_canvas_item_given_type(atom_changed, atom_other, bt, root);
   }
   void change_bond_order(const lig_build::atom_t &atom_changed,
			  const lig_build::atom_t &atom_other,
			  GooCanvasItem *root) {
      change_bond_order(atom_changed, atom_other, 0, root);
   } 
   void change_bond_order(const lig_build::atom_t &atom_changed,
			  const lig_build::atom_t &atom_other,
			  bool allow_triple_toggle,
			  GooCanvasItem *root) {
      lig_build:: atom_t at_1 = atom_changed;
      lig_build:: atom_t at_2 = atom_other;
      lig_build::bond_t::bond_type_t bt = get_bond_type();
      if (bt == lig_build::bond_t::SINGLE_BOND) { 
	 if (allow_triple_toggle)
	    bt = lig_build::bond_t::TRIPLE_BOND;
	 else 
	    bt = lig_build::bond_t::DOUBLE_BOND;
      } else { 
	 if (bt == lig_build::bond_t::DOUBLE_BOND)
	    if (allow_triple_toggle)
	       bt = lig_build::bond_t::TRIPLE_BOND;
	    else
	       bt = lig_build::bond_t::SINGLE_BOND;
	 else 
	    if (bt == lig_build::bond_t::TRIPLE_BOND)
	       bt = lig_build::bond_t::DOUBLE_BOND;
	    else
	       if (bt == lig_build::bond_t::IN_BOND) {
		  // std::cout << " add a reverse direction here " << std::endl;
		  std::swap(atom_1, atom_2);
		  std::swap(at_1, at_2);
		  bt = lig_build::bond_t::OUT_BOND;
	       } else {
		  if (bt == lig_build::bond_t::OUT_BOND) { 
		     bt = lig_build::bond_t::IN_BOND;
		  }
	       }
      }
      
//       std::cout << "changing bond type from " << get_bond_type() << " to "
// 		<< bt << std::endl;
      set_bond_type(bt);
      make_new_canvas_item_given_type(at_1, at_2, bt, root);
   }
   void close(GooCanvasItem *root) {
      // std::cout << " closing sub-class bond" << std::endl;
      lig_build::bond_t::close();
      update_canvas_item(NULL, root);
   }
   int mmdb_bond_type() const {
      int mmdb_bt = 1;
      switch (get_bond_type()) { 
      case SINGLE_BOND:
	 mmdb_bt = 1;
	 break;
      case IN_BOND:
	 mmdb_bt = 1;
	 break;
      case OUT_BOND:
	 mmdb_bt = 1;
	 break;
      case DOUBLE_BOND:
	 mmdb_bt = 2;
	 break;
      case TRIPLE_BOND:
	 mmdb_bt = 3;
	 break;
      case BOND_UNDEFINED:
	 mmdb_bt = UNASSIGNED_INDEX;
	 break;
      }
      return mmdb_bt;
   }
   void add_centre(const lig_build::pos_t &centre_in) {
      set_centre_pos(centre_in);
   }

   int get_other_index(int atom_index) const {
      int idx = get_atom_1_index();
      if (idx == atom_index)
	 idx = get_atom_2_index();
      return idx;
   }

};


// ====================================================================
//                     widgeted_molecule_t
// ====================================================================

#define MAX_SEARCH_DEPTH 9

class widgeted_molecule_t : public lig_build::molecule_t<widgeted_atom_t, widgeted_bond_t> {
   
private:
   std::string group;
   void init() {
      mol_in_max_y = 0;
      mol_in_min_y = 0;
      scale_correction.first = 0;
      scale_correction.second = 1;
   }
   bool member(const int &ind, const std::vector<int> &no_pass_atoms) const {
      bool found = 0;
      for (unsigned int i=0; i<no_pass_atoms.size(); i++) { 
	 if (no_pass_atoms[i] == ind) {
	    found = 1;
	    break;
	 }
      }
      return found;
   } 
   // Return a vector of bonds.  If empty, then it didn't find self.
   // 
   std::pair<bool, std::vector<int> >
   found_self_through_bonds(int atom_index_start, int atom_index_other) const;
   std::pair<bool, std::vector<int> >
   find_bonded_atoms_with_no_pass(int atom_index_start,
				  int atom_index_other, // must pass through this
				  int this_atom_index,
				  const std::vector<int> &no_pass_atoms,
				  int depth) const;
   void debug_pass_atoms(int atom_index, int this_atom_index, 
			 int depth,  const std::vector<int> &local_no_pass_atoms) const;
   std::pair<bool, double>
   get_scale_correction(const lig_build::molfile_molecule_t &mol_in) const;
   int get_number_of_atom_including_hydrogens() const;


public:
   widgeted_molecule_t() { init(); }
   widgeted_molecule_t(const lig_build::molfile_molecule_t &mol_in);

   // return 0 as first if not highlighting a bond
   std::pair<bool, widgeted_bond_t> highlighted_bond_p(int x, int y) const;

   // return -1 as the atom index if not highlighting an atom.
   std::pair<int, widgeted_atom_t> highlighted_atom_p(int x, int y) const;
   bool write_mdl_molfile(const std::string &file_name) const;
   bool write_minimal_cif_file(const std::string &file_name) const;
   bool close_bond(int ib, GooCanvasItem *root, bool handle_post_delete_stray_atoms_flag);
   bool close_atom(int iat, GooCanvasItem *root);
   std::vector<int> get_unconnected_atoms() const;

   // don't count closed bonds.
   std::vector<int> bonds_having_atom_with_atom_index(int test_atom_index) const;
   bool operator==(const widgeted_molecule_t &mol_other) const;
   int n_stray_atoms() const; // unbonded atoms
   std::vector<int> stray_atoms() const;
   void translate(const lig_build::pos_t &delta); // move the atoms
   lig_build::pos_t input_coords_to_canvas_coords(const clipper::Coord_orth &in) const;
      

   // make private when bug is fixed.
   lig_build::pos_t centre_correction;
   std::pair<bool, double> scale_correction;
   double mol_in_min_y;
   double mol_in_max_y;
   
};


// ====================================================================
//                     lbg_info_t
// ====================================================================

class lbg_info_t {

public:
   class highlight_data_t {
      int n_atoms_;
      lig_build::pos_t pos_1_;
      lig_build::pos_t pos_2_;
      std::pair<int,int> bond_indices; // the atom indices of the bond.
      int atom_index; // the index of the atom in the molecule (single atom highlighting)
      bool has_ring_centre_flag;
      lig_build::pos_t ring_centre;
      GooCanvasItem *highlight_widget;
      lig_build::polygon_position_info_t get_new_polygon_centre_using_2_atoms(int n_edges,
								const double &radius) const;
      lig_build::polygon_position_info_t
      get_new_polygon_centre_using_1_atom(int n_edges,
					  bool spiro,
					  const double &radius_std,
					  const double &radius_corr,
					  const widgeted_molecule_t &mol) const;
   public:
      highlight_data_t(GooCanvasItem *w_in,
		       std::pair<int, int> bond_indices_in,
		       const lig_build::pos_t &p1,
		       const lig_build::pos_t &p2) {
	 n_atoms_ = 2;
	 pos_1_ = p1;
	 pos_2_ = p2;
	 highlight_widget = w_in;
	 atom_index = UNASSIGNED_INDEX; // no value
	 bond_indices = bond_indices_in;
      }
      highlight_data_t(GooCanvasItem *w_in,
		       const lig_build::pos_t &p, int index_in) {
	 n_atoms_ = 1;
	 pos_1_ = p;
	 highlight_widget = w_in;
	 atom_index = index_in;
	 bond_indices = std::pair<int, int> (UNASSIGNED_INDEX, UNASSIGNED_INDEX);
      }

      highlight_data_t() {
	 highlight_widget = NULL;
	 n_atoms_ = 0;
	 atom_index = UNASSIGNED_INDEX; // unset
	 bond_indices = std::pair<int, int> (UNASSIGNED_INDEX, UNASSIGNED_INDEX);
      }
      bool has_contents() const {
	 if (highlight_widget)
	    return 1;
	 else
	    return 0;
      }
      bool single_atom() const {
	 return (n_atoms_ == 1);
      }
      int get_atom_index() const { return atom_index; }
      lig_build::pos_t get_atom_1_pos() const {
	 return pos_1_;
      } 
      void clear(GooCanvasItem *root) {
	 n_atoms_ = 0;
	 if (highlight_widget) {
	    gint child_index = goo_canvas_item_find_child(root, highlight_widget);
	    if (child_index != -1) {
	       goo_canvas_item_remove_child(root, child_index);
	       highlight_widget = NULL;
	    }
	 } else {
	    std::cout << "in clear() NULL highlight_widget" << std::endl;
	 }
      }
      lig_build::polygon_position_info_t
      get_new_polygon_centre(int n_edges,
			     bool spiro_flag,
			     const double &radius_standard,
			     const double &radius_corrected,
			     const widgeted_molecule_t &mol) const;
      std::pair<int, int> get_bond_indices() const { return bond_indices; }
   }; // finish highlight_data_t class

   class match_results_t {
   public:
      bool success;
      std::string name;
      std::string comp_id;
      CResidue *res;
      // clipper::RTop_orth
      match_results_t(const std::string &comp_id_in, const std::string &name_in, CResidue *res_in) {
	 name = name_in;
	 comp_id = comp_id_in;
	 res = res_in;
	 if (res_in)
	    success = 1;
	 else
	    success = 0;
      }
   };

   class residue_circle_t {
   public:
      double pos_x;
      double pos_y;
      double pos_z;
      std::string residue_type;
      std::string residue_label;
      residue_circle_t(const double &x_in, const double &y_in, const double &z_in,
		       const std::string &type_in,
		       const std::string &label_in) {
	 pos_x = x_in;
	 pos_y = y_in;
	 pos_z = z_in;
	 residue_type = type_in;
	 residue_label = label_in;
      }
   };
			  
   
private:
   bool try_stamp_polygon(int n_edges, int x_pos_centre, int y_pos_centre,
			  bool is_spiro, bool is_aromatic);
   void stamp_polygon_anywhere(int n_edges, int x_pos_centre, int y_pos_centre,
			       bool is_aromatic, GooCanvasItem *root);
   std::vector<int> stamp_polygon(int n_edges, lig_build::polygon_position_info_t ppi,
				  bool aromatic_flag, GooCanvasItem *root);
   std::vector<int> try_stamp_polygon_using_highlighted_data(int n_edges,
							     bool spiro_flag,
							     bool aromomatic_flag,
							     GooCanvasItem *root);
   bool try_stamp_hexagon_aromatic(int x_pos_centre, int y_pos_centre, bool shift_is_pressed);
   std::vector<widgeted_molecule_t> previous_molecules;
   int save_molecule_index;
   bool in_delete_mode_;
   highlight_data_t highlight_data;
   bool is_atom_element(int addition_mode) const;
   bool is_bond(int addition_mode) const;
   bool try_change_to_element(int addition_element_mode); // check for highlighted atom;
   bool try_add_or_modify_bond(int canvas_addition_mode, int x, int y); //  ditto.
   bool add_bond_to_atom(int atom_index, int canvas_addition_mode);
   void add_bond_to_atom_with_0_neighbours(int atom_index, int canvas_addition_mode);
   void add_bond_to_atom_with_1_neighbour(int atom_index, int canvas_addition_mode,
					  int bond_index);
   void add_bond_to_atom_with_2_neighbours(int atom_index, int canvas_addition_mode,
					   const std::vector<int> &bond_indices);
   void add_bond_to_atom_with_3_neighbours(int atom_index, int canvas_addition_mode,
					   const std::vector<int> &bond_indices);
   std::string to_element(int addition_mode) const;
   std::string font_colour(int addition_element_mode) const;
   std::string font_colour(const std::string &ele) const;
   lig_build::bond_t::bond_type_t addition_mode_to_bond_type(int canvas_addition_mode) const;
   void try_stamp_bond_anywhere(int canvas_addition_mode, int x_mouse, int y_mouse); // always modifies.
   bool change_atom_element(int atom_index, std::string new_element, std::string fc);
   void change_atom_name_maybe(int atom_index);
   lig_build::pos_t mouse_at_click;
   void save_molecule();
   std::string mdl_file_name; // for save function.
   void add_search_combobox_text() const;
   bool make_saves_mutex;
   double canvas_scale; 
   void init_internal() {
      in_delete_mode_ = 0;
      save_molecule_index = UNASSIGNED_INDEX;
      make_saves_mutex = 1; // allow saves initially
      search_similarity = 0.93;
      coot_mdl_time = 0;
      canvas_scale = 1.0;
      canvas_drag_offset = lig_build::pos_t(0,0);
   }
   
   // return a status and a vector of atoms (bonded to atom_index) having
   // only one bond.
   // 
   std::pair<bool, std::vector<int> > 
   have_2_stubs_attached_to_atom(int atom_index, const std::vector<int> &bond_indices) const;
   void squeeze_in_a_4th_bond(int atom_index, int canvas_addition_mode,
			      const std::vector<int> &bond_indices);
   std::vector<double>
   get_angles(int atom_index, const std::vector<int> &bond_indices) const;
   lig_build::pos_t  new_pos_by_bisection(int atom_index,
					  const std::vector<int> &bond_indices,
					  const std::vector<double> &angles,
					  GooCanvasItem *root) const;
   bool all_closed_rings(int atom_index, const std::vector<int> &bond_indices) const;
   std::vector<lig_build::pos_t>
   get_centres_from_bond_indices(const std::vector<int> &bond_indices) const;
   lig_build::pos_t get_new_pos_not_towards_ring_centres(int atom_index,
							 const std::vector<int> &bond_indices) const;

   // sbase functions
   CSBase *SBase;
   float search_similarity;
   int init_sbase(const std::string &sbase_monomer_dir_in);
   std::vector<match_results_t> compare_vs_sbase(CGraph *graph1,
						 float similarity,
						 int n_vertices) const;
   // return 0 on strangeness, to pass in search.
   // 

   int get_min_match(const int &n1) const {
      int most_1 = int (search_similarity * float(n1));
      return most_1;
   }

   // not used.
   int get_min_match(const int &n1, const int &n2) const {
      int r = 0;
      int most_1 = int (search_similarity * float(n1));
      int most_2 = int (search_similarity * float(n2));
      if ((n2>=most_1) && (n1>=most_2))  {
	 r = (most_2 > most_1) ? most_2 : most_1;
      }
      return r;
   }
   match_results_t residue_from_best_match(CGraph &graph_1, CGraph &graph_2,
					   CGraphMatch &match, int n_match, 
					   CSBStructure *SBS) const;
   void display_search_results(const std::vector<lbg_info_t::match_results_t> v) const;
   void orthogonalise_2_bonds(int atom_index,
			      const std::vector<int> &attached_bonds,
			      const std::vector<int> &bond_indices);

   void translate_molecule(const lig_build::pos_t &delta);
   void clear_canvas();
   void add_residue_circle(const residue_circle_t &rc, const lig_build::pos_t &p);

   std::pair<std::string, std::string>
   get_residue_circle_colour(const std::string &residue_type) const;
   lig_build::pos_t canvas_drag_offset;
   
public:
   lbg_info_t(GtkWidget *canvas_in) {
      canvas = canvas_in;
      init_internal();
   }
   lbg_info_t() { init_internal(); }
   // toggle button modes, mutually exclusive
   enum { NONE, TRIANLE, SQUARE, PENTAGON, HEXAGON, HEXAGON_AROMATIC, HEPTAGON, OCTAGON,
	  ATOM_C, ATOM_N, ATOM_O, ATOM_S, ATOM_P, ATOM_F, ATOM_CL, ATOM_I, ATOM_BR, ATOM_X,
	  CHARGE, ADD_SINGLE_BOND, ADD_DOUBLE_BOND, ADD_TRIPLE_BOND, ADD_STEREO_OUT_BOND,
	  DELETE_MODE};
   void init(GtkBuilder *builder);
   GtkWidget *about_dialog; 
   GtkWidget *search_combobox;
   GtkWidget *open_dialog;
   GtkWidget *save_as_dialog;
   GtkWidget *lbg_export_as_pdf_dialog;
   GtkWidget *lbg_export_as_png_dialog;
   GtkWidget *lbg_sbase_search_results_dialog;
   GtkWidget *lbg_sbase_search_results_vbox;
   GtkWidget *canvas;
   std::map<std::string, GtkToggleToolButton *> widget_names;
   widgeted_molecule_t mol;
   int canvas_addition_mode;
   bool save_togglebutton_widgets(GtkBuilder *builder);
   void handle_item_add(GdkEventButton *event);
   void handle_item_delete(GdkEventButton *event);
   void untoggle_others_except(GtkToggleToolButton *button_toggled_on);
   bool item_highlight_maybe(int x, int y);
   void highlight_bond(const lig_build::bond_t &bond, bool delete_mode);
   void highlight_atom(const lig_build::atom_t &atom, int atom_index, bool delete_mode);
   void remove_bond_and_atom_highlighting();
   void set_in_delete_mode(bool v) {
      in_delete_mode_ = v;
   }
   bool in_delete_mode_p() const { return in_delete_mode_; }
   double radius(int n_edges) const; // depends on zoom? (for future).
   void clear();
   std::string get_stroke_colour(int i, int n) const;
   void drag_canvas(int mouse_x, int mouse_y);
   void write_pdf(const std::string &file_name) const;
   void write_png(const std::string &file_name) const;
   void set_mouse_pos_at_click(int xpos, int ypos) {
      mouse_at_click = lig_build::pos_t(double(xpos), double(ypos));
   }
   void render_from_molecule(const widgeted_molecule_t &mol_in);
   void undo();
   void search() const;
   void import(const lig_build::molfile_molecule_t &mol_in, const std::string &filename);
   static void on_sbase_search_result_button_clicked(GtkButton *button, gpointer user_data);
   static gboolean watch_for_mdl_from_coot(gpointer user_data);
   time_t coot_mdl_time;
   void read_draw_residues(const std::string &file_name);
   std::vector<residue_circle_t> read_residues(const std::string &file_name) const;
   void draw_residue_circles(const std::vector<residue_circle_t> &v);
};

#endif // LBG_HH
