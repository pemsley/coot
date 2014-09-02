/* lbg/wmolecule.hh
 * 
 * Author: Paul Emsley
 * Copyright 2010, 2011 by The University of Oxford
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

#ifndef WMOLECULE_HH
#define WMOLECULE_HH

#include "lidia-core/lbg-shared.hh"

// #define dark "#111111"

class solvent_accessible_atom_t {
public:
   std::string atom_name;
   clipper::Coord_orth pt;
   double solvent_accessibility;
   std::vector<coot::bash_distance_t> bash_distances;
   solvent_accessible_atom_t(const std::string &at,
			     const clipper::Coord_orth &pt_in,
			     double sa) {
      atom_name = at;
      pt = pt_in;
      solvent_accessibility = sa;
   }
   solvent_accessible_atom_t() { }
   void add_bash_dist(double d) {
      bash_distances.push_back(coot::bash_distance_t(d));
   }
   void add_unlimited() {
      coot::bash_distance_t bd;      
      bash_distances.push_back(bd);
   }
};


class ligand_layout_graphic_primitives {
public:
   ligand_layout_graphic_primitives() {
      dark = "#111111";
   }
   std::string dark;
   int i;
   virtual GooCanvasItem *wrap_goo_canvas_group_new (GooCanvasItem *root,
					     const std::string &stroke_colour) const {
      // need fill-colour too?
      return goo_canvas_group_new(root, "stroke-color", stroke_colour.c_str(), NULL);
   }

   virtual GooCanvasItem *wrap_goo_canvas_text_new(GooCanvasItem *group,
					   const std::string &text,
					   double x_pos, double y_pos, 
					   int something,
					   GtkAnchorType anchor_type,
					   const std::string &font,
					   const std::string &fill_colour) const {
      return goo_canvas_text_new(group,
				 text.c_str(),
				 x_pos, y_pos, 
				 -1,
				 anchor_type,
				 "font", font.c_str(),
				 "fill_color", fill_colour.c_str(),
				 NULL);
   }
   virtual GooCanvasItem *wrap_goo_canvas_polyline_new_line(GooCanvasItem *root,
						    double pos_1_x, double pos_1_y,
						    double pos_2_x, double pos_2_y,
						    const std::string &key="stroke-color",
						    const std::string &value="#111111") const { 
      GooCanvasItem *item = 
 	 goo_canvas_polyline_new_line(root, 
 				      pos_1_x, pos_1_y,
 				      pos_2_x, pos_2_y,
 				      key.c_str(), value.c_str(),
 				      NULL);
      return item;
   }
   virtual GooCanvasItem *wrap_goo_canvas_polyline_new_line(GooCanvasItem *root,
							    double pos_1_x, double pos_1_y,
							    double pos_2_x, double pos_2_y,
							    double pos_3_x, double pos_3_y,
							    double pos_4_x, double pos_4_y,
							    const std::string &stroke_colour,
							    const std::string &fill_colour) const {
      
      GooCanvasItem *item = goo_canvas_polyline_new_line(root,
							 TRUE, 4,
							 pos_1_x, pos_1_y,
							 pos_2_x, pos_2_y,
							 pos_3_x, pos_3_y,
							 pos_4_x, pos_4_y,
							 "stroke-color", stroke_colour.c_str(),
							 "fill-color", fill_colour.c_str(),
							 NULL);
      return item;
   }

   virtual GooCanvasItem *
   wrap_goo_canvas_polyline_new(GooCanvasItem *root,
				double sharp_point_2_x, double sharp_point_2_y, 
				double sharp_point_1_x, double sharp_point_1_y, 
				double short_edge_pt_1_x, double short_edge_pt_1_y,
				double short_edge_pt_2_x, double short_edge_pt_2_y,
				std::string fc, std::string sc) const {

      GooCanvasItem *item = 
	 goo_canvas_polyline_new(root, 
				 TRUE, 4,
				 sharp_point_2_x, sharp_point_2_y, 
				 sharp_point_1_x, sharp_point_1_y, 
				 short_edge_pt_1_x, short_edge_pt_1_y,
				 short_edge_pt_2_x, short_edge_pt_2_y,
				 "fill-color", fc.c_str(),
				 "stroke-color", sc.c_str(),
				 NULL);
      return item;
   } 
   
   void wrap_goo_canvas_item_rotate(GooCanvasItem *ci,
					 double degrees, double cx, double cy) const {
      goo_canvas_item_rotate(ci, degrees, cx, cy);
   } 
};

// ====================================================================
//                     widgeted_atom_t
// ====================================================================

class widgeted_atom_t : public lig_build::atom_t , ligand_layout_graphic_primitives {
   std::string font_colour;
   double solvent_accessibility;
   GooCanvasItem *ci;
   // std::string atom_name; // typically names from a PDB file. 20111229 base class now
   void clear(GooCanvasItem *root) {
      gint child_index = goo_canvas_item_find_child(root, ci);
      if (child_index != -1) {
	 goo_canvas_item_remove_child(root, child_index);
      }
      ci = NULL;
   }
   void blank_text(GooCanvasItem *root) {
      gint child_index = goo_canvas_item_find_child(root, ci);
      if (child_index != -1) {
	 // set text to ""
	 g_object_set(ci, "visibility", GOO_CANVAS_ITEM_INVISIBLE, NULL);
	 // gtk_widget_hide(GTK_WIDGET(ci)); invalid cast
      } else {
	 // std::cout << "debug::       blank_text() ci not found in children of root "<< std::endl;
      } 
   }

   // then general form, use this not the below 2 (which are used by
   // this function)
   // 
   GooCanvasItem *make_canvas_text_item(const lig_build::atom_id_info_t &atom_id_info_in,
					const std::string &fc,
					GooCanvasItem *root);

public:

   widgeted_atom_t(lig_build::atom_t &at_in, GooCanvasItem *ci_in) : lig_build::atom_t(at_in) {
      ci = ci_in;
      font_colour = "yellow";
      solvent_accessibility = -1;
   }
   widgeted_atom_t(lig_build::pos_t pos_in,
		   std::string ele_in,
		   int charge_in,
		   GooCanvasItem *ci_in) :
      lig_build::atom_t(pos_in, ele_in, charge_in) {
      ci = ci_in;
      font_colour = "hotpink";
      solvent_accessibility = -1;
   }
   GooCanvasItem *get_canvas_item() const { return ci; }
   void update_canvas_item(GooCanvasItem *new_item, GooCanvasItem *root) {
      blank_text(root);
      ci = new_item;
   }
   bool update_atom_id_maybe(const lig_build::atom_id_info_t &atom_id_info_in,
			     GooCanvasItem *root) {
      return update_atom_id_maybe(atom_id_info_in, font_colour, root);
   }
   bool update_atom_id_maybe(const lig_build::atom_id_info_t &atom_id_info_in,
			     const std::string &fc,
			     GooCanvasItem *root) {
      bool changed_status = 0;
      font_colour = fc;
      GooCanvasItem *text_item = NULL;
      if (! is_closed()) {
	 if (atom_id_info_in.atom_id != get_atom_id()) {
	    changed_status = set_atom_id(atom_id_info_in.atom_id);
	    if (changed_status) {
	       text_item = make_canvas_text_item(atom_id_info_in, fc, root);
	       update_canvas_item(text_item, root);
	    }
	 }
      } else {
	 update_canvas_item(text_item, root); // close atom, replace with null.
      } 
      return changed_status;
   }
   void update_atom_id_forced(const lig_build::atom_id_info_t &atom_id_info_in,
			      const std::string &fc, 
			      GooCanvasItem *root) {
      set_atom_id(atom_id_info_in.atom_id);
      GooCanvasItem *text_item = make_canvas_text_item(atom_id_info_in, fc, root);
      update_canvas_item(text_item, root);
   }
   
   void update_atom_id_forced(const lig_build::atom_id_info_t &atom_id_info_in,
			      GooCanvasItem *root) {
      update_atom_id_forced(atom_id_info_in, font_colour, root);
   }
   void add_solvent_accessibility(double sa) {
      solvent_accessibility = sa;
   }
   double get_solvent_accessibility() const { return solvent_accessibility; } // negative for none.
   void close(GooCanvasItem *root) {
      // std::cout << " closing subclass atom" << std::endl;
      lig_build::atom_t::close();
      update_canvas_item(NULL, root);
   }
   void set_atom_name(const std::string atom_name_in) {
      atom_name = atom_name_in;
   }
//    std::string get_atom_name() const {
//       return atom_name;
//    } 
   std::vector<coot::bash_distance_t> bash_distances;
};

// ====================================================================
//                     widgeted_bond_t
// ====================================================================

static gboolean
on_wmolecule_key_press_event (GooCanvasItem *item,
			      GooCanvasItem *target,
			      GdkEventKey *event,
			      gpointer data)
{
  gchar *id = 0;
  // id = g_object_get_data (G_OBJECT (item), "id");
  g_print ("%s received key-press event\n", id ? id : "unknown");
  return FALSE;
}


class widgeted_bond_t : public lig_build::bond_t, ligand_layout_graphic_primitives {
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
      if (atom_first.atom_id != "C") { 
	 shorten_first = 1;
      } 
      if (atom_second.atom_id != "C") { 
	 shorten_second = 1;
      }
      lig_build::pos_t pos_1 =  atom_first.atom_position;
      lig_build::pos_t pos_2 = atom_second.atom_position;
      ci = canvas_item_for_bond(pos_1, pos_2, shorten_first, shorten_second, bt, root);

      // std::cout << "construct_internal() made ci " << ci << std::endl;
      
      g_signal_connect (ci, "key_press_event",
			G_CALLBACK (on_wmolecule_key_press_event), NULL);
      
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
      if (atom_changed.atom_id != "C")
	 shorten_first = 1;
      if (atom_other.atom_id != "C")
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
   
   void update_canvas_item(GooCanvasItem *new_item, GooCanvasItem *root) {
      clear(root);
      ci = new_item;
   }

   void rotate_canvas_item(gdouble cx, gdouble cy, gdouble degrees) {
      wrap_goo_canvas_item_rotate(ci, degrees, cx, cy);
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
      // std::cout << "change_bond_order " << atom_changed << " " << atom_other << std::endl;
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
      case AROMATIC_BOND:
      case DELOC_ONE_AND_A_HALF:
      case SINGLE_OR_DOUBLE:
      case SINGLE_OR_AROMATIC:
      case DOUBLE_OR_AROMATIC:
      case BOND_ANY:
      case BOND_UNDEFINED:
	 mmdb_bt = UNASSIGNED_INDEX;
	 break;
      }
      return mmdb_bt;
   }

}; // end of widgeted_bond_t

// trivial container for a (copy of an) atom an its ring centre (if
// it has one)
class widgeted_atom_ring_centre_info_t {
public:
   widgeted_atom_t atom;
   bool has_ring_centre_flag;
   lig_build::pos_t ring_centre;
   widgeted_atom_ring_centre_info_t(const widgeted_atom_t &at) : atom(at) {
      has_ring_centre_flag = 0;
   }
   void add_ring_centre(const lig_build::pos_t &pos) {
      ring_centre = pos;
      has_ring_centre_flag = 1;
   } 
};
std::ostream& operator<<(std::ostream &s, widgeted_atom_ring_centre_info_t wa);



// ====================================================================
//                     widgeted_molecule_t
// ====================================================================


class widgeted_molecule_t : public lig_build::molecule_t<widgeted_atom_t, widgeted_bond_t> {
   
private:
   std::string group;
   void init() {
      mol_in_max_y = 0;
      mol_in_min_y = 0;
      scale_correction.first = 0;
      scale_correction.second = 1;
      // have_cached_bond_ring_centres_flag = 0; in base class now
   }
   // Return a vector of bonds.  If empty, then it didn't find self.
   //

   // 20111229
//    std::pair<bool, std::vector<int> >
//    found_self_through_bonds(int atom_index_start, int atom_index_other) const;
//    std::pair<bool, std::vector<int> >
//    find_bonded_atoms_with_no_pass(int atom_index_start,
// 				  int atom_index_other, // must pass through this
// 				  int this_atom_index,
// 				  const std::vector<int> &no_pass_atoms,
// 				  int depth) const;
   
   // void debug_pass_atoms(int atom_index, int this_atom_index, 
   // int depth,  const std::vector<int> &local_no_pass_atoms) const;
   
   std::pair<bool, double>
      get_scale_correction(const lig_build::molfile_molecule_t &mol_in) const;
   
   // int get_number_of_atom_including_hydrogens() const;
   // return negative if not solvent accessibility available.
   double get_solvent_accessibility(const clipper::Coord_orth &pt,
				    const std::vector<solvent_accessible_atom_t> &sa) const;

   // This uses mmdb and clipper, so is not in lig-build.hh
   // 
   std::string get_atom_name(const clipper::Coord_orth &pt, CMMDBManager *mol) const;

public:
   widgeted_molecule_t() { init(); }
   widgeted_molecule_t(const lig_build::molfile_molecule_t &mol_in, CMMDBManager *pdb_mol);

   // return 0 as first if not highlighting a bond
   std::pair<bool, widgeted_bond_t> highlighted_bond_p(int x, int y) const;

   // return -1 as the atom index if not highlighting an atom.
   std::pair<int, widgeted_atom_t> highlighted_atom_p(int x, int y) const;
   bool write_mdl_molfile(const std::string &file_name) const;
   bool write_minimal_cif_file(const std::string &file_name) const;
   bool close_bond(int ib, GooCanvasItem *root, bool handle_post_delete_stray_atoms_flag);
   bool close_atom(int iat, GooCanvasItem *root);
   // std::vector<int> get_unconnected_atoms() const; 20111229 base class now

   // don't count closed bonds.
   // 
   // std::vector<int> bonds_having_atom_with_atom_index(int test_atom_index) const;

   // bool operator==(const widgeted_molecule_t &mol_other) const;
   // int n_stray_atoms() const; // unbonded atoms
   // std::vector<int> stray_atoms() const;
   // void translate(const lig_build::pos_t &delta); // move the atoms 20111229 base class now

   lig_build::pos_t input_coords_to_canvas_coords(const clipper::Coord_orth &in) const;
      

   // make private when bug is fixed.
   lig_build::pos_t centre_correction;
   std::pair<bool, double> scale_correction;
   double mol_in_min_y;
   double mol_in_max_y;

   // how much are the atoms of this moleclue scaled up (c.f. vs 1.5
   // units for a bond length) and what is the centre (in atom coords)
   // of the atom.
   // 
   // return a number greater than 0 when we have enough bonds and
   // atoms to determine this (we need this function because a
   // widgeted molecule has atoms in canvas coords and when we make
   // and rdkit molecule from a widgeted_molecule, we should put the
   // atoms back on sensible molecule scale. see lbg_info_t::rdkit_mol()).
   // 
   std::pair<double, lig_build::pos_t> current_scale_and_centre() const;

   void map_solvent_accessibilities_to_atoms(std::vector<solvent_accessible_atom_t> solvent_accessible_atoms);

   // can throw an exception
   // 
   lig_build::pos_t get_atom_canvas_position(const std::string &atom_name) const;

   // can throw an exception (no atoms)
   // 
   // lig_build::pos_t get_ligand_centre() const; // 20111229 base class now

   // If the ligand is not created in lbg, (as is the case if it gets
   // its ligand from Coot/Prodrg-FLAT) then these are not read and
   // correct ring centres, they are just points that are on the
   // correct side of the ring (so that the double bond can be drawn
   // on the correct side).  However, these control points are often
   // at ring centres and can be used to reject residue cirlces in the
   // position refinement.  In my tests so far there are considerably
   // more control points than genuine ring centres.
   //
   // in coot/src/lbg-graph.hh/cc, there is an algorithm to find
   // aromatic ring centres.  Maybe we should use that (extending it
   // to find ring non-aromatic centres) instead of the algorithm now
   // in get_ring_centres().
   // 
   // This function can throw an exception (no bonds).
   // 
   // (not const because it caches the return value)
   // 
   // std::vector<lig_build::pos_t> get_ring_centres(); 20111229 base class now

// base class now   
//    // can throw an exception (no atoms)
//    // 
//    lig_build::pos_t get_ring_centre(const std::vector<std::string> &ring_atom_names) const;

//    // can throw an exception (no rings with this atom)
//    //
//    lig_build::pos_t get_ring_centre(const widgeted_atom_ring_centre_info_t &atom) const;


//    // can throw an exception (no atoms) - top-left (small small)
//    // bottom-right (high high)
//    //
//    std::pair<lig_build::pos_t, lig_build::pos_t> ligand_extents() const;

//    int n_open_bonds() const;

//    bool is_close_to_non_last_atom(const lig_build::pos_t &test_post) const;

   void delete_hydrogens(GooCanvasItem *root);

};

// -----------------------------------------------------------------
//                   chirality
// -----------------------------------------------------------------

class topological_equivalence_t {

   // internal copy of input params.
   // 
   std::vector<widgeted_atom_t> atoms;
   std::vector<widgeted_bond_t> bonds;
   
   std::vector<bool> unique; // is the atom index marked as unique? initially all 0.
   std::vector<int> isn;     // invariant-sequence-numbers
   std::map<std::string, std::vector<int> > atom_map;

   // the number of different EC values in the molecules.
   int n_extended_connectivity(const std::vector<long int> &equivalent_classes) const;

   bool continue_ec_calculations_p(const std::vector<long int> &curr_eqv, 
				   const std::vector<long int> &prev_eqv);
   
   // fiddles with unique, return true if at least one unique was assigned.
   bool assign_uniques(const std::vector<long int> &extended_connectivity);

   void assign_invariant_sequence_number(const std::vector<long int> &curr_ec);

   // return the next isn index to be used.
   int assign_invariant_sequence_number(const std::vector<long int> &curr_ec,
					const std::vector<std::pair<int, int> > &atom_index,
					int next_index);

   bool atoms_have_unassigned_isn_p() const;

   // return a flag to let us know that it was done.
   bool mark_isn(int atom_index, int i_s_n); 

   // old
   // fiddles with unique
   bool identified_unique_p(const std::vector<long int> &curr_eqv, 
			    const std::vector<long int> &prev_eqv);

   std::vector<long int> assign_initial_topo_indices();
   std::vector<long int> assign_topo_indices(const std::vector<long int> &prev_eqv,
					     int round);

   // Return a list of atom indices that are connected to 3 or 4 other
   // atoms (return the indices of those other atoms too.
   // 
   std::vector<std::pair<int, std::vector<int> > > tetrahedral_atoms() const;

      
public:
   topological_equivalence_t(const std::vector<widgeted_atom_t> &atoms,
			     const std::vector<widgeted_bond_t> &bonds);

   std::vector<std::string> chiral_centres() const;
};


#endif // WMOLECULE_HH
