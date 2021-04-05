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

#include <string>
#include <goocanvas.h>
#include <clipper/core/coords.h>

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
                                                   GooCanvasAnchorType anchor_type,
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

      // For reasons not clear to me, the default line width is 2 when we get here.
      // We can force line width here (1.5 is nice) but I'd rather put it on the canvas
      // root and let this item inherit it from there
      
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

   // for 5 points darts - this is ugly
   virtual GooCanvasItem *
   wrap_goo_canvas_polyline_new(GooCanvasItem *root,
				double sharp_point_2_x,   double sharp_point_2_y, 
				double sharp_point_1_x,   double sharp_point_1_y, 
				double in_point_x,        double in_point_y, 
				double short_edge_pt_1_x, double short_edge_pt_1_y,
				double short_edge_pt_2_x, double short_edge_pt_2_y,
				std::string fc, std::string sc) const {

      GooCanvasItem *item = 
	 goo_canvas_polyline_new(root, 
				 TRUE, 5,
				 sharp_point_2_x, sharp_point_2_y, 
				 sharp_point_1_x, sharp_point_1_y, 
				 in_point_x, in_point_y, 
				 short_edge_pt_1_x, short_edge_pt_1_y,
				 short_edge_pt_2_x, short_edge_pt_2_y,
				 "fill-color", fc.c_str(),
				 "stroke-color", sc.c_str(),
				 NULL);
      return item;
   }

   virtual GooCanvasItem *
   wrap_goo_canvas_polyline_new_vp(GooCanvasItem *root,
				   const std::vector<lig_build::pos_t> &pts,
				   std::string fc,
				   std::string sc) const {

      GooCanvasItem *item = NULL;
      if (pts.size() == 6) {
	 item = goo_canvas_polyline_new(root, 
					TRUE, 6,
					pts[0].x, pts[0].y,
					pts[1].x, pts[1].y,
					pts[2].x, pts[2].y,
					pts[3].x, pts[3].y,
					pts[4].x, pts[4].y,
					pts[5].x, pts[5].y,
					"fill-color", fc.c_str(),
					"stroke-color", sc.c_str(),
					NULL);
      } else {
	 if (pts.size() == 5) {
	    item = goo_canvas_polyline_new(root, 
					   TRUE, 5,
					   pts[0].x, pts[0].y,
					   pts[1].x, pts[1].y,
					   pts[2].x, pts[2].y,
					   pts[3].x, pts[3].y,
					   pts[4].x, pts[4].y,
					   "fill-color", fc.c_str(),
					   "stroke-color", sc.c_str(),
					   NULL);
	 }
      }
      return item;
   }
   
   
   void wrap_goo_canvas_item_rotate(GooCanvasItem *ci,
					 double degrees, double cx, double cy) const {
      goo_canvas_item_rotate(ci, degrees, cx, cy);
   } 
};

#include "w-atom.hh"

#include "w-bond.hh"


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
   std::string get_atom_name(const clipper::Coord_orth &pt, mmdb::Manager *mol) const;

public:
   widgeted_molecule_t() { init(); }
   widgeted_molecule_t(const lig_build::molfile_molecule_t &mol_in, mmdb::Manager *pdb_mol);
   virtual ~widgeted_molecule_t();

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

   // Moved down
   // to draw double bonds without centre correctly (and below)
   // std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> >
   // make_other_connections_to_first_atom_info(unsigned int bond_index) const;
   // to draw wedge bonds correctly
   //    std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> >
   // make_other_connections_to_second_atom_info(unsigned int bond_index) const;
   
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

   enum {X_AXIS, Y_AXIS}; 
   void flip(int axis); // X or Y
   void rotate_z(double angle);  // in degrees

};

#include "topological-equivalence.hh"

#endif // WMOLECULE_HH
