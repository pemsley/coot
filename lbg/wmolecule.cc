/* lbg/wmolecule.cc
 * 
 * Author: Paul Emsley
 * Copyright 2010, 2011 by The University of Oxford
 * Copyright 2013, 2016 by Medical Research Council
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

#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "lbg.hh"

widgeted_molecule_t::~widgeted_molecule_t() {}

template<class widgeted_atom_t, class widgeted_bond_t> lig_build::molecule_t<widgeted_atom_t, widgeted_bond_t>::~molecule_t() {}

// Don't forget that this function will be used in
// render_from_molecule, which will add canvas item.
//
// Perhaps this should be in the base class, as it doesn't use widgets
// at all.
// 
widgeted_molecule_t::widgeted_molecule_t(const lig_build::molfile_molecule_t &mol_in,
					 mmdb::Manager *pdb_mol) { 

   bool debug_local = false;

   // the input coordinates are not necessarily centred on (0,0), so
   // let's find the centre of the input molecule first.
   //
   double centre_x = 0;
   double centre_y = 0;
   double sum_x = 0;
   double sum_y = 0;
   mol_in_min_y =  9999999;
   mol_in_max_y = -9999999;
   
   for (unsigned int iat=0; iat<mol_in.atoms.size(); iat++) {
      sum_x += mol_in.atoms[iat].atom_position.x();
      sum_y += mol_in.atoms[iat].atom_position.y();
      if (mol_in.atoms[iat].atom_position.y() > mol_in_max_y)
	 mol_in_max_y = mol_in.atoms[iat].atom_position.y();
      if (mol_in.atoms[iat].atom_position.y() < mol_in_min_y)
	 mol_in_min_y = mol_in.atoms[iat].atom_position.y();
   }

   // Sometimes (e.g. drugbank sdf files) the bonds are short (they
   // seem to be 1.0, or thereabouts), the atoms of the molecule need
   // to be scaled up.  So, to test if the molecule has short bonds,
   // calculate the median length of the bonds (not including bonds to
   // hydrogens).

   scale_correction = get_scale_correction(mol_in);
   
   if (mol_in.atoms.size() > 0) {
      centre_x = sum_x/double(mol_in.atoms.size());
      centre_y = sum_y/double(mol_in.atoms.size());

      centre_correction = lig_build::pos_t(centre_x, centre_y);

      if (debug_local)
	 std::cout << "::::::::::::::: y stats: extents " << mol_in_min_y << " "
		   << mol_in_max_y << " centre correction: " << centre_correction.y
		   << std::endl;
   
      for (unsigned int iat=0; iat<mol_in.atoms.size(); iat++) {

	 clipper::Coord_orth pt = mol_in.atoms[iat].atom_position;
	 lig_build::pos_t pos = input_coords_to_canvas_coords(pt);
	 std::string element = mol_in.atoms[iat].element;
	 int charge = mol_in.atoms[iat].formal_charge;
	 GooCanvasItem *ci = NULL;
	 widgeted_atom_t at(pos, element, charge, ci);

	 std::string current_atom_name = mol_in.atoms[iat].name;

	 // if the atom name was not set already, try to set it.
	 // 
	 // if (current_atom_name == "") {
	 if (current_atom_name.length() == 1) {

	    std::string atom_name = get_atom_name(pt, pdb_mol);
	    // and if atom_name is not "", then set atom name of at.
	    if (atom_name != "") {
	       if (debug_local) 
		  std::cout << "Hoorah setting atom at " << pt.format() << " with name :"
			    << atom_name << ":" << std::endl;
	       at.set_atom_name(atom_name);
	    } else {
	       if (0) 
		  std::cout << "Boo!! atom at " << pt.format() << " with no name"
			    << std::endl;
	    }
	 } else {
	    at.set_atom_name(current_atom_name);
	 } 
	 
	 if (0)
	    std::cout << "Element " << element << " at " << pos << " "
		      << mol_in_min_y << " "
		      << mol_in_max_y << std::endl;
	 // std::cout << "::::::: pushing back widgeted atom " << at
	 // << " with name :"  << at.get_atom_name() << ":" << std::endl;
	 atoms.push_back(at);
      }

      for (unsigned int ib=0; ib<mol_in.bonds.size(); ib++) {
	 int index_1 = mol_in.bonds[ib].index_1;
	 int index_2 = mol_in.bonds[ib].index_2;
	 lig_build::bond_t::bond_type_t bt = mol_in.bonds[ib].bond_type;
	 GooCanvasItem *ci = NULL;
	 bool shorten_first  = false;
	 bool shorten_second = false;
	 if (mol_in.atoms[index_1].element != "C")
	    shorten_first = true;
	 if (mol_in.atoms[index_2].element != "C")
	    shorten_second = true;

	 std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > empty;
	 widgeted_bond_t bond(index_1, index_2, atoms[index_1], atoms[index_2],
			      shorten_first, shorten_second, bt, empty, empty, ci);
	 bonds.push_back(bond);
      }

      assign_ring_centres();
   }
}

GooCanvasItem *
widgeted_atom_t::make_canvas_text_item(const lig_build::atom_id_info_t &atom_id_info_in,
				       const std::string &fc,
				       GooCanvasItem *root) {

   GooCanvasItem *group = NULL;
   if (atom_id_info_in.atom_id != "C") {
      group = wrap_goo_canvas_group_new (root, fc);
      // run through each of the offsets (i.e. each letter)
      for (unsigned int i=0; i<atom_id_info_in.n_offsets(); i++) {
	 double x_o = -5; 
	 double y_o =  8;
	 if (atom_id_info_in[i].text_pos_offset == lig_build::offset_text_t::UP)
	    y_o += 12;
	 if (atom_id_info_in[i].text_pos_offset == lig_build::offset_text_t::DOWN)
	    y_o += -12;

	 std::string font = "Sans 10";

	 if (atom_id_info_in.size_hint == -1)
	    font = "Sans 6";
	 
	 double x_pos = atom_position.x + atom_id_info_in.offsets[i].tweak.x + x_o;
	 double y_pos = atom_position.y + atom_id_info_in.offsets[i].tweak.y + y_o;

	 // hacketty hack because smaller font go down and to the left relatively
	 // (left is good because atom name are wider) - but let's compensate
	 // for the downness
	 if (atom_id_info_in.size_hint == -1) {
	    y_pos -= 3;
	 }
	 
	 // subscripts are lower and smaller font
	 if (atom_id_info_in.offsets[i].subscript) { 
	    font = "Sans 8";
	    y_pos += 3;
	 }
	 if (atom_id_info_in.offsets[i].superscript) { 
	    font = "Sans 8"; // 6 is too small for O-
	    y_pos -= 6;
	 }

	 if (0)
	    std::cout << "Rendering :" << atom_id_info_in[i].text << ": with tweak "
		      << atom_id_info_in[i].tweak << " and y-offset due to UP/DOWN "
		      << y_o << std::endl;
	 
	 GooCanvasItem *item =
	    wrap_goo_canvas_text_new(group,
				     atom_id_info_in.offsets[i].text.c_str(),
				     x_pos, y_pos, 
				     -1,
				     // GTK_ANCHOR_CENTER,
				     GOO_CANVAS_ANCHOR_SW,
				     font, fc);
      }
   }
   return group;
}



std::pair<bool, double>
widgeted_molecule_t::get_scale_correction(const lig_build::molfile_molecule_t &mol_in) const {
   return mol_in.get_scale_correction();
}


lig_build::pos_t
widgeted_molecule_t::input_coords_to_canvas_coords(const clipper::Coord_orth &pos_in) const {

   // convert from JME-style (top canvas atoms have big Y) to internal
   // coordinates (top canvas has small Y).

   if (0)
      std::cout << "   scale_correction  " << scale_correction.first << " " << scale_correction.second
		<< " centre_correction " << centre_correction << " mol_in_max_y " << mol_in_max_y
		<< " mol_in_min_y " << mol_in_min_y << std::endl;

   // scale_correction is typically 0.9 - and is: how much do I need
   // to scale up the coordinates of this molecule so that the average
   // bond length is 1.0?

   double x =   scale_correction.second * (pos_in.x() - centre_correction.x) * SINGLE_BOND_CANVAS_LENGTH;
   double y = - scale_correction.second * (pos_in.y() - centre_correction.y) * SINGLE_BOND_CANVAS_LENGTH;

   // double y_offset = 60 + scale_correction.second * (centre_correction.y - mol_in_min_y) * 20;
   double y_offset = 40 + scale_correction.second * (centre_correction.y - mol_in_min_y) * 30;

   if (0)
      std::cout << "input_coords_to_canvas_coords()" << pos_in.format()
		<< " y_offset: " << y_offset << " from " << scale_correction.first << " "
		<< scale_correction.second << " and mol_in_min_y "
		<< mol_in_min_y << std::endl;
   
   x += 300;
   y += y_offset;

   lig_build::pos_t pos(x,y);
   return pos;
}


// how much are the atoms of this moleclue scaled up (c.f. vs 1.3
// units for a bond length) and what is the centre (in atom coords)
// of the atom.
// 
std::pair<double, lig_build::pos_t>
widgeted_molecule_t::current_scale_and_centre() const {

   double scale = -1;
   lig_build::pos_t centre;
   
   int n_atoms = 0;
   double sum_x = 0.0;
   double sum_y = 0.0;
   for (unsigned int iat=0; iat<atoms.size(); iat++) { 
      if (! atoms[iat].is_closed()) {
	 sum_x += atoms[iat].atom_position.x;
	 sum_y += atoms[iat].atom_position.y;
	 n_atoms++;
      } 
   }

   if (n_atoms > 0) {
      double inv_n = 1.0/double(n_atoms);
      centre = lig_build::pos_t(sum_x * inv_n, sum_y * inv_n);

      double sum_bond_lengths = 0.0;
      int n_bond_lengths = 0; 
      for (unsigned int ibond=0; ibond<bonds.size(); ibond++) {
	 if (! bonds[ibond].is_closed()) {
	    lig_build::pos_t p1 = atoms[bonds[ibond].get_atom_1_index()].atom_position;
	    lig_build::pos_t p2 = atoms[bonds[ibond].get_atom_2_index()].atom_position;
	    double b = lig_build::pos_t::length(p2,p1);
	    sum_bond_lengths += b;
	    n_bond_lengths++;
	 }
      }

      if (n_bond_lengths) {
	 scale = sum_bond_lengths/(1.5*double(n_bond_lengths));
      }
   }

   return std::pair<double, lig_build::pos_t> (scale, centre);
} 



// void
// widgeted_molecule_t::translate(const lig_build::pos_t &delta) {

//    for (unsigned int iat=0; iat<atoms.size(); iat++) {
//       atoms[iat].atom_position += delta;
//    }
//    for (unsigned int ib=0; ib<bonds.size(); ib++) {
//       if (bonds[ib].have_centre_pos()) {
// 	 bonds[ib].move_centre_pos(delta);
//       }
//    }
// }


// All bonds are made this way
//
// simple (non-ring system) bonds.
//
// now deals with stereo_out/wedged/OUT_BOND
// 
GooCanvasItem *
widgeted_bond_t::canvas_item_for_bond(const lig_build::atom_t &at_1,
				      const lig_build::atom_t &at_2,
				      bool shorten_first,
				      bool shorten_second,
				      bond_type_t bt,
				      const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
				      const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
				      GooCanvasItem *root) const {

   // We can shorten the bonds to the atoms by different amounts, eg. N+ = 0
   // needs to be shorted asymmetrically.
   // 
   double shorten_fraction_1 = 0.72; // was 0.76
   double shorten_fraction_2 = 0.72;
   
   lig_build::pos_t pos_1_in = at_1.atom_position;
   lig_build::pos_t pos_2_in = at_2.atom_position;

   // Note: deltas from the south-east direction neighbour have a delta of ~ [20,20].
   // 
   // 20160628: When the bond comes in from the right and we have a Cl we need extra shortening.
   //           To do that we need to be passed the atom info (now done)
   //           These tweaks may need extension in future.
   //
   if (bt == SINGLE_BOND) {
      lig_build::pos_t delta;
      if (at_1.element == "Cl" ||
	  at_2.element == "Cl" ||
	  at_1.element == "Br" ||
	  at_2.element == "Br") {
	 if (at_1.element == "Cl" || at_1.element == "Br") {
	    // delta is the difference vector from the Cl to the other atom
	    delta = at_2.atom_position - at_1.atom_position;
	 }
	 if (at_2.element == "Cl" || at_2.element == "Br") {
	    delta = at_1.atom_position - at_2.atom_position;
	 }
	 if (delta.x > 10) {
	    // shorten both, but the flag is not set for the bond to the (presumably) C.
	    shorten_fraction_1 -= 0.13 * delta.x/26.8;
	    shorten_fraction_2 -= 0.13 * delta.x/26.8;
	 }
      }
   }

   if (at_1.element == "N" || at_2.element == "N") {
      if (at_1.element == "N") {
	 if (at_1.charge == 1) {
	    // N+ : shorten most if bond comes in from the NE (delta is [X, -Y])
	    //
	    lig_build::pos_t delta = at_2.atom_position - at_1.atom_position;
	    double f = delta.x - delta.y; // from 0 to 40;
	    double theta = delta.theta();
	    // I want sc_1 to maximize at theta = -45 degrees (bond from NE corner)
	    double sc_1 = 0.5 * (1.0 + cos(theta - M_PI_4)); // 0 -> 1
	    sc_1 *= sc_1; // sharpen
	    double sc_2 = 0.3 * sc_1;
	    shorten_fraction_1 -= + (sc_2 - 0.1);
	    if (false)
	       std::cout << "N is at_1"
			 << " delta: " << delta
			 << " theta " << (180/M_PI) * theta
			 << " delta-theta " << (180 / M_PI) * (theta - M_PI_4)
			 << " sc_1 " << sc_1
			 << std::endl;
	 }
      }
      if (at_2.element == "N") {
	 if (at_2.charge == 1) {
	    // N+ : shorten if bond comes in from the NE (as above)
	    lig_build::pos_t delta = at_1.atom_position - at_2.atom_position;
	    double theta = delta.theta();
	    double sc_1 = 0.5 * (1.0 + cos(theta - M_PI_4)); // 0 -> 1
	    sc_1 *= sc_1;
	    double sc_2 = 0.3 * sc_1;
	    shorten_fraction_2 -= + (sc_2 - 0.1);
	    if (false)
	       std::cout << "N is at_2"
			 << " delta: " << delta
			 << " theta " << (180/M_PI) * theta
			 << " delta-theta " << (180 / M_PI) * (theta - M_PI_4)
			 << " sc_1 " << sc_1
			 << std::endl;
	 }
      }
   }

   lig_build::pos_t pos_1 = pos_1_in;
   lig_build::pos_t pos_2 = pos_2_in;
   //
   // fraction_point() returns a point that is (say) 0.8 of the way
   // from p1 (first arg) to p2 (second arg).
   // 
   if (shorten_first)
      pos_1 = lig_build::pos_t::fraction_point(pos_2_in, pos_1_in, shorten_fraction_1);
   if (shorten_second)
      pos_2 = lig_build::pos_t::fraction_point(pos_1_in, pos_2_in, shorten_fraction_2);

   GooCanvasItem *ci = NULL;
   
   switch (bt) {
   case SINGLE_BOND:
      // Add new cases, a bit of a hack of course.
   case SINGLE_OR_DOUBLE:
   case SINGLE_OR_AROMATIC:
   case AROMATIC_BOND:        // this should not happen 
   case DELOC_ONE_AND_A_HALF: // this should not happen either
   case BOND_ANY:
      {
	 bool at_1_is_singleton = false;
	 bool at_2_is_singleton = false;
	 if (other_connections_to_first_atom.size()  == 0) at_1_is_singleton = true;
	 if (other_connections_to_second_atom.size() == 0) at_2_is_singleton = true;
	 std::pair<lig_build::pos_t, lig_build::pos_t> bp =
	    coords_for_single_bond(at_1, at_2, at_1_is_singleton, at_2_is_singleton);
	 ci = wrap_goo_canvas_polyline_new_line(root,
						pos_1.x, pos_1.y,
						pos_2.x, pos_2.y,
						"stroke-color", dark);
      }
      break;
   case DOUBLE_BOND:
   case DOUBLE_OR_AROMATIC:
      {
	 if (have_centre_pos()) {

	    // we want to draw this sort of double bond
	    // for double bonds that have atoms (other than the bond atoms) connect to each
	    // atom of the bond.
	    // i.e. we should draw a simple/symmetric double bond (canvas_item_double_bond) 
	    // only if one or other (or both) of the atoms is not connected to other atoms.
	    //
	    // And to make the decision, this function needs to be passed
	    // other_connections_to_first_atom.
	    //
	    ci = canvas_item_double_with_shortened_side_bond(pos_1, pos_2, root);
	 } else {
	    
	    ci = canvas_item_double_bond(pos_1, pos_2,
					 other_connections_to_first_atom,
					 other_connections_to_second_atom, root);
	 } 
      }
      break;
   case TRIPLE_BOND:
      { 
	 GooCanvasItem *group = wrap_goo_canvas_group_new (root, dark);
      
	 lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
	 lig_build::pos_t buv_90 = buv.rotate(90);
	 double small = 4;
	 lig_build::pos_t p1 = pos_1 + buv_90 * small;
	 lig_build::pos_t p2 = pos_2 + buv_90 * small;
	 lig_build::pos_t p3 = pos_1;
	 lig_build::pos_t p4 = pos_2;
	 lig_build::pos_t p5 = pos_1 - buv_90 * small;
	 lig_build::pos_t p6 = pos_2 - buv_90 * small;
	 GooCanvasItem *ci_1 = wrap_goo_canvas_polyline_new_line(group,
								 p1.x, p1.y,
								 p2.x, p2.y);
	 GooCanvasItem *ci_2 = wrap_goo_canvas_polyline_new_line(group,
								 p3.x, p3.y,
								 p4.x, p4.y);
	 GooCanvasItem *ci_3 = wrap_goo_canvas_polyline_new_line(group,
								 p5.x, p5.y,
								 p6.x, p6.y);
	 ci = group;
      }
      break;
   case IN_BOND:
      ci = make_wedge_bond_item(pos_1, pos_2, bt, other_connections_to_second_atom, root);
      break;
   case OUT_BOND:
      // to be like MarvinSketch we need to know what is connected to the non-pointy end
      ci = make_wedge_bond_item(pos_1, pos_2, bt, other_connections_to_second_atom, root);
      break;
   case BOND_UNDEFINED:
      break;
   }
   return ci;
}

// pos_1 and pos_2 are pre-shortened as needed.
//
GooCanvasItem *
widgeted_bond_t::canvas_item_double_bond(const lig_build::pos_t &pos_1,
					 const lig_build::pos_t &pos_2,
					 const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
					 const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
					 GooCanvasItem *root) const {

   if ((other_connections_to_second_atom.size() == 0) ||
       (other_connections_to_first_atom.size()  == 0)) {
      return canvas_item_double_bond_simple(pos_1, pos_2, root);

   } else {
   
      // somewhat like shortened side bond

      GooCanvasItem *group = wrap_goo_canvas_group_new (root, dark);
      GooCanvasItem *ci_1 = wrap_goo_canvas_polyline_new_line(group, pos_1.x, pos_1.y, pos_2.x, pos_2.y);
      double bond_length = lig_build::pos_t::length(pos_2, pos_1); // shortened possibly.
      lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
      lig_build::pos_t buv_90 = buv.rotate(90);

      // if the bond is cis, we want to be on the inside of that, if it's trans, it doesn't matter.
      // So does pos_1 and pos_2 have neighbours that are on the same side as each other?

      double delta_x = pos_2.x-pos_1.x;
      if (delta_x == 0) delta_x = 0.01; // don't divide by zero
      double m = (pos_2.y-pos_1.y)/delta_x;
      double c = - m * pos_2.x + pos_2.y;

      bool done = false; // true when pos_neighb_ref is set
      lig_build::pos_t pos_neighb_ref;
      
      for (unsigned int ib1=0; ib1<other_connections_to_first_atom.size(); ib1++) {
	 for (unsigned int ib2=0; ib2<other_connections_to_second_atom.size(); ib2++) {
	    const lig_build::pos_t &p_1 = other_connections_to_first_atom[ib1].first.atom_position;
	    const lig_build::pos_t &p_2 = other_connections_to_second_atom[ib2].first.atom_position;

	    double r1 = m * p_1.x - p_1.y + c;  // ax + by + c > 0 ? where a = m, b = -1
	    double r2 = m * p_2.x - p_2.y + c;
	    double r1r2 = r1 * r2;
	    if (r1r2 > 0) { // on same side

	       pos_neighb_ref = p_1;
	       done = true;
	    }
	    if (done) break;
	 }
	 if (done) break;
      }
      if (! done) {
	 // we had a trans bond (only) - either side will do
	 if (other_connections_to_first_atom.size()) {
	    pos_neighb_ref = other_connections_to_first_atom[0].first.atom_position;
	    done = true;
	 } else {
	    if (other_connections_to_second_atom.size()) {
	       pos_neighb_ref = other_connections_to_second_atom[0].first.atom_position;
	       done = true;
	    }
	 }
      }

      if (done) {
	 double nice_dist = 5.0;
	 lig_build::pos_t pos_offset_bond_start_t1 = pos_1 + buv_90 * nice_dist;
	 lig_build::pos_t pos_offset_bond_start_t2 = pos_1 - buv_90 * nice_dist;
	 double d1 = lig_build::pos_t::length(pos_neighb_ref, pos_offset_bond_start_t1);
	 double d2 = lig_build::pos_t::length(pos_neighb_ref, pos_offset_bond_start_t2);

	 lig_build::pos_t sp = pos_offset_bond_start_t1;
	 if (d2 < d1)
	    sp = pos_offset_bond_start_t2;
	 lig_build::pos_t ep = sp + buv * bond_length;

	 lig_build::pos_t cutened_inner_start_point = lig_build::pos_t::fraction_point(sp, ep, 0.14);
	 lig_build::pos_t cutened_inner_end_point   = lig_build::pos_t::fraction_point(sp, ep, 0.86);
	       
	 GooCanvasItem *ci_2 =
	    wrap_goo_canvas_polyline_new_line(group,
					      cutened_inner_start_point.x, 
					      cutened_inner_start_point.y, 
					      cutened_inner_end_point.x, 
					      cutened_inner_end_point.y);
      }
      return group;
   }
}
					



// symmetric (both offset, no shortening of a side-bond)
//
GooCanvasItem *
widgeted_bond_t::canvas_item_double_bond_simple(const lig_build::pos_t &pos_1,
						const lig_build::pos_t &pos_2,
						GooCanvasItem *root) const { 
	 
   GooCanvasItem *group = wrap_goo_canvas_group_new (root, dark);
						     
   lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
   lig_build::pos_t buv_90 = buv.rotate(90);

   double small = 2;
   lig_build::pos_t p1 = pos_1 + buv_90 * small;
   lig_build::pos_t p2 = pos_2 + buv_90 * small;
   lig_build::pos_t p3 = pos_1 - buv_90 * small;
   lig_build::pos_t p4 = pos_2 - buv_90 * small;
   GooCanvasItem *ci_1 = wrap_goo_canvas_polyline_new_line(group,
						      p1.x, p1.y,
						      p2.x, p2.y);
   GooCanvasItem *ci_2 = wrap_goo_canvas_polyline_new_line(group,
							   p3.x, p3.y,
							   p4.x, p4.y);
   return group;
}


// for rings
GooCanvasItem *
widgeted_bond_t::canvas_item_double_with_shortened_side_bond(const lig_build::pos_t &pos_1,
							     const lig_build::pos_t &pos_2,
							     GooCanvasItem *root) const {

   GooCanvasItem *group = wrap_goo_canvas_group_new (root, dark);

   GooCanvasItem *ci_1 = wrap_goo_canvas_polyline_new_line(group,
							   pos_1.x, pos_1.y,
							   pos_2.x, pos_2.y);

   lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
   lig_build::pos_t buv_90 = buv.rotate(90);

   // Which side of the pos_1 -> pos_2 vector shall we put this bond?
   // 
   // So create a T piece, and measure the distance to the centre
   // point, if we are on the inside then the distance to the centre
   // will be shorter.
   //
   double nice_dist = 5.0;
   double bond_length = lig_build::pos_t::length(pos_2, pos_1); // shortened possibly.
   lig_build::pos_t test_pt_1 = pos_1 + buv_90 * nice_dist;
   lig_build::pos_t test_pt_2 = pos_1 - buv_90 * nice_dist;
   double d_1 = lig_build::pos_t::length(test_pt_1, centre_pos());
   double d_2 = lig_build::pos_t::length(test_pt_2, centre_pos());

   lig_build::pos_t inner_start_point = test_pt_1;
   if (d_2 < d_1)
      inner_start_point = test_pt_2;

   lig_build::pos_t inner_end_point = inner_start_point + buv * bond_length;

   // These values (0.1, 0.9) are approximately OK for 6-membered rings
   // (but create bonds that are slightly too long, in fact)
   // but definately need to be shorter for 3,4 membered rings.
   //
   // To get it right, we need to know the positions of the atoms that are connected
   // to the atoms of the double bonds (or slightly kludgey), we could work it out
   // by knowing how many bonds/atoms in the ring of which this bond is a member.
   // 
   lig_build::pos_t cutened_inner_start_point =
      lig_build::pos_t::fraction_point(inner_start_point, inner_end_point, 0.14); // was 0.1
   lig_build::pos_t cutened_inner_end_point =
      lig_build::pos_t::fraction_point(inner_start_point, inner_end_point, 0.86); // was 0.9
   
   GooCanvasItem *ci_2 =
      wrap_goo_canvas_polyline_new_line(group,
					cutened_inner_start_point.x, 
					cutened_inner_start_point.y, 
					cutened_inner_end_point.x, 
					cutened_inner_end_point.y);

      return group;
}

GooCanvasItem *
widgeted_bond_t::make_wedge_bond_item(const lig_build::pos_t &pos_1,
				      const lig_build::pos_t &pos_2,
				      const lig_build::bond_t::bond_type_t &bt,
				      const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
				      GooCanvasItem *root) const {

   GooCanvasItem *item = NULL;

   if (bt == lig_build::bond_t::OUT_BOND)
      item = make_wedge_out_bond_item(pos_1, pos_2, other_connections_to_second_atom, root);
   if (bt == lig_build::bond_t::IN_BOND)
      item = make_wedge_in_bond_item(pos_1, pos_2, root);

   return item;
}

// pos_1 is the position of the atom at the sharp point (the chiral centre).
// pos_2 is the position of the other atom.
GooCanvasItem *
widgeted_bond_t::make_wedge_out_bond_item(const lig_build::pos_t &pos_1,
					  const lig_build::pos_t &pos_2,
					  const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
					  GooCanvasItem *root) const {


   if (false) {
      std::cout << "   debug:: in make_wedge_out_bond_item() n-other-connections-to-second "
		<< other_connections_to_second_atom.size() << std::endl;
      for (unsigned int i=0; i<other_connections_to_second_atom.size(); i++) {
	 std::cout << "  "
		   << other_connections_to_second_atom[i].first << " "
		   << other_connections_to_second_atom[i].second << std::endl;
      }
   }

   GooCanvasItem *item = NULL; // updated as return value
   
   if (other_connections_to_second_atom.size() > 0) {
      item = make_sheared_or_darted_wedge_bond(pos_1, pos_2, other_connections_to_second_atom, root);
   } else {
   
      // A filled shape (quadralateral (almost triangle)).

      lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
      lig_build::pos_t buv_90 = buv.rotate(90);
      // How long is this bond?  The width of the fat (i.e. non-pointy) end should be
      // proprotional to the length of the bond.
      //
      double l = lig_build::pos_t::length(pos_1, pos_2);
      double bond_length_ratio = l/SINGLE_BOND_CANVAS_LENGTH;

      lig_build::pos_t short_edge_pt_1 = pos_2 + buv_90 * 3 * bond_length_ratio;
      lig_build::pos_t short_edge_pt_2 = pos_2 - buv_90 * 3 * bond_length_ratio;
   
      // the line width means that the sharp angle at pos_1 here results
      // in a few pixels beyond the pos_1, so artificially shorten it a
      // tiny amount.
      //
      // Also, make it a quadralateral, with the sharp points very close,
      // this make the spike go away.
      //
      // 
      // lig_build::pos_t sharp_point = lig_build::pos_t::fraction_point(pos_1, pos_2, 0.11);
      // lig_build::pos_t sharp_point = lig_build::pos_t::fraction_point(pos_1, pos_2, 0.07);
      lig_build::pos_t sharp_point = lig_build::pos_t::fraction_point(pos_1, pos_2, 0.04);
   
      lig_build::pos_t sharp_point_1 = sharp_point + buv_90 * 0.03;
      lig_build::pos_t sharp_point_2 = sharp_point - buv_90 * 0.03;

      item = wrap_goo_canvas_polyline_new(root, 
				      sharp_point_2.x, sharp_point_2.y,
				      sharp_point_1.x, sharp_point_1.y,
				      short_edge_pt_1.x, short_edge_pt_1.y,
				      short_edge_pt_2.x, short_edge_pt_2.y,
				      dark, dark);
   }
   return item;
}


// pos_1 is the position of the atom at the sharp point (the chiral centre).
// pos_2 is the position of the other atom.
GooCanvasItem *
widgeted_bond_t::make_sheared_or_darted_wedge_bond(const lig_build::pos_t &pos_1,
						   const lig_build::pos_t &pos_2,
						   const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
						   GooCanvasItem *root) const {
   GooCanvasItem *item = NULL;

   if (other_connections_to_second_atom.size() > 0) {

      if (other_connections_to_second_atom.size() == 1) {

	 const lig_build::pos_t  &third_atom_pos = other_connections_to_second_atom[0].first.atom_position;
	 const lig_build::bond_t &third_bond     = other_connections_to_second_atom[0].second;

	 lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
	 lig_build::pos_t buv_90 = buv.rotate(90);
	 lig_build::pos_t sharp_point = lig_build::pos_t::fraction_point(pos_1, pos_2, 0.04);
   
	 lig_build::pos_t sharp_point_1 = sharp_point + buv_90 * 0.03;
	 lig_build::pos_t sharp_point_2 = sharp_point - buv_90 * 0.03;

	 lig_build::pos_t bfrom3rd = pos_2 - third_atom_pos;
	 lig_build::pos_t bond_from_3rd_atom_extension   = pos_2 + bfrom3rd*0.1;
	 lig_build::pos_t bond_from_3rd_atom_contraction = pos_2 - bfrom3rd*0.18;

	 if (third_bond.get_bond_type() == lig_build::bond_t::DOUBLE_BOND) {
	    // we need to make this shorter.
	    bond_from_3rd_atom_extension   -= buv * 2.6;
	    bond_from_3rd_atom_contraction -= buv * 2.6;
	 }

	 if (false) {
	    std::cout << " buv             " << buv << std::endl;
	    std::cout << " pos_2           " << pos_2 << std::endl;
	    std::cout << " 3rd atom pos    " << third_atom_pos << std::endl;
	    std::cout << " bfrom3rd        " << bfrom3rd << std::endl;
	    std::cout << " sheared points: " << sharp_point_2 << std::endl;
	    std::cout << "                 " << sharp_point_1 << std::endl;
	    std::cout << "                 " << bond_from_3rd_atom_extension   << std::endl;
	    std::cout << "                 " << bond_from_3rd_atom_contraction << std::endl;
	 }

	 item = wrap_goo_canvas_polyline_new(root,
					     sharp_point_2.x, sharp_point_2.y,
					     sharp_point_1.x, sharp_point_1.y,
					     bond_from_3rd_atom_extension.x,
					     bond_from_3rd_atom_extension.y,
					     bond_from_3rd_atom_contraction.x,
					     bond_from_3rd_atom_contraction.y,
					     dark, dark);

      } else {

	 // make a dart (there are 2 third atoms)
	 
	 const lig_build::pos_t  &third_atom_1_pos = other_connections_to_second_atom[0].first.atom_position;
	 const lig_build::pos_t  &third_atom_2_pos = other_connections_to_second_atom[1].first.atom_position;

	 lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
	 lig_build::pos_t buv_90 = buv.rotate(90);
	 lig_build::pos_t sharp_point = lig_build::pos_t::fraction_point(pos_1, pos_2, 0.04);
	 lig_build::pos_t sharp_point_1 = sharp_point + buv_90 * 0.03;
	 lig_build::pos_t sharp_point_2 = sharp_point - buv_90 * 0.03;
   
	 lig_build::pos_t bfrom3rd_1 = pos_2 - third_atom_1_pos;
	 lig_build::pos_t bfrom3rd_2 = pos_2 - third_atom_2_pos;
	 lig_build::pos_t bond_from_3rd_atom_1_contraction = pos_2 - bfrom3rd_1*0.15;
	 lig_build::pos_t bond_from_3rd_atom_2_contraction = pos_2 - bfrom3rd_2*0.15;

	 std::vector<lig_build::pos_t> pts;
	 pts.push_back(sharp_point_2);
	 pts.push_back(sharp_point_1);
	 pts.push_back(bond_from_3rd_atom_1_contraction);
	 pts.push_back(pos_2);
	 pts.push_back(bond_from_3rd_atom_2_contraction);
	 item = wrap_goo_canvas_polyline_new_vp(root, pts, dark, dark);
      }
   }
   return item;
}


// pos_1 is at the sharp end.  Into the page.  A series of lines.
// 
GooCanvasItem *
widgeted_bond_t::make_wedge_in_bond_item(const lig_build::pos_t &pos_1,
					 const lig_build::pos_t &pos_2,
					 GooCanvasItem *root) const {
   
   GooCanvasItem *group = wrap_goo_canvas_group_new (root, dark);
   lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
   lig_build::pos_t buv_90 = buv.rotate(90);
   int n_lines = 5;

   // How long is this bond?  The width of the fat (i.e. non-pointy) end should be
   // proprotional to the length of the bond.
   //
   double l = lig_build::pos_t::length(pos_1, pos_2);
   double bond_length_ratio = l/SINGLE_BOND_CANVAS_LENGTH;

   for (int i=1; i<=n_lines; i++) {
      // then centre point of the line, some way along the pos_1 -> pos_2 vector;
      double len = double(i) * 1.0;
      double frac = (double(i)- 0.3)/double(n_lines);
      lig_build::pos_t fp = lig_build::pos_t::fraction_point(pos_1, pos_2, frac);
      lig_build::pos_t p1 = fp + buv_90 * len * bond_length_ratio;
      lig_build::pos_t p2 = fp - buv_90 * len * bond_length_ratio;
      GooCanvasItem *ci_1 = wrap_goo_canvas_polyline_new_line(group,
							 p1.x, p1.y,
							 p2.x, p2.y);
   }
   return group;
}


// // to draw wedge bonds correctly
// std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> >
// widgeted_molecule_t::make_other_connections_to_second_atom_info(unsigned int bond_index) const {

//    std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > v;
//    int atom_chiral_idx = bonds[bond_index].get_atom_1_index();
//    int atom_other_idx  = bonds[bond_index].get_atom_2_index();

//    for (unsigned int ibond=0; ibond<bonds.size(); ibond++) {
//       if (ibond != bond_index) {
// 	 int at_1_idx = bonds[ibond].get_atom_1_index();
// 	 int at_2_idx = bonds[ibond].get_atom_2_index();
// 	 if (at_1_idx == atom_other_idx) {
// 	    if (at_2_idx != atom_chiral_idx) { // should always be
// 	       std::pair<lig_build::atom_t, lig_build::bond_t> p(atoms[at_2_idx], bonds[ibond]);
// 	       v.push_back(p);
// 	    }
// 	 }
// 	 if (at_2_idx == atom_other_idx) {
// 	    if (at_1_idx != atom_chiral_idx) {
// 	       std::pair<lig_build::atom_t, lig_build::bond_t> p(atoms[at_1_idx], bonds[ibond]);
// 	       v.push_back(p);
// 	    }
// 	 }
//       }
//    }

//    // std::cout << "from make_other_connections_to_second_atom_info() returning v of size "
//    // << v.size() << std::endl;

//    return v;
// }

// std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> >
// widgeted_molecule_t::make_other_connections_to_first_atom_info(unsigned int bond_index) const {

//    std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > v;
//    int atom_1_ref_idx = bonds[bond_index].get_atom_1_index();
//    int atom_2_ref_idx = bonds[bond_index].get_atom_2_index();

//    for (unsigned int ibond=0; ibond<bonds.size(); ibond++) {
//       if (ibond != bond_index) {
// 	 int at_1_idx = bonds[ibond].get_atom_1_index();
// 	 int at_2_idx = bonds[ibond].get_atom_2_index();
// 	 if (at_1_idx == atom_1_ref_idx) {
// 	    if (at_2_idx != atom_2_ref_idx) {
// 	       std::pair<lig_build::atom_t, lig_build::bond_t> p(atoms[at_2_idx], bonds[ibond]);
// 	       v.push_back(p);
// 	    }
// 	 }
// 	 if (at_2_idx == atom_1_ref_idx) {
// 	    if (at_1_idx != atom_2_ref_idx) {
// 	       std::pair<lig_build::atom_t, lig_build::bond_t> p(atoms[at_1_idx], bonds[ibond]);
// 	       v.push_back(p);
// 	    }
// 	 }
//       }
//    }
//    return v;
// }


std::pair<bool, widgeted_bond_t>
widgeted_molecule_t::highlighted_bond_p(int x, int y) const {

   widgeted_bond_t b;

   std::pair<bool, widgeted_bond_t> r(0, b);
   double dx(x);
   double dy(y);
   
   for (unsigned int i=0; i<bonds.size(); i++) {
      if (! bonds[i].is_closed()) { 
	 if (bonds[i].over_bond(dx, dy,
				atoms[bonds[i].get_atom_1_index()],
				atoms[bonds[i].get_atom_2_index()])) {
	    r.first = 1;
	    r.second = bonds[i];
	    break;
	 }
      }
   }
   return r;
}

// return -1 as the atom index if not highlighting an atom.
std::pair<int, widgeted_atom_t>
widgeted_molecule_t::highlighted_atom_p(int x, int y) const {

   GooCanvasItem *ci = NULL;
   lig_build::atom_t at(lig_build::pos_t() , "", 0);
   widgeted_atom_t a(at, ci);
   std::pair<int, widgeted_atom_t> r(UNASSIGNED_INDEX, a);
   double dx(x);
   double dy(y);
   for (unsigned int i=0; i<atoms.size(); i++) {
      if (! atoms[i].is_closed()) { 
	 if (atoms[i].over_atom(dx, dy)) {
	    r.first = i;
	    r.second = atoms[i];
	    break;
	 }
      }
   }
   return r;
}


// int
// widgeted_molecule_t::n_stray_atoms() const {  // unbonded atoms
//    return stray_atoms().size();
// }

// std::vector<int>
// widgeted_molecule_t::stray_atoms() const {

//    std::vector<int> strays;
   
//    bool found[atoms.size()];

//    for (unsigned int i=0; i<atoms.size(); i++)
//       found[i] = 0;

//    for (unsigned int ib=0; ib<bonds.size(); ib++) { 
//       int iat_1 = bonds[ib].get_atom_1_index();
//       int iat_2 = bonds[ib].get_atom_2_index();
//       if (! atoms[iat_1].is_closed())
// 	 found[iat_1] = 1;
//       if (! atoms[iat_2].is_closed())
// 	 found[iat_2] = 1;
//    }

//    for (unsigned int i=0; i<atoms.size(); i++)
//       if (! found[i])
// 	 strays.push_back(i);

//    return strays;
// }

std::ostream&
operator<<(std::ostream &s, widgeted_atom_ring_centre_info_t wa) {
   s << wa.atom << " ring-centre: " << wa.has_ring_centre_flag;
   if (wa.has_ring_centre_flag) {
      s << " " << wa.ring_centre;
   }
   return s;
}



// bool
// widgeted_molecule_t::operator==(const widgeted_molecule_t &mol_other) const {

//    bool status = 0;

//    // need to check that bonds are the same (atom indexing can be
//    // different) and also stray atoms need to be checked after bonds.
//    //
//    int n_bond_hits = 0;

//    if (mol_other.bonds.size() != bonds.size())
//       return status;
//    for (unsigned int ib=0; ib<bonds.size(); ib++) {
//       lig_build::atom_t atom_i_1 = atoms[bonds[ib].get_atom_1_index()];
//       lig_build::atom_t atom_i_2 = atoms[bonds[ib].get_atom_2_index()];
//       int n_hits = 0;
//       for (unsigned int jb=0; jb<mol_other.bonds.size(); jb++) {
// 	 lig_build::atom_t atom_j_1 = mol_other.atoms[bonds[ib].get_atom_1_index()];
// 	 lig_build::atom_t atom_j_2 = mol_other.atoms[bonds[ib].get_atom_2_index()];
// 	 if (atom_i_1 == atom_j_1)
// 	    if (atom_i_2 == atom_j_2)
// 	       n_hits++;
//       }
//       if (n_hits == 1) {
// 	 n_bond_hits++;
//       }
//    }

//    if (n_bond_hits == bonds.size()) {

//       // So the bonds were the same, now the strays...

//       if (mol_other.n_stray_atoms() == n_stray_atoms()) {
// 	 std::vector<int> i_stray_atoms = stray_atoms();
// 	 std::vector<int> j_stray_atoms = mol_other.stray_atoms();
// 	 int n_stray_hits = 0;
// 	 for (unsigned int i=0; i<i_stray_atoms.size(); i++) {
// 	    for (unsigned int j=0; j<j_stray_atoms.size(); j++) {
// 	       if (atoms[i_stray_atoms[i]] == mol_other.atoms[j_stray_atoms[j]])
// 		  n_stray_hits++;
// 	    }
// 	 }
// 	 if (n_stray_hits == n_stray_atoms())
// 	    status = 1;
//       }

//    }
//    return status;
// }


// std::vector<int>
// widgeted_molecule_t::bonds_having_atom_with_atom_index(int test_atom_index) const {

//    std::vector<int> v;
//    std::vector<int> vb =  bond_indices_with_atom_index(test_atom_index);

//    for (unsigned int iv=0; iv<vb.size(); iv++) {
//       if (! bonds[vb[iv]].is_closed())
// 	 v.push_back(vb[iv]);
//    }

//    return v;
// } 



// // return a vector of the atom indices of unconnected atoms
// //
// std::vector<int>
// widgeted_molecule_t::get_unconnected_atoms() const {

//    std::vector<int> v;
//    for (unsigned int iat=0; iat<atoms.size(); iat++) {
//       if (! atoms[iat].is_closed()) { 
// 	 bool in_a_bond = 0;
// 	 for (unsigned int ib=0; ib<bonds.size(); ib++) { 
// 	    if (! bonds[ib].is_closed()) {
// 	       if (bonds[ib].get_atom_1_index() == iat)
// 		  in_a_bond = 1;
// 	       if (bonds[ib].get_atom_2_index() == iat)
// 		  in_a_bond = 1;
// 	    }
// 	    if (in_a_bond)
// 	       break;
// 	 }
// 	 if (! in_a_bond)
// 	    v.push_back(iat);
//       }
//    }
//    return v;
// }

// // can throw an exception (no atoms)
// // 
// lig_build::pos_t
// widgeted_molecule_t::get_ligand_centre() const {

//    lig_build::pos_t centre(0,0);

//    if (atoms.size() == 0) {
//       std::string message("No atoms in ligand");
//       throw std::runtime_error(message);
//    } else {
//       lig_build::pos_t centre_sum(0,0);
//       for (unsigned int iat=0; iat<atoms.size(); iat++) {
// 	 centre_sum += atoms[iat].atom_position;
//       }
//       if (atoms.size() > 0)
// 	 centre = centre_sum * (1.0/double(atoms.size()));
//    }
//    return centre;
// }

// // can throw an exception (no atoms)
// // 
// lig_build::pos_t
// widgeted_molecule_t::get_ring_centre(const std::vector<std::string> &ring_atom_names) const {

//    lig_build::pos_t positions_sum(0,0);
//    int n_found = 0;
//    for (unsigned int jat=0; jat<ring_atom_names.size(); jat++) {
//       for (unsigned int iat=0; iat<atoms.size(); iat++) {
// 	 if (ring_atom_names[jat] == atoms[iat].get_atom_name()) {
// 	    positions_sum += atoms[iat].atom_position;
// 	    n_found++;
// 	    break;
// 	 }
//       }
//    }
//    if (n_found == 0) {
//       std::string mess = "No ring atom names found in ligand!";
//       throw(std::runtime_error(mess));
//    }
//    lig_build::pos_t centre = positions_sum * (1.0/double(n_found));
//    return centre;
// }

// // can throw an exception (no rings with this atom)
// //
// lig_build::pos_t
// widgeted_molecule_t::get_ring_centre(const widgeted_atom_ring_centre_info_t &atom) const {

//    lig_build::pos_t position(0,0);
//    bool found = 0;

//    for (unsigned int ibond=0; ibond<bonds.size(); ibond++) { 
//       if ((atoms[bonds[ibond].get_atom_1_index()] == atom.atom) ||
// 	  (atoms[bonds[ibond].get_atom_2_index()] == atom.atom)) {
// 	 if (bonds[ibond].have_centre_pos()) {
// 	    position = bonds[ibond].centre_pos();
// 	    found = 1;
// 	 }
//       }
//       if (found)
// 	 break;
//    }

//    if (! found) {
//       std::string mess("No atom ");
//       mess += atom.atom.get_atom_name();
//       mess += " found to be in a ring";
//       throw(std::runtime_error(mess));
//    }
//    return position;
// }


// // can throw an exception (no bonds)
// //
// // not const because it now caches the return value;
// //
// std::vector<lig_build::pos_t>
// widgeted_molecule_t::get_ring_centres() {

//    if (have_cached_bond_ring_centres_flag) {
//       return cached_bond_ring_centres;
//    } else { 
//       std::vector<lig_build::pos_t> v;
//       for (unsigned int ib=0; ib<bonds.size(); ib++) {
// 	 if (bonds[ib].have_centre_pos()) {
// 	    lig_build::pos_t rc = bonds[ib].centre_pos();
// 	    bool found = 0;
// 	    for (unsigned int i=0; i<v.size(); i++) {
// 	       if (v[i].near_point(rc, 7)) { // 7 seems good, others tested.
// 		  found = 1;
// 		  break;
// 	       }
// 	    }
// 	    if (! found)
// 	       v.push_back(rc);
// 	 }
//       }
//       //    std::cout << "get_ring_centres returns\n";
//       //    for (unsigned int iv=0; iv<v.size(); iv++) {
//       //       std::cout << "   "  << iv << " " << v[iv] << "\n";
//       //    }
//       cached_bond_ring_centres = v;
//       have_cached_bond_ring_centres_flag = 1;
//       return v;
//    }
// }


// return was-really-closed status
bool
widgeted_molecule_t::close_bond(int ib, GooCanvasItem *root,
				bool handle_post_delete_stray_atoms_flag) {

   // on killing a bond, an N at one end of the bond may need its atom_id
   // changed to NH or so.

   bool status = false;
   int n_bonds = bonds.size();
   if ((ib >= 0) && (ib<n_bonds)) {
      int ind_1 = bonds[ib].get_atom_1_index();
      int ind_2 = bonds[ib].get_atom_2_index();
      bonds[ib].close(root);
      status = true;

      // ind_1
      std::string ele = atoms[ind_1].element;
      std::vector<unsigned int> local_bonds = bonds_having_atom_with_atom_index(ind_1);
      bool gl_flag = false; // not a GL render engine
      lig_build::atom_id_info_t new_atom_id =
	 make_atom_id_by_using_bonds(ind_1, ele, local_bonds, gl_flag);
      atoms[ind_1].update_atom_id_maybe(new_atom_id, root);
      // ind_2
      ele = atoms[ind_2].element;
      local_bonds = bonds_having_atom_with_atom_index(ind_2);
      new_atom_id = make_atom_id_by_using_bonds(ind_2, ele, local_bonds, gl_flag);
      atoms[ind_2].update_atom_id_maybe(new_atom_id, root);
   
      if (handle_post_delete_stray_atoms_flag) {
	 // are there any atoms that now don't have bonds?
	 std::vector<unsigned int> stray_atoms = get_unconnected_atoms();
	 std::vector<unsigned int> bonds; // empty
	 // std::cout << "got " << stray_atoms.size() << " stray atoms " << std::endl;
	 if (stray_atoms.size()) {
	    for (unsigned int istray=0; istray<stray_atoms.size(); istray++) {
	       const std::string &ele = atoms[stray_atoms[istray]].element;
	       lig_build::atom_id_info_t atom_id_info =
		  make_atom_id_by_using_bonds(ind_1, ele, bonds, gl_flag);
	       atoms[stray_atoms[istray]].update_atom_id_forced(atom_id_info, root);
	    }
	 }
      }
   }
   return status;
}

bool
widgeted_molecule_t::close_atom(int iat, GooCanvasItem *root) {

   std::vector<unsigned int> bds = bonds_having_atom_with_atom_index(iat);
   
   bool status = 0;
   if ((iat >= 0) && (iat<int(atoms.size()))) {
      atoms[iat].close(root);
      status = 1;
   }

   // erase the bonds that were attached to atom iat:
   //
   // the atoms attached at the other end of the bonds may need to
   // change (say from N to NH).  Make a list of the atoms affected by
   // deleting the bond(s).

   std::vector<unsigned int> affected_neighbour_atoms;

   for (unsigned int i=0; i<bds.size(); i++) {
      // std::cout << "close_atom() means closing bond " << bds[i] << std::endl;
      int ind1 = bonds[bds[i]].get_atom_1_index();
      int ind2 = bonds[bds[i]].get_atom_2_index();
      if (ind1 != iat) affected_neighbour_atoms.push_back(ind1);
      if (ind2 != iat) affected_neighbour_atoms.push_back(ind2);
      close_bond(bds[i], root, 1);
   }
   return status;
}








bool
widgeted_molecule_t::write_mdl_molfile(const std::string &file_name) const {

   bool status = 0;

   std::ofstream of(file_name.c_str());

   // we need to convert between atoms vector (which have closed
   // atoms) and mdl atoms, which are ordered and start from 1.
   // 
   int mdl_atoms[bonds.size()];
   int current_index = 1;
   int n_non_closed_atoms = 0;
   int n_non_closed_bonds = 0;
   for (unsigned int iat=0; iat<atoms.size(); iat++) {
      if (atoms[iat].is_closed()) { 
	 mdl_atoms[iat] = UNASSIGNED_INDEX;
      } else { 
	 mdl_atoms[iat] = current_index;
	 current_index++;
	 n_non_closed_atoms++;
      }
   }
   for (unsigned int ib=0; ib<bonds.size(); ib++) {
      if (! bonds[ib].is_closed()) {
	 n_non_closed_bonds++;
      }
   }
   

   if (! of) {
      std::cout << "WARNING:: Cannot open file :" << file_name << ": to write molfile"
		<< std::endl;
   } else {

      // 3 lines of header
      of << "# Molecule sketched in LIDIA (part of Coot)\n";
      of << "#\n";
      of << "#\n";
      
//       timeval start_time;
//       int time_status = gettimeofday(&start_time, NULL);
//       tm *lt = localtime(start_time);
      
      // next the counts line:
      of.width(3);

      // is this a chiral molecule?
      bool chiral = false;
      // yes, if any of the bonds are stereo
      for (unsigned int ib=0; ib<bonds.size(); ib++) {
	 const widgeted_bond_t &bond = bonds[ib];
	 if (! bond.is_closed()) {
	    if (bond.get_bond_type() == lig_build::bond_t::IN_BOND || 
		bond.get_bond_type() == lig_build::bond_t::OUT_BOND) { 
	       chiral = true;
	       break;
	    }
	 }
      }
     
      of << n_non_closed_atoms;
      of.width(3);
      of << n_non_closed_bonds;
      of << "  1   ";
      of.width(3);
      of << chiral;
      of << "  0" << "            " // obsolete entries
	 << "999" << std::endl;

      // let's calculate the centre of the atoms and take that away
      // from all the atom coordinates.
      //
      lig_build::pos_t centre_sum(0,0);
      lig_build::pos_t centre(0,0);
      for (unsigned int iat=0; iat<atoms.size(); iat++) {
	 centre_sum += atoms[iat].atom_position;
      }
      if (atoms.size() > 0)
	 centre = centre_sum * (1.0/double(atoms.size()));

      
      // atom table:
      for (unsigned int iat=0; iat<atoms.size(); iat++) {
	 if (! atoms[iat].is_closed()) {
	    // shall we invert the y coordinate? Yes.  Internal
	    // coordinates have the atoms top of the canvas with
	    // lowest Y values.  The JME has them with highest Y
	    // values, so swap on output.
	    of << std::setiosflags(std::ios::fixed) << std::setw(10) << std::setprecision(4)
	       << (atoms[iat].atom_position.x - centre.x) * 1.54/SINGLE_BOND_CANVAS_LENGTH; // file output
	    of << std::setiosflags(std::ios::fixed) << std::setw(10) << std::setprecision(4)
	       << (- atoms[iat].atom_position.y + centre.y) * 1.54/SINGLE_BOND_CANVAS_LENGTH;
	    double z = 0.0;
	    of << std::setiosflags(std::ios::fixed) << std::setw(10) << std::setprecision(4)
	       << z;
	    // of.width(3);
	    std::string ele = atoms[iat].element;
	    if (ele.length() == 1) { 
	       ele = " " + atoms[iat].element + " ";
	    } else {
	       ele = " " + ele;
	    }
	    of << ele;
	    int mass_diff = 0;
	    of.width(3);
	    of << mass_diff;

	    int file_charge = 0;
	    // int charge = 0;
	    // from ctfile.pdf:
	    // 1 means charge +3
	    // 2 means charge +2
	    // 3 means charge +1
	    // 4 = doublet radical
	    // 5 is charge -1
	    // 6 is charge -2
	    // 7 is charge -3
	    if (atoms[iat].charge == 3)
	       file_charge = 1;
	    if (atoms[iat].charge == 2)
	       file_charge = 1;
	    if (atoms[iat].charge == 1)
	       file_charge = 3;
	    if (atoms[iat].charge == -1)
	       file_charge = 5;
	    if (atoms[iat].charge == -1)
	       file_charge = 5;
	    if (atoms[iat].charge == -2)
	       file_charge = 6;
	    if (atoms[iat].charge == -7)
	       file_charge = 7;
	    of.width(3);
	    of << file_charge;
	    
	    int stereo_parity = 0;
	    of.width(3);
	    of << stereo_parity;
	    int hydrogen_count = 0; 
	    of.width(3);
	    of << hydrogen_count;
	    int stereo_care_box = 0; 
	    of.width(3);
	    of << stereo_care_box;
	    int valence = 0; 
	    of.width(3);
	    of << valence;
	    int H0_designator = 0; 
	    of.width(3);
	    of << H0_designator;
	    of << "      ";
	    of.width(3);
	    of << mdl_atoms[iat]; // maybe
	    int inversion_flag = 0;
	    of.width(3);
	    of << inversion_flag;
	    int exact_charge_flag = 0;
	    of.width(3);
	    of << exact_charge_flag;
	    of << "\n";
	 }
      }

      // bond table
      for (unsigned int ib=0; ib<bonds.size(); ib++) {
	 const widgeted_bond_t &bond = bonds[ib];
	 if (! bond.is_closed()) {

	    // bond stereo, for wedge bonds. 0=Not-Stereo, 1=Up, 6=Down,
	    // 3=Cis ot trans (either).
	    // The wedge (narrow) end of the stereo bond is the first atom.
	    // IN_BOND is Down and OUT_BOND is Up
	    // 
	    int bond_stereo = 0;
	    if (bond.get_bond_type() == lig_build::bond_t::IN_BOND)
	       bond_stereo = 6;
	    if (bond.get_bond_type() == lig_build::bond_t::OUT_BOND)
	       bond_stereo = 1;
	    
	    int idx_1 = mdl_atoms[bonds[ib].get_atom_1_index()];
	    int idx_2 = mdl_atoms[bonds[ib].get_atom_2_index()];
	    
	    of.width(3);
	    of << idx_1;
	    of.width(3);
	    of << idx_2;
	    int bond_type = 1;
	    if (bonds[ib].get_bond_type() == lig_build::bond_t::DOUBLE_BOND)
	       bond_type = 2;
	    if (bonds[ib].get_bond_type() == lig_build::bond_t::TRIPLE_BOND)
	       bond_type = 3;
	    if (bonds[ib].get_bond_type() == lig_build::bond_t::AROMATIC_BOND)
	       bond_type = 4;
	    if (bonds[ib].get_bond_type() == lig_build::bond_t::SINGLE_OR_DOUBLE)
	       bond_type = 5;
	    if (bonds[ib].get_bond_type() == lig_build::bond_t::SINGLE_OR_AROMATIC)
	       bond_type = 6;
	    if (bonds[ib].get_bond_type() == lig_build::bond_t::DOUBLE_OR_AROMATIC)
	       bond_type = 7;
	    if (bonds[ib].get_bond_type() == lig_build::bond_t::BOND_ANY)
	       bond_type = 8;
	    of.width(3);
	    of << bond_type;

	    // write the bond stereo now
	    of.width(3);
	    of << bond_stereo;
	    
	    of << "   ";
	    int bond_topology = 0; // 0 = either, 1 = ring, 2 = chain
	    of.width(3);
	    of << bond_topology;
	    int reacting_centre_status = 0;
	    of.width(3);
	    of << reacting_centre_status;
	    of << "\n";
	 }
      }

      // end
      of << "M  END\n";
      // of << "$$$$\n"; // not needed for mol2 files.
   }
   of.close();
   
   return status;
}


bool
widgeted_molecule_t::write_minimal_cif_file(const std::string &file_name) const {

   bool status = 0;

   mmdb::mmcif::File *mmCIFFile = new mmdb::mmcif::File();
   mmdb::mmcif::Data   *mmCIFData = NULL;
   mmdb::mmcif::Struct *mmCIFStruct;
   char S[2000];

   int rc = mmCIFFile->AddCIFData("comp_list");
   mmCIFData = mmCIFFile->GetCIFData("comp_list");
   rc = mmCIFData->AddStructure ("_chem_comp", mmCIFStruct);
   // std::cout << "rc on AddStructure returned " << rc << std::endl;
   if (rc!=mmdb::mmcif::CIFRC_Ok && rc!=mmdb::mmcif::CIFRC_Created)  {
      // badness!
      std::cout << "rc not mmdb::mmcif::CIFRC_Ok " << rc << std::endl;
      printf ( " **** error: attempt to retrieve Loop as a Structure.\n" );
      if (!mmCIFStruct)  {
         printf ( " **** error: mmCIFStruct is NULL - report as a bug\n" );
      }
   } else {

      status = 1;
      std::string comp_id = "DRG";
      std::string three_letter_code = "DRG";
      std::string name = "ligand";
      int n_atoms_all = get_number_of_atoms_including_hydrogens();
      int n_non_H_atoms = atoms.size();
      std::string description_level = "M";

      mmdb::mmcif::Loop *mmCIFLoop = new mmdb::mmcif::Loop; // 20100212

      rc = mmCIFData->AddLoop("_chem_comp", mmCIFLoop);
      int i=0;
      const char *s;
      mmCIFLoop->PutString(comp_id.c_str(), "comp_id", i);
      s = three_letter_code.c_str();
      mmCIFLoop->PutString(s, "three_letter_code", i);
      s = name.c_str();
      mmCIFLoop->PutString(s, "name", i);
      s =  group.c_str();
      mmCIFLoop->PutString(s, "group", i);
      mmCIFLoop->PutInteger(n_atoms_all, "number_atoms_all", i);
      mmCIFLoop->PutInteger(n_non_H_atoms, "number_atoms_nh", i);
      mmCIFLoop->PutString(description_level.c_str(), "description_level", i);

      std::string comp_record = "comp_list";
      mmCIFData->PutDataName(comp_record.c_str()); // 'data_' record

      // atom loop

      std::string comp_monomer_name = "comp_";
      comp_monomer_name += comp_id;
      rc = mmCIFFile->AddCIFData(comp_monomer_name.c_str());
      mmCIFData = mmCIFFile->GetCIFData(comp_monomer_name.c_str());
      rc = mmCIFData->AddLoop("_chem_comp_atom", mmCIFLoop);

      if (rc == mmdb::mmcif::CIFRC_Ok || rc == mmdb::mmcif::CIFRC_Created) {
         for (unsigned int i=0; i<atoms.size(); i++) {

            mmCIFLoop->PutString(comp_id.c_str(), "comp_id", i);

            std::string ss = atoms[i].get_atom_id().c_str();
            mmCIFLoop->PutString(ss.c_str(), "atom_id", i);

            ss = atoms[i].element;
            mmCIFLoop->PutString(ss.c_str(), "type_symbol", i);
         }
      }

      // bond loop

      rc = mmCIFData->AddLoop("_chem_comp_bond", mmCIFLoop);
      if (rc == mmdb::mmcif::CIFRC_Ok || rc == mmdb::mmcif::CIFRC_Created) {
         // std::cout << " number of bonds: " << bond_restraint.size() << std::endl;
         for (unsigned int i=0; i<bonds.size(); i++) {
            // std::cout << "ading bond number " << i << std::endl;
            mmCIFLoop->PutString(comp_id.c_str(), "comp_id", i);
            std::string atom_name_1 = atoms[bonds[i].get_atom_1_index()].get_atom_id();
            mmCIFLoop->PutString(atom_name_1.c_str(), "atom_id_1", i);
            std::string atom_name_2 = atoms[bonds[i].get_atom_2_index()].get_atom_id();
            mmCIFLoop->PutString(atom_name_2.c_str(), "atom_id_1", i);

            std::string bond_type = "single";
            if (bonds[i].get_bond_type() == lig_build::bond_t::DOUBLE_BOND)
               bond_type = "double";
            if (bonds[i].get_bond_type() == lig_build::bond_t::TRIPLE_BOND)
               bond_type = "triple";
            mmCIFLoop->PutString(bond_type.c_str(), "type", i);
         }
      }
      mmCIFFile->WriteMMCIFFile(file_name.c_str());
   }
   delete mmCIFFile; // deletes all its attributes too.
   return status;
}

// void
// widgeted_molecule_t::add_solvent_accesibilities_to_atoms(std::vector<solvent_accessible_atom_t> solvent_accessible_atoms) {

//    for (unsigned int i=0; i<atoms.size(); i++) { 
//       for (unsigned int j=0; j<solvent_accessible_atoms.size(); j++) { 
// 	 if (solvent_accessible_atoms[j].atom_name == atoms[i].get_atom_name()) {
// 	    atoms[i].add_solvent_accessibility(solvent_accessible_atoms[j].solvent_accessibility);
// 	 }
//       }
//    }

// } 



// return negative if not solvent accessibility available.
double
widgeted_molecule_t::get_solvent_accessibility(const clipper::Coord_orth &pt,
					       const std::vector<solvent_accessible_atom_t> &sav) const {

   double sa = -1;
   double dsq = 0.03 * 0.03;
   
   for (unsigned int i=0; i<sav.size(); i++) { 
      if ((sav[i].pt - pt).lengthsq() < dsq) {
	 sa = sav[i].solvent_accessibility;
	 break;
      }
   }

   return sa;
}

std::string
widgeted_molecule_t::get_atom_name(const clipper::Coord_orth &pt, mmdb::Manager *mol) const {

   std::string atom_name;
   double close_2 = 0.01 * 0.01;

   if (! mol) {
      // we don't need to be always told this - this *is* the case if
      // we read a straight mol file.
      // 
      // std::cout << "Null molecule in get_atom_name " << std::endl;

   } else { 

      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::Residue *residue_p;
	 mmdb::Atom *at;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    int n_atoms = residue_p->GetNumberOfAtoms();
	    for (int iat=0; iat<n_atoms; iat++) {
	       at = residue_p->GetAtom(iat);
	       clipper::Coord_orth atom_pos(at->x, at->y, at->z);
	       double d = (atom_pos-pt).lengthsq();
	       if (d < close_2) {
		  atom_name = at->name;
		  break;
	       } 
	    }
	 }
      }
   }
   return atom_name;
} 

// solvent accessibility and bash distances in fact.
// 
void
widgeted_molecule_t::map_solvent_accessibilities_to_atoms(std::vector<solvent_accessible_atom_t> solvent_accessible_atoms) {

   if (0)
      std::cout << "in map_solvent_accessibilities_to_atoms() " << atoms.size() << " atoms "
		<< "and " << solvent_accessible_atoms.size() << " solvent atoms "
		<< std::endl;
   
   for (unsigned int i=0; i<atoms.size(); i++) {
      for (unsigned int j=0; j<solvent_accessible_atoms.size(); j++) {
	 if (0) 
	    std::cout << "   for solvent accessibility" << i << "  " << j
		      << "  mol atom name :"
		      << atoms[i].get_atom_name()
		      <<  ": :" << solvent_accessible_atoms[j].atom_name << ":" << std::endl;
	 if (atoms[i].get_atom_name() == solvent_accessible_atoms[j].atom_name) {
	    atoms[i].add_solvent_accessibility(solvent_accessible_atoms[j].solvent_accessibility);
	    if (0)
	       std::cout << "transfering " << solvent_accessible_atoms[j].bash_distances.size()
			 << " bash distances for atom " << atoms[i].get_atom_name() << std::endl;
	    atoms[i].bash_distances = solvent_accessible_atoms[j].bash_distances;
	    break;
	 } 
      }
   }
} 

// can throw a runtime_error exception.
// 
lig_build::pos_t
widgeted_molecule_t::get_atom_canvas_position(const std::string &atom_name) const {

   lig_build::pos_t p;
   bool ifound = 0;
   for (unsigned int i=0; i<atoms.size(); i++) { 
      if (atoms[i].get_atom_name() == atom_name) {
	 p = atoms[i].atom_position;
	 return p;
      }
   }

   if (! ifound) {
      std::string mess = "No atom name \"";
      mess += atom_name;
      mess += "\" found in ligand";
      throw std::runtime_error(mess);
   }
   return p;
}


// // can throw an exception (no atoms) return pseudo points top-left
// // (small small) bottom-right (high high).
// //
// std::pair<lig_build::pos_t, lig_build::pos_t>
// widgeted_molecule_t::ligand_extents() const {

//    lig_build::pos_t top_left;
//    lig_build::pos_t bottom_right;

//    double mol_min_x =  9999999;
//    double mol_max_x = -9999999;
//    double mol_min_y =  9999999;
//    double mol_max_y = -9999999;

//    if (atoms.size()) { 
      
//       for (unsigned int iat=0; iat<atoms.size(); iat++) {
// 	 if (atoms[iat].atom_position.x > mol_max_x)
// 	    mol_max_x = atoms[iat].atom_position.x;
// 	 if (atoms[iat].atom_position.x < mol_min_x)
// 	    mol_min_x = atoms[iat].atom_position.x;
// 	 if (atoms[iat].atom_position.y > mol_max_y)
// 	    mol_max_y = atoms[iat].atom_position.y;
// 	 if (atoms[iat].atom_position.y < mol_min_y)
// 	    mol_min_y = atoms[iat].atom_position.y;
//       }
//       top_left     = lig_build::pos_t(mol_min_x, mol_min_y);
//       bottom_right = lig_build::pos_t(mol_max_x, mol_max_y);
//    } else {
//       std::string mess = "WARNING:: no atoms in ligand_extents()";
//       throw std::runtime_error(mess);
//    }
   
//    return std::pair<lig_build::pos_t, lig_build::pos_t> (top_left, bottom_right);

// }



// int
// widgeted_molecule_t::n_open_bonds() const {

//    int n_bonds = 0;
//    for (unsigned int i=0; i<bonds.size(); i++) { 
//       if (! bonds[i].is_closed())
// 	 n_bonds++;
//    }
//    return n_bonds;
// }


// bool
// widgeted_molecule_t::is_close_to_non_last_atom(const lig_build::pos_t &test_pos) const {

//    bool close = 0;
//    int n_atoms_for_test = atoms.size() - 1; 
//    for (int iat=0; iat<n_atoms_for_test; iat++) {
//       if (! atoms[iat].is_closed()) {
// 	 if (atoms[iat].atom_position.near_point(test_pos, 2.1)) {
// 	    close = 1;
// 	    break;
// 	 }
//       }
//    }
//    return close;
// }

void
widgeted_molecule_t::delete_hydrogens(GooCanvasItem *root) {

   for (unsigned int iat=0; iat<atoms.size(); iat++) {

   // for (unsigned int iat=0; iat<28; iat++) {
      if (atoms[iat].element == "H") {
	 // std::cout << "closing atom number " << iat << std::endl;
	 close_atom(iat, root);
      } 
   } 
}

// X or Y
void
widgeted_molecule_t::flip(int axis) {

   std::pair<double, lig_build::pos_t> s_c = current_scale_and_centre();
   const lig_build::pos_t &centre = s_c.second;
   for (unsigned int iat=0; iat<atoms.size(); iat++) { 
      if (! atoms[iat].is_closed()) {
	 if (axis == X_AXIS) {
	    double y_new = 2 * centre.y - atoms[iat].atom_position.y;
	    atoms[iat].atom_position.y = y_new;
	 }
	 if (axis == Y_AXIS) {
	    double x_new = 2 * centre.x - atoms[iat].atom_position.x;
	    atoms[iat].atom_position.x = x_new;
	 }
      }
   }

   for (unsigned int ibond=0; ibond<bonds.size(); ibond++) {
      lig_build::bond_t &bond = bonds[ibond];
      if (! bond.is_closed()) {

	 // when we invert the positions of the atom, then the centre info
	 // (for double bonds) becomes invalid, so we need to reset the bond
	 // to a copy of the bond without centre info.
	 //
	 if (bond.get_bond_type() == lig_build::bond_t::IN_BOND ||
	     bond.get_bond_type() == lig_build::bond_t::OUT_BOND) {
	    lig_build::bond_t::bond_type_t other_dir = lig_build::bond_t::IN_BOND;
	    if (bond.get_bond_type() == lig_build::bond_t::IN_BOND)
	       other_dir = lig_build::bond_t::OUT_BOND;
	    if (bond.get_bond_type() == lig_build::bond_t::OUT_BOND)
	       other_dir = lig_build::bond_t::IN_BOND;
	    bond.set_bond_type(other_dir);
	 }

	 bool shorten_first  = false;
	 bool shorten_second = false;
	 if (bond.get_bond_type() == lig_build::bond_t::DOUBLE_BOND) {

	    std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> >
	       other_connections_to_first_atom =
	       make_other_connections_to_first_atom_info(ibond);
	    std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> >
	       other_connections_to_second_atom =
	       make_other_connections_to_second_atom_info(ibond);
	    widgeted_bond_t new_bond(bond.get_atom_1_index(), bond.get_atom_2_index(),
				     atoms[bond.get_atom_1_index()], atoms[bond.get_atom_2_index()],
				     shorten_first, shorten_second,
				     lig_build::bond_t::DOUBLE_BOND,
				     other_connections_to_first_atom,
				     other_connections_to_second_atom,
				     NULL);
	    bonds[ibond] = new_bond;
	 }
      }
   }
   assign_ring_centres(true);
}

void
widgeted_molecule_t::rotate_z(double angle) {

   std::pair<double, lig_build::pos_t> s_c = current_scale_and_centre();
   const lig_build::pos_t &centre = s_c.second;
   for (unsigned int iat=0; iat<atoms.size(); iat++) {
      if (! atoms[iat].is_closed()) {
	 atoms[iat].atom_position = atoms[iat].atom_position.rotate_about(centre, angle);
      }
   }
   assign_ring_centres(true);
}


topological_equivalence_t::topological_equivalence_t(const std::vector<widgeted_atom_t> &atoms_in,
						     const std::vector<widgeted_bond_t> &bonds_in) {

   atoms = atoms_in;
   bonds = bonds_in;

   unique.resize(atoms.size(), false);

   for (unsigned int iat=0; iat<atoms.size(); iat++)
      atom_map[atoms[iat].element].push_back(iat);
      

   // The Morgan algorithm - or something like it.
   
   std::vector<long int> prev_eqv = assign_initial_topo_indices();
   std::vector<long int> curr_eqv;
   int round = 0;

   bool eq_changed = true; // fake initial value
   bool uniques_assigned = assign_uniques(prev_eqv);

   while (eq_changed || uniques_assigned) {

      std::cout << "::::::::::::::::::::::::::: round: "
		<< round << " ::::::::::::::::::::::::::::::::::::" << std::endl;

      curr_eqv = assign_topo_indices(prev_eqv, round);

      // should we do this only if eq_changed is true?
      uniques_assigned = assign_uniques(curr_eqv);
      
      eq_changed = continue_ec_calculations_p(curr_eqv, prev_eqv);
      std::cout << ":::::::::::: eq_changed: " << eq_changed
		<< "    uniques_assigned: " << uniques_assigned<< std::endl;

      // for next round
      round++;
      prev_eqv = curr_eqv;
   }

   // OK, when we get here, unique is set and the curr_eqv tells us
   // what are the equivalent atoms.
   //

   // EC numbers:
   std::cout << "--- final EC numbers -------- " << std::endl;
   for (unsigned int iec=0; iec<prev_eqv.size(); iec++)
      std::cout << "  atom_index: " << iec << " ec-final: "
		<< prev_eqv[iec] << std::endl;

   std::cout << "--- final uniques -------- " << std::endl;
   for (unsigned int iun=0; iun<unique.size(); iun++) {
      if (unique[iun])
	 std::cout << "  atom index " << iun << " marked as unique" << std::endl;
   }
   
   std::cout << "--- final equivalents -------- " << std::endl;
   std::map<std::string, std::vector<int> >::const_iterator it;
   for (it=atom_map.begin(); it!=atom_map.end(); it++) {

      // it->second points to a vector of atoms.  If they have the
      // same prev_eqv value then they are equivalent
      
      for (unsigned int iat=0; iat<it->second.size(); iat++) {

	 // is there another atom that matches the prev_eqv of this atom
	 int atom_current_equiv = prev_eqv[it->second[iat]];
	 for (unsigned int ivi=0; ivi<it->second.size(); ivi++) {
	    if (!unique[it->second[ivi]] && !unique[it->second[iat]]) { 
	       if (it->second[ivi] != it->second[iat]) {
		  if (prev_eqv[it->second[ivi]] == atom_current_equiv) {
		     std::cout << "  atom index " << it->second[ivi]
			       << " equivalent to atom index " << it->second[iat]
			       << "  " << prev_eqv[it->second[ivi]]
			       << std::endl;
		  }
	       }
	    }
	 }
      }
   }
   std::cout << "--------------------- " << std::endl;

   assign_invariant_sequence_number(curr_eqv);
} 

// Have more atoms become topologically equivalent by this iteration?
//
// (If no, then end the loop in the calling function)
// 
bool
topological_equivalence_t::continue_ec_calculations_p(const std::vector<long int> &curr_eqv,
						      const std::vector<long int> &prev_eqv) {

   int ec_curr = n_extended_connectivity(curr_eqv);
   int ec_prev = n_extended_connectivity(prev_eqv);

   std::cout << "in continue_ec_calculations_p: comparing current ec " << ec_curr
	     << " with previous ec " << ec_prev << std::endl;

   // we want to continue if the number of ec from current is more
   // than from previous round.
   
   if (ec_curr <= ec_prev)
      return false;
   else
      return true;

}

// Have more atoms become topologically equivalent by this iteration?
//
// (If no, then end the loop in the calling function)
// 
bool
topological_equivalence_t::identified_unique_p(const std::vector<long int> &curr_eqv, 
					       const std::vector<long int> &prev_eqv) {

   bool r = 0; // not changed.

   // sanity checks
   // 
   if (atoms.size() == 0) {
      std::cout << "identified_unique_p: no atoms" << std::endl;
      return r;
   }
   if (curr_eqv.size() == 0) {
      std::cout << "identified_unique_p: no current equivalences" << std::endl;
      return r;
   }
   if (curr_eqv.size() == 0) {
      std::cout << "identified_unique_p: no previous equivalences" << std::endl;
      return r;
   }
   if (curr_eqv.size() == 0) {
      std::cout << "identified_unique_p: no previous equivalences and current equivalences"
		<< " not of equal size"  << std::endl;
      return r;
   }
      
   std::map<std::string, std::vector<int> >::const_iterator it;
   for (it=atom_map.begin(); it!=atom_map.end(); it++) {
      
      // it points to a vector of atoms that have the same element.
      // Do they have new topological equivalence?
      //
      // We test that by looking at the number of different
      // topo_indices for the set of atoms (with the same element).
      //
      std::map<long int, std::vector<int> > topo_indices_curr;
      std::map<long int, std::vector<int> > topo_indices_prev;

      for (unsigned int ivi=0; ivi<it->second.size(); ivi++) {
	 if (! unique[it->second[ivi]]) { 
	    topo_indices_prev[curr_eqv[it->second[ivi]]].push_back(it->second[ivi]);
	    topo_indices_curr[prev_eqv[it->second[ivi]]].push_back(it->second[ivi]);
	    std::cout << "for topo_indices_curr added " << it->second[ivi]
		      << " for topo map key " << prev_eqv[it->second[ivi]] << std::endl;
	 } else {
	    std::cout << "topo_indices setting blocked because "
		      << it->second[ivi] << " was marked as unique"
		      << std::endl; 
	 }
      }

      // now mark the uniques
      std::map<long int, std::vector<int> >::const_iterator it_topo;

      // debug
//       for (it_topo=topo_indices_curr.begin(); it_topo!=topo_indices_curr.end(); it_topo++) {
// 	 std::cout << "debug:: key: " << it_topo->first << std::endl;
// 	 for (unsigned int jj=0; jj<it_topo->second.size(); jj++)
// 	    std::cout << "debug::  vector values   " << it_topo->second[jj] << std::endl;
//       }
      
      for (it_topo=topo_indices_curr.begin(); it_topo!=topo_indices_curr.end(); it_topo++) {
	 std::cout << "in marking uniques, it->second has size " << it_topo->second.size()
		   << std::endl;
	 if (it_topo->second.size() == 1) {
	    int unique_index = it_topo->second[0];
	    std::cout << "........ marking " << unique_index << " as unique "
		      << std::endl;
	    unique[unique_index] = true;
	    // r = 1; // we make a unique
	 }
      }

       std::cout << "for equivalence changed test: comparing  n_prev: "
		<< topo_indices_prev.size() << " and n_curr: " << topo_indices_curr.size()
 		<< std::endl;
	 
       if (topo_indices_curr.size() != topo_indices_prev.size()) {
	 r = 1;
 	 break;
       }
   } 
   return r; 
}


std::vector<long int>
topological_equivalence_t::assign_initial_topo_indices() {
   
   std::vector<long int> r(atoms.size());

   for (unsigned int ibond=0; ibond<bonds.size(); ibond++) {
      r[bonds[ibond].get_atom_1_index()]++;
      r[bonds[ibond].get_atom_2_index()]++;
   }

   for (unsigned int i=0; i<r.size(); i++)
      std::cout << "initial topo index: " << i << "  " << r[i] << std::endl;
      
   return r;
}

std::vector<long int>
topological_equivalence_t::assign_topo_indices(const std::vector<long int> &prev_eqv,
					       int round) {

   std::vector<long int> r(prev_eqv.size(), 0);

   // which atoms are bonded to which other atoms?
   std::vector<std::vector<int> > atom_bonds(atoms.size()); // atom index vector
   
   for (unsigned int ibond=0; ibond<bonds.size(); ibond++) {
      atom_bonds[bonds[ibond].get_atom_1_index()].push_back(bonds[ibond].get_atom_2_index());
      atom_bonds[bonds[ibond].get_atom_2_index()].push_back(bonds[ibond].get_atom_1_index());
   }

   for (unsigned int iat=0; iat<atoms.size(); iat++) {
      for (unsigned int ii=0; ii<atom_bonds[iat].size(); ii++) {
// 	 std::cout << "   index iat: "<< iat << " ii: " << ii << " of " << atom_bonds[iat].size()
// 		   << "    adding prev_eqv[" << atom_bonds[iat][ii] << "]="
// 		   << prev_eqv[atom_bonds[iat][ii]] << " to "
// 		   << r[iat] << std::endl;
	 r[iat] += prev_eqv[atom_bonds[iat][ii]];
      }
   }
   for (unsigned int i=0; i<r.size(); i++)
      std::cout << "   round " << round << " assigned topo index: "
		<< i << "  " << r[i] << std::endl;
   return r;
}

bool
topological_equivalence_t::assign_uniques(const std::vector<long int> &extended_connectivity) {

   bool r = 0;

   std::map<std::string, std::vector<int> >::const_iterator it;
   for (it=atom_map.begin(); it!=atom_map.end(); it++) {
      
      // it points to a vector of atoms indices that have the same element.
      //
      std::map<long int, std::vector<int> > topo_indices;
   
      for (unsigned int ivi=0; ivi<it->second.size(); ivi++) {
	 if (! unique[it->second[ivi]]) {
	    topo_indices[extended_connectivity[it->second[ivi]]].push_back(it->second[ivi]);
	 }
      }

      // were there any vectors that had size 1?
      std::map<long int, std::vector<int> >::const_iterator it_topo;
      for (it_topo=topo_indices.begin(); it_topo!=topo_indices.end(); it_topo++) {
	 if (it_topo->second.size() == 1) {
	    int unique_index = it_topo->second[0];
	    std::cout << "........ marking " << unique_index << " as unique "
		      << std::endl;
	    unique[unique_index] = true;
	    r = 1;
	 }
      }
   }
   return r; 
} 

// there is an EC value for each atom.  How many different EC values are there?
int
topological_equivalence_t::n_extended_connectivity(const std::vector<long int> &extended_connectivity) const {

   int r = 0; 

   std::map<std::string, std::vector<int> >::const_iterator it;
   for (it=atom_map.begin(); it!=atom_map.end(); it++) {
      
      // it points to a vector of atoms indices that have the same element.
      //
      std::map<long int, std::vector<int> > topo_indices;

      for (unsigned int ivi=0; ivi<it->second.size(); ivi++) {
	 if (! unique[it->second[ivi]]) {
	    topo_indices[extended_connectivity[it->second[ivi]]].push_back(it->second[ivi]);
	 }
      }
      r += topo_indices.size();
   }

   return r;
   
} 

// return a vector the same size as atoms, with the invariant sequence numbers
void
topological_equivalence_t::assign_invariant_sequence_number(const std::vector<long int> &extended_connectivity) {

   isn.resize(atoms.size(), 0);
   std::map<long int, std::vector<int> > topo_indices;
   std::map<long int, std::vector<int> >::const_reverse_iterator it_topo;
   std::map<long int, std::vector<int> >::const_reverse_iterator start_it_topo  = topo_indices.rbegin();
   std::map<long int, std::vector<int> >::const_reverse_iterator end_it_topo = topo_indices.rend();

   for (unsigned int iat=0; iat<extended_connectivity.size(); iat++)
      topo_indices[extended_connectivity[iat]].push_back(iat);

   // recall that std::maps are auto-sorted (forward low to high)

   
   // so now iterate through topo_indices, hight to low
   // 
   int next_index = 1;
   for (it_topo=start_it_topo; it_topo!=end_it_topo; it_topo++) { 

      const std::vector<int> &atom_index_deepest = it_topo->second;
      unsigned int first_atom_index(atom_index_deepest[0]);
   
      if (atom_index_deepest.size() == 1) {

	 bool done = mark_isn(first_atom_index, next_index);
	 if (done)
	    next_index++;
      
      } else {

	 std::vector<std::pair<int, int> > connected_indices;
	 for (unsigned int iat=0; iat<atom_index_deepest.size(); iat++) {
	    std::pair<int, int> p(atom_index_deepest[iat], lig_build::bond_t::BOND_UNDEFINED);
	    connected_indices.push_back(p);
	 }
	 next_index = assign_invariant_sequence_number(extended_connectivity,
						       connected_indices, next_index);
      } 
      
      // what atoms are connected to first_atom_index?
      // make a vector of them.
      std::vector<std::pair<int, int> > connected_indices;
      for (unsigned int ibond=0; ibond<bonds.size(); ibond++) { 
	 if (bonds[ibond].get_atom_1_index() == first_atom_index) { 
	    std::pair<int, int> p(bonds[ibond].get_atom_2_index(), ibond);
	    connected_indices.push_back(p);
	 }
	 if (bonds[ibond].get_atom_2_index() == first_atom_index) {
	    std::pair<int, int> p(bonds[ibond].get_atom_1_index(), ibond);
	    connected_indices.push_back(p);
	 }
      }
      next_index = assign_invariant_sequence_number(extended_connectivity,
						    connected_indices,
						    next_index);
   }

   std::cout << "-------------- isns ------------------" << std::endl;
   for (unsigned int i_isn=0; i_isn<isn.size(); i_isn++) { 
      std::cout << "   atom-index:  " << i_isn << " isn: " << isn[i_isn] << std::endl;
   }
}

// assign the isns to atom indices of atom_index
// 
int
topological_equivalence_t::assign_invariant_sequence_number(const std::vector<long int> &extended_connectivity,
							    const std::vector<std::pair<int, int> > &atom_and_bond_index,
							    int next_index) {


   // which of the atom indices have the highest ec?  Sort them
   //
   std::map<long int, std::vector<std::pair<int, int> > > ec_map;
   // fill the map
   for (unsigned int ivi=0; ivi<atom_and_bond_index.size(); ivi++) {
      std::pair<int, int> ai_pair = atom_and_bond_index[ivi];
      ec_map[extended_connectivity[ai_pair.first]].push_back(ai_pair);
   }

   // ec_map is sorted lowest first, run through that list backwards
   // (incrementing the counter, oh yes).

   std::map<long int, std::vector<std::pair<int, int> > >::const_reverse_iterator it_ec;
   std::map<long int, std::vector<std::pair<int, int> > >::const_reverse_iterator start_val = ec_map.rbegin();
   std::map<long int, std::vector<std::pair<int, int> > >::const_reverse_iterator end_val   = ec_map.rend();

   for (it_ec=start_val; it_ec!=end_val; it_ec++) {
      const std::vector<std::pair<int, int> > &atom_indices_inner = it_ec->second;

      if (atom_indices_inner.size() == 1) {
	 bool done = mark_isn(atom_indices_inner[0].first, next_index);
	 if (done)
	    next_index++;
      } else { 
	 // botheration, we have to assign equivalent atoms.  Just
	 // do it "randomly", i.e. by order for now.
	 //
	 // Should test on bond order
	 for (unsigned int ivi=0; ivi<atom_indices_inner.size(); ivi++) {
	    // use the bond?
	    const widgeted_bond_t &bond = bonds[atom_indices_inner[ivi].second];
	    bool done = mark_isn(atom_indices_inner[ivi].first, next_index);
	    if (done)
	       next_index++;
	 } 
      }
   }
   return next_index;
   
}

bool
topological_equivalence_t::atoms_have_unassigned_isn_p() const {

   int r = 0;  // non unassigned initially.
   
   for (unsigned int i=0; i<isn.size(); i++) { 
      if (!isn[i]) {
	 r = 0;
	 break;
      }
   }
   return r;
} 


bool
topological_equivalence_t::mark_isn(int atom_index, int i_s_n) {

   bool done = false;
   if (isn[atom_index] == 0) { 
      isn[atom_index] = i_s_n;
      done = true;
//       std::cout << "        marked " << atom_index << " as invarient-sequence-number "
// 		<< i_s_n << std::endl;
   } else {
//       std::cout << "      atom with index " << atom_index << " already given isn "
// 		<< isn[atom_index] << " failed to slot in " << i_s_n << std::endl;
   } 
   return done;
} 


std::vector<std::string>
topological_equivalence_t::chiral_centres() const {

   std::vector<std::string> v;

   return v;

} 

// Return a list of atom indices that are connected to 3 or 4 other
// atoms (return the indices of those other atoms too.
// 
std::vector<std::pair<int, std::vector<int> > >
topological_equivalence_t::tetrahedral_atoms() const {

   std::vector<std::pair<int, std::vector<int> > > v;
   
   return v;
}
