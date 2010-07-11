/* lbg/wmolecule.cc
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

#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "lbg.hh"

// Don't forget that this function will be used in
// render_from_molecule, which will add canvas item.
//
// Perhaps this should be in the base class, as it doesn't use widgets
// at all.
// 
widgeted_molecule_t::widgeted_molecule_t(const lig_build::molfile_molecule_t &mol_in,
					 CMMDBManager *pdb_mol) { 


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

      std::cout << "::::::::::::::: y stats: extents " << mol_in_min_y << " "
		<< mol_in_max_y << " centre correction: " << centre_correction.y
		<< std::endl;

   
      for (unsigned int iat=0; iat<mol_in.atoms.size(); iat++) {

	 clipper::Coord_orth pt = mol_in.atoms[iat].atom_position;
	 lig_build::pos_t pos = input_coords_to_canvas_coords(pt);
	 std::string element = mol_in.atoms[iat].element;
	 int charge = 0;
	 GooCanvasItem *ci = NULL;
	 widgeted_atom_t at(pos, element, charge, ci);

	 std::string atom_name = get_atom_name(pt, pdb_mol);
	 // and if atom_name is not "", then set atom name of at.
	 if (atom_name != "") {
	    if (0) 
	       std::cout << "Hoorah setting atom at " << pt.format() << " with name :"
			 << atom_name << ":" << std::endl;
	    at.set_atom_name(atom_name);
	 } else { 
	    std::cout << "Boo!! atom at " << pt.format() << " with no name"
		      << std::endl;
	 } 
	 
	 if (0)
	    std::cout << "Element " << element << " at " << pos << " "
		      << mol_in_min_y << " "
		      << mol_in_max_y << std::endl;
	 atoms.push_back(at);
      }

      for (unsigned int ib=0; ib<mol_in.bonds.size(); ib++) {
	 int index_1 = mol_in.bonds[ib].index_1;
	 int index_2 = mol_in.bonds[ib].index_2;
	 lig_build::bond_t::bond_type_t bt = mol_in.bonds[ib].bond_type;
	 GooCanvasItem *ci = NULL;
	 widgeted_bond_t bond(index_1, index_2, atoms[index_1], atoms[index_2], bt, ci);
	 bonds.push_back(bond);
      }


      // add ring centres
      //
      bool debug = 0;
      for (unsigned int ib=0; ib<bonds.size(); ib++) {
	 if (! bonds[ib].have_centre_pos()) { 
	    int atom_index = bonds[ib].get_atom_1_index();
	    int atom_index_other = bonds[ib].get_atom_2_index();
	    if (debug) 
	       std::cout << "=============== checking ring for atom index "
			 << atom_index << " ===============" << std::endl;
	    // path must pass through atom_index_other
	    std::pair<bool, std::vector<int> > found =
	       found_self_through_bonds(atom_index, atom_index_other);
	    if (debug) { 
	       std::cout << "-- constructor of widgeted_bond_t atom " << atom_index
			 << " other bond index (not tested) " << bonds[ib].get_atom_2_index()
			 << ", found status "
			 << found.first << " (";
	       for (unsigned int i=0; i<found.second.size(); i++)
		  std::cout << found.second[i] << " ";
	       std::cout << ")\n"
			 << std::endl;
	    }
	    if (found.first) {
	       if (found.second.size() > 0) { 
		  lig_build::pos_t centre_pos_sum;
		  std::string centre_pos_atoms_string;
		  for (unsigned int i_ring_atom_list=0;
		       i_ring_atom_list<found.second.size();
		       i_ring_atom_list++) {
		     centre_pos_sum += atoms[found.second[i_ring_atom_list]].atom_position;
		     // these for debugging.
		     centre_pos_atoms_string += coot::util::int_to_string(found.second[i_ring_atom_list]);
		     centre_pos_atoms_string += " ";
		  }
		  lig_build::pos_t centre_pos = centre_pos_sum * (1.0/double(found.second.size()));
		  bonds[ib].add_centre(centre_pos);
		  if (debug)
		     std::cout << "adding centre at " << centre_pos
			       << " generated from (" << centre_pos_atoms_string << ")" 
			       << " for bond " << ib
			       << " connecting " << bonds[ib].get_atom_1_index() << " to "
			       << bonds[ib].get_atom_2_index() << std::endl;
	       }
	    }
	 }
      }
   }
}

std::pair<bool, double>
widgeted_molecule_t::get_scale_correction(const lig_build::molfile_molecule_t &mol_in) const {

   bool status = 0;
   double scale = 1.0;
   std::vector<double> bond_lengths; // process this to scale up the mol file molecule if
                  		     // needed.
   
   for (unsigned int i=0; i<mol_in.bonds.size(); i++) {
      int index_1 = mol_in.bonds[i].index_1;
      int index_2 = mol_in.bonds[i].index_2;
      if (mol_in.atoms[index_1].element != "H") { 
	 if (mol_in.atoms[index_2].element != "H") { 
	    double l =
	       clipper::Coord_orth::length(mol_in.atoms[index_1].atom_position,
					   mol_in.atoms[index_2].atom_position);
	    bond_lengths.push_back(l);
	 }
      }
   }

   if (bond_lengths.size() > 0) {
      status = 1;
      std::sort(bond_lengths.begin(), bond_lengths.end());
      int index = bond_lengths.size()/2;
      double bll = bond_lengths[index];
      scale = 1.54/bll; // sqrt(1.54) is 1.24
   }

   return std::pair<bool, double> (status, scale);
}


lig_build::pos_t
widgeted_molecule_t::input_coords_to_canvas_coords(const clipper::Coord_orth &pos_in) const {

   // convert from JME-style (top canvas atoms have big Y) to internal
   // coordinates (top canvas has small Y).

   if (0) 
      std::cout << "   scale_correction  " << scale_correction.first << " " << scale_correction.second
		<< " centre_correction " << centre_correction << " mol_in_max_y " << mol_in_max_y
		<< " mol_in_min_y " << mol_in_min_y << std::endl;

   double x =   scale_correction.second * (pos_in.x() - centre_correction.x) * SINGLE_BOND_CANVAS_LENGTH/1.3;
   double y = - scale_correction.second * (pos_in.y() - centre_correction.y) * SINGLE_BOND_CANVAS_LENGTH/1.3;

   // double y_offset = 60 + scale_correction.second * (centre_correction.y - mol_in_min_y) * 20;
   double y_offset = 110 + scale_correction.second * (centre_correction.y - mol_in_min_y) * 30;
   
   x += 300;
   y += y_offset;

   return lig_build::pos_t(x,y);
}


void
widgeted_molecule_t::translate(const lig_build::pos_t &delta) {

   for (unsigned int iat=0; iat<atoms.size(); iat++) {
      atoms[iat].atom_position += delta;
   }
   for (unsigned int ib=0; ib<bonds.size(); ib++) {
      if (bonds[ib].have_centre_pos()) {
	 bonds[ib].move_centre_pos(delta);
      }
   }
}


// simple (non-ring system) bonds.
//
// now deals with stereo_out/wedged/OUT_BOND
// 
GooCanvasItem *
widgeted_bond_t::canvas_item_for_bond(const lig_build::pos_t &pos_1_raw,
				      const lig_build::pos_t &pos_2_raw,
				      bool shorten_first,
				      bool shorten_second,
				      bond_type_t bt,
				      GooCanvasItem *root) const {

   double shorten_fraction = 0.76;
   
   lig_build::pos_t pos_1 = pos_1_raw;
   lig_build::pos_t pos_2 = pos_2_raw;

   // fraction_point() returns a point that is (say) 0.8 of the way
   // from p1 (first arg) to p2 (second arg).
   // 
   if (shorten_first)
      pos_1 = lig_build::pos_t::fraction_point(pos_2_raw, pos_1_raw, shorten_fraction);
   if (shorten_second)
      pos_2 = lig_build::pos_t::fraction_point(pos_1_raw, pos_2_raw, shorten_fraction);


   GooCanvasItem *ci = NULL;
   switch (bt) {
   case SINGLE_BOND:
      ci = goo_canvas_polyline_new_line(root,
					pos_1.x, pos_1.y,
					pos_2.x, pos_2.y,
					"stroke-color", dark,
					NULL);
      break;
   case DOUBLE_BOND:
      {
	 if (have_centre_pos()) {
	    ci = canvas_item_double_aromatic_bond(pos_1, pos_2, root);
	 } else {
	    ci = canvas_item_double_bond(pos_1, pos_2, root);
	 } 
      }
      break;
   case TRIPLE_BOND:
      { 
	 GooCanvasItem *group = goo_canvas_group_new (root, "stroke-color", dark,
						      NULL);
      
	 lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
	 lig_build::pos_t buv_90 = buv.rotate(90);
	 double small_t = 4;
	 lig_build::pos_t p1 = pos_1 + buv_90 * small_t;
	 lig_build::pos_t p2 = pos_2 + buv_90 * small_t;
	 lig_build::pos_t p3 = pos_1;
	 lig_build::pos_t p4 = pos_2;
	 lig_build::pos_t p5 = pos_1 - buv_90 * small_t;
	 lig_build::pos_t p6 = pos_2 - buv_90 * small_t;
	 GooCanvasItem *ci_1 = goo_canvas_polyline_new_line(group,
							    p1.x, p1.y,
							    p2.x, p2.y, NULL);
	 GooCanvasItem *ci_2 = goo_canvas_polyline_new_line(group,
							    p3.x, p3.y,
							    p4.x, p4.y,
							    NULL);
	 GooCanvasItem *ci_3 = goo_canvas_polyline_new_line(group,
							    p5.x, p5.y,
							    p6.x, p6.y,
							    NULL);
	 ci = group;
      }
      break;
   case IN_BOND:
      ci = make_wedge_bond_item(pos_1, pos_2, bt, root);
      break;
   case OUT_BOND:
      ci = make_wedge_bond_item(pos_1, pos_2, bt, root);
      break;
   case BOND_UNDEFINED:
      break;
   }
   return ci;
}


GooCanvasItem *
widgeted_bond_t::canvas_item_double_bond(const lig_build::pos_t &pos_1,
					 const lig_build::pos_t &pos_2,
					 GooCanvasItem *root) const { 
	 
   GooCanvasItem *group = goo_canvas_group_new (root, "stroke-color", dark,
						NULL);
      
   lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
   lig_build::pos_t buv_90 = buv.rotate(90);

   double small_t = 2;
   lig_build::pos_t p1 = pos_1 + buv_90 * small_t;
   lig_build::pos_t p2 = pos_2 + buv_90 * small_t;
   lig_build::pos_t p3 = pos_1 - buv_90 * small_t;
   lig_build::pos_t p4 = pos_2 - buv_90 * small_t;
   GooCanvasItem *ci_1 = goo_canvas_polyline_new_line(group,
						      p1.x, p1.y,
						      p2.x, p2.y, NULL);
   GooCanvasItem *ci_2 = goo_canvas_polyline_new_line(group,
						      p3.x, p3.y,
						      p4.x, p4.y,
						      NULL);
   return group;
}


// Note of course, that a shortened aromatic double bond is probably
// wrong.  Only carbons (AFAICS) will have aromatic double bonds.
// 
GooCanvasItem *
widgeted_bond_t::canvas_item_double_aromatic_bond(const lig_build::pos_t &pos_1,
						  const lig_build::pos_t &pos_2,
						  GooCanvasItem *root) const {

   GooCanvasItem *group = goo_canvas_group_new (root, "stroke-color", dark,
						NULL);
   GooCanvasItem *ci_1 = goo_canvas_polyline_new_line(group,
						      pos_1.x, pos_1.y,
						      pos_2.x, pos_2.y, NULL);

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

   lig_build::pos_t cutened_inner_start_point =
      lig_build::pos_t::fraction_point(inner_start_point, inner_end_point, 0.1);
   lig_build::pos_t cutened_inner_end_point =
      lig_build::pos_t::fraction_point(inner_start_point, inner_end_point, 0.9);
   
   GooCanvasItem *ci_2 =
      goo_canvas_polyline_new_line(group,
				   cutened_inner_start_point.x, 
				   cutened_inner_start_point.y, 
				   cutened_inner_end_point.x, 
				   cutened_inner_end_point.y, 
				   NULL);
   return group;
}

GooCanvasItem *
widgeted_bond_t::make_wedge_bond_item(const lig_build::pos_t &pos_1,
				      const lig_build::pos_t &pos_2,
				      const lig_build::bond_t::bond_type_t &bt,
				      GooCanvasItem *root) const {
   
   GooCanvasItem *item = NULL;

   if (bt == lig_build::bond_t::OUT_BOND)
      item = make_wedge_out_bond_item(pos_1, pos_2, root);
   if (bt == lig_build::bond_t::IN_BOND)
      item = make_wedge_in_bond_item(pos_1, pos_2, root);

   return item;
}

GooCanvasItem *
widgeted_bond_t::make_wedge_out_bond_item(const lig_build::pos_t &pos_1,
					  const lig_build::pos_t &pos_2,
					  GooCanvasItem *root) const {

   lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
   lig_build::pos_t buv_90 = buv.rotate(90);
   lig_build::pos_t short_edge_pt_1 = pos_2 + buv_90 * 4;
   lig_build::pos_t short_edge_pt_2 = pos_2 - buv_90 * 4;

   // the line width means that the sharp angle an pos_1 here results
   // in a few pixels beyond the pos_1, so artificially shorten it a
   // tiny amount.
   // 
   lig_build::pos_t sharp_point = lig_build::pos_t::fraction_point(pos_1, pos_2, 0.11);
   
   GooCanvasItem *item =
      goo_canvas_polyline_new(root, TRUE, 3,
			      sharp_point.x, sharp_point.y, 
			      short_edge_pt_1.x, short_edge_pt_1.y,
			      short_edge_pt_2.x, short_edge_pt_2.y,
			      "stroke-color", dark,
			      "fill-color", dark,
			      NULL);
   return item;
}

// pos_1 is at the sharp end.  Into the page.  A series of lines.
// 
GooCanvasItem *
widgeted_bond_t::make_wedge_in_bond_item(const lig_build::pos_t &pos_1,
					 const lig_build::pos_t &pos_2,
					 GooCanvasItem *root) const {
   
   GooCanvasItem *group = goo_canvas_group_new (root, "stroke-color", dark,
						NULL);
   lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
   lig_build::pos_t buv_90 = buv.rotate(90);
   int n_lines = 5;
   for (unsigned int i=1; i<=n_lines; i++) {
      // then centre point of the line, some way along the pos_1 -> pos_2 vector;
      double len = double(i) * 1.4;
      double frac = (double(i)- 0.3)/double(n_lines);
      lig_build::pos_t fp = lig_build::pos_t::fraction_point(pos_1, pos_2, frac);
      lig_build::pos_t p1 = fp + buv_90 * len;
      lig_build::pos_t p2 = fp - buv_90 * len;
      GooCanvasItem *ci_1 = goo_canvas_polyline_new_line(group,
							 p1.x, p1.y,
							 p2.x, p2.y, NULL);
   }
   return group;
}


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


int
widgeted_molecule_t::n_stray_atoms() const {  // unbonded atoms

   return stray_atoms().size();
   
}

std::vector<int>
widgeted_molecule_t::stray_atoms() const {

   std::vector<int> strays;
   
   bool found[atoms.size()];

   for (unsigned int i=0; i<atoms.size(); i++)
      found[i] = 0;

   for (unsigned int ib=0; ib<bonds.size(); ib++) { 
      int iat_1 = bonds[ib].get_atom_1_index();
      int iat_2 = bonds[ib].get_atom_2_index();
      if (! atoms[iat_1].is_closed())
	 found[iat_1] = 1;
      if (! atoms[iat_2].is_closed())
	 found[iat_2] = 1;
   }

   for (unsigned int i=0; i<atoms.size(); i++)
      if (! found[i])
	 strays.push_back(i);

   return strays;
}


bool
widgeted_molecule_t::operator==(const widgeted_molecule_t &mol_other) const {

   bool status = 0;

   // need to check that bonds are the same (atom indexing can be
   // different) and also stray atoms need to be checked after bonds.
   //
   int n_bond_hits = 0;

   if (mol_other.bonds.size() != bonds.size())
      return status;
   for (unsigned int ib=0; ib<bonds.size(); ib++) {
      lig_build::atom_t atom_i_1 = atoms[bonds[ib].get_atom_1_index()];
      lig_build::atom_t atom_i_2 = atoms[bonds[ib].get_atom_2_index()];
      int n_hits = 0;
      for (unsigned int jb=0; jb<mol_other.bonds.size(); jb++) {
	 lig_build::atom_t atom_j_1 = mol_other.atoms[bonds[ib].get_atom_1_index()];
	 lig_build::atom_t atom_j_2 = mol_other.atoms[bonds[ib].get_atom_2_index()];
	 if (atom_i_1 == atom_j_1)
	    if (atom_i_2 == atom_j_2)
	       n_hits++;
      }
      if (n_hits == 1) {
	 n_bond_hits++;
      }
   }

   if (n_bond_hits == bonds.size()) {

      // So the bonds were the same, now the strays...

      if (mol_other.n_stray_atoms() == n_stray_atoms()) {
	 std::vector<int> i_stray_atoms = stray_atoms();
	 std::vector<int> j_stray_atoms = mol_other.stray_atoms();
	 int n_stray_hits = 0;
	 for (unsigned int i=0; i<i_stray_atoms.size(); i++) {
	    for (unsigned int j=0; j<j_stray_atoms.size(); j++) {
	       if (atoms[i_stray_atoms[i]] == mol_other.atoms[j_stray_atoms[j]])
		  n_stray_hits++;
	    }
	 }
	 if (n_stray_hits == n_stray_atoms())
	    status = 1;
      }

   }
   return status;
}


std::vector<int>
widgeted_molecule_t::bonds_having_atom_with_atom_index(int test_atom_index) const {

   std::vector<int> v;
   std::vector<int> vb =  bond_indices_with_atom_index(test_atom_index);

   for (unsigned int iv=0; iv<vb.size(); iv++) {
      if (! bonds[vb[iv]].is_closed())
	 v.push_back(vb[iv]);
   }

   return v;
} 



// return a vector of the atom indices of unconnected atoms
//
std::vector<int>
widgeted_molecule_t::get_unconnected_atoms() const {

   std::vector<int> v;
   for (unsigned int iat=0; iat<atoms.size(); iat++) {
      if (! atoms[iat].is_closed()) { 
	 bool in_a_bond = 0;
	 for (unsigned int ib=0; ib<bonds.size(); ib++) { 
	    if (! bonds[ib].is_closed()) {
	       if (bonds[ib].get_atom_1_index() == iat)
		  in_a_bond = 1;
	       if (bonds[ib].get_atom_2_index() == iat)
		  in_a_bond = 1;
	    }
	    if (in_a_bond)
	       break;
	 }
	 if (! in_a_bond)
	    v.push_back(iat);
      }
   }
   return v;
}

// can throw an exception (no atoms)
// 
lig_build::pos_t
widgeted_molecule_t::get_ligand_centre() const {

   lig_build::pos_t centre(0,0);

   if (atoms.size() == 0) {
      std::string message("No atoms in ligand");
      throw std::runtime_error(message);
   } else {
      lig_build::pos_t centre_sum(0,0);
      for (unsigned int iat=0; iat<atoms.size(); iat++) {
	 centre_sum += atoms[iat].atom_position;
      }
      if (atoms.size() > 0)
	 centre = centre_sum * (1.0/double(atoms.size()));
   }
   return centre;
}

// can throw an exception (no atoms)
// 
lig_build::pos_t
widgeted_molecule_t::get_ring_centre(const std::vector<std::string> &ring_atom_names) const {

   lig_build::pos_t positions_sum(0,0);
   int n_found = 0;
   for (unsigned int jat=0; jat<ring_atom_names.size(); jat++) {
      for (unsigned int iat=0; iat<atoms.size(); iat++) {
	 if (ring_atom_names[jat] == atoms[iat].get_atom_name()) {
	    positions_sum += atoms[iat].atom_position;
	    n_found++;
	    break;
	 }
      }
   }
   if (n_found == 0) {
      std::string mess = "No ring atom names found in ligand!";
      throw(std::runtime_error(mess));
   }
   lig_build::pos_t centre = positions_sum * (1.0/double(n_found));
   return centre;
}

// can throw an exception (no rings with this atom)
//
lig_build::pos_t
widgeted_molecule_t::get_ring_centre(const widgeted_atom_ring_centre_info_t &atom) const {

   lig_build::pos_t position(0,0);
   bool found = 0;

   for (unsigned int ibond=0; ibond<bonds.size(); ibond++) { 
      if ((atoms[bonds[ibond].get_atom_1_index()] == atom.atom) ||
	  (atoms[bonds[ibond].get_atom_2_index()] == atom.atom)) {
	 if (bonds[ibond].have_centre_pos()) {
	    position = bonds[ibond].centre_pos();
	    found = 1;
	 }
      }
      if (found)
	 break;
   }

   if (! found) {
      std::string mess("No atom ");
      mess += atom.atom.get_atom_name();
      mess += " found to be in a ring";
      throw(std::runtime_error(mess));
   }
   return position;
}


// can throw an exception (no bonds)
//
// not const because it now caches the return value;
//
std::vector<lig_build::pos_t>
widgeted_molecule_t::get_ring_centres() {

   if (have_cached_bond_ring_centres_flag) {
      return cached_bond_ring_centres;
   } else { 
      std::vector<lig_build::pos_t> v;
      for (unsigned int ib=0; ib<bonds.size(); ib++) {
	 if (bonds[ib].have_centre_pos()) {
	    lig_build::pos_t rc = bonds[ib].centre_pos();
	    bool found = 0;
	    for (unsigned int i=0; i<v.size(); i++) {
	       if (v[i].near_point(rc, 7)) { // 7 seems good, others tested.
		  found = 1;
		  break;
	       }
	    }
	    if (! found)
	       v.push_back(rc);
	 }
      }
      //    std::cout << "get_ring_centres returns\n";
      //    for (unsigned int iv=0; iv<v.size(); iv++) {
      //       std::cout << "   "  << iv << " " << v[iv] << "\n";
      //    }
      cached_bond_ring_centres = v;
      have_cached_bond_ring_centres_flag = 1;
      return v;
   }
}



bool
widgeted_molecule_t::close_bond(int ib, GooCanvasItem *root,
				bool handle_post_delete_stray_atoms_flag) {

   // on killing a bond, an N at one end of the bond may need its atom_id
   // changed to NH or so.

   bool status = 0;
   if ((ib >= 0) && (ib<bonds.size())) {
      int ind_1 = bonds[ib].get_atom_1_index();
      int ind_2 = bonds[ib].get_atom_2_index();
      bonds[ib].close(root);
      status = 1;

      // ind_1
      std::string ele = atoms[ind_1].element;
      std::vector<int> local_bonds = bonds_having_atom_with_atom_index(ind_1);
      std::string new_atom_id = make_atom_id_by_using_bonds(ele, local_bonds);
      atoms[ind_1].update_atom_id_maybe(new_atom_id, root);
      // ind_2
      ele = atoms[ind_2].element;
      local_bonds = bonds_having_atom_with_atom_index(ind_2);
      new_atom_id = make_atom_id_by_using_bonds(ele, local_bonds);
      atoms[ind_2].update_atom_id_maybe(new_atom_id, root);
   
      if (handle_post_delete_stray_atoms_flag) {
	 // are there any atoms that now don't have bonds?
	 std::vector<int> stray_atoms = get_unconnected_atoms();
	 std::vector<int> bonds; // empty
	 // std::cout << "got " << stray_atoms.size() << " stray atoms " << std::endl;
	 if (stray_atoms.size()) {
	    for (unsigned int istray=0; istray<stray_atoms.size(); istray++) { 
	       std::string atom_id =
		  make_atom_id_by_using_bonds(atoms[stray_atoms[istray]].element, bonds);
	       atoms[stray_atoms[istray]].update_atom_id_forced(atom_id, root);
	    }
	 } 
      }
   }
   return status;
}

bool
widgeted_molecule_t::close_atom(int iat, GooCanvasItem *root) {

   std::vector<int> bds = bonds_having_atom_with_atom_index(iat);
   
   bool status = 0;
   if ((iat >= 0) && (iat<atoms.size())) {
      atoms[iat].close(root);
      status = 1;
   }

   // erase the bonds that were attached to atom iat:
   //
   // the atoms attached at the other end of the bonds may need to
   // change (say from N to NH).  Make a list of the atoms affected by
   // deleting the bond(s).

   std::vector<int> affected_neighbour_atoms;

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


// atom_index_other is the index of the other atom, the path back to
// atom_index_start must pass through atom_index_other.
// 
std::pair<bool, std::vector<int> >
widgeted_molecule_t::found_self_through_bonds(int atom_index_start,
					      int atom_index_other) const {

   std::vector<int> empty_no_pass_atoms;
   empty_no_pass_atoms.push_back(atom_index_start);
   std::pair<bool, std::vector<int> > r =
      find_bonded_atoms_with_no_pass(atom_index_start, atom_index_start, atom_index_other, empty_no_pass_atoms,
				     MAX_SEARCH_DEPTH);
   return r;
}


// Can I find start_atom_index from this_atom_index running through
// bond tree until given depth?
//
// Must pass through atom_index_other for good solution.
// 
std::pair<bool, std::vector<int> >
widgeted_molecule_t::find_bonded_atoms_with_no_pass(int start_atom_index,
						    int atom_index_other,
						    int this_atom_index,
						    const std::vector<int> &no_pass_atoms,
						    int depth) const {

   std::vector<int> atoms_bonded_to_atom_index_start;
   std::vector<int> local_no_pass_atoms = no_pass_atoms;
   if (depth == 0) {
      std::vector<int> empty;
      return std::pair<bool, std::vector<int> > (0, empty);
   } else {

      // get a list of all the atoms that are bonded to this atom not
      // in the no_pass list
      
      for (unsigned int i=0; i<bonds.size(); i++) { 
	 if (bonds[i].get_atom_1_index() == this_atom_index) {
	    int idx = bonds[i].get_atom_2_index();
	    if (idx == start_atom_index)
	       if (depth < (MAX_SEARCH_DEPTH-1)) {
		  if (member(atom_index_other, local_no_pass_atoms)) {
		     // if this_atom_index was not alread in
		     // local_no_pass_atoms, then add it.
		     bool ifound = 0; 
		     for (unsigned int ilnp=0; ilnp<local_no_pass_atoms.size(); ilnp++) {
			if (local_no_pass_atoms[ilnp] == this_atom_index) {
			   ifound = 1;
			   break;
			}
		     }
		     if (! ifound)
			local_no_pass_atoms.push_back(this_atom_index);
		     // debug_pass_atoms(start_atom_index, this_atom_index, depth, local_no_pass_atoms);
		     return std::pair<bool, std::vector<int> > (1, local_no_pass_atoms);
		  }
	       }
	    if (! member(idx, local_no_pass_atoms)) {
	       atoms_bonded_to_atom_index_start.push_back(idx);
	       // add this_atom_index to local_no_pass_atoms if it is not already a member.
	       bool found_this_atom_in_no_pass_atoms = 0;
	       for (unsigned int inp=0; inp<local_no_pass_atoms.size(); inp++) { 
		  if (local_no_pass_atoms[inp] == this_atom_index) {
		     found_this_atom_in_no_pass_atoms = 1;
		     break;
		  }
	       }
	       if (! found_this_atom_in_no_pass_atoms)
		  local_no_pass_atoms.push_back(this_atom_index);
	    } 
	 }
	 if (bonds[i].get_atom_2_index() == this_atom_index) {
	    int idx = bonds[i].get_atom_1_index();
	    if (idx == start_atom_index)
	       if (depth < (MAX_SEARCH_DEPTH-1)) {
		  if (member(atom_index_other, local_no_pass_atoms)) { 
		     local_no_pass_atoms.push_back(this_atom_index);
		     // if this_atom_index was not alread in
		     // local_no_pass_atoms, then add it.
		     bool ifound = 0; 
		     for (unsigned int ilnp=0; ilnp<local_no_pass_atoms.size(); ilnp++) {
			if (local_no_pass_atoms[ilnp] == this_atom_index) {
			   ifound = 1;
			   break;
			}
		     }
		     if (! ifound)
			local_no_pass_atoms.push_back(this_atom_index);
		     // debug_pass_atoms(start_atom_index, this_atom_index, depth, local_no_pass_atoms);
		     return std::pair<bool, std::vector<int> > (1, local_no_pass_atoms);
		  }
	       }
	    if (! member(idx, local_no_pass_atoms)) {
	       atoms_bonded_to_atom_index_start.push_back(idx);
	       // add this_atom_index to local_no_pass_atoms if it is not already a member.
	       bool found_this_atom_in_no_pass_atoms = 0;
	       for (unsigned int inp=0; inp<local_no_pass_atoms.size(); inp++) { 
		  if (local_no_pass_atoms[inp] == this_atom_index) {
		     found_this_atom_in_no_pass_atoms = 1;
		     break;
		  }
	       }
	       if (! found_this_atom_in_no_pass_atoms)
		  local_no_pass_atoms.push_back(this_atom_index);
	    } 
	 }
      }

      if (0) {  // debug;
	 std::cout << "     atom index " << this_atom_index << " has "
		   << atoms_bonded_to_atom_index_start.size()
		   << " connected atoms not in the no-pass list (";
	 for (unsigned int i=0; i<local_no_pass_atoms.size(); i++)
	    std::cout << local_no_pass_atoms[i] << " ";
	 std::cout << ")\n";
      }
	 
      for (unsigned int iat=0; iat<atoms_bonded_to_atom_index_start.size(); iat++) { 
	 std::pair<bool, std::vector<int> > r =
	    find_bonded_atoms_with_no_pass(start_atom_index,
					   atom_index_other,
					   atoms_bonded_to_atom_index_start[iat],
					   local_no_pass_atoms, depth-1);
	 if (r.first) {
	    // std::cout << "    passing on success...." << std::endl;
	    // debug_pass_atoms(start_atom_index, this_atom_index, depth, r.second);
	    return r;
	 }
      }
      
   } // end depth test

   if (0) { 
      std::cout << "returning 0 with this depth " << depth << " no-pass-list: (";
      for (unsigned int i=0; i<local_no_pass_atoms.size(); i++)
	 std::cout << local_no_pass_atoms[i] << " ";
      std::cout << ")\n";
   }
   std::vector<int> empty;
   return std::pair<bool, std::vector<int> > (0, empty);
}

void
widgeted_molecule_t::debug_pass_atoms(int atom_index_start, int this_atom_index, int depth,
				      const std::vector<int> &local_no_pass_atoms) const {

   std::cout << "    found atom index " << atom_index_start << " from this atom: " << this_atom_index
	     << ", at depth " << depth << " no-pass-atoms: (";
   for (unsigned int i=0; i<local_no_pass_atoms.size(); i++) { 
      std::cout << local_no_pass_atoms[i] << " ";
   }
   std::cout << ")" << std::endl;
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
      std::cout << "WARNING:: Cannot open file " << file_name << " to write molfile" << std::endl;
   } else {

      // 3 lines of header
      of << "# Molecule sketched in Liebig (part of Coot)\n";
      of << "#\n";
      of << "#\n";
      
//       timeval start_time;
//       int time_status = gettimeofday(&start_time, NULL);
//       tm *lt = localtime(start_time);
      
      // next the counts line:
      of.width(3);
      bool chiral = 0;
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
	    of << setiosflags(std::ios::fixed) << std::setw(10) << std::setprecision(4)
	       << (atoms[iat].atom_position.x - centre.x) * 1.33/SINGLE_BOND_CANVAS_LENGTH;
	    of << setiosflags(std::ios::fixed) << std::setw(10) << std::setprecision(4)
	       << (- atoms[iat].atom_position.y + centre.y) * 1.33/SINGLE_BOND_CANVAS_LENGTH;
	    double z = 0.0;
	    of << setiosflags(std::ios::fixed) << std::setw(10) << std::setprecision(4)
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
	    int charge = 0;
	    of.width(3);
	    of << charge;
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
	 if (! bonds[ib].is_closed()) { 
	    of.width(3);
	    of << mdl_atoms[bonds[ib].get_atom_1_index()];
	    of.width(3);
	    of << mdl_atoms[bonds[ib].get_atom_2_index()];
	    int bond_type = 1;
	    if (bonds[ib].get_bond_type() == lig_build::bond_t::DOUBLE_BOND)
	       bond_type = 2;
	    if (bonds[ib].get_bond_type() == lig_build::bond_t::TRIPLE_BOND)
	       bond_type = 3;
	    of.width(3);
	    of << bond_type;
	    int bond_stereo = 0; // fixme
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
      of << "$$$$\n";
 }
   return status;
}

int
widgeted_molecule_t::get_number_of_atom_including_hydrogens() const {

   return atoms.size();

}


bool
widgeted_molecule_t::write_minimal_cif_file(const std::string &file_name) const {

   bool status = 0;

   PCMMCIFFile mmCIFFile = new CMMCIFFile();

   PCMMCIFData   mmCIFData = NULL;
   PCMMCIFStruct mmCIFStruct;
   char S[2000];

   int rc = mmCIFFile->AddMMCIFData("comp_list");
   mmCIFData = mmCIFFile->GetCIFData("comp_list");
   rc = mmCIFData->AddStructure ("_chem_comp", mmCIFStruct);
   // std::cout << "rc on AddStructure returned " << rc << std::endl;
   if (rc!=CIFRC_Ok && rc!=CIFRC_Created)  {
      // badness!
      std::cout << "rc not CIFRC_Ok " << rc << std::endl;
      printf ( " **** error: attempt to retrieve Loop as a Structure.\n" );
      if (!mmCIFStruct)  {
         printf ( " **** error: mmCIFStruct is NULL - report as a bug\n" );
      }
   } else {

      status = 1;
      std::string comp_id = "DRG";
      std::string three_letter_code = "DRG";
      std::string name = "ligand";
      int n_atoms_all = get_number_of_atom_including_hydrogens();
      int n_non_H_atoms = atoms.size();
      std::string description_level = "M";

      PCMMCIFLoop mmCIFLoop = new CMMCIFLoop; // 20100212

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
      rc = mmCIFFile->AddMMCIFData(comp_monomer_name.c_str());
      mmCIFData = mmCIFFile->GetCIFData(comp_monomer_name.c_str());
      rc = mmCIFData->AddLoop("_chem_comp_atom", mmCIFLoop);

      if (rc == CIFRC_Ok || rc == CIFRC_Created) {
         for (int i=0; i<atoms.size(); i++) {

            mmCIFLoop->PutString(comp_id.c_str(), "comp_id", i);

            std::string ss = atoms[i].get_atom_id().c_str();
            mmCIFLoop->PutString(ss.c_str(), "atom_id", i);

            ss = atoms[i].element;
            mmCIFLoop->PutString(ss.c_str(), "type_symbol", i);
         }
      }

      // bond loop

      rc = mmCIFData->AddLoop("_chem_comp_bond", mmCIFLoop);
      if (rc == CIFRC_Ok || rc == CIFRC_Created) {
         // std::cout << " number of bonds: " << bond_restraint.size() << std::endl;
         for (int i=0; i<bonds.size(); i++) {
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
widgeted_molecule_t::get_atom_name(const clipper::Coord_orth &pt, CMMDBManager *mol) const {

   std::string atom_name;
   double close_2 = 0.01 * 0.01;

   if (! mol) {
      // we don't need to be always told this - this *is* the case if
      // we read a straight mol file.
      // 
      // std::cout << "Null molecule in get_atom_name " << std::endl;

   } else { 

      int imod = 1;
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 CResidue *residue_p;
	 CAtom *at;
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

   std::cout << "in map_solvent_accessibilities_to_atoms() " << atoms.size() << " atoms "
	     << "and " << solvent_accessible_atoms.size() << " solvent atoms " << std::endl;
   for (unsigned int i=0; i<atoms.size(); i++) {
      for (unsigned int j=0; j<solvent_accessible_atoms.size(); j++) {
	 if (0) 
	    std::cout << "   for solvent accessibility" << i << "  " << j << "  :"
		      << atoms[i].get_atom_name()
		      <<  ": :" << solvent_accessible_atoms[j].atom_name << ":" << std::endl;
	 if (atoms[i].get_atom_name() == solvent_accessible_atoms[j].atom_name) {
	    atoms[i].add_solvent_accessibility(solvent_accessible_atoms[j].solvent_accessibility);
	    atoms[i].bash_distances = solvent_accessible_atoms[j].bash_distances;
	    break;
	 } 
      }
   }
} 

// can throw an exception
// 
lig_build::pos_t
widgeted_molecule_t::get_atom_canvas_position(const std::string &atom_name) const {

   lig_build::pos_t p;
   bool ifound = 0;
   for (unsigned int i=0; i<atoms.size(); i++) { 
      if (atoms[i].get_atom_name() == atom_name) {
	 p = atoms[i].atom_position;
	 ifound = 1;
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


// can throw an exception (no atoms) return pseudo points top-left
// (small small) bottom-right (high high).
//
std::pair<lig_build::pos_t, lig_build::pos_t>
widgeted_molecule_t::ligand_extents() const {

   lig_build::pos_t top_left;
   lig_build::pos_t bottom_right;

   double mol_min_x =  9999999;
   double mol_max_x = -9999999;
   double mol_min_y =  9999999;
   double mol_max_y = -9999999;

   if (atoms.size()) { 
      
      for (unsigned int iat=0; iat<atoms.size(); iat++) {
	 if (atoms[iat].atom_position.x > mol_max_x)
	    mol_max_x = atoms[iat].atom_position.x;
	 if (atoms[iat].atom_position.x < mol_min_x)
	    mol_min_x = atoms[iat].atom_position.x;
	 if (atoms[iat].atom_position.y > mol_max_y)
	    mol_max_y = atoms[iat].atom_position.y;
	 if (atoms[iat].atom_position.y < mol_min_y)
	    mol_min_y = atoms[iat].atom_position.y;
      }
      top_left     = lig_build::pos_t(mol_min_x, mol_min_y);
      bottom_right = lig_build::pos_t(mol_max_x, mol_max_y);
   } else {
      std::string mess = "WARNING:: no atoms in ligand_extents()";
      throw std::runtime_error(mess);
   }
   
   return std::pair<lig_build::pos_t, lig_build::pos_t> (top_left, bottom_right);

}



int
widgeted_molecule_t::n_open_bonds() const {

   int n_bonds = 0;
   for (unsigned int i=0; i<bonds.size(); i++) { 
      if (! bonds[i].is_closed())
	 n_bonds++;
   }
   return n_bonds;
}


bool
widgeted_molecule_t::is_close_to_non_last_atom(const lig_build::pos_t &test_pos) const {

   bool close = 0;
   int n_atoms_for_test = atoms.size() - 1; 
   for (int iat=0; iat<n_atoms_for_test; iat++) {
      if (! atoms[iat].is_closed()) {
	 if (atoms[iat].atom_position.near_point(test_pos, 2.1)) {
	    close = 1;
	    break;
	 }
      }
   }
   return close;
}

