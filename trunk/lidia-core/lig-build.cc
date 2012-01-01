/* lbg/lig-build.cc
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

#include <iostream>

#include "lig-build.hh"


// is x,y (from the mouse pointer position) over this bond?
bool
lig_build::bond_t::over_bond(double x_in, double y_in,
			     const lig_build::atom_t &atom_1_at,
			     const lig_build::atom_t &atom_2_at) const {

   bool status = 0;
   lig_build::pos_t pos_in(x_in, y_in);
   for(double icp=0.25; icp<=0.75; icp+=0.1) {
      lig_build::pos_t test_pt =
	 lig_build::pos_t::fraction_point(atom_1_at.atom_position,
						    atom_2_at.atom_position, icp);
      if (test_pt.close_point(pos_in)) {
	 status = 1;
	 break;
      }
   }
   return status;
}

//       static
std::vector<std::pair<lig_build::pos_t, lig_build::pos_t> >
lig_build::pos_t::make_wedge_in_bond(const pos_t &pos_1, const pos_t &pos_2) {
   
   std::vector<std::pair<pos_t, pos_t> > lines;
   pos_t buv = (pos_2-pos_1).unit_vector();
   pos_t buv_90 = buv.rotate(90);
   int n_lines = 5;
   double bond_length = pos_t::length(pos_2, pos_1);
   double length_scale = 0.03 * bond_length;
   for (unsigned int i=1; i<=n_lines; i++) {
      // then centre point of the line, some way along the pos_1 -> pos_2 vector;
      double len = double(i) * length_scale;
      double frac = (double(i)- 0.3)/double(n_lines);
      pos_t fp = pos_t::fraction_point(pos_1, pos_2, frac);
      pos_t p1 = fp + buv_90 * len;
      pos_t p2 = fp - buv_90 * len;
      lines.push_back(std::pair<pos_t, pos_t>(p1,p2));
   }
   return lines;
}

// static
std::vector<lig_build::pos_t>
lig_build::pos_t::make_wedge_out_bond(const pos_t &pos_1, const pos_t &pos_2) {

   // double sharp_point_sharp_factor = 0.1; // for canvas
   double sharp_point_sharp_factor = 0.02;
   std::vector<pos_t> v;
   pos_t buv = (pos_2-pos_1).unit_vector();
   pos_t buv_90 = buv.rotate(90);
   double bond_length = pos_t::length(pos_2, pos_1);
   double length_scale = 0.1 * bond_length;
   pos_t short_edge_pt_1 = pos_2 + buv_90 * length_scale;
   pos_t short_edge_pt_2 = pos_2 - buv_90 * length_scale;
   // the line width means that the sharp angle an pos_1 here results
   // in a few pixels beyond the pos_1, so artificially shorten it a
   // tiny amount.
   //
   // Also, make it a quadralateral, with the sharp points very close,
   // this make the spike go away.
   //
   // 
   // pos_t sharp_point = pos_t::fraction_point(pos_1, pos_2, 0.11);
   // pos_t sharp_point = pos_t::fraction_point(pos_1, pos_2, 0.07);
   pos_t sharp_point = pos_t::fraction_point(pos_1, pos_2, 0.03);
   
   pos_t sharp_point_1 = sharp_point + buv_90 * sharp_point_sharp_factor;
   pos_t sharp_point_2 = sharp_point - buv_90 * sharp_point_sharp_factor;
   v.push_back(sharp_point_2);
   v.push_back(sharp_point_1);
   v.push_back(short_edge_pt_1);
   v.push_back(short_edge_pt_2);
   return v;
}

// return angle to X axis in degrees
double
lig_build::pos_t::axis_orientation() const {
   double angle = atan2(y,x)/DEG_TO_RAD;
   return angle;
}

bool
lig_build::pos_t::non_zero() const {
   if (fabs(x) > 0.00001)
      return true;
   if (fabs(y) > 0.00001)
      return true;
   return false;
} 

bool
lig_build::pos_t::operator==(const pos_t &pos) const {
   if (fabs(pos.x-x) < 0.00001) {
      if (fabs(pos.x-x) < 0.00001) {
	 return 1;
      } else {
	 return 0;
      }
   } else {
      return 0;
   } 
}


std::pair<lig_build::pos_t, lig_build::pos_t>
lig_build::bond_t::make_double_aromatic_short_stick(const pos_t &pos_1, const pos_t &pos_2) const {

   std::pair<lig_build::pos_t, lig_build::pos_t> p;
   lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
   lig_build::pos_t buv_90 = buv.rotate(90);

   // Which side of the pos_1 -> pos_2 vector shall we put this bond?
   // 
   // So create a T piece, and measure the distance to the centre
   // point, if we are on the inside then the distance to the centre
   // will be shorter.
   //
   double bond_length = lig_build::pos_t::length(pos_2, pos_1); // shortened possibly.
   double nice_dist = bond_length * 0.18;
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

   return std::pair<lig_build::pos_t, lig_build::pos_t>(cutened_inner_start_point,
							cutened_inner_end_point);
}

std::pair<std::pair<lig_build::pos_t, lig_build::pos_t>, std::pair<lig_build::pos_t, lig_build::pos_t> >
lig_build::bond_t::make_double_bond(const pos_t &pos_1, const pos_t &pos_2) const {

   lig_build::pos_t buv = (pos_2-pos_1).unit_vector();
   lig_build::pos_t buv_90 = buv.rotate(90);

   double small = lig_build::pos_t::length(pos_1, pos_2) * 0.1;
   lig_build::pos_t p1 = pos_1 + buv_90 * small;
   lig_build::pos_t p2 = pos_2 + buv_90 * small;
   lig_build::pos_t p3 = pos_1 - buv_90 * small;
   lig_build::pos_t p4 = pos_2 - buv_90 * small;

   std::pair<lig_build::pos_t, lig_build::pos_t> pair1(p1, p2);
   std::pair<lig_build::pos_t, lig_build::pos_t> pair2(p3, p4);
   return std::pair<std::pair<lig_build::pos_t, lig_build::pos_t>, std::pair<lig_build::pos_t, lig_build::pos_t> > (pair1, pair2);
} 



std::ostream&
lig_build::operator<<(std::ostream &s, lig_build::atom_t atom) {
   s << "[ atom at " << atom.atom_position.x << "," << atom.atom_position.y
     << " ele: " << atom.element << " charge: " << atom.charge << "]";
   return s;
}


std::ostream&
lig_build::operator<<(std::ostream &s, lig_build::bond_t bond) {
   s << "bond " << bond.atom_1 << " to " << bond.atom_2 << " type "
     << bond.bond_type;
   return s;
}

std::ostream&
lig_build::operator<<(std::ostream &s, const pos_t &p) {
   s << "[pos " << p.x << " " << p.y << "]";
   return s;
}

std::ostream&
lig_build::operator<<(std::ostream &s, atom_id_info_t a) {
   s << "atom_id_info: :" << a.atom_id << ": with " << a.offsets.size() << " offsets\n";
   for (unsigned int io=0; io<a.offsets.size(); io++) { 
      s << a[io];
   }
   return s;
}

std::ostream&
lig_build::operator<<(std::ostream &s, offset_text_t ot) {

      s << "    :" << ot.text << ": here/up/down: " << ot.text_pos_offset << " tweak: "
	<< ot.tweak << "\n";
      return s;
}

std::ostream&
lig_build::operator<<(std::ostream &s, atom_ring_centre_info_t wa) {
   s << wa.atom << " ring-centre: " << wa.has_ring_centre_flag;
   if (wa.has_ring_centre_flag) {
      s << " " << wa.ring_centre;
   }
   return s;
} 
