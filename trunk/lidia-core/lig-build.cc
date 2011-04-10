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
