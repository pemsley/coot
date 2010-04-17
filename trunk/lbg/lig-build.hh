/* lbg/lig-build.hh
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

#ifndef LIG_BUILD_HH
#define LIG_BUILD_HH

#include <math.h>  // for fabs, cos, sin
#include <string>
#include <vector>


#define DEG_TO_RAD .01745327 // defined elsewhere maybe.

enum { UNASSIGNED_INDEX = -1 };

namespace lig_build {

   class pos_t {
   public:
      double x;
      double y;
      pos_t(double x_in, double y_in) {
	 x = x_in;
	 y = y_in;
      }
      pos_t() {
	 x = -1;
	 y = -1;
      }
      pos_t unit_vector() const {
	 double l_sqrd = x*x + y*y;
	 double l = sqrt(l_sqrd);
	 return pos_t(x/l, y/l);
      }
      pos_t operator+(const pos_t &p) const {
	 return pos_t(x+p.x, y+p.y);
      }
      pos_t operator-(const pos_t &p) const {
	 return pos_t(x-p.x, y-p.y);
      }
      void operator+=(const pos_t &p) {
	 x += p.x;
	 y += p.y;
      }
      void operator-=(const pos_t &p) {
	 x -= p.x;
	 y -= p.y;
      }
      bool close_point(const pos_t &other) const {
	 bool status = 0;
	 double small = 3;
	 return near_point(other, small);
      }
      bool near_point(const pos_t &other, double small) const {
	 bool status = 0;
	 if (fabs(x-other.x) < small) { 
	    if (fabs(y-other.y) < small) {
	       status = 1;
	    }
	 }
	 return status; 
      }
      static pos_t fraction_point(const pos_t &first,
				  const pos_t &second,
				  double frac) {
	 double d_x = second.x - first.x;
	 double d_y = second.y - first.y;
	 return pos_t(first.x + frac * d_x,
				first.y + frac * d_y);
      }
      static pos_t mid_point(const pos_t &first,
				       const pos_t &second) {
	 return fraction_point(first, second, 0.5);
      }
      // angle in degrees
      pos_t rotate(double angle) const {
	 double theta = angle * DEG_TO_RAD;
	 double sin_theta = sin(theta);
	 double cos_theta = cos(theta);
	 double new_x = x * cos_theta - y * sin_theta;
	 double new_y = x * sin_theta + y * cos_theta;
	 return pos_t(new_x, new_y);
      }
      pos_t operator*(float sc) {
	 return pos_t(x*sc, y*sc);
      }
      double length() const {
	 return sqrt(x*x + y*y);
      } 
      static double length(const pos_t &p1, const pos_t &p2) {
	 double a = (p2.x - p1.x);
	 double b = (p2.y - p1.y);
	 double c = a*a + b*b;
	 if (c < 0.0) c = 0.0;
	 return sqrt(c);
      }
      // angle between this points (relative to the origin of course).
      static double angle(const pos_t &a, const pos_t &b) {
	 double a_dot_b = a.x * b.x + a.y * b.y;
	 double ab = a.length() * b.length();
	 // std::cout << a << " dot " << b << " is " << a_dot_b << std::endl;
	 double cos_theta = a_dot_b/ab;
	 // 	 std::cout << " |a| is " << a.length() << " |b| is " << b.length()
	 //                << std::endl;
	 // 	 std::cout << a << " cos_theta is " << a_dot_b << "/" << ab << " = "
	 // 		   << cos_theta << std::endl;
	 return acos(cos_theta)/DEG_TO_RAD;
      }
      static double cross(const pos_t &a, const pos_t &b) {
	 return (a.x*b.x + a.y*b.y);
      }
      // return angle to X axis in degrees
      double axis_orientation() const {
	 double angle = atan2(y,x)/DEG_TO_RAD;
	 return angle;
      }
      bool operator==(const pos_t &pos) const {
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
      friend std::ostream& operator<<(std::ostream &s, const pos_t &p);
   };
   std::ostream& operator<<(std::ostream &s, const pos_t &p);

   // -----------------------------------------------------------------
   //                   atom_t
   // -----------------------------------------------------------------

   // Note!
   //
   // element is the element, C, N, O etc
   //
   // atom_id is the text on the screen, e.g. NH, C+, OH etc.
   //
   // atom_name is the matching name from the PDB file (if any),
   // e.g. " C12", " OD2" etc.  This can be "" if there was no atom
   // name assigned.
   
   
   class atom_t {
      bool is_closed_; // because we can't delete atoms, we can only close
                       // them (deleting/erase()ing atoms would mess up
		       // the atom indexes in the bonds descriptions.
   public:
      pos_t atom_position;
      std::string element;
      std::string atom_id;
      int charge;
      bool aromatic;
      atom_t(pos_t pos_in, std::string ele_in, int charge_in) {
	 atom_position = pos_in;
	 atom_id = ele_in;
	 element = ele_in;
	 charge = charge_in;
	 is_closed_ = 0;
      }
      bool over_atom(const double &x_in, const double &y_in) const {
	 pos_t mouse(x_in, y_in);
	 double d = pos_t::length(mouse, atom_position);
	 if (d < 5)
	    return 1;
	 else
	    return 0;
      }
      // return the changed status (if input element is the same as
      // internal element, change status is 0).
      bool change_element(const std::string &ele_in) {
	 bool ch_status = 0;
	 if (element != ele_in) { 
	    element = ele_in;
	    ch_status = 1;
	 }
	 return ch_status;
      }
      void set_aromatic(bool flag) {
	 aromatic = flag;
      }
      bool atomatic_p() const {
	 return aromatic;
      }
      std::string get_atom_id() const { return atom_id; }
      // return atom-id-was-changed status
      bool set_atom_id(const std::string &atom_id_in ) {
	 bool changed_status = 0;
	 if (atom_id != atom_id_in) { 
	    changed_status = 1;
	    atom_id = atom_id_in;
	 }
	 return changed_status;
      }
      friend std::ostream& operator<<(std::ostream &s, atom_t);
      bool operator==(const atom_t &at) const {
	 bool status = 0;
	 if (at.atom_position == atom_position)
	    if (at.element == element)
	       if (at.atom_id == atom_id)
		  if (at.charge == charge)
		     if (at.aromatic == aromatic)
			status = 1;
	 return status;
      }
      bool is_closed() const { return is_closed_; }
      void close() {
	 // std::cout << " closing base class atom" << std::endl;
	 is_closed_ = 1;
      }
   };
   
   std::ostream& operator<<(std::ostream &s, atom_t);
   

   // -----------------------------------------------------------------
   //                   bond_t
   // -----------------------------------------------------------------
   class bond_t {
   public:
      enum bond_type_t { BOND_UNDEFINED=100, SINGLE_BOND=101, DOUBLE_BOND=102, TRIPLE_BOND=103,
			 IN_BOND=104, OUT_BOND=105 };
   protected: // atom_1 and atom_2 get swaped when turning an IN_BOND to an OUT_BOND;
      int atom_1;
      int atom_2;
   private:
      bond_type_t bond_type;
      pos_t centre_pos_; // the position of the polygen of
				  // which this bond is a part
      bool have_centre_pos_;  // was the bond from a polygon or just an external bond?
      bool is_closed_;
   public:
      bond_t() {
	 atom_1 = UNASSIGNED_INDEX;
	 atom_2 = UNASSIGNED_INDEX;
	 have_centre_pos_ = 0;
	 is_closed_ = 0;
	 bond_type = BOND_UNDEFINED;
      }
      bond_t(int first, int second, bond_type_t bt) {
	 atom_1 = first;
	 atom_2 = second;
	 bond_type = bt;
	 have_centre_pos_ = 0;
	 is_closed_ = 0;
      }
      bond_t(int first, int second, pos_t centre_pos_in, bond_type_t bt) {
	 atom_1 = first;
	 atom_2 = second;
	 bond_type = bt;
	 have_centre_pos_ = 1;
	 centre_pos_ = centre_pos_in;
	 is_closed_ = 0;
      }
      // mouse is hovering over bond?
      bool over_bond(double x, double y,
		     const atom_t &atom_1_at, const atom_t &atom_2_at) const;
      friend std::ostream& operator<<(std::ostream &s, bond_t);
      int get_atom_1_index() const {
	 return atom_1;
      }
      int get_atom_2_index() const {
	 return atom_2;
      }
      bond_type_t get_bond_type() const { return bond_type; }
      void set_bond_type(bond_type_t bt) { bond_type = bt; }
      bool have_centre_pos() const { return have_centre_pos_; }
      // post-hoc bond modification
      void set_centre_pos(const pos_t &pos_in) {
	 centre_pos_ = pos_in;
	 have_centre_pos_ = 1;
      } 
      pos_t centre_pos() const { return centre_pos_; }
      bool operator==(const bond_t &bond_in) const {
	 bool status = 0;
	 if (bond_in.get_bond_type() == bond_type)
	    status = 1;
	 return status;
      }
      void move_centre_pos(const pos_t &delta) {
	 centre_pos_ += delta;
      }
      void close() {
	 // std::cout << " closing base class bond" << std::endl;
	 is_closed_ = 1;
      }
      bool is_closed() const { return is_closed_; }
   };
   std::ostream& operator<<(std::ostream &s, bond_t);

   // -----------------------------------------------------------------
   //                   template molecule_t
   // -----------------------------------------------------------------

   template<class Ta, class Tb> class molecule_t {
   private:
      std::pair<bool,int> checked_add(const Ta &at) {
	 // check that there are no atoms that are close already positioned
	 int new_atom_index = UNASSIGNED_INDEX;
	 bool is_new = 0;
	 for (unsigned int iat=0; iat<atoms.size(); iat++) {
	    if (! atoms[iat].is_closed()) { 
	       if (atoms[iat].atom_position.near_point(at.atom_position, 2)) {
		  new_atom_index = iat; // an old atom
		  break;
	       } 
	    }
	 }

	 if (new_atom_index == UNASSIGNED_INDEX) { // not previously set
	    atoms.push_back(at);
	    new_atom_index = atoms.size() - 1;
	    is_new = 1;
	 }
	 return std::pair<bool, int> (is_new, new_atom_index);
      }
//       bool equivalent_bonds(<class Ta, class Tb> molecule_t) const {
// 	 bool status = 0;
// 	 return status;
//       } 
   public:
      std::vector<Ta> atoms;
      std::vector<Tb> bonds;
      
      // Return new atom index (int) and whether or not the atom was
      // added (1) or returned the index of an extant atom (0).
      // 
      virtual std::pair<bool, int> add_atom(const Ta &at) {
	 return checked_add(at);
      }
       int add_bond(const Tb &bond) {
	  bonds.push_back(bond);
	  return bonds.size() -1;
      }
      bool is_empty() const {
	 bool status = 1;
	 if (atoms.size())
	    status = 0;
	 if (bonds.size())
	    status = 0;
	 return status;
      }
      void clear() {
	 atoms.clear();
	 bonds.clear();
      }

      // "Put the template class method definitions in the .h file" - Thanks, Kevin!
      // 
      virtual std::vector<Tb> bonds_with_vertex(const pos_t &pos) const {
	 std::vector<Tb> rv;
	 for (unsigned int i=0; i<bonds.size(); i++) {
	    if (atoms[bonds[i].get_atom_1_index()].atom_position.near_point(pos, 1))
	       rv.push_back(bonds[i]);
	    if (atoms[bonds[i].get_atom_2_index()].atom_position.near_point(pos, 1))
	       rv.push_back(bonds[i]);
	 }
	 return rv;
      }

      // We dont want a copy of the bond, we want a reference to the
      // bond (so that it can be manipulated).
      // 
      std::vector<int> bond_indices_with_atom_index(int test_atom_index) const {
	 std::vector<int> rv;
	 if (!atoms[test_atom_index].is_closed()) { 
	    for (unsigned int i=0; i<bonds.size(); i++) {
	       if (bonds[i].get_atom_1_index() == test_atom_index) { 
		  if (! atoms[bonds[i].get_atom_2_index()].is_closed())
		     rv.push_back(i);
	       } else {
		  if (bonds[i].get_atom_2_index() == test_atom_index)
		     if (! atoms[bonds[i].get_atom_1_index()].is_closed())
			rv.push_back(i);
	       }
	    }
	 }
	 return rv;
      }

      // we have 2 indices, what is the bond that has these 2 atom indices?
      int get_bond_index(int atom_index_in_1, int atom_index_in_2) const {
	 int bi = UNASSIGNED_INDEX;
	 if (atom_index_in_1 != atom_index_in_2) { 
	    for (unsigned int i=0; i<bonds.size(); i++) { 
	       if ((bonds[i].get_atom_1_index() == atom_index_in_1) ||
		   (bonds[i].get_atom_2_index() == atom_index_in_1)) {
		  if ((bonds[i].get_atom_1_index() == atom_index_in_2) ||
		      (bonds[i].get_atom_2_index() == atom_index_in_2)) {
		     bi = i;
		     break;
		  }
	       }
	    }
	 }
	 return bi;
      } 

      // return the closest approach of test_pos to the atoms in the
      // molecule (except dont test against avoid_atom_index).
      //
      // return a status too saying that the distance was set.
      // 
      std::pair<bool, double>
      dist_to_other_atoms_except(const std::vector<int> &avoid_atom_index,
				 const pos_t &test_pos) const {
	 double dist_closest = 99999999999.9;
	 bool set_status = 0;
	 for (unsigned int iat=0; iat<atoms.size(); iat++) {
	    bool ok_atom = 1; 
	    for (unsigned int j=0; j<avoid_atom_index.size(); j++) {
	       if (iat == avoid_atom_index[j]) {
		  ok_atom = 0;
		  break;
	       }
	    }
	    if (ok_atom) { 
	       double d = pos_t::length(test_pos, atoms[iat].atom_position);
	       if (d < dist_closest) {
		  set_status = 1;
		  dist_closest = d;
	       }
	    }
	 }
	 return std::pair<bool, double> (set_status, dist_closest);
      }

      //
      std::string
      make_atom_id_by_using_bonds(const std::string &ele,
				    const std::vector<int> &bond_indices) const {
   
	 std::string atom_id = ele;
   
	 // Have we added an NH2, an NH or an N (for example)
	 //
	 int sum_neigb_bond_order = 0;
	 for (unsigned int ib=0; ib<bond_indices.size(); ib++) {
	    if (bonds[bond_indices[ib]].get_bond_type() == lig_build::bond_t::SINGLE_BOND)
	       sum_neigb_bond_order += 1;
	    if (bonds[bond_indices[ib]].get_bond_type() == lig_build::bond_t::IN_BOND)
	       sum_neigb_bond_order += 1;
	    if (bonds[bond_indices[ib]].get_bond_type() == lig_build::bond_t::OUT_BOND)
	       sum_neigb_bond_order += 1;
	    if (bonds[bond_indices[ib]].get_bond_type() == lig_build::bond_t::DOUBLE_BOND)
	       sum_neigb_bond_order += 2;
	    if (bonds[bond_indices[ib]].get_bond_type() == lig_build::bond_t::TRIPLE_BOND)
	       sum_neigb_bond_order += 3;
	 }

	 if (0) 
	    std::cout << "  in make_atom_id_by_using_bonds() "<< ele << " sum: "
		      << sum_neigb_bond_order << " " << std::endl;

	 if (ele == "N") {
	    if (sum_neigb_bond_order == 5)
	       atom_id = "N+2";
	    if (sum_neigb_bond_order == 4)
	       atom_id = "N+";
	    if (sum_neigb_bond_order == 3)
	       atom_id = "N";
	    if (sum_neigb_bond_order == 2)
	       atom_id = "NH";
	    if (sum_neigb_bond_order == 1)
	       atom_id = "NH2";
	    if (sum_neigb_bond_order == 0)
	       atom_id = "NH3";
	 }
							   
	 if (ele == "O") {
	    if (sum_neigb_bond_order == 3) { 
	       atom_id = "0+";
	    }
	    if (sum_neigb_bond_order == 2)
	       atom_id = "O";
	    if (sum_neigb_bond_order == 1)
	       atom_id = "OH";
	    if (sum_neigb_bond_order == 0)
	       atom_id = "OH2";
	 }

	 if (ele == "S") {
	    if (sum_neigb_bond_order == 1) {
	       atom_id = "SH";
	    }
	    if (sum_neigb_bond_order == 0) {
	       atom_id = "SH2";
	    }
	 }

	 if (ele == "C") {
	    if (sum_neigb_bond_order == 0) {
	       atom_id = "CH4";
	    }
	    if (sum_neigb_bond_order == 5) {
	       atom_id = "C+";
	    }
	 }

	 return atom_id;
      }

      // write out the atom and bond tables:
      // 
      void debug() const {
	 for (unsigned int i=0; i<atoms.size(); i++) {
	    std::cout << "Atom " << i << ": " << atoms[i].element << " "
		      << atoms[i].atom_id << " at "
		      << atoms[i].atom_position << std::endl;
	 }
	 for (unsigned int i=0; i<bonds.size(); i++) { 
	    std::cout << "Bond " << i << " atom indices: " << bonds[i].get_atom_1_index()
		      << " " << bonds[i].get_atom_2_index();
	    if (bonds[i].have_centre_pos())
	       std::cout << " centre_pos: " << bonds[i].centre_pos();
	    std::cout << " type " << bonds[i].get_bond_type() << std::endl;
	 }
      } 
      
   }; // end of molecule_t class
   

   class base_molecule_t : public molecule_t<atom_t, bond_t> {};

   class polygon_position_info_t {
   public:
      bool apply_internal_angle_offset_flag;
      pos_t pos;
      double angle_offset;
      bool can_stamp; // we are not overlaying an already existing ring (say).

      polygon_position_info_t() {
	 can_stamp = 0; // unset internals.
      } 
      polygon_position_info_t(double x_in, double y_in, double angle_in) {
	 pos.x = x_in;
	 pos.y = y_in;
	 angle_offset = angle_in;
	 apply_internal_angle_offset_flag = 1;
	 can_stamp = 1;
      }
      polygon_position_info_t(pos_t pos_in, double angle_in) {
	 pos = pos_in;
	 angle_offset = angle_in;
	 apply_internal_angle_offset_flag = 1;
	 can_stamp = 1;
      }
      polygon_position_info_t(pos_t pos_in, double angle_in, bool io) {
	 pos = pos_in;
	 angle_offset = angle_in;
	 apply_internal_angle_offset_flag = io;
	 can_stamp = 1;
      }
   };
}

#endif // LIG_BUILD_HH
