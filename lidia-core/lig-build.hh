/* lbg/lig-build.hh
 * 
 * Author: Paul Emsley
 * Copyright 2010 by The University of Oxford
 * Copyright 2015, 2016 by Medical Research Council
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
#include <iostream>
#include <algorithm>
#include <set>
#include <stdexcept>

#define MAX_SEARCH_DEPTH 9

#define DEG_TO_RAD .01745327 // defined elsewhere maybe.

enum { UNASSIGNED_INDEX = -1 };

namespace lig_build {

   enum { UNASSIGNED_BOND_INDEX = 2147483640 };

   // -----------------------------------------------------------------
   //                   pos_t
   // -----------------------------------------------------------------
   // 
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
	 double small_bit = 3;
	 return near_point(other, small_bit);
      }
      bool near_point(const pos_t &other, double small_bit) const {
	 bool status = 0;
	 if (fabs(x-other.x) < small_bit) { 
	    if (fabs(y-other.y) < small_bit) {
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
      static pos_t mid_point(const pos_t &first, const pos_t &second) {
	 return fraction_point(first, second, 0.5);
      }
      // angle in degrees
      pos_t rotate(double angle) const {
	 double theta_l = angle * DEG_TO_RAD;
	 double sin_theta = sin(theta_l);
	 double cos_theta = cos(theta_l);
	 double new_x = x * cos_theta - y * sin_theta;
	 double new_y = x * sin_theta + y * cos_theta;
	 return pos_t(new_x, new_y);
      }

      pos_t rotate_about(double x_cen, double y_cen, double angle) {
	 double theta_l = angle * DEG_TO_RAD;
	 double sin_theta = sin(theta_l);
	 double cos_theta = cos(theta_l);
	 double new_x = x_cen + (x - x_cen) * cos_theta - (y - y_cen) * sin_theta;
	 double new_y = y_cen + (x - x_cen) * sin_theta + (y - y_cen) * cos_theta;
	 return pos_t(new_x, new_y);
      }
      
      pos_t rotate_about(const pos_t &cen, double angle) {
	 double theta_l = angle * DEG_TO_RAD;
	 double sin_theta = sin(theta_l);
	 double cos_theta = cos(theta_l);
	 double new_x = cen.x + (x - cen.x) * cos_theta - (y - cen.y) * sin_theta;
	 double new_y = cen.y + (x - cen.x) * sin_theta + (y - cen.y) * cos_theta;
	 return pos_t(new_x, new_y);
      }
      
      pos_t operator*(float sc) const {
	 return pos_t(x*sc, y*sc);
      }
      double length() const {
	 return sqrt(x*x + y*y);
      } 
      double lengthsq() const {
	 return (x*x + y*y);
      }
      double theta() const {
	 // down the canvas is positve Y, so we need to correct for that
	 return atan2(-y,x);
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
      static double dot(const pos_t &a, const pos_t &b) {
	 return (a.x*b.x + a.y*b.y);
      }
      static std::vector<std::pair<pos_t, pos_t> > make_wedge_in_bond(const pos_t &pos_1,
								      const pos_t &pos_2);
      static std::vector<pos_t> make_wedge_out_bond(const pos_t &pos_1,
						    const pos_t &pos_2);

      // return angle to X axis in degrees
      double axis_orientation() const;
      bool non_zero() const;
      bool operator==(const pos_t &pos) const;
      friend std::ostream& operator<<(std::ostream &s, const pos_t &p);
   };
   std::ostream& operator<<(std::ostream &s, const pos_t &p);

   // -----------------------------------------------------------------
   //                   offset text
   // -----------------------------------------------------------------
   // 
   // 20110410 this is so that NH and OH can be typeset in various
   // ways, e.g. HO or with the N appearing above the H.  It depends
   // on the bonds to which this atom is attached.  So instead of
   // returning a simple string, as we did for
   // make_atom_id_by_using_bonds, we return a a container for vector
   // of offsets;
   //
   class offset_text_t {
   public:
      enum text_pos_offset_t { HERE=0, UP=-1, DOWN=1};
      offset_text_t(const std::string &text_in) {
	 text = text_in;
	 text_pos_offset = HERE;
	 tweak = pos_t(0,0);
	 subscript = false;
	 superscript = false;
      }
      offset_text_t(const std::string &text_in, text_pos_offset_t text_pos_offset_in) {
	 text = text_in;
	 text_pos_offset = text_pos_offset_in;
	 tweak = pos_t(0,0);
	 subscript = false;
	 superscript = false;
      }
      std::string text;
      text_pos_offset_t text_pos_offset;
      pos_t tweak;
      bool subscript;
      bool superscript;
      friend std::ostream& operator<<(std::ostream &s, offset_text_t a);
   };
   std::ostream& operator<<(std::ostream &s, offset_text_t a);

   class atom_id_info_t {
   public:

      // simple case
      atom_id_info_t(const std::string &atom_id_in) : atom_id(atom_id_in) {
	 offsets.push_back(offset_text_t(atom_id_in));
	 size_hint = 0;
      }

      // make a superscript for the formal charge
      atom_id_info_t(const std::string &atom_id_in, int formal_charge) : atom_id(atom_id_in){
	 offsets.push_back(offset_text_t(atom_id_in));
	 if (formal_charge != 0) { 
	    offset_text_t superscript("");
	    // if (formal_charge == -1) superscript=offset_text_t("-"); // try U+2012 or U+2013
	    if (formal_charge == -1) superscript=offset_text_t("−");
	    if (formal_charge == -2) superscript=offset_text_t("2-");
	    if (formal_charge == +1) superscript=offset_text_t("+");
	    if (formal_charge == +2) superscript=offset_text_t("2+");
	    superscript.superscript = true;
	    superscript.tweak = pos_t(8,0);
	    offsets.push_back(superscript);
	 } else {
	    // as above (simple constructor)
	    atom_id = atom_id_in;
	    offsets.push_back(offset_text_t(atom_id_in));
	 }
	 size_hint = 0;
      }

      // atom_id_info to be filled in by function that knows about
      // bonds by adding offsets.
      atom_id_info_t() {  size_hint = 0; }

      // convenience constructor (for NH2, OH2 typically)
      atom_id_info_t(const std::string &front,
		     const std::string &subscripted_text) {
	 offset_text_t ot1(front);
	 offset_text_t ot2(subscripted_text);
	 atom_id = front + subscripted_text;
	 ot1.tweak = pos_t(0,0);
	 ot2.tweak = pos_t(front.length()*10,0);
	 ot2.subscript = true;
	 offsets.push_back(ot1);
	 offsets.push_back(ot2);
	 size_hint = 0;
      } 

      std::vector<offset_text_t> offsets;
      int size_hint; // a hint so that atom names can be rendered smaller than normal
                     // -1 means smaller
      
      void set_atom_id(const std::string &atom_id_in) {
	 atom_id = atom_id_in;
      }
      
      const offset_text_t &operator[](const unsigned int &indx) const {
	 return offsets[indx];
      }
      void add(const offset_text_t &off) {
	 offsets.push_back(off);
      } 
      unsigned int n_offsets() const { return offsets.size(); } 
      
      std::string atom_id; // what we used to return, this is used for
			   // testing against what we had so we know
			   // if we need to change the id.
      friend std::ostream& operator<<(std::ostream &s, atom_id_info_t a);
   };
   std::ostream& operator<<(std::ostream &s, atom_id_info_t a);

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
      std::string atom_id;
      std::string element;
      std::string atom_name; // PDB atom names typically
      int charge;
      bool aromatic;
      atom_t(const pos_t &pos_in, const std::string &ele_in, int charge_in) :
         atom_position(pos_in), atom_id(ele_in), element(ele_in), charge(charge_in) {
	 is_closed_ = 0;
         aromatic = false;
      }
      ~atom_t() {}
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
      std::string get_atom_name() const {
	 return atom_name;
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
      // return the angle in degrees relative to the horizontal axis.
      // Recall that molecules are drawn upside down (inverted Y)
      // these are molecule coordinates though, not canvas coordinates.
      static double angle_of_bond(const atom_t &at_1, const atom_t &at_2) {
	 pos_t delta = at_2.atom_position - at_1.atom_position;
	 double theta = atan2(delta.y, delta.x);
	 return theta;
      }
   };
   
   std::ostream& operator<<(std::ostream &s, atom_t);
   
   // trivial container for a (copy of an) atom an its ring centre (if
   // it has one)
   class atom_ring_centre_info_t {
   public:
      atom_t atom;
      bool has_ring_centre_flag;
      pos_t ring_centre;
      atom_ring_centre_info_t(const atom_t &at) : atom(at) {
	 has_ring_centre_flag = 0;
      }
      void add_ring_centre(const pos_t &pos) {
	 ring_centre = pos;
	 has_ring_centre_flag = 1;
      }
   };
   std::ostream& operator<<(std::ostream &s, atom_ring_centre_info_t wa);

   // -----------------------------------------------------------------
   //                   bond_t
   // -----------------------------------------------------------------
   class bond_t {
   public:
      // IN_BOND and OUT_BOND are a type of single bond.
      // 
      enum bond_type_t { BOND_UNDEFINED=100, SINGLE_BOND=101, DOUBLE_BOND=102,
			 TRIPLE_BOND=103, AROMATIC_BOND=4, IN_BOND=104, OUT_BOND=105,
			 SINGLE_OR_DOUBLE=5, SINGLE_OR_AROMATIC=6,
			 DOUBLE_OR_AROMATIC=7,
			 DELOC_ONE_AND_A_HALF=8,
			 BOND_ANY=9 };
   protected: // atom_1 and atom_2 get swapped when turning an IN_BOND to an OUT_BOND;
      unsigned int atom_1;
      unsigned int atom_2;
   private:
      bond_type_t bond_type;
      pos_t centre_pos_; // the position of the polygen (i.e. ring of atoms) of
	   	         // which this bond is a part
      bool have_centre_pos_;  // was the bond from a polygon or just an external bond?
      int n_ring_atoms_; // ring double bond shortening factor depends on this
      bool is_closed_;
   public:
      bond_t() {
	 atom_1 = UNASSIGNED_INDEX;
	 atom_2 = UNASSIGNED_INDEX;
	 have_centre_pos_ = false;
	 is_closed_ = 0;
	 bond_type = BOND_UNDEFINED;
      }
      bond_t(unsigned int first, unsigned int second, bond_type_t bt) {
	 atom_1 = first;
	 atom_2 = second;
	 bond_type = bt;
	 have_centre_pos_ = false;
	 is_closed_ = 0;
      }
      bond_t(unsigned int first, unsigned int second, pos_t centre_pos_in, bond_type_t bt) {
	 atom_1 = first;
	 atom_2 = second;
	 bond_type = bt;
	 have_centre_pos_ = true;
	 centre_pos_ = centre_pos_in;
	 is_closed_ = 0;
      }
      ~bond_t() {}
      // mouse is hovering over bond?
      bool over_bond(double x, double y,
		     const atom_t &atom_1_at, const atom_t &atom_2_at) const;
      friend std::ostream& operator<<(std::ostream &s, bond_t);
      unsigned int get_atom_1_index() const {
	 return atom_1;
      }
      unsigned int get_atom_2_index() const {
	 return atom_2;
      }
      bond_type_t get_bond_type() const { return bond_type; }
      void set_bond_type(bond_type_t bt) { bond_type = bt; }
      bool have_centre_pos() const { return have_centre_pos_; }
      // post-hoc bond modification
      void set_centre_pos(const pos_t &pos_in, int n_ring_atoms_in) {
	 centre_pos_ = pos_in;
	 have_centre_pos_ = true;
	 n_ring_atoms_ = n_ring_atoms_in;
      }
      pos_t centre_pos() const { return centre_pos_; }
      void add_centre(const pos_t &centre_in, int n_ring_atoms_in) {
	 set_centre_pos(centre_in, n_ring_atoms_in);
      }
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
      bool matches_indices(const bond_t &test_bond) const {
	 if (get_atom_1_index() == test_bond.get_atom_1_index()) {
	    if (get_atom_2_index() == test_bond.get_atom_2_index()) {
	       return true;
	    } else {
	       return false;
	    }
	 } else {
	    return false;
	 } 
      }

      unsigned int get_other_index(const unsigned int &atom_index) const {
	 unsigned int idx = get_atom_1_index();
	 if (idx == atom_index)
	    idx = get_atom_2_index();
	 return idx;
      }

      // for double bonds.  We need more info to get this right.
      // Need connections to at_2
      //
      // pass const atom_t &at_1, const atom_t &at_2,
      bool draw_as_asymmetry_shortened(unsigned int n_neigh_1,
				       unsigned int n_neigh_2) const {
	 bool state = false;
	 if (bond_type == DOUBLE_BOND) {
	    if (have_centre_pos()) {
	       state = true;
	    } else {
	       // non-ring double bond
	       unsigned int n_C = 0;

	       // if (at_1.element == "C") n_C++;
	       // if (at_2.element == "C") n_C++;

	       bool do_symmetric_case = false;
	       if (n_C == 0) {
		  do_symmetric_case = true;
	       } else {
		  // more than 1 Carbon
		  if (n_C == 2) {
		     if (n_neigh_1 == 0)
			if (n_neigh_2 == 0)
			   do_symmetric_case = true;
		  }
		  if (n_C == 1) {
		     // if (at_1.element == "C")
		     if (n_neigh_1 > 1)
			do_symmetric_case = true;
		     // if (at_2.element == "C")
		     if (n_neigh_2 > 1)
			do_symmetric_case = true;
		  }
	       }
	       state = ! do_symmetric_case;
	    }
	 }
	 return state;
      }

      std::pair<pos_t, pos_t> coords_for_single_bond(const atom_t &at_1,
						     const atom_t &at_2,
						     bool at_1_is_singleton,
						     bool at_2_is_singleton) const {
	 pos_t pos_1_in = at_1.atom_position;
	 pos_t pos_2_in = at_2.atom_position;
	 pos_t pos_1 = pos_1_in;
	 pos_t pos_2 = pos_2_in;
	 bool shorten_first  = false;
	 bool shorten_second = false;

	 if (at_1.element != "C") {
	    shorten_first = true;
	 } else {
	    // Does this C have one bond?
	    if (at_1_is_singleton)
	       shorten_first = true;
	 }
	 if (at_2.element != "C") {
	    shorten_second = true;
	 } else {
	    // Does this C have one bond?
	    if (at_2_is_singleton)
	       shorten_second = true;
	 }
	 if (shorten_first) {
	    double shorten_fraction = 0.76; // standard
	    if (at_1.element == "Cl" || at_1.element == "Br" ||
		at_1.element == "As") {
	       double theta = atom_t::angle_of_bond(at_1, at_2); // radians
	       // std::cout << "got angle for " << at_1 << " " << at_2 << " " << a << std::endl;
	       // Actually, it would be better to have a more severe cut
	       // when the bond comes from the right (because of the l of the Cl)
	       // i.e. when theta is beteen -90 and 90 (-M_PI_4 -> M_PI_4)
	       double sc_1 = 0.5 * (1.0 + cos(theta - M_PI_4)); // 0 -> 1, max at the corners
	       shorten_fraction -= sc_1 * 0.25;
	    }
	    if (at_1.element == "N") {
	       double theta = atom_t::angle_of_bond(at_1, at_2);
	       double sc_1 = 0.5 * (1.0 + cos(theta - M_PI_4));
	       shorten_fraction -= sc_1 * 0.05; // need optimization
	    }

	    pos_1 = pos_t::fraction_point(pos_2_in, pos_1_in, shorten_fraction);

	 }
	 if (shorten_second) {
	    double shorten_fraction = 0.74; // standard
	    if (at_2.element == "Cl" || at_2.element == "Br" ||
		at_2.element == "As") {
	       double theta = atom_t::angle_of_bond(at_2, at_1); // radians
	       // std::cout << "got angle for " << at_1 << " " << at_2 << " " << a << std::endl;
	       double sc_1 = 0.5 * (1.0 + cos(theta - M_PI_4)); // 0 -> 1, max at the corners
	       shorten_fraction -= sc_1 * 0.25;
	    }
	    if (at_2.element == "N") {
	       double theta = atom_t::angle_of_bond(at_2, at_1);
	       double sc_1 = 0.5 * (1.0 + cos(theta - M_PI_4));
	       shorten_fraction -= sc_1 * 0.05;
	    }
	    pos_2 = pos_t::fraction_point(pos_1_in, pos_2_in, shorten_fraction);
	 }

	 return std::pair<pos_t, pos_t> (pos_1, pos_2);
      }

      std::vector<pos_t>
      coords_for_sheared_or_darted_wedge_bond(const pos_t &pos_1,
					      const pos_t &pos_2,
					      const std::vector<std::pair<atom_t, bond_t> > &other_connections_to_second_atom) const {

	 std::vector<pos_t> pts;

	 if (other_connections_to_second_atom.size() > 0) {

	    if (other_connections_to_second_atom.size() == 1) {

	       const pos_t  &third_atom_pos = other_connections_to_second_atom[0].first.atom_position;
	       const bond_t &third_bond     = other_connections_to_second_atom[0].second;

	       if (third_bond.bond_type == TRIPLE_BOND) {
		  // make an ordinary wedge (but shortened)
		  pos_t pos_2_cut = pos_t::fraction_point(pos_2, pos_1, 0.16);
		  pts = pos_t::make_wedge_out_bond(pos_1, pos_2_cut);

	       } else {

		  pos_t b = pos_2-pos_1;
		  pos_t buv = b.unit_vector();
		  pos_t buv_90 = buv.rotate(90);
		  pos_t sharp_point = pos_t::fraction_point(pos_1, pos_2, 0.04);

		  pos_t sharp_point_1 = sharp_point + buv_90 * 0.03; // was 0.03
		  pos_t sharp_point_2 = sharp_point - buv_90 * 0.03; // ditto

		  pos_t bfrom3rd = pos_2 - third_atom_pos;
		  pos_t bond_from_3rd_atom_extension   = pos_2 + bfrom3rd*0.08;  // was 0.1
		  pos_t bond_from_3rd_atom_contraction = pos_2 - bfrom3rd*0.16;  // was 0.18

		  if (third_bond.get_bond_type() == bond_t::DOUBLE_BOND) {
		     // do we need to make this shorter? (Looks OK in cairo-molecule)
		     // bond_from_3rd_atom_extension   -= b * 0.05;
		     // bond_from_3rd_atom_contraction -= b * 0.05;
		  }

		  if (false) {
		     std::cout << " pos_1           " << pos_1 << std::endl;
		     std::cout << " pos_2           " << pos_2 << std::endl;
		     std::cout << " buv             " << buv << std::endl;
		     std::cout << " 3rd atom pos    " << third_atom_pos << std::endl;
		     std::cout << " bfrom3rd        " << bfrom3rd << std::endl;
		     std::cout << " sheared points: " << sharp_point_2 << std::endl;
		     std::cout << "                 " << sharp_point_1 << std::endl;
		     std::cout << "                 " << bond_from_3rd_atom_extension   << std::endl;
		     std::cout << "                 " << bond_from_3rd_atom_contraction << std::endl;
		  }

		  // is bfrom3rd_2 in the same direction as buv_90?
		  // If not, then we need to swap around sharp_point_1 and sharp_point_2.
		  //
		  double dp = pos_t::dot(bfrom3rd, buv_90);
		  if (dp < 0)
		     std::swap(sharp_point_1, sharp_point_2);

		  pts.push_back(sharp_point_2);
		  pts.push_back(sharp_point_1);
		  pts.push_back(bond_from_3rd_atom_extension);
		  pts.push_back(bond_from_3rd_atom_contraction);
	       }

	    } else {

	       // make a dart (there are 2 third atoms)
	 
	       const lig_build::pos_t  &third_atom_1_pos = other_connections_to_second_atom[0].first.atom_position;
	       const lig_build::pos_t  &third_atom_2_pos = other_connections_to_second_atom[1].first.atom_position;

	       const bond_t &third_bond_1     = other_connections_to_second_atom[0].second;
	       const bond_t &third_bond_2     = other_connections_to_second_atom[1].second;
	       lig_build::pos_t b = pos_2-pos_1;
	       lig_build::pos_t buv = b.unit_vector();
	       lig_build::pos_t buv_90 = buv.rotate(90);
	       lig_build::pos_t sharp_point = lig_build::pos_t::fraction_point(pos_1, pos_2, 0.04);
	       lig_build::pos_t sharp_point_1 = sharp_point + buv_90 * 0.03;
	       lig_build::pos_t sharp_point_2 = sharp_point - buv_90 * 0.03;

	       lig_build::pos_t bfrom3rd_1 = pos_2 - third_atom_1_pos;
	       lig_build::pos_t bfrom3rd_2 = pos_2 - third_atom_2_pos;
	       lig_build::pos_t bond_from_3rd_atom_1_contraction = pos_2 - bfrom3rd_1*0.15;
	       lig_build::pos_t bond_from_3rd_atom_2_contraction = pos_2 - bfrom3rd_2*0.15;

	       // is bfrom3rd_2 in the same direction as buv_90?
	       // If not, then we need to swap around sharp_point_1 and sharp_point_2.
	       //
	       double dp = pos_t::dot(bfrom3rd_2, buv_90);
	       if (dp < 0)
		  std::swap(sharp_point_1, sharp_point_2);

	       // for double bonds, we don't want to have the dart point along the bond,
	       // we want it to be on the line of the inner bond
	       //
	       // but we don't want to do this for double bonds in rings
	       // (because those double bonds are not draw with symmetric offset)
	       // in fact, there are non-ring double bonds drawn with offset
	       // Hmmm... How do we know if the third bonds are drawn with offset?
	       //
	       // atom_t third_bond_1_atom_1 = atoms[third_bond_1.get_atom_1_index()];
	       // atom_t third_bond_1_atom_2 = atoms[third_bond_1.get_atom_2_index()];
	       // atom_t third_bond_2_atom_1 = atoms[third_bond_2.get_atom_1_index()];
	       // atom_t third_bond_2_atom_2 = atoms[third_bond_2.get_atom_2_index()];
	       // I don't know what is first or second - baah - needs some
	       // clear thinking.
	       // std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > oc =
	       // make_other_connections_to_second_atom_info(bond_index);
               unsigned int n_1 = other_connections_to_second_atom.size(); // yes TIWIM
	       unsigned int n_2 = 1;
	       if (third_bond_1.get_bond_type() == bond_t::DOUBLE_BOND)
		  if (! third_bond_1.draw_as_asymmetry_shortened(n_1, n_2))
		     // change by half of the separation between double bonds
		     bond_from_3rd_atom_1_contraction -= b * 0.1;
	       if (third_bond_2.get_bond_type() == bond_t::DOUBLE_BOND)
		  if (! third_bond_2.draw_as_asymmetry_shortened(n_1, n_2))
		     // change by half of the separation between double bonds
		     bond_from_3rd_atom_2_contraction -= b * 0.1;

	       pts.push_back(sharp_point_2);
	       pts.push_back(sharp_point_1);
	       pts.push_back(bond_from_3rd_atom_1_contraction);
	       pts.push_back(pos_2);
	       pts.push_back(bond_from_3rd_atom_2_contraction);

	    }
	 }
	 return pts;
      }
      // (the long/normal stick for this is just pos_1_in to pos_2_in)
      std::pair<pos_t, pos_t>
      make_double_aromatic_short_stick(const pos_t &pos_1_in, const pos_t &pos_2_in,
				       bool shorten_first,
				       bool shorten_second) const;

      // symmetric offset version
      std::pair<std::pair<pos_t, pos_t>, std::pair<pos_t, pos_t> >
			       make_double_bond(const pos_t &pos_1, const pos_t &pos_2,
						bool shorten_first, bool shorten_second) const;

      std::pair<std::pair<pos_t, pos_t>, std::pair<pos_t, pos_t> >
			       make_double_bond(const pos_t &pos_1, const pos_t &pos_2,
						bool shorten_first, bool shorten_second,
						const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
						const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom) const;

   };
   std::ostream& operator<<(std::ostream &s, bond_t);


   // -----------------------------------------------------------------
   //                   template molecule_t
   // -----------------------------------------------------------------

   // from sub-class back down here.
   // 20111229

   template<class Ta, class Tb> class molecule_t {
   private:
      std::pair<bool,int> checked_add(const Ta &at) {
	 // check that there are no atoms that are close already positioned
	 int new_atom_index = UNASSIGNED_INDEX;
	 bool is_new = 0;
	 for (unsigned int iat=0; iat<atoms.size(); iat++) {
	    if (! atoms[iat].is_closed()) {

	       // FIXME (rethink)
	       // when we import 3D coords, 2 is not a good distance.  This whole
	       // "checked_add" thing should be reconsidered.  Why did I need it?
	       // 
	       // if (atoms[iat].atom_position.near_point(at.atom_position, 2)) {
	       if (atoms[iat].atom_position.near_point(at.atom_position, 0.01)) {
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
      bool member(const unsigned int &ind, const std::vector<unsigned int> &no_pass_atoms) const {
	 bool found = false;
	 for (unsigned int i=0; i<no_pass_atoms.size(); i++) { 
	    if (no_pass_atoms[i] == ind) {
	       found = true;
	       break;
	    }
	 }
	 return found;
      } 

      // this helps find the fist ring containing start_atom_index that goes through atom_index_other
      // 
      std::pair<bool, std::vector<unsigned int> >
      find_bonded_atoms_with_no_pass(unsigned int start_atom_index,
				     unsigned int atom_index_other,
				     unsigned int this_atom_index,
				     const std::vector<unsigned int> &no_pass_atoms,
				     unsigned int depth) const {
	 
	 std::vector<unsigned int> atoms_bonded_to_atom_index_start;
	 std::vector<unsigned int> local_no_pass_atoms = no_pass_atoms;
	 if (depth == 0) {
	    std::vector<unsigned int> empty;
	    return std::pair<bool, std::vector<unsigned int> > (0, empty);
	 } else {

	    // get a list of all the atoms that are bonded to this atom not
	    // in the no_pass list
      
	    for (unsigned int i=0; i<bonds.size(); i++) { 
	       if (bonds[i].get_atom_1_index() == this_atom_index) {
		  unsigned int idx = bonds[i].get_atom_2_index();
		  if (idx == start_atom_index)
		     if (depth < (MAX_SEARCH_DEPTH-1)) {
			if (member(atom_index_other, local_no_pass_atoms)) {
			   // if this_atom_index was not already in
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
			   return std::pair<bool, std::vector<unsigned int> > (1, local_no_pass_atoms);
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
		  unsigned int idx = bonds[i].get_atom_1_index();
		  if (idx == start_atom_index)
		     if (depth < (MAX_SEARCH_DEPTH-1)) {
			if (member(atom_index_other, local_no_pass_atoms)) { 
			   local_no_pass_atoms.push_back(this_atom_index);
			   // if this_atom_index was not alread in
			   // local_no_pass_atoms, then add it.
			   bool ifound = false; 
			   for (unsigned int ilnp=0; ilnp<local_no_pass_atoms.size(); ilnp++) {
			      if (local_no_pass_atoms[ilnp] == this_atom_index) {
				 ifound = true;
				 break;
			      }
			   }
			   if (! ifound)
			      local_no_pass_atoms.push_back(this_atom_index);
			   // debug_pass_atoms(start_atom_index, this_atom_index, depth, local_no_pass_atoms);
			   return std::pair<bool, std::vector<unsigned int> > (1, local_no_pass_atoms);
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
	       std::pair<bool, std::vector<unsigned int> > r =
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
	 std::vector<unsigned int> empty;
	 return std::pair<bool, std::vector<unsigned int> > (0, empty);
      }

      std::vector<std::set<unsigned int> >
      find_rings_including_atom_internal(unsigned int start_atom_index,
					 unsigned int atom_index_other,
					 unsigned int this_atom_index,
					 const std::set<unsigned int> &no_pass_atoms,
					 unsigned int depth) const {

	 std::vector<std::set<unsigned int> > v;

	 std::set<unsigned int> atoms_bonded_to_this_atom;
	 std::set<unsigned int> local_no_pass_atoms = no_pass_atoms;

	 if (depth == 0) {
	    std::vector<std::set<unsigned int> > empty;
	    return empty;
	 } else {

	    // construct a list of all the atoms that are bonded to this atom not
	    // in the no_pass list

	    for (unsigned int i=0; i<bonds.size(); i++) {
	       if (bonds[i].get_atom_1_index() == this_atom_index) {
		  unsigned int idx = bonds[i].get_atom_2_index();
		  if (idx == start_atom_index) {
		     if (depth < (MAX_SEARCH_DEPTH-1)) {

			if (local_no_pass_atoms.find(atom_index_other) != local_no_pass_atoms.end()) {
			   local_no_pass_atoms.insert(this_atom_index);
			   v.push_back(local_no_pass_atoms);
			}
		     }
		  }
		  if (local_no_pass_atoms.find(idx) == local_no_pass_atoms.end()) {
		     atoms_bonded_to_this_atom.insert(idx);
		     local_no_pass_atoms.insert(this_atom_index);
		  }
	       }

	       if (bonds[i].get_atom_2_index() == this_atom_index) {
		  unsigned int idx = bonds[i].get_atom_1_index();
		  if (idx == start_atom_index) {
		     if (depth < (MAX_SEARCH_DEPTH-1)) {

			if (local_no_pass_atoms.find(atom_index_other) != local_no_pass_atoms.end()) {
			   local_no_pass_atoms.insert(this_atom_index);
			   std::vector<std::set<unsigned int> > f;
			   v.push_back(local_no_pass_atoms);
			}
		     }
		  }
		  if (local_no_pass_atoms.find(idx) == local_no_pass_atoms.end()) {
		     atoms_bonded_to_this_atom.insert(idx);
		     local_no_pass_atoms.insert(this_atom_index);
		  }
	       }
	    } // end the bond loop

	    std::set<unsigned int>::const_iterator it;
	    for (it  = atoms_bonded_to_this_atom.begin();
		 it != atoms_bonded_to_this_atom.end();
		 ++it) {
	       std::vector<std::set<unsigned int> > r =
		  find_rings_including_atom_internal(start_atom_index,
						     atom_index_other,
						     *it, local_no_pass_atoms, depth-1);
	       // merge r onto v
	       for (unsigned int jj=0; jj<r.size(); jj++)
		  v.push_back(r[jj]);
	    }
	 }
	 return v;
      }
      

      std::vector<std::set<unsigned int> >
      find_rings_including_atom_simple_internal(unsigned int start_atom_index,
						unsigned int this_atom_index,
						const std::set<unsigned int> &no_pass_atoms,
						unsigned int depth) const {

	 std::vector<std::set<unsigned int> > v;

	 std::set<unsigned int> atoms_bonded_to_this_atom;
	 std::set<unsigned int> local_no_pass_atoms = no_pass_atoms;

	 if (depth == 0) {
	    std::vector<std::set<unsigned int> > empty;
	    return empty;
	 } else {
      
	    // construct a list of all the atoms that are bonded to this atom not
	    // in the no_pass list
      
	    for (unsigned int i=0; i<bonds.size(); i++) { 
	       if (bonds[i].get_atom_1_index() == this_atom_index) {
		  unsigned int idx = bonds[i].get_atom_2_index();
		  if (idx == start_atom_index) {
		     if (depth < (MAX_SEARCH_DEPTH-1)) {

			local_no_pass_atoms.insert(this_atom_index);
			v.push_back(local_no_pass_atoms);
		     }
		  }
		  if (local_no_pass_atoms.find(idx) == local_no_pass_atoms.end()) {
		     atoms_bonded_to_this_atom.insert(idx);
		     local_no_pass_atoms.insert(this_atom_index);
		  }
	       }

	       if (bonds[i].get_atom_2_index() == this_atom_index) {
		  unsigned int idx = bonds[i].get_atom_1_index();
		  if (idx == start_atom_index) {
		     if (depth < (MAX_SEARCH_DEPTH-1)) {

			local_no_pass_atoms.insert(this_atom_index);
			v.push_back(local_no_pass_atoms);
		     }
		  }
		  if (local_no_pass_atoms.find(idx) == local_no_pass_atoms.end()) {
		     atoms_bonded_to_this_atom.insert(idx);
		     local_no_pass_atoms.insert(this_atom_index);
		  }
	       }
	    } // end the bond loop

	    std::set<unsigned int>::const_iterator it;
	    for (it  = atoms_bonded_to_this_atom.begin();
		 it != atoms_bonded_to_this_atom.end();
		 it++) {
	       std::vector<std::set<unsigned int> > r =
		  find_rings_including_atom_simple_internal(start_atom_index, *it,
							    local_no_pass_atoms, depth-1);
	       // merge r onto v
	       for (unsigned int jj=0; jj<r.size(); jj++)
		  v.push_back(r[jj]);
	    }
	 }
	 return v;
      }

      // can throw an exception (no bonds)
      //
      // not const because it now caches the return value;
      //
      bool have_cached_bond_ring_centres_flag;
      //
      // ring centres should incorporate the number of atoms contributing to the
      // ring, i.e. std::vector<std::pair<unsigned int, pos_t> >
      // so that the ring centre info can be used to calculate the correct
      // position of the inner bond of a ring double bond
      //
      std::vector<pos_t> cached_bond_ring_centres;

   public:
      molecule_t() {
	 have_cached_bond_ring_centres_flag = false;
      }
      std::vector<Ta> atoms;
      std::vector<Tb> bonds;

      virtual ~molecule_t() = 0;

      // Return new atom index (int) and whether or not the atom was
      // added (1) or returned the index of an extant atom (0).
      // 
      virtual std::pair<bool, int> add_atom(const Ta &at) {
	 return checked_add(at);
      }

      // Add a bond. But only add a bond if there is not already a
      // bond between the given atom indices in the bond list.
      //
      // In that case, return -1, otherwise return the index of the
      // newly added bond.
      // 
      int add_bond(const Tb &bond) {
	 
	 bool matches_indices = 0;
	 for (unsigned int i=0; i<bonds.size(); i++) { 
	    if (bonds[i].matches_indices(bond)) {
	       matches_indices = 1;
	       break;
	    }
	 }
	 if (! matches_indices) { 
	    bonds.push_back(bond);
	    return bonds.size() -1;
	 } else {
	    return -1;
	 } 
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

      std::pair<bool, std::vector<unsigned int> >
      found_self_through_bonds(unsigned int atom_index_start,
			       unsigned int atom_index_other) const {
	 std::vector<unsigned int> empty_no_pass_atoms;
	 empty_no_pass_atoms.push_back(atom_index_start);
	 std::pair<bool, std::vector<unsigned int> > r =
	    find_bonded_atoms_with_no_pass(atom_index_start, atom_index_start,
					   atom_index_other, empty_no_pass_atoms,
					   MAX_SEARCH_DEPTH);
	 return r;
      }

      // use like find_bonded_atoms_with_no_pass
      std::vector<std::set<unsigned int> > rings_including_atom(unsigned int atom_index_start,
								unsigned int atom_index_other) const {

         std::set<unsigned int> empty_no_pass_atoms;
         empty_no_pass_atoms.insert(atom_index_start);
         // use like find_bonded_atoms_with_no_pass
         std::vector<std::set<unsigned int> > v =
	    find_rings_including_atom_internal(atom_index_start, atom_index_start,
             	                               atom_index_other, empty_no_pass_atoms, MAX_SEARCH_DEPTH);
         return v;
      }

      // as above - but we don't have a bond through which this ring must pass.
      std::vector<std::set<unsigned int> > rings_including_atom(unsigned int atom_index_start) const {

         std::set<unsigned int> empty_no_pass_atoms;
         empty_no_pass_atoms.insert(atom_index_start);
         // use like find_bonded_atoms_with_no_pass
         std::vector<std::set<unsigned int> > v =
	    find_rings_including_atom_simple_internal(atom_index_start, atom_index_start,
						      empty_no_pass_atoms, MAX_SEARCH_DEPTH);
         return v;
      }

      bool in_ring_p(unsigned int atom_index) const {
	 bool in_ring = false;
	 std::vector<std::set<unsigned int> > v = rings_including_atom(atom_index);
	 for (unsigned int i=0; i<v.size(); i++) {
	    if (v[i].find(atom_index) != v[i].end()) {
	       in_ring = true;
	       break;
	    }
	 }
	 return in_ring;
      }

      // return a vector of bond indices
      //
      std::vector<unsigned int> ring_atoms_to_ring_bonds(const std::set<unsigned int> &atom_indices) const {
	 std::vector<unsigned int> ring_bonds;
	 for (std::size_t iat=0; iat<atom_indices.size(); iat++) {
	    // what are the bonds that contain iat?
	    std::vector<unsigned int> bds = bonds_having_atom_with_atom_index(iat);
	    for (std::size_t ibd=0; ibd<bds.size(); ibd++) {
	       int idx_other = bonds[bds[ibd]].get_other_index(iat);
	       if (std::find(atom_indices.begin(),
			     atom_indices.end(), idx_other) != atom_indices.end()) {
		  if (std::find(ring_bonds.begin(), ring_bonds.end(), ibd) == ring_bonds.end())
		     ring_bonds.push_back(ibd);
	       }
	    }
	 }
	 return ring_bonds;
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

      // can throw a runtime_error exception
      Ta &operator[](unsigned int idx) {
	 if (idx >= atoms.size()) {
	    std::string s = "molecule_t atom index error: ";
	    throw std::runtime_error(s);
	 } else {
	    return atoms[idx];
	 }
      }

      // We dont want a copy of the bond, we want a reference to the
      // bond (so that it can be manipulated).
      // 
      std::vector<unsigned int> bond_indices_with_atom_index(unsigned int test_atom_index) const {
	 std::vector<unsigned int> rv;
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
      unsigned int get_bond_index(unsigned int atom_index_in_1, unsigned int atom_index_in_2) const {
	 unsigned int bi = UNASSIGNED_BOND_INDEX;
	 if (atom_index_in_1 != atom_index_in_2) {
	    for (unsigned int i=0; i<bonds.size(); i++) {
	       if (! bonds[i].is_closed()) {
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
	 }
	 return bi;
      }

      void debug_pass_atoms(unsigned int atom_index_start, unsigned int this_atom_index,
			    unsigned int depth,
			    const std::vector<unsigned int> &local_no_pass_atoms) const {

	 std::cout << "    found atom index " << atom_index_start << " from this atom: "
		   << this_atom_index
		   << ", at depth " << depth << " no-pass-atoms: (";
	 for (unsigned int i=0; i<local_no_pass_atoms.size(); i++) { 
	    std::cout << local_no_pass_atoms[i] << " ";
	 }
	 std::cout << ")" << std::endl;
      }

      void debug_pass_atoms(unsigned int atom_index_start, unsigned int this_atom_index,
			    unsigned int depth,
			    const std::set<unsigned int> &local_no_pass_atoms) const {

	 std::cout << "    debug_pass_atoms(): found atom index "
		   << atom_index_start << " from this atom: "
		   << this_atom_index
		   << ", at depth " << depth << " no-pass-atoms: (";
	 std::set<unsigned int>::const_iterator it;
	 for (it=local_no_pass_atoms.begin(); it!=local_no_pass_atoms.end(); it++)
	    std::cout << *it << " ";
	 std::cout << ")" << std::endl;
      }
      

      unsigned int get_number_of_atoms_including_hydrogens() const {
	 return atoms.size();
      }

      void assign_ring_centres(bool force=false) {
	 bool debug = false;

	 // first cache the ring info for the atoms.  The ring info should ideally/in future contain
	 // information about the number (and placement?) of the double bonds in the ring (see comments
	 // for favourite_ring_id().

	 for (unsigned int ib=0; ib<bonds.size(); ib++) {
	    if (! bonds[ib].have_centre_pos() || force) {
	       unsigned int atom_index     = bonds[ib].get_atom_1_index();
	       unsigned int atom_idx_other = bonds[ib].get_atom_2_index();
	       if (debug)
		  std::cout << "=============== checking bond " << ib
			    << " for rings for atom index "
			    << atom_index << " ===============" << std::endl;

	       // all the rings for atom_index that pass throught atom_idx_other (because
	       // for a double bond connected to a fused ring atom we want the ring that contains
	       // this bond (not the other ring).
	       // 
	       std::vector<std::set<unsigned int> > rings = rings_including_atom(atom_index, atom_idx_other);

	       if (debug) {
		  std::cout << "   constructor of widgeted_bond_t atom " << atom_index
			    << " other bond index (not tested) " << bonds[ib].get_atom_2_index()
			    << ", found n_rings: " << rings.size() << std::endl;
		  for (unsigned int ir=0; ir<rings.size(); ir++) {
		     std::cout << "   ring " << ir << " :: ";
		     std::set<unsigned int>::const_iterator it;
		     for (it=rings[ir].begin(); it != rings[ir].end(); it ++)
			std::cout << " " << *it;
		     std::cout << std::endl;
		  }

	       }
	       if (rings.size() > 0) {
		  lig_build::pos_t centre_pos_sum;
		  std::string centre_pos_atoms_string;
		  unsigned int fav_ring_id = favourite_ring_id(rings);
		  std::set<unsigned int>::const_iterator it;
		  for (it=rings[fav_ring_id].begin(); it != rings[fav_ring_id].end(); it ++) {
		     centre_pos_sum += atoms[*it].atom_position;
		  }
		  lig_build::pos_t centre_pos =
		     centre_pos_sum * (1.0/double(rings[fav_ring_id].size()));
		  bonds[ib].add_centre(centre_pos, rings[fav_ring_id].size());
		  if (debug) {
		     std::cout << "   adding centre at " << centre_pos
			       << " generated from (";
		     std::cout << "   (";
		     for (it=rings[fav_ring_id].begin(); it != rings[fav_ring_id].end(); it ++)
			std::cout << " " << *it;
		     std::cout << ")" 
			       << " for bond " << ib
			       << " connecting " << bonds[ib].get_atom_1_index() << " to "
			       << bonds[ib].get_atom_2_index() << std::endl;
		  }
	       }
	    }
	 }
      }

      // Say double bond is in a 5-membered and a 6-membered ring.
      // On which ring should a double bond be placed?
      // The 6-membered.  Or failing that it has a 6-membered ring, the biggest ring.
      //
      // IUPAC Guidelines for chemical structure diagrams says that it should be offset
      // to the inside of the ring that contains the most double bonds.
      // In 6-membered rings or alternating single and double bonds, it is especially prefered
      // if the bond is on the inside of this ring (in preference to the ring with the
      // largest number of double bonds). e.g. tetralin.
      // 
      // calling function should ensure that there should be at least one ring passed in rings.
      // 
      unsigned int favourite_ring_id(const std::vector<std::set<unsigned int> > &rings) const {

	 unsigned int id = 0;
	 unsigned int n_best = 0;

	 for (unsigned int i=0; i<rings.size(); i++) {
	    std::vector<unsigned int> ring_bonds = ring_atoms_to_ring_bonds(rings[i]);
	    unsigned int n_bonds_double = 0;
	    for (std::size_t ib=0; ib<ring_bonds.size(); ib++) {
	       if (false)
		  std::cout << "in favourite_ring_id() iring " << i << " ib " << ib
			    << " ring_bonds[ " << ib << "] = " << ring_bonds[ib]
			    << " " << bonds[ring_bonds[ib]] << std::endl;
	       if (bonds[ring_bonds[ib]].get_bond_type() == bond_t::DOUBLE_BOND ||
		   bonds[ring_bonds[ib]].get_bond_type() == bond_t::AROMATIC_BOND ||
		   bonds[ring_bonds[ib]].get_bond_type() == bond_t::DOUBLE_OR_AROMATIC ||
		   bonds[ring_bonds[ib]].get_bond_type() == bond_t::SINGLE_OR_AROMATIC) {
		  n_bonds_double++;
	       }
	    }
	    if (n_bonds_double > n_best) {
	       n_best = n_bonds_double;
	       id = i;
	    }
	 }
	 return id;
      }

      // 
      std::vector<pos_t> get_ring_centres() {

	 if (have_cached_bond_ring_centres_flag) {
	    return cached_bond_ring_centres;
	 } else { 
	    std::vector<pos_t> v;
	    for (unsigned int ib=0; ib<bonds.size(); ib++) {
	       if (bonds[ib].have_centre_pos()) {
		  pos_t rc = bonds[ib].centre_pos();
		  bool found = 0;
		  for (unsigned int i=0; i<v.size(); i++) {
		     if (v[i].near_point(rc, 7)) { // "7"? is that gocanvas coords? FIXME
			found = 1;
			break;
		     }
		  }
		  if (! found)
		     v.push_back(rc);
	       }
	    }
// 	    std::cout << "DEBUG:: get_ring_centres() returns\n";
// 	    for (unsigned int iv=0; iv<v.size(); iv++) {
// 	       std::cout << "   "  << iv << " " << v[iv] << "\n";
// 	    }
	    cached_bond_ring_centres = v;
	    have_cached_bond_ring_centres_flag = true;
	    return v;
	 }
      }

      // can throw an exception (no atoms)
      // 
      pos_t
      get_ring_centre(const std::vector<std::string> &ring_atom_names) const {

	 pos_t positions_sum(0,0);
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
	 pos_t centre = positions_sum * (1.0/double(n_found));
	 return centre;
      }

      // can throw an exception (no rings with this atom)
      //
      pos_t
      get_ring_centre(const atom_ring_centre_info_t &rc_atom) const {

	 pos_t position(0,0);
	 bool found = 0;

	 for (unsigned int ibond=0; ibond<bonds.size(); ibond++) { 
	    if ((atoms[bonds[ibond].get_atom_1_index()] == rc_atom.atom) ||
		(atoms[bonds[ibond].get_atom_2_index()] == rc_atom.atom)) {
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
	    mess += rc_atom.atom.get_atom_name();
	    mess += " found to be in a ring";
	    throw(std::runtime_error(mess));
	 }
	 return position;
      }


      // can throw an exception (no atoms) return pseudo points top-left
      // (small small) bottom-right (high high).
      //
      // return coords for top_left, bottom_right
      // 
      std::pair<pos_t, pos_t>
      ligand_extents() const {

	 pos_t top_left;
	 pos_t bottom_right;

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
	    top_left     = pos_t(mol_min_x, mol_min_y);
	    bottom_right = pos_t(mol_max_x, mol_max_y);
	 } else {
	    std::string mess = "WARNING:: no atoms in ligand_extents()";
	    throw std::runtime_error(mess);
	 }
	 return std::pair<pos_t, pos_t> (top_left, bottom_right);
      }


      int n_open_bonds() const {

	 int n_bonds = 0;
	 for (unsigned int i=0; i<bonds.size(); i++) { 
	    if (! bonds[i].is_closed())
	       n_bonds++;
	 }
	 return n_bonds;
      }
      
      bool is_close_to_non_last_atom(const pos_t &test_pos) const {

	 bool close = false;
	 int n_atoms_for_test = atoms.size() - 1; 
	 for (int iat=0; iat<n_atoms_for_test; iat++) {
	    if (! atoms[iat].is_closed()) {
	       if (atoms[iat].atom_position.near_point(test_pos, 2.1)) {
		  close = true;
		  break;
	       }
	    }
	 }
	 return close;
      }

      void delete_hydrogens() {
	 for (unsigned int iat=0; iat<atoms.size(); iat++) { 
	    if (atoms[iat].element == "H") {
	       std::vector<unsigned int> bds = bonds_having_atom_with_atom_index(iat);
	       atoms[iat].close();
	       for (unsigned int i=0; i<bds.size(); i++)
		  bonds[bds[i]].close();
	    }
	 }
      }

      std::vector<unsigned int> bonds_having_atom_with_atom_index(unsigned int test_atom_index) const {

	 std::vector<unsigned int> v;
	 std::vector<unsigned int> vb =  bond_indices_with_atom_index(test_atom_index);
	 
	 for (unsigned int iv=0; iv<vb.size(); iv++) {
	    if (! bonds[vb[iv]].is_closed())
	       v.push_back(vb[iv]);
	 }
	 
	 return v;
      }

      // return a vector of the atom indices of unconnected atoms
      //
      std::vector<unsigned int> get_unconnected_atoms() const {

	 std::vector<unsigned int> v;
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
      pos_t get_ligand_centre() const {

	 pos_t centre(0,0);

	 if (atoms.size() == 0) {
	    std::string message("No atoms in ligand");
	    throw std::runtime_error(message);
	 } else {
	    pos_t centre_sum(0,0);
	    for (unsigned int iat=0; iat<atoms.size(); iat++) {
	       centre_sum += atoms[iat].atom_position;
	    }
	    if (atoms.size() > 0)
	       centre = centre_sum * (1.0/double(atoms.size()));
	 }
	 return centre;
      }

      double median_bond_length() const {
	 double bl = -1;
	 std::vector<double> bls;
	 bls.reserve(bonds.size());
	 for (unsigned int ib=0; ib<bonds.size(); ib++) {
	    int idx_1 = bonds[ib].get_atom_1_index();
	    int idx_2 = bonds[ib].get_atom_2_index();
	    if ((idx_1 != UNASSIGNED_INDEX) && (idx_2 != UNASSIGNED_INDEX)) {
	       const pos_t &pos_1 = atoms[idx_1].atom_position;
	       const pos_t &pos_2 = atoms[idx_2].atom_position;
	       pos_t delta = pos_2 - pos_1;
	       bls.push_back(delta.length());
	    }
	 }
	 if (bls.size() > 0) {
	    std::sort(bls.begin(), bls.end());
	    unsigned int idx = bls.size() / 2;
	    bl = bls[idx];
	 }
	 return bl;
      }

      bool
      operator==(const molecule_t &mol_other) const {

	 bool status = 0;

	 // need to check that bonds are the same (atom indexing can be
	 // different) and also stray atoms need to be checked after bonds.
	 //
	 int n_bond_hits = 0;

	 if (mol_other.bonds.size() != bonds.size())
	    return status;
	 for (unsigned int ib=0; ib<bonds.size(); ib++) {
	    atom_t atom_i_1 = atoms[bonds[ib].get_atom_1_index()];
	    atom_t atom_i_2 = atoms[bonds[ib].get_atom_2_index()];
	    int n_hits = 0;
	    for (unsigned int jb=0; jb<mol_other.bonds.size(); jb++) {
	       atom_t atom_j_1 = mol_other.atoms[bonds[ib].get_atom_1_index()];
	       atom_t atom_j_2 = mol_other.atoms[bonds[ib].get_atom_2_index()];
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
	       std::vector<unsigned int> i_stray_atoms = stray_atoms();
	       std::vector<unsigned int> j_stray_atoms = mol_other.stray_atoms();
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


      unsigned int n_stray_atoms() const {  // unbonded atoms
	 return stray_atoms().size();
      }

      unsigned int num_atoms() const {
	 int n = 0;
	 for (unsigned int i=0; i<atoms.size(); i++) {
	    if (!atoms[i].is_closed())
	       n++;
	 }
	 return n;
      }
      
      unsigned int num_bonds_for_atom(int atom_index) {
	 unsigned n = 0;
	 for (unsigned int ib=0; ib<bonds.size(); ib++) {
	    if (! bonds[ib].is_closed()) {
	       int iat_1 = bonds[ib].get_atom_1_index();
	       int iat_2 = bonds[ib].get_atom_2_index();
	       if (! atoms[iat_1].is_closed())
		  if (iat_1 == atom_index)
		     n++;
	       if (! atoms[iat_2].is_closed())
		  if (iat_2 == atom_index)
		     n++;
	    }
	 }
	 return n;
      }
      
      std::vector<unsigned int> stray_atoms() const {
	 std::vector<unsigned int> strays;
	 std::vector<bool> found(atoms.size(), false);
	 for (unsigned int ib=0; ib<bonds.size(); ib++) { 
	    int iat_1 = bonds[ib].get_atom_1_index();
	    int iat_2 = bonds[ib].get_atom_2_index();
	    if (! atoms[iat_1].is_closed())
	       found[iat_1] = true;
	    if (! atoms[iat_2].is_closed())
	       found[iat_2] = true;
	 }
	 for (unsigned int i=0; i<atoms.size(); i++)
	    if (! atoms[i].is_closed())
	       if (! found[i])
		  strays.push_back(i);
	 return strays;
      }

      void translate(const pos_t &delta) {

	 for (unsigned int iat=0; iat<atoms.size(); iat++) {
	    atoms[iat].atom_position += delta;
	 }
	 for (unsigned int ib=0; ib<bonds.size(); ib++) {
	    if (bonds[ib].have_centre_pos()) {
	       bonds[ib].move_centre_pos(delta);
	    }
	 }
      }


      // return the closest approach of test_pos to the atoms in the
      // molecule (except dont test against avoid_atom_index).
      //
      // return a status too saying that the distance was set.
      // 
      std::pair<bool, double>
      dist_to_other_atoms_except(const std::vector<unsigned int> &avoid_atom_index,
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

      // the sum of the neighbour vectors from (relative to) the central atom
      //
      pos_t get_sum_delta_neighbours(unsigned int atom_index,
				     const std::vector<unsigned int> &bond_indices) const {
	 pos_t sum_delta(0,0);
	 for (unsigned int ibond=0; ibond<bond_indices.size(); ibond++) {
	    int idx_other = bonds[bond_indices[ibond]].get_other_index(atom_index);
	    pos_t delta = atoms[idx_other].atom_position - atoms[atom_index].atom_position;
	    sum_delta += delta;
	 }
	 return sum_delta;
      }

      pos_t get_sum_delta_neighbours(const atom_t &atom) const {
	 pos_t sum_delta(0,0);
	 for (unsigned int ibond=0; ibond<bonds.size(); ibond++) {
	    int idx_1 = bonds[ibond].get_atom_1_index();
	    int idx_2 = bonds[ibond].get_atom_2_index();
	    if (atoms[idx_1] == atom) {
	       // std::cout << "found atom idx1 " << atom << " in atoms: " << idx_1 << std::endl;
	       sum_delta += atoms[idx_2].atom_position - atom.atom_position;
	    }
	    if (atoms[idx_2] == atom) {
	       // std::cout << "found atom idx2 " << atom << " in atoms: " << idx_1 << std::endl;
	       sum_delta += atoms[idx_1].atom_position - atom.atom_position;
	    }
	 }
	 return sum_delta;
      }

      // to draw double bonds without centre correctly (and below)
      std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> >
      make_other_connections_to_first_atom_info(unsigned int bond_index) const {
	 std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > v;
	 int atom_1_ref_idx = bonds[bond_index].get_atom_1_index();
	 int atom_2_ref_idx = bonds[bond_index].get_atom_2_index();

	 for (unsigned int ibond=0; ibond<bonds.size(); ibond++) {
	    if (ibond != bond_index) {
	       int at_1_idx = bonds[ibond].get_atom_1_index();
	       int at_2_idx = bonds[ibond].get_atom_2_index();
	       if (at_1_idx == atom_1_ref_idx) {
		  if (at_2_idx != atom_2_ref_idx) {
		     std::pair<lig_build::atom_t, lig_build::bond_t> p(atoms[at_2_idx], bonds[ibond]);
		     v.push_back(p);
		  }
	       }
	       if (at_2_idx == atom_1_ref_idx) {
		  if (at_1_idx != atom_2_ref_idx) {
		     std::pair<lig_build::atom_t, lig_build::bond_t> p(atoms[at_1_idx], bonds[ibond]);
		     v.push_back(p);
		  }
	       }
	    }
	 }
	 return v;
      }

      // to draw wedge bonds correctly
      std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> >
      make_other_connections_to_second_atom_info(unsigned int bond_index) const {
	 std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > v;
	 int atom_chiral_idx = bonds[bond_index].get_atom_1_index();
	 int atom_other_idx  = bonds[bond_index].get_atom_2_index();

	 for (unsigned int ibond=0; ibond<bonds.size(); ibond++) {
	    if (ibond != bond_index) {
	       int at_1_idx = bonds[ibond].get_atom_1_index();
	       int at_2_idx = bonds[ibond].get_atom_2_index();
	       if (at_1_idx == atom_other_idx) {
		  if (at_2_idx != atom_chiral_idx) { // should always be
		     std::pair<lig_build::atom_t, lig_build::bond_t> p(atoms[at_2_idx], bonds[ibond]);
		     v.push_back(p);
		  }
	       }
	       if (at_2_idx == atom_other_idx) {
		  if (at_1_idx != atom_chiral_idx) {
		     std::pair<lig_build::atom_t, lig_build::bond_t> p(atoms[at_1_idx], bonds[ibond]);
		     v.push_back(p);
		  }
	       }
	    }
	 }

	 // std::cout << "from make_other_connections_to_second_atom_info() returning v of size "
	 // << v.size() << std::endl;

	 return v;
      }
      
      std::pair<bool, bool> shorten_flags(unsigned int bond_idx) const {
	 bool shorten_first  = false;
	 bool shorten_second = false;
	 if (bond_idx < bonds.size()) {
	    unsigned int idx_1 = bonds[bond_idx].get_atom_1_index();
	    unsigned int idx_2 = bonds[bond_idx].get_atom_2_index();
	    if (atoms[idx_1].element != "C") {
	       shorten_first = true;
	    } else {
	       // Does this C have one bond?
	       std::vector<unsigned int> v = bonds_having_atom_with_atom_index(idx_1);
	       if (v.size() == 1)
		  shorten_first = true;
	    }
	    if (atoms[idx_2].element != "C") {
	       shorten_second = true;
	    } else {
	       std::vector<unsigned int> v = bonds_having_atom_with_atom_index(idx_2);
	       if (v.size() == 1)
		  shorten_second = true;
	    }
	 }
	 return std::pair<bool, bool> (shorten_first, shorten_second);
      }

      // using bonds and charge
      // 
      atom_id_info_t
      make_atom_id_by_using_bonds(unsigned int atom_index,
				  const std::string &ele,
				  const std::vector<unsigned int> &bond_indices,
				  bool simple_gl_render) const {

	 // The offset tweaks depend on the way the text is positioned on the canvas.
	 // presumably simple OpenGL is anchored bottom left.
	 //
	 // So, for example, whereas before "OH" needed an x-tweak
	 // with simple_gl_render, it does not.
	 
	 std::string atom_id = ele;
	 int charge = atoms[atom_index].charge;
   
	 // Have we added an NH2, an NH or an N (for example)
	 //
	 int sum_neigb_bond_order = 0;
	 // This is a hack to deal with aromatic bonds (non-integers - doh!)
	 int n_aromatic = 0;

	 // 20160701: This was commented out.  I wonder why.  I have
	 //           now uncommented it.
	 //           It was commented out because when we have charge the
	 //           sum_neigb_bond_order gets confusing, e.g. N+H4 for
	 //           sum_neigb_bond_order == 1 and
	 //           N+H3 for sum_neigb_bond_order == 2.
	 // sum_neigb_bond_order += charge;
	 
	 for (unsigned int ib=0; ib<bond_indices.size(); ib++) {
	    if (bonds[bond_indices[ib]].get_bond_type() == bond_t::SINGLE_BOND)
	       sum_neigb_bond_order += 1;
	    if (bonds[bond_indices[ib]].get_bond_type() == bond_t::IN_BOND)
	       sum_neigb_bond_order += 1;
	    if (bonds[bond_indices[ib]].get_bond_type() == bond_t::OUT_BOND)
	       sum_neigb_bond_order += 1;
	    if (bonds[bond_indices[ib]].get_bond_type() == bond_t::DOUBLE_BOND)
	       sum_neigb_bond_order += 2;
	    if (bonds[bond_indices[ib]].get_bond_type() == bond_t::TRIPLE_BOND)
	       sum_neigb_bond_order += 3;
	    if (bonds[bond_indices[ib]].get_bond_type() == bond_t::AROMATIC_BOND)
	       n_aromatic++;
	 }

	 if (n_aromatic == 3) // unusual
	    sum_neigb_bond_order += 3;
	 if (n_aromatic == 2) // pyridine
	    sum_neigb_bond_order += 3;
	 if (n_aromatic == 1) // weird
	    sum_neigb_bond_order += 2;

	 if (ele == "N") {

	    if (charge == -2) { // weird
	       atom_id = "N-2";
	    }

	    if (charge == -1) { // weird
	       atom_id = "N-";
	    }
	    
	    if (charge == 0) {
	       if (sum_neigb_bond_order == 0)
		  atom_id = "NH3";
	       if (sum_neigb_bond_order == 1)
		  atom_id = "NH2";
	       if (sum_neigb_bond_order == 2)
		  atom_id = "NH";
	       if (sum_neigb_bond_order == 3)
		  atom_id = "N";
	       if (sum_neigb_bond_order == 4)
		  atom_id = "N";
	    }
	    if (charge == +1) { // very typical
	       if (sum_neigb_bond_order == 0)
		  atom_id = "NH4+"; // ammonium
	       if (sum_neigb_bond_order == 1)
		  atom_id = "NH3+"; // lys
	       if (sum_neigb_bond_order == 2)
		  atom_id = "NH2+";
	       if (sum_neigb_bond_order == 3)
		  atom_id = "NH+";
	       if (sum_neigb_bond_order == 4)
		  atom_id = "N+";
	       if (sum_neigb_bond_order == 5)
		  atom_id = "N+"; // weird
	    }
	    if (charge == +2) { // all weird
	       if (sum_neigb_bond_order == 0)
		  atom_id = "N+2H4";
	       if (sum_neigb_bond_order == 1)
		  atom_id = "N+2H3";
	       if (sum_neigb_bond_order == 2)
		  atom_id = "N2H2+";
	       if (sum_neigb_bond_order == 3)
		  atom_id = "NH+2";
	       if (sum_neigb_bond_order == 4)
		  atom_id = "N+2";
	    }
	 }
							   
	 if (ele == "O") {

	    if (charge == 2) {
	       if (sum_neigb_bond_order == 1)
		  atom_id = "O+2H3";
	       if (sum_neigb_bond_order == 2)
		  atom_id = "O+2H2";
	    }

	    if (charge == 1) {
	       if (sum_neigb_bond_order == 1)
		  atom_id = "O+H2";
	       if (sum_neigb_bond_order == 2)
		  atom_id = "O+H";
	    }
	    
	    if (charge == 0) {
	       if (sum_neigb_bond_order == 0)
		  atom_id = "OH2";
	       if (sum_neigb_bond_order == 1)
		  atom_id = "OH";
	       if (sum_neigb_bond_order == 2)
		  atom_id = "O";
	       if (sum_neigb_bond_order == 3)
		  atom_id = "O"; // doesn't make sense
	    }

	    if (charge == -1) {
	       atom_id = "O-";
	    }

	    if (charge == -2) {
	       atom_id = "O-2";
	    }
	 }

	 if (ele == "S") {

	    // 20160701 need charge-specific atom_ids here?
	    
	    if (sum_neigb_bond_order == 0)
	       atom_id = "SH2";
	    if (sum_neigb_bond_order == 1)
	       atom_id = "SH";
	    if (sum_neigb_bond_order == 2)
	       atom_id = "S";
	    if (sum_neigb_bond_order == 3)
	       atom_id = "S";
	 }

	 if (ele == "C") {
	    if (sum_neigb_bond_order == 0) {
	       return atom_id_info_t("CH", "4");
	    }
	    if (charge == 1) {
	       if (sum_neigb_bond_order == 5) {
		  atom_id = "C="; // Hmm... carbon can go either way
	       }
	    }

	    if (charge == 0) {

	       if (bond_indices.size() > 1) {
		  atom_id = "C";
	       } else {
		  std::string h_count = "3";
		  if (sum_neigb_bond_order == 2)
		     h_count = "2";

		  pos_t sum_delta = get_sum_delta_neighbours(atom_index, bond_indices);
		  double theta = (180.0/M_PI) * atan2(sum_delta.y, sum_delta.x);
		  if (! (theta > -95.0 && theta < 95.0)) {
		     // if (sum_delta.x < 3.2) { // prefer CH3 to H3C when (nearly) vertical. (old)
		     atom_id_info_t id("CH");
		     offset_text_t ot(h_count);
		     ot.tweak = pos_t(18, 0);
		     ot.subscript = true;
		     id.add(ot);
		     return id;
		  } else {
		     atom_id_info_t id;
		     offset_text_t otH("H");
		     offset_text_t ot2(h_count);
		     offset_text_t otN("C");
		     ot2.subscript = true;
		     otH.tweak = pos_t(-16, 0);
		     ot2.tweak = pos_t(-8, 0);
		     otN.tweak = pos_t(0, 0);
		     id.add(otH);
		     id.add(ot2);
		     id.add(otN);
		     return id;
		  }
	       }
	    }
	 }

	 // ----------------- OK, atom_id is set, how do we type set that?

	 if (atom_id != "NH" && atom_id != "OH" && atom_id != "SH") {

	    if (atom_id == "NH2") {
	       // Might be typeset as "NH2" or "H2N", let's investigate

	       if (bond_indices.size() != 1) {
		  // strange case - I don't know what else to do..
		  return atom_id_info_t("NH", "2");
		  
	       } else {

		  // If the bond is on the left of the text, we're good.
		  // If it comes from the right, then we need to write
		  // H2N, with the 2 subscripted.
		  // 
		  pos_t sum_delta = get_sum_delta_neighbours(atom_index, bond_indices);
		  double theta = (180.0/M_PI) * atan2(sum_delta.y, sum_delta.x);
		  // std::cout << "::: here for NH2 with sum_delta " << sum_delta
		  // << " theta " << theta << std::endl;
		  // old (dependent on molecule scaling)
		  // if (sum_delta.x < 3.2) { // prefer NH2 to H2N when (nearly) vertical.

		  if (theta > -95.0 && theta < 95.0) {
		     // more tricky case then...
		     atom_id_info_t id;
		     offset_text_t otH("H");
		     offset_text_t ot2("2");
		     offset_text_t otN("N");
		     ot2.subscript = true;
		     otH.tweak = pos_t(-16, 0);
		     ot2.tweak = pos_t(-7, 0);
		     otN.tweak = pos_t(0, 0);
		     id.set_atom_id("NH2");
		     id.add(otH);
		     id.add(ot2);
		     id.add(otN);
		     return id;
		  } else {
		     return atom_id_info_t("NH", "2");
		  }
	       }

	    } else {

	       // not NH2...
	       // 
	       // still in the [not NH, OH, SH] block.

	       if (atom_id == "N+") {
		  // nitro

		  atom_id_info_t atom_id_info("N", 1);
		  // offset_text_t n("N");
		  // atom_id_info.add(n);
		  // offset_text_t pl("+");
		  // pl.superscript = true;
		  // pl.tweak = pos_t(6, 0);
		  // atom_id_info.add(pl);
		  return atom_id_info;

	       } else {

		  if (atom_id == "NH+") {

		     pos_t sum_delta = get_sum_delta_neighbours(atom_index, bond_indices);
		     if (bond_indices.size() == 1) {
			if (sum_delta.x < 0) {
			   atom_id_info_t id("N", 1);
			   offset_text_t h("H");
			   h.tweak = pos_t(14, 0);
			   id.add(h);
			   return id;
			} else {
			   atom_id_info_t id;
			   offset_text_t h("H");
			   offset_text_t n("N");
			   offset_text_t plus("+");
			   h.tweak = pos_t(-10, 0);
			   plus.tweak = pos_t(8, 0);
			   plus.superscript = true;
			   id.add(h);
			   id.add(n);
			   id.add(plus);
			   return id;
			}
		     } else {
			if (bond_indices.size() == 2) {

			   if (fabs(sum_delta.y) > fabs(sum_delta.x)) {

			      if (sum_delta.y > 0) {

				 //    \   /
				 //      N+
				 //      H
				 //

				 atom_id_info_t id("N", 1);
				 offset_text_t h("H", offset_text_t::DOWN);
				 id.add(h);
				 return id;
				 
			      } else {

				 //      H
				 //      N+
				 //    /   \           .[not-multi-line comment]
				 //
				 atom_id_info_t id("N", 1);
				 offset_text_t h("H", offset_text_t::UP);
				 id.add(h);
				 return id;

			      }
			   } else {
			      // normal sideways

			      if (sum_delta.x < 0) {
				 atom_id_info_t id("N", 1);
				 offset_text_t h("H");
				 h.tweak = pos_t(14, 0);
				 id.add(h);
				 return id;
			      } else {
				 atom_id_info_t id;
				 offset_text_t h("H");
				 offset_text_t n("N");
				 offset_text_t plus("+");
				 h.tweak = pos_t(-10, 0);
				 plus.tweak = pos_t(8, 0);
				 plus.superscript = true;
				 id.add(h);
				 id.add(n);
				 id.add(plus);
				 return id;
			      }
			   }
			   
			} else {
			   // number of bond_indices is not 2.
			   // what is the madness?
			   return atom_id_info_t("NH+");
			}
		     } 
		  } else {
		     
		     if (atom_id == "NH2+") {

			pos_t sum_delta = get_sum_delta_neighbours(atom_index, bond_indices);

			if (bond_indices.size() == 2) {

			   if (fabs(sum_delta.y) > fabs(sum_delta.x)) {

			      if (sum_delta.y > 0) {

				 //      H2
				 //      N+
				 //    /   \           .[not multi-line comment]
				 //
				 atom_id_info_t id("N", 1);
				 offset_text_t h("H", offset_text_t::DOWN); // weird
				 offset_text_t t("2", offset_text_t::DOWN);
				 t.subscript = true;
				 t.tweak = pos_t(10,-4);
				 id.add(h);
				 id.add(t);
				 return id;

			      } else {
		  
				 //    \   /
				 //      N+
				 //      H2
				 //
				 atom_id_info_t id("N", 1);
				 offset_text_t h("H", offset_text_t::UP);
				 offset_text_t t("2", offset_text_t::UP);
				 t.subscript = true;
				 t.tweak = pos_t(9,-2);
				 id.add(h);
				 id.add(t);
				 return id;

			      }
			   } else {
			      // normal sideways of N+H2

			      if (sum_delta.x < 0) {
				 atom_id_info_t id("N"); // No charge set here
				 offset_text_t plus("+"); // the plus goes after the H
				 offset_text_t h("H");
				 offset_text_t t("2");
				 h.tweak = pos_t(10, 0);
				 t.tweak = pos_t(18, 0);
				 plus.tweak = pos_t(18, 0);
				 t.subscript = true;
				 plus.superscript = true;
				 id.add(h);
				 id.add(t);
				 id.add(plus);
				 return id;
			      } else {
				 atom_id_info_t id;
				 offset_text_t h("H");
				 offset_text_t t("2");
				 offset_text_t n("N");
				 offset_text_t plus("+");
				 h.tweak = pos_t(-16, 0);
				 t.tweak = pos_t(-7, 0);
				 t.subscript = true;
				 plus.tweak = pos_t(8, 0);
				 plus.superscript = true;
				 id.add(h);
				 id.add(n);
				 id.add(t);
				 id.add(plus);
				 return id;
			      }
			   }
			   
			} else {

			   // sum bond_indices is 1 (well, not 2)

			   if (sum_delta.x < 0) {
			      // bond from the left
			      atom_id_info_t id("N", 1);
			      offset_text_t h("H");
			      h.tweak = pos_t(15, 0);
			      offset_text_t two("2");
			      two.subscript = true;
			      two.tweak = pos_t(24, 0);
			      id.add(h);
			      id.add(two);
			      return id;
			   } else {
			      // bond from the right
			      atom_id_info_t id("N", 1);
			      offset_text_t h("H");
			      offset_text_t two("2");
			      h.tweak = pos_t(-14, 0);
			      two.tweak = pos_t(-6, 0);
			      two.subscript = true;
			      id.add(h);
			      id.add(two);
			      return id;
			   }
			}

		     } else {

			if (atom_id == "O-") {
			   return atom_id_info_t("O", -1);

			} else {
			   if (atom_id == "O-2") {
			      return atom_id_info_t("O", -2);

			   } else {

			      if (atom_id == "OH2") {
			
				 return atom_id_info_t("OH", "2");

			      } else {
				 
				 if (atom_id == "N-2") {
				    return atom_id_info_t("N", -2);

				 } else {
				    if (atom_id == "N-") {
				       return atom_id_info_t("N", -1);
			   
				    } else {
				       atom_id_info_t simple(atom_id);
				       return simple;
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 } else {

	    // NH or OH or SH
	    
	    pos_t sum_delta = get_sum_delta_neighbours(atom_index, bond_indices);
	    atom_id_info_t atom_id_info;
	    
	    if (bond_indices.size() == 1) {
	       // OH, HO, NH, HN
	       if (sum_delta.x > 0) {

		  atom_id_info.set_atom_id(atom_id);
		  std::string txt = "HO";
		  if (ele == "N")
		     txt = "HN";
		  if (ele == "S")
		     txt = "HS";
		  offset_text_t ot(txt);
		  ot.tweak = pos_t(-10, 0);
		  atom_id_info.add(ot);
	       } else {
		  // simple
		  atom_id_info = atom_id_info_t(atom_id);

		  // We don't need to do this with a simple GL renderer
		  if (! simple_gl_render) 
		     // add a tweak to that, push it to the right a tiny
		     // bit (we want the O (or N) at the end of the bond)
		     atom_id_info.offsets.back().tweak += pos_t(0,0);
	       } 
	    }

	    if (bond_indices.size() == 2) {

	       // Does this happen for S?

	       // Add a tweak factor (1.3) to prefer horizontal orientation.
	       if (fabs(sum_delta.y) > fabs(sum_delta.x) * 1.3) { 
		  
		  //    \   /
		  //      N
		  //      H
		  //
		  if (sum_delta.y > 0) {
		     offset_text_t n(ele);
		     offset_text_t h("H", offset_text_t::DOWN);
		     atom_id_info.add(n);
		     atom_id_info.add(h);
		  } else { 
		  
		  
		  //      H
		  //      N
		  //    /   \  .
		  //
		     offset_text_t n(ele);
		     offset_text_t h("H", offset_text_t::UP);
		     atom_id_info.add(n);
		     atom_id_info.add(h);
		  } 
	       } else {

		  // these are not simple, they need a position tweak.
		  if (sum_delta.x > 0.0) {
		     // H pokes to the left
		     atom_id_info = atom_id_info_t();
		     offset_text_t n(std::string("H") + ele);
		     n.tweak = pos_t(-10,0);
		     atom_id_info.add(n);
		  } else {
		     atom_id_info = atom_id_info_t();
		     offset_text_t n(ele + "H");
		     n.tweak = pos_t(0,0);
		     atom_id_info.add(n);
		  } 
	       } 
	    }

// 	    std::cout << "----------- make_atom_id_by_using_bonds() returning-------------- "
// 		      << atom_id_info.offsets.size() << " offsets "
// 		      << std::endl;
// 	    for (unsigned int ioff=0; ioff<atom_id_info.offsets.size(); ioff++) { 
// 	       std::cout << "   :"
// 			 << atom_id_info.offsets[ioff].text << ": here/up/down: "
// 			 << atom_id_info.offsets[ioff].text_pos_offset << " tweak: "
// 			 << atom_id_info.offsets[ioff].tweak << std::endl;
// 	    }
// 	    std::cout << std::endl;
	    
	    return atom_id_info;
	 }
      } // end of make_atom_id_by_using_bonds()

      bool delete_bond_between(unsigned int idx_1, unsigned int idx_2) {
	 bool status = false;

	 for (unsigned int ibond=0; ibond<bonds.size(); ibond++) { 
	    unsigned int index_1 = bonds[ibond].get_atom_1_index();
	    unsigned int index_2 = bonds[ibond].get_atom_2_index();
	    if (idx_1 == index_1) { 
	       if (idx_2 == index_2) {
		  bonds.erase(bonds.begin()+ibond);
		  status = true;
		  break;
	       }
	    }
	 }
	 return status;
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
