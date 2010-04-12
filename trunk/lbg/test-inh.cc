

#include <math.h>  // for fabs, cos, sin
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>


#define DEG_TO_RAD .01745327 // defined elsewhere maybe.

namespace lig_build {

   class atom_position_t {
   public:
      double x;
      double y;
      atom_position_t(double x_in, double y_in) {
	 x = x_in;
	 y = y_in;
      }
      atom_position_t() {
	 x = -1;
	 y = -1;
      }
      atom_position_t unit_vector() const {
	 double l_sqrd = x*x + y*y;
	 double l = sqrt(l_sqrd);
	 return atom_position_t(x/l, y/l);
      }
      atom_position_t operator+(const atom_position_t &p) const {
	 return atom_position_t(x+p.x, y+p.y);
      }
      atom_position_t operator-(const atom_position_t &p) const {
	 return atom_position_t(x-p.x, y-p.y);
      }
      void operator+=(const atom_position_t &p) {
	 x += p.x;
	 y += p.y;
      }
      void operator-=(const atom_position_t &p) {
	 x -= p.x;
	 y -= p.y;
      }
      bool close_point(const atom_position_t &other) const {
	 bool status = 0;
	 double small = 3;
	 return near_point(other, small);
      }
      bool near_point(const atom_position_t &other, double small) const {
	 bool status = 0;
	 if (fabs(x-other.x) < small) { 
	    if (fabs(y-other.y) < small) {
	       status = 1;
	    }
	 }
	 return status; 
      }
      static atom_position_t fraction_point(const atom_position_t &first,
					    const atom_position_t &second,
					    double frac) {
	 double d_x = second.x - first.x;
	 double d_y = second.y - first.y;
	 return atom_position_t(first.x + frac * d_x,
				first.y + frac * d_y);
      }
      static atom_position_t mid_point(const atom_position_t &first,
				       const atom_position_t &second) {
	 return fraction_point(first, second, 0.5);
      }
      // angle in degrees
      atom_position_t rotate(double angle) const {
	 double theta = angle * DEG_TO_RAD;
	 double sin_theta = sin(theta);
	 double cos_theta = cos(theta);
	 double new_x = x * cos_theta - y * sin_theta;
	 double new_y = x * sin_theta + y * cos_theta;
	 return atom_position_t(new_x, new_y);
      }
      atom_position_t operator*(float sc) {
	 return atom_position_t(x*sc, y*sc);
      }
      double length() const {
	 return sqrt(x*x + y*y);
      } 
      static double length(const atom_position_t &p1, const atom_position_t &p2) {
	 double a = (p2.x - p1.x);
	 double b = (p2.y - p1.y);
	 double c = a*a + b*b;
	 if (c < 0.0) c = 0.0;
	 return sqrt(c);
      }
      // angle between this points (relative to the origin of course).
      static double angle(const atom_position_t &a, const atom_position_t &b) {
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
      // return angle to X axis in degrees
      double axis_orientation() const {
	 double angle = atan2(y,x)/DEG_TO_RAD;
	 return angle;
      }
      bool operator==(const atom_position_t &pos) {
	 if (fabs(pos.x-x) < 0.00001) {
	    
	    if (fabs(pos.x-x) < 0.00001) {
	    } else {
	       return 0;
	    }
	 } else {
	    return 0;
	 } 
      }
   };

   class atom_t {
   public:
      atom_position_t atom_position;
      std::string element;
      std::string name;
      int charge;
      atom_t(atom_position_t pos_in, std::string ele_in, int charge_in) {
	 atom_position = pos_in;
	 name = ele_in;
	 element = ele_in;
	 charge = charge_in;
      }
      bool over_atom(const double &x_in, const double &y_in) const {
	 atom_position_t mouse(x_in, y_in);
	 double d = atom_position_t::length(mouse, atom_position);
	 if (d < 5)
	    return 1;
	 else
	    return 0;
      }
   };

   class bond_t {
   public:
      enum bond_type_t { SINGLE_BOND, DOUBLE_BOND, TRIPLE_BOND, IN_BOND, OUT_BOND };
   private: 
      int atom_1;
      int atom_2;
      bond_type_t bond_type;
      atom_position_t centre_pos; // the position of the polygen of
				  // which this bond is a part
      bool have_centre_pos;  // was the bond from a polygon or just an external bond?
   public:
      bond_t(int first, int second, bond_type_t bt) {
	 atom_1 = first;
	 atom_2 = second;
	 bond_type = bt;
	 have_centre_pos = 0;
      }
      bond_t(int first, int second, atom_position_t centre_pos_in, bond_type_t bt) {
	 atom_1 = first;
	 atom_2 = second;
	 bond_type = bt;
	 have_centre_pos = 1;
	 centre_pos = centre_pos_in;
      }
      // mouse is hovering over bond?
      bool over_bond(double x, double y,
		     const atom_t &atom_1_at, const atom_t &atom_2_at) const;
      int get_atom_1_index() const {
	 return atom_1;
      }
      int get_atom_2_index() const {
	 return atom_2;
      }
   };

   template<class Ta, class Tb> class molecule_t {
   public: 
      std::vector<Ta> atoms;
      std::vector<Tb> bonds;
      virtual Tb get_bond(int bond_index) const {
	 return bonds[bond_index];
      }
      virtual void add_bond(Tb b) {
	 bonds.push_back(b);
      }
   };

   // ----------------------------------------------------------------------------------
   //                     sub-classes
   // ----------------------------------------------------------------------------------

   class w_bond_t : public bond_t {
   public:
      double canvas_item;
      w_bond_t(bond_type_t bt) : bond_t(1,2, bt) {
	 canvas_item = 0;
      }
   };

   class w_atom_t : public atom_t {
   public:
      double canvas_item;
      w_atom_t(atom_position_t pos_in, std::string ele_in, int charge_in) : atom_t(pos_in, ele_in, charge_in) {
	 canvas_item = 0;
      }
   };

   class w_molecule_t : public molecule_t<w_atom_t, w_bond_t> {
   };

}


int main(int argc, char **argv) {

   lig_build::bond_t::bond_type_t bt = lig_build::bond_t::SINGLE_BOND;
   lig_build::w_bond_t wb(bt);

   lig_build::w_molecule_t mol;

   mol.add_bond(wb);

   lig_build::bond_t wb_1 = mol.get_bond(0);

   double x = 800000.0/81.0;
   std::cout << std::setprecision(10);
   std::cout << setiosflags(std::ios::fixed) << std::setw(10) << std::setprecision(4) << x << ":" << std::endl;
      /// std::cout << x << ":" << std::endl;
   return 0;

   return 0;

}

