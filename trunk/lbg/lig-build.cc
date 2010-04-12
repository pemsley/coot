
#include "lbg.hh"


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

