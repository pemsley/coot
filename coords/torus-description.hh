

#ifndef COORDS_TORUS_DESCRIPTION_HH
#define COORDS_TORUS_DESCRIPTION_HH

namespace coot {

   // For OpenGL solid model, it is much better to draw a torus than a
   // set of short sticks (particularly for the ring representing
   // aromaticity).  So now (20100831 Bond_lines_container contains a
   // number of torus descriptions).
   //
   class torus_description_t {
   public:
      double inner_radius;
      double outer_radius;
      int n_sides;
      int n_rings;
      clipper::Coord_orth centre;
      clipper::Coord_orth normal;
      torus_description_t(const clipper::Coord_orth &pt,
			  const clipper::Coord_orth &normal_in,
			  double ir1, double ir2, int n1, int n2) {
	 inner_radius = ir1;
	 outer_radius = ir2;
	 n_sides = n1;
	 n_rings = n2;
	 centre = pt;
	 normal = normal_in;
      }
   };

}

#endif // COORDS_TORUS_DESCRIPTION_HH
