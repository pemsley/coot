
#include "gl-matrix.h"

#include <ostream>

class quaternion {

   double q[4];
public:
   quaternion(const double &q0,
	      const double &q1,
	      const double &q2,
	      const double &q3) {

      q[0] = q0;
      q[1] = q1;
      q[2] = q2;
      q[3] = q3;
   }
   // co_uv is unit vector, phi in radians
   quaternion(const clipper::Coord_orth &co_uv);
   quaternion add(const quaternion &q_in);
   GL_matrix build_matrix();

   static quaternion slerp(const quaternion &start, const quaternion &end, double fraction);

   static double dot(const quaternion &start_uv, const quaternion &end_uv);

   friend std::ostream &operator<<(std::ostream &s, const quaternion &q);
   double len3() const;
   double length() const;
   double length_sq() const;
   clipper::Coord_orth first3() const { return clipper::Coord_orth(q[0], q[1], q[2]); }

};
std::ostream &operator<<(std::ostream &s, const quaternion &q);


