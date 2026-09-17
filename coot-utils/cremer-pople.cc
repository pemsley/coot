#include <cmath>
#include <limits>
#include "cremer-pople.hh"

coot::cremer_pople_t::cremer_pople_t(const std::vector<clipper::Coord_orth> &ring,
                                     const clipper::Coord_orth *up_reference) {

   filled = false;
   n_ring = static_cast<unsigned int>(ring.size());
   Q = 0.0; q2 = 0.0; q3 = 0.0;
   theta = std::numeric_limits<double>::quiet_NaN();
   phi = std::numeric_limits<double>::quiet_NaN();
   normal = clipper::Coord_orth(0,0,0);

   if (n_ring != 5 && n_ring != 6) return;
   const double N = static_cast<double>(n_ring);

   // 1. Centre on the ring centroid.
   clipper::Coord_orth centre(0,0,0);
   for (const auto &p : ring) centre += p;
   centre = clipper::Coord_orth(centre.x()/N, centre.y()/N, centre.z()/N);

   std::vector<clipper::Coord_orth> R;
   for (const auto &p : ring) R.push_back(p - centre);

   // 2. The Cremer-Pople mean plane. NOT an LSQ plane: the CP plane is defined
   //    by sum(z) = 0 AND sum(z cos) = 0 AND sum(z sin) = 0, which LSQ does not
   //    satisfy. R1 x R2 gives a normal obeying all three by construction.
   clipper::Coord_orth R1(0,0,0), R2(0,0,0);
   for (unsigned int j=0; j<n_ring; j++) {
      double a = 2.0 * M_PI * static_cast<double>(j) / N;
      R1 += clipper::Coord_orth(R[j].x()*sin(a), R[j].y()*sin(a), R[j].z()*sin(a));
      R2 += clipper::Coord_orth(R[j].x()*cos(a), R[j].y()*cos(a), R[j].z()*cos(a));
   }
   clipper::Coord_orth n(clipper::Coord_orth::cross(R1, R2));
   double n_len = std::sqrt(n.lengthsq());
   if (n_len < 1e-12) return;   // degenerate (collinear) ring
   n = clipper::Coord_orth(n.x()/n_len, n.y()/n_len, n.z()/n_len);

   // 3. Face selection. Note that reversing the ring traversal direction flips
   //    R1 (sine weights are antisymmetric) and hence flips n -- so under the
   //    right-hand-rule default, reversal maps theta to 180-theta. An
   //    up_reference pins the face from geometry instead, making theta
   //    reversal-invariant and letting enantiomers superimpose.
   if (up_reference) {
      clipper::Coord_orth u = *up_reference - centre;
      if ((u.x()*n.x() + u.y()*n.y() + u.z()*n.z()) < 0.0)
         n = clipper::Coord_orth(-n.x(), -n.y(), -n.z());
   }
   normal = n;

   // 4. Perpendicular displacements.
   z.clear();
   for (unsigned int j=0; j<n_ring; j++)
      z.push_back(R[j].x()*n.x() + R[j].y()*n.y() + R[j].z()*n.z());

   // 5. The m=2 component (present for both N=5 and N=6).
   double C = 0.0, S = 0.0;
   for (unsigned int j=0; j<n_ring; j++) {
      double a = 4.0 * M_PI * static_cast<double>(j) / N;   // m=2 -> 2*(2 pi j / N)
      C += z[j] * cos(a);
      S += z[j] * sin(a);
   }
   double k2 = std::sqrt(2.0/N);
   q2  = k2 * std::sqrt(C*C + S*S);
   phi = atan2(-S, C);
   if (phi < 0.0) phi += 2.0 * M_PI;

   if (n_ring == 6) {
      // 6. The m=3 component exists only for N=6. Its (-1)^j parity is why an
      //    odd rotation of the start atom flips q3, and so swaps 4C1 and 1C4.
      double A = 0.0;
      for (unsigned int j=0; j<6; j++) A += z[j] * ((j % 2 == 0) ? 1.0 : -1.0);
      q3 = std::sqrt(1.0/6.0) * A;
      Q  = std::sqrt(q2*q2 + q3*q3);
      theta = atan2(q2, q3);          // q2 >= 0, so theta is in [0, pi]
   } else {
      q3 = 0.0;
      Q  = q2;                       // no m=3 term: a 5-ring lives on a circle
      // theta stays NaN
   }

   filled = true;
}
