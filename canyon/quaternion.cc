/* canyon/quaternion.cc
 * 
 * Copyright 2014 by Medical Research Council
 * Author: Paul Emsley
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

#include "quaternion.hh"

quaternion::quaternion(const clipper::Coord_orth &co_uv) {

   for (unsigned int i=0; i<3; i++)
      q[i] = co_uv[i];
   q[3] = 0;
} 

std::ostream&
operator<<(std::ostream &s, const quaternion &q) {
   s << q.q[0] << " " << q.q[1] << " " << q.q[2] << " " << q.q[3];
   return s;
} 


// static
quaternion
quaternion::slerp(const quaternion &q1, const quaternion &q2, double fraction) { 

   double dp = quaternion::dot(q1, q2);
   double omega = acos(dp);

   double one_over_sin_omega = 1.0/sin(omega);
   double frac1 = sin((1-fraction)*omega) * one_over_sin_omega;
   double frac2 = sin(fraction*omega) * one_over_sin_omega;
   
   return quaternion(frac1*q1.q[0] + frac2*q2.q[0],
		     frac1*q1.q[1] + frac2*q2.q[1],
		     frac1*q1.q[2] + frac2*q2.q[2],
		     frac1*q1.q[3] + frac2*q2.q[3]);
}


// static
double
quaternion::dot(const quaternion &start_uv, const quaternion &end_uv) {

   return
      start_uv.q[0] * end_uv.q[0] +
      start_uv.q[1] * end_uv.q[1] + 
      start_uv.q[2] * end_uv.q[2] + 
      start_uv.q[3] * end_uv.q[3];

} 


double
quaternion::len3() const {
   return sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
}

double
quaternion::length() const {
   double l = length_sq();
   if (l < 0)
      l = 0;
   return sqrt(l);
}

double
quaternion::length_sq() const {
   return q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
} 

