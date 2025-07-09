/*
 * src/matrix-utils.cc
 *
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include "matrix-utils.hh"

clipper::Mat33<double> glm_to_mat33(const glm::mat4 &quat_mat) {

   clipper::Mat33<double> m(quat_mat[0][0], quat_mat[0][1], quat_mat[0][2],
                            quat_mat[1][0], quat_mat[1][1], quat_mat[1][2],
                            quat_mat[2][0], quat_mat[2][1], quat_mat[2][2]);
   return m;

}

glm::quat coot_quaternion_to_glm(const coot::util::quaternion &q) {

   return glm::quat(q.q3, q.q0, q.q1, q.q2);

}

coot::util::quaternion glm_to_coot_quaternion(const glm::quat &q) {

   coot::util::quaternion qq(q[3], q[0], q[1], q[2]);
   return qq;
}

// return quaternion as glm format (w (x,y,z))
// from https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2015/01/matrix-to-quat.pdf
glm::quat matrix_to_quaternion(const clipper::Mat33<double> &m) {

   glm::quat q;
   double t;
   double m00=m(0,0);
   double m01=m(0,1);
   double m02=m(0,2);
   double m10=m(1,0);
   double m11=m(1,1);
   double m12=m(1,2);
   double m20=m(2,0);
   double m21=m(2,1);
   double m22=m(2,2);

   if (m22 < 0) {
      if (m00 > m11) {
         t = 1 + m00 - m11 - m22;
         q = glm::quat(m12-m21, t, m01+m10, m20+m02);
      } else {
         t = 1 - m00 + m11 - m22;
         q = glm::quat(m20-m02, m01+m10, t, m12+m21);
      }
   } else {
      if (m00 < -m11) {
         t = 1 - m00 - m11 + m22;
         q = glm::quat(m01-m10, m20+m02, m12+m21, t);
      } else {
         t = 1 + m00 + m11 + m22;
         q = glm::quat(t, m12-m21, m20-m02, m01-m10);
      }
   }
   q *= 0.5 / sqrt(t);
   return q;
}
