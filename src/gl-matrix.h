/* src/gl-matrix.h
 * 
 * Copyright 2002,  by The University of York
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

#ifndef HAVE_GL_MATRIX_H
#define HAVE_GL_MATRIX_H

#include <ostream>

enum { TRANSPOSE }; 

#ifdef HAVE_GSL
#include "gsl/gsl_linalg.h"
#endif

#include <clipper/core/coords.h>
#include "coords/Cartesian.h"

class GL_matrix {

   float mat[16];

 public:

   GL_matrix();
   GL_matrix(float m11, float m12, float m13,
	     float m21, float m22, float m23,
	     float m31, float m32, float m33);
   GL_matrix(const clipper::Mat33<double> &m);
   // We want an orientation (doesn't matter which) that is along normal.
   GL_matrix(const clipper::Coord_orth &normal);

   GL_matrix(float quat[4]);
   
   GL_matrix( GL_matrix, int manip);  

   clipper::Mat33<double> to_clipper_mat() const;

   float* operator()() const;

   const float* get() const;

   void rotate_X(float angle);
   void rotate_Y(float angle);
   void rotate_Z(float angle);


   // return a "this is a useful matrix" flag, so that we don't draw
   // elipsoids for atoms with non-positive definite U matrices.
   std::pair<bool,GL_matrix> cholesky() const; 

   GL_matrix transpose() const; 

   float cholesky_non_diag(const GL_matrix &l, int i, int j) const; 

   float cholesky_diag(const GL_matrix &l, int i) const; 


   float squared(float a) const { return a*a; }

   void from_quaternion(float quat[4]); // constructor really, if 
   // we could pass a quaternion, not just a float[4], I would fix
   // it now.

   GL_matrix mat_mult(GL_matrix in) const; 

   coot::Cartesian mult(const coot::Cartesian &in) const;

   void print_matrix() const;

   // irow and icol: 0 -> 2 inclusive
   float matrix_element(int icol, int irow) const;

   friend std::ostream& operator<<(std::ostream&, const GL_matrix &m);

};
std::ostream& operator<<(std::ostream &s, const GL_matrix &m);

void my_aniso_error_handler (const char * reason,
			     const char * file,
			     int line,
			     int gsl_errno);

#endif // HAVE_GL_MATRIX_H
