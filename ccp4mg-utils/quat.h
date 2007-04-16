//   CCP4 Molecular Graphics Program
//
//   Copyright 2004 The University of York
//   Author: Stuart McNicholas and Liz Potterton
//
//   This program is free software and is distributed under the terms
//   and conditions of the CCP4 licence agreement as `Part 0' (Annex 2)
//   software, which is version 2.1 of the GNU Lesser General Public
//   Licence (LGPL) with the following additional clause:
//
//      `You may also combine or link a "work that uses the Library"
//      to produce a work containing portions of the Library, and
//      distribute that work under terms of your choice, provided that
//      you give prominent notice with each copy of the work that the
//      specified version of the Library is used in it, and that you
//      include or provide public access to the complete corresponding
//      machine-readable source code for the Library including whatever
//      changes were used in the work. (i.e. If you make changes to the
//      Library you must distribute those, but you do not need to
//      distribute source or object code to those portions of the work
//      not covered by this licence.)'
//
//   Note that this clause grants an additional right and does not
//   impose any additional restriction, and so does not affect
//   compatibility with the GNU General Public Licence (GPL). If you
//   wish to negotiate other terms, please contact the maintainer.
//   You can redistribute it and/or modify the program under the terms
//   of the GNU Lesser General Public License as published by the Free
//   Software Foundation; either version 2.1 of the License, or (at
//   your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//   Lesser General Public License for more details.
//
//   You should have received a copy of the CCP4 licence and/or GNU
//   Lesser General Public License along with this program; if not,
//   write to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4
//   4AD, UK. The GNU Lesser General Public can also be obtained by
//   writing to the Free Software Foundation, Inc., 51 Franklin
//   Street, Fifth Floor, Boston, MA 02110-1301, USA


#ifndef _QUAT_
#define _QUAT_

#include <vector>
#include "cartesian.h"

class matrix;

class Quat{
   std::vector<double> dval;
  public:
   void Setdval(const double *din);
   void Setdval(const std::vector<double> &din);
   const double* Getdval(void) const;
   Quat(double x=0,double y=0,double z=0,int wi=-1,double angle=0);
   Quat(const Cartesian& axis,int wi=-1,double angle=0);
   Quat(const Cartesian &n, const Cartesian &up);
   ~Quat();
   void reset(void);
   void set(double x,double y,double z);
   Quat(const Quat &qfrom);
   void seta(double angle_in,double x_in,double y_in,double z_in);
   void postMult(const Quat &q);
   void multAndSet(const Quat &quat1,const Quat &quat2);
   void normalize(void);
   matrix getInvMatrix(void) const;
   matrix getMatrix(void) const;
   void rotate_about_axes(const Cartesian& xaxis, const Cartesian& yaxis, const Cartesian& zaxis, double *drot);
   friend std::ostream& operator<<(std::ostream &c, Quat a);
   Quat operator+(const Quat &) const;
   Quat operator-(const Quat &) const; 
   Quat operator-() const; 
   void operator+=(const Quat &);
   void operator-=(const Quat &);
   friend Quat operator* (double val,   const Quat&);
   friend Quat operator* (const Quat &a, double val);
   friend Quat operator/ (const Quat &a, double val);
 
   Quat& operator*= (double val);
   Quat& operator/= (double val);

};
#endif
