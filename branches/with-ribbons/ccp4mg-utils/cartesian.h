//   CCP4 Molecular Graphics Program
//
//   Copyright 2004 The University of York
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
//   write to the CCP4 Secretary, Daresbury Laboratory, Warrington
//   WA4 4AD, UK. The GNU Lesser General Public can also be obtained
//   writing to the Free Software Foundation, Inc., 51 Franklin
//   Street, Fifth Floor, Boston, MA 02110-1301, USA


#ifndef _CCP4MG_CART_
#define _CCP4MG_CART_
#include <vector>

#include "matrix.h"

class Cartesian{
  double x;
  double y;
  double z;
  double a;
 public:
  Cartesian();
  Cartesian(double *coords_in);
  Cartesian(const std::vector<double> &coords_in);
  Cartesian(double x_in, double y_in, double z_in, double a_in=1.0);
  double *getxyza(void) const;
  void setxyza(double *coords_in);
  std::vector<double> getxyza_vec(void) const;
  void setxyza_vec(const std::vector<double> &coords_in);
  inline double get_x(void) const {return x;};
  inline double get_y(void) const {return y;};
  inline double get_z(void) const {return z;};
  inline double get_a(void) const {return a;};
  inline void set_x(double x_in){x=x_in;};
  inline void set_y(double y_in){y=y_in;};
  inline void set_z(double z_in){z=z_in;};
  inline void set_a(double a_in){a=a_in;};
  static Cartesian CrossProduct(const Cartesian &a, const Cartesian &b);
  static double DotProduct(const Cartesian &a, const Cartesian &b);
  static Cartesian MidPoint(const Cartesian &v1, const Cartesian &v2);
  static Cartesian MidPoint(const std::vector<Cartesian> &v);
  static std::vector<Cartesian> PrincipalComponentAnalysis(const std::vector<Cartesian> &carts);
  void normalize(double radius=1.0);
  Cartesian operator+(const Cartesian &) const;
  Cartesian operator-(const Cartesian &) const; 
  Cartesian operator-() const; 
  void operator+=(const Cartesian &);
  void operator-=(const Cartesian &);
  friend Cartesian operator* (double val,   const Cartesian&);
  friend Cartesian operator* (const Cartesian &a, double val);
  friend Cartesian operator/ (const Cartesian &a, double val);

  friend Cartesian operator* (matrix objrotmat, const Cartesian &prim);

  Cartesian& operator*= (double val);
  Cartesian& operator/= (double val);
  double length() const;
  double *to_dp() const;

  friend std::ostream& operator<<(std::ostream &c, const Cartesian &a);
  friend std::istream& operator>>(std::istream &c, Cartesian &a);

  void Scale(double a, double b, double c);

};

double Angle(const Cartesian &A, const Cartesian &B, const Cartesian &C);


#endif
