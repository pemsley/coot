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
//   write to the CCP4 Secretary, Daresbury Laboratory, Warrington
//   WA4 4AD, UK. The GNU Lesser General Public can also be obtained
//   writing to the Free Software Foundation, Inc., 51 Franklin
//   Street, Fifth Floor, Boston, MA 02110-1301, USA


#ifndef _CCP4MG_GEOMUTIL_
#define _CCP4MG_GEOMUTIL_
#include "cartesian.h"
#include "quat.h"
#include <string>
Cartesian GetCartFrom3Carts(const Cartesian &Atom1, double blength, const Cartesian &Atom2, double angle1, const Cartesian &Atom3, double angle2, int chiral=0);
double LineLength(const Cartesian &at1,  const Cartesian &at2);
double DihedralAngle(const Cartesian &at1,  const Cartesian &at2,  const Cartesian &at3,  const Cartesian &at4);
std::vector<double> DistanceBetweenTwoLines(const Cartesian &a, const Cartesian &b, const Cartesian &c, const Cartesian &d);
std::vector<double> DistanceBetweenPointAndLine(const Cartesian &ls, const Cartesian &le, const Cartesian &p);
Cartesian PointAtWhichTangentToOneLineIntersectsAnotherLine(const Cartesian &p, const Cartesian &ls, const Cartesian &le,  const Cartesian &ols,  const Cartesian &ole);
std::vector<double> LeastSquares2D(const std::vector<Cartesian> &p);
std::vector<Cartesian> LeastSquaresOrtho3D(const std::vector<Cartesian> &p);
Cartesian Mean(const std::vector<Cartesian> &v);
Quat GetStandardRotation(const std::string &s);
#endif
