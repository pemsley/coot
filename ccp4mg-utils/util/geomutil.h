/*
     util/geomutil.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
*/


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
std::vector<double> LeastSquaresPolyFit(const std::vector<double> &xs, std::vector<double> &ys, const int order=1);
std::vector<double> LeastSquaresQuadraticFit3D(const std::vector<Cartesian> &carts);
std::vector<double> LeastSquaresCubicFit3D(const std::vector<Cartesian> &carts);
#endif
