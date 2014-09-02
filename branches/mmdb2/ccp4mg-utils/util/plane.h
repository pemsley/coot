/*
     util/plane.h: CCP4MG Molecular Graphics Program
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


#ifndef _CCP4MG_PLANE_
#define _CCP4MG_PLANE_

#include <vector>
#include "cartesian.h"
#include "matrix.h"

class Plane{
  double A;
  double B;
  double C;
  double D;
 public:
  Plane();
  ~Plane();
  Plane(double *p1, double *p2, double *p3);
  Plane(const Cartesian &p1, const Cartesian &p2, const Cartesian &p3);
  Plane(double a, double b, double c, double d);
  Plane(const Cartesian &p0);
  void CalculateFromVertices(double *p1, double *p2, double *p3);
  void SetABCD(double a, double b, double c, double d);
  double get_A() const;
  double get_B() const;
  double get_C() const;
  double get_D() const;
  void set_A(double a_in);
  void set_B(double b_in);
  void set_C(double c_in);
  void set_D(double d_in);
  void Normalize(void) ;
  std::vector <Cartesian> find_points_on_plane(void) const;
  Cartesian line_plane_intersection(const Cartesian &p0, const Cartesian &p1) const;
  matrix GetProjectionMatrix() const;
  Cartesian get_normal() const;
  friend std::ostream& operator<<(std::ostream &c, Plane a);
  void Print(void){ std::cout << *this;};
};

#endif
