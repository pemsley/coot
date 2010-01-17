/*
     util/cartesian.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009 University of York

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


#ifndef _CCP4MG_CART_
#define _CCP4MG_CART_
#include <vector>

#include "matrix.h"

class Cartesian{
  double x;
  double y;
  double z;
  double a;
  static std::vector<Cartesian> PrincipalComponentAnalysis(const std::vector<double> &X, const std::vector<double> &Y, const std::vector<double> &Z);
 public:
  Cartesian();
  Cartesian(const double *coords_in);
  Cartesian(const std::vector<double> &coords_in);
  Cartesian(double x_in, double y_in, double z_in, double a_in=1.0);
  double *getxyza(void) const;
  void setxyza(const double *coords_in);
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
  static std::vector<Cartesian> PrincipalComponentAnalysis(const std::vector<Cartesian>::iterator &i1, const std::vector<Cartesian>::iterator &i2);
  void normalize(double radius=1.0);
  Cartesian operator+(const Cartesian &) const;
  Cartesian operator-(const Cartesian &) const; 
  Cartesian operator-() const; 
  void operator+=(const Cartesian &);
  void operator-=(const Cartesian &);
  friend Cartesian operator* (double val,   const Cartesian&);
  friend Cartesian operator* (const Cartesian &a, double val);
  friend Cartesian operator/ (const Cartesian &a, double val);

  friend Cartesian operator* (const matrix &objrotmat, const Cartesian &prim);

  Cartesian& operator*= (double val);
  Cartesian& operator/= (double val);
  double length() const;
  double *to_dp() const;

  friend std::ostream& operator<<(std::ostream &c, const Cartesian &a);
  friend std::istream& operator>>(std::istream &c, Cartesian &a);

  void Scale(double a, double b, double c);
  void Print(void){ std::cout << *this;};

};

double Angle(const Cartesian &A, const Cartesian &B, const Cartesian &C);


#endif
