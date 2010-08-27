/*
     util/quat.h: CCP4MG Molecular Graphics Program
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


#ifndef _QUAT_
#define _QUAT_

#include <vector>
#include "cartesian.h"

class matrix;

class Quat{
   std::vector<double> dval;
   static Quat getSi(const Quat &q, const Quat &qp1, const Quat &qm1);
  public:
   void Setdval(const double *din);
   void Setdval(const std::vector<double> &din);
   const double* Getdval(void) const;
   Quat(double x=0,double y=0,double z=0,int wi=-1,double angle=0);
   Quat(const Cartesian& axis,int wi=-1,double angle=0);
   Quat(const Cartesian &n, const Cartesian &up);
   Quat(const matrix &m);
   Quat Inverse() const;
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
   void Print(void){ std::cout << *this;};
   Quat operator+(const Quat &) const;
   Quat operator-(const Quat &) const; 
   Quat operator-() const; 
   void operator+=(const Quat &);
   void operator-=(const Quat &);
   friend Quat operator* (double val,   const Quat&);
   friend Quat operator* (const Quat &a, double val);
   friend Quat operator/ (const Quat &a, double val);
   friend Quat operator/ (double val, const Quat &a);
 
   Quat& operator*= (double val);
   Quat& operator/= (double val);

   static double Distance(const Quat &q1,const Quat &q2);
   static Quat Slerp(const Quat &q1,const Quat &q2, double h);

   static std::vector<Quat> Squad(const std::vector<Quat> &qs, unsigned accu);

   Quat operator/(const Quat &) const;
   void operator/=(const Quat &);
   Quat operator*(const Quat &) const;
   void operator*=(const Quat &);
   Quat Conjugate() const;
   static double DotProduct(const Quat &q1,const Quat &q2);
   static Quat OuterProduct(const Quat &q1,const Quat &q2);
   Quat Log() const;
   Quat Exp() const;
   Quat Power(double p) const;
   double Norm() const;

};
#endif
