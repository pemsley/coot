/*
     util/plane.cc: CCP4MG Molecular Graphics Program
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


#include "plane.h"
#include <stdio.h>
#include <vector>
#include <math.h>

std::ostream& operator<<(std::ostream& c, Plane a){
  c << a.get_A() << " " << a.get_B() << " " << a.get_C() << " " << a.get_D();
  return c;
}

matrix Plane::GetProjectionMatrix() const{
  matrix projmatrix(4,4);
  projmatrix.Zero();
  projmatrix(3,3) =  1.0;
  projmatrix(0,0) =  get_B() * get_B() + get_C() * get_C();
  projmatrix(0,1) = -get_A() * get_B();
  projmatrix(0,2) = -get_A() * get_C();
  projmatrix(1,0) = -get_A() * get_B();
  projmatrix(2,0) = -get_A() * get_C();
  projmatrix(1,2) = -get_B() * get_C();
  projmatrix(2,1) = -get_B() * get_C();
  projmatrix(1,1) =  get_A() * get_A() + get_C() * get_C();
  projmatrix(2,2) =  get_A() * get_A() + get_B() * get_B();
  return projmatrix;
}

Plane::Plane(const Cartesian &p0){
  set_A(p0.get_x());
  set_B(p0.get_y());
  set_C(p0.get_z());
  set_D(0.0);
}

void Plane::Normalize(void){
  Cartesian n = get_normal();
  double fac = n.length();
  n.normalize();
  set_A(n.get_x());
  set_B(n.get_y());
  set_C(n.get_z());
  set_D(get_D()/fac);
}

Plane::Plane(){
}

Plane::Plane(double *p1, double *p2, double *p3){
  CalculateFromVertices(p1,p2,p3);
}

Plane::Plane(const Cartesian &p1, const Cartesian &p2, const Cartesian &p3){
  double *p1_dp = p1.to_dp();
  double *p2_dp = p2.to_dp();
  double *p3_dp = p3.to_dp();
  CalculateFromVertices(p1_dp,p2_dp,p3_dp);
  delete [] p1_dp;
  delete [] p2_dp;
  delete [] p3_dp;
}

void Plane::CalculateFromVertices(double *p1, double *p2, double *p3){
  A = p1[1] * (p2[2] - p3[2]) + p2[1] * (p3[2] - p1[2]) + p3[1] * (p1[2] - p2[2]);
  B = p1[2] * (p2[0] - p3[0]) + p2[2] * (p3[0] - p1[0]) + p3[2] * (p1[0] - p2[0]);
  C = p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]); 
  D = -(p1[0] * (p2[1] * p3[2] - p3[1] * p2[2]) + p2[0] * (p3[1] * p1[2] - p1[1] * p3[2]) + p3[0] * (p1[1] * p2[2] - p2[1] * p1[2]));
}

void Plane::set_A(double a_in){
  A = a_in;
}
void Plane::set_B(double b_in){
  B = b_in;
}
void Plane::set_C(double c_in){
  C = c_in;
}
void Plane::set_D(double d_in){
  D = d_in;
}

void Plane::SetABCD(double a, double b, double c, double d){
  A = a; B = b; C = c; D = d;
}

Plane::Plane(double a, double b, double c, double d){
  A = a; B = b; C = c; D = d;
}

Plane::~Plane(){
}

double Plane::get_A(void) const{
  return A;
}

double Plane::get_B(void) const{
  return B;
}

double Plane::get_C(void) const{
  return C;
}

double Plane::get_D(void) const{
  return D;
}

Cartesian Plane::line_plane_intersection(const Cartesian &p0, const Cartesian &p1) const{
  Cartesian intersec;

  Cartesian p = p1-p0;

  double num = get_A()*p0.get_x() + get_B()*p0.get_y() + get_C()*p0.get_z() + get_D();

  double den = get_A()*p.get_x() + get_B()*p.get_y() + get_C()*p.get_z();

  if(fabs(den) < 0.00000001){
     return Cartesian(1e+7,1e+7,1e+7,1e+7);// This is dodgy, better way?
  }
  double u = -num/den;
  intersec = p0 + u*(p1-p0);
  intersec.set_a(1.0);
  return intersec;
}

std::vector <Cartesian> Plane::find_points_on_plane() const{
  // Find three (totally arbitrary!) points on this plane.
  std::vector<Cartesian> retval;

  Cartesian norm  = get_normal();

  Cartesian p0 = line_plane_intersection(Cartesian(0,0,0),norm);

  Cartesian xaxis(1,0,0);
  Cartesian yaxis(0,1,0);
  Cartesian zaxis(0,0,1);

  norm.normalize();

  double ax = fabs(norm.DotProduct(norm,xaxis));
  double ay = fabs(norm.DotProduct(norm,yaxis));
  double az = fabs(norm.DotProduct(norm,zaxis));

  double maxa = ax;
  Cartesian axis = xaxis;
  if(ay < maxa){
       maxa = ay;
       axis = yaxis;
  }
  if(az < maxa){
       maxa = az;
       axis = zaxis;
  }

  Cartesian v = norm.CrossProduct(norm,axis);


  Cartesian p1 = p0 + v;
  Cartesian v2 = norm.CrossProduct(norm,v);
  Cartesian p2 = p0 + v2;

  double opfac = 1.0/sqrt(get_A()*get_A() +
          get_B()*get_B() +
          get_C()*get_C());

  Plane new_plane = Plane(p0,p1,p2);
  double npfac = 1.0/sqrt(new_plane.get_A()*new_plane.get_A() +
          new_plane.get_B()*new_plane.get_B() +
          new_plane.get_C()*new_plane.get_C());

  if(opfac/npfac < 0.0){
    retval.push_back(Cartesian(p1));
    retval.push_back(Cartesian(p0));
    retval.push_back(Cartesian(p2));
  }else{
    retval.push_back(Cartesian(p0));
    retval.push_back(Cartesian(p1));
    retval.push_back(Cartesian(p2));
  }

  return retval;
}

Cartesian Plane::get_normal() const{
  return Cartesian(A,B,C);
}
