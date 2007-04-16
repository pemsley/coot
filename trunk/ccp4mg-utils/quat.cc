//   CCP4 Molecular Graphics Program
//
//   Copyright 2004 The University of York
//   Author: Stuart McNicholas and Liz Potterton
//
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


#include <math.h>
#include <iostream>
#include "quat.h"
#include "matrix.h"
#include "cartesian.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

std::ostream& operator<<(std::ostream& c, Quat a){
  c << a.dval[0] << " " << a.dval[1] << " " << a.dval[2] << " " << a.dval[3];
  return c;
}

const double* Quat::Getdval(void) const{
  double *retdval = new double[4];
  retdval[0] = dval[0];
  retdval[1] = dval[1];
  retdval[2] = dval[2];
  retdval[3] = dval[3];
  return retdval;
}

Quat::~Quat(){
  dval.clear();
}

Quat::Quat(const Cartesian& axis,int wi,double angle){
  *this = Quat(axis.get_x(),axis.get_y(),axis.get_z(),wi,angle);
}

Quat::Quat(double x, double y, double z, int wi, double angle){
  dval.clear();
  dval.push_back(0.0);
  dval.push_back(0.0);
  dval.push_back(0.0);
  dval.push_back(0.0);

 dval[0] = 1.0;
 dval[1] = 0.0;
 dval[2] = 0.0;
 dval[3] = 0.0;

 if(wi == 0)
   set(x,y,z);
 if(wi == 1)
   seta(angle,x,y,z);
}

void Quat::reset(void){

 dval[0] = 1.0;
 dval[1] = 0.0;
 dval[2] = 0.0;
 dval[3] = 0.0;

}

void Quat::set(double x,double y,double z){

  Quat xQ = Quat(1.0, 0.0, 0.0, 1, x);
  Quat yQ = Quat(0.0, 1.0, 0.0, 1, y);
  Quat zQ = Quat(0.0, 0.0, 1.0, 1, z);

  Setdval(xQ.dval);
  postMult(yQ);
  postMult(zQ);

}

void Quat::seta(double angle_in, double x_in, double y_in, double z_in){

  double angle = angle_in * M_PI/180.0;
  double radius = sqrt(x_in*x_in + y_in*y_in +z_in*z_in);
  if(radius < 0.000000001){
     std::cout << "zero length in Quat::seta\n";
     return;
  }
  double x = x_in / radius;
  double y = y_in / radius;
  double z = z_in / radius;

  dval[0] = cos(angle/2.0);
  dval[1] = x*sin(angle/2.0);
  dval[2] = y*sin(angle/2.0);
  dval[3] = z*sin(angle/2.0);
}

void Quat::postMult(const Quat &q){

  Quat temp;
  temp.Setdval(dval);

  multAndSet(temp,q);

}


void Quat::multAndSet(const Quat &quat1,const Quat &quat2){
  double dval_new[4];
  std::vector<double> quat1dval = quat1.dval;
  std::vector<double> quat2dval = quat2.dval;
  dval_new[0] = quat2dval[0] * quat1dval[0] - quat2dval[1] * quat1dval[1]
                                               - quat2dval[2] * quat1dval[2]
                                               - quat2dval[3] * quat1dval[3];

  dval_new[1] = quat2dval[0] * quat1dval[1] + quat2dval[1] * quat1dval[0]
                                               + quat2dval[2] * quat1dval[3]
                                               - quat2dval[3] * quat1dval[2];

  dval_new[2] = quat2dval[0] * quat1dval[2] - quat2dval[1] * quat1dval[3]
                                               + quat2dval[2] * quat1dval[0]
                                               + quat2dval[3] * quat1dval[1];

  dval_new[3] = quat2dval[0] * quat1dval[3] + quat2dval[1] * quat1dval[2]
                                               - quat2dval[2] * quat1dval[1]
                                               + quat2dval[3] * quat1dval[0];
  Setdval(dval_new);
}

void Quat::normalize(void){
   double d = sqrt(dval[0]*dval[0]+dval[1]*dval[1]+dval[2]*dval[2]+dval[3]*dval[3]);

   if (fabs(d)<1.0e-12) {
       std::cout << "zero length vector in Quat::normalize\n";
      return;
   }
   double radius = 1.0;

   dval[0] /= (d/radius);
   dval[1] /= (d/radius);
   dval[2] /= (d/radius);
   dval[3] /= (d/radius);

}

matrix Quat::getInvMatrix(void) const{

  Quat quat;
  quat.dval[0] = dval[0];
  quat.dval[1] = -dval[1];
  quat.dval[2] = -dval[2];
  quat.dval[3] = -dval[3];
  quat.normalize();
  return quat.getMatrix();
}

void Quat::rotate_about_axes(const Cartesian& xaxis, const Cartesian& yaxis, const Cartesian& zaxis, double *drot){
  Quat rotx(xaxis.get_x(),xaxis.get_y(),xaxis.get_z(),1,drot[0]);
  Quat roty(yaxis.get_x(),yaxis.get_y(),yaxis.get_z(),1,drot[1]);
  Quat rotz(zaxis.get_x(),zaxis.get_y(),zaxis.get_z(),1,drot[2]);
  postMult(rotx);
  postMult(roty);
  postMult(rotz);
}

matrix Quat::getMatrix(void) const{

  Quat quat = *this;
  quat.normalize();

  double w = quat.dval[0];
  double x = quat.dval[1];
  double y = quat.dval[2];
  double z = quat.dval[3];

  double xx = x*x;
  double yy = y*y;
  double zz = z*z;

  matrix M(4,4);

  M(0,0) = 1.0 - 2.0 * ( yy + zz );
  M(0,1) = 2.0 * ( x * y + w * z );
  M(0,2) = 2.0 * ( x * z - w * y );
  M(0,3) = 0.0;

  M(1,0) = 2.0 * ( x * y - w * z );
  M(1,1) = 1.0 - 2.0 * ( xx + zz );
  M(1,2) = 2.0 * ( y * z + w * x );
  M(1,3) = 0.0;

  M(2,0) = 2.0 * ( x * z + w * y );
  M(2,1) = 2.0 * ( y * z - w * x );
  M(2,2) = 1.0 - 2.0 * ( xx + yy );
  M(2,3) = 0.0;

  M(3,0) = 0.0;
  M(3,1) = 0.0;
  M(3,2) = 0.0;
  M(3,3) = 1.0;

  return M;
}

Quat::Quat(const Cartesian &n, const Cartesian &up){
  /*
   * Quaternion formed from description of scene by outward 'z' vector n
   * and up vector.
   */

  double rotation_angle = acos(Cartesian::DotProduct(n,Cartesian(0,0,1,1)));
  if(fabs(rotation_angle)>1e-7){
    Cartesian rotation_vector = Cartesian::CrossProduct(n,Cartesian(0,0,1,1));
    *this = Quat(rotation_vector,1,rotation_angle * 180.0/M_PI);
  } else {
    *this = Quat(n,1,0.0);
  }

  double up_rotation_angle = acos(Cartesian::DotProduct(up,getMatrix()*Cartesian(0,1,0,1)));
  if(fabs(up_rotation_angle)>1e-7){
    Cartesian up_rotation_vector = Cartesian::CrossProduct(up,getMatrix()*Cartesian(0,1,0,1));
    postMult(Quat(up_rotation_vector,1,up_rotation_angle * 180.0/M_PI));
  }

}

void Quat::Setdval(const std::vector<double> &din){
  dval = din;
}

void Quat::Setdval(const double *din){
  dval[0] = din[0];
  dval[1] = din[1];
  dval[2] = din[2];
  dval[3] = din[3];
}

Quat::Quat(const Quat &qin){
  dval.clear();
  dval.push_back(0.0);
  dval.push_back(0.0);
  dval.push_back(0.0);
  dval.push_back(0.0);
  Setdval(qin.dval);
}

Quat Quat::operator-(const Quat &qin) const {

  return (*this) + -qin;

}

Quat Quat::operator+(const Quat &qin) const {
  std::vector<double> d1 = dval;
  std::vector<double> d2 = qin.dval;
  std::vector<double> dval_new;
  dval_new.push_back(d1[0]+d2[0]);
  dval_new.push_back(d1[1]+d2[1]);
  dval_new.push_back(d1[2]+d2[2]);
  dval_new.push_back(d1[3]+d2[3]);
  Quat qnew;
  qnew.Setdval(dval_new);
  qnew.normalize();
  return qnew;
}

void Quat::operator+=(const Quat &qin){
  *this = *this + qin;
}

void Quat::operator-=(const Quat &qin){
  *this = *this - qin;
}

Quat Quat::operator-() const {
  return -1.0 * (*this);
}

Quat operator*(double val, const Quat &q){
  std::vector<double> d = q.dval;
  std::vector<double> dval_new;
  dval_new.push_back(d[0]*val);
  dval_new.push_back(d[1]*val);
  dval_new.push_back(d[2]*val);
  dval_new.push_back(d[3]*val);
  Quat out;
  out.Setdval(dval_new);
  return out;
}

Quat operator*(const Quat &q, double val){
  return val * q;
}

Quat operator/(const Quat &q, double val){
  return 1.0/val * q;
}

Quat& Quat::operator*= (double val){
  *this = *this * val;
  return *this;
}

Quat& Quat::operator/= (double val){
  *this = *this / val;
  return *this;
}
