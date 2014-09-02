/*
     util/quat.cc: CCP4MG Molecular Graphics Program
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

Quat::Quat(const matrix &m_in){
  matrix m = m_in.Transpose(); // I should fix code below so that this is not necessary ...
  double trace = m(0,0) + m(1,1) + m(2,2) + m(3,3);
  dval.push_back(0.0);
  dval.push_back(0.0);
  dval.push_back(0.0);
  dval.push_back(0.0);

  if( trace > 1e-7 ) {
    double s = 0.5 / sqrt(trace);
    dval[0] = 0.25 / s;
    dval[1] = ( m(2,1)- m(1,2)) * s;
    dval[2] = ( m(0,2)- m(2,0)) * s;
    dval[3] = ( m(1,0)- m(0,1)) * s;
  } else {
    if ( m(0,0)> m(1,1)&& m(0,0)> m(2,2)) {
      double s = 2.0 * sqrt( 1.0 + m(0,0)- m(1,1)- m(2,2));
      dval[1] = 0.25 * s;
      dval[2] = (m(0,1)+ m(1,0)) / s;
      dval[3] = (m(0,2)+ m(2,0)) / s;
      dval[0] = (m(1,2)- m(2,1)) / s;

    } else if (m(1,1)> m(2,2)) {
      double s = 2.0 * sqrt( 1.0 + m(1,1)- m(0,0)- m(2,2));
      dval[1] = (m(0,1)+ m(1,0)) / s;
      dval[2] = 0.25 * s;
      dval[3] = (m(1,2)+ m(2,1)) / s;
      dval[0] = (m(0,2)- m(2,0)) / s;
    } else {
      double s = 2.0 * sqrt( 1.0 + m(2,2)- m(0,0)- m(1,1));
      dval[1] = (m(0,2)+ m(2,0)) / s;
      dval[2] = (m(1,2)+ m(2,1)) / s;
      dval[3] = 0.25 * s;
      dval[0] = (m(0,1)- m(1,0)) / s;
    }
  }
}

Quat::Quat(const Cartesian &n, const Cartesian &up){
  /*
   * Quaternion formed from description of scene by outward 'z' vector n
   * and up vector.
   */

  double rotation_angle = acos(Cartesian::DotProduct(n,Cartesian(0,0,1,1)));
  if(fabs(rotation_angle)>1e-7){
    Cartesian rotation_vector = Cartesian::CrossProduct(n,Cartesian(0,0,1,1));
    rotation_vector.normalize();
    *this = Quat(rotation_vector,1,rotation_angle * 180.0/M_PI);
  } else {
    *this = Quat(n,1,0.0);
  }

  double up_rotation_angle = acos(Cartesian::DotProduct(up,getMatrix()*Cartesian(0,1,0,1)));
  if(fabs(up_rotation_angle)>1e-7){
    Cartesian cross = getInvMatrix()*Cartesian::CrossProduct(getMatrix()*Cartesian(0,1,0,1),up);
    Cartesian up_rotation_vector(0,0,1);
    if(cross.get_z()<0.0)
      postMult(Quat(up_rotation_vector,1,up_rotation_angle * 180.0/M_PI));
    else
      postMult(Quat(up_rotation_vector,1,-up_rotation_angle * 180.0/M_PI));
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
  //qnew.normalize();
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

Quat Quat::Inverse() const {
  double modqsq = dval[0] * dval[0] + dval[1] * dval[1] + dval[2] * dval[2] + dval[3] * dval[3];
  Quat conjugate;
  conjugate.dval[0] = dval[0];
  conjugate.dval[1] = -dval[1];
  conjugate.dval[2] = -dval[2];
  conjugate.dval[3] = -dval[3];
  return conjugate/modqsq;
}

Quat operator/ (double val, const Quat &a){
  return val * a.Inverse();
}

double Quat::Distance(const Quat &q1,const Quat &q2){
  Quat diff = q1-q2;
  //double q1l = q1.dval[0] * q1.dval[0] + q1.dval[1] * q1.dval[1] + q1.dval[2] * q1.dval[2] + q1.dval[3] * q1.dval[3];
  //double q2l = q2.dval[0] * q2.dval[0] + q2.dval[1] * q2.dval[1] + q2.dval[2] * q2.dval[2] + q2.dval[3] * q2.dval[3];
  //double diffl = diff.dval[0] * diff.dval[0] + diff.dval[1] * diff.dval[1] + diff.dval[2] * diff.dval[2] + diff.dval[3] * diff.dval[3];
  //std::cout << "Q1: " << q1.dval[0] << " "  << q1.dval[1] << " " << q1.dval[2] << " " << q1.dval[3] << "(" << q1l << ")\n";
  //std::cout << "Q2: " << q2.dval[0] << " "  << q2.dval[1] << " " << q2.dval[2] << " " << q2.dval[3] << "(" << q2l << ")\n";
  //std::cout << "DI: " << diff.dval[0] << " "  << diff.dval[1] << " " << diff.dval[2] << " " << diff.dval[3] << "(" << diffl << ")\n";
  double lengthsq = diff.dval[0] * diff.dval[0] + diff.dval[1] * diff.dval[1] + diff.dval[2] * diff.dval[2] + diff.dval[3] * diff.dval[3];
  return sqrt(lengthsq);
}

Quat Quat::Conjugate() const {
  Quat conjugate;
  conjugate.dval[0] = dval[0];
  conjugate.dval[1] = -dval[1];
  conjugate.dval[2] = -dval[2];
  conjugate.dval[3] = -dval[3];
  return conjugate;
}

Quat Quat::operator*(const Quat &qin) const {
  double dval_new[4];
  std::vector<double> quat1dval = this->dval;
  std::vector<double> quat2dval = qin.dval;
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
  Quat qnew;
  qnew.Setdval(dval_new);
  return qnew;
}

void Quat::operator*=(const Quat &qin){
  *this = *this * qin;
}

Quat Quat::operator/(const Quat &qin) const {

  Quat qinv = qin.Inverse();
  return *this * qinv;

}

void Quat::operator/=(const Quat &qin){
  *this = *this / qin;
}

double Quat::DotProduct(const Quat &q1,const Quat &q2){
  return q1.dval[0] * q2.dval[0] + q1.dval[1] * q2.dval[1] + q1.dval[2] * q2.dval[2] + q1.dval[3] * q2.dval[3];
}

Quat Quat::OuterProduct(const Quat &q1,const Quat &q2){
  return (q1.Conjugate()*q2 - q2.Conjugate()*q1)/2.0;
}

Quat Quat::Log() const{
  /* Only valid for unit quaternion. */
  double dval_new[4];
  dval_new[0] = 0.0;
  double theta = acos(dval[0]);
  dval_new[1] = dval[1] * theta/sin(theta);
  dval_new[2] = dval[2] * theta/sin(theta);
  dval_new[3] = dval[3] * theta/sin(theta);
  Quat qnew;
  qnew.Setdval(dval_new);
  return qnew;
}

double Quat::Norm() const {
  double modsq = dval[0] * dval[0] + dval[1] * dval[1] + dval[2] * dval[2] + dval[3] * dval[3];
  return sqrt(modsq);
}

Quat Quat::Exp() const{

  Quat qnew;
  double dval_new[4];

  double l = sqrt(dval[1]*dval[1]+dval[2]*dval[2]+dval[3]*dval[3]);

  if(fabs(l)<1e-7){
    //dval_new[0] = exp(dval_new[0]); // Old surely wrong
    dval_new[0] = exp(dval[0]);
    qnew.Setdval(dval_new);
    return qnew;
  }
    
  double s = sin(l);
  double c = cos(l);
  double e = exp(dval[0]);

  double t = e * s / l;

  dval_new[0] = e*c;
  dval_new[1] = t * dval[1];
  dval_new[2] = t * dval[2];
  dval_new[3] = t * dval[3];

  qnew.Setdval(dval_new);
  return qnew;
}

Quat Quat::Power(double p) const{
  return (p*Log()).Exp();
}

Quat Quat::Slerp(const Quat &q1,const Quat &q2, double h){
  double cosw = DotProduct(q1,q2);
  if(cosw>1.0) cosw = 1.0;
  if(cosw<-1.0) cosw = -1.0;

  double omega = acos(cosw);

  return (q1*sin((1.0-h)*omega) + q2*sin(h*omega))/sin(omega);
}

Quat Quat::getSi(const Quat &q, const Quat &qp1, const Quat &qm1){
  Quat si = q*(-((q.Inverse()*(qp1)).Log()+((q.Inverse())*(qm1)).Log())/4.0).Exp();
  return si;
}

std::vector<Quat> Quat::Squad(const std::vector<Quat> &qs, unsigned accu){
  std::vector<Quat> slerp;
  std::vector<Quat>::const_iterator qiter = qs.begin();
  while(qiter!=qs.end()-1){
    for(unsigned i=0;i<accu;i++){
      double h = double(i)/double(accu);
      Quat si;
      Quat sip1;
      if(qiter==qs.begin())
        si = *qiter;
      else
        si = getSi(*qiter,*(qiter+1),*(qiter-1));
      if(qiter==qs.end()-2)
        sip1 = *(qs.end()-1);
      else
        sip1 = getSi(*(qiter+1),*(qiter+2),*qiter);
      slerp.push_back(Slerp(Slerp(*qiter,*(qiter+1),h),Slerp(si,sip1,h),2*h*(1-h)));
    }
    qiter++;
  }
  return slerp;
}

