/*
     util/cartesian.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York
     Copyright (C) 2012 STFC

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


#include "cartesian.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "matrix.h"
#include <stdlib.h>

std::istream& operator>>(std::istream& c, Cartesian &a){

  c >> a.x >> a.y >> a.z;

  return c;
}

std::ostream& operator<<(std::ostream& c, const Cartesian &a){

  c.flags(std::ios::fixed|std::ios::right);
  c.precision(6);

  c << std::setw(10) << a.get_x() << " " << std::setw(10) << a.get_y() << " " << std::setw(10) << a.get_z() << " " << std::setw(10) << a.get_a();
  
  return c;
}

double Angle(const Cartesian &A, const Cartesian &B, const Cartesian &C){

  Cartesian BA = B-A;
  Cartesian CA = C-A;
  Cartesian CB = C-B;

  double ab = BA.length();
  double ac = CA.length();
  double bc = CB.length();

  double absq = ab*ab;
  double acsq = ac*ac;
  double bcsq = bc*bc;

  return  acos((bcsq + absq - acsq)/(2*bc*ab));
}

void Cartesian::operator*=(const std::vector<double> &objrotmat) {
  double result[4];
  double input[4] = {x, y, z, a};

  for(unsigned i=0;i<4;i++){
    result[i] = 0.0;
    for(unsigned j=0;j<4;j++){
      result[i] += input[j]*objrotmat[i*4+j];
    }
  }
  x = result[0];
  y = result[1];
  z = result[2];
  a = result[3];
}

Cartesian operator*(const std::vector<double> &objrotmat, const Cartesian &prim){
  double result[4];
  double input[4] = {prim.x, prim.y, prim.z, prim.a};

  for(unsigned i=0;i<4;i++){
    result[i] = 0.0;
    for(unsigned j=0;j<4;j++){
      result[i] += input[j]*objrotmat[i*4+j];
    }
  }
  
  return Cartesian(result);
}

Cartesian operator*(const matrix &objrotmat, const Cartesian &prim){
  double result[4];
  double input[4] = {prim.x, prim.y, prim.z, prim.a};

  for(unsigned i=0;i<objrotmat.get_columns();i++){
    result[i] = 0.0;
    for(unsigned j=0;j<objrotmat.get_rows();j++){
      result[i] += input[j]*objrotmat(i,j);
    }
  }
  
  return Cartesian(result);
}


Cartesian operator* (double val, const Cartesian &a){
  
  Cartesian c(a.x*val,a.y*val,a.z*val);
  return c;
}

Cartesian operator* (const Cartesian &a, double val){
  
  Cartesian c(a.x*val,a.y*val,a.z*val);
  return c;
}

Cartesian operator/ (const Cartesian &a, double val){
  
  Cartesian c(a.x/val,a.y/val,a.z/val);
  return c;
}

Cartesian& Cartesian::operator*= (double val){
  
  //*this = *this * val;
  x *= val;
  y *= val;
  z *= val;
  return *this;
}

Cartesian& Cartesian::operator/= (double val){
  
  //*this = *this / val;
  x /= val;
  y /= val;
  z /= val;
  return *this;
}

Cartesian::Cartesian(){

  x = y = z = 0.0;
  a = 1.0;
}

Cartesian::Cartesian(const double *coords_in){
  x = coords_in[0];
  y = coords_in[1];
  z = coords_in[2];
  a = 1.0;
}

Cartesian::Cartesian(const std::vector<double> &coords_in){
  x = coords_in[0];
  y = coords_in[1];
  z = coords_in[2];
  if(coords_in.size()>3)
    a = coords_in[3];
  else
    a = 1.0;
  
}

Cartesian::Cartesian(double x_in, double y_in, double z_in, double a_in){
  x = x_in; y = y_in; z= z_in; a = a_in;
}

double *Cartesian::getxyza(void) const{
  double *result = new double[4];
  result[0] = x;
  result[1] = y;
  result[2] = z;
  result[3] = a;
  return result;
}

void Cartesian::setxyza(const double *coords_in){
  x = coords_in[0];
  y = coords_in[1];
  z = coords_in[2];
  a = coords_in[3];
}  

std::vector<double> Cartesian::getxyza_vec(void) const{
  std::vector <double>result;
  result.push_back(x);
  result.push_back(y);
  result.push_back(z);
  result.push_back(a);
  return result;
}

void Cartesian::setxyza_vec(const std::vector<double> &coords_in){
  x = coords_in[0];
  y = coords_in[1];
  z = coords_in[2];
  a = coords_in[3];
}

double Covariance(const std::vector<double> &X, const std::vector<double> &Y, const double Xmean, const double Ymean){

  if(X.size()!=Y.size()){
    std::cerr << "Error in calculating covarience " << X.size() << " != " << Y.size() << "\n";
    return 1e-7;
  }
  if(X.size()<2){
    std::cerr << "Error in calculating covarience, size of problem < 2\n";
  }


  double cov = 0.0;

  for(unsigned i=0;i<X.size();i++){
    cov += (X[i]-Xmean)*(Y[i]-Ymean);
  }
  cov /= (X.size()-1);

  return cov;
}

std::vector<Cartesian> Cartesian::PrincipalComponentAnalysis(const std::vector<Cartesian>::iterator &i1, const std::vector<Cartesian>::iterator &i2){
  std::vector<double> X;
  std::vector<double> Y;
  std::vector<double> Z;
  std::vector<Cartesian>::iterator iter=i1;
  while(iter!=i2+1){
    X.push_back(iter->x);
    Y.push_back(iter->y);
    Z.push_back(iter->z);
    iter++;
  }
  return PrincipalComponentAnalysis(X,Y,Z);
}

std::vector<Cartesian> Cartesian::PrincipalComponentAnalysis(const std::vector<Cartesian> &carts){
  std::vector<double> X;
  std::vector<double> Y;
  std::vector<double> Z;
  for(unsigned i=0;i<carts.size();i++){
    X.push_back(carts[i].x);
    Y.push_back(carts[i].y);
    Z.push_back(carts[i].z);
  }
  return PrincipalComponentAnalysis(X,Y,Z);
}
  
std::vector<Cartesian> Cartesian::PrincipalComponentAnalysis(const std::vector<double> &X, const std::vector<double> &Y, const std::vector<double> &Z){

  double Xmean = 0.0;
  double Ymean = 0.0;
  double Zmean = 0.0;

  for(unsigned i=0;i<X.size();i++){
    Xmean += X[i];
    Ymean += Y[i];
    Zmean += Z[i];
  }
  Xmean /= X.size();
  Ymean /= X.size();
  Zmean /= Z.size();
  

  std::vector<double> cov;
  double covXX = Covariance(X,X,Xmean,Xmean);
  double covYY = Covariance(Y,Y,Ymean,Ymean);
  double covZZ = Covariance(Z,Z,Zmean,Zmean);
  double covXY = Covariance(X,Y,Xmean,Ymean);
  double covXZ = Covariance(X,Z,Xmean,Zmean);
  double covYZ = Covariance(Y,Z,Xmean,Zmean);
  cov.push_back(covXX);
  cov.push_back(covXY);
  cov.push_back(covXZ);
  cov.push_back(covXY);
  cov.push_back(covYY);
  cov.push_back(covYZ);
  cov.push_back(covXZ);
  cov.push_back(covYZ);
  cov.push_back(covZZ);
  
  matrix covmat(3,3,cov);
  std::vector<matrix> eigen = covmat.Eigen();
  eigen = eigen[0].SortEigenvalues(eigen); // Sort should be static, no?

  //std::cout << "Eigenvalues and vectors\n";
  //std::cout << eigen[0] << "\n";
  //std::cout << eigen[1] << "\n";
  //std::cout << "\n";

  std::vector<Cartesian> pca(5);
  pca[2].set_x(eigen[0](0,0));
  pca[2].set_y(eigen[0](1,0));
  pca[2].set_z(eigen[0](2,0));
  pca[1].set_x(eigen[0](0,1));
  pca[1].set_y(eigen[0](1,1));
  pca[1].set_z(eigen[0](2,1));
  pca[0].set_x(eigen[0](0,2));
  pca[0].set_y(eigen[0](1,2));
  pca[0].set_z(eigen[0](2,2));

  pca[3].set_x(Xmean);
  pca[3].set_y(Ymean);
  pca[3].set_z(Zmean);

  pca[4].set_x(eigen[1](2,0));
  pca[4].set_y(eigen[1](1,0));
  pca[4].set_z(eigen[1](0,0));

  //std::cout << pca[0] << "\n";
  //std::cout << pca[1] << "\n";
  //std::cout << pca[2] << "\n";
  //std::cout << "\n";
  
  return pca;
}

Cartesian Cartesian::CrossProduct(const Cartesian &v1, const Cartesian &v2){
  Cartesian v;
  v.x = v1.y * v2.z - v2.y * v1.z;
  v.y = v1.z * v2.x - v2.z * v1.x;
  v.z = v1.x * v2.y - v2.x * v1.y;
  return v;
}

Cartesian Cartesian::MidPoint(const std::vector<Cartesian> &v1){
  Cartesian v;
  if(v1.size()==0) return v;
  for(unsigned i=0;i<v1.size();i++)
    v += v1[i];
  v /= v1.size();
  return v;
}

Cartesian Cartesian::MidPoint(const Cartesian &v1, const Cartesian &v2){
  Cartesian v;
  v.x = v1.x+(v2.x-v1.x)/2;
  v.y = v1.y+(v2.y-v1.y)/2;
  v.z = v1.z+(v2.z-v1.z)/2;
  return v;
}


double Cartesian::DotProduct(const Cartesian &v1, const Cartesian &v2){
  double v;
  
  v  = v1.x * v2.x;
  v += v1.y * v2.y;
  v += v1.z * v2.z;
  return v;
}

void Cartesian::normalize(double radius){
   double d = sqrt(x*x+y*y+z*z);

   if (fabs(d)<1.0e-12) {
           std::cerr << "zero length vector in Cartesian::normalize" << *this << "\n";
      //abort();
      return;
   }

   double rdivd = radius/d;
   x *= rdivd;
   y *= rdivd;
   z *= rdivd;

}

Cartesian Cartesian::operator+(const Cartesian &in1) const {
   
   Cartesian out(x+in1.x, y+in1.y, z+in1.z);

   return out;

} 

Cartesian Cartesian::operator-() const {

   Cartesian out(-x, -y, -z);

   return out;
}

Cartesian Cartesian::operator-(const Cartesian &in1) const {

   Cartesian out(x-in1.x, y-in1.y, z-in1.z);

   return out;
}

void Cartesian::operator+=(const Cartesian &in) {

   x += in.x;
   y += in.y;
   z += in.z;

   // return *this;

}

void Cartesian::operator-=(const Cartesian &in) {

   x -= in.x;
   y -= in.y;
   z -= in.z;

}

double DotCart(const Cartesian &v1, const Cartesian &v2){
  double v;
  v  = v1.get_x() * v2.get_x();
  v += v1.get_y() * v2.get_y();
  v += v1.get_z() * v2.get_z();
  return v;
}

double Cartesian::length(void) const{
  return  sqrt(x*x+y*y+z*z);
}

double* Cartesian::to_dp() const{
  double *dp = new double[4];
  dp[0] = x;
  dp[1] = y;
  dp[2] = z;
  dp[3] = a;
  return dp;
}

void Cartesian::Scale(double a, double b, double c){
  x *= a;
  y *= b;
  z *= c;
}

void Cartesian::setMultiplyAndAdd(double m1, const Cartesian &c1, double m2, const Cartesian &c2){
  x = m1 * c1.x + m2 * c2.x;
  y = m1 * c1.y + m2 * c2.y;
  z = m1 * c1.z + m2 * c2.z;
}

bool Cartesian::CheckDistanceRangeMin(const Cartesian &v1, const Cartesian &v2, double minv){
  double d2 = (v1.get_x()-v2.get_x())*(v1.get_x()-v2.get_x()) + (v1.get_y()-v2.get_y())*(v1.get_y()-v2.get_y()) + (v1.get_z()-v2.get_z())*(v1.get_z()-v2.get_z());
  double minv2 = minv*minv;
  if(d2>minv2) return true;
  return false;
}

bool Cartesian::CheckDistanceRange(const Cartesian &v1, const Cartesian &v2, double minv, double maxv){
  double d2 = (v1.get_x()-v2.get_x())*(v1.get_x()-v2.get_x()) + (v1.get_y()-v2.get_y())*(v1.get_y()-v2.get_y()) + (v1.get_z()-v2.get_z())*(v1.get_z()-v2.get_z());
  double minv2 = minv*minv;
  double maxv2 = maxv*maxv;
  if((d2>minv2)&&(d2<maxv2)) return true;
  return false;
}
