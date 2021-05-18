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


#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <algorithm>
#include "matrix.h"
#include "cartesian.h"
#include "quat.h"

Cartesian RotMxV(double *objrotmat, const Cartesian &prim){
  double result[4];
  double input[4] = {prim.get_x(), prim.get_y(), prim.get_z(), prim.get_a()};

  for(int i=0;i<4;i++){
    result[i] = 0.0;
    for(int j=0;j<4;j++){
      result[i] += input[j]*objrotmat[i*4+j];
    }
  }
  
  return Cartesian(result);
}

matrix matrix::BlockMatrix(const std::vector<std::vector<matrix> > &blocks){

  if(blocks.size()==0||blocks[0].size()==0)
    return matrix();

  unsigned total_cols = 0;
  unsigned total_rows = 0;

  for(unsigned ii=0;ii<blocks.size();ii++)
    total_rows += blocks[ii][0].get_rows();

  for(unsigned jj=0;jj<blocks[0].size();jj++)
    total_cols += blocks[0][jj].get_columns();

  matrix a(total_rows,total_cols);

  unsigned ioff = 0;
  for(unsigned ii=0;ii<blocks.size();ii++){
    unsigned joff = 0;
    for(unsigned jj=0;jj<blocks[0].size();jj++){
      for(unsigned i=0;i<blocks[ii][jj].get_rows();i++){
	for(unsigned j=0;j<blocks[ii][jj].get_columns();j++){
	  a(ioff+i,joff+j) = blocks[ii][jj](i,j);
	}
      }
      joff += blocks[ii][jj].get_columns();
    }
    ioff += blocks[ii][0].get_rows();
  }

  return a;

}

void matrix::SwitchRows(const unsigned int &i1, const unsigned int &i2){

  std::vector<double> row1 = mat[i1];
  std::vector<double> row2 = mat[i2];

  mat[i1] = row2;
  mat[i2] = row1;

}

matrix matrix::TriangularMatrix() const {
  /* Use Gaussian elimination to create upper triangular matrix */
  if(get_rows()<1||get_columns()<1||get_rows()!=get_columns()){
    std::cout << "Cannot calculate triangular matrix of non-square matrix\n";
    std::cout << get_rows() << "\n";
    std::cout << get_columns()<< "\n";
    return *this;
  }

  matrix b = *this;

  unsigned int i=0;
  unsigned int j=0;

  int iswitch = 1;
  double mult = 1.0;

  while (i<b.get_rows()&&j<get_columns()) { 
    // Find pivot in column j, starting in row i:
    double max_val = b(i,j);
    unsigned int max_ind = i;
    for(unsigned int k=i+1;k<b.get_rows();k++){
      if(fabs(b(k,j))>fabs(max_val)){
        max_val = b(k,j);
        max_ind = k;
      }
    }
    mult *= max_val;

    if(fabs(max_val)>1e-7) {
      if(i!=max_ind){
        b.SwitchRows(i,max_ind);
        iswitch *= -1;
      }
      for(unsigned int k=0;k<b.get_columns();k++) 
        b(i,k) /= max_val;
      for(unsigned int u = i+1;u<b.get_rows();u++){
	double buj = b(u,j);
        for(unsigned int k=0;k<b.get_columns();k++){
          b(u,k) += -buj * b(i,k);
	}
      }
      i++;
    }
    j++;
  }

  for(unsigned int k=0;k<b.get_columns();k++) 
    b(0,k) *= mult *iswitch;

  return b;

}

matrix matrix::SolveLinearEquations(const matrix &a, const matrix &b){
  std::vector<int> perm;
  int parity;

  matrix lu = matrix::LUDecomposition(a,perm,parity);
  return matrix::LUSubstitution(lu,b,perm);

}

matrix matrix::LUSubstitution(const matrix &a, const matrix &b, const std::vector<int> &perm){
  /* 
     Performs substitution to solve for x in Ax=b, where A is a compund
     LU decomposed matrix.
  */

  matrix y=b;

  for(unsigned i=0;i<a.get_rows();i++){
    double sum = y(i,0);
    for(unsigned int j=0;j<i;j++)
      sum -= a(i,j) * y(j,0);
    y(i,0) = sum;
  }

  for(int i=a.get_rows()-1;i>=0;i--){
    double sum = y(i,0);
    for(unsigned j=i+1;j<a.get_columns();j++)
      sum -= a(i,j) * y(j,0);
    y(i,0) = sum/a(i,i);
  }

  matrix xx = y;

  for(unsigned i=0;i<xx.get_rows();i++){
    y(perm[i],0) = xx(i,0);
  }

  return y;
}

matrix matrix::LUMult(const matrix &a, const matrix &b, const std::vector<int> &perm){

  /* 
     Performs matrix * vector multiplication A*b, where A is a compound
     LU decomposed matrix. Illustartive example.
  */

  matrix bb = a.GetLowerTriangle(1) * a.GetUpperTriangle() * b;
  matrix c(b.get_rows(),1);
  for(unsigned i=0;i<b.get_rows();i++){
    c(perm[i],0) = bb(i,0);
  }

  return c;

}

matrix matrix::GetUpperTriangle() const {

  matrix c(get_rows(),get_columns());

  for(unsigned i=0;i<get_rows();i++)
    for(unsigned j=i;j<get_columns();j++)
      c(i,j) = mat[i][j];

  return c;

}

matrix matrix::GetLowerTriangle(int lu) const {

  matrix c(get_rows(),get_columns());

  for(unsigned i=0;i<get_rows();i++){
    if(lu) 
      c(i,i) = 1.0;
    else
      c(i,i) = mat[i][i];
    for(unsigned j=0;j<i;j++)
      c(i,j) = mat[i][j];
  }

  return c;

}

matrix matrix::LUDecomposition(const matrix &a, std::vector<int> &perm, int &parity){


  matrix lu = a;
  return lu;
    
}

matrix matrix::TriangularSolveForward(const matrix &a, const matrix &b){

  /*
    Solve a set of equations Ax = b where x and b are vectors and
     A is a lower triangular matrix.
  */
  matrix x(b.get_rows(),1);

  for(unsigned i=0;i<b.get_rows();i++){
    x(i,0) = -b(i,0);
    for(unsigned j=0;j<i;j++){
      x(i,0) += a(i,j) * x(j,0);
    }
    x(i,0) /= -a(i,i);
  }
  return x;

}

matrix matrix::TriangularSolveBack(const matrix &a, const matrix &b){

  /*
    Solve a set of equations Ax = b where x and b are vectors and
     A is an upper triangular matrix.
  */
  matrix x(b.get_rows(),1);

  for(int i=b.get_rows()-1;i>=0;i--){
    x(i,0) = -b(i,0);
    for(unsigned j=i+1;j<b.get_rows();j++){
      x(i,0) += a(i,j) * x(j,0);
    }
    x(i,0) /= -a(i,i);
  }
  return x;

}

matrix matrix::Inverse() const {

  /*
  std::vector<int> perm;
  int parity;
  matrix col(get_rows(),1);
  matrix inv(get_rows(),get_columns());

  matrix clu = LUDecomposition(*this,perm,parity);

  for(unsigned i=0;i<get_rows();i++){
    for(unsigned j=0;j<get_columns();j++)
      col(j,0) = 0.0;
    col(i,0) = 1.0;
    col = LUSubstitution(clu,col,perm);
    for(unsigned j=0;j<get_columns();j++)
      inv(j,i) = col(j,0);
  }

  return inv;
  */

  /* Use Gaussian elimination to create upper triangular matrix */
  if(get_rows()<1||get_columns()<1||get_rows()!=get_columns()){
    std::cout << "Cannot calculate inverse of non-square matrix\n";
    std::cout << get_rows() << "\n";
    std::cout << get_columns()<< "\n";
    return *this;
  }

  matrix b(get_rows(),get_columns()*2);

  for(unsigned int i=0;i<b.get_rows();i++){
    for(unsigned int j=0;j<b.get_columns();j++){
      if(j<get_columns()){
	b(i,j) = mat[i][j];
      }else{
	if(j==i+get_columns())
	  b(i,j) = 1.0;
	else
	  b(i,j) = 0.0;
      }
    }
  }

  unsigned int i=0;
  unsigned int j=0;

  int iswitch = 1;
  double mult = 1.0;

  while (i<b.get_rows()&&j<get_columns()) { 
    // Find pivot in column j, starting in row i:
    double max_val = fabs(b(i,j));
    unsigned int max_ind = i;
    for(unsigned int k=i+1;k<b.get_rows();k++){
      if(fabs(b(k,j))>fabs(max_val)){
        max_val = b(k,j);
        max_ind = k;
      }
    }
    mult *= max_val;

    if(fabs(max_val)>1e-7) {
      if(i!=max_ind){
        b.SwitchRows(i,max_ind);
        iswitch *= -1;
      }
      for(unsigned int k=0;k<b.get_columns();k++) 
        b(i,k) /= max_val;
      for(unsigned int u = i+1;u<b.get_rows();u++){
	double buj = b(u,j);
        for(unsigned int k=0;k<b.get_columns();k++){
          b(u,k) += -buj * b(i,k);
	}
      }
      i++;
    }
    j++;
  }

  //std::cout << "Upper triangle matrix:\n" << b << "\n"; std::cout.flush();

  // Now reduce to left hand side to unit matrix

  for(unsigned ii=1;ii<get_rows();ii++){
    for(unsigned i=0;i<get_rows()-ii;i++){
      //std::cout << "i, ii: " << i << " " << ii << "\n"; std::cout.flush();
      double val = b(i,i+ii)/b(i+ii,i+ii); // if b(i+1,i+1) matrix is singular
      for(unsigned int k=0;k<b.get_columns();k++){
	b(i,k) -= val * b(i+ii,k);
      }
      //std::cout << "New matrix:\n" << b << "\n"; std::cout.flush();
    }
    //std::cout << "ii loop matrix:\n" << b << "\n"; std::cout.flush();
  }
  //std::cout << "Unit matrix:\n" << b << "\n"; std::cout.flush();

  //for(unsigned int k=0;k<b.get_columns();k++) 
  //  b(0,k) *= mult *iswitch;

  matrix c(get_rows(),get_columns());

  for(unsigned int i=0;i<c.get_rows();i++){
    for(unsigned int j=0;j<c.get_columns();j++){
      c(i,j) = b(i,j+get_columns());
    }
  }

  //std::cout << "c * *this" << "\n";  
  //std::cout << c * *this << "\n";
  //std::cout << "*this * c" << "\n";  
  //std::cout << *this * c << "\n";

  return c;

}

double matrix::Determinant() const {
  if(get_rows()<1||get_columns()<1||get_rows()!=get_columns()){
    std::cout << "Cannot calculate determinant of non-square matrix\n";
    return 0.0;
  }
  matrix b = TriangularMatrix();

  double det = b(0,0);
  for(unsigned int i = 1; i<b.get_columns(); i++)
    det *= b(i,i);

  return det;
}

double matrix::Minor(const matrix &a, const unsigned int &row, const unsigned int &col){
  return MinorMatrix(a,row,col).Determinant();
}

matrix matrix::MinorMatrix(const matrix &a, const unsigned int &row, const unsigned int &col){
  if(a.get_rows()<1||a.get_columns()<1||row>a.get_rows()-1||col>a.get_columns()-1)
    return matrix(1,1);

  matrix minormat(a.get_rows()-1,a.get_columns()-1);
  int ii,jj;

  for(unsigned int i=0;i<a.get_rows();i++){
    if(i!=row) {
      for(unsigned int j=0;j<a.get_columns();j++){
         if(j!=col){
           if(i<row)
             ii = i;
           else
             ii = i-1;
           if(j<col)
             jj = j;
           else
             jj = j-1;
           minormat(ii,jj) = a(i,j);
        }
      }
    }
  }
  return minormat;
}

double* matrix::to_dp(){

   double *result = 0;

   if (! mat.empty()) {
      result = new double[get_columns()*get_rows()];
      int l=0;

      for(unsigned int i=0;i<get_rows();i++){
         for(unsigned int j=0;j<get_columns();j++){
            result[l++] = mat[i][j];
         }
      }
   }
   return result;
}

unsigned int matrix::get_rows() const{
  return mat.size();
}

unsigned int matrix::get_columns() const{
  if(mat.size()>0)
    return mat[0].size();
  return 0;
}

// CONSTRUCTORS

void matrix::Zero(){
  for(unsigned int i=0;i<get_rows();i++){
    for(unsigned int j=0;j<get_columns();j++){
      mat[i][j]=0.0;
    }
  }
}

matrix& matrix::operator=(const matrix &a){
  //std::cout << "matrix::operator=(const matrix &a)\n"; std::cout.flush();
  mat = a.mat;
  return *this;
}

matrix::matrix(const Quat &q){
  //std::cout << "matrix::matrix(const Quat &q)\n"; std::cout.flush();
  *this = q.getMatrix();
}

matrix::matrix(const matrix &a){
  //std::cout << "matrix::matrix(const matrix &a)\n"; std::cout.flush();
  mat = a.mat;
}

matrix::~matrix(){
  //std::cout << "Deleting matrix " << get_rows() << ", " << get_columns() << "\n"; std::cout.flush();
  for(unsigned int i=0;i<get_rows();i++)
    mat[i].clear();
  mat.clear();
}

matrix::matrix(unsigned int x,unsigned int y){
  //std::cout << "matrix::matrix(unsigned int x,unsigned int y)\n"; std::cout.flush();

  mat = std::vector<std::vector<double> >(x);
  for(unsigned int i=0;i<x;i++)
    mat[i] = std::vector<double>(y);
  Zero();
}

matrix::matrix(unsigned int x,unsigned int y, double *array){
  //std::cout << "matrix::matrix(unsigned int x,unsigned int y, double *array)\n"; std::cout.flush();

  mat = std::vector<std::vector<double> >(x);
  for(unsigned int i=0;i<x;i++)
    mat[i] = std::vector<double>(y);

  for(unsigned int i=0;i<get_rows();i++){
    for(unsigned int j=0;j<get_columns();j++){
      mat[i][j]=*array++;
    } 
  }
}
matrix::matrix(unsigned int x,unsigned int y, const std::vector<double> &mat_in){

  //std::cout << "matrix::matrix(unsigned int x,unsigned int y, const std::vector<double> &mat_in)\n"; std::cout.flush();
  mat = std::vector<std::vector<double> >(x);
  for(unsigned int i=0;i<x;i++)
    mat[i] = std::vector<double>(y);
  std::vector<double>::const_iterator k = mat_in.begin();

  for(unsigned int i=0;i<get_rows();i++){
    for(unsigned int j=0;j<get_columns();j++){
      mat[i][j]=*k++;
    } 
  }
}

matrix::matrix(unsigned int x,double(*fun)(int,int,int)){
  //std::cout << "matrix::matrix(unsigned int x,double(*fun)(int,int,int))\n"; std::cout.flush();
  
  mat = std::vector<std::vector<double> >(x);
  for(unsigned int i=0;i<x;i++)
    mat[i] = std::vector<double>(x);

  for(unsigned int i=0;i<get_rows();i++){
    for(unsigned int j=0;j<get_columns();j++){
      mat[i][j]=fun(i,j,get_columns());
    } 
  }
}

matrix matrix::operator+ (const matrix &b) const {

  if((get_rows() != b.get_columns())||(get_rows()!=b.get_columns())){
    std::cerr << "Matrices are not same size in addition!\n";
    return matrix(0,1);
  }
  
  matrix c(get_rows(),get_columns());

  for(unsigned int i=0;i<get_rows();i++){
    for(unsigned int j=0;j<get_columns();j++){
     c(i,j) =  mat[i][j] + b(i,j);
    }
  }
  return c;
}

matrix matrix::operator- (const matrix &b) const {
    return *this+(-b);
}

matrix matrix::operator- () const {
  
  matrix c(get_rows(),get_columns());

  for(unsigned int i=0;i<get_rows();i++){
    for(unsigned int j=0;j<get_columns();j++){
      c.mat[i][j] = -mat[i][j];
    }
  }
  return c;
}

void matrix::operator+= (const matrix &a){
  
  *this = *this + a;
}

void matrix::operator-= (const matrix &a){
  
  *this = *this - a;
}

void matrix::operator*= (const matrix &a){
  
  *this = *this * a;
}

void matrix::operator/= (const matrix &a){
  
  *this = *this / a;
}

matrix& matrix::operator*= (double val){
  
  *this = *this * val;
  return *this;
}

matrix& matrix::operator/= (double val){
  
  *this = *this / val;
  return *this;
}

matrix operator/ (double val, const matrix &a){
  return val * a.Inverse();
}

matrix operator* (double val, const matrix &a){

  matrix c(a.get_rows(),a.get_columns());

  for(unsigned int i=0;i<a.get_rows();i++){
    for(unsigned int j=0;j<a.get_columns();j++){
      c.mat[i][j] = val*a.mat[i][j];
    }
  }
  return c;

}

matrix operator/ (const matrix &a, double val){
  return (1.0/val)*a;
}

matrix operator* (const matrix &a, double val){
  return val*a;
}

matrix matrix::operator/ (const matrix &b) const {
  return *this * b.Inverse();
}

matrix matrix::operator* (const matrix &b) const {

  if(get_columns()!=b.get_rows()){
     std::cerr << "Matrices cannot be multiplied\n";
     return matrix(0,1);
  }

  matrix c(get_rows(),b.get_columns());

  for(unsigned int k=0;k<get_rows();k++){
      for(unsigned int s=0;s<b.get_columns();s++){
	  for(unsigned int i=0;i<get_columns();i++){
	      c.mat[k][s] += mat[k][i] * b.mat[i][s];
	  }
      }
  }
  return c;
}

double matrix::Trace() const{

  double sum=0.0;
  
  if(get_columns()!=get_rows()){
     std::cerr << "Cannot compute trace of non-square matrix\n";
     return 0;
  }
  
  for(unsigned int i=0;i<get_columns();i++){
      sum += mat[i][i];
  }
  return sum;
}


matrix matrix::Transpose() const{

    matrix c(get_columns(),get_rows());

    for(unsigned int i=0;i<get_rows();i++){
      for(unsigned int j=0;j<get_columns();j++){
	c.mat[j][i] = mat[i][j];
      }
    }
    return c;
}


matrix pow(const matrix &a, double p){

    matrix c=a;
    std::vector<matrix> eigen = c.Eigen();
    matrix evals = eigen[1];
    matrix d=eigen[0];

    for(unsigned int i=0;i<c.get_rows();i++)
      for(unsigned int j=0;j<c.get_columns();j++)
	c.mat[i][j] = i==j?(pow(evals(i,0),p)):0.0;

    return d*c*d.Transpose();
}

// Apply arbitrary function to matrix. *Any* function which takes one argument, sin, cos, exp, etc.

matrix matrix::fun(const matrix &a, const std::string &fun) const {

  if(fun=="acos")
    return a.fun(a,(double (*)(double))acos);
  if(fun=="atan")
    return a.fun(a,(double (*)(double))atan);
  if(fun=="asin")
    return a.fun(a,(double (*)(double))asin);
  if(fun=="cos")
    return a.fun(a,(double (*)(double))cos);
  if(fun=="tan")
    return a.fun(a,(double (*)(double))tan);
  if(fun=="sin")
    return a.fun(a,(double (*)(double))sin);
#ifndef WIN32
  if(fun=="acosh")
    return a.fun(a,(double (*)(double))acosh);
  if(fun=="atanh")
    return a.fun(a,(double (*)(double))atanh);
  if(fun=="asinh")
    return a.fun(a,(double (*)(double))asinh);
#endif
  if(fun=="cosh")
    return a.fun(a,(double (*)(double))cosh);
  if(fun=="tanh")
    return a.fun(a,(double (*)(double))tanh);
  if(fun=="sinh")
    return a.fun(a,(double (*)(double))sinh);
  if(fun=="exp")
    return a.fun(a,(double (*)(double))exp);
  if(fun=="log")
    return a.fun(a,(double (*)(double))log);
  if(fun=="sqrt")
    return a.fun(a,(double (*)(double))sqrt);

  return a;
}

matrix matrix::fun(const matrix &a, double(*f)(double)) const {

    matrix c=a;
    std::vector<matrix> eigen = c.Eigen();
    matrix evals = eigen[1];
    matrix d=eigen[0];

    for(unsigned int i=0;i<d.get_rows();i++)
      for(unsigned int j=0;j<d.get_columns();j++)
	c.mat[i][j] = i==j?(f(evals(i,0))):0.0;

    return d*c*d.Transpose();
}

// Delta Function

double kdelta(int i,int j,int k){
  return i==j?1.0:0.0;
}

double imag2(int i,int j,int k){
  
  if(i==0&&j==1) {
    return 1;
  }else if(i==1&&j==0) {
    return -1;
  }
  return 0;
}

#ifndef M_PI
#define M_PI 3.14159265358979323844
#endif

#define SWAP(g,h) {y=(g);(g)=(h);(h)=y;}

int imax(int i,int j){
  return i>j?i:j;
}

double dsign(double z, double p){
  return z*(fabs(p)/p);
}

#define TOL 1.0e-16

std::vector<matrix> matrix::SortEigenvalues(const std::vector<matrix>&eigen) const{

  matrix vectors = eigen[0];
  matrix evals = eigen[1];
  matrix b = eigen[2];

  std::vector<std::pair<double,int> > sort_map;

  for(unsigned i=0;i<evals.get_rows();i++){
    sort_map.push_back(std::pair<double,int>(evals(i,0),i));
  }
  
  std::sort(sort_map.begin(),sort_map.end());

  for(unsigned i=0;i<vectors.get_rows();i++){
    evals(i,0) = eigen[1](sort_map[i].second,0);
    b(i,i)     = eigen[2](sort_map[i].second,sort_map[i].second);
    for(unsigned j=0;j<vectors.get_columns();j++){
      vectors(i,j) = eigen[0](i,sort_map[j].second);
    }
  }

  std::vector<matrix> sorted;
  sorted.push_back(vectors);
  sorted.push_back(evals);
  sorted.push_back(b);
  return sorted;

}

std::vector<matrix> matrix::Eigen() const{

  std::vector<matrix> eigenprob;

  // Whether to calc vectors or not, this will
  // be an option in class definition.

  unsigned int i, j, p=0, q=0, iter=0, k, s;
  double maxval=0.0, theta, bip, biq, cpp, cqq, cpq, cqp;
  matrix b=*this;
  matrix vectors(get_rows(),kdelta);
  matrix tmpvectors(get_rows(),get_columns());
  matrix r(get_rows(),kdelta);

  if(get_rows()!=get_columns()){
     std::cerr << "Cannot calculate eigenvalues of non-square matrix!\n";
     return eigenprob;
  }
  
  while(iter<get_columns()*get_columns()*100){
    for(i=1;i<b.get_rows();i++){
      for(j=0;j<i;j++){
	if(fabs(b(i,j)-b(j,i))>1e-12){
	  std::cerr << "Cannot calculate eigenvalues of non-symmetric matrix using Eigen!\n";
	  return eigenprob;
	}
	if(fabs(b(i,j))>fabs(maxval)) {
	  maxval=b(i,j);
	  p=i;
	  q=j;
	}
      }
    }
    
    if(fabs(maxval)<TOL) {
      matrix evals(get_columns(),1);
      for(i=0;i<get_columns();i++) evals(i,0) = b(i,i);
      eigenprob.push_back(vectors);
      eigenprob.push_back(evals);
      eigenprob.push_back(b);
      return eigenprob;
    }
    
    if((fabs(b(p,p)-b(q,q)))>TOL){
      theta = (atan(2*b(p,q)/(b(p,p)-b(q,q))))/2;
    } else {
      theta = b(p,q)/fabs(b(p,q)) * M_PI/4;
    }
    
    
    r(p,p) = cos(theta);
    r(p,q) = -sin(theta);
    r(q,p) = sin(theta);
    r(q,q) = cos(theta);
    
    tmpvectors = vectors;
    for(k=0;k<r.get_rows();k++){
      for(s=0;s<r.get_columns();s++){
	if((s==p)||(s==q)){
	  tmpvectors(k,s) = 0.0;
	  for(i=0;i<r.get_columns();i++){
	    tmpvectors(k,s) += vectors(k,i) * r(i,s);
	  }
	}
      }
    }
    vectors = tmpvectors;

    r(p,p) = 1.0;
    r(q,q) = 1.0;
    r(p,q) = 0.0;
    r(q,p) = 0.0;

    for(i=0;i<b.get_rows();i++){
      if((i != p) && (i != q)){
	bip = b(i,p)*cos(theta)+b(i,q)*sin(theta);
	biq =-b(i,p)*sin(theta)+b(i,q)*cos(theta);
	b(i,p) = b(p,i) = bip;
	b(i,q) = b(q,i) = biq;
      }
    }

    cpp =   b(p,p)*cos(theta)*cos(theta)
      + 2*b(p,q)*cos(theta)*sin(theta)
      +   b(q,q)*sin(theta)*sin(theta);

    cqq =   b(p,p)*sin(theta)*sin(theta)
      - 2*b(p,q)*cos(theta)*sin(theta)
      +   b(q,q)*cos(theta)*cos(theta);

    cpq = cqp = (b(q,q)-b(p,p))*cos(theta)*sin(theta)
      + b(p,q)*(cos(theta)*cos(theta) - sin(theta)*sin(theta)); 
    
    b(p,p) = cpp;
    b(q,q) = cqq;
    b(p,q) = cpq;
    b(q,p) = cqp;

    maxval=0.0;
    iter++;
  }
  std::cout << "Diagonalization was probably unsuccessful!\n";
  std::cout << "Diagonalization used " << iter << " iterations\n";

  return eigenprob;

}

matrix matrix::DirSum(const matrix &a, const matrix &b){

  unsigned int i,j;

  matrix c(a.get_rows()+b.get_rows(),a.get_columns()+b.get_columns());

  for(i=0;i<a.get_rows();i++){
    for(j=0;j<a.get_columns();j++){
      c(i,j)=a(i,j);
    }
  }
  for(i=0;i<b.get_columns();i++){
    for(j=0;j<b.get_rows();j++){
      c(i+a.get_rows(),j+a.get_columns())=b(i,j);
    }
  }
  
  return c;
}

/*
matrix matrix::DirProd(const matrix &a, const matrix &b){
    
    int i,j,k,l,ij,kl,aij,bkl;
    
    matrix c(a.get_rows()*b.get_rows(),a.get_columns()*b.get_columns());

    aij=0;
    for(i=0;i<a.get_rows();i++){
      for(j=0;j<a.get_columns();j++){
	ij = b.get_columns()*b.get_rows()*a.get_columns()*i + b.get_columns()*j;
	bkl=0;
	for(k=0;k<b.get_rows();k++){
	  for(l=0;l<b.get_columns();l++){
	    kl = a.get_columns()*b.get_columns()*k + l;

	    c.mat[ij+kl] = a.mat[aij]*b.mat[bkl];

	    bkl++;
	  }
	}
	aij++;
      }
    }
    
    return c;
}
*/
std::ostream& operator<<(std::ostream& c,matrix a){

  c.flags(std::ios::fixed|std::ios::right);
  c.precision(6);

  for(unsigned int i=0;i<a.get_rows();i++){
    for(unsigned int j=0;j<a.get_columns();j++){
      c << std::setw(10) << a(i,j) << " ";
    }
    c << "\n";
  }
  
  return c;
}
