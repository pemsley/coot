/*
     util/matrix.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2011 University of York

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


#include <iostream>
#include <iomanip>
#include <math.h>
#ifdef _DO_TIMINGS_
#include <time.h>
#endif
#include <vector>
#include <algorithm>
#include "matrix.h"
#include "cartesian.h"
#include "quat.h"

#define PIOVER4 M_PI/4 ;

/* 

 Need to implement:
   Templatized versions
     o float/double.
     o complex
*/

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

matrix matrix::SolveLinearEquations(const matrix &a, const matrix &b){
  if(a.get_rows()!=a.get_columns()||a.get_rows()!=b.get_rows()||b.get_columns()!=1){
    std::cerr << "Impossible linear equations\n";
    return a;
  }
  matrix gin = matrix(a.get_rows(),a.get_columns()+1);
  for(unsigned i=0;i<a.get_rows();i++){
    for(unsigned j=0;j<a.get_columns();j++){
      gin(i,j) = a(i,j);
    }
    gin(i,a.get_columns()) = b(i,0);
  }
  //std::cout << "gin:\n" << gin << "\n";
  matrix gauss = gin.GaussianElimination();
  //std::cout << "Gauss\n" << gauss << "\n";
  matrix gout = matrix(a.get_rows(),a.get_columns());
  matrix bout = matrix(a.get_rows(),1);
  for(unsigned i=0;i<a.get_rows();i++){
    for(unsigned j=0;j<a.get_columns();j++){
      gout(i,j) = gauss(i,j);
    }
    bout(i,0) = gauss(i,a.get_columns());
  }
  //std::cout << "gout:\n" << gout << "\n";
  //std::cout << "bout:\n" << bout << "\n";

  matrix x = TriangularSolveBack(gout,bout)[0];
  //std::cout << "x: " << x << "\n";
  return x;

}

matrix matrix::GaussianElimination(bool unitDiag) const {
  matrix g = *this;

  int i=0;
  int j=0;
  int k;

  int m = get_rows();
  int n = get_columns();

  int maxi;

  while(i<m&&j<n){
    maxi = i;
    for(k=i+1;k<m;k++){
      if(fabs(g(k,j)) > fabs(g(maxi,j)))
        maxi = k;
    }
    //std::cout<< "Pivot: " << maxi << "\n";
    if(fabs(g(maxi,j))>1e-12){
      //swap rows i and maxi, but do not change the value of i
      g.SwitchRows(i,maxi);
      //std::cout << "Swapped:\n" << g << "\n\n";

      //Now A[i,j] will contain the old value of A[maxi,j].
      //divide each entry in row i by A[i,j]
      double multt = g(i,j);
      //std:: cout << "mult t: " << multt << " " << i << " " << j << "\n";
      for(int l=0;l<n&&fabs(multt)>1e-12;l++){
         g(i,l) = g(i,l)/multt;
      }
      //std::cout << "Divided:\n" << g << "\n\n";
      //Now A[i,j] will have the value 1.
      //subtract A[u,j] * row i from row u
      for(int u=i+1;u<m;u++){
        double t = g(u,j);
        //std:: cout << "mult: " << t << "\n";
        for(int l=0;l<n;l++){
          //std::cout << "orig: " << g(u,l) << "\n";
          g(u,l) -= t*g(i,l);
          //std::cout << "new: " << g(u,l) << "\n";
        }
        //Now A[u,j] will be 0, since A[u,j] - A[i,j] * A[u,j] = A[u,j] - 1 * A[u,j] = 0.
      }
      //std::cout << "Minussed:\n" << g << "\n\n";
      if(!unitDiag){
        for(int l=0;l<n;l++){
          g(i,l) = g(i,l)*multt;
        }
      }
      i = i + 1;
    }
    j = j + 1;
  }

  return g;
}

matrix matrix::GaussJordanElimination() const {
  matrix g=GaussianElimination(true);
  //std::cout << "gout:\n" << g << "\n";
  for(unsigned i=0;i<get_rows()-1;i++){
    for(unsigned k=i+1;k<get_rows();k++){
      double mult = g(i,k);
      //std::cout << i << " " << k << " " << get_columns() << " " << mult << "\n";
      for(unsigned j=i+1;j<g.get_columns();j++){
        g(i,j) = g(i,j) - mult*g(k,j);
      }
      //std::cout << g << "\n";
    }
  }
  return g;
}

std::vector<matrix> matrix::TriangularSolveForward(const matrix &a, const matrix &b, unsigned nsolutions){

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
    if(fabs(a(i,i))>1e-12)
    x(i,0) /= -a(i,i);
  }
  std::vector<matrix> solutions;
  solutions.push_back(x);
  return solutions;

}

std::vector<matrix> matrix::TriangularSolveBack(const matrix &a, const matrix &b, unsigned nsolutions){

  /*
    Solve a set of equations Ax = b where x and b are vectors and
     A is an upper triangular matrix.
  */

  std::vector<matrix> solutions;
  matrix gotten(a.get_rows(),a.get_rows());

  matrix x(b.get_rows(),1);

  if(b.get_rows()==0||a.get_rows()==0||a.get_rows()!=b.get_rows()) return solutions;
  bool homo = true;
  for(unsigned i=0;i<b.get_rows();i++){
    if(fabs(b(i,0))>1e-6) homo = false;
  }

  //Should now loop over nsolutions if homogeneous;

  for(unsigned isol=0;isol<nsolutions;isol++){
    x.Zero();
  bool doneSomething = false;
  for(int i=b.get_rows()-1;i>=0;i--){
    // This attempts to find non-trivial solutions, this may not be very general.
    // Needs testing!.
    if(homo&&!doneSomething&&(i!=(int)(b.get_rows()-1))&&(fabs(a(i,i))>1e-6)&&(fabs(a(i+1,i+1))<1e-6)){
            for(unsigned j=i+1;j<a.get_columns();j++){
               if(fabs(a(i,j))>1e-6&&!gotten(i,j)){
                 if(fabs(x(j,0))<1e-6){
                   //std::cout << "Setting " << j << " for row " << i << "\n";
                   x(j,0) = 1;
                   doneSomething = true;
                   gotten(i,j) = 1;
                   break;
                 }
               }
            }
    }
    if(homo&&!doneSomething&&(i!=(int)(b.get_rows()-1))&&!a.GetRow(i).isNull()){
        //std::cout << "Do something for row " << i << ", solution: " << isol << "\n";
            int first = -1;
            int second = -1;
            for(unsigned j=i+1;j<a.get_columns()&&second==-1;j++){

               if(fabs(a(i,j))>1e-6&&!gotten(i,j)){
                 if(first==-1) first = j;
                 else second = j;
               }

            }
            if(first!=-1&&second!=-1){
              //std::cout << first << " " << second << "\n";
              x(second,0) = 1;
              x(first,0) = -a(i,second)/a(i,first);
              doneSomething = true;
              gotten(i,second) = 1;
            }
    }
    x(i,0) = -b(i,0);
    for(unsigned j=i+1;j<b.get_rows();j++){
      x(i,0) += a(i,j) * x(j,0);
    }
    if(fabs(a(i,i))>1e-12)
    x(i,0) /= -a(i,i);
  }
  //std::cout << x << "\n";
  //std::cout << (a*x).isNull() << "\n";
  if(homo&&!doneSomething&&x.isNull()){
    //std::cout << "Hmm, solution " << isol << " is Null\n";
    // Try some unit vectors.
    for(unsigned ii=0;ii<a.get_rows();ii++){
      bool doneThisLine = false;
      for(unsigned jj=ii;jj<a.get_rows();jj++){
         if(gotten(ii,jj)>0){
           doneThisLine = true;
           break;
         }
      }
      if(!doneThisLine) {
        //std::cout << "Row " << ii << " available\n";
        x(ii,0) = 1;
        if((a*x).isNull()){
           gotten(ii,ii) = 1;
           //std::cout << "Row " << ii << " works\n";
           break;
        } else x.Zero();
      }
    }
    
  }
  solutions.push_back(x);
  }

  return solutions;

}

double matrix::FrobeniusNorm() const {
    std::vector<matrix> svd = SVD();
    matrix s = svd[1];
    unsigned mind = s.get_rows();
    if(s.get_columns()<mind) mind = s.get_columns();
    double norm = 0.0;
    for(unsigned i=0;i<mind;i++)
      norm += s(i,i);
  
    return sqrt(norm);
}

matrix matrix::Inverse() const {

  if(get_rows()<1||get_columns()<1||get_rows()!=get_columns()){
    std::vector<matrix> svd = SVD();
    matrix ut = svd[2];
    matrix v = svd[0];
    matrix s = svd[1];
    
    unsigned mind = s.get_rows();
    if(s.get_columns()<mind) mind = s.get_columns();
    for(unsigned i=0;i<mind;i++){
      if(fabs(s(i,i))>1e-12) s(i,i) = 1./s(i,i);
    }

    matrix mpinv = (v*s*ut).Transpose();

    return mpinv;
  }
  matrix gin = matrix(get_rows(),get_columns()*2);
  for(unsigned i=0;i<get_rows();i++){
    for(unsigned j=0;j<get_columns();j++){
      gin(i,j) = (*this)(i,j);
    }
    gin(i,i+get_columns()) = 1.0;
  }
  //std::cout << "gin:\n" << gin << "\n";
  matrix g=gin.GaussJordanElimination();
  matrix inv(get_columns(),get_rows());
  for(unsigned i=0;i<get_rows();i++){
    for(unsigned j=0;j<get_columns();j++){
      inv(i,j) = g(i,j+get_columns());
    }
  }
  //std::cout << g << "\n";
  return inv;
}

double matrix::Determinant() const {
  if(get_rows()<1||get_columns()<1||get_rows()!=get_columns()){
    std::cout << "Cannot calculate determinant of non-square matrix\n";
    return 0.0;
  }
  matrix gauss = GaussianElimination();
  //std::cout << "Gauss\n" << gauss << "\n";
  double det=1.0;
  for(unsigned i=0;i<get_rows();i++)
    det *= gauss(i,i);

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

  unsigned rows = a.get_rows();
  unsigned cols = a.get_columns();
  for(unsigned int i=0;i<rows;i++){
    if(i!=row) {
      for(unsigned int j=0;j<cols;j++){
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

  double *result = new double[get_columns()*get_rows()];
  int l=0;

  for(unsigned int i=0;i<get_rows();i++){
    for(unsigned int j=0;j<get_columns();j++){
      result[l++] = mat[i][j];
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
  unsigned cols = get_columns();
  for(unsigned int i=0;i<get_rows();i++){
    for(unsigned int j=0;j<cols;j++){
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
    mat[i].resize(y);
  // this is not necessary is it?
  //Zero();
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

  if((get_rows() != b.get_rows())||(get_columns()!=b.get_columns())){
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

void matrix::MultiplyAndSetLeft (const matrix &b) {
}
void matrix::MultiplyAndSetRight (const matrix &b) {
  unsigned int myRows = get_rows();
  unsigned int myCols = get_columns();
  unsigned int bCols = b.get_columns();

  std::vector<std::vector<double> > bmat = b.Transpose().mat;
  std::vector<double> matknew(myRows);
  for(unsigned int k=0;k<myRows;k++){
     std::vector<double> matk = mat[k];
     for(unsigned int s=0;s<bCols;s++){
        double sum = 0;
        std::vector<double> bmats = bmat[s];
        for(unsigned int i=0;i<myCols;i++){
           sum += matk[i] * bmats[i];
        }
        matknew[s] = sum;
    }
    mat[k] = matknew;
  }
}

matrix matrix::operator* (const matrix &b) const {

  if(get_columns()!=b.get_rows()){
     std::cerr << "Matrices cannot be multiplied\n";
     return matrix(0,1);
  }

  matrix c(get_rows(),b.get_columns());

  unsigned int myRows = get_rows();
  unsigned int myCols = get_columns();
  unsigned int bCols = b.get_columns();

  std::vector<std::vector<double> > bmat = b.Transpose().mat;

  for(unsigned int k=0;k<myRows;k++){
     std::vector<double> matk = mat[k];
     for(unsigned int s=0;s<bCols;s++){
        double sum = 0;
        std::vector<double> bmats = bmat[s];
        for(unsigned int i=0;i<myCols;i++){
           sum += matk[i] * bmats[i];
        }
        c.mat[k][s] = sum;
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
#ifndef _WIN32
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

#define TOL 1.0e-12

std::vector<double> matrix::matrixVector() {
  std::vector<double> matout;
  for(unsigned i=0;i<mat.size();i++){
    for(unsigned j=0;j<mat[i].size();j++){
      matout.push_back(mat[j][i]);
    }
  }
  return matout;
}
    


std::vector<matrix> matrix::SortEigenvalues(const std::vector<matrix>&eigen,bool descending){

  std::vector<matrix> sorted;

  if(eigen.size()==0){
    std::cerr << "Cannot sort eigenvalues!!\n";
    std::cerr << "Null eigenproblem\n";
    return sorted;
  }

  matrix vectors = eigen[0];
  matrix evals = eigen[1];
  matrix b = eigen[2];

  std::vector<std::pair<double,int> > sort_map;

  for(unsigned i=0;i<evals.get_rows();i++){
    sort_map.push_back(std::pair<double,int>(evals(i,0),i));
  }
  
  std::sort(sort_map.begin(),sort_map.end());
  if(descending) std::reverse(sort_map.begin(),sort_map.end());

  for(unsigned i=0;i<vectors.get_rows();i++){
    evals(i,0) = eigen[1](sort_map[i].second,0);
    b(i,i)     = eigen[2](sort_map[i].second,sort_map[i].second);
    for(unsigned j=0;j<vectors.get_columns();j++){
      vectors(i,j) = eigen[0](i,sort_map[j].second);
    }
  }

  sorted.push_back(vectors);
  sorted.push_back(evals);
  sorted.push_back(b);
  return sorted;

}

std::vector<matrix> GetTheBlocks(const matrix &a){
	std::vector<matrix> qrorig = a.QRDecomposition();
	std::vector<matrix> blocks;
	matrix q = qrorig[0];
	std::vector<std::vector<int> > pairs(q.get_rows());
	std::vector<unsigned> gone;
	for(unsigned i=0;i<q.get_rows();i++){
		if(std::find(gone.begin(),gone.end(),i)==gone.end()){
			for(unsigned j=i+1;j<q.get_rows();j++){
				if(fabs(q(i,j))>1e-7){
					gone.push_back(j);
					pairs[i].push_back(j);
				}
			}
			matrix b(pairs[i].size()+1,pairs[i].size()+1);
			b(0,0) = a(i,i);
			for(unsigned k=0;k<pairs[i].size();k++){
				b(k+1,k+1) = a(pairs[i][k],pairs[i][k]);
				b(k,k+1) = a(i,pairs[i][k]);
				b(k+1,k) = a(pairs[i][k],i);
			}
			blocks.push_back(b);
		}
	}
	return blocks;
}

std::vector<double>  SearchForDegenerateEigenValues(const matrix &a){
	std::vector<double> dgenEvals;
	for(unsigned i=0;i<a.get_rows();i++){
		for(unsigned j=i+1;j<a.get_rows();j++){
			matrix m(2,2);
			m(0,0) = a(i,i);
			m(1,1) = a(j,j);
			m(0,1) = a(i,j);
			m(1,0) = a(j,i);
			//std::cout << m << "\n";
			double A = 1.0;
			double B = -m.Trace();
			double C = m.Determinant();
			//std::cout << B*B << " " << 4*A*C << "\n";
			if(B*B-4*A*C<-1e-7){
				//std::cout << "BAD no real sub-matrix roots\n";
			} else {
				//std::cout << "OK real sub-matrix roots\n";
				if(fabs(B)<1e-7){
					//std::cout << "Degenerate pairs\n";
					//std::cout << i << "  " << j << " " << sqrt(B*B-4*A*C)/2 << " " << -sqrt(B*B-4*A*C)/2 << "\n";
					dgenEvals.push_back(sqrt(-C));
					dgenEvals.push_back(-sqrt(-C));
				}else if(fabs(B*B-4*A*C)<1e-7){
					//std::cout << "Degenerate pairs\n";
					//std::cout << i << "  " << j << " " << -B/2 << "\n";
					dgenEvals.push_back(-B/2);
					dgenEvals.push_back(-B/2);
				}
			}
		}
	}
	return dgenEvals;
}

matrix matrix::QREigen() const {

	matrix qall = matrix(get_rows(),get_columns());
	
	matrix anew = *this;
	
	for(unsigned i=0;i<get_columns();i++)
		qall(i,i) = 1.0;

	/* This should be convergant iterations! */
	bool converged = false;
	//for(unsigned i=0;i<10;i++){
	int iter = 0;
	double diff = 0;
	int i0 = -1;
	int j0 = -1;
	//for(unsigned i=0;i<50;i++){
	int max_iter = 500;
	while(!converged&&iter<max_iter){
		std::vector<matrix> qr = anew.QRDecomposition();
		matrix q = qr[0];
		matrix r = qr[1];
		matrix anewTmp = r*q;
		converged = true;
		for(unsigned i=0;i<get_rows()&&converged;i++){
			for(unsigned j=0;j<get_columns()&&converged;j++){
				if(fabs(fabs(anewTmp(i,j))-fabs(anew(i,j)))>1e-7){
					converged = false;
					diff = fabs(anewTmp(i,j))-fabs(anew(i,j));
					i0 = i;
					j0 = j;
					break;
				}
			}
		}
		iter++;
		//std::cout << iter << " " << diff << "\n"; std::cout.flush();
		/*
		if(iter==max_iter){
			std::cerr << "QR iterations failed\n";
			std::cerr << "Old RQ\n" << anew << "\n";
			std::cerr << "New RQ\n" << anewTmp << "\n";
			std::cerr << "diff " << diff << ": " << i0 << " " << j0 << "\n";
		}
		*/
		anew = anewTmp;
		qall = qall *q;
	}

	matrix evalsMat = anew;

	if(!converged){
		// Search for blocks, in fact, it is probably a good idea to do this anyway...
		std::vector<matrix> blocks = GetTheBlocks(*this);
		// Then for each block:
		int iblock = 0;
		for(unsigned k=0;k<blocks.size();k++){
			//std::cout << "Block " << k << "\n" << blocks[k] << "\n";
			std::vector<double> dgenEvals = SearchForDegenerateEigenValues(blocks[k]);
			if(dgenEvals.size()>0){
				//std::cout << "Degen:\n";
				for(unsigned l=0;l<dgenEvals.size();l++){
					//std::cout << dgenEvals[l] << "\n";
					evalsMat(l+iblock,l+iblock) = dgenEvals[l];
				}
			}else if(blocks[k].get_rows()==1){
					evalsMat(iblock,iblock) = blocks[k](0,0);
			}else if(blocks.size()>1){
				// This had better converge!
				matrix nonDegen = blocks[k].QREigen();
				//std::cout << "Non-degen\n" << nonDegen << "\n";
				for(unsigned l=0;l<nonDegen.get_rows();l++){
					evalsMat(l+iblock,l+iblock) = nonDegen(l,l);
				}
			} else {
				std::cerr << "QR eigenproblem failed to converge\n";
				std::cerr << "There may be no real solution to this eigenproblem\n";
				return matrix(1,1);
			}
			iblock += blocks[k].get_rows();
		}
	} 
	
	return evalsMat;
}


std::vector<matrix> matrix::GetEvecsAndEvalsFromEvals(const matrix& evalsMat_in) const {
          matrix evalsMat = evalsMat_in;
	  matrix evecsMat(get_rows(),get_rows());
	  matrix evecsMatNew(get_rows(),get_rows());

	  matrix evalsVec = matrix(evalsMat.get_rows(),1);
	  matrix null = matrix(get_rows(),1);

	  std::vector<std::vector<int> > mult;
          mult.push_back(std::vector<int>());

	  std::vector<double> evsort;
	  for(unsigned ii=0;ii<get_rows();ii++)
             evsort.push_back(evalsMat(ii,ii));
          std::sort(evsort.begin(),evsort.end());
	  for(unsigned ii=0;ii<get_rows();ii++)
             evalsMat(ii,ii) = evsort[ii];
	  for(unsigned ii=0;ii<get_rows();ii++){
		  evalsVec(ii,0) = evalsMat(ii,ii);
		  if(ii>0&&fabs(evalsMat(ii,ii)-evalsMat(ii-1,ii-1))>1e-8){
                    mult.push_back(std::vector<int>());
		  }
	          mult.back().push_back(ii);
	  }
          int irow = 0;
	  for(unsigned im=0;im<mult.size();im++){
		  matrix charDet = (*this)-evalsVec(mult[im][0],0)*matrix(get_rows(),kdelta);
		  matrix g = charDet.GaussianElimination();
		  std::vector<matrix> evecs = TriangularSolveBack(g,null,mult[im].size());
                  for(unsigned iv=0;iv<evecs.size();iv++){
                    evecs[iv] = evecs[iv].Orthonormalize();
		    evecsMatNew.SetColumn(irow,evecs[iv]);
                    irow++;
		  }
          }
	  std::vector<matrix> qre;
	  qre.push_back(evecsMatNew);
	  qre.push_back(evalsVec);
	  qre.push_back(evalsMat);
          return qre;
}

std::vector<matrix> matrix::Eigen() const{

  std::vector<matrix> eigenprob;

  // Whether to calc vectors or not, this will
  // be an option in class definition.

  if(get_rows()!=get_columns()){
     std::cerr << "Cannot calculate eigenvalues of non-square matrix!\n";
     return eigenprob;
  }
  //return QREigen();
  
  unsigned int i, j, p=0, q=0, iter=0, k;
  double amaxval=0.0, maxval=0.0, theta, bip, biq, cpp, cqq, cpq, cqp;
  matrix b=*this;
  matrix vectors(get_rows(),kdelta);
  matrix tmpvectors(get_rows(),get_columns());
  unsigned myCols = get_columns();
  unsigned myRows = get_rows();

  std::vector<int> maxvals(get_rows());

  for(i=1;i<myRows;i++){
    for(j=0;j<i;j++){
       if(fabs(b(i,j)-b(j,i))>1e-12){
          matrix evalsMat = QREigen();
          if(evalsMat.get_rows()==1) return eigenprob;
          std::vector<matrix> qre = GetEvecsAndEvalsFromEvals(evalsMat);
          return qre;
      }
    }
  }

/*
#ifdef _DO_TIMINGS_
  clock_t t1hess =clock();
#endif
  matrix hess=Hessenberg();
  b = hess;
#ifdef _DO_TIMINGS_
  clock_t t2hess =clock();
  std::cout << "Hessenberg " << ((t2hess-t1hess)*1000.0/CLOCKS_PER_SEC)/1000.0 << "\n";
#endif
*/


  double costheta;
  double sintheta;

  double veckp;
  double veckq;

  for(i=1;i<myRows;i++){
    amaxval=0.0, maxval=0.0;
    for(j=0;j<i;j++){
      if(b(i,j)>amaxval||-b(i,j)>amaxval) {
        maxvals[i]=j;
        maxval=b(i,j);
        if(maxval<0)
          amaxval = -maxval;
        else
          amaxval = maxval;
      }
    }
  }

#ifdef _DO_TIMINGS_
  clock_t t1=0,t2=0,t3=0,t1s,t2s,t3s;
#endif

  amaxval=0.0, maxval=0.0;
  while(iter<myCols*myCols*100){
#ifdef _DO_TIMINGS_
    t1s = clock();
#endif
    /*
    for(i=1;i<myRows;i++){
      for(j=0;j<i;j++){
	if(b(i,j)>amaxval||-b(i,j)>amaxval) {
	  maxval=b(i,j);
	  if(maxval<0)
		  amaxval =- maxval;
	  else
		  amaxval = maxval;
          //amaxval = fabs(maxval);
	  p=i;
	  q=j;
	}
      }
    }
    */

    for(i=1;i<myRows;i++){
      j=maxvals[i];
	if(b(i,j)>amaxval||-b(i,j)>amaxval) {
	  maxval=b(i,j);
	  if(maxval<0)
		  amaxval = -maxval;
	  else
		  amaxval = maxval;
          //amaxval = fabs(maxval);
	  p=i;
	  q=j;
	}
    }

#ifdef _DO_TIMINGS_
    t1 += clock() - t1s;
    
    t2s = clock();
#endif
    if(fabs(maxval)<TOL) {
      matrix evals(myCols,1);
      for(i=0;i<myCols;i++) evals(i,0) = b(i,i);
      eigenprob.push_back(vectors.Transpose());
      eigenprob.push_back(evals);
      eigenprob.push_back(b);
#ifdef _DO_TIMINGS_
      std::cout << "Time for section1 " << ((t1)*1000.0/CLOCKS_PER_SEC)/1000.0 << "\n";
      std::cout << "Time for section2 " << ((t2)*1000.0/CLOCKS_PER_SEC)/1000.0 << "\n";
      std::cout << "Time for section3 " << ((t3)*1000.0/CLOCKS_PER_SEC)/1000.0 << "\n";
      std::cout << "Iterations " <<  iter << "\n";
#endif
      return eigenprob;
    }
    
    if((fabs(b(p,p)-b(q,q)))>TOL){
      theta = (atan(2*b(p,q)/(b(p,p)-b(q,q))))*0.5;
    } else {
      theta = PIOVER4;
    }
    
    
    costheta = cos(theta);
    sintheta = sin(theta);
 
    for(k=0;k<myRows;k++){
      veckq = vectors(q,k);
      veckp = vectors(p,k);
      vectors(p,k) = veckq * sintheta + veckp * costheta;
      vectors(q,k) = veckq * costheta - veckp * sintheta;
    }

#ifdef _DO_TIMINGS_
    t2 += clock() - t2s;
    t3s = clock();
#endif
    for(i=0;i<q;i++){
        bip = b(p,i);
        biq = b(q,i);
	b(p,i) = bip*costheta+biq*sintheta;
	b(q,i) = -bip*sintheta+biq*costheta;
    }
    for(i=q+1;i<p;i++){
        bip = b(p,i);
        biq = b(i,q);
	b(p,i) = bip*costheta+biq*sintheta;
	b(i,q) = -bip*sintheta+biq*costheta;
    }
    for(i=p+1;i<myRows;i++){
        bip = b(i,p);
        biq = b(i,q);
	b(i,p) = bip*costheta+biq*sintheta;
	b(i,q) = -bip*sintheta+biq*costheta;
    }

    cpp =   b(p,p)*costheta*costheta
      + 2*b(p,q)*costheta*sintheta
      +   b(q,q)*sintheta*sintheta;

    cqq =   b(p,p)*sintheta*sintheta
      - 2*b(p,q)*costheta*sintheta
      +   b(q,q)*costheta*costheta;

    cpq = cqp = (b(q,q)-b(p,p))*costheta*sintheta
      + b(p,q)*(costheta*costheta - sintheta*sintheta); 
    
    b(p,p) = cpp;
    b(q,q) = cqq;
    b(p,q) = cpq;
    b(q,p) = cqp;
#ifdef _DO_TIMINGS_
    t3 += clock() - t3s;
    t1s = clock();
#endif

      amaxval=0.0, maxval=0.0;
      i=p;
        for(j=0;j<i;j++){
          if(b(i,j)>amaxval||-b(i,j)>amaxval) {
            maxvals[i]=j;
            maxval=b(i,j);
            if(maxval<0)
              amaxval = -maxval;
            else
              amaxval = maxval;
          }
        }
      amaxval=0.0, maxval=0.0;
      i=q;
        for(j=0;j<i;j++){
          if(b(i,j)>amaxval||-b(i,j)>amaxval) {
            maxvals[i]=j;
            maxval=b(i,j);
            if(maxval<0)
              amaxval = -maxval;
            else
              amaxval = maxval;
          }
        }
#ifdef _DO_TIMINGS_
    t1 += clock() - t1s;
#endif

    maxval=0.0;
    amaxval=0.0;

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

matrix matrix::Cholesky() const {
	/* Returns the lower triangle A in  U = AA^T */
	matrix a(get_rows(),get_columns());
	if(get_rows()!=get_columns()){
		std::cerr << "Cannot do Cholesky decomposition of non-square matrix\n";
		return a;
	}

	matrix u = *this;

	for(unsigned int i=0;i<u.get_rows();i++){
		for(unsigned int j=0;j<i;j++){
			double sum = u(i,j);
			for(unsigned k=0;k<i;k++) sum -= a(i,k)*a(j,k);
			a(i,j) = sum/a(j,j);
		}
		double sum = u(i,i);
		for(unsigned k=0;k<i;k++){
			sum -= a(i,k)*a(i,k);
		}
		a(i,i) = sqrt(sum);
	}

	/*
	a(0,0) = sqrt(u(0,0));

	if(get_rows()!=3)
		return a;

	if(fabs(a(0,0))<1e-7)
		a(1,0) = 0.0;	
	else
		a(1,0) = u(0,1)/a(0,0);

	if(fabs(a(0,0))<1e-7)
		a(2,0) = 0.0;
	else
		a(2,0) = u(0,2)/a(0,0);

	a(1,1) = sqrt(u(1,1)-a(1,0)*a(1,0));

	if(fabs(a(1,1))<1e-7)
		a(2,1) = 0.0;
	else
		a(2,1) = (u(2,1)-a(2,0)*a(1,0))/a(1,1);

	a(2,2) = sqrt(u(2,2)-a(2,0)*a(2,0)-a(2,1)*a(2,1));
	*/


	return a;
}

Cartesian matrix::GetRotationAxis() const {
  /* Returns the rotational aspect of a 4x4 matrix.
     Should be extended to 3x3 and 2x2, I guess.
   */
  Quat q(*this);
  Cartesian rotAxis = q.GetRotationAxis();
  return rotAxis;
}

matrix GetP(const matrix &v, unsigned m, unsigned n){
	matrix ident = matrix(m,n);
	//ident.Zero();
	unsigned maxv = m;
	if(n>m) maxv = n;
	for(unsigned i=0;i<maxv;i++)
		ident(i,i) = 1.0;
	matrix p = ident - (2.0*v*(v.Transpose()));
	return p;
}
	
matrix GetV(const matrix &a, unsigned i, double alpha, double r){
	matrix v = matrix(a.get_rows(),1);
	for(unsigned k=0;k<i+1;k++)
		v(k,0) = 0.0;
	if(fabs(r)>1e-12){
	v(i+1,0) = (a(i+1,i)-alpha)/(2*r);
	for(unsigned k=i+2;k<a.get_rows();k++)
		v(k,0) = a(k,i)/(2*r);
	}
	return v;
}

double sgn(double x){
  return fabs(x)/x;
}

std::vector<double> GetAlpahRHess(const matrix &a, int i){
	double s = 0;
	for(unsigned j=i+1;j<a.get_rows();j++)
		s = s + a(j,i) * a(j,i);

	double alpha;
	if(fabs(a(i+1,i))>1e-12)
		alpha = -sgn(a(i+1,i)) * sqrt(s);
	else
		alpha = sqrt(s);
	double r = sqrt(0.5*(alpha*alpha-alpha*a(i+1,i)));
	std::vector<double> ar;
	ar.push_back(alpha);
	ar.push_back(r);
	return ar;
}

matrix GetXQR(const matrix &a, int i){
	matrix x = matrix(a.get_rows(),1);
	for(unsigned k=0;k<a.get_rows();k++)
		x(k,0) = a(k,i);
	return x;
}

double GetNorm(const matrix &x){
	double s = 0;
	for(unsigned j=0;j<x.get_rows();j++)
		s = s + x(j,0) * x(j,0);
	return sqrt(s);
}

double GetAlpahQRFromX(const matrix &x){
	double s = 0;
	for(unsigned j=0;j<x.get_rows();j++)
		s = s + x(j,0) * x(j,0);
	return sqrt(s);
}

matrix GetUQR(const matrix &x, double alpha, int i){
	matrix e = matrix(x.get_rows(),1);
	//e.Zero();
	e(i,0) = 1;
	return x - e*alpha;
}


std::vector<matrix> matrix::QRDecomposition() const {
	std::vector<matrix> QR;
	matrix anew = *this;
	matrix atmp = *this;
	
	std::vector<matrix> qs;
	for(unsigned i=0;i<get_columns()-1;i++){
		matrix x = GetXQR(anew,0);
		double alpha = GetAlpahQRFromX(x);
		matrix u = GetUQR(x,alpha,0);
		matrix v = u;
		double norm = GetNorm(u);
	       	if(fabs(norm)>1e-12) v /= norm;
		matrix q1 = GetP(v,anew.get_rows(),anew.get_columns());
		matrix q2;
		if(i>0){
			matrix ident = matrix(i,kdelta);
			q2 = DirSum(ident,q1);
		}else{
			q2 = q1;
		}
		qs.push_back(q2.Transpose());
		atmp = q2 * atmp;
		matrix m11 = anew.MinorMatrix(atmp,0,0);
		for(unsigned k=0;k<i;k++)
			m11 = anew.MinorMatrix(m11,0,0);
		anew = m11;
	}

	matrix q = qs[0];
	for(unsigned k=1;k<qs.size();k++){
		//q = q * qs[k];	
		q.MultiplyAndSetRight(qs[k]);
        }
	matrix r = q.Transpose()*(*this);

	QR.push_back(q);
	QR.push_back(r);
	return QR;
}

matrix matrix::Hessenberg() const {
	matrix anew = *this;

	double alpha,r;
	std::vector<double> ar;
	for(unsigned i=0;i<anew.get_rows()-1;i++){
		ar = GetAlpahRHess(anew,i);
		alpha = ar[0];
		r = ar[1];
		if(fabs(r)>1e-12){
			matrix v = GetV(anew,i,alpha,r);
			matrix p = GetP(v,anew.get_rows(),anew.get_columns());
			anew = p*anew*p;
		}
	}

	return anew;
}

matrix matrix::GetRow(int i) const {
	matrix a(1,get_columns());
	for(unsigned j=0;j<get_columns();j++){
		a(0,j) = (*this)(i,j);
	}
	return a;
}

matrix matrix::GetColumn(int j) const {
	matrix a(get_rows(),1);
	for(unsigned i=0;i<get_rows();i++){
		a(i,0) = (*this)(i,j);
	}
	return a;
}

void matrix::SetColumn(int j, const matrix& colVec){
	for(unsigned i=0;i<get_rows();i++){
		(*this)(i,j) = colVec(i,0);
	}
}

void matrix::SetRow(int i, const matrix& rowVec){
	for(unsigned j=0;j<get_columns();j++){
		(*this)(i,j) = rowVec(0,j);
	}
}

bool matrix::isNull() const {
	double rms = 0.0;
	for(unsigned i=0;i<get_rows();i++){
		for(unsigned j=0;j<get_columns();j++){
			rms += (*this)(i,j) * (*this)(i,j);
		}
	}
	if(rms<1e-10) return true;
	return false;
}

std::vector<matrix> matrix::Bidiagonalize() const {
	//Golub-Kahan-Lanczos Bidiagonalization Procedure
	//http://www.cs.utk.edu/~dongarra/etemplates/node198.html
	
	matrix a = matrix(get_rows(),get_columns());
	matrix v = matrix(get_columns(),get_columns());
	matrix u = matrix(get_rows(),get_rows());
	matrix vk = matrix(get_columns(),1);
	matrix uk = matrix(get_rows(),1);
	//vk.Zero();
	//v.Zero();
	vk(0,0) = 1.0;
	v(0,0) = 1.0;

	double bk = 0;

	for(unsigned k=0;k<get_columns()-1;k++){
		uk = (*this)*vk - uk*bk;
		double ak = GetNorm(uk);
		//std::cout << "ak" << ": " << ak << "\n";
		if(fabs(ak)>1e-12) uk = uk/ak;
		//std::cout << "vk*ak" << ": " << vk*ak << "\n";
		//std::cout << "Transpose()" << ": " << Transpose() << "\n";
		//std::cout << "uk" << ": " << uk << "\n";
		vk = (Transpose())*uk  - vk*ak;
		//std::cout << "(Transpose())*uk" << ": " << (Transpose())*uk << "\n";
		//std::cout << "vk" << ": " << vk << "\n";
		bk = GetNorm(vk);
		//std::cout << "bk" << ": " << bk << "\n";
		if(fabs(bk)>1e-12) vk = vk/bk;
		if(k<u.get_columns()){
			for(unsigned i=0;i<u.get_rows();i++)
				u(i,k) = uk(i,0);
		}
		for(unsigned i=0;i<v.get_rows();i++)
			v(i,k+1) = vk(i,0);
		a(k,k) = ak;
		if(k<a.get_columns()-1) a(k,k+1) = bk;
	}

	/*
	print "UU*"; sys.stdout.flush()
	(u*u.Transpose()).Print(); streams.FlushCStdOut(); print; sys.stdout.flush()
	print; sys.stdout.flush()

	print "VV*"; sys.stdout.flush()
	(v*v.Transpose()).Print(); streams.FlushCStdOut(); print; sys.stdout.flush()
	print; sys.stdout.flush()

	print "U*MV"; sys.stdout.flush()
	(u.Transpose()*m*v).Print(); streams.FlushCStdOut(); print; sys.stdout.flush()
	print; sys.stdout.flush()

	print "ATA"; sys.stdout.flush()
	(a.Transpose()*a).Print(); streams.FlushCStdOut(); print; sys.stdout.flush()
	print; sys.stdout.flush()

	print "A"; sys.stdout.flush()
	(a).Print(); streams.FlushCStdOut(); print; sys.stdout.flush()
	print; sys.stdout.flush()
	*/
	std::vector<matrix> uav;
	uav.push_back(u);
	uav.push_back(a);
	uav.push_back(v);
	return uav;
}

matrix matrix::Orthonormalize() const {
  // Takes a matrix which contains a set of column vectors
  // and returns a matrix containg orthonormalized column vectors.

  matrix ortho = *this;

  ortho.SetColumn(0,GetColumn(0));
  for(unsigned j=0;j<get_columns();j++){
    matrix vj = GetColumn(j);
    if(vj.isNull()) {
       for(unsigned knorm=0;knorm<vj.get_rows();knorm++){
          bool have_norm = true;
          vj(knorm,0) = 1.0;
          for(unsigned jnorm=0;jnorm<j;jnorm++){
             double n1 = (ortho.GetColumn(jnorm).Transpose()*vj)(0,0);
             if(fabs(n1)>1e-12){
                have_norm = false;
             }
          }
          if(!have_norm)
           vj.Zero();
          else
           break;
       }
    }
    matrix uj = vj;
    for(unsigned i=0;i<j;i++){
      //double proj = (vj.Transpose()*ortho.GetColumn(i))(0,0);
      double n1 = (ortho.GetColumn(i).Transpose()*ortho.GetColumn(i))(0,0);
      if(fabs(n1)>1e-12){
        matrix projuv = (vj.Transpose()*ortho.GetColumn(i))(0,0)/n1 * ortho.GetColumn(i);
        uj -= projuv;
      }
    }
    double n = GetNorm(uj);
    if(fabs(n)>1e-12) uj = uj/n;
    ortho.SetColumn(j,uj);
  }

  return ortho;
}

std::vector<matrix> matrix::SVD() const {

	// Ideally we would want to do bidiagonalization first.
	// This algorithm does not yet deal with singular values of zero.

	matrix vin = Transpose() * (*this);

	std::vector<matrix> veig = vin.Eigen();
	veig = SortEigenvalues(veig,true);

	matrix vt = veig[0].Transpose();
	//std::cout << "VT(c++)\n" << vt << "\n";

	matrix s = matrix(get_rows(),get_columns());
	matrix sinv = matrix(vt.get_rows(),vt.get_rows());

	unsigned mind = get_rows();
	if(s.get_columns()<mind) mind = get_columns();

	bool inorm = false;
	for(unsigned i=0;i<mind;i++){
		if(fabs(veig[1](i,0))> 1e-12){
			sinv(i,i) = 1.0/sqrt(veig[1](i,0));
			s(i,i) = sqrt(veig[1](i,0));
		} else{
		// We have to use Gram-Schmidt to orthonormalize!
		inorm = true;
		}
	}

	matrix u = (*this) * (vt.Transpose()) * sinv;
	matrix unew = matrix(u.get_rows(),u.get_rows());
	for(unsigned i=0;i<u.get_rows();i++){
		for(unsigned j=0;j<u.get_rows();j++){
			unew(i,j) = u(i,j);
		}
	}
	if(inorm)
		u = unew.Orthonormalize();
	else
		u = unew;
	//std::cout << "U(c++)\n" << u << "\n";
	std::vector<matrix> usvt;
	usvt.push_back(u);
	usvt.push_back(s);
	usvt.push_back(vt);
	return usvt;
}

std::vector<matrix> matrix::LUDecomposition() const {

	matrix u = *this;

	std::vector<matrix> ls;
	std::vector<matrix> ps;

	matrix p;
	matrix l0;
	matrix l;

	std::vector<matrix> lup;

	for(unsigned n=0;n<get_rows()-1;n++){
		p = matrix(get_rows(),get_rows());
		//p.Zero();
		for(unsigned i=0;i<get_rows();i++)
			p(i,i) = 1;

		l = matrix(get_rows(),get_rows());
		//l.Zero();

		for(unsigned i=0;i<get_rows();i++)
			l(i,i) = 1;

		if(fabs(u(n,n))< 1e-7){

			bool success=false;
			for(unsigned i=n+1;i<get_rows();i++){
				if(fabs(u(i,n)) > 1e-7){
					u.SwitchRows(i,n);
					p.SwitchRows(i,n);
					success=true;
					break;
				}
			}
			if(!success){
				// No pivot, so give up at this point.

				p = ps[0];
				l0 = p * ls[0];

				for(unsigned i=1;i<ls.size();i++){
					l0 = l0 * ps[i]* ls[i];
					p = ps[i] * p;
				}

				l0 = p * l0;

				lup.push_back(l0);
				lup.push_back(u);
				lup.push_back(p);
				return lup;
			}
		}

		for(unsigned i=n+1;i<get_rows();i++){
			l(i,n) = -(u(i,n)/u(n,n));
		}

		u = l * u;
		ls.push_back(l.Inverse());
		ps.push_back(p);
	}

	p = ps[0];
	l0 = p * ls[0];

	for(unsigned i=1;i<ls.size();i++){
		l0 = l0 * ps[i]* ls[i];
		p = ps[i] * p;
	}

	l0 = p * l0;

	lup.push_back(l0);
	lup.push_back(u);
	lup.push_back(p);
	
	return lup;

}

