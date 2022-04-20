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


#ifndef _CCP4MG_MATRIX_
#define _CCP4MG_MATRIX_

#include <vector>
#include <string>

class Quat;

class matrix {
  std::vector<std::vector<double> > mat;

public:
   matrix(){/*std::cout << "matrix default constructor\n"; std::cout.flush();*/};
   ~matrix();
   matrix(const matrix &a);
   matrix(const Quat &q);
   matrix& operator=(const matrix &a);
   matrix(int unsigned x,unsigned int y);
   matrix(int unsigned x,unsigned int y, double *array);
   matrix(unsigned int x,unsigned int y, const std::vector<double> &mat_in);
   explicit matrix(const std::vector<std::vector<double> > &mat_in) : mat(mat_in) {};
   matrix(unsigned int x,double(*fun)(int i,int j,int k));

   double *to_dp();
   void Zero();

  unsigned int get_rows() const;
  unsigned int  get_columns() const;

  double Trace() const;
  matrix Transpose() const;

  friend matrix pow(const matrix &a, double power);
  // Pointer to function, for sin,cos, tan, exp, etc.,
  matrix fun(const matrix &a,double(*fun)(double)) const;
  matrix fun(const matrix &a, const std::string &fun) const;
  

  // Operators 
  
  // Binary +-*/ 
  matrix operator+ (const matrix &b) const ;
  matrix operator- (const matrix &b) const ;
  matrix operator* (const matrix &b) const ;
  matrix operator/ (const matrix &b) const ;
  
  // Combine assignment and operation
  void operator+= (const matrix &a);
  void operator-= (const matrix &a);
  void operator*= (const matrix &a);
  void operator/= (const matrix &a);

  matrix& operator*= (double val);
  matrix& operator/= (double val);
  
  // Unary -
  matrix operator- () const;
  
  // Multiply/divide by constant and vv.
  friend matrix operator* (double val,   const matrix&);
  friend matrix operator/ (double val,   const matrix&);
  friend matrix operator* (const matrix &a, double val);
  friend matrix operator/ (const matrix &a, double val);
  
  //Matrix Stuff
  matrix DirSum(const matrix &a, const matrix &b);
  //matrix DirProd(const matrix &a, const matrix &b);
  std::vector<matrix> Eigen() const;
  std::vector<matrix> SortEigenvalues(const std::vector<matrix>&eigen) const;

  static matrix MinorMatrix(const matrix &a, const unsigned int &row, const unsigned int &col);
  static double Minor(const matrix &a, const unsigned int &row, const unsigned int &col);
  static matrix TriangularSolveForward(const matrix &a, const matrix &b);
  static matrix TriangularSolveBack(const matrix &a, const matrix &b);
  static matrix LUDecomposition(const matrix &a, std::vector<int> &perm, int &parity);
  static matrix LUSubstitution(const matrix &a, const matrix &b, const std::vector<int> &perm);
  static matrix SolveLinearEquations(const matrix &a, const matrix &b);
  static matrix LUMult(const matrix &a, const matrix &b, const std::vector<int> &perm);

  double Determinant() const;
  matrix TriangularMatrix() const;
  matrix Inverse() const;
  matrix BlockMatrix(const std::vector<std::vector<matrix> > &blocks);
  matrix GetUpperTriangle() const;
  matrix GetLowerTriangle(int lu=0) const;

  void SwitchRows(const unsigned int &i1, const unsigned int &i2);

  // overloaded round brackets for indices
  double &operator()(int i,int j);
  double const &operator()(int i,int j) const;

  // stream operator
  friend std::ostream& operator<<(std::ostream &c, matrix a);

};

// overloaded round brackets for indices
inline double &matrix::operator()(int i,int j){
  return mat[i][j];
}

inline const double &matrix::operator()(int i,int j) const{
  return mat[i][j];
}

double kdelta(int i,int j,int k);

#endif
