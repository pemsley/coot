/*
     mmut/mmut_nma.h: CCP4MG Molecular Graphics Program
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

#ifndef _CCP4_MMUT_NMA_H_
#define _CCP4_MMUT_NMA_H_

enum {MMUT_ANM,MMUT_GNM};

#include <string>
#include <vector>
#include <matrix.h>
#include <cartesian.h>

class NormalModeAnalysis {
  matrix hessian;
  std::vector<matrix> eigen;
  int _type;
  matrix hinv;
  void NMA_GNM(const std::vector<Cartesian> &carts, const double cutoff);
  void NMA_ANM(const std::vector<Cartesian> &carts, const double cutoff);
 public:
  NormalModeAnalysis(const std::vector<Cartesian> &carts, const int type=MMUT_GNM, const double cutoff=-1.0);
  //NormalModeAnalysis(const std::vector<Cartesian> &carts, const std::string &type="GNM", const double cutoff=-1.0);
  std::vector<matrix> GetEigen() const {return eigen;};
  // 3N-5 (for ANM) vectors of n steps. We calculate a complete cycle of the mode. This is first approxn.
  std::vector<std::vector<Cartesian> > GetDisplacements(const int nsteps) const;
  matrix GetHessian() const {return hessian;};
  matrix GetInverseHessian() const;
  matrix GetDiagonalInverseHessian() const;
  void StoreInverseHessian();
  int GetType() const {return _type;};
  std::vector<double>  GetBValues() const ;
  std::vector<std::vector<double> > GetModeShapes(const double gammainv=1.0) const ;
};
#endif //_CCP4_MMUT_NMA_H_
