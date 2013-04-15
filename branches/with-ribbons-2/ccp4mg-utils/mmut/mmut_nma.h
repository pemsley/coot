/*
     mmut/mmut_nma.h: CCP4MG Molecular Graphics Program
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

#ifndef _CCP4_MMUT_NMA_H_
#define _CCP4_MMUT_NMA_H_

enum {MMUT_ANM,MMUT_GNM,MMUT_NMA_NONE};

#include <string>
#include <vector>
#include <matrix.h>
#include <cartesian.h>

class NormalModeDisplacements {
  std::vector<std::vector<std::vector<Cartesian> > > displacements;
 public:
  NormalModeDisplacements();
  const std::vector<Cartesian>& getDisplacements(unsigned mode, int step) const;
  const std::vector<std::vector<Cartesian> >& getDisplacements(unsigned mode) const;
  void addMode(){displacements.push_back(std::vector<std::vector<Cartesian> >());};
  int addDisplacements(unsigned mode, const std::vector<Cartesian> &disp);
  unsigned getNumberOfModes();
  void clear() {displacements.clear();} ;
};

class NormalModeAnalysis {
  matrix hessian;
  std::vector<matrix> eigen;
  int _type;
  matrix hinv;
  void NMA_GNM(const std::vector<Cartesian> &carts, const double cutoff);
  void NMA_ANM(const std::vector<Cartesian> &carts, const double cutoff);
 public:
  void Calculate(const std::vector<Cartesian> &carts, const int type=MMUT_GNM, const double cutoff=-1.0);
  NormalModeAnalysis(const std::vector<Cartesian> &carts, const int type=MMUT_GNM, const double cutoff=-1.0);
  NormalModeAnalysis();
  //NormalModeAnalysis(const std::vector<Cartesian> &carts, const std::string &type="GNM", const double cutoff=-1.0);
  std::vector<matrix> GetEigen() const {return eigen;};
  // 3N-5 (for ANM) vectors of n steps. We calculate a complete cycle of the mode. This is first approxn.
  NormalModeDisplacements GetDisplacements(const int nsteps) const;
  matrix GetHessian() const {return hessian;};
  matrix GetInverseHessian() const;
  matrix GetDiagonalInverseHessian() const;
  void StoreInverseHessian();
  int GetType() const {return _type;};
  std::vector<double>  GetBValues() const ;
  std::vector<std::vector<double> > GetModeShapes(const double gammainv=1.0) const ;
  std::vector<std::vector<Cartesian> > GetModeShapesAsCartesians(const double gammainv) const ;
  matrix GetCorrelations(const double gammainv) const;
};
#endif //_CCP4_MMUT_NMA_H_
