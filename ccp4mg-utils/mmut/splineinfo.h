/*
     mmut/splineinfo.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009 University of York

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

#ifndef _CCP4MG_SPLINEINFO_
#define _CCP4MG_SPLINEINFO_

#include "atom_util.h"

struct SplineInfo {
  std::vector<std::vector<Cartesian> > splines;
  std::vector<std::vector<Cartesian> > n1_splines;
  std::vector<std::vector<Cartesian> > n2_splines;
  std::vector<std::vector<std::vector<int> > > secstr_indices;
  std::vector<std::vector<Cartesian> > colours;
  std::vector<std::vector<Cartesian> > nasplines;
  std::vector<std::vector<Cartesian> > n1_nasplines;
  std::vector<std::vector<Cartesian> > n2_nasplines;
  std::vector<std::vector<Cartesian> > nacolours;
};

SplineInfo GetSplineInfo(CMMANManager *molH, int atom_selHnd, AtomColourVector *atm_col_vect, int spline_accu , int udd_chain = -1, int udd_CA = -1, int flatten_beta_sheet =0, int flatten_loop=0, int smooth_helix=0);
int GetCAFromSelection(CMMANManager *molH, int atom_selHnd_in);

#endif
