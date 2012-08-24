/*
     mmut/mmut_util.h: CCP4MG Molecular Graphics Program
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


#ifndef __MMUT_Util__
#define __MMUT_Util__

#include <mmdb_manager.h>
#include <vector>
#include <cartesian.h>


inline void setrealp(realtype *a, int i, realtype val){

  a[i] = val;

}

inline void setrealpp(realtype **a, int i, int j, realtype val){

  a[i][j] = val;

}

inline void setintp(int *a, int i, int val){

  a[i] = val;

}

inline void setintpp(int **a, int i, int j, int val){

  a[i][j] = val;

}

inline realtype getrealp(realtype *a, int i){

  return a[i];

}

inline realtype getrealpp(realtype **a, int i, int j){

  return a[i][j];

}

inline int getintp(int *a, int i){

  return a[i];

}

inline int getintpp(int **a, int i, int j){

  return a[i][j];

}
void printrealp(realtype *a, int size);
void printrealpp(realtype **a, int rows, int cols);
void printintp(int *a, int size);
void printintpp(int **a, int rows, int cols);

std::vector<Cartesian> PPCAtomsToCartesians(int natoms, PPCAtom atoms);

#endif
