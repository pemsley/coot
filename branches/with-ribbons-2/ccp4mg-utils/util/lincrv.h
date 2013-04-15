/*
     util/lincrv.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York

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


#ifndef __LINCRV__
#define __LINCRV__
#include "cartesian.h"
#include <vector>

/****** lincrv.h ******/
/* Ken Shoemake, 1994 */

typedef struct CART3D {
  double x;
  double y;
  double z;
} CART3D ;

void DialASpline(double t, const std::vector<double> &a,  const std::vector<Cartesian> &p, int Cn, int interp, std::vector<Cartesian> &output, const int idx, std::vector<CART3D> &work);
#endif

