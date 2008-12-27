/*
     pygl/bezier.cc: CCP4MG Molecular Graphics Program
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

#include <vector>
#include <fstream>
#include <cartesian.h>

std::vector<Cartesian> BezierCurve(const std::vector<Cartesian> &carts, const unsigned int accu){
  std::vector<Cartesian> spline;
  unsigned nsteps = accu * (carts.size()-1);
  double tstep = 1.0 / double(nsteps);

  double t = 0.0;
  for (unsigned k = 0; k <= nsteps; k++){
   std::vector<Cartesian> cartsi = carts;

   for (int j = int(carts.size()-1); j > 0; j--)
    for (int i = 0; i < j; i++)
     cartsi[i] = (1-t)* cartsi[i] + t*cartsi[i+1];

   spline.push_back(cartsi[0]);
   t += tstep;
  }

  return spline;

}
