/*
     mmut/mmut_morph.h: CCP4MG Molecular Graphics Program
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
#include <cartesian.h>
#include <mgtree.h>

class CMMUTMorph {
 public:
  std::vector<std::vector<Cartesian> > MorphCartesians(const std::vector<Cartesian> &a, const std::vector<Cartesian> &b, int nsteps);
  std::vector<Tree> MorphTrees(const Tree &a, const Tree &b, int nsteps=10);
};
