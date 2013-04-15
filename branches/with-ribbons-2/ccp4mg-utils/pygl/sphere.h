/*
     pygl/sphere.h: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2005 University of York, CCLRC

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


#ifndef __SPHERE__
#define __SPHERE__

void sphere(int accu, float radius, bool force_dl=false);

#define X .525731112119133606 
#define Z .850650808352039932

static GLfloat icosa[12][3] = {    
  {-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},    
  {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},    
  {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0} 
};

static GLuint tindices[20][3] = { 
  {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},    
  {8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},    
  {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6}, 
  {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11} };

#endif
