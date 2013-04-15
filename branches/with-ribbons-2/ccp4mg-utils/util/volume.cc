/*
     util/volume.cc: CCP4MG Molecular Graphics Program
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


#include "volume.h"
#include "plane.h"
#include "cartesian.h"

std::ostream& operator<<(std::ostream& c, Volume a){
  for(int i=0;i<a.GetNumberOfPlanes();i++)
    c << "Plane(" << i << "): " << a.GetPlane(i);
  return c;
}

Volume::Volume(){
}

Volume::~Volume(){
}

void Volume::AddPlane(Plane plane){
   planes.push_back(plane);
}

int Volume::PointInVolume(double *point)const {
  Cartesian cart_point(point);
  return PointInVolume(cart_point);
}

int Volume::PointInVolume(const Cartesian &p) const{
  int involume = 1;
  for(int i=0;i<GetNumberOfPlanes();i++){
    Plane plane = GetPlane(i);
    Cartesian point = plane.find_points_on_plane()[0];
    Cartesian n = plane.get_normal();
    Cartesian p2prim = point-p;
    n.normalize();
    p2prim.normalize();
    if(n.DotProduct(n,p2prim)>0.0){
      involume = 0;
    }
  }
  return involume;
}

int Volume::GetNumberOfPlanes(void) const {
  return planes.size();
}
