/*
     pygl/billboard.cc: CCP4MG Molecular Graphics Program
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

#ifdef _WIN32
#include <windows.h>
#endif

#include <string>
#include <vector>
#include "cdisplayobject.h"
#include "cprimitive.h"
#include "cartesian.h"

void AddBillBoardImage(Displayobject &obj, std::string filename, std::vector <Cartesian> carts){
  BillBoard *bill;
  bill = new BillBoard(carts,carts[0],filename.c_str());
  obj.add_primitive(bill);
}

void AddBillBoardImage(Displayobject &obj, const char* filename, double *pcarts){
  std::vector <Cartesian> carts;
  carts.push_back(Cartesian(pcarts[0],pcarts[1],pcarts[2]));
  carts.push_back(Cartesian(pcarts[3],pcarts[4],pcarts[5]));
  AddBillBoardImage(obj,std::string(filename),carts);
}
