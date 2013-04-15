/*
     mmut/mmut_util.cc: CCP4MG Molecular Graphics Program
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


#include <iostream>
#include <iomanip>
#include "mmut_util.h"
#include <mmdb/mmdb_manager.h>

using namespace std;

void printrealp(realtype *a, int size){
 
  int i,k;

  k = 0;
  while(k<size){
    for(i=0;i<6&&k<size;i++,k++){
      cout << setw(10) << a[k];      
    }
    cout << endl;
  }
  cout << endl;

}

void printrealpp(realtype **a, int rows, int cols){
 
  int i,k,l;

  k = 0;
  for(k=0;k<rows;k++){
    cout << "Row: " << k << endl;
    l = 0;
    while(l<cols){
      for(i=0;i<6&&l<cols;i++,l++){
	cout << setw(10) << a[k][l];      
      }
      cout << endl;
    }
    cout << endl;
  }
  cout << endl;

}

void printintp(int *a, int size){
 
  int i,k;

  k = 0;
  while(k<size){
    for(i=0;i<6&&k<size;i++,k++){
      cout << setw(10) << a[k];      
    }
    cout << endl;
  }
  cout << endl;

}

void printintpp(int **a, int rows, int cols){

  int i,k,l;

  k = 0;
  for(k=0;k<rows;k++){
    cout << "Row: " << k << endl;
    l = 0;
    while(l<cols){
      for(i=0;i<6&&l<cols;i++,l++){
	cout << setw(10) << a[k][l];      
      }
      cout << endl;
    }
    cout << endl;
  }
  cout << endl;

}

std::vector<Cartesian> PPCAtomsToCartesians(int natoms, PPCAtom atoms){

  std::vector<Cartesian> carts(natoms);

  for(int i=0;i<natoms;i++){
    carts[i] = Cartesian(atoms[i]->x,atoms[i]->y,atoms[i]->z);
  }

  return carts;

}
