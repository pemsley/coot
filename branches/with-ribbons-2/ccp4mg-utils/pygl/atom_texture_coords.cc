/*
     pygl/atom_texture_coords.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2011 University of York

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
#include "atom_texture_coords.h"

std::map<std::string,std::vector<double> > TexCoords::coords;

std::vector<double> TexCoords::GetCoords(const std::string &name) {
  std::map<std::string, std::vector<double> >::iterator it = coords.find(name);
  if(it==coords.end()){
    it = coords.find(" X");
    if(it!=coords.end()){
      //std::cout << "Return X for " << name << "\n";
      return coords[" X"];
    }
    std::vector<double> c;
    c.push_back(0.875);
    c.push_back(1.0);
    c.push_back(0.875);
    c.push_back(0.9995);
    //std::cout << "Return def for " << name << "\n";
    return c;
  }
  //std::cout << "Return name for " << name << "\n";
  return coords[name];
}

void TexCoords::AddCoords(const std::string &name, double tsx, double tex, double tsy, double tey){
  std::vector<double> c;
  c.push_back(tsx);
  c.push_back(tex);
  c.push_back(tsy);
  c.push_back(tey);
  coords[name] = c;
}

void TexCoords::AddCoords(const std::string &name, const std::vector<double>&c){
  coords[name] = c;
}
