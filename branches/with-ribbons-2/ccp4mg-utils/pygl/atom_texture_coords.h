/*
     pygl/atom_texture_coords.h: CCP4MG Molecular Graphics Program
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


#ifndef _CCP4MG_TEXCOORDS_
#define _CCP4MG_TEXCOORDS_


#include <map>
#include <vector>
#include <string>

class TexCoords {
  static std::map<std::string,std::vector<double> > coords;
 public:
  TexCoords(){};
  ~TexCoords(){};

  static std::vector<double> GetCoords(const std::string &name);
  static void AddCoords(const std::string &name, double tsx, double tex, double tsy, double tey);
  static void AddCoords(const std::string &name, const std::vector<double>&c);
};

#endif
