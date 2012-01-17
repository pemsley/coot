/*
     util/CParamsManager.h: CCP4MG Molecular Graphics Program
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
#ifndef _CCP4MG_CPARAMSMANAGER_
#define _CCP4MG_CPARAMSMANAGER_

#include <map>
#include <string>

class CParamsManager {
 private:
  std::map < std::string,int > Ints;
  std::map < std::string, float > Floats;
 public:
  CParamsManager();
  ~CParamsManager();
  void SetFloat ( const std::string &key, const float value ){ Floats[key] = value;}
  void SetInt ( const std::string &key, const int value ) {Ints[key]=value;  }
  float GetFloat ( const std::string &key ) const;
  int GetInt ( const std::string &key )const;
};

#endif
