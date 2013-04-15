/*
     util/mgutil.cc: CCP4MG Molecular Graphics Program
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


#include <string>
#include <stdio.h>
#include <ctype.h>

/*
std::string IntToString(int fn_isize){

  char *tmpstr = new char[20];
  sprintf(tmpstr,"%d",fn_isize);
  std::string fn_size = std::string(tmpstr);
  delete [] tmpstr;
  return fn_size;
}
*/

std::string ToUpper(const std::string &in){
  std::string upper;
  std::string::const_iterator s=in.begin();
  while(s!=in.end()){
    upper.push_back(toupper(*s));
    s++;
  } 
  return upper;
}

std::string IntToString(int fn_isize,const char *format){

  char *tmpstr = new char[20];
  sprintf(tmpstr,format,fn_isize);
  std::string fn_size = std::string(tmpstr);
  delete [] tmpstr;
  return fn_size;
}

std::string FloatToString(float fn_isize,const char *format){

  char *tmpstr = new char[20];
  sprintf(tmpstr,format,fn_isize);
  std::string fn_size = std::string(tmpstr);
  delete [] tmpstr;
  return fn_size;
}

/*
std::string FloatToString(float fn_isize){

  char *tmpstr = new char[20];
  sprintf(tmpstr,"%f",fn_isize);
  std::string fn_size = std::string(tmpstr);
  delete [] tmpstr;
  return fn_size;
}
*/

