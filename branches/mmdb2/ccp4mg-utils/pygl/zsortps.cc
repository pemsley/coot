/*
     pygl/zsortps.cc: CCP4MG Molecular Graphics Program
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


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <algorithm>
#ifdef __GNUC__
#if (__GNUC__==2)
#include <streambuf.h>
#else
#include <streambuf>
#endif
#else
#include <streambuf>
#endif

class valueline{
 public:
  double val;
  std::string str;
  valueline(double val_in, std::string str_in){val = val_in; str = str_in;};
};

class sort_valueline{
  public:
    int operator()(const valueline &p1, const valueline &p2) const
       { return p1.val < p2.val; }
};


void zsortps(const std::string &filename, int fog){

  std::string s="";
  std::vector<std::string> oh,ot;
  char c;
  std::vector<valueline> zslines;
  std::ifstream str(filename.c_str());
  //std::ifstream str2(filename.c_str()); // If only seekg worked.
  int insortlines = 0;

  /*
  double fog_start = 1;
  double fog_end = 1000;
  while(str2.get(c)){
     s += c;
     if(c=='\n'){
       if(s.find("FOG_START:")!=std::string::npos){
	 std::string sub = s.substr(s.find("FOG_START:")+10);
	 std::istringstream istr(sub);
	 istr.get(c);
	 istr >> fog_start; 
         std::cout << "fog_start: " << fog_start << "\n";
       }else if(s.find("FOG_END:")!=std::string::npos){
	 std::string sub = s.substr(s.find("FOG_END:")+7);
	 std::istringstream istr(sub);
	 istr.get(c);
	 istr >> fog_end; 
         std::cout << "fog_end: " << fog_end << "\n";
       }
       s = "";
     }
  }
  str2.close();
  */

  while(str.get(c)){
     s += c;
     if(c=='\n'){
       size_t t = s.find("ZVALUE");
       if(t!=std::string::npos){
	 insortlines = 1;
	 std::string sub = s.substr(s.find("ZVALUE")+8);
         zslines.push_back(valueline(atof(sub.c_str()),s));
       }else{
         size_t t = s.find("TEXT");
         if(t!=std::string::npos){
	     ot.push_back(s);
	 }else{
	   if(!insortlines)
             oh.push_back(s);
	   else
             ot.push_back(s);
	 }
       }
       s = "";
     }
  }
  if(zslines.size()==0)
     return;
  std::sort(zslines.begin(),zslines.end(),sort_valueline());

  str.close();
#ifdef __GNUC__
#if (__GNUC__==2)
  std::ofstream ostr(filename.c_str(),_IO_TRUNC);
#else
  std::ofstream ostr(filename.c_str(),std::ios_base::trunc);
#endif
#else
  std::ofstream ostr(filename.c_str(),std::ios_base::trunc);
#endif

  for(unsigned int i = 0; i < oh.size(); i++)
    ostr << oh[i];

  if(fog){
    double range = zslines[zslines.size()-1].val - zslines[0].val;
    for(unsigned int i = 0; i < zslines.size(); i++){
      double fog_scale = (zslines[i].val - zslines[0].val)/range;
      //std::cout << fog_scale << " " << (fog_end-zslines[i].val)/(fog_end-fog_start) << "\n";
      size_t t = zslines[i].str.find("srgb");
      if(t!=std::string::npos&&t>5){
	size_t g_pos = zslines[i].str.rfind(" ",t-2);
	if(g_pos!=std::string::npos&&g_pos>3){
	  size_t b_pos = zslines[i].str.rfind(" ",g_pos-2);
	  if(b_pos!=std::string::npos&&b_pos>1){
	    size_t r_pos = zslines[i].str.rfind(" ",b_pos-2);
	    if(r_pos==std::string::npos){
	      r_pos = 0;
	    }
	    std::istringstream istr(zslines[i].str.substr(r_pos));
	    double r; istr >> r; 
	    double g; istr >> g; 
	    double b; istr >> b;
	    ostr << zslines[i].str.substr(0,r_pos) << " " << fog_scale*r+(1-fog_scale) << " " << fog_scale*g+(1-fog_scale) << " " << fog_scale*b+(1-fog_scale) << " " << zslines[i].str.substr(t);
	  }
	}
      }
    }
  }else{
    for(unsigned int i = 0; i < zslines.size(); i++)
      ostr << zslines[i].str;
  }

  for(unsigned int i = 0; i < ot.size(); i++)
    ostr << ot[i];

  //std::cout filename << " has " << oh.size()+ot.size()+zslines.size() << " lines\n";
  //std::cout filename << " has " << zslines.size() << " lines with %ZVALUE:\n";

}
