/*
     util/rgbreps.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York
     Copyright (C) 2012 STFC

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
#include <map>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>

#include "rgbreps.h"

using namespace std;

unsigned int RGBReps::GetNumberOfColours(){
  return colours.size();
}

std::string RGBReps::GetNameOfColour(unsigned int i){
  /* Return name of colour as in order of input file */
  if(i>=colours.size()){
    //cerr << "Warning colour out of range, returning default " << i << " " << GetNumberOfColours() << "\n";
    return std::string("default");
  }
  return colours_name_map[i];
}

int  RGBReps::GetColourNumber(std::string name ) {
  for (int i=0;i<int(colours_name_map.size());i++) {
    if (!name.compare(colours_name_map[i])) {
      return i;
    }
  }
  return -1;
}

std::vector<int>  RGBReps::GetColourNumbers(std::vector<std::string> names ) {
  std::vector<int> codes(names.size());
  for (int n=0;n<int(names.size());n++) {
    //cout << "n,names[n] " << n << " " << names[n];
    for (int i=0;i<int(colours_name_map.size());i++) {
      if (names[n]==colours_name_map[i]) {
        codes[n]=i;
        //cout << " code " << codes[n] << endl;
        break;
      }
    }
  }
  return codes;
}

std::vector<double> RGBReps::GetColour(unsigned int i){
  /* Return colours in order of input file */
  if(i>=colours.size()){
    //cerr << "RGBReps::GetColour Warning colour out of range, returning default " << i << " " << GetNumberOfColours() << "\n";
    std::vector<double> col(4); 
    col[0] = col[1] = col[2] = col[3] = 1.0;
    return col;
  }

  return colours[i];
}

std::vector<std::vector<double> > RGBReps::GetColourDefinitions(std::vector<int>codes,int mode) { 
  int ii;
  std::vector <std::vector<double> > colvect;
  std::vector<double> col; 
  int n_codes =  codes.size();
  for (int i=0;i< n_codes;i++) {
      col.clear();
      ii = codes[i];
      if ( mode == HSVCOLOURREP ) {
         if (ii>=int(hsv.size())){
	   col.push_back(0.0);
	   col.push_back(0.0);
	   col.push_back(1.0);
	   col.push_back(1.0);
         } else {
           col.push_back(hsv[ii][0]);
           col.push_back(hsv[ii][1]);
           col.push_back(hsv[ii][2]);
           col.push_back(1.0);
         }
      } else {
         if (ii>=int(colours.size())){
	   col.push_back(1.0);
	   col.push_back(1.0);
	   col.push_back(1.0);
	   col.push_back(1.0);
         } else {
           col.push_back(colours[ii][0]);
           col.push_back(colours[ii][1]);
           col.push_back(colours[ii][2]);
           col.push_back(1.0);
         }
      }
      //std::cout << i << " " << ii << " " << array[i][0] << " " << array[i][1] << " " << array[i][2] << endl;
      //std::cout << "GetColourDefinitions " << col[0] << " " << col[1] << " " << col[2] << " " << std::endl; 
    colvect.push_back(col);
  }
  return colvect;
}

double* RGBReps::GetColourP(unsigned int i){
  /* Return colours in order of input file */
  double* col = new double[4]; 
  if(i>=colours.size()){
    //cerr << "RGBReps::GetColourP Warning colour out of range, returning default " << i << " " << GetNumberOfColours() << "\n";
    col[0] = col[1] = col[2] = col[3] = 1.0;
    return col;
  }
  col[0] = colours[i][0];
  col[1] = colours[i][1];
  col[2] = colours[i][2];
  col[3] = colours[i][3];
  return col;
}

std::vector<double> RGBReps::GetColourHSV(unsigned int i){
  /* Return colours in order of input file */
  std::vector <double> col(4); 
  if(i>=colours.size()){
    //cerr << "RGBReps::GetColourP Warning colour out of range, returning default " << i << " " << GetNumberOfColours() << "\n";
    col[0] = col[1] = col[2] = col[3] = 1.0;
    return col;
  }
  col[0] = hsv[i][0];
  col[1] = hsv[i][1];
  col[2] = hsv[i][2];
  col[3] = hsv[i][3];
  return col;
}

std::vector<std::vector<double> > RGBReps::GetColours(){
  /* Return colours in order of input file */
  return colours;
}

std::vector<double> RGBReps::GetColour(const std::string &name){
  /* Return colours BY NAME */
  
  if(string_colours[name].size()==0){
    if (name == "complement" ) {
      std::vector<double> col(4);
      col = complement_colour; 
      return col;      
    } else {
      //cerr << "RGBReps::GetColour(const std::string &name) Warning unknown colour name, returning default " << name << " " << GetNumberOfColours() << "\n";
      std::vector<double> col(4); 
      col[0] = col[1] = col[2] = col[3] = 1.0;
      return col;
    }
  }

  return string_colours[name];
}

RGBReps::RGBReps() {
  std::vector<double> comp_col(4);
  for (int n=0;n<4;n++)comp_col[n] = 1.0;
  RGBReps::SetBackgroundColour(comp_col);
}  

RGBReps::RGBReps(const char *infile){
  std::vector<double> comp_col(4);
  for (int n=0;n<4;n++)comp_col[n] = 1.0;
  RGBReps::SetBackgroundColour(comp_col);
  SetFile(infile);
}
RGBReps::RGBReps(std::string &infile){
  std::vector<double> comp_col(4);
  for (int n=0;n<4;n++)comp_col[n] = 1.0;
  RGBReps::SetBackgroundColour(comp_col);
  SetFile(infile);
}

RGBReps::RGBReps(std::ifstream &infile){
  std::vector<double> comp_col(4);
  for (int n=0;n<4;n++)comp_col[n] = 1.0;
  RGBReps::SetBackgroundColour(comp_col);
  SetFile(infile);
}

void RGBReps::SetFile(const char *infile){
  std::ifstream s(infile);
  RGBReps::BuildRGBReps(s);
}
void RGBReps::SetFile(std::string &infile){
  std::ifstream s(infile.c_str());
  RGBReps::BuildRGBReps(s);
}

void RGBReps::SetFile(std::ifstream &infile){
  RGBReps::BuildRGBReps(infile);
}

void RGBReps::SetBackgroundColour(std::string colour) {
  std::vector<double> comp_col(4);
  std::vector<double> col(4);
  col = RGBReps::GetColour(colour);
  for (int n=0;n<3;n++)comp_col[n] = 1.0 - col[n];
  comp_col[3] = 1.0;
  complement_colour = comp_col;
}
void RGBReps::SetBackgroundColour(std::vector<double> col) {
  std::vector<double> comp_col(4);
  for (int n=0;n<3;n++)comp_col[n] = 1.0 - col[n];
  comp_col[3] = 1.0;
  complement_colour = comp_col;
}

void RGBReps::DeleteColour(string colour){
  if(string_colours[colour].size()==0){
    //cerr << "RGBReps::DeleteColour Warning unknown colour name, deleting nothing\n";
    return;
  }

  string_colours.erase(colour);
  
  std::vector<std::string> colours_name_map_tmp;
  std::vector<std::vector<double> > colours_tmp;

  vector<std::string>::const_iterator col_name_iter=colours_name_map.begin();
  std::vector<std::vector<double> >::const_iterator col_iter = colours.begin();

  while(col_name_iter!=colours_name_map.end()){
    if(colour!=*col_name_iter){
       colours_name_map_tmp.push_back(*col_name_iter);
       colours_tmp.push_back(*col_iter);
    }
    col_name_iter++;
    col_iter++;
  }

  colours_name_map = colours_name_map_tmp;
  colours = colours_tmp;

}

void RGBReps::ChangeColour(string colour, double r, double g, double b, double a){
  if(string_colours[colour].size()!=0){
    AddColour(colour,r,g,b,a);
  }else{
    cerr << "Colour " <<  colour << " does not exist cannot change\n";
  }
}

void RGBReps::AddColour(string colour, double r, double g, double b, double a){
  vector<double> rgba(4);
  vector<double> hsva;
  rgba[0] = r;
  rgba[1] = g;
  rgba[2] = b;
  rgba[3] = a;
 
  hsva = rgbtohsv(rgba);
  //std::cout << "rgba " << rgba[0] << " "  << rgba[1] << " "   << rgba[2] << endl;
  //std::cout << "hsva " << hsva[0] << " "  << hsva[1] << " "   << hsva[2] << endl;

  if(string_colours[colour].size()==0){
    colours.push_back(rgba);
    hsv.push_back(hsva);
    //std::cout << " " << rgba[0] << " " <<  rgba[1] << " " <<  rgba[2] << " " << std::endl;
    //std::cout << "hsva " << hsva[0] << " " << hsva[1] << " " <<  hsva[2] << " " << std::endl;
    colours_name_map.push_back(colour);
  }else{
    int i = 0;
    vector<std::string>::const_iterator col_name_iter=colours_name_map.begin();
    while(col_name_iter!=colours_name_map.end()){
      if(colour==*col_name_iter){
	colours[i] = rgba; 
	hsv[i] = hsva;
      }
      col_name_iter++; i++;
    }
  }
  string_colours[colour] = rgba;
}

void RGBReps::BuildRGBReps(ifstream &infile){

  string colour;
  string red, green, blue, alpha;
  char buf[1024];
  string word;

  colours.clear();
  string_colours.clear();

  vector<double> rgba(4);

  while(infile.getline(buf,1024)){
    word = string(buf);
    colour = word.substr(0,word.find(','));
    word = word.substr(word.find(',')+1);
    red = word.substr(0,word.find(','));
    word = word.substr(word.find(',')+1);
    green = word.substr(0,word.find(','));
    word = word.substr(word.find(',')+1);
    blue = word.substr(0,word.find(','));
    word = word.substr(word.find(',')+1);
    alpha = word.substr(0,word.find(','));
    
    rgba[0] = atof(red.c_str());
    rgba[1] = atof(green.c_str());
    rgba[2] = atof(blue.c_str());
    rgba[3] = atof(alpha.c_str());
    
    colours.push_back(rgba);
    vector<double> hsva = rgbtohsv(rgba);
    hsv.push_back(hsva);
    string_colours[colour] = rgba;
    colours_name_map.push_back(colour);

  }

  if(string_colours["default"].size() == 0){
    rgba[0] = rgba[1] = rgba[2] = rgba[3] = 1.0;
    colours.push_back(rgba);
    string_colours["default"]  = rgba;
    colours_name_map.push_back("default");
  }

}

//-------------------------------------------------------------------------
std::vector<double> RGBReps::rgbtohsv ( std::vector<double> rgb) {
  //-------------------------------------------------------------------------
  double s,h,v;
  vector<double> hsv(4);

  hsv[3] = rgb[3];

  double maxrgb = max(rgb[0],rgb[1]);
         maxrgb = max(maxrgb,rgb[2]);
  double minrgb = min(rgb[0],rgb[1]);
         minrgb = min(minrgb,rgb[2]);
  
  v = maxrgb;
  double delta = maxrgb - minrgb;

  if (fabs(delta)<1e-6) {
    hsv[0] = 0;
    hsv[1] = 0;
    hsv[2] = v;
    return hsv;
  }
 
  if ( fabs(maxrgb) > 1e-6 ) {
    s = delta / maxrgb;
  } else {
    hsv[0] = 0;
    hsv[1] = -1;
    hsv[2] = v;
    return hsv;
  }

  if ( fabs(rgb[0]-maxrgb) < 1e-6 ) {
    h = float(rgb[1] - rgb[2]) /delta;
  } else if (fabs(rgb[1] - maxrgb)< 1e-6 ) {
    h= 2 + (float( rgb[2] - rgb[0] )/delta);
  } else {
    h = 4 + (float( rgb[0] - rgb[1] )  / delta);
  }

  h= h * 60;
  if (h < 0 ) h = h + 360;

  hsv[0] = h;
  hsv[1] = s;
  hsv[2] = v;
  return hsv;

}

//-----------------------------------------------------------------------
std::vector<double> RGBReps::hsvtorgb (std::vector<double> hsv) {
//-----------------------------------------------------------------------
  double r,g,b;
  std::vector<double> rgb(4);


  rgb[3] = hsv[3];


  if ( fabs(hsv[0]) < 1e-6 && fabs(hsv[1]) < 1e-6 && fabs(hsv[2]) < 1e-6 ) {
    rgb[0] = 1.0; 
    rgb[1] = 1.0; 
    rgb[2] = 1.0; 
    return rgb;
  }
  if ( fabs(hsv[1]) < 1e-6 ) {
    rgb[0] = rgb[1] = rgb[2] = hsv[2];
    return rgb;
  }

  double h = hsv[0] / 60.0;
  int i = int(floor(hsv[0]/60.0)); // Need to round up
  double f = h - i;
  double p = hsv[2] * (1.0 - hsv[1]);
  double q = hsv[2] * (1.0 - hsv[1] * f);
  double t = hsv[2] * (1.0 - hsv[1] * (1.0 - f));
  double v = hsv[2];

  switch (i) {
    case 0 :
        r = v;
        g  = t;
        b  = p;
        break;
    case  1 : 
        r  = q;
        g  = v;
        b  = p;
        break;
    case 2:
        r = p;
        g = v;
        b = t;
        break;
    case 3 :
        r  = p;
        g = q;
        b = v;
        break;
    case 4 :
        r = t;
        g = p;
        b = v;
        break;
    default :
        r = v;
        g = p;
        b = q;
        break;
    }
    rgb[0] = r;
    rgb[1] = g;
    rgb[2] = b;

    return rgb;
}

//-------------------------------------------------------------------------
double *RGBReps::rgbtohsv ( double *rgb) {
//-------------------------------------------------------------------------
  std::vector<double> rgbvec(4);
  rgbvec[0] = rgb[0];
  rgbvec[1] = rgb[1];
  rgbvec[2] = rgb[2];
  std::vector<double> hsvvec = rgbtohsv(rgbvec);
  double * hsv = new double[4];
  hsv[0] = hsvvec[0];
  hsv[1] = hsvvec[1];
  hsv[2] = hsvvec[2];
  hsv[3] = rgb[3];
  return hsv;

}

//-----------------------------------------------------------------------
double *RGBReps::hsvtorgb (double *hsv) {
//-----------------------------------------------------------------------
  std::vector<double> hsvvec(4);
  hsvvec[0] = hsv[0];
  hsvvec[1] = hsv[1];
  hsvvec[2] = hsv[2];
  std::vector<double> rgbvec = hsvtorgb(hsvvec);
  double *rgb = new double[4];
  rgb[0] = rgbvec[0];
  rgb[1] = rgbvec[1];
  rgb[2] = rgbvec[2];
  rgb[3] = hsv[3];
  return rgb;
}


//-----------------------------------------------------------------------
void RGBReps::clear(){
//-----------------------------------------------------------------------
  colours_name_map.clear();
  string_colours.clear();
  colours.clear();
  hsv.clear();
}

std::vector<std::string> RGBReps::colours_name_map;
std::map<std::string,std::vector<double> > RGBReps::string_colours;
std::vector<std::vector<double> > RGBReps::colours;
std::vector<std::vector<double> > RGBReps::hsv;
std::vector<double> RGBReps::complement_colour;
