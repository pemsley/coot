//   CCP4 Molecular Graphics Program
//
//   Copyright 2004 The University of York
//   Author: Stuart McNicholas and Liz Potterton
//
//   This program is free software and is distributed under the terms
//   and conditions of the CCP4 licence agreement as `Part 0' (Annex 2)
//   software, which is version 2.1 of the GNU Lesser General Public
//   Licence (LGPL) with the following additional clause:
//
//      `You may also combine or link a "work that uses the Library"
//      to produce a work containing portions of the Library, and
//      distribute that work under terms of your choice, provided that
//      you give prominent notice with each copy of the work that the
//      specified version of the Library is used in it, and that you
//      include or provide public access to the complete corresponding
//      machine-readable source code for the Library including whatever
//      changes were used in the work. (i.e. If you make changes to the
//      Library you must distribute those, but you do not need to
//      distribute source or object code to those portions of the work
//      not covered by this licence.)'
//
//   Note that this clause grants an additional right and does not
//   impose any additional restriction, and so does not affect
//   compatibility with the GNU General Public Licence (GPL). If you
//   wish to negotiate other terms, please contact the maintainer.
//   You can redistribute it and/or modify the program under the terms
//   of the GNU Lesser General Public License as published by the Free
//   Software Foundation; either version 2.1 of the License, or (at
//   your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//   Lesser General Public License for more details.
//
//   You should have received a copy of the CCP4 licence and/or GNU
//   Lesser General Public License along with this program; if not,
//   write to the CCP4 Secretary, Daresbury Laboratory, Warrington
//   WA4 4AD, UK. The GNU Lesser General Public can also be obtained
//   by writing to the Free Software Foundation, Inc., 59 Temple Place,
//   Suite 330, Boston, MA 02111-1307 USA


#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>

#include "rgbreps.h"


unsigned int RGBReps::GetNumberOfColours(){
  return colours.size();
}

std::string RGBReps::GetNameOfColour(unsigned int i){
  /* Return name of colour as in order of input file */
  if(i>=colours.size()){
     std::cerr << "ERROR:: colour out of range in GetNameOfColour, returning default "
	       << i << " " << GetNumberOfColours() << "\n";
    return std::string("default");
  }
  return colours_name_map[i];
}

int  RGBReps::GetColourNumber(std::string name ) {
  for (int i=0;i<colours_name_map.size();i++) {
    if (!name.compare(colours_name_map[i])) {
      return i;
    }
  }
  return -1;
}

std::vector<int>  RGBReps::GetColourNumbers(std::vector<std::string> names ) {
  std::vector<int> codes(names.size());
  for (unsigned int n=0;n<names.size();n++) {
    //cout << "n,names[n] " << n << " " << names[n];
    for (unsigned int i=0; i<colours_name_map.size(); i++) {
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
     std::cerr << "ERROR in RGBReps::GetColour Warning colour out of range, returning default "
	       << i << " " << GetNumberOfColours() << "\n";
     std::vector<double> col(4); 
     col[0] = col[1] = col[2] = col[3] = 255;
    return col;
  }

  return colours[i];
}

std::vector<std::vector<double> > RGBReps::GetColourDefinitions(std::vector<int>codes,int mode) { 
  int ii;
  std::vector <std::vector<double> > colvect;
  std::vector<double> col; 
  int n_codes =  codes.size();
  for (unsigned int i=0;i< n_codes;i++) {
      col.clear();
      ii = codes[i];
      if (ii>colours.size()){
	col.push_back(255.0);
	col.push_back(255.0);
	col.push_back(255.0);
	col.push_back(255.0);
      } else if ( mode == HSVCOLOURREP ) {
         col.push_back(hsv[ii][0]);
         col.push_back(hsv[ii][1]);
         col.push_back(hsv[ii][2]);
         col.push_back(255.0);
      } else {
         col.push_back(colours[ii][0]);
         col.push_back(colours[ii][1]);
         col.push_back(colours[ii][2]);
         col.push_back(255.0);
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
     std::cerr << "RGBReps::GetColourP Warning colour out of range, returning default " << i << " " << GetNumberOfColours() << "\n";
    col[0] = col[1] = col[2] = col[3] = 255;
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
     std::cerr << "RGBReps::GetColourHSV Warning colour out of range, returning default "
	       << i << " " << GetNumberOfColours() << "\n";
     col[0] = col[1] = col[2] = col[3] = 255;
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
    std::cerr << "RGBReps::GetColour(const std::string &name) Warning unknown colour name, returning default " << name << " " << GetNumberOfColours() << "\n";
    std::vector<double> col(4); 
    col[0] = col[1] = col[2] = col[3] = 255;
    return col;
  }

  return string_colours[name];
}

RGBReps::RGBReps(const char *infile){
  SetFile(infile);
}
RGBReps::RGBReps(std::string &infile){
  SetFile(infile);
}

RGBReps::RGBReps(std::ifstream &infile){
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

void RGBReps::DeleteColour(std::string colour){
  if(string_colours[colour].size()==0){
    std::cerr << "RGBReps::DeleteColour Warning unknown colour name, deleting nothing\n";
    return;
  }

  string_colours.erase(colour);
  
  std::vector<std::string> colours_name_map_tmp;
  std::vector<std::vector<double> > colours_tmp;

  std::vector<std::string>::const_iterator col_name_iter=colours_name_map.begin();
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

void RGBReps::ChangeColour(std::string colour, double r, double g, double b, double a){
  if(string_colours[colour].size()!=0){
    AddColour(colour,r,g,b,a);
  }else{
    std::cerr << "Colour " <<  colour << " does not exist cannot change\n";
  }
}

void RGBReps::AddColour(std::string colour, double r, double g, double b, double a){
  std::vector<double> rgba(4);
  std::vector<double> hsva;
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
    colours_name_map.push_back(colour);
  }else{
    int i = 0;
    std::vector<std::string>::const_iterator col_name_iter=colours_name_map.begin();
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

void RGBReps::BuildRGBReps(std::ifstream &infile){

  std::string colour;
  std::string red, green, blue, alpha;
  char buf[1024];
  std::string word;

  colours.clear();
  string_colours.clear();

  std::vector<double> rgba(4);

  while(infile.getline(buf,1024)){
    word = std::string(buf);
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
    std::vector<double> hsva = rgbtohsv(rgba);
    hsv.push_back(hsva);
    string_colours[colour] = rgba;
    colours_name_map.push_back(colour);

  }

  if(string_colours["default"].size() == 0){
    rgba[0] = rgba[1] = rgba[2] = rgba[3] = 255.0;
    colours.push_back(rgba);
    string_colours["default"]  = rgba;
    colours_name_map.push_back("default");
  }

}

//-------------------------------------------------------------------------
std::vector<double> RGBReps::rgbtohsv ( std::vector<double> rgb) {
  //-------------------------------------------------------------------------
  double s,h,v;
  std::vector<double> hsv(4);

  hsv[3] = rgb[3];

  double maxrgb = std::max(rgb[0],rgb[1]);
         maxrgb = std::max(maxrgb,rgb[2]);
  double minrgb = std::min(rgb[0],rgb[1]);
         minrgb = std::min(minrgb,rgb[2]);
  
  v = maxrgb;
  double delta = maxrgb - minrgb;

  if (delta == 0) {
    hsv[0] = 0;
    hsv[1] = 0;
    hsv[2] = floor(v);
    return hsv;
  }
 
  if ( maxrgb > 0 ) {
    s = delta * 255.0 / maxrgb;
  } else {
    hsv[0] = 0;
    hsv[1] = -1;
    hsv[2] = floor(v);
    return hsv;
  }
 
  if ( rgb[0] == maxrgb ) {
    h = float(rgb[1] - rgb[2]) /delta;
  } else if (rgb[1] == maxrgb ) {
    h= 2 + (float( rgb[2] - rgb[0] )/delta);
  } else {
    h = 4 + (float( rgb[0] - rgb[1] )  / delta);
  }

  h= h * 60;
  if (h < 0 ) h = h + 360;

  hsv[0] = floor(h);
  hsv[1] = floor(s);
  hsv[2] = floor(v);
  return hsv;

}

//-----------------------------------------------------------------------
std::vector<double> RGBReps::hsvtorgb (std::vector<double> hsv) {
  //-----------------------------------------------------------------------
  double r,g,b;
  std::vector<double> rgb(4);

  rgb[3] = hsv[3];

  if ( hsv[0] == 0 && hsv[1] == 0 && hsv[2] == 0 ) {
    rgb[0] =255; 
    rgb[1] =255; 
    rgb[2] =255; 
    return rgb;
  }
  if ( hsv[1] == 0 ) {
    rgb[0] = rgb[1] = rgb[2] = floor(hsv[2]);
    return rgb;
  }

  double h = hsv[0] / 60.0;
  int i = int(floor(hsv[0]/60.0)); // Need to round up
  double f = h - i;
  double p = (float(hsv[2])/255.0) * (1.0 - float(hsv[1])/255.0);
  double q = (float(hsv[2])/255.0) * (1.0 - (float(hsv[1])/255.0) * f);
  double t = (float(hsv[2])/255.0) * (1.0 - (float(hsv[1])/255.0) * (1.0 - f));
  double v = hsv[2];

  switch (i) {
    case 0 :
        r = v;
        g  = t * 255;
        b  = p * 255;
        break;
    case  1 : 
        r  = q * 255;
        g  = v;
        b  = p * 255;
        break;
    case 2:
        r = p * 255;
        g = v;
        b = t * 255;
        break;
    case 3 :
        r  = p * 255;
        g = q * 255;
        b = v;
        break;
    case 4 :
        r = t * 255;
        g = p * 255;
        b = v;
        break;
    default :
        r = v;
        g = p * 255;
        b = q * 255;
        break;
    }
    rgb[0] = floor(r);
    rgb[1] = floor(g);
    rgb[2] = floor(b);
    return rgb;
}

//-------------------------------------------------------------------------
double *RGBReps::rgbtohsv ( double *rgb) {
  //-------------------------------------------------------------------------
  double s,h,v;
  //vector<double> hsv(4);
  double *hsv = new double[4];
  hsv[3] = rgb[3];

  double maxrgb = std::max(rgb[0],rgb[1]);
         maxrgb = std::max(maxrgb,rgb[2]);
  double minrgb = std::min(rgb[0],rgb[1]);
         minrgb = std::min(minrgb,rgb[2]);
  
  v = maxrgb;
  double delta = maxrgb - minrgb;

  if (delta == 0) {
    hsv[0] = 0;
    hsv[1] = 0;
    hsv[2] = floor(v);
    return hsv;
  }
 
  if ( maxrgb > 0 ) {
    s = delta * 255.0 / maxrgb;
  } else {
    hsv[0] = 0;
    hsv[1] = -1;
    hsv[2] = floor(v);
    return hsv;
  }
 
  if ( rgb[0] == maxrgb ) {
    h = float(rgb[1] - rgb[2]) /delta;
  } else if (rgb[1] == maxrgb ) {
    h= 2 + (float( rgb[2] - rgb[0] )/delta);
  } else {
    h = 4 + (float( rgb[0] - rgb[1] )  / delta);
  }

  h= h * 60;
  if (h < 0 ) h = h + 360;

  hsv[0] = floor(h);
  hsv[1] = floor(s);
  hsv[2] = floor(v);
  return hsv;

}

//-----------------------------------------------------------------------
double *RGBReps::hsvtorgb (double *hsv) {
  //-----------------------------------------------------------------------
  double r,g,b;
  //std::vector<double> rgb(4);
  double * rgb = new double [4];
  rgb[3] = hsv[3];

  if ( hsv[0] == 0 && hsv[1] == 0 && hsv[2] == 0 ) {
    rgb[0] =255; 
    rgb[1] =255; 
    rgb[2] =255; 
    return rgb;
  }
  if ( hsv[1] == 0 ) {
    rgb[0] = rgb[1] = rgb[2] = floor(hsv[2]);
    return rgb;
  }

  double h = hsv[0] / 60.0;
  int i = int(ceil(hsv[0]/60.0)); // Need to round up
  double f = h - i;
  double p = (float(hsv[2])/255.0) * (1.0 - float(hsv[1])/255.0);
  double q = (float(hsv[2])/255.0) * (1.0 - (float(hsv[1])/255.0) * f);
  double t = (float(hsv[2])/255.0) * (1.0 - (float(hsv[1])/255.0) * (1.0 - f));
  double v = hsv[2];

  switch (i) {
    case 0 :
        r = v;
        g  = t * 255;
        b  = p * 255;
        break;
    case  1 : 
        r  = q * 255;
        g  = v;
        b  = p * 255;
        break;
    case 2:
        r = p * 255;
        g = v;
        b = t * 255;
        break;
    case 3 :
        r  = p * 255;
        g = q * 255;
        b = v;
        break;
    case 4 :
        r = t * 255;
        g = p * 255;
        b = v;
        break;
    default :
        r = v;
        g = p * 255;
        b = q * 255;
        break;
    }
    rgb[0] = floor(r);
    rgb[1] = floor(g);
    rgb[2] = floor(b);
    return rgb;
}



std::vector<std::string> RGBReps::colours_name_map;
std::map<std::string,std::vector<double> > RGBReps::string_colours;
std::vector<std::vector<double> > RGBReps::colours;
std::vector<std::vector<double> > RGBReps::hsv;
