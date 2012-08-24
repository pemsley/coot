/*
     pygl/rgbreps.h: CCP4MG Molecular Graphics Program
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


#ifndef _CCP4MG_RGBREPS_
#define _CCP4MG_RGBREPS_


#include <map>
#include <vector>
#include <string>
#include <fstream>

enum ColourRepresentation { NOCOLOURREP , RGBCOLOURREP , HSVCOLOURREP };

class RGBReps {
  static std::vector<std::string> colours_name_map;
  static std::map<std::string,std::vector<double> > string_colours;
  static std::vector<std::vector<double> > colours;
  static std::vector<std::vector<double> >hsv;
  static std::vector<double> complement_colour;

  static void BuildRGBReps(std::ifstream &infile);
 public:
  RGBReps(std::ifstream &infile);
  RGBReps(std::string &infile);
  RGBReps(const char *infile);
  RGBReps();
  ~RGBReps(){};

  static void SetFile(std::ifstream &infile);
  static void SetFile(std::string &infile);
  static void SetFile(const char *infile);

  static std::vector<std::vector<double> > GetColours();
  static std::vector<double> GetColour(unsigned int i);
  static double* GetColourP(unsigned int i);
  static std::vector<double> GetColourHSV(unsigned int i);
  static std::vector<std::vector<double> > GetColourDefinitions( std::vector<int> codes,int mode);
  static std::string GetNameOfColour(unsigned int i);
  static std::vector<double> GetColour(const std::string &name);
  static unsigned int GetNumberOfColours();
  static void AddColour(std::string colour, double r, double g, double b, double a=255.0);
  static void ChangeColour(std::string colour, double r, double g, double b, double a=255.0);
  static void DeleteColour(std::string colour);
  static int GetColourNumber(std::string name ); 
  static std::vector<int> GetColourNumbers(std::vector<std::string> names ); 
  static std::vector<double> hsvtorgb (std::vector<double> hsv);
  static std::vector<double> rgbtohsv (std::vector<double> rgb );
  static double* hsvtorgb (double *hsv);
  static double* rgbtohsv (double *rgb );
  static void SetBackgroundColour(std::string colour);
  static void SetBackgroundColour(std::vector<double> col);
};

#endif
